import ipywidgets as widgets
from traitlets import Unicode, observe, List, HasTraits, Instance, Integer
from IPython.display import display, Javascript
from io import StringIO
import pandas as pd
import uuid
import time
import os
from .auxfuncs import *

class StyleSheet(widgets.HTML):
    def __init__(self, file):
        super().__init__()
        if os.path.exists(file):
            with open(file, 'r') as f:
                css = f.read()
                self.value = f'<style>{css}</style>'
        else:
            print(f"Warning: {file} not found.")

class SpinnerWidget(widgets.HTML):
    """A simple spinner icon that can be used to indicate loading or processing."""
    def __init__(self, icon, color='auto', style=None):
        super().__init__()
        self.add_class('spinner-widget')
        style = style or {}
        styling = ';'.join(f'{k}: {v}' for k, v in style.items()) + ';'
        icon_string = f'<i class="fa fa-{icon}" style="color: {color}; transform: rotate(0deg); {styling}"></i>'
        self.value = icon_string

    def spin(self):
        val = self.value
        rot = parse.search('rotate({rot:d}deg)', val)['rot']
        new_rot = (rot + 45) % 360
        self.value = val.replace(f'rotate({rot:d}deg)', f'rotate({new_rot:d}deg)')

    def reset(self):
        val = self.value
        rot = parse.search('rotate({rot:d}deg)', val)['rot']
        self.value = val.replace(f'rotate({rot:d}deg)', 'rotate(0deg)')

    @property
    def icon(self):
        return parse.search('class="fa fa-{icon}"', self.value)['icon']
    @icon.setter
    def icon(self, new_icon):
        self.value = self.value.replace(f'fa fa-{self.icon}', f'fa fa-{new_icon}')

    @property
    def color(self):
        return parse.search('color: {color};', self.value)['color']
    @color.setter
    def color(self, new_color):
        self.value = self.value.replace(f'color: {self.color};', f'color: {new_color};')


class EditableLabel(widgets.VBox):
    """A label that becomes editable when double-clicked."""
    
    value = Unicode('').tag(sync=True)
    
    def __init__(self, value='', **kwargs):
        self._last_click_time = 0
        self._double_click_threshold = 0.5  # seconds
        self.selected = False
        
        # Create widgets
        self.label = widgets.Button(
            description=value,
            disabled=False,
            button_style='',
            layout=widgets.Layout(border='none', width='auto')
        )
        self.editor = widgets.Text(value=value, layout=widgets.Layout(display='none'))
        self.value = value
        
        # Setup events
        self.label.on_click(self._handle_click)
        self.editor.on_submit(self._finish_edit)
        
        super().__init__([self.label, self.editor], **kwargs)
        self.add_class('EditableLabel')
    
    def _handle_click(self, button):
        """Handle single/double click detection."""
        current_time = time.time()
        
        if current_time - self._last_click_time < self._double_click_threshold:
            # Double click detected
            self._start_edit()
        else:
            self.selected = not self.selected
        
        self._last_click_time = current_time
    
    def _start_edit(self):
        """Switch to edit mode."""
        self.label.layout.display = 'none'
        self.editor.layout.display = 'block'
        self.editor.value = self.value
    
    def _finish_edit(self, change=None):
        """Switch back to label mode and update value."""
        self.value = self.editor.value
        self.label.description = self.value
        self.label.layout.display = 'block'
        self.editor.layout.display = 'none'
    
    @observe('value')
    def _on_value_change(self, change):
        """Update label when value changes."""
        self.label.description = change["new"]
        self.editor.value = change["new"]

class HTMLDataFrame(widgets.HTML):
    def __init__(self, df, **kwargs):
        """Create an HTML representation of a DataFrame."""
        super().__init__(**kwargs)
        self.df = df.reset_index(drop=True)
        self.value = self.df.to_html(index=False)
        self.add_class('HTMLDataFrame')
        self.add_class('jp-RenderedHTMLCommon')

class EditableTable(widgets.VBox, HasTraits):
    style = """
        <style>
        .EditableTable {
            width: fit-content;
        }
        .EditableTable .jp-RenderedHTMLCommon {
            padding: 0;
        }
        .EditableTable .dataframe {
            margin: 0;
        }
        .EditableTable .highlighted-row {
            background-color: var(--jp-brand-color1);
        }
        </style>
    """

    def __init__(self, df: pd.DataFrame, max_height='350px', **kwargs):
        self.df = df.reset_index(drop=True)
        self.sid = str(uuid.uuid4()).replace('-', '_')
        self.tid = f"editable_table_{self.sid}"
        self.tid_class = f"table_export_{self.sid}"
        self.container_class = f"EditableTable-{self.sid}"  # Unique container class

        # Create editable HTML table
        html_table = self.df.to_html(table_id=self.tid, **kwargs)
        html_table = html_table.replace('<td>', '<td contenteditable="true">')
        html_table += self.style
        self.initial_table = html_table
        self.highlighted_row = None

        # Create widgets
        self.html_widget = widgets.HTML(value=html_table)
        self.html_widget.add_class('jp-RenderedHTMLCommon')

        self.export_widget = widgets.Textarea(layout={'display': 'none'})
        self.export_widget.add_class(self.tid_class)

        # JavaScript output widget
        self.js_output = widgets.Output(layout={'display': 'none'})

        @log_callback_errors
        def sync_table(change):
            self._value = self.export_widget.value

        self.export_widget.observe(sync_table, names='value')

        super().__init__([self.html_widget, self.export_widget, self.js_output])
        self.layout.max_height = max_height
        self.add_class('EditableTable')
        self.add_class(self.container_class)  # Add unique class for targeting
        
        # Store the current value privately
        self._value = html_table
        
        # Set up event listener and mutation observer on this specific container
        self._setup_container_functionality()

    def _setup_container_functionality(self):
        """Set up both event listener and mutation observer for this container"""
        js = f"""
            function scrollToElementInContainer(element, container) {{
                    // Calculate the position of the element relative to the container
                    var elementRect = element.getBoundingClientRect();
                    var containerRect = container.getBoundingClientRect();
                    
                    // Calculate the scroll position needed to center the element
                    var relativeTop = elementRect.top - containerRect.top + container.scrollTop;
                    var centerPosition = relativeTop - (container.clientHeight / 2) + (element.offsetHeight / 2);
                    
                    // Smooth scroll to the calculated position
                    container.scrollTo({{
                        top: Math.max(0, centerPosition),
                        behavior: 'smooth'
                    }});
                }}
            setTimeout(function() {{
                // Find this specific container by unique class
                var container = document.querySelector('.{self.container_class}');
                
                if (!container) {{
                    console.warn("EditableTable: Container not found for {self.container_class}");
                    return;
                }}
                
                console.log("Setting up functionality for {self.container_class}");
                
                // 1. Set up focusout event listener for content sync
                container.addEventListener('focusout', function(event) {{
                    var target = event.target;
                    
                    // Check if this focusout is from a contenteditable cell
                    if (!target.contentEditable || target.contentEditable !== 'true') {{
                        return;
                    }}
                    
                    // Find the table within this container
                    var table = container.querySelector('table[id="{self.tid}"]');
                    if (!table) {{
                        return;
                    }}
                    
                    // Find the textarea within this container
                    var textarea = container.querySelector('.{self.tid_class} textarea');
                    if (!textarea) {{
                        console.warn("EditableTable: Textarea not found in {self.container_class}");
                        return;
                    }}
                    
                    // Sync the table content to textarea
                    setTimeout(function() {{
                        textarea.value = table.outerHTML;
                        textarea.dispatchEvent(new Event('input', {{ bubbles: true }}));
                        console.log("Table content synced for {self.tid}");
                    }}, 10);
                    
                }}, true); // Use capture phase
                
                // 2. Set up MutationObserver for highlighted row scrolling
                var observer = new MutationObserver(function(mutations) {{
                    mutations.forEach(function(mutation) {{
                        // Check for class changes on table rows
                        if (mutation.type === 'attributes' && mutation.attributeName === 'class') {{
                            var target = mutation.target;
                            
                            // Check if this is a table row with highlighted-row class
                            if (target.tagName === 'TR' && target.classList.contains('highlighted-row')) {{
                                console.log("Highlighted row detected, scrolling to view");
                                scrollToElementInContainer(target, container);
                            }}
                        }}
                        
                        // Also check for newly added nodes that might have highlighted-row class
                        if (mutation.type === 'childList') {{
                            mutation.addedNodes.forEach(function(node) {{
                                if (node.nodeType === Node.ELEMENT_NODE) {{
                                    // Check if the added node itself is a highlighted row
                                    if (node.tagName === 'TR' && node.classList.contains('highlighted-row')) {{
                                        console.log("New highlighted row added, scrolling to view");
                                        setTimeout(function() {{
                                            scrollToElementInContainer(node, container);
                                        }}, 100); // Small delay to ensure DOM is settled
                                    }}
                                    
                                    // Check for highlighted rows within the added node
                                    var highlightedRows = node.querySelectorAll && node.querySelectorAll('tr.highlighted-row');
                                    if (highlightedRows && highlightedRows.length > 0) {{
                                        console.log("Highlighted rows found in new content, scrolling to first one");
                                        setTimeout(function() {{
                                            scrollToElementInContainer(highlightedRows[0], container);
                                        }}, 100);
                                    }}
                                }}
                            }});
                        }}
                    }});
                }});
                
                // Start observing the container for changes
                observer.observe(container, {{
                    attributes: true,           // Watch for attribute changes
                    attributeFilter: ['class'], // Only watch class attribute changes
                    childList: true,           // Watch for added/removed nodes
                    subtree: true              // Watch all descendants
                }});
                
                // Store the observer reference for potential cleanup
                container.highlightObserver = observer;
                
                console.log("Event listener and MutationObserver attached to {self.container_class}");
                
            }}, 100); // Small delay to ensure DOM is ready
        """
        
        with self.js_output:
            display(Javascript(js))

    @property
    def value(self):
        return getattr(self, '_value', self.html_widget.value)
    
    @value.setter
    def value(self, new_value):
        self._value = new_value
        self.html_widget.value = new_value

    def reset_html(self):
        self.value = ''
        self.value = self.initial_table

    def highlight_row(self, idx):
        new = self.value
        new = new.replace('<tr class="highlighted-row">', '<tr>')
        if idx is not None:
            new = replace_nth(new, '<tr>', '<tr class="highlighted-row">', idx+1)
        self.value = new
        self.highlighted_row = idx
        # The MutationObserver will automatically handle scrolling!

    def get_dataframe(self):
        """Returns the current state of the table as a DataFrame"""
        try:
            html = StringIO(self.value)
            return pd.read_html(html)[0]
        except Exception as e:
            print("⚠️ Failed to parse table HTML:", e)
            return self.df.copy()