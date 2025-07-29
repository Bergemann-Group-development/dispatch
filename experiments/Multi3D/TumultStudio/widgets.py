import ipywidgets as widgets
from traitlets import Unicode, observe, List, HasTraits, Instance, Integer
from IPython.display import display, Javascript
from io import StringIO
import pandas as pd
import uuid
import time
import os
from .auxfuncs import replace_nth

class StyleSheet(widgets.HTML):
    def __init__(self, file):
        super().__init__()
        if os.path.exists(file):
            with open(file, 'r') as f:
                css = f.read()
                self.value = f'<style>{css}</style>'
        else:
            print(f"Warning: {file} not found.")

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


class EditableList(widgets.VBox, HasTraits):
    """A list widget that displays objects as EditableLabel widgets."""
    
    objects = List()
    selected = Instance(widgets.HBox, allow_none=True)
    
    def __init__(self, objects=None, **kwargs):
        super().__init__(**kwargs)
        self.add_class('EditableList')

        self.selected = None
        self.list_items = []  # Stores HBox containers with label + delete button
        self.label_to_object = {}  # Maps EditableLabel to original object
        self.objects = objects or []
        
        # Initial population
        self.update_list()
    
    @observe('objects')
    def _on_objects_change(self, change):
        """Update the list when objects change."""
        self.update_list_incremental()
    
    def update_list(self):
        """Full update - recreate all widgets (for initial setup)."""
        # Clear existing widgets
        self.children = []
        self.list_items.clear()
        self.label_to_object.clear()
        
        # Create new widgets for all objects
        for obj in self.objects:
            self._create_item_widget(obj)
        
        # Update children
        self.children = self.list_items
    
    def update_list_incremental(self):
        """Incremental update - only modify what changed."""
        # Get current objects mapped by id/reference
        current_objects = {id(obj): obj for obj in self.objects}
        existing_objects = {id(obj): obj for obj in self.label_to_object.values()}
        
        # Find objects to remove
        objects_to_remove = set(existing_objects.keys()) - set(current_objects.keys())
        
        # Find objects to add
        objects_to_add = set(current_objects.keys()) - set(existing_objects.keys())
        
        # Remove widgets for deleted objects
        for obj_id in objects_to_remove:
            obj = existing_objects[obj_id]
            self._remove_item_widget(obj)
        
        # Add widgets for new objects
        for obj_id in objects_to_add:
            obj = current_objects[obj_id]
            self._create_item_widget(obj)
        
        # Update children order to match objects order
        self._reorder_children()
        
        # Update children
        self.children = self.list_items
    
    def _create_item_widget(self, obj):
        """Create a new item widget for an object."""
        # Create label widget
        label_widget = EditableLabel(value=obj.label)
        
        # Create delete button
        delete_button = widgets.Button(
            description='×',
            tooltip='Delete item'
        )
        delete_button.add_class('del-btn')
        
        # Create container for label + delete button
        item_container = widgets.HBox([
            label_widget,
            delete_button
        ])
        item_container.add_class('row-box')
        item_container.object = obj
        
        # Store mapping
        self.label_to_object[label_widget] = obj
        self.list_items.append(item_container)
        
        # Setup bindings
        self._setup_binding(label_widget, obj)
        self._setup_delete_button(delete_button, obj)

        def select_container(btn):
            self.selected = item_container

        def delete_selection(btn):
            if self.selected is item_container:
                self.selected = None

        label_widget.label.on_click(select_container)
        delete_button.on_click(delete_selection)

    @observe('selected')
    def _update_selected(self, change):
        selected = change['new']
        for cnt in self.list_items:
            if cnt is selected:
                cnt.add_class('selected')
            else:
                cnt.remove_class('selected')
    
    def _remove_item_widget(self, obj):
        """Remove the widget for an object."""
        # Find the label widget for this object
        label_widget = None
        for label, mapped_obj in self.label_to_object.items():
            if mapped_obj is obj:
                label_widget = label
                break
        
        if label_widget:
            # Find the container holding this label
            container_to_remove = None
            for container in self.list_items:
                if container.children[0] is label_widget:  # First child is the label
                    container_to_remove = container
                    break
            
            if container_to_remove:
                # Remove from list and mapping
                self.list_items.remove(container_to_remove)
                del self.label_to_object[label_widget]
    
    def _reorder_children(self):
        """Reorder children to match the objects order."""
        # Create a mapping from object to container
        obj_to_container = {}
        for container in self.list_items:
            label_widget = container.children[0]  # First child is the label
            obj = self.label_to_object[label_widget]
            obj_to_container[obj] = container
        
        # Reorder list_items to match objects order
        self.list_items = [obj_to_container[obj] for obj in self.objects if obj in obj_to_container]
    
    def _setup_binding(self, label_widget, obj):
        """Setup two-way binding between EditableLabel and object."""
        
        def on_label_change(change):
            # Update object when label changes
            obj.label = change['new']
        
        def on_object_change(change):
            # Update label when object changes
            label_widget.value = change['new']
        
        # Bind label changes to object
        label_widget.observe(on_label_change, 'description')
        
        # Bind object changes to label (if object supports observation)
        if hasattr(obj, 'observe'):
            obj.observe(on_object_change, 'label')
    
    def _setup_delete_button(self, delete_button, obj):
        """Setup delete button click handler."""
        
        def on_delete_click(button):
            
            # Call object's delete method if available
            if hasattr(obj, 'delete') and callable(getattr(obj, 'delete')):
                try:
                    obj.delete()
                except Exception as e:
                    print(f"Error calling delete method on object: {e}")
            
            # Remove object from list
            self.remove_object(obj)
        
        delete_button.on_click(on_delete_click)
    
    def add_object(self, obj, label=None):
        """Add a new object to the list."""
        if label: obj.__dict__['label'] = label
        self.objects = self.objects + [obj]  # Trigger trait change
    
    def remove_object(self, obj):
        """Remove an object from the list."""
        if obj in self.objects:
            new_objects = [o for o in self.objects if o != obj]
            self.objects = new_objects
    
    def clear(self):
        """Clear all objects from the list."""
        self.objects = []
    
    def get_object_by_label(self, label_widget):
        """Get the original object associated with a label widget."""
        return self.label_to_object.get(label_widget)
    
    def set_objects(self, objects_list):
        """Set the list of objects."""
        self.objects = objects_list



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

        # Create widgets
        self.html_widget = widgets.HTML(value=html_table)
        self.html_widget.add_class('jp-RenderedHTMLCommon')

        self.export_widget = widgets.Textarea(layout={'display': 'none'})
        self.export_widget.add_class(self.tid_class)

        # JavaScript output widget
        self.js_output = widgets.Output(layout={'display': 'none'})

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
        new = replace_nth(new, '<tr>', '<tr class="highlighted-row">', idx+1)
        self.value = new
        # The MutationObserver will automatically handle scrolling!

    def get_dataframe(self):
        """Returns the current state of the table as a DataFrame"""
        try:
            html = StringIO(self.value)
            return pd.read_html(html)[0]
        except Exception as e:
            print("⚠️ Failed to parse table HTML:", e)
            return self.df.copy()