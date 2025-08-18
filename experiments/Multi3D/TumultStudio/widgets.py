import ipywidgets as widgets
from traitlets import Unicode, observe, List, HasTraits, Instance, Integer
from IPython.display import display, Javascript
from io import StringIO
import pandas as pd
import uuid
import time
import os
from .auxfuncs import *
import logging

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
    
    def __init__(self, objects=None, buttons=None, **kwargs):
        super().__init__(**kwargs)
        self.add_class('EditableList')
        self.buttons = buttons  # List of buttons to show for each item

        self.selected = None
        self.list_items = []  # Stores HBox containers with label + delete button
        self.label_to_object = {}  # Maps EditableLabel to original object
        self.object_to_container = {} # Maps object to its container widget
        self.container_pool = []  # Pool of reusable containers
        self.objects = objects or []
        self._observers = {}  # Track observers for cleanup
        
        # Initial population
        self.update_list()
    
    @observe('objects')
    def _on_objects_change(self, change):
        """Update the list when objects change."""
        self.update_list_incremental()
    
    def update_list(self):
        """Full update - recreate all widgets (for initial setup)."""
        # Clear existing widgets but keep containers in pool
        self.children = []
        self._clear_all_bindings()
        self.list_items.clear()
        
        # Ensure we have enough containers
        self._ensure_container_pool_size(len(self.objects))
        
        # Assign objects to containers
        for i, obj in enumerate(self.objects):
            self._assign_object_to_container(obj, self.container_pool[i])
        
        # Update list_items and children
        self.list_items = self.container_pool[:len(self.objects)]
        self.children = self.list_items
    
    @log_callback_errors
    def update_list_incremental(self):
        """Incremental update - reuse existing containers."""
        # Clear all current bindings
        self._clear_all_bindings()
        
        # Ensure we have enough containers for the new object count
        self._ensure_container_pool_size(len(self.objects))
        
        # Assign new objects to containers (reusing existing ones)
        for i, obj in enumerate(self.objects):
            container = self.container_pool[i]
            self._assign_object_to_container(obj, container)
        
        # Update list_items to match current objects
        self.list_items = self.container_pool[:len(self.objects)]
        
        # Update children
        self.children = self.list_items
        
        logging.info(f'Updated list with {len(self.objects)} objects using {len(self.list_items)} containers')
    
    def _ensure_container_pool_size(self, needed_size):
        """Ensure the container pool has at least needed_size containers."""
        while len(self.container_pool) < needed_size:
            container = self._create_empty_container()
            self.container_pool.append(container)
    
    def _create_empty_container(self):
        """Create a new empty container widget."""
        # Create label widget
        label_widget = EditableLabel(value="")
        
        # Create container for label + buttons
        item_container = widgets.HBox([label_widget])
        item_container.add_class('row-box')
        item_container.label = label_widget

        # Add buttons based on configuration
        hide_button = widgets.Button(
            description='',
            icon='eye',
            tooltip='Hide/Display item'
        )
        hide_button.add_class('hide-btn')
        item_container.children += (hide_button,)
        item_container.hide = hide_button

        delete_button = widgets.Button(
            description='×',
            tooltip='Delete item'
        )
        delete_button.add_class('del-btn')
        item_container.children += (delete_button,)
        item_container.delete = delete_button
        
        return item_container
    
    def _assign_object_to_container(self, obj, container):
        """Assign an object to a specific container and set up bindings."""
        # Update label
        container.label.value = getattr(obj, 'label', str(obj))
        
        # Store mappings
        container.object = obj
        self.label_to_object[container.label] = obj
        self.object_to_container[obj] = container
        
        # Setup bindings
        self._setup_binding(container.label, obj)
        has_hide = hasattr(obj, 'hide') and callable(obj.hide)
        has_delete = hasattr(obj, 'delete') and callable(obj.delete)
        container.hide.layout.display = 'block' if has_hide else 'none'
        container.delete.layout.display = 'block' if has_delete else 'none'
        
        # Setup button handlers
        if has_hide:
            self._setup_hide_button(container.hide, obj)
        if has_delete:
            self._setup_delete_button(container.delete, obj)

        # Setup selection handlers
        def select_container(btn):
            self.selected = container

        def delete_selection(btn):
            if self.selected is container:
                self.selected = None

        container.label.label.on_click(select_container)
        if hasattr(container, 'delete'):
            container.delete.on_click(delete_selection)

    def _clear_all_bindings(self):
        """Clear all current object-to-container bindings and observers."""
        # Ensure _observers exists
        if not hasattr(self, '_observers'):
            self._observers = {}

        # Clean up all observers
        for obj_id in list(self._observers.keys()):
            self._cleanup_object_observers_by_id(obj_id)
        
        # Clear mappings
        self.label_to_object.clear()
        self.object_to_container.clear()
        
        # Clear object references from containers
        for container in self.container_pool:
            if hasattr(container, 'object'):
                delattr(container, 'object')

    @observe('selected')
    def _update_selected(self, change):
        selected = change['new']
        for cnt in self.list_items:
            if cnt is selected:
                cnt.add_class('selected')
            else:
                cnt.remove_class('selected')
    
    def _register_observer(self, target, handler, name):
        """Register an observer for later cleanup."""
        obj_id = id(target)
        if obj_id not in self._observers:
            self._observers[obj_id] = []
        self._observers[obj_id].append({
            'target': target,
            'handler': handler,
            'name': name
        })

    def _cleanup_object_observers_by_id(self, obj_id):
        """Clean up observers for a specific object ID."""
        if obj_id in self._observers:
            for observer_info in self._observers[obj_id]:
                try:
                    if hasattr(observer_info['target'], 'unobserve'):
                        observer_info['target'].unobserve(
                            observer_info['handler'], 
                            observer_info['name']
                        )
                except Exception:
                    pass  # Ignore cleanup errors
            del self._observers[obj_id]

    def _cleanup_object_observers(self, obj):
        """Clean up observers for a specific object."""
        self._cleanup_object_observers_by_id(id(obj))
    
    def _setup_binding(self, label_widget, obj):
        """Setup two-way binding between EditableLabel and object."""
        
        @log_callback_errors
        def on_label_change(change):
            # Update object when label changes
            obj.label = change['new']
        
        @log_callback_errors
        def on_object_change(change):
            # Update label when object changes
            label_widget.value = change['new']
        
        # Bind label changes to object
        label_widget.observe(on_label_change, 'description')
        self._register_observer(label_widget, on_label_change, 'description')
        
        # Bind object changes to label (if object supports observation)
        if hasattr(obj, 'observe'):
            obj.observe(on_object_change, 'label')
            self._register_observer(obj, on_object_change, 'label')
    
    def _setup_delete_button(self, delete_button, obj):
        """Setup delete button click handler."""
        @log_callback_errors
        def on_delete_click(button):
            logging.info(f"Deleting object: {getattr(obj, 'label', str(obj))}")
            if hasattr(obj, 'delete') and callable(getattr(obj, 'delete')):
                try:
                    obj.delete()
                except Exception as e:
                    print(f"Error calling delete method on object: {e}")
            
            # Remove object from list
            self.remove_object(obj)
        
        # Clear any existing handlers to avoid duplicates
        delete_button._click_handlers.callbacks.clear()
        delete_button.on_click(on_delete_click)

    def _setup_hide_button(self, hide_button, obj):
        """Setup hide button click handler."""
        @log_callback_errors
        def on_hide_click(button):
            if button.icon == 'eye':
                if hasattr(obj, 'hide') and callable(getattr(obj, 'hide')):
                    obj.hide()
            else:
                if hasattr(obj, 'show') and callable(getattr(obj, 'show')):
                    obj.show()
            button.icon = 'eye-slash' if button.icon == 'eye' else 'eye'
        
        # Clear any existing handlers to avoid duplicates
        hide_button._click_handlers.callbacks.clear()
        hide_button.on_click(on_hide_click)
    
    def add_object(self, obj, label=None):
        """Add a new object to the list."""
        if label: 
            obj.__dict__['label'] = label
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
    
    def get_container_by_object(self, obj):
        """Get the container widget for a specific object."""
        return self.object_to_container.get(obj)
    
    def replace_object(self, old_obj, new_obj):
        """Replace an object in the list, triggering proper updates."""
        try:
            idx = self.objects.index(old_obj)
            new_objects = list(self.objects)
            new_objects[idx] = new_obj
            self.objects = new_objects  # Trigger trait change
        except ValueError:
            raise ValueError(f"Object {old_obj} not found in list")
    
    def set_objects(self, objects_list):
        """Set the list of objects."""
        self.objects = objects_list

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