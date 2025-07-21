import ipywidgets as widgets
from traitlets import Unicode, observe, List, HasTraits, Instance
import time
import os

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