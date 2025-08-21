from traitlets import Unicode, observe, List, HasTraits, Instance, Integer
import ipywidgets as widgets
import logging
from .auxfuncs import *
from .widgets import *

class BrowserItem(widgets.VBox):
    browser = None

    def __init__(self, obj, parent, label=None, deletable=True):
        super().__init__()
        self.obj = obj  # The wrapped object
        self.parent = parent  # Reference to the parent browser widget
        self._label = label
        self.deletable = deletable
        self._init_widgets()
        self._setup_bindings()
        self.add_class('browser-item')


    def _init_widgets(self):
        editable_label = EditableLabel(value='')
        icon = SpinnerWidget(icon='spinner', style={'opacity': '0.3'})

        caret_button = widgets.Button(
            description='',
            icon='caret-down',
            tooltip='Expand/Collapse item'
        )
        caret_button.add_class('caret-btn')

        # Add buttons based on configuration
        hide_button = widgets.Button(
            description='',
            icon='eye',
            tooltip='Hide/Display item'
        )
        hide_button.add_class('hide-btn')

        delete_button = widgets.Button(
            description='×',
            tooltip='Delete item'
        )
        delete_button.add_class('del-btn')

        progress = widgets.IntProgress(description='', bar_style='info')
        progress.layout.display = 'none'

        self.editable_label = editable_label
        self.icon = icon
        self.caret = caret_button
        self.hide = hide_button
        self.delete = delete_button
        self.progress = progress
        self.row = widgets.HBox([caret_button, icon, editable_label, hide_button, delete_button])
        self.row.add_class('item-row')
        self.row.layout.padding = f'0 5px 0 {self.parent.level*10}px'  # Indent based on level
        self.children = [self.row, self.progress]

    def _setup_bindings(self):
        self._setup_label()
        self._setup_icon()
        self._setup_delete_button()
        self._setup_hide_button()
        self._setup_caret_button()

    def _setup_label(self):
        label = self._label or getattr(self.obj, 'label', str(self.obj))
        self.label = label

    def _setup_icon(self):
        self.icon.icon = getattr(self.obj, 'icon', 'spinner')
        self.icon.color = getattr(self.obj, 'icon_color', 'auto')
    
    def _setup_delete_button(self):
        """Setup delete button click handler."""
        obj = self.obj
        delete_button = self.delete
        delete_button._click_handlers.callbacks.clear()
        has_delete = hasattr(obj, 'delete') and callable(obj.delete)
        self.delete.layout.display = 'block' if self.deletable else 'none'
        if has_delete and self.deletable:
            @log_callback_errors
            def on_delete_click(button):
                obj.delete()
            
            delete_button.on_click(on_delete_click)

        @log_callback_errors
        def delete_item(btn):
            if self.parent.master.selected is self:
                self.parent.master.selected = None
            # Remove object from list
            self.parent.remove_object(obj)

        delete_button.on_click(delete_item)
        
    def _setup_hide_button(self):
        """Setup hide button click handler."""
        obj = self.obj
        hide_button = self.hide
        hide_button._click_handlers.callbacks.clear()
        has_hide = hasattr(obj, 'hide') and callable(obj.hide)
        has_show = hasattr(obj, 'show') and callable(obj.show)
        if has_hide and has_show:
            hide_button.layout.display = 'block'
            @log_callback_errors
            def on_hide_click(button):
                was_hidden = button.icon == 'eye-slash'
                if was_hidden:
                    obj.show()
                else:
                    obj.hide()
                button.icon = 'eye' if was_hidden else 'eye-slash'
            
            hide_button.on_click(on_hide_click)
        else:
            hide_button.layout.display = 'none'

    def _setup_caret_button(self):
        """Setup caret button click handler for expandable items."""
        caret_button = self.caret
        caret_button._click_handlers.callbacks.clear()
        has_browser = isinstance(self.browser, BrowserWidget)
        caret_button.disabled = not has_browser
        if has_browser:
            caret_button.icon = 'caret-down'
            caret_button.tooltip = 'Expand/Collapse children'
            @log_callback_errors
            def on_caret_click(button):
                was_collapsed = button.icon == 'caret-right'
                if was_collapsed:
                    button.icon = 'caret-down'
                    self.browser.layout.display = 'block'
                else:
                    button.icon = 'caret-right'
                    self.browser.layout.display = 'none'
                
            caret_button.on_click(on_caret_click)
        else:
            caret_button.icon = 'none'
            caret_button.tooltip = ''

    def replace_object(self, new_obj):
        self.obj = new_obj
        self._setup_bindings()
        self.parent.notify_change({'name': 'items',
                                   'old': self.parent.items,
                                   'new': self.parent.items,
                                   'type': 'change'})

    def _init_child_browser(self):
        master = self.parent.master
        level = self.parent.level + 1
        self.browser = BrowserWidget(master=master, level=level)
        self.children = [self.row, self.progress, self.browser]
        self._setup_caret_button()

        @log_callback_errors
        def on_child_browser_change(change):
            logging.info(f"Child browser changed: {change}")
            items = change['new']
            if len(items) == 0:
                self.children = [self.row]
                self.browser = None
                self._setup_caret_button()

        self.browser.observe(on_child_browser_change, names='items')
    
    def add_child(self, obj, **kwargs):
        if self.browser is None:
            self._init_child_browser()

        item = self.browser.add_object(obj, **kwargs)
        return item

    def remove_child(self, obj):
        if self.browser is None:
            logging.warning("No editable list found to remove object from.")
        else:
            self.browser.remove_object(obj)
    
    @property
    def label(self):
        """Access the original wrapped object"""
        return self.editable_label.value
    
    @label.setter
    def label(self, value):
        """Set the label of the wrapped object"""
        self.editable_label.value = value
        # if hasattr(self.obj, 'label'):
        #     self.obj.label = value


class BrowserWidget(widgets.VBox, HasTraits):
    """A list widget that displays objects as EditableLabel widgets."""
    
    items = List()
    selected = Instance(widgets.Box, allow_none=True)
    
    def __init__(self, master=None, level=0, **kwargs):
        super().__init__(**kwargs)
        self.add_class('browser-widget')
        self.master = master or self
        self.level = level  # Level in the hierarchy, used for indentation
        self.selected = None
        self.items = []

        self.master.observe(self._update_selected, names='selected')
    
    @observe('items')
    def _on_items_change(self, change):
        """Update the list when objects change."""
        # Update children
        self.object_to_item = {id(item.obj): item for item in self.items}
        self.children = self.items

    def _create_new_item(self, obj, **kwargs):
        """Create a new BrowserItem for the given object."""
        item = BrowserItem(obj, self, **kwargs)

        # Setup selection handlers
        def select_container(btn):
            self.master.selected = item

        item.editable_label.label.on_click(select_container)
        return item

    def _update_selected(self, change):
        selected = change['new']
        for item in self.items:
            if item is selected:
                item.row.add_class('selected')
            else:
                item.row.remove_class('selected')
    
    def add_object(self, obj, **kwargs):
        """Add a new object to the list."""
        # Wrap in BrowserItem for consistency
        item = self._create_new_item(obj, **kwargs)
        self.items = self.items + [item]  # Trigger trait change
        return item

    def remove_object(self, obj):
        """Remove an object from the list."""
        if obj in self.objects:
            new_items = [i for i in self.items if i.obj != obj]
            self.items = new_items
        else:
            logging.warning(f"Object {obj} not found in browser items.")

    @property
    def objects(self):
        return [i.obj for i in self.items]
    
    def clear(self):
        """Clear all objects from the list."""
        self.items = []
    
    def get_item_by_object(self, obj):
        """Get the item widget for a specific object."""
        return self.object_to_item.get(id(obj))
    
    def get_unique_label(self, label):
        labels = [i.label for i in self.items]
        new_label = get_unique_name(label, labels, '_{:03d}')
        return new_label

    def get_objects_by_type(self, typ, exclude=None):
        objects = [o for o in self.objects if isinstance(o, typ)]
        if exclude:
            objects = [o for o in objects if not isinstance(o, exclude)]
        return objects