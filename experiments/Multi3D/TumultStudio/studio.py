from .auxfuncs import *
from .ServerManager import *
from .LineMask import *
from .widgets import *
from .synth_form import synth_form, launch_synth
from .fit_form import fit_form, fit_synth
from .settings import settings_form, settings_page, saved_settings
from bokeh.models import *
from bokeh.layouts import *
from bokeh.events import *
from bokeh.plotting import figure
import numpy as np
import m3dis
from ipyformkit.core import load_stylesheets
from traitlets import HasTraits

class TumultStudio(HasTraits):
    def __init__(self, source):
        self.source = source
        self.xdata = source.data['x'] # original data
        self.ydata = source.data['y'] # original data
        self.nav_source = ColumnDataSource(dict(x=self.xdata, y=self.ydata))
        self.nlam = len(self.xdata)
        self.max_points = 1000
        self.browser = EditableList()
        self.browser.add_object(self.source, label='Observations')

        title = widgets.Label('Controls:')
        self.controls = widgets.HBox([title])
        self.controls.add_class('controls')
        
        self.editor = widgets.VBox([self.controls])
        self.editor.add_class('editor')
        
        self.right_sidebar = widgets.VBox([])
        self.right_sidebar.add_class('sidebar')
        self.right_sidebar.add_class('right-sidebar')
        
        title = widgets.Label('Objects:')
        self.left_sidebar = widgets.VBox([title])
        self.left_sidebar.add_class('sidebar')
        self.left_sidebar.add_class('left-sidebar')
        
        self.content = widgets.HBox([self.left_sidebar, self.editor, self.right_sidebar])
        self.content.add_class('content')
        
        self.header_title = widgets.Label('TUMULT Studio')
        self.settings_button = widgets.Button(icon='gear')
        self.header = widgets.HBox([self.header_title, self.settings_button])
        self.header.add_class('header')

        self.status = widgets.Label('Status')
        self.status.add_class('status')
        self.footer = widgets.HBox([self.status])
        self.footer.add_class('footer')

        self._gui = [settings_page, self.header, self.content, self.footer]
        self.gui = widgets.VBox(self._gui)
        self.gui.add_class('tum-gui')
        
        self.stylesheets = [StyleSheet(module_dir + 'style.css'),
                            StyleSheet(module_dir + 'widgets.css'),
                            *load_stylesheets()]

    def launch(self, xmin=None, xmax=None, height=500):
        self.init_settings()
        self.init_controls()
        self.init_sidebars()
        canvas = self.init_canvas(xmin=xmin, xmax=xmax)
        self.editor.children += (canvas,)
        self.gui.layout.height = f'{height}px'
        
        return display(*self.stylesheets, self.gui)

    
    def init_settings(self):
        def open_settings(btn):
            settings_page.add_class('in-front')

        def close_settings(btn):
            settings_page.remove_class('in-front')

        self.settings_button.on_click(open_settings)
        settings_page.close_button.on_click(close_settings)

    
    def init_sidebars(self):
        self.init_synth_form()
        self.init_fit_form()
        self.left_sidebar.children += (self.browser,)
        self.right_sidebar.children = (synth_form.vbox,)

        def set_right_sidebar(change):
            obj = self.browser.selected.object
            typ = type(obj)
            if typ==m3dis.m3dis:
                lm = obj.line_mask
                xmin = lm.source.data['x'][lm.islc[0]]
                xmax = lm.source.data['x'][lm.islc[-1]]
                reset_xaxis(self.fig, [xmin, xmax])
                reset_xaxis(self.nav, [xmin, xmax], margin=2)
                if hasattr(obj, 'fit_params'):
                    fit_form.set_values(obj.fit_params)
                self.right_sidebar.children = (fit_form.vbox,)
            elif typ==LineMask:
                lm = obj
                xmin = lm.source.data['x'][lm.islc[0]]
                xmax = lm.source.data['x'][lm.islc[-1]]
                reset_xaxis(self.fig, [xmin, xmax])
                reset_xaxis(self.nav, [xmin, xmax], margin=2)
                synth_form.widgets_dict['line_mask'].wid.value = lm
                self.right_sidebar.children = (synth_form.vbox,)
            elif typ==ColumnDataSource:
                reset_xaxis(self.fig, [self.saved_xmin, self.saved_xmax])
                reset_xaxis(self.nav, [self.xdata[0], self.xdata[-1]])
                self.right_sidebar.children = (synth_form.vbox,)
                
        self.browser.observe(set_right_sidebar, names='selected')
        
                
    def init_controls(self):
        #self.init_norm_buttons()
        self.init_line_mask_button()
        
        
    def init_canvas(self, xmin=None, xmax=None):
        x0 = self.xdata[0]
        x1 = self.xdata[-1]
        
        # Create plot with box select tool AND tap tool
        self.fig = init_main_canvas(sizing_mode='stretch_both',
                               tools="xbox_select,tap,xbox_zoom")
        
        self.nav = init_nav_canvas(sizing_mode='stretch_width', height=100,
                              tools="xbox_select,xwheel_zoom,reset",
                              toolbar_sticky=False)
        
        self.nav.line('x', 'y', source=self.nav_source, color='dimgrey')
        self.nav.x_range = Range1d(x0, x1, bounds=(x0, x1))
        
        self.fig.line('x', 'y', source=self.source, color='lightgrey')
        self.fig.scatter('x', 'y', size=3, source=self.source, color='dimgrey', 
                    selection_color='black', nonselection_color='lightgrey', nonselection_alpha=1)

        def yrange_update_callback(attr, new, old):
            slc = get_crop_range(self.xdata, self.fig.x_range)
            yvis = self.ydata[slc]
            if len(yvis)>0: reset_yaxis(self.fig, yvis, on_server=True)

        self.fig.x_range.on_change('start', yrange_update_callback)
        self.fig.x_range.on_change('end', yrange_update_callback)
        
        if xmin: self.fig.x_range.start = xmin
        if xmax: self.fig.x_range.end   = xmax
        self.saved_xmin = self.fig.x_range.start
        self.saved_xmax = self.fig.x_range.end
            
        @log_callback_errors
        def nav_drag_callback(event):
            geo = event.geometry
            xmin = geo['x0']
            xmax = geo['x1']
                
            self.saved_xmin = xmin
            self.saved_xmax = xmax
            reset_xaxis(self.fig, [xmin, xmax], on_server=True)
            reset_xaxis(self.nav, [xmin, xmax], on_server=True, margin=2)
            
        self.nav.on_event(SelectionGeometry, nav_drag_callback)
        
        # Canvas layout
        bk_canvas = column(
            self.nav,
            self.fig,
            sizing_mode='stretch_both',
        )
        
        # Create main layout with plot on left, sidebar on right
        canvas = iframe_server(bk_canvas)
        canvas.add_class('canvas')
        return canvas

        
    def init_norm_buttons(self):
        norm_button = widgets.Button(description='Normalise')
        norm_reset_button = widgets.Button(description='Reset Normalisation')
        norm_button.layout.display = 'none'
        norm_reset_button.layout.display = 'none'

        @log_callback_errors
        def selection_callback(attr, old, new):
            nothing_selected = len(self.source.selected.indices)==0
            display = 'none' if nothing_selected else 'block'
            norm_button.layout.display = display
        
        @log_callback_errors
        def norm_clicked(btn):
            islc = np.array(self.source.selected.indices)
            selection = self.source.data['x'][islc]
            xmin = np.min(selection)
            xmax = np.max(selection)
            newy = renorm(self.xdata, self.ydata, xmin, xmax)
            update_cds(self.source, 'y', newy)
            norm_reset_button.layout.display = 'block'
            reset_yaxis(self.fig, self.source.data['y'][islc])
            
        @log_callback_errors
        def norm_reset_clicked(btn):
            update_cds(self.source, 'y', self.ydata)
            reset_yaxis(self.fig, self.ydata)
            btn.layout.display = 'none'

        self.source.selected.on_change('indices', selection_callback)
        norm_button.on_click(norm_clicked)
        norm_reset_button.on_click(norm_reset_clicked)
        self.controls.children += (norm_button, norm_reset_button)

    def init_line_mask_button(self):
        new_line_mask_button = widgets.Button(description='New Line Mask')
        new_line_mask_button.layout.display = 'none'
        
        @log_callback_errors
        def selection_callback(attr, old, new):
            nothing_selected = len(self.source.selected.indices)==0
            display = 'none' if nothing_selected else 'block'
            new_line_mask_button.layout.display = display
            
        @log_callback_errors
        def new_line_mask(btn):
            line_masks = self._get_objects_by_type(LineMask)
            i = find_smallest_missing([lm.ic for lm in line_masks])
            label = f'LineMask_{i:03d}'
            new_mask = LineMask(self.source, label=label)
            new_mask.ic = i
            self.browser.add_object(new_mask)
            palette = ColorCycle(saved_settings['bokeh palette'])
            new_mask.draw_boxes(self.nav, self.fig, line_color=None, fill_color=palette[i], fill_alpha=0.1)
        
        self.source.selected.on_change('indices', selection_callback)
        new_line_mask_button.on_click(new_line_mask)
        self.controls.children += (new_line_mask_button,)

    
    def init_synth_form(self):
        launch_button = synth_form.widgets_dict['launch'].wid
        launch_button.on_click(self.plot_synth)
    
        def update_line_mask_dropdown(change):
            line_masks = self._get_objects_by_type(LineMask)
            wid = synth_form.widgets_dict['line_mask'].wid
            wid.options = tuple(line_masks)
    
        self.browser.observe(update_line_mask_dropdown, names='objects')

    
    def init_fit_form(self):
        fit_button = fit_form.widgets_dict['fit'].wid
        fit_button.on_click(self.fit_synth)

    
    @log_callback_errors
    def plot_synth(self, btn):
        new_run = launch_synth()
        synth_runs = self._get_objects_by_type(m3dis.m3dis)
        i = find_smallest_missing([sr.ic for sr in synth_runs])
        new_run.label = f'SynthRun_{i:03d}'
        new_run.ic = i
        
        self.browser.add_object(new_run)
        x, y = new_run.crop(LTE=True)
        new_run.synth = ColumnDataSource(dict(x=x, y=y))
    
        def draw_synth():
            palette = ColorCycle(saved_settings['bokeh palette'])
            new_run.glyph = self.fig.line('x', 'y', source=new_run.synth, color=palette[i])
    
        def delete_synth():
            def update():
                new_run.glyph.visible = False
            server_manager.next_tick_update(update)
    
        new_run.delete = delete_synth
        server_manager.next_tick_update(draw_synth)

    
    @log_callback_errors
    def fit_synth(self, btn):
        run = self.browser.selected.object
        xfit, yfit = fit_synth(run)

        def update():
            run.synth.data['x'] = xfit
            run.synth.data['y'] = yfit

        server_manager.next_tick_update(update)
        

    def _get_objects_by_type(self, typ):
        objects = [o for o in self.browser.objects if type(o)==typ]
        return objects
        
        

def init_main_canvas(xlabel='Wavelength [Å]', ylabel='Flux', **kwargs):
    p = figure(output_backend='svg',
               x_range=[0, 1],
               y_range=[-0.05, 1.05],
               **kwargs)
    p.add_layout(LinearAxis(), 'above')
    p.add_layout(LinearAxis(), 'right')
    # Optional: show ticks but not labels on the top/right
    p.axis[0].major_label_text_font_size = '0pt'  # Top
    p.axis[3].major_label_text_font_size = '0pt'  # Right
    p.xaxis[1].axis_label = xlabel
    p.yaxis[0].axis_label = ylabel
    p.toolbar.logo = None
    tools = {type(tool):tool for tool in p.toolbar.tools}
    p.toolbar.active_drag = tools[BoxSelectTool]
    return p

def init_nav_canvas(**kwargs):
    p = figure(toolbar_location='above',
               x_axis_location='above',
               **kwargs)
    p.yaxis.visible = False
    p.xaxis.axis_line_color = None
    p.xaxis.major_tick_line_color = None
    p.xaxis.minor_tick_line_color = None
    p.toolbar.logo = None
    tools = {type(tool):tool for tool in p.toolbar.tools}
    p.toolbar.active_scroll = tools[WheelZoomTool]
    return p

def reset_xaxis(p, xx, margin=0, on_server=False):
    xmin = np.min(xx)
    xmax = np.max(xx)
    d = (xmax-xmin)*margin

    def update():
        p.x_range.start = xmin-d
        p.x_range.end   = xmax+d

    if on_server:
        update()
    else:
        server_manager.next_tick_update(update)
        
def reset_yaxis(p, yy, margin=0.05, on_server=False):
    ymin = np.min(yy)
    ymax = np.max(yy)
    d = (ymax-ymin)*margin

    def update():
        p.y_range.start = ymin-d
        p.y_range.end   = ymax+d

    if on_server:
        update()
    else:
        server_manager.next_tick_update(update)