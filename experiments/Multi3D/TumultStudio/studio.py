from .auxfuncs import *
from . import ServerManager
from .ServerManager import reset_xaxis, reset_yaxis
from .LineMask import *
from .widgets import *
from .synth_form import synth_form, launch_synth, monitor_synth
from .fit_form import fit_form, set_fit_form, fit_synth
from .rerun_form import rerun_form, rerun_synth
from .settings import settings_form, settings_page, saved_settings
from .llist_annotations import plot_line_annotations, remove_line_annotations, get_llist_table
from bokeh.models import *
from bokeh.layouts import *
from bokeh.events import *
from bokeh.plotting import figure
from functools import wraps
import numpy as np
import tumult
import time
from traitlets import HasTraits

class TumultStudio(HasTraits):
    def __init__(self, source):
        self.source = source
        self.xdata = source.data['x'] # original data
        self.ydata = source.data['y'] # original data
        self.nav_source = ColumnDataSource(dict(x=[], y=[]))
        self.fig_source = ColumnDataSource(dict(x=[], y=[]))
        self.synth_source = ColumnDataSource(dict(xs=[], ys=[], colors=[], labels=[], alphas=[],
                                                  vmic=[], FeH=[], atmos=[], run_ids=[]))
        self.llist_source = ColumnDataSource(dict(xs=[], ys=[], x=[], y=[], text=[]))
        self.nlam = len(self.xdata)
        self.max_points = 1000
        self.browser = EditableList()
        self.browser.add_object(self.source, label='Observations')

        title = widgets.Label('Controls:')
        self.controls = widgets.HBox([title])
        self.controls.add_class('controls')
        
        self.editor = widgets.VBox([self.controls])
        self.editor.add_class('editor')
        
        self.right_sidebar = widgets.Tab([])
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
        self.progress = widgets.IntProgress(description=' ', bar_style='info')
        self.status.add_class('status')
        self.footer = widgets.HBox([self.status, self.progress])
        self.footer.add_class('footer')

        self._gui = [settings_page, self.header, self.content, self.footer]
        
        self.stylesheets = [*synth_form.load_stylesheets(),
                            StyleSheet(module_dir + 'style.css'),
                            StyleSheet(module_dir + 'widgets.css'),
                            synth_form.load_tooltip_script()]
        
        # Load image from file
        with open(module_dir + 'SplashScreen.png', 'rb') as f:
            image_data = f.read()

        # Create splash screen with image
        splash_image = widgets.Image(value=image_data, format='png')
        splash_image.add_class('splash-image')

        splash_text = widgets.Label('TumultStudio')
        splash_text.add_class('splash-text')

        self.gui = widgets.VBox([])
        #self.gui.children = [splash_image, splash_text]
        self.gui.add_class('tum-gui')

    def launch(self, xmin=None, xmax=None, height=500):
        ServerManager.create_global()
        self.gui.layout.height = f'{height}px'
        display(self.gui, *self.stylesheets)
        self.init_settings()
        self.init_controls()
        self.init_sidebars()
        canvas = self.init_canvas(xmin=xmin, xmax=xmax)
        self.editor.children += (canvas,)
        self.gui.children = self._gui
    
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
        self.init_rerun_form()
        self.left_sidebar.children += (self.browser,)
        self.right_sidebar.children = (synth_form.vbox,)
        self.right_sidebar.titles = ['Synthesis Settings']

        def set_right_sidebar(change):
            sel = self.browser.selected
            if not sel:
                reset_xaxis(self.fig, [self.saved_xmin, self.saved_xmax])
                reset_xaxis(self.nav, [self.xdata[0], self.xdata[-1]])
                self.right_sidebar.children = (synth_form.vbox,)
                self.right_sidebar.titles = ['Synthesis Settings']
            else:
                obj = sel.object
                if isinstance(obj, tumult.tumult):
                    run = obj
                    lm = run.line_mask
                    ids = self.synth_source.data['run_ids']
                    i = ids.index(id(run))
                    reorder_data(self.synth_source, i, len(ids))  # move to back to display on top
                    reset_xaxis(self.fig, [lm.xmin, lm.xmax])
                    reset_xaxis(self.nav, [lm.xmin, lm.xmax], margin=2)
                    set_fit_form(run)
                    llist_vbox = widgets.VBox([run.llist_html, rerun_form.vbox])
                    llist_vbox.add_class('llist-vbox')
                    self.right_sidebar.children = (fit_form.vbox, llist_vbox)
                    self.right_sidebar.titles = ['Fit Settings', 'Line List']
                elif isinstance(obj, LineMask):
                    lm = obj
                    reset_xaxis(self.fig, [lm.xmin, lm.xmax])
                    reset_xaxis(self.nav, [lm.xmin, lm.xmax], margin=2)
                    synth_form.widgets_dict['line_mask'].wid.value = lm
                    self.right_sidebar.children = (synth_form.vbox,)
                    self.right_sidebar.titles = ['Synthesis Settings']
                elif isinstance(obj, ColumnDataSource):
                    reset_xaxis(self.fig, [self.saved_xmin, self.saved_xmax])
                    reset_xaxis(self.nav, [self.xdata[0], self.xdata[-1]])
                    self.right_sidebar.children = (synth_form.vbox,)
                    self.right_sidebar.titles = ['Synthesis Settings']
            
            ServerManager.server.next_tick_update(self.yrange_update_callback)
            ServerManager.server.next_tick_update(self.update_annotations)
                
        self.browser.observe(set_right_sidebar, names='selected')
        
                
    def init_controls(self):
        self.init_line_mask_button()
        
        
    def init_canvas(self, xmin=None, xmax=None):
        # Create plot with box select tool AND tap tool
        self.init_main_canvas(xmin, xmax, sizing_mode='stretch_both',
                              tools="xbox_select,tap,xbox_zoom,pan,save")
        
        self.init_nav_canvas(sizing_mode='stretch_width', height=100,
                             tools="xbox_select,xwheel_zoom,reset",
                             toolbar_sticky=False)
        
        self.nav.line('x', 'y', source=self.nav_source, color='dimgrey')

        self.obs_line = self.fig.line('x', 'y', source=self.fig_source, color='lightgrey')
        self.obs_scatter = self.fig.scatter('x', 'y', size=4, source=self.fig_source, color='dimgrey', 
                    selection_color='black', nonselection_color='lightgrey', nonselection_alpha=1)
        
        self.synth_line = self.fig.multi_line('xs', 'ys', line_color='colors', legend_field='labels',
                                              source=self.synth_source, line_width=3, line_alpha='alphas',
                                              nonselection_line_alpha=1.0,
                                              hover_line_width=5)
        
        self.llist_line = self.fig.multi_line('xs', 'ys', line_color='black',
                                              source=self.llist_source, line_width=0.8,
                                              nonselection_line_alpha=1.0,
                                              hover_line_width=3, )
        
        self.llist_labels = LabelSet(x='x', y='y', text='text',
                                     text_font_size='10pt',
                                     text_font='times',
                                     text_align='left', text_baseline='middle',
                                     text_color='black',
                                     source=self.llist_source, angle=np.pi/2)
        
        # TOOLTIPS = [
        #     ('λ', '$x{0.000}'),
        #     ('[Fe/H]', '@FeH{0.00}'),
        #     ('vmic', '@vmic{0.00}'),
        # ]
        
        hover = HoverTool(renderers=[self.llist_line, self.synth_line],
                  tooltips=None, mode='mouse')  # You can add tooltips if needed
        
        self.fig.add_layout(self.llist_labels)
        self.fig.add_tools(hover)
        self.fig.legend.location = "bottom_left"
            
        @log_callback_errors
        def nav_drag_callback(event):
            geo = event.geometry
            xmin = geo['x0']
            xmax = geo['x1']
                
            self.saved_xmin = xmin
            self.saved_xmax = xmax
            reset_xaxis(self.fig, [xmin, xmax])
            reset_xaxis(self.nav, [xmin, xmax], margin=2)
            
        self.nav.on_event(SelectionGeometry, nav_drag_callback)

        def synth_selection(attr, old, new):
            selected = new
            if not selected: return
            idx = selected[0] # last clicked line
            rid = self.synth_source.data['run_ids'][idx]
            item = [item for item in self.browser.list_items if id(item.object)==rid]
            if item:
                self.browser.selected = item[0]
            else:
                self.browser.selected = None

        self.synth_source.selected.on_change('indices', synth_selection)

        def synth_line_selection(attr, old, new):
            selected = new
            if not selected: return
            idx = selected[0] # last clicked line
            run = self.browser.selected.object
            llist = run.llist
            lam0 = self.llist_source.data['xs'][idx][0]
            i = llist.reset_index().index[llist.lam0==lam0][0]
            run.llist_html.highlight_row(i)
            self.right_sidebar.selected_index = 1  # switch to Line List tab

        self.llist_source.selected.on_change('indices', synth_line_selection)
        
        # Canvas layout
        bk_canvas = column(
            self.nav,
            self.fig,
            sizing_mode='stretch_both',
        )
        
        # Create main layout with plot on left, sidebar on right
        canvas = ServerManager.iframe_server(bk_canvas)
        canvas.add_class('canvas')
        return canvas

    def init_line_mask_button(self):
        new_line_mask_button = widgets.Button(description='New Line Mask')
        new_line_mask_button.layout.display = 'none'
        
        @log_callback_errors
        def selection_callback(attr, old, new):
            nothing_selected = len(self.fig_source.selected.indices)==0
            display = 'none' if nothing_selected else 'block'
            new_line_mask_button.layout.display = display
            
        @log_callback_errors
        def new_line_mask(btn):
            line_masks = self._get_objects_by_type(LineMask)
            i = find_smallest_missing([lm.ic for lm in line_masks])
            label = f'LineMask_{i:03d}'
            new_mask = LineMask(self.fig_source, self.source.data, label=label)
            new_mask.ic = i
            self.browser.add_object(new_mask)
            color = ColorCycle(saved_settings['bokeh palette'])[i]
            new_mask.draw_boxes(self.nav, self.fig, line_color=None, fill_color=color, fill_alpha=0.1)
        
        self.fig_source.selected.on_change('indices', selection_callback)
        new_line_mask_button.on_click(new_line_mask)
        self.controls.children += (new_line_mask_button,)

    
    def init_synth_form(self):
        def launch_on_click(btn):
            self.plot_synth(launch_synth)

        launch_button = synth_form.widgets_dict['launch'].wid
        launch_button.on_click(launch_on_click)
    
        def update_line_mask_dropdown(change):
            line_masks = self._get_objects_by_type(LineMask)
            wid = synth_form.widgets_dict['line_mask'].wid
            if line_masks:
                wid.options = tuple(line_masks)
                wid.value = line_masks[-1]
            else:
                wid.options = ('no line masks yet',)
    
        self.browser.observe(update_line_mask_dropdown, names='objects')
        update_line_mask_dropdown(None)

    
    def init_fit_form(self):
        fit_button = fit_form.widgets_dict['apply'].wid
        fit_button.on_click(self.fit_synth)


    def init_rerun_form(self):
        def rerun_on_click(btn):
            run = self.browser.selected.object
            self.plot_synth(rerun_synth, run)
        
        rerun_button = rerun_form.widgets_dict['rerun'].wid
        rerun_button.on_click(rerun_on_click)

    
    @log_callback_errors
    def plot_synth(self, launch_function, *args):
        synth_runs = self._get_objects_by_type(tumult.tumult)
        i = find_smallest_missing([sr.ic for sr in synth_runs])
        color = ColorCycle(saved_settings['bokeh palette'])[i]
        name = f'SynthRun_{i:03d}'

        self.progress.description = f'Running {name}:'
        new_run = launch_function(name, *args, progress=self.progress)
        new_run.ic = i
        
        self.progress.description = f'Plotting {name}:'
        self.browser.add_object(new_run)
        widget_label = self.browser.list_items[-1].children[0].label

        if hasattr(new_run, 'fit_params'):
            crop_params = new_run.fit_params
        else:
            crop_params = {}

        x, y = new_run.crop(LTE=True, **crop_params)

        @update_data_source(self.synth_source)
        def draw_synth(data):
            data['xs'].append(x)
            data['ys'].append(y)
            data['colors'].append(color)
            data['labels'].append(new_run.label)
            data['vmic'].append(new_run.nml['atmos_params']['vmic'])
            data['FeH'].append(new_run.nml['atmos_params']['FeH'])
            data['atmos'].append(new_run.nml['atmos_params']['atmos_file'])
            data['run_ids'].append(id(new_run))
            data['alphas'].append(1.0)

        draw_synth()
        self.set_synth_controls(new_run, widget_label)
        self.progress.value = self.progress.max
        self.progress.description = f'Done!'

    
    @log_callback_errors
    def fit_synth(self, btn):
        old_run = self.browser.selected.object
        i = self.synth_source.data['run_ids'].index(id(old_run))
        self.progress.description = f'Fitting {old_run.label}:'
        new_run, xfit, yfit = fit_synth(old_run, self.fig_source, self.progress)

        @update_data_source(self.synth_source)
        def update_synth(data):
            data['xs'][i] = xfit
            data['ys'][i] = yfit
            data['run_ids'][i] = id(new_run)
            data['labels'][i] = new_run.label

            if new_run is not old_run:
                self.browser.remove_object(old_run)
                self.browser.add_object(new_run)
                self.browser.selected = self.browser.list_items[-1]  # select the new run
                widget_label = self.browser.selected.children[0].label
                self.set_synth_controls(new_run, widget_label)

        update_synth()
        self.progress.value = self.progress.max
        self.progress.description = f'Done!'

    def set_synth_controls(self, run, widget_label=None):
        rid = id(run)

        @update_data_source(self.synth_source)
        def delete_synth(data):
            i = data['run_ids'].index(rid)
            for key in data:
                data[key].pop(i)
        
        @update_data_source(self.synth_source)
        def hide_synth(data):
            i = data['run_ids'].index(rid)
            data['alphas'][i] = 0.0

        @update_data_source(self.synth_source)
        def show_synth(data):
            i = data['run_ids'].index(rid)
            data['alphas'][i] = 1.0

        run.delete = delete_synth
        run.hide = hide_synth
        run.show = show_synth

        if widget_label:
            def update_label(change):
                @update_data_source(self.synth_source)
                def update(data):
                    i = data['run_ids'].index(rid)
                    data['labels'][i] = change['new']
                update()

            widget_label.observe(update_label, names='description')


    def _get_objects_by_type(self, typ):
        objects = [o for o in self.browser.objects if isinstance(o, typ)]
        return objects
        
        
    def init_main_canvas(self, xmin, xmax, xlabel='Wavelength [Å]', ylabel='Flux', **kwargs):
        p = figure(output_backend='svg',
                x_range=[0, 1],
                y_range=[-0.05, 1.05],
                **kwargs)
        
        self.fig = p
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

        def xrange_update_callback(attr, new, old):
            slc = get_crop_range(self.xdata, p.x_range, nmax=1000, margin=0.01)
            self.fig_source.data = dict(x=self.xdata[slc], y=self.ydata[slc])
            if hasattr(self, 'fig_slc'):
                off = self.fig_slc.start - slc.start
                idx = np.array(self.fig_source.selected.indices) + off
                self.fig_source.selected.indices = list(idx)
            self.fig_slc = slc

        p.x_range.on_change('start', xrange_update_callback)
        p.x_range.on_change('end', xrange_update_callback)
        p.x_range.on_change('start', self.yrange_update_callback)
        p.x_range.on_change('end', self.yrange_update_callback)
        p.x_range.on_change('start', self.update_annotations)
        p.x_range.on_change('end', self.update_annotations)
        
        if xmin: p.x_range.start = xmin
        if xmax: p.x_range.end   = xmax
        self.saved_xmin = p.x_range.start
        self.saved_xmax = p.x_range.end
    
    def yrange_update_callback(self, attr, old, new):
        yvis = self.fig_source.data['y']
        for i, x in enumerate(self.synth_source.data['xs']):
            slc = get_crop_range(x, self.fig.x_range)
            yvis = np.concatenate((yvis, self.synth_source.data['ys'][i][slc]))
        if len(yvis)>0:
            reset_yaxis(self.fig, yvis)

    def update_annotations(self, attr, old, new):
        remove_line_annotations(self.llist_source)
        selected = self.browser.selected
        if selected:
            obj = selected.object
            if isinstance(obj, tumult.tumult):
                plot_line_annotations(self.fig, obj, self.llist_source)

    def init_nav_canvas(self, **kwargs):
        x0 = self.xdata[0]
        x1 = self.xdata[-1]
        p = figure(output_backend='svg',
                   toolbar_location='above',
                   x_axis_location='above',
                   **kwargs)
        self.nav = p
        p.yaxis.visible = False
        p.xaxis.axis_line_color = None
        p.xaxis.major_tick_line_color = None
        p.xaxis.minor_tick_line_color = None
        p.toolbar.logo = None
        tools = {type(tool):tool for tool in p.toolbar.tools}
        p.toolbar.active_scroll = tools[WheelZoomTool]

        def xrange_update_callback(attr, new, old):
            slc = get_crop_range(self.xdata, p.x_range, nmax=10000)
            self.nav_source.data = dict(x=self.xdata[slc], y=self.ydata[slc])

        p.x_range = Range1d(0, 1, bounds=(x0, x1))
        p.x_range.on_change('start', xrange_update_callback)
        p.x_range.on_change('end', xrange_update_callback)
        p.x_range.start = x0
        p.x_range.end = x1
        return p

def reorder_data(source, iold, inew):
    @update_data_source(source)
    def reorder(data):
        for key in data:
            item = data[key].pop(iold)
            data[key].insert(inew, item)
    reorder()


def update_data_source(source):
    def decorator(update_func):
        @wraps(update_func)
        def wrapper():
            @log_callback_errors
            def update():
                current_data = source.data.copy()
                # Call the specific update function
                update_func(current_data)
                source.data = current_data
                        
            ServerManager.server.next_tick_update(update)
        return wrapper
    return decorator