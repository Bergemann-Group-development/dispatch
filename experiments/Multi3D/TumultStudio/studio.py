from .auxfuncs import *
from .ServerManager import *
from .LineMask import *
from .widgets import *
from .synth_form import synth_form, launch_synth, monitor_synth
from .fit_form import fit_form, fit_synth
from .settings import settings_form, settings_page, saved_settings
from .llist_annotations import plot_line_annotations, remove_line_annotations, get_llist_table
from bokeh.models import *
from bokeh.layouts import *
from bokeh.events import *
from bokeh.plotting import figure
import numpy as np
import m3dis
import time
from traitlets import HasTraits

class TumultStudio(HasTraits):
    def __init__(self, source):
        self.source = source
        self.xdata = source.data['x'] # original data
        self.ydata = source.data['y'] # original data
        self.nav_source = ColumnDataSource(dict(x=[], y=[]))
        self.fig_source = ColumnDataSource(dict(x=[], y=[]))
        self.synth_source = ColumnDataSource(dict(xs=[], ys=[], colors=[], labels=[], run_ids=[]))
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
        self.gui = widgets.VBox(self._gui)
        self.gui.add_class('tum-gui')
        
        self.stylesheets = [*synth_form.load_stylesheets(),
                            StyleSheet(module_dir + 'style.css'),
                            StyleSheet(module_dir + 'widgets.css'),
                            synth_form.load_tooltip_script()]

    def launch(self, xmin=None, xmax=None, height=500):
        self.init_settings()
        self.init_controls()
        self.init_sidebars()
        canvas = self.init_canvas(xmin=xmin, xmax=xmax)
        self.editor.children += (canvas,)
        self.gui.layout.height = f'{height}px'
        
        return display(self.gui, *self.stylesheets)

    
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
        self.right_sidebar.titles = ['Synthesis Settings']
        self.rerun_button = widgets.Button(description='Rerun')

        def set_right_sidebar(change):
            sel = self.browser.selected
            if not sel:
                reset_xaxis(self.fig, [self.saved_xmin, self.saved_xmax])
                reset_xaxis(self.nav, [self.xdata[0], self.xdata[-1]])
                self.right_sidebar.children = (synth_form.vbox,)
                self.right_sidebar.titles = ['Synthesis Settings']
            else:
                obj = sel.object
                typ = type(obj)
                if typ==m3dis.m3dis:
                    lm = obj.line_mask
                    reset_xaxis(self.fig, [lm.xmin, lm.xmax])
                    reset_xaxis(self.nav, [lm.xmin, lm.xmax], margin=2)
                    if hasattr(obj, 'fit_params'):
                        fit_form.set_values(obj.fit_params)
                    llist_vbox = widgets.VBox([obj.llist_html, self.rerun_button])
                    llist_vbox.add_class('llist-vbox')
                    self.right_sidebar.children = (fit_form.vbox, llist_vbox)
                    self.right_sidebar.titles = ['Fit Settings', 'Line List']
                elif typ==LineMask:
                    lm = obj
                    reset_xaxis(self.fig, [lm.xmin, lm.xmax], margin=0.05)
                    reset_xaxis(self.nav, [lm.xmin, lm.xmax], margin=2)
                    synth_form.widgets_dict['line_mask'].wid.value = lm
                    self.right_sidebar.children = (synth_form.vbox,)
                    self.right_sidebar.titles = ['Synthesis Settings']
                elif typ==ColumnDataSource:
                    reset_xaxis(self.fig, [self.saved_xmin, self.saved_xmax])
                    reset_xaxis(self.nav, [self.xdata[0], self.xdata[-1]])
                    self.right_sidebar.children = (synth_form.vbox,)
                    self.right_sidebar.titles = ['Synthesis Settings']
                
        self.browser.observe(set_right_sidebar, names='selected')
        
                
    def init_controls(self):
        self.init_line_mask_button()
        
        
    def init_canvas(self, xmin=None, xmax=None):
        # Create plot with box select tool AND tap tool
        self.fig = self.init_main_canvas(xmin, xmax, sizing_mode='stretch_both',
                               tools="xbox_select,tap,xbox_zoom,pan")
        
        self.nav = self.init_nav_canvas(sizing_mode='stretch_width', height=100,
                              tools="xbox_select,xwheel_zoom,reset",
                              toolbar_sticky=False)
        
        self.nav.line('x', 'y', source=self.nav_source, color='dimgrey')

        self.fig.line('x', 'y', source=self.fig_source, color='lightgrey')
        self.fig.scatter('x', 'y', size=4, source=self.fig_source, color='dimgrey', 
                    selection_color='black', nonselection_color='lightgrey', nonselection_alpha=1)
        
        self.synth_line = self.fig.multi_line('xs', 'ys', line_color='colors', legend_field='labels',
                                              source=self.synth_source, line_width=3,
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
            reset_xaxis(self.fig, [xmin, xmax], on_server=True)
            reset_xaxis(self.nav, [xmin, xmax], on_server=True, margin=2)
            
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
            frame = run.llist_html.html_widget
            new = frame.value
            new = new.replace('<tr class="highlighted-row">', '<tr>')
            new = replace_nth(new, '<tr>', '<tr class="highlighted-row">', i+1)
            frame.value = new

        self.llist_source.selected.on_change('indices', synth_line_selection)
        
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
        launch_button = synth_form.widgets_dict['launch'].wid
        launch_button.on_click(self.plot_synth)
    
        def update_line_mask_dropdown(change):
            line_masks = self._get_objects_by_type(LineMask)
            wid = synth_form.widgets_dict['line_mask'].wid
            if line_masks:
                wid.options = tuple(line_masks)
                wid.value = line_masks[-1]
            else:
                wid.options = ('no line masks yet',)
    
        self.browser.observe(update_line_mask_dropdown, names='objects')

    
    def init_fit_form(self):
        fit_button = fit_form.widgets_dict['fit'].wid
        fit_button.on_click(self.fit_synth)

    
    @log_callback_errors
    def plot_synth(self, btn):
        synth_runs = self._get_objects_by_type(m3dis.m3dis)
        i = find_smallest_missing([sr.ic for sr in synth_runs])
        color = ColorCycle(saved_settings['bokeh palette'])[i]
        name = f'SynthRun_{i:03d}'
        self.progress.description = f'Running {name}:'

        thread = run_in_thread(launch_synth)(name)

        while thread.is_alive():
            monitor_synth(name, self.progress)
            time.sleep(0.25)

        new_run = thread.out['result']
        new_run.label = name
        new_run.ic = i
        
        self.progress.description = f'Plotting {name}:'
        self.browser.add_object(new_run)
        widget_label = self.browser.list_items[-1].children[0].label
        x, y = new_run.crop(LTE=True)
        rid = id(new_run)
    
        @log_callback_errors
        def draw_synth():
            current_data = self.synth_source.data.copy()
            current_data['xs'].append(x)
            current_data['ys'].append(y)
            current_data['colors'].append(color)
            current_data['labels'].append(new_run.label)
            current_data['run_ids'].append(rid)
            self.synth_source.data = current_data
            self.progress.description = f'Done!'
            self.progress.value = self.progress.max

        @log_callback_errors
        def delete_synth():
            def update():
                i = self.synth_source.data['run_ids'].index(rid)
                current_data = self.synth_source.data.copy()
                for key in current_data:
                    current_data[key].pop(i)
                self.synth_source.data = current_data
            server_manager.next_tick_update(update)

        @log_callback_errors
        def update_label(change):
            def update():
                i = self.synth_source.data['run_ids'].index(rid)
                current_data = self.synth_source.data.copy()
                current_data['labels'][i] = change['new']
                self.synth_source.data = current_data
            server_manager.next_tick_update(update)
    
        new_run.delete = delete_synth
        server_manager.next_tick_update(draw_synth)
        widget_label.observe(update_label, names='description')

    
    @log_callback_errors
    def fit_synth(self, btn):
        run = self.browser.selected.object
        i = self.synth_source.data['run_ids'].index(id(run))
        xfit, yfit = fit_synth(run, self.fig_source)

        def update():
            current_data = self.synth_source.data.copy()
            current_data['xs'][i] = xfit
            current_data['ys'][i] = yfit
            self.synth_source.data = current_data

        server_manager.next_tick_update(update)
        

    def _get_objects_by_type(self, typ):
        objects = [o for o in self.browser.objects if type(o)==typ]
        return objects
        
        
    def init_main_canvas(self, xmin, xmax, xlabel='Wavelength [Å]', ylabel='Flux', **kwargs):
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

        def xrange_update_callback(attr, new, old):
            slc = get_crop_range(self.xdata, p.x_range, nmax=1000)
            self.fig_source.data = dict(x=self.xdata[slc], y=self.ydata[slc])

        def yrange_update_callback(attr, new, old):
            yvis = self.fig_source.data['y']
            for i, x in enumerate(self.synth_source.data['xs']):
                slc = get_crop_range(x, p.x_range)
                yvis = np.concatenate((yvis, self.synth_source.data['ys'][i][slc]))
            if len(yvis)>0: reset_yaxis(p, yvis, on_server=True)
            remove_line_annotations(self.llist_source)
            selected = self.browser.selected
            if selected:
                obj = selected.object
                if type(obj)==m3dis.m3dis:
                    plot_line_annotations(p, obj, self.llist_source)

        p.x_range.on_change('start', xrange_update_callback)
        p.x_range.on_change('end', xrange_update_callback)
        p.x_range.on_change('start', yrange_update_callback)
        p.x_range.on_change('end', yrange_update_callback)
        
        if xmin: p.x_range.start = xmin
        if xmax: p.x_range.end   = xmax
        self.saved_xmin = p.x_range.start
        self.saved_xmax = p.x_range.end
        return p

    def init_nav_canvas(self, **kwargs):
        x0 = self.xdata[0]
        x1 = self.xdata[-1]
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

        def xrange_update_callback(attr, new, old):
            slc = get_crop_range(self.xdata, p.x_range, nmax=10000)
            self.nav_source.data = dict(x=self.xdata[slc], y=self.ydata[slc])

        p.x_range = Range1d(0, 1, bounds=(x0, x1))
        p.x_range.on_change('start', xrange_update_callback)
        p.x_range.on_change('end', xrange_update_callback)
        p.x_range.start = x0
        p.x_range.end = x1
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