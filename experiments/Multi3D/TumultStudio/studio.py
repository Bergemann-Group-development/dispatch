from .auxfuncs import *
from . import ServerManager
from .ServerManager import *
from .LineMask import *
from .widgets import *
from .browser import BrowserWidget
from .synth import StudioSynth
from .StudioObs import StudioObs
from types import SimpleNamespace
from .synth_form import synth_form, init_synth_form
from .fit_form import fit_form, set_fit_form, init_fit_form
from .rerun_form import rerun_form, init_rerun_form
from .atmos_grid_form import init_interp_atmos_form
from .norm_form import *
from .departure_form import init_departure_form
from .StudioDepartures import StudioDeparture, init_departure_plot_form
from .atmos_form import init_import_atmos_button
from .StudioAtmos import StudioAtmos, init_atmos_plot_form
from .StudioGrid import StudioGrid
from .settings import settings_form, settings_page, saved_settings
from .threading import execute_in_thread, profile_checkbox, get_profile_results
from bokeh.models import *
from bokeh.models import ColumnDataSource as CDS
from bokeh.layouts import *
from bokeh.events import *
from bokeh.plotting import figure
import numpy as np
from traitlets import HasTraits
import pathlib

first_time = True

class TumultStudio(HasTraits):
    nav_source = CDS(dict(x=[], y=[]))
    fig_source = CDS(dict(x=[], y=[]))
    fig_line_source = CDS(dict(x=[], y=[]))
    synth_source = CDS(dict(xs=[], ys=[], colors=[], labels=[], alphas=[],
                                    vmic=[], FeH=[], atmos=[], run_ids=[]))
    llist_source = CDS(dict(xs=[], ys=[], x=[], y=[], text=[]))
    cntm_source = CDS(dict(x=[], y=[]))
    spline_source = CDS(dict(x=[], y=[]))
    atmos_scatter_source = CDS(dict(x=[], y=[]))
    atmos_line_source = CDS(dict(x=[], y=[]))
    atmos_interp_source = CDS(dict(x=[], y=[]))
    
    synth_folder = pathlib.Path(saved_settings['synthesis folder'])
    synth_folder.mkdir(parents=True, exist_ok=True)
    synth_folder = create_temp_dir_for_kernel(synth_folder)

    def __init__(self, **kwargs):
        self.obs = StudioObs(**kwargs)
        self.browser = BrowserWidget()
        self.browser.add_object(self.obs, label='Observations', deletable=False)

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
        self.header = widgets.HBox([self.header_title, profile_checkbox, self.settings_button])
        self.header.add_class('header')

        self.status = widgets.Label('Status')
        self.progress = widgets.IntProgress(description=' ', bar_style='info')
        self.status.add_class('status')
        self.footer = widgets.HBox([self.status, self.progress])
        self.footer.add_class('footer')

        self._gui = [settings_page, self.header, self.content, self.footer]
        
        self.stylesheets = [*synth_form.load_stylesheets(),
                            StyleSheet(module_dir / 'style.css'),
                            StyleSheet(module_dir / 'widgets.css'),
                            synth_form.load_tooltip_script()]
        
        # # Load image from file
        # with open(module_dir / 'SplashScreen.png', 'rb') as f:
        #     image_data = f.read()

        # # Create splash screen with image
        # splash_image = widgets.Image(value=image_data, format='png')
        # splash_image.add_class('splash-image')

        # splash_text = widgets.Label('TumultStudio')
        # splash_text.add_class('splash-text')

        self.gui = widgets.VBox([])
        #self.gui.children = [splash_image, splash_text]
        self.gui.add_class('tum-gui')

    def launch(self, xmin=None, xmax=None, height=500):
        #global first_time
        #if first_time:
        #    first_time = False
        #else:
        #    raise RuntimeError("Studio can only be launched once per kernel session.")
        
        ServerManager.create_global()
        self.gui.layout.height = f'{height}px'
        self.init_settings()
        canvas = self.init_canvas(xmin=xmin, xmax=xmax)
        self.init_controls()
        self.init_sidebars()
        self.init_browser()
        self.editor.children += (canvas,)
        self.gui.children = self._gui
        display(self.gui, *self.stylesheets)
    
    def init_settings(self):
        def open_settings(btn):
            settings_page.add_class('in-front')

        def close_settings(btn):
            settings_page.remove_class('in-front')

        self.settings_button.on_click(open_settings)
        settings_page.close_button.on_click(close_settings)

    def init_sidebars(self):
        init_synth_form(self)
        init_fit_form(self)
        init_rerun_form(self)
        init_norm_form(self)
        init_atmos_plot_form(self)
        init_departure_form(self)
        init_interp_atmos_form(self)
        init_departure_plot_form(self)
        self.left_sidebar.children += (self.browser,)
        self.right_sidebar.children = (synth_form.vbox,)
        self.right_sidebar.titles = ['Synthesis Settings']
                
    def init_controls(self):
        init_norm_checkbox(self)
        init_import_atmos_button(self)
        init_line_mask_button(self)
        
    def init_browser(self):
        def on_browser_selection(change):
            obj = self.current

            if obj is None:
                self.right_sidebar.children = (synth_form.vbox,)
                self.right_sidebar.titles = ['Synthesis Settings']
            elif hasattr(obj, 'chosen'):
                def update_canvas():
                    show_atmos = isinstance(obj, (StudioAtmos, StudioDeparture, StudioGrid))
                    self.atmos_fig.visible = show_atmos
                    self.nav.visible = not show_atmos
                    self.fig.visible = not show_atmos
                
                ServerManager.server.next_tick_update(update_canvas)

                d = obj.chosen(self)
                if d is not None:
                    self.right_sidebar.children = d['tabs']
                    self.right_sidebar.titles = d['titles']

                    if 'xmin' in d and 'xmax' in d:
                        reset_xaxis(self.fig, [d['xmin'], d['xmax']])
                        reset_xaxis(self.nav, [d['xmin'], d['xmax']], margin=2)
            
            self.reset_yaxis_and_annotations()
                
        self.browser.observe(on_browser_selection, names='selected')
        self.browser.selected = self.browser.items[0]
        
        
    def init_canvas(self, xmin=None, xmax=None):
        # Create plot with box select tool AND tap tool
        self.init_main_canvas(xmin, xmax, sizing_mode='stretch_both',
                              tools="xbox_select,tap,xbox_zoom,pan,save,reset,undo")
        
        self.init_nav_canvas(sizing_mode='stretch_width', height=100,
                             tools="xbox_select,xwheel_zoom,reset",
                             toolbar_sticky=False)
        
        self.init_atmos_canvas(sizing_mode='stretch_both')
        
        self.nav.line('x', 'y', source=self.nav_source, color='dimgrey')

        self.obs_line = self.fig.line('x', 'y', source=self.fig_line_source, color='lightgrey')
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
        self.fig.add_layout(self.llist_labels)
        
        if not self.obs.normalised:
            self.spline_line = self.fig.line('x', 'y', source=self.cntm_source, color='red', alpha=0.5)
            self.spline_knots = self.fig.scatter('x', 'y', source=self.spline_source, size=6, color='red')
            draw_tool = PointDrawTool(renderers=[self.spline_knots])
            self.fig.add_tools(draw_tool)

        self.atmos_scatter = self.atmos_fig.scatter('x', 'y', source=self.atmos_scatter_source, size=4,
                                                    color='black', nonselection_fill_color='lightgray', nonselection_alpha=1,
                                                    nonselection_line_color='none', selection_color='cornflowerblue')

        self.atmos_line = self.atmos_fig.line('x', 'y', source=self.atmos_line_source, color='red')
        self.atmos_interp_scatter = self.atmos_fig.scatter('x', 'y', source=self.atmos_interp_source, size=15,
                                                           color='red', alpha=0.5)
        
        # TOOLTIPS = [
        #     ('λ', '$x{0.000}'),
        #     ('[Fe/H]', '@FeH{0.00}'),
        #     ('vmic', '@vmic{0.00}'),
        # ]
        
        hover = HoverTool(renderers=[self.llist_line, self.synth_line],
                  tooltips=None, mode='mouse')  # You can add tooltips if needed
        
        self.fig.add_tools(hover)
        self.fig.legend.location = "bottom_left"
        self.fig.legend.spacing = 0

        def synth_selection(attr, old, new):
            selected = new
            if not selected: return
            idx = selected[0] # last clicked line
            rid = self.synth_source.data['run_ids'][idx]
            item = [item for item in self.browser.items if id(item.obj)==rid]
            if item:
                self.browser.selected = item[0]
            else:
                self.browser.selected = None

        self.synth_source.selected.on_change('indices', synth_selection)

        def synth_line_selection(attr, old, new):
            selected = new
            if not selected: return
            idx = selected[0] # last clicked line
            run = self.current
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
            self.atmos_fig,
            sizing_mode='stretch_both',
        )
        
        # Create main layout with plot on left, sidebar on right
        canvas = ServerManager.iframe_server(bk_canvas)
        canvas.add_class('canvas')
        return canvas


    @log_callback_errors
    def plot_synth(self, studio, name, launch_function, *args, runner=None):
        run = launch_function(self, name, *args, runner=runner)
        if run:
            crop_params = getattr(run, 'fit_params', {})
            x, y = run.crop(**crop_params)

            # run data is always normalised, so we need to transform
            # to view coordinates if the view is set to unnormalised
            y = conditional_denormalisation(x, y)

            @update_source_data(self.synth_source)
            def draw_synth(data):
                color_cylce = ColorCycle(saved_settings['bokeh palette'])
                color = color_cylce.get_new_color(data['colors'])
                data['xs'].append(x)
                data['ys'].append(y)
                data['colors'].append(color)
                data['labels'].append(run.label)
                data['vmic'].append(run.nml['atmos_params']['vmic'])
                data['FeH'].append(run.nml['atmos_params']['FeH'])
                data['atmos'].append(run.nml['atmos_params']['atmos_file'])
                data['run_ids'].append(id(run))
                data['alphas'].append(1.0)

            draw_synth()
            #self.set_synth_controls(run)
            return run
        

    def set_synth_controls(self, run):
        rid = id(run)
        item = self.browser.get_item_by_object(run)
        widget_label = item.label.label

        def update_label(change):
            @update_source_data(self.synth_source)
            def update(data):
                i = data['run_ids'].index(rid)
                data['labels'][i] = change['new']
            update()

        widget_label.observe(update_label, names='description')
        
    def init_main_canvas(self, xmin=None, xmax=None, xlabel='Wavelength [Å]', ylabel='Flux', **kwargs):
        p = figure(output_backend='svg',
                x_range=(0,1), # for some reason this is crucial
                y_range=(-0.5, 1.05), # for some reason this is crucial
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
            self.update_data_views()

        p.x_range.on_change('start', xrange_update_callback)
        p.x_range.on_change('end', xrange_update_callback)
        p.x_range.on_change('start', self.yrange_update_callback)
        p.x_range.on_change('end', self.yrange_update_callback)
        p.x_range.on_change('start', self.update_annotations)
        p.x_range.on_change('end', self.update_annotations)

        xmin = xmin or self.xdata[0]
        xmax = xmax or self.xdata[-1]
        p.x_range.start = xmin
        p.x_range.end   = xmax

    def update_data_views(self):
        slc1 = get_crop_range(self.xdata, self.fig.x_range, nmax=10000, margin=0.01)
        slc2 = get_crop_range(self.xdata, self.fig.x_range, nmax=1000, margin=0.01)

        # ydata is always unnormalised, so we need to normalise if view is normed
        y1 = conditional_normalisation(self.xdata[slc1], self.ydata[slc1])
        y2 = conditional_normalisation(self.xdata[slc2], self.ydata[slc2])
        self.fig_line_source.data = dict(x=self.xdata[slc1], y=y1)
        self.fig_source.data = dict(x=self.xdata[slc2], y=y2)
        self.cntm_source.data = dict(x=self.xdata[slc1], y=self.cntm[slc1])
        self.fig_source.selected.indices = []

    def reset_yaxis_and_annotations(self):
        logging.info('Updating yrange and line list annotations')
        ServerManager.server.next_tick_update(self.yrange_update_callback)
        ServerManager.server.next_tick_update(self.update_annotations)
    
    def yrange_update_callback(self, attr, old, new):
        yvis = self.fig_source.data['y']
        xs = self.synth_source.data['xs']
        ys = self.synth_source.data['ys']
        
        for x, y in zip(xs, ys):
            slc = get_crop_range(x, self.fig.x_range)
            yvis = np.concatenate((yvis, y[slc]))

        if len(yvis) > 0:
            reset_yaxis(self.fig, yvis)

    def update_annotations(self, attr, old, new):
        logging.info('Updating line list annotations')
        clear_data(self.llist_source)
        obj = self.current
        if isinstance(obj, StudioSynth) and hasattr(obj, 'llist'):
            logging.info(f'Plotting line list for {obj.label}')
            fig_width = max(self.fig.inner_width, self.atmos_fig.inner_width)
            ymin, ymax = obj.plot_llist(self.fig, self.llist_source, fig_width)
            if ymin and ymax:
                reset_yaxis(self.fig, [ymin, ymax])

    def init_nav_canvas(self, **kwargs):
        x0 = self.xdata[0]
        x1 = self.xdata[-1]
        p = figure(output_backend='svg',
                   toolbar_location='above',
                   x_axis_location='above',
                   x_range=(0,1), # for some reason this is crucial
                   y_range=(-0.5, 1.05), # for some reason this is crucial
                   **kwargs)
        self.nav = p
        p.yaxis.visible = False
        p.xaxis.axis_line_color = None
        p.xaxis.major_tick_line_color = None
        p.xaxis.minor_tick_line_color = None
        p.toolbar.logo = None
        tools = {type(tool):tool for tool in p.toolbar.tools}
        p.toolbar.active_scroll = tools[WheelZoomTool]

        def nav_xrange_update_callback(attr, new, old):
            slc = get_crop_range(self.xdata, p.x_range, nmax=10000)
            self.nav_source.data = dict(x=self.xdata[slc], y=self.ydata[slc])
            reset_yaxis(self.nav, self.ydata[slc])

        p.x_range = Range1d(0, 1, bounds=(x0, x1))
        p.x_range.on_change('start', nav_xrange_update_callback)
        p.x_range.on_change('end', nav_xrange_update_callback)
        p.x_range.start = x0
        p.x_range.end = x1
            
        @log_callback_errors
        def nav_drag_callback(event):
            geo = event.geometry
            xmin = geo['x0']
            xmax = geo['x1']
            reset_xaxis(self.nav, [xmin, xmax], margin=0.05)
            reset_xaxis(self.fig, [xmin, xmax])
            
        p.on_event(SelectionGeometry, nav_drag_callback)

    def init_atmos_canvas(self, **kwargs):
        p = figure(output_backend='webgl', **kwargs)
        self.atmos_fig = p
        p.toolbar.logo = None

    @property
    def xdata(self):
        return self.obs.xdata
    
    @property
    def ydata(self):
        return self.obs.ydata
    
    @property
    def cntm(self):
        return self.obs.cntm
    
    @cntm.setter
    def cntm(self, value):
        self.obs.cntm = value
    
    @property
    def current(self):
        return self.browser.selected.obj if self.browser.selected else None

    def get_profile_results(self):
        return get_profile_results()
    
    def add_object(self, name, func, *args, parent=None, **kwargs):
        # first we create a dummy object and add it to the browser
        dummy = SimpleNamespace()
        uname = self.browser.get_unique_label(name)

        if parent:
            item = parent.add_child(dummy, label=uname)
        else:
            item = self.browser.add_object(dummy, label=name)

        # then we define the async callback that will run the function
        @execute_in_thread(cancel=item.delete, spinner=item.icon, progress=item.progress)
        def async_callback(execution):
            obj = func(self, uname, *args, runner=execution, **kwargs)
            if obj and not execution.cancelled:
                # if the function returned an object, we replace the dummy with it
                item.replace_object(obj)
            else:
                # if the function did not return an object, we remove the dummy
                self.browser.remove_object(dummy)

        async_callback()