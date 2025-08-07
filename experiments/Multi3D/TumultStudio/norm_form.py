from .auxfuncs import *
from .ServerManager import *
from .settings import saved_settings
import ipyformkit as ifk
from tumult.linelist import LineList
from tumult.auxfuncs import pchip_mono
import copy
from scipy.signal import argrelextrema
from scipy.interpolate import make_splrep
import ipywidgets as widgets
from .LineMask import LineMask

fields = {'method': ('spline',),
          'peak distance': 5.0,
          'Normalise': None,
         }

norm_form = ifk.Form(fields)
norm_form.vbox.add_class('norm-form')
norm_button = norm_form.widgets_dict['Normalise'].wid
norm_checkbox = widgets.Checkbox(value=False, description='View normalised', indent=False)
norm_checkbox.add_class('norm-checkbox')
cntm_spline = None
    
def conditional_normalisation(x, y):
    if cntm_spline is None: return y
    if norm_checkbox.value:
        return y / cntm_spline(x)
    else:
        return y
    
def conditional_denormalisation(x, y):
    if cntm_spline is None: return y
    if not norm_checkbox.value:
        return y * cntm_spline(x)
    else:
        return y

def toggle_normalisation(x, y):
    y = conditional_normalisation(x, y)
    y = conditional_denormalisation(x, y)
    return y

def init_norm_checkbox(studio):
    if studio.normalised: return

    @log_callback_errors
    def on_norm_toggle(change):
        # A couple of things have to be updated
        # 1. The viewed slice of the spectrum, i.e. fig_source and fig_line_source
        # 2. The synthetic spectra, i.e. synth_source
        # 3. The annotations, which will simply be reset

        def norm_data(data):
          if 'x' in data:
            data['y'] = toggle_normalisation(data['x'], data['y'])
          elif 'xs' in data:
            for i in range(len(data['xs'])):
                data['ys'][i] = toggle_normalisation(data['xs'][i], data['ys'][i])

          # we do this here so it is run on the server
          norm = change['new']
          studio.spline_line.visible = not norm
          studio.spline_knots.visible = not norm

        update_data_source(studio.fig_source)(norm_data)()
        update_data_source(studio.fig_line_source)(norm_data)()
        update_data_source(studio.synth_source)(norm_data)()
        studio.reset_yaxis_and_annotations()

    norm_checkbox.observe(on_norm_toggle, names='value')
    studio.controls.children += (norm_checkbox,)

def init_norm_form(studio):
    if studio.normalised: return

    @log_callback_errors
    def on_norm_click(btn):
        values = norm_form.check_and_return_values()
        if values['method'] == 'spline':
            spline_normalisation(studio)

    @log_callback_errors
    def update_spline(attr, old, new):
        knots = studio.spline_source.data
        global cntm_spline
        order = np.argsort(knots['x'])
        kx = np.array(knots['x'])[order]
        ky = np.array(knots['y'])[order]
        cntm_spline = pchip_mono(kx, ky)
        studio.cntm = cntm_spline(studio.xdata)
        line_masks = studio._get_objects_by_type(LineMask)
        for lm in line_masks:
            lm.normalise(cntm_spline)
        studio.update_data_views()

    @update_data_source(studio.spline_source)
    def init_spline_knots(data):
        data['x'] = [studio.xdata[0], studio.xdata[-1]]
        data['y'] = [1, 1]

    norm_button.on_click(on_norm_click)
    studio.spline_source.on_change('data', update_spline)
    init_spline_knots()

def spline_normalisation(studio):
    values = norm_form.check_and_return_values()
    peak_distance = values['peak distance']
    dx = np.mean(np.diff(studio.xdata))
    nstep = int(peak_distance/dx)
    ex = argrelextrema(studio.ydata, np.greater_equal, order=nstep)[0]

    @update_data_source(studio.spline_source)
    def update_spline_knots(data):
        data['x'] = studio.xdata[ex]
        data['y'] = studio.ydata[ex]

    update_spline_knots()

    