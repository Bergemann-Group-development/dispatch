import ipywidgets as widgets
from tumult.auxfuncs import get_ion_weight, EK
from tumult.linelist import LineList
from .settings import saved_settings
from .ServerManager import reset_xaxis, update_source_data, clear_data
from . import ServerManager
from .auxfuncs import *
from bokeh.models import ColumnDataSource as CDS
import numpy as np
import pandas as pd
import ipyformkit as ifk

fields = {('sort by', 'window size [Å]'):
          (('wavelength', 'approx. strength (descending)', 'loggf (descending)', 'evi (ascending)'), 5)}

check = {'window size [Å]': lambda d: d['window size [Å]'] > 0}
inputs = ifk.Form(fields, check=check)
sort_wid = inputs.widgets_dict['sort by'].wid
range_wid = inputs.widgets_dict['window size [Å]'].wid

ion_wid = widgets.Select(options=[])
line_wid = widgets.Select(option=[])
LineNavigator_btn = widgets.Button(description='Line Navigator')
ion_wid.layout.width = 'fit-content'
line_wid.layout.width = '100%'
selects = widgets.HBox([ion_wid, line_wid])

# Keep track of the current timeout callback
current_timeout_callback = None
dummy_T = 5000

def init_LineNavigator(obs, studio):
  master = LineList.from_tumult_format(saved_settings['master linelist'])
  master['strength'] = master.lam0 * np.exp(-master.evi*EK/dummy_T) * 10**master.loggf

  indicator_source = CDS(data=dict(x=[], y=[]))
  indicator = studio.fig.line('x', 'y', source=indicator_source, line_width=0.5, color='red')

  ions = master.index.get_level_values(0).unique()
  weights = [get_ion_weight(ion) for ion in ions]
  ions = ions[np.argsort(weights)]
  ion_wid.options = ions

  @log_callback_errors
  def ion_selected(change):
    if change['new'] is None: return
    ion = ion_wid.value
    df = master.loc[ion]

    df = group_by_distance(df, 'lam0', 0.02)
    df = df.groupby(['evi', 'group']).apply(superline, include_groups=False)

    xmin = obs.xdata[0]
    xmax = obs.xdata[-1]
    mask = (df.lam0 >= xmin) & (df.lam0 <= xmax)
    df = df[mask]

    sort = sort_wid.value
    if sort.startswith('loggf'):
      df = df.sort_values('loggf', ascending=False)
    elif sort.startswith('evi'):
      df = df.sort_values('evi', ascending=True)
    elif sort.startswith('approx.'):
      df = df.sort_values('strength', ascending=False)
    else:
      df = df.sort_values('lam0', ascending=True)

    line_wid.options = df.lam0.map('{:.3f}'.format)

  @log_callback_errors
  def line_selected(change):
    if line_wid.value is None: return
    lam0 = float(line_wid.value)
    dlam = range_wid.value
    xmin = lam0 - dlam/2
    xmax = lam0 + dlam/2

    def update():
      reset_xaxis(studio.fig, [xmin, xmax])
      reset_xaxis(studio.nav, [xmin, xmax], margin=2)
      show_vertical_line_with_timeout(studio.fig, indicator_source, lam0)

    ServerManager.server.next_tick_update(update)

  ion_wid.observe(ion_selected, names='value')
  sort_wid.observe(ion_selected, names='value')
  line_wid.observe(line_selected, names='value')
  range_wid.observe(line_selected, names='value')

  LineNavigator = widgets.VBox([selects, inputs.vbox])
  LineNavigator.add_class('line-navigator')
  LineNavigator.indicator = indicator
  LineNavigator.source = indicator_source
  obs.master = master
  return LineNavigator


def show_vertical_line_with_timeout(p, source, xpos):
    """Show vertical line and set/reset timeout to clear it"""
    global current_timeout_callback
    doc = ServerManager.server.doc
    
    # Cancel existing timeout if there is one
    if current_timeout_callback is not None:
        doc.remove_timeout_callback(current_timeout_callback)
    
    # Show the line
    ymin = p.y_range.start
    ymax = p.y_range.end

    @update_source_data(source)
    def update_indicator_source(data):
        data['x'] = [xpos, xpos]
        data['y'] = [ymin, ymax]
    update_indicator_source()
    
    def clear_line():
      global current_timeout_callback
      clear_data(source)
      current_timeout_callback = None
    
    # Schedule removal after 1 second (1000ms)
    current_timeout_callback = doc.add_timeout_callback(clear_line, 1000)

def superline(group):
    w = group['strength']
    return pd.Series({
        'lam0': np.average(group['lam0'], weights=w),
        'loggf': np.log10(np.sum(10**group['loggf'])),
        'strength': np.sum(group['strength']),
    })

def group_by_distance(df, column, max_distance):
    """Group rows where consecutive values are within max_distance of each other"""
    df_sorted = df.sort_values(column).reset_index(drop=True)
    
    # Calculate differences between consecutive values
    diffs = df_sorted[column].diff().fillna(0)
    
    # Create group boundaries where difference > max_distance
    group_breaks = (diffs > max_distance).cumsum()
    
    df_sorted['group'] = group_breaks
    return df_sorted