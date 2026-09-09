"""Legacy Bokeh plotting helpers, kept on the GUI side of the package boundary."""
from bokeh.io import curdoc
from bokeh.themes import Theme
from bokeh.plotting import figure, show
from bokeh.io import output_notebook, push_notebook
from bokeh.models import LinearAxis, Range1d
from IPython import get_ipython
from IPython.display import display
from inspect import signature
from bokeh.resources import INLINE

custom_theme = Theme(json={
    'attrs': {
        'figure': {
            'width': 450,
            'height': 333,
            'background_fill_color': 'white',
            'outline_line_color': None,
        },
        'Axis': {
            'major_tick_line_color': 'black',
            'minor_tick_line_color': 'black',
            'major_tick_in': 10,
            'minor_tick_in': 5,
            'major_tick_out': 0,
            'minor_tick_out': 0,
            'major_tick_line_width': 0.8,
            'minor_tick_line_width': 0.6,
            'axis_line_color': 'black',
            'major_label_text_font_size': '10pt',
            'major_label_text_font': 'sans-serif',
            'axis_label_text_font_size': '10pt',
            'axis_label_text_font_style': 'normal',  # <-- This makes axis labels non-italic
        },
        'Legend': {
            'background_fill_alpha': 0.0,
            'border_line_width': 0.0,
            'label_text_font_size': '10pt',
            'label_text_font': 'sans-serif',
            'spacing': 10,
        },
        'Title': {
            'text_font_size': '12pt',
            'text_font': 'sans-serif',
        },
        'ColorBar': {
            'title_text_font_size': '12pt',
            'major_label_text_font_size': '10pt',
        },
        'Grid': {
            'grid_line_color': None
        },
        'Line': {
            'line_width': 2  # Set your desired default width here
        }
    },
    'palette': ['#006BA4', '#FF800E', '#ABABAB', '#595959', '#5F9ED1', 
                '#C85200', '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']
})

curdoc().theme = custom_theme
palette = custom_theme._json['palette']

ip = get_ipython()

class bokeh_wrapper:
  def __init__(self, tools='pan, reset', **kwargs):
    output_notebook(resources=INLINE, hide_banner=True)
    self.p = figure(tools=tools, output_backend='svg', **kwargs)
    self.p.add_layout(LinearAxis(), 'above')
    self.p.add_layout(LinearAxis(), 'right')
    # Optional: show ticks but not labels on the top/right
    self.p.axis[0].major_label_text_font_size = '0pt'  # Top
    self.p.axis[3].major_label_text_font_size = '0pt'  # Right
    # Optional: fine-tune the others if needed
    self.p.min_border = 10
    self.p.min_border_bottom = 30
    self.p.min_border_left = 30
    self.i_color = 0
    self.count = ip.execution_count
    self.kwargs = kwargs

  def figure(self, tools='pan, reset', **kwargs):
    if self.count != ip.execution_count:
      self.__init__(tools=tools, **self.kwargs)
    return self.p

  def plot(self, *args, lw=2, ls='solid', ms=None, color=None, tools='pan, reset', **kwargs):
    if self.count != ip.execution_count:
      self.__init__(tools=tools, **self.kwargs)

    if not color:
      color = palette[self.i_color]
      self.i_color = (self.i_color+1) % len(palette)

    if 'label' in kwargs:
      if kwargs['label']: kwargs['legend_label'] = kwargs['label']
      del kwargs['label']

    if lw:
      self.p.line(*args, color=color, line_width=lw, line_dash=ls, **kwargs)
    if ms:
      self.p.scatter(*args, size=ms, color=color, **kwargs)

    return self
  
  def xlabel(self, label):
    self.p.xaxis[-1].axis_label = label
    return self

  def ylabel(self, label):
    self.p.yaxis[0].axis_label = label
    return self

  def xlim(self, xmin, xmax):
    self.p.x_range = Range1d(start=xmin, end=xmax)
    return self

  def ylim(self, ymin, ymax):
    self.p.y_range = Range1d(start=ymin, end=ymax)
    return self
  
  def show(self):
    show(self.p, notebook_handle=True)

  def _repr_html_(self):
    self.show()

bkh = bokeh_wrapper()
bkh.plot_kwargs = signature(bkh.plot).parameters
