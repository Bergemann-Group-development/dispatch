import plotly.graph_objects as go
import plotly.io as pio
from IPython import get_ipython
from inspect import signature

# Define the custom template
custom_template = go.layout.Template(
    layout=dict(
        font=dict(
            family="sans-serif",
            size=10,
            color="black"
        ),
        xaxis=dict(
            tickfont=dict(size=12, color="black"),
            mirror="allticks",  # adds top axis line and ticks
            ticks="inside",     # makes the ticks point inward
            showline=True,
            linecolor="black",
            automargin=True,
            showgrid=False,
            ticklen=8,
            minor=dict({
              'ticks':'inside',
              'ticklen':4,
            })
        ),
        yaxis=dict(
            tickfont=dict(size=12, color="black"),
            mirror="allticks",  # adds top axis line and ticks
            ticks="inside",     # makes the ticks point inward
            showline=True,
            linecolor="black",
            automargin=True,
            showgrid=False,
            ticklen=8,
            minor=dict({
              'ticks':'inside',
              'ticklen':4,
            })
        ),
        legend=dict(
            x=0.95,  # Horizontal position (0 is left, 1 is right)
            y=0.05,  # Vertical position (0 is bottom, 1 is top)
            bgcolor="rgba(255, 255, 255, 0.5)",  # Optional: Set a semi-transparent background
            bordercolor="black",  # Optional: Add a border color
            borderwidth=0,        # Optional: Set border width
            xanchor="right",  # Anchor the legend box's right edge to the x position
            yanchor="bottom",  # Anchor the legend box's right edge to the x position
        ),
        width=500,  # width in pixels (approx. 5 inches at 100 dpi)
        height=375, # height in pixels (approx. 3.75 inches at 100 dpi)
        margin=dict(
            l=30,  # Left margin
            r=20,  # Right margin
            t=30,  # Top margin
            b=30   # Bottom margin
        ),
    )
)

# Register and activate the template
pio.templates["bergemann2020"] = custom_template
pio.templates.default = "bergemann2020"

palette = ['#006BA4', '#FF800E', '#ABABAB', '#595959', '#5F9ED1', 
           '#C85200', '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']

ip = get_ipython()

class plotly_wrapper:
  def __init__(self, width=450, height=333, **kwargs):
    self.p = go.Figure()
    self.p = self.p.update_layout(width=width, height=height, dragmode='pan', **kwargs)
    self.i_color = 0
    self.count = ip.execution_count
    kwargs['width'] = width
    kwargs['height'] = height
    self.kwargs = kwargs

  def figure(self, tools='pan, reset', **kwargs):
    if self.count != ip.execution_count:
      self.__init__(**self.kwargs)
    return self.p

  def plot(self, x, y, lw=2, line=None, marker=None, mode=None, color=None, label=None, **kwargs):
    if self.count != ip.execution_count:
      self.__init__(**self.kwargs)

    if not color:
      color = palette[self.i_color]
      self.i_color = (self.i_color+1) % len(palette)

    self.p.add_trace(go.Scatter(x=x, y=y, line=line, marker=marker, mode=mode, name=label, **kwargs))
    return self
  
  def xlabel(self, label):
    self.p.update_layout(xaxis_title=label)
    return self

  def ylabel(self, label):
    self.p.update_layout(yaxis_title=label)
    return self

  def xlim(self, xmin, xmax):
    self.p.update_layout(xaxis=dict(range=[xmin, xmax]))
    return self

  def ylim(self, ymin, ymax):
    self.p.update_layout(yaxis=dict(range=[ymin, ymax]))
    return self
  
  def show(self):
    self.p.show()

  def _repr_html_(self):
    self.show()

ptl = plotly_wrapper()
ptl.plot_kwargs = signature(ptl.plot).parameters