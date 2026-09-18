"""The plot palette and Bokeh figure styling shared across features.

The colours mirror the page's CSS tokens (STYLE_REFERENCE.md §2): the dark surfaces, the
one blue accent, and the line palette the spectra are drawn in. Every Bokeh figure in the
app (the canvas, the filter view, the diagnostic panels, the curve of growth) is styled by
:func:`style_dark` so it reads as one surface.
"""

BG = "#1e1e1e"
SURFACE = "#2a2a2a"
BORDER = "#3a3a3a"
GRID = "#2c2c2c"
MUTED = "#888888"
ACCENT = "#4a90e2"
# distinct from ACCENT so the abundance best-fit synthetic reads clearly against the (ACCENT)
# EW hatch it is drawn over
SYNTH_FIT = "#3ec46d"
LINE = "#cfd2d6"

# fixed line-colour palette; PALETTE[0] == LINE is the default (gray) object colour,
# PALETTE[1:] are the non-gray colours assigned to tumult runs. The first run colour is a
# cornflower blue kept distinct from ACCENT (#4a90e2) so a run never matches the UI accent.
PALETTE = ["#cfd2d6", "#7b8cde", "#50c878", "#e0a030", "#e05c5c", "#b084e8", "#34c9a0"]
RUN_COLORS = PALETTE[1:]

# objects are distinguished by dash; runs by colour (then dash once colours run out)
DASHES = ["solid", "dashed", "dotted", "dotdash", "dashdot"]
# a Bokeh line dash -> the CSS border style the sidebar swatch draws it with
DASH_CSS = {"solid": "solid", "dashed": "dashed", "dotted": "dotted",
            "dotdash": "dashed", "dashdot": "dotted"}


def style_dark(p):
    """Give a Bokeh figure the app's dark theme and no toolbar."""
    p.background_fill_color = BG
    p.border_fill_color = BG  # dark plot surround to match the darkened work surface
    p.outline_line_color = BORDER
    p.grid.grid_line_color = GRID
    p.axis.axis_line_color = BORDER
    p.axis.major_tick_line_color = BORDER
    p.axis.minor_tick_line_color = None
    p.axis.major_label_text_color = MUTED
    p.axis.axis_label_text_color = MUTED
    # no Bokeh toolbar anywhere: the tools stay usable (wheel zoom, pan, tap, hover, the
    # active drag the tools switch) and the page's own toolbar holds every button
    p.toolbar_location = None
