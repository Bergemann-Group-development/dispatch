import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from inspect import signature
import numpy as np
import pandas as pd
from scipy.ndimage import convolve
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

plt.style.use('default')
plt.style.use(os.path.join(current_dir, 'Bergemann2020.mplstyle'))

plot_kwargs = signature(plt.plot).parameters
line2d_kwargs = signature(matplotlib.lines.Line2D).parameters
artist_properties = {
    'alpha', 'zorder', 'label', 'visible', 'animated', 'clip_on', 
    'clip_box', 'clip_path', 'transform', 'agg_filter', 'figure'}
combined_kwargs = {**plot_kwargs, **line2d_kwargs}
plt_kwargs = set(combined_kwargs.keys()).union(artist_properties)

#=====================================================================
def set_title(title, default=None):
    if title is None:
        title = default

    if title is not None:
        plt.title(title)

#=====================================================================
def set_legend(**kwargs):
    if 'label' in kwargs:
        plt.legend(frameon=False)

#=====================================================================
def subplot_labels(ax=None, xlabel=None, ylabel=None, title=None):
    if ax == None: ax = plt.gca()
    sub = ax.get_subplotspec()
    grid = sub.get_gridspec()
    if xlabel is not None:
        if sub.rowspan.stop == grid._nrows:
            ax.set_xlabel(xlabel)
    if ylabel is not None:
        if sub.colspan.start == 0:
            ax.set_ylabel(ylabel)
    if title is not None:
        if sub.rowspan.start == 0:
            ax.set_title(title)


#=====================================================================
def adjust_max_xticks(ax=None):
    ax = ax or plt.gca()
    transform = ax.transData.inverted()
    xticklabels = ax.get_xticklabels()
    widths = [label.get_window_extent().transformed(transform).width for label in xticklabels]
    axwidth = ax.get_window_extent().transformed(transform).width
    if np.sum(widths) > axwidth:
        n = np.floor(axwidth / np.mean(widths))
        ax.xaxis.set_major_locator(plt.MaxNLocator(n))

#=====================================================================
def comparison_scatter(xx, yy, marker='.', rasterized=True, **kwargs):
  p = plt.scatter(xx, yy, marker=marker, rasterized=rasterized, **kwargs)
  xmin = min(np.nanmin(xx), np.nanmin(yy))
  xmax = max(np.nanmax(xx), np.nanmax(yy))
  plt.plot([xmin, xmax], [xmin, xmax], zorder=-10, ls='--', color='lightgrey')
  return p

#=====================================================================
def scatter_density(xx, yy, radius=0, cmap='binary', vmin=0, vmax=None, aspect='auto',
            xmin=None, xmax=None, ymin=None, ymax=None, interpolation='nearest', **kwargs):
  #-----------------------------------------------------------------
  
  mask = np.isnan(xx) | np.isnan(yy) | np.isinf(xx) | np.isinf(yy)
  xx = xx[~mask]
  yy = yy[~mask]
  
  if xmin is None: xmin = xx.min()
  if xmax is None: xmax = xx.max()
  if ymin is None: ymin = yy.min()
  if ymax is None: ymax = yy.max()

  fig = plt.gcf()
  nx, ny = fig.get_size_inches() * fig.get_dpi()

  # Get the bounding box of the axes in figure-relative coordinates
  bbox = plt.gca().get_position()

  # Compute the size of the axes in inches
  nx = bbox.width * nx
  ny = bbox.height * ny

  nx = np.ceil(nx).astype(int)
  ny = np.ceil(ny).astype(int)
  dx = (xmax - xmin) / (nx-1)
  dy = (ymax - ymin) / (ny-1)
  
  ix = (xx - xmin) / dx
  iy = (yy - ymin) / dy
  iix = np.floor(ix).astype(int)
  iiy = np.floor(iy).astype(int)
  mask = (iix >= 0) & (iix < nx-1) & (iiy >= 0) & (iiy < ny-1)
  ix = ix[mask]
  iy = iy[mask]
  iix = iix[mask]
  iiy = iiy[mask]
  fx = ix - iix
  fy = iy - iiy

  # Image
  image = np.zeros((nx, ny))  # Replace with your actual image
  np.add.at(image, (iix-1, iiy-1), (1-fx)*(1-fy))
  np.add.at(image, (iix, iiy-1), fx*(1-fy))
  np.add.at(image, (iix-1, iiy), (1-fx)*fy)
  np.add.at(image, (iix, iiy), fx*fy)

  if radius > 0:
    # Circle convolution kernel
    x = np.arange(-radius, radius+1)
    y = np.arange(-radius, radius+1)
    x, y = np.meshgrid(x, y)
    circle = np.sqrt(x**2 + y**2)
    circle = np.interp(circle.ravel(), [radius, radius+1], [1, 0]).reshape(circle.shape)
    circle = circle / circle.sum()
    image = convolve(image, circle)

  with np.errstate(divide='ignore', invalid='ignore'):
    image = np.log(image)

  if vmax is None: vmax = image.max()
  return plt.imshow(image.T, origin='lower', extent=(xmin, xmax, ymin, ymax),
                    vmin=vmin, vmax=vmax, aspect=aspect, cmap=cmap, 
                    interpolation=interpolation, **kwargs)

#=====================================================================
def alpha_cmap(cmap_name, clip=None):
    # Load the original colormap and extract RGB values
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0, 1, cmap.N))
    
    # Separate RGB and calculate alpha based on lightness
    rgba_colors = np.zeros_like(colors)
    rgba_colors[:, :3] = colors[:, :3]  # Copy the RGB values
    alphas = 1 - colors[:, :3].mean(axis=1)  # Higher lightness -> lower alpha

    # Ensure alpha is in [0, 1] range
    if clip is not None:
      alphas = np.clip(alphas, clip, 1.0)  # Preventing 0 alpha for visibility
    
    # Adjust RGB based on alpha to blend correctly over white
    rgba_colors[:, 3] = alphas # Set the alpha channel
    alphas = alphas[:, None]
    rgba_colors[:, 0:3] = (colors[:, 0:3] - 1 + alphas) / alphas

    # Clip RGB to stay within valid range [0, 1]
    rgba_colors[:, :3] = np.clip(rgba_colors[:, :3], 0, 1)

    # Create a new colormap with the adjusted RGBA values
    new_cmap = LinearSegmentedColormap.from_list('custom_map', rgba_colors)
    return new_cmap

#=====================================================================
def create_centered_colormap(cmap_name, centers):
    # Get the colormap from the specified name
    cmap = plt.get_cmap(cmap_name)
    
    # Define the boundaries for the levels
    bounds = np.array(centers)
    diff = np.diff(centers)
    bounds = np.concatenate([bounds - diff[0]/2, [bounds[-1] + diff[-1]/2]])
    num_levels = len(bounds)
    
    # Create the normalization
    norm = BoundaryNorm(bounds, num_levels)
    
    # Create the stepwise colormap
    stepwise_cmap = ListedColormap(cmap(np.linspace(0, 1, num_levels)))
    
    return stepwise_cmap, norm

#=====================================================================
def plot_arrow(start, end, color=None, linewidth=None, markersize=None):
    """
    Plots an arrow from `start` to `end` with a line and a rotated 'v' marker.

    Parameters:
    - start: tuple (x0, y0), starting point of the arrow.
    - end: tuple (x1, y1), ending point of the arrow.
    - color: color of the line and marker.
    - linewidth: width of the line.
    - markersize: size of the 'v' marker.
    """
    # Calculate the direction vector in data coordinates
    dx, dy = end[0] - start[0], end[1] - start[1]
    
    # Get axis data limits and aspect ratio correction
    ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    data_ratio = ax.get_data_ratio()  # Ratio of y-unit to x-unit

    # Scale dy by the data ratio to account for aspect ratio
    dy_corrected = dy / data_ratio

    # Calculate the angle in degrees from the corrected vector
    angle = np.degrees(np.arctan2(dy_corrected, dx)) - 90

    # Plot the line segment
    p, = plt.plot([start[0], end[0]], [start[1], end[1]], color=color, linewidth=linewidth)
    plt.plot(end[0], end[1], marker=(3, 0, angle), color=p.get_color(), markersize=markersize)

#=====================================================================
def plot_line_annotations(xx, labels, weights=None, ha='center', va='bottom',
                      rotation=90, fontfamily='serif', color='k', xmin=None, xmax=None,
                      ymin=None, ymax=None, linewidth=0.8, ntop=-1, focus='strongest', **kwargs):
    #-----------------------------------------------------------------
    ax = plt.gca()
    xmn, xmx = ax.get_xlim()
    ymn, ymx = ax.get_ylim()
    if xmin is None: xmin = xmn
    if xmax is None: xmax = xmx
    if ymin is None: ymin = ymn
    if ymax is None: ymax = ymx
    h = 0.1 * (ymax - ymin)

    mask = (xx > xmin) & (xx < xmax)
    xx = xx[mask]
    labels = np.array(labels)[mask]
    if weights is not None: weights = weights[mask]

    plt.gcf().draw_without_rendering()
    t = ax.text(xmin, ymin, 'Mg IV', ha=ha, va=va, rotation=rotation, fontfamily=fontfamily, **kwargs)
    bbox = t.get_window_extent().transformed(ax.transData.inverted())
    t.remove()

    if rotation == 90:
        pad = 1.4
    else:
        pad = 1.2

    dxmin = bbox.width*pad
    newn = int((xmax - xmin) / dxmin)

    if ntop == -1: ntop = newn
    if ntop is not None and weights is not None:
        
        if focus == 'strongest':
            order = np.argsort(weights)[::-1][:ntop]
            mask = np.isin(np.arange(xx.size), order)
        elif focus == 'bins':
            # 1. Define bin edges and assign indices
            edges = np.histogram_bin_edges(xx, bins=ntop)
            indices = np.digitize(xx, edges)
            indices[indices == ntop + 1] = ntop  # Include the rightmost edge in the final bin

            # 2. Find the maximum weights value within each bin
            max_vals = np.full(ntop + 2, -np.inf)
            np.maximum.at(max_vals, indices, weights)
            mask = weights == max_vals[indices]

        labels = labels[mask]
        xx = xx[mask]

    n = xx.size

    newn = max(newn, n)
    nn = np.ones(n)

    newx = np.linspace(xmin, xmax, newn+2)[1:-1]

    xlist = list(xx)
    binlist = list(newx)

    newx = []
    for x in xlist:
        i = np.argmin(np.abs(binlist - x))
        newx.append(binlist[i])
        del binlist[i]

    oldx = xx
    newx = np.sort(newx)

    xl = np.vstack([newx, newx, oldx, oldx])
    yl = ymax*nn + np.vstack([h*nn, 0.8*h*nn, 0.2*h*nn, 0*nn])

    plt.plot(xl, yl, color=color, linewidth=linewidth)

    for x, label in zip(newx, labels):
        plt.text(x, ymax+1.2*h, label, ha=ha, va=va, rotation=rotation,
                  fontfamily=fontfamily, color=color, **kwargs)

    if rotation == 90:
        pad = 1.5
    else:
        pad = 2.5

    plt.ylim(top=ymax+1.2*h + bbox.height*pad)


def quickplot(x, y, ax=None, color=None, trendline=None, xerr=None, yerr=None,
              xline=None, yline=None, xlabel=None, ylabel=None, margins=None,
              xlog=False, ylog=False, trend_color=None, trend_alpha=1, 
              **kwargs):

  # Error propagation for log scaling
  if xlog:
      xerr = np.asarray(xerr) / x if xerr is not None else None
      x = np.log(x)

  if ylog:
      yerr = np.asarray(yerr) / y if yerr is not None else None
      y = np.log(y)

  ax = ax or plt.gca()
  pt = ax.errorbar(x, y, xerr=xerr, yerr=yerr, color=color, **kwargs)

  if trendline:
    c = trend_color or pt.lines[0]._color
    w = 1/yerr**2 if yerr is not None else None
    p = np.polyfit(x, y, w=w, deg=trendline)
    newx = np.linspace(np.min(x), np.max(x), 200)
    f = np.polyval(p, newx)
    ax.plot(newx, f, color=c, alpha=trend_alpha, ls='--')

  xline = xline or []
  yline = yline or []
  xline = [xline] if not hasattr(xline, '__iter__') else xline
  yline = [yline] if not hasattr(yline, '__iter__') else yline

  [ax.axhline(y, color='lightgrey', ls='--') for y in yline]
  [ax.axvline(x, color='lightgrey', ls='--') for x in xline]

  if xlabel: ax.set_xlabel(xlabel)
  if ylabel: ax.set_ylabel(ylabel)
  if margins: ax.margins(*margins)

  return pt
