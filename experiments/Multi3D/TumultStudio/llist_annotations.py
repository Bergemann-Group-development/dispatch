from .auxfuncs import *
from .widgets import EditableTable
import ipywidgets as widgets
import numpy as np

fmt = {
    'lam0': lambda x: f'{x:.2f}',
    'loggf': lambda x: f'{x:.3f}',
    'evi': lambda x: f'{x:.3f}',
    'ga': lambda x: f'{x:.2e}'}

head = ['ion', 'iso', 'lam0', 'loggf', 'evi', 'ga', 'gv', 'gq']

#=====================================================================
def get_llist_table(run):
    #-----------------------------------------------------------------
    """
    Returns a table of the line list with formatted values.
    """
    html_widget = EditableTable(run.llist, columns=head, index=False, formatters=fmt)
    html_widget.layout.max_height = 'none'
    html_widget.layout.width = '100%'
    return html_widget

#=====================================================================
def remove_line_annotations(source):
    #-----------------------------------------------------------------
    cdata = source.data.copy()
    for key in source.data:
        cdata[key] = []
    source.data = cdata

#=====================================================================
def plot_line_annotations(p, run, source):
    #-----------------------------------------------------------------
    xx = run.llist.lam0.values
    labels = run.llist.ion.values
    weights = run.llist.strength.values

    pxmin = p.x_range.start
    pxmax = p.x_range.end
    pymin = p.y_range.start
    pymax = max(1.0, p.y_range.end)

    rxmin = run.lam[0]
    rxmax = run.lam[-1]

    h = 0.1 * (pymax - pymin)
    pmask = (xx > pxmin) & (xx < pxmax)
    rmask = (rxmin < xx) & (xx < rxmax)
    mask = pmask & rmask
    if not np.any(mask):
        remove_line_annotations(source)
        return

    xx = xx[mask]
    labels = np.array(labels)[mask]
    weights = weights[mask]

    lxmin = np.min(xx)
    lxmax = np.max(xx)
    wratio = (lxmax - lxmin) / (pxmax - pxmin)
    wpix = p.inner_width * wratio
    ntop = int(wpix / 20)  # 50 pixels per label

    arr = labels
    change_indices = np.where(arr[1:] != arr[:-1])[0] + 1
    starts = np.concatenate(([0], change_indices))
    ends = np.concatenate((change_indices, [len(arr)]))
    wts = np.array([np.sum(weights[s:e]) for s, e in zip(starts, ends)])

    # select the top ntop labels based on their weights
    tops = np.argsort(wts)[::-1][:ntop]
    mask = np.isin(np.arange(wts.size), tops)

    starts = starts[mask]
    ends = ends[mask]
    reduced_labels = arr[starts]
    label_repeats = ends - starts

    # remove consecutive duplicate labels
    oldx = []
    wxx = []
    for s, e in zip(starts, ends):
        oldx = np.concatenate((oldx, xx[s:e]))
        w = weights[s:e] / np.sum(weights[s:e])
        wxx.append(np.sum(xx[s:e]*w))

    n = reduced_labels.size
    newn = max(ntop, n)
    xbins = np.linspace(lxmin, lxmax, newn+2)[1:-1]
    xbins = list(xbins)

    newx = []
    for x in wxx:
        i = np.argmin(np.abs(np.array(xbins) - x))
        newx.append(xbins[i])
        del xbins[i] # remove the used bin

    newx = np.sort(newx)
    newx = np.repeat(newx, label_repeats)
    labels = np.repeat(reduced_labels, label_repeats)
    idx = np.where(newx[1:] == newx[:-1])[0] + 1
    labels[idx] = ''

    nn = np.ones(newx.size)
    xl = np.vstack([oldx, oldx, newx, newx])
    yl = pymax*nn + h*np.vstack([0*nn, 0.2*nn, 0.8*nn, 1*nn])

    current_data = source.data.copy()
    current_data['xs'] += [x for x in xl.T]
    current_data['ys'] += [y for y in yl.T]
    current_data['x'] += [x[-1] for x in xl.T]
    current_data['y'] += [y[-1]+0.1*h for y in yl.T]
    current_data['text'] += labels.tolist()
    source.data = current_data
    
    p.y_range.end = pymax + 3*h