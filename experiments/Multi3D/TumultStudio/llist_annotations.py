from .auxfuncs import *
from .widgets import EditableTable
import ipywidgets as widgets
import numpy as np
import logging
from .ServerManager import update_source_data, clear_data

fmt = {
    'lam0': lambda x: f'{x:.3f}',
    'loggf': lambda x: f'{x:.3f}',
    'evi': lambda x: f'{x:.3f}',
    'ga': lambda x: f'{x:.2e}'}

head = ['ion', 'lam0', 'loggf', 'evi', 'ga', 'gv', 'gq', 'i', 'j', 'iso']

#=====================================================================
def get_llist_table(llist):
    #-----------------------------------------------------------------
    """
    Returns a table of the line list with formatted values.
    """
    html_widget = EditableTable(llist, columns=head, index=False, formatters=fmt)
    html_widget.layout.max_height = 'none'
    html_widget.layout.width = '100%'
    return html_widget

#=====================================================================
def plot_line_annotations(run, p, source, fig_width):
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
        clear_data(source)
        logging.info("No line annotations in the current range.")
        return None, None

    xx = xx[mask]
    labels = np.array(labels)[mask]
    weights = weights[mask]

    lxmin = np.min(xx)
    lxmax = np.max(xx)
    wratio = (lxmax - lxmin) / (pxmax - pxmin)
    wpix = fig_width * wratio
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

    @update_source_data(source)
    def update_llist_source(data):
        data['xs'] += [x for x in xl.T]
        data['ys'] += [y for y in yl.T]
        data['x'] += [x[-1] for x in xl.T]
        data['y'] += [y[-1]+0.1*h for y in yl.T]
        data['text'] += labels.tolist()
    update_llist_source()

    return pymin, pymax + 3*h