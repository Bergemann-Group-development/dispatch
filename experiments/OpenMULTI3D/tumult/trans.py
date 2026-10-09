"""Bound-bound transitions of a model atom, and their spectra in a run."""
import numpy as np
from .const import CC
from .plotting_scripts import set_legend, set_title
from scipy.signal import argrelextrema


class Transition(object):
    """Bound-bound transition ``kr`` (from 0) of an atom, row ``kr`` of its ``bbt``
    table, read through ``rows`` (``atom.transitions``): the table's columns
    (``lam0``, ``loggf``, ``i``, ``j``, ``nq``, ...) read as attributes. In a run,
    ``ib``/``ir`` index its wavelength range on the run's grid
    (:meth:`tumult.atom.Atom.load_sp`)."""
    def __init__(self, rows, kr):
        self.rows = rows
        self.kr = kr

    def __getattr__(self, attr):
        rows = self.__dict__.get('rows')       # absent while a copy is being built
        if rows is not None and attr in rows.columns:
            return rows.columns[attr][self.kr]
        raise AttributeError(f'Transition has no attribute {attr!r}')

    def __repr__(self):
        return self.atom.bbt.iloc[self.kr].__repr__()

    @property
    def atom(self):
        return self.rows.atom

    @property
    def lower(self):
        """The lower level, a :class:`tumult.atom.Level`."""
        return self.rows.levels[int(self.i)]

    @property
    def upper(self):
        """The upper level, a :class:`tumult.atom.Level`."""
        return self.rows.levels[int(self.j)]

    @property
    def run(self):
        return self.atom.run

    @property
    def title(self):
        return 'Line = {:.2f} [Å]'.format(self.lam0)

    #=====================================================================
    def crop(self, **kwargs):
        #-----------------------------------------------------------------
        spec = self.run.crop(line=self, **kwargs)
        return spec
    
    #=====================================================================
    def get_clv(self, dlam=None, **kwargs):
        """Equivalent width [mÅ] at each μ of the run (:meth:`Run.get_clv`), over the
        line's own wavelength range or, with ``dlam``, over lam0 ± dlam [Å]."""
        #-----------------------------------------------------------------
        if dlam is None:
            kwargs['line'] = self
        else:
            kwargs['window'] = (self.lam0 - dlam, self.lam0 + dlam)
        xx, yy = self.run.get_clv(**kwargs)
        return xx, yy

    #=====================================================================
    def plot_clv(self, dlam=None, LTE=False, title=None, **kwargs):
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        xx, yy = self.get_clv(dlam=dlam, LTE=LTE)

        if not 'linestyle' in kwargs:
            kwargs['linestyle'] = 'none'

        if not 'marker' in kwargs:
            kwargs['marker'] = '.'

        
        plt.xlabel(r'μ')
        plt.ylabel('Equivalent width [mÅ]')

        set_title(title, default=self.title+', ' + self.atom.get_title())
        set_legend(**kwargs)

        return plt.plot(xx, yy, **kwargs)

    #=====================================================================
    def get_xx(self, dlam, velocity_units=False, **kwargs):
        #-----------------------------------------------------------------

        xx = self.run.lam
        if dlam is not None and dlam is not False:
            xx = xx - self.lam0

        xx, mask = self.run.get_xx(xx, slc=slice(self.ib, self.ir), **kwargs)

        if velocity_units:
            xx = xx / self.lam0
            if dlam is None or dlam is False:
                xx = xx - 1
            xx = xx * CC * 1e-5

        return xx, mask
        
    #=====================================================================
    def calc_blueshift(self, dvel=4, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(velocity_units=True, norm=True, **kwargs)
        mask = (xx > -dvel) & (xx < dvel)
        xx = xx[mask]
        yy = yy[mask]
        p = np.polyfit(xx, yy, 2)
        pdiff = np.polyder(p)
        poly_root = np.roots(pdiff)
        blues = poly_root[0]
        depth = 1 - np.polyval(p, blues)
        return blues, depth

    #=====================================================================
    def plot_resolved(self, ang=0, norm=True, dlam=None, velocity_units=False, 
                vdop=None, vshift=None, LTE=False, title=None,
                color='k', linewidth=1, alpha=0.2, step=1, **kwargs):
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        if not self.run.nml['output_params']['save_resolved']:
            print('Resolved profiles were not saved!')
            return None

        xx, mask = self.get_xx(dlam, velocity_units=velocity_units)
        self.run.set_xaxis(xx, dlam=dlam, velocity_units=velocity_units)
        
        cnt = self.run.c3[...,ang]

        if LTE:
            yy = self.run.i3_lte[...,ang]
        else:
            yy = self.run.i3[...,ang]

        nnu = xx.size
        yy  = yy [mask, ::step, ::step]
        cnt = cnt[mask, ::step, ::step]

        if norm:
            yy = yy / np.mean(cnt, axis=(1,2))[:,None,None]
        
        plt.plot(xx, yy.reshape((nnu, -1)), color=color, linewidth=linewidth,
                 alpha=alpha, **kwargs)
        
        set_title(title, default=self.title+', '+self.atom.get_title())
        set_legend(**kwargs)
        
        if 'label' in kwargs:
            #remove duplicates
            handles, labels = plt.gca().get_legend_handles_labels()
            newLabels, newHandles = [], []
            for handle, label in zip(handles, labels):
                if label not in newLabels:
                    newLabels.append(label)
                    newHandles.append(handle)
            
            plt.legend(newHandles, newLabels, frameon=False)
        
        return xx, yy
    
    #=====================================================================
    def get_bisector(self, n=100, frac=0.99, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(**kwargs)
        ex = argrelextrema(yy, np.less_equal)[0]
        i0 = ex[np.argmin(np.abs(ex - self.i0 + self.ib))]
        ex = argrelextrema(yy, np.greater_equal)[0]
        i = np.searchsorted(ex, i0)
        ib = ex[i-1]
        ir = ex[i]
        ymin = np.max(yy[i0-1:i0+2])
        ymax = np.min([yy[ib], yy[ir]])
        ymax = ymin + (ymax - ymin)*frac
        newy = np.linspace(ymax, ymin, n)
        lwing = np.interp(newy, yy[i0:ib:-1], xx[i0:ib:-1])
        rwing = np.interp(newy, yy[i0:ir], xx[i0:ir])
        newx = 0.5*(lwing + rwing)
        return newx, newy