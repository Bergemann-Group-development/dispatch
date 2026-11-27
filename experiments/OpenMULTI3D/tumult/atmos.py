"""A model atmosphere on a mesh: the one a run saved (:meth:`AtmosClass.from_run`) or a
model file read directly (:meth:`AtmosClass.from_file`), with plots against log tau_500
and height.

Arrays have shape (nx, ny, nz) in 3D and (nz,) in 1D, the deepest point first: ``temp``
[K], ``rho`` [g/cm3], ``ne`` [1/cm3], ``vx, vy, vz`` [cm/s], ``vmic`` [km/s], ``chi``
the extinction at 500 nm [1/cm], ``tau`` and ``ltau`` (log10) tau_500. ``xx``, ``yy``
and ``zz`` are the axes in Mm, ``zz`` increasing outward, and ``height`` is ``zz`` above
the mean height of tau_500 = 1.

The plot methods draw into ``ax``, by default the current axes, and label the axes they
draw into while those have no label.
"""
import numpy as np
from .const import BK, CC, UU
from .auxfuncs import (atomic_weights, crop_inputs, interp_tau,
                       near_equal_split)
import hopacity
from hopacity import optical_depth
import pandas as pd
from .plotting_scripts import scatter_density
import pathlib

#: symbol and unit of each variable, for the axis labels
LABELS = {'temp': ('T', 'K'), 'rho': (r'$\rho$', r'g cm$^{-3}$'), 'ne': (r'$n_e$', r'cm$^{-3}$'),
          'pe': (r'$P_e$', r'dyn cm$^{-2}$'), 'vx': (r'$v_x$', r'cm s$^{-1}$'),
          'vy': (r'$v_y$', r'cm s$^{-1}$'), 'vz': (r'$v_z$', r'cm s$^{-1}$'),
          'vmic': (r'$v_{\rm mic}$', r'km s$^{-1}$'), 'chi': (r'$\chi_{500}$', r'cm$^{-1}$'),
          'height': ('height', 'Mm')}
LTAU_LABEL = r'$\log_{10}\tau_{500}$'


def read_atmosphere(path, atmos_format='unknown', feh=None):
    """The model-atmosphere file ``path`` (a 3D snapshot or a 1D model) as an
    :class:`AtmosClass`, the kind of object a run's ``atmos`` is, without a run
    (:meth:`AtmosClass.from_file`)."""
    return AtmosClass.from_file(path, atmos_format, feh)


def _axes(ax):
    if ax is None:
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        ax = plt.gca()
    return ax


def _label(ax, xlabel=None, ylabel=None):
    """Label the axes of ``ax`` that have no label yet."""
    if xlabel and not ax.get_xlabel():
        ax.set_xlabel(xlabel)
    if ylabel and not ax.get_ylabel():
        ax.set_ylabel(ylabel)


def _var_label(var, suffix=''):
    """'T [K]' for 'temp', with ``suffix`` after the symbol; None for an array."""
    if not isinstance(var, str) or var not in LABELS:
        return None
    symbol, unit = LABELS[var]
    return f'{symbol}{suffix} [{unit}]'


class AtmosClass(object):
    @property
    def pe(self):
        if hasattr(self, '_pe'):
            return self._pe
        return self.ne * BK * self.temp

    @pe.setter
    def pe(self, value):
        self._pe = value

    #=====================================================================
    @classmethod
    def from_df(cls, df):
        self = cls()
        for key in ['zz', 'pe', 'ne', 'rho', 'temp', 'vx', 'vy', 'vz', 'chi', 'vmic', 'tau', 'ltau']:
            if key in df.columns:
                setattr(self, key, np.array(df[key]))

        return self

    #=====================================================================
    @classmethod
    def from_run(cls, run):
        """The atmosphere ``run`` saved on its mesh."""
        self = cls()
        self.run = run
        atmos_params = run.nml['atmos_params']
        saved_atmos = 'save' in atmos_params['atmos_format'].lower()

        atmos_folder = pathlib.Path(atmos_params['atmos_file']) if saved_atmos else run.folder
        sfolder = atmos_folder / 'save'

        dims = atmos_params['dims']

        with open(atmos_folder / 'atmos_mesh.txt', 'r') as mesh:
            if atmos_params['amr']:
                self.dz = []
                nnz = []
                for _ in range(dims):
                    nx, ny, nz = [int(i) for i in mesh.readline().split()]
                    dx, dy, dz = [float(i) for i in mesh.readline().split()]
                    self.dz = np.concatenate([self.dz, np.ones(nz)*dz])
                    nnz.append(nz)
            else:
                nx, ny, nz = [int(i) for i in mesh.readline().split()]
                dx, dy, dz = [float(i) for i in mesh.readline().split()]
                self.dz = np.ones(nz)*dz
                nnz = near_equal_split(nz, dims)

            self.dx = dx
            self.dy = dy
            self.xx = np.arange(nx)*dx
            self.yy = np.arange(ny)*dy
            self.zz = np.cumsum(self.dz)
            self.atmos_cuts = np.cumsum(nnz)

        self.dim = '3D' if nx*ny > 1 else '1D'

        if (sfolder / 'atmos_meta.txt').exists():
            mem = self.run.read_patch_save('atmos', folder=sfolder, zdim=2)
            mem = mem[..., 0]
            if self.dim == '1D': mem = mem.squeeze()
            self.ne   = np.array(mem[...,0])
            self.rho  = np.array(mem[...,1])
            self.temp = np.array(mem[...,2])
            self.vx   = np.array(mem[...,3]) * CC
            self.vy   = np.array(mem[...,4]) * CC
            self.vz   = np.array(mem[...,5]) * CC
            self.chi  = np.array(mem[...,6])
            self.vmic = np.array(mem[...,7])
            self.tau  = optical_depth(self.chi, np.log(self.rho), self.zz*1e8)
            self.ltau = np.log10(self.tau)
            self._set_height()

        return self

    #=====================================================================
    @classmethod
    def from_file(cls, path, atmos_format='unknown', feh=None):
        """A model-atmosphere file at its own resolution, read as dispatch.x reads it
        (:func:`hopacity.read_atmosphere`), without a run.

        ``atmos_format`` takes the &atmos_params names; 'unknown' infers it from the
        name. tau_500 is integrated from HOP's continuum opacity at the composition
        ``feh`` (:mod:`tumult.avg3d`), by default the file's own [Fe/H]; a file that
        does not state it (Stagger2) needs ``feh``. The ghost zones are dropped, and the
        axes start at zero.
        """
        from . import avg3d
        snap = hopacity.read_atmosphere(str(path), atmos_format)
        feh = snap.feh if feh is None else feh
        if feh is None:
            raise ValueError(f'{path} does not state its [Fe/H]: pass feh')
        table = avg3d.OpacityTable.for_snapshots([snap], feh=feh)
        cols = avg3d.columns(snap, table)
        nx, ny, nz = snap.shape
        self = cls()
        self.dim = '3D' if nx*ny > 1 else '1D'

        def on_mesh(a):     # (columns, depth) top first -> the mesh, deepest first
            a = np.asarray(a).reshape(nx, ny, -1)[..., ::-1]
            return a.squeeze() if self.dim == '1D' else a

        self.temp, self.rho, self.pe = on_mesh(cols.temp), on_mesh(cols.rho), on_mesh(cols.pe)
        self.ne = self.pe / (BK * self.temp)
        self.vx, self.vy, self.vz = on_mesh(cols.vx), on_mesh(cols.vy), on_mesh(cols.vz)
        self.ltau = on_mesh(cols.ltau500)
        self.tau = 10**self.ltau
        self.chi = table.chi500(np.log(self.temp), np.log(self.rho))
        inner = slice(snap.gz, nz - snap.gz)
        vmic = snap.vmic[..., inner].astype(float)
        lrho = snap.lrho[..., inner]
        if lrho[..., 0].mean() > lrho[..., -1].mean():          # the file has its top last
            vmic = vmic[..., ::-1]
        self.vmic = on_mesh(vmic)
        z = np.asarray(cols.z[0])[::-1] * 1e-8                  # height [Mm], deepest first
        self.zz = z - z[0]
        self.dz = np.gradient(self.zz)
        self.xx = (snap.xx - snap.xx[0]).astype(float) * 1e-8
        self.yy = (snap.yy - snap.yy[0]).astype(float) * 1e-8
        self.dx = float(self.xx[1]) if nx > 1 else 0.0
        self.dy = float(self.yy[1]) if ny > 1 else 0.0
        self.atmos_cuts = np.array([], int)
        self._set_height()
        return self

    def _set_height(self):
        """``z0``, the mean height of the tau_500 = 1 surface, and ``height`` above it."""
        shape = (len(self.xx), len(self.yy), 1) if self.dim == '3D' else (1,)
        zz = np.tile(self.zz, shape)
        self.z0 = np.mean(self.get_ltau_surface(0, zz))
        self.height = self.zz - self.z0

    #=====================================================================
    def get_toth(self):
        abundances = self.run.abundances
        xmy = np.array([10**(abundances[key]-12)*atomic_weights[key]
                        for key in abundances.keys()])
        grph = np.sum(xmy) * UU
        toth = self.rho / grph
        return toth

    #=====================================================================
    def get_variable(self, var):
        if type(var) is str:
            var = getattr(self, var)
        return var

    #=====================================================================
    def plot_slice(self, var, ix=None, iy=None, iz=None, shading='nearest',
                   cmap='hot', rasterized=True, logNorm=False, contrast=False,
                   show_cuts=False, title=None, tau_contours=False, ax=None):
        """A cut through the cube: vertical at column ``ix`` or row ``iy``, horizontal
        at depth index ``iz``, heights relative to the mean height of tau_500 = 1."""
        import matplotlib
        from matplotlib.colors import LogNorm
        ax = _axes(ax)
        var = self.get_variable(var)

        if self.dim == '3D':
            ltau = np.mean(self.ltau, axis=(0,1))
        else:
            ltau = self.ltau

        step = -1 if ltau[0] > ltau[-1] else 1
        z0 = np.interp(0, ltau[::step], self.zz[::step])

        xx = self.xx
        yy = self.yy
        zz = self.zz - z0

        if ix is not None:
            slc = (ix, slice(None), slice(None))
            X, Y = np.meshgrid(yy, zz)
            tle = f'X = {xx[ix]:.1f} Mm'
            _label(ax, 'Y [Mm]', 'Z [Mm]')
        elif iy is not None:
            slc = (slice(None), iy, slice(None))
            X, Y = np.meshgrid(xx, zz)
            tle = f'Y = {yy[iy]:.1f} Mm'
            _label(ax, 'X [Mm]', 'Z [Mm]')
        elif iz is not None:
            slc = (slice(None), slice(None), iz)
            X, Y = np.meshgrid(xx, yy)
            tle = f'Z = {zz[iz]:.1f} Mm'
            _label(ax, 'X [Mm]', 'Y [Mm]')

        Z = var[slc]

        if iz is None and contrast:
            Z = Z / np.mean(var, axis=(0,1))

        norm = [None, LogNorm()][logNorm]
        title = [title, tle][title is None]

        pmsh = ax.pcolormesh(X, Y, Z.T, shading=shading, cmap=cmap, rasterized=rasterized, norm=norm)
        ax.set_title(title)

        if iz is None:
            if show_cuts:
                [ax.axhline(zz[i-1], lw=0.5, color='k') for i in self.atmos_cuts]

            if tau_contours:
                ax2 = ax.twinx()
                Z = self.ltau[slc]
                levels = np.arange(int(np.min(Z)), int(np.max(Z))+1)
                cm = matplotlib.colormaps[cmap]
                cmap2 = ['binary', 'binary_r'][cm(1) < cm(0)]
                ax2.contour(X, Y, Z.T, levels=levels, cmap=cmap2, linewidths=1, linestyles=':')
                ax2.set_ylabel(LTAU_LABEL)
                zlevels = np.interp(levels, Z[-1,:][::step], zz[::step])
                ax2.set_yticks(zlevels)
                ax2.set_yticklabels(levels)

        return pmsh

    #=====================================================================
    def plot_ltau_slice(self, var, ltau=0, cmap='hot', shading='nearest', rasterized=True,
                        colorbar=True, ax=None, **kwargs):
        """Map of ``var`` on the surface log10 tau_500 = ``ltau`` of a 3D atmosphere,
        with a labelled colour bar unless ``colorbar`` is False. Returns the mesh."""
        if self.dim != '3D':
            raise ValueError('plot_ltau_slice maps a 3D atmosphere')
        ax = _axes(ax)
        surf = self.get_ltau_surface(ltau, var)[..., 0]
        mesh = ax.pcolormesh(self.xx, self.yy, surf.T, shading=shading, cmap=cmap,
                             rasterized=rasterized, **kwargs)
        ax.set_aspect('equal')
        _label(ax, 'X [Mm]', 'Y [Mm]')
        if colorbar:
            label = _var_label(var, rf' at $\log_{{10}}\tau_{{500}}={ltau:g}$')
            ax.figure.colorbar(mesh, ax=ax, label=label)
        return mesh

    #=====================================================================
    def interp_surf(self, xvar, yvar, newx=None, log=False, power=None, zaxis=-1):
        xvar = self.get_variable(xvar)
        yvar = self.get_variable(yvar)
        if zaxis != -1: yvar = np.moveaxis(yvar, zaxis, -1)
        step = -1 if xvar[...,0].ravel()[0] > xvar[...,-1].ravel()[0] else 1
        xx = xvar[...,::step].copy()
        yy = yvar[...,::step].copy()
        if log: yy = np.log(yy)
        if power: yy = yy**power
        if newx is None:
            axis = tuple(i for i in range(xx.ndim-1))
            newx = np.mean(xvar, axis=axis)
        newy = interp_tau(newx, xx, yy)
        if power: newy = newy**(1/power)
        if log: newy = np.exp(newy)
        if zaxis != -1: newy = np.moveaxis(newy, -1, zaxis)
        return newx, newy

    #=====================================================================
    def plot_geo_cube(self, var, marker=',', linestyle='none', rasterized=True,
                      xmin=None, xmax=None, ymin=None, ymax=None, ax=None, **kwargs):
        """Every point of ``var`` against the height."""
        ax = _axes(ax)
        label = _var_label(var)
        var = self.get_variable(var)
        xx = (self.height * np.ones_like(self.ltau)).ravel()
        yy = var.ravel()
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        ax.set_xlim(np.nanmin(xx), np.nanmax(xx))
        ax.set_ylim(np.nanmin(yy), np.nanmax(yy))
        _label(ax, _var_label('height'), label)
        return ax.plot(xx, yy, marker=marker, linestyle=linestyle,
                       rasterized=rasterized, **kwargs)

    #=====================================================================
    def plot_ltau_cube(self, var, marker=',', linestyle='none', rasterized=True,
                       xmin=None, xmax=None, ymin=None, ymax=None, density=False,
                       cmap='GnBu', ax=None, **kwargs):
        """Every point of ``var`` against log10 tau_500, or with ``density`` the
        logarithm of their number as an image (:func:`scatter_density`)."""
        ax = _axes(ax)
        label = _var_label(var)
        var = self.get_variable(var)
        xx = self.ltau.ravel()
        yy = var.ravel()
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        ax.set_xlim(np.nanmin(xx), np.nanmax(xx))
        ax.set_ylim(np.nanmin(yy), np.nanmax(yy))
        _label(ax, LTAU_LABEL, label)

        if density:
            return scatter_density(xx, yy, cmap=cmap, xmin=xmin, xmax=xmax, ymin=ymin,
                                   ymax=ymax, ax=ax, **kwargs)
        return ax.plot(xx, yy, marker=marker, linestyle=linestyle,
                       rasterized=rasterized, **kwargs)

    #=====================================================================
    def get_surf_mean(self, xvar, yvar, newx=None, log=False, power=None, axis=-1, **kwargs):
        xvar = self.get_variable(xvar)
        yvar = self.get_variable(yvar)

        if self.dim == '3D':
            newx, yy3d = self.interp_surf(xvar, yvar, newx=newx, log=log, power=power, zaxis=axis, **kwargs)
            axis = tuple(np.delete(np.arange(yy3d.ndim), axis))
            newy = advanced_mean(yy3d, axis=axis, log=log, power=power)
        else:
            newx = newx or xvar
            step = -1 if xvar[0] > xvar[-1] else 1
            newy = np.interp(newx, xvar[::step], yvar[::step])

        return newx, newy

    #=====================================================================
    def plot_ltau_mean(self, var, log=False, newx=None, xmin=None, xmax=None,
                       ymin=None, ymax=None, power=None, ax=None, **kwargs):
        """The mean of ``var`` on surfaces of constant tau_500 (geometric with ``log``,
        of ``var**power`` with ``power``) against log10 tau_500."""
        ax = _axes(ax)
        xx, yy = self.get_surf_mean('ltau', var, log=log, newx=newx, power=power)
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        _label(ax, LTAU_LABEL, _var_label(var))
        return ax.plot(xx, yy, **kwargs)

    #=====================================================================
    def get_geo_mean(self, var, newx=None, log=False, power=None):
        if hasattr(self, 'runs'):
            xx = np.mean([run.height for run in self.runs], axis=0)
            yy = [run.get_geo_mean(var, newx=xx, log=log, power=power) for run in self.runs]
            yy = advanced_mean(yy, axis=0, log=log, power=power)
            return xx, yy
        else:
            yy = self.get_variable(var)

        xx = newx if newx else self.height

        if self.dim == '3D':
            yy = advanced_mean(yy, axis=(0,1), log=log, power=power)

        if xx is not self.height:
            yy = np.interp(xx, self.height, yy)

        return xx, yy

    #=====================================================================
    def plot_geo_mean(self, var, log=False, xmin=None, xmax=None,
                      ymin=None, ymax=None, ax=None, **kwargs):
        """The horizontal mean of ``var`` against the height."""
        ax = _axes(ax)
        xx, yy = self.get_geo_mean(var, log=log)
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        _label(ax, _var_label('height'), _var_label(var))
        return ax.plot(xx, yy, **kwargs)

    #=====================================================================
    def get_ltau_surface(self, ltau, var, log=False, zaxis=-1):
        if hasattr(self, 'runs'):
            yy = np.array([run.get_ltau_surface(ltau, var, log=log, zaxis=zaxis) for run in self.runs])
        else:
            var = self.get_variable(var)

            if self.dim == '3D':
                newx = np.array([ltau])
                _, yy = self.interp_surf(self.ltau, var, newx=newx, log=log, zaxis=zaxis)
            else:
                ltau1d = self.ltau
                step = -1 if ltau1d[0] > ltau1d[-1] else 1
                yy = np.interp(ltau, ltau1d[::step], var[::step])

        return yy

    #=====================================================================
    def plot_ltau_hist(self, ltau, var, log=False, power=None, zaxis=-1, show_mean=False,
                       color=None, ax=None, **kwargs):
        """Histogram of ``var`` over the columns on the surface log10 tau_500 = ``ltau``,
        with its mean dashed if ``show_mean``."""
        ax = _axes(ax)
        yy = self.get_ltau_surface(ltau, var, log=log, zaxis=zaxis)

        if show_mean:
            mean = advanced_mean(yy, log=log, power=power)
            ax.axvline(mean, ls='--', color=color)

        _label(ax, _var_label(var, rf' at $\log_{{10}}\tau_{{500}}={ltau:g}$'), 'Columns')
        return ax.hist(yy.ravel(), color=color, **kwargs)

    #=====================================================================
    def plot_cube_slice(self, cube, ix=None, iy=None, iz=None,
                        ang=None, rasterized=True, shading='nearest',
                        show_labels=True, cmap='hot', rot90=0, ax=None, **kwargs):
        """A cut through ``cube`` on this mesh (a 2D ``cube`` is a horizontal map),
        foreshortened for the ray ``ang`` of the run."""
        ax = _axes(ax)
        if iz is not None:
            xx = self.xx + self.dx/2
            yy = self.yy + self.dy/2
            Z = cube[:,:,iz].T
            xlabel = 'X [Mm]'
            ylabel = 'Y [Mm]'
        elif iy is not None:
            xx = self.xx + self.dx/2
            yy = self.zz - self.z0
            Z = cube[:,iy,:].T
            xlabel = 'X [Mm]'
            ylabel = 'Z [Mm]'
        elif ix is not None:
            xx = self.yy + self.dy/2
            yy = self.zz - self.z0
            Z = cube[ix,:,:].T
            xlabel = 'Y [Mm]'
            ylabel = 'Z [Mm]'
        elif cube.ndim == 2:
            xx = self.xx + self.dx/2
            yy = self.yy + self.dy/2
            Z = cube.T
            xlabel = 'X [Mm]'
            ylabel = 'Y [Mm]'

        if ang is None:
            mux, muy, muz = (0.0, 0.0, 1.0)
        else:
            mux, muy, muz = self.run.vec[ang]

        for _ in range(rot90):
            tmp = -xx
            xx = yy
            yy = tmp

            tmp = -mux
            mux = muy
            muy = tmp
            Z = np.rot90(Z, k=1)

        X, Y = np.meshgrid(xx, yy)
        ax.pcolormesh(X, Y, Z, shading=shading, rasterized=rasterized, cmap=cmap, **kwargs)

        xfac = np.sqrt( (1 - mux**2)**2 + (mux*muy)**2 + (mux*muz)**2 )
        yfac = np.sqrt( (1 - muy**2)**2 + (mux*muy)**2 + (muy*muz)**2 )
        ax.set_aspect(yfac/xfac)

        if show_labels:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if ang is not None:
                ax.set_title(f'({mux:.2f},{muy:.2f},{muz:.2f})')
        else:
            ax.set_xticks([])
            ax.set_yticks([])

    #=====================================================================
    def to_dispatch1d_atmos(self, file, xvar='ltau', vmic=None, FeH=''):
        _, ltau = self.get_surf_mean(xvar, 'ltau', log=False)
        _, temp = self.get_surf_mean(xvar, 'temp', log=True)
        _, ne   = self.get_surf_mean(xvar, 'ne',   log=True)
        _, rho  = self.get_surf_mean(xvar, 'rho',  log=True)

        if vmic is None:
            vx = self.get_surf_mean(xvar, 'vx')[1]**2
            vy = self.get_surf_mean(xvar, 'vy')[1]**2
            vz = self.get_surf_mean(xvar, 'vz')[1]**2
            vx2 = self.get_surf_mean(xvar, 'vx', power=2)[1]**2
            vy2 = self.get_surf_mean(xvar, 'vy', power=2)[1]**2
            vz2 = self.get_surf_mean(xvar, 'vz', power=2)[1]**2
            vmic = 1e-5 * np.sqrt(((vx2 - vx) + (vy2 - vy) + (vz2 - vz)) / 3)

        df = pd.DataFrame()
        df['lTau'] = ltau
        df['Temp'] = temp
        df['Pe'] = ne*temp*BK
        df['Rho'] = rho
        df['Vturb'] = vmic

        file = pathlib.Path(file)
        file.parent.mkdir(parents=True, exist_ok=True)
        f = open(file, 'w+')
        f.write(f'{file.with_suffix("").name:}\n')
        f.write(f'{len(df):} {FeH}\n')
        text = df[::-1].to_string(index=False, float_format='%.4e')
        text = '*' + text[1:]
        f.write(text)
        f.close()

      #=================================================================
    def to_multi1d_atmos(self, folder, name, logg, abund=None, FeH=0, vz=0.0):
        df = pd.DataFrame()
        df['ltau'] = self.ltau
        df['Temp'] = self.temp
        df['Ne'] = self.ne
        df['vz'] = vz
        df['vturb'] = self.vmic
        df = df[::-1].reset_index(drop=True)
        folder = pathlib.Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        file = folder / f'atmos.{name}'
        f = open(file, 'w+')
        f.write(f'{name}\n')
        f.write('TAU\n')
        f.write('* LG G\n')
        f.write(f' {logg:.2f}\n')
        f.write('* NDEP\n')
        f.write(f' {len(self.ltau):d}\n')
        text = df.to_string(index=False, float_format='%.4e')
        text = '*' + text[1:]
        f.write(text)
        f.close()
        file = folder / f'dscale.{name}'
        f = open(file, 'w+')
        f.write(f'{name}\n')
        f.write('TAU\n')
        f.write(f' {len(df.ltau):d} {df.ltau[0]:.4e}\n')
        f.write(df.to_string(index=False, columns=['ltau'], float_format='%.4e', header=None))
        f.close()
        if abund is not None:
            file = folder / f'abund.{name}'
            f = open(file, 'w+')
            for eid, abnd in abund.items():
                f.write(f'{eid.upper():2} {abnd+FeH:8.4f}\n')
            f.close()
  
        #=================================================================
    @classmethod
    def combine(cls, atmos_list, weights=None):
      #-----------------------------------------------------------------
        self = cls()
        if weights is not None:
                taumin = np.min([a.ltau for a in atmos_list])
                taumax = np.max([a.ltau for a in atmos_list])
                ltau = np.linspace(taumin, taumax, 128)

                def average_atmos(var):
                        return np.sum(self.get_ltau_surface(ltau, var).T * weights, axis=-1)

                self.ne = average_atmos('ne')
                self.rho = average_atmos('rho')
                self.temp = average_atmos('temp')
                self.vx = average_atmos('vx')
                self.vy = average_atmos('vy')
                self.vz = average_atmos('vz')
                self.chi = average_atmos('chi')
                self.tau = 10**ltau
                self.ltau = ltau

        else:
                self.ne = np.array([a.ne for a in atmos_list])
                self.rho = np.array([a.rho for a in atmos_list])
                self.temp = np.array([a.temp for a in atmos_list])
                self.vx = np.array([a.vx for a in atmos_list])
                self.vy = np.array([a.vy for a in atmos_list])
                self.vz = np.array([a.vz for a in atmos_list])
                self.chi = np.array([a.chi for a in atmos_list])
                self.tau = np.array([a.tau for a in atmos_list])
                self.ltau = np.array([a.ltau for a in atmos_list])

        return self
  
def advanced_mean(values, log=False, power=None, **kwargs):
      if log: values = np.log(values)
      if power: values = values**power
      mean_values = np.mean(values, **kwargs)
      if power: mean_values = mean_values**(1/power)
      if log: mean_values = np.exp(mean_values)
      return mean_values