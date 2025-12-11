import numpy as np
from .auxfuncs import *
from .plotting_scripts import *
import pathlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

class AtmosClass(object):
  @property
  def pe(self):
    if hasattr(self, '_pe'):
        return self._pe
    else:
        return self.ne * BK * self.temp
  
  @pe.setter
  def pe(self, value):
    self._pe = value

  #=====================================================================
  @classmethod
  def from_run(cls, run):
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
        if self.dim=='1D': mem = mem.squeeze()
        self.ne   = np.array(mem[...,0])
        self.rho  = np.array(mem[...,1])
        self.temp = np.array(mem[...,2])
        self.vx   = np.array(mem[...,3]) * CC
        self.vy   = np.array(mem[...,4]) * CC
        self.vz   = np.array(mem[...,5]) * CC
        self.chi  = np.array(mem[...,6])
        self.vmic = np.array(mem[...,7])
        self.tau  = compute_tau_scale(self.chi, self.zz*1e8, axis=-1)
        self.ltau = np.log10(self.tau)

        zz = np.tile(self.zz, (nx, ny, 1)).squeeze()
        self.z0 = np.mean(self.get_ltau_surface(0, zz))
        self.height = self.zz - self.z0

    return self

  #=====================================================================
  def get_toth(self):
    #-----------------------------------------------------------------
    abundances = self.run.abundances
    xmy = np.array([10**(abundances[key]-12)*atomic_weights[key]
                    for key in abundances.keys()])
    grph = np.sum(xmy) * UU
    toth = self.rho / grph
    return toth
  
  #=====================================================================
  def get_variable(self, var):
    #-----------------------------------------------------------------
    if type(var) is str:
        var = getattr(self, var)
    return var
  
  #=================================================================    
  def plot_slice(self, var, ix=None, iy=None, iz=None, shading='nearest',
              cmap='hot', rasterized=True, logNorm=False, contrast=False,
              show_cuts=False, title=None, tau_contours=False):
    #-----------------------------------------------------------------
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
        plt.xlabel('Y [Mm]')
        plt.ylabel('Z [Mm]')
    elif iy is not None:
        slc = (slice(None), iy, slice(None))
        X, Y = np.meshgrid(xx, zz)
        tle = f'Y = {yy[iy]:.1f} Mm'
        plt.xlabel('X [Mm]')
        plt.ylabel('Z [Mm]')
    elif iz is not None:
        slc = (slice(None), slice(None), iz)
        X, Y = np.meshgrid(xx, yy)
        tle = f'Z = {zz[iz]:.1f} Mm'
        plt.xlabel('X [Mm]')
        plt.ylabel('Y [Mm]')

    Z = var[slc]

    if iz is None and contrast:
        Z = Z / np.mean(var, axis=(0,1))

    norm = [None, LogNorm()][logNorm]
    title = [title, tle] [title is None]

    pmsh = plt.pcolormesh(X, Y, Z.T, shading=shading, cmap=cmap, rasterized=rasterized, norm=norm)
    plt.title(title)

    if iz is None:
        if show_cuts:
            [plt.axhline(zz[i-1], lw=0.5, color='k') for i in self.atmos_cuts]

        if tau_contours:
            ax = plt.gca()
            ax2 = ax.twinx()
            Z = self.ltau[slc]
            levels = np.arange(int(np.min(Z)), int(np.max(Z))+1)
            cm = plt.cm.__dict__[cmap]
            cmap2 = ['binary', 'binary_r'][cm(1) < cm(0)]
            ax2.contour(X, Y, Z.T, levels=levels, cmap=cmap2, linewidths=1, linestyles=':')
            ax2.set_ylabel(r'$\log \tau_{500}$')
            zlevels = np.interp(levels, Z[-1,:][::step], zz[::step])
            ax2.set_yticks(zlevels)
            ax2.set_yticklabels(levels)

    return pmsh

  #=====================================================================
  def interp_ltau_surf(self, var, newx=None, log=False, zaxis=-1):
    #-----------------------------------------------------------------
    var = self.get_variable(var)
    if zaxis != -1: var = np.moveaxis(var, zaxis, -1)
    ltau = self.ltau
    step = -1 if ltau[...,0].ravel()[0] > ltau[...,-1].ravel()[0] else 1
    xx = ltau[...,::step].copy()
    yy = var[...,::step].copy()
    if log: yy = np.log(yy)
    if newx is None:
        axis = tuple(i for i in range(xx.ndim-1))
        newx = np.mean(xx, axis=axis)
    newy = interp_tau(newx, xx, yy)
    if log: newy = np.exp(newy)
    if zaxis != -1: newy = np.moveaxis(newy, -1, zaxis)
    return newx, newy
  
  #=====================================================================
  def plot_geo_cube(self, var, marker=',', linestyle='none', rasterized=True,
                      xmin=None, xmax=None, ymin=None, ymax=None, **kwargs):
    #-----------------------------------------------------------------
    var = self.get_variable(var)
    xx = (self.height * np.ones_like(self.ltau)).ravel()
    yy = var.ravel()
    if xmin is None: xmin = np.nanmin(xx)
    if xmax is None: xmax = np.nanmax(xx)
    if ymin is None: ymin = np.nanmin(yy)
    if ymax is None: ymax = np.nanmax(yy)
    mask = (xx > xmin) * (xx < xmax) * (yy > ymin) * (yy < ymax)
    xx = xx[mask]
    yy = yy[mask]
    plt.xlim(np.nanmin(xx), np.nanmax(xx))
    plt.ylim(np.nanmin(yy), np.nanmax(yy))
    return plt.plot(xx, yy, marker=marker, linestyle=linestyle,
                    rasterized=rasterized, **kwargs)
  
  #=====================================================================
  def plot_ltau_cube(self, var, marker=',', linestyle='none', rasterized=True,
                      xmin=None, xmax=None, ymin=None, ymax=None, density=False,
                      cmap='GnBu', **kwargs):
    #-----------------------------------------------------------------
    var = self.get_variable(var)
    xx = self.ltau.ravel()
    yy = var.ravel()
    xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
    plt.xlim(np.nanmin(xx), np.nanmax(xx))
    plt.ylim(np.nanmin(yy), np.nanmax(yy))

    if density:
        return scatter_density(xx, yy, cmap=cmap, **kwargs)
    else:
        return plt.plot(xx, yy, marker=marker, linestyle=linestyle,
                        rasterized=rasterized, **kwargs)

  #=====================================================================
  def get_ltau_mean(self, var, log=False, newx=None, **kwargs):
    #-----------------------------------------------------------------
    var = self.get_variable(var)
    if self.dim == '3D':
        newx, yy3d = self.interp_ltau_surf(var, log=log, newx=newx, **kwargs)
        if log: yy3d = np.log(yy3d)
        axis = tuple(i for i in range(yy3d.ndim-1))
        yy1d = np.mean(yy3d, axis=axis)
        if log: yy1d = np.exp(yy1d)
    else:
        if newx is None: newx = self.ltau[::-1]
        step = -1 if self.ltau[0] > self.ltau[-1] else 1
        yy1d = np.interp(newx, self.ltau[::step], var[::step])

    return newx, yy1d
  
  #=====================================================================
  def plot_ltau_mean(self, var, log=False, newx=None, xmin=None, xmax=None,
                      ymin=None, ymax=None, **kwargs):
    #-----------------------------------------------------------------
    xx, yy = self.get_ltau_mean(var, log=log, newx=newx)
    xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
    return plt.plot(xx, yy, **kwargs)
  
  #=====================================================================
  def get_geo_mean(self, var, newx=None, log=False):
    if hasattr(self, 'runs'):
        xx = np.mean([run.height for run in self.runs], axis=0)
        yy = np.mean([run.get_geo_mean(var, newx=xx) for run in self.runs], axis=0)
        return xx, yy
    else:
      yy = self.get_variable(var)

    xx = newx if newx else self.height

    if self.dim == '3D':
        if log: yy = np.log(yy)
        yy = np.mean(yy, axis=(0,1))
        if log: yy = np.exp(yy)
        
    if xx is not self.height:
        yy = np.interp(xx, self.height, yy)

    return xx, yy
  
  #=====================================================================
  def plot_geo_mean(self, var, log=False, xmin=None, xmax=None,
                      ymin=None, ymax=None, **kwargs):
    #-----------------------------------------------------------------
    xx, yy = self.get_geo_mean(var, log=log)
    xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
    return plt.plot(xx, yy, **kwargs)
  
  #=================================================================
  def get_ltau_surface(self, ltau, var, log=False, zaxis=-1):
    #-----------------------------------------------------------------
    if hasattr(self, 'runs'):
        yy = np.array([run.get_ltau_surface(ltau, var, log=log, zaxis=zaxis) for run in self.runs])
    else:
        var = self.get_variable(var)

        if self.dim == '3D':
            newx = np.array([ltau])
            _, yy = self.interp_ltau_surf(var, newx=newx, log=log, zaxis=zaxis)
        else:
            ltau1d = self.ltau
            step = -1 if ltau1d[0] > ltau1d[-1] else 1
            yy = np.interp(ltau, ltau1d[::step], var[::step])
            
    return yy

  #=================================================================
  def plot_ltau_hist(self, ltau, var, log=False, zaxis=-1, show_mean=False, color=None, **kwargs):
    #-----------------------------------------------------------------
    yy = self.get_ltau_surface(ltau, var, log=log, zaxis=zaxis)

    if show_mean:
        if log:
            mean = 10**np.mean(np.log10(yy))
        else:
            mean = np.mean(yy)
            
        plt.axvline(mean, ls='--', color=color)

    return plt.hist(yy.ravel(), color=color, **kwargs)
  
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