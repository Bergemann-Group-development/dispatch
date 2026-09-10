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
    def from_df(cls, df):
        self = cls()
        for key in ['zz', 'pe', 'ne', 'rho', 'temp', 'vx', 'vy', 'vz', 'chi', 'vmic', 'tau', 'ltau']:
            if key in df.columns:
                    setattr(self, key, np.array(df[key]))

        return self

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
    def interp_surf(self, xvar, yvar, newx=None, log=False, power=None, zaxis=-1):
      #-----------------------------------------------------------------
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
    def get_surf_mean(self, xvar, yvar, newx=None, log=False, power=None, axis=-1, **kwargs):
      #-----------------------------------------------------------------
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
                                            ymin=None, ymax=None, power=None, **kwargs):
                          #-----------------------------------------------------------------
        xx, yy = self.get_surf_mean('ltau', var, log=log, newx=newx, power=power)
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        return plt.plot(xx, yy, **kwargs)
  
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
                        _, yy = self.interp_surf(self.ltau, var, newx=newx, log=log, zaxis=zaxis)
                else:
                        ltau1d = self.ltau
                        step = -1 if ltau1d[0] > ltau1d[-1] else 1
                        yy = np.interp(ltau, ltau1d[::step], var[::step])
            
        return yy

      #=================================================================
    def plot_ltau_hist(self, ltau, var, log=False, power=None, zaxis=-1, show_mean=False, color=None, **kwargs):
      #-----------------------------------------------------------------
        yy = self.get_ltau_surface(ltau, var, log=log, zaxis=zaxis)

        if show_mean:
                mean = advanced_mean(yy, log=log, power=power)
                plt.axvline(mean, ls='--', color=color)

        return plt.hist(yy.ravel(), color=color, **kwargs)
  
      #=====================================================================
    def plot_cube_slice(self, cube, ix=None, iy=None, iz=None,
                                            ang=None, rasterized=True, shading='nearest', 
                                            show_labels=True, cmap='hot', rot90=0, **kwargs):
                            #-----------------------------------------------------------------
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

            for i in range(rot90):
                    tmp = -xx
                    xx = yy
                    yy = tmp

                    tmp = -mux
                    mux = muy
                    muy = tmp
                    Z = np.rot90(Z, k=1)

            X, Y = np.meshgrid(xx, yy)
            plt.pcolormesh(X, Y, Z, shading=shading, rasterized=rasterized, cmap=cmap, **kwargs)
            ax = plt.gca()

            xfac = np.sqrt( (1 - mux**2)**2 + (mux*muy)**2 + (mux*muz)**2 )
            yfac = np.sqrt( (1 - muy**2)**2 + (mux*muy)**2 + (muy*muz)**2 )
            ax.set_aspect(yfac/xfac)

            if show_labels:
                    plt.xlabel(xlabel)
                    plt.ylabel(ylabel)
                    if ang is not None:
                            plt.title(f'({mux:.2f},{muy:.2f},{muz:.2f})')
            else:
                    plt.xticks([])
                    plt.yticks([])

            #=================================================================
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