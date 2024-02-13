from .atom import atom_class
from .trans import trans_class
from .auxfuncs import vacuum2obs, continuum, freq_quad, convol, SBOL, CC, interp_tau, check_file, compute_tau_scale
import numpy as np
import pandas as pd
import dask.array as da
import matplotlib.pyplot as plt
import glob
import f90nml
from scipy.signal import argrelextrema
from scipy.io import FortranFile
from matplotlib.collections import LineCollection
from matplotlib.ticker import MultipleLocator
import os

def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=32, log_file=None):
    cwd = os.getcwd()
    check_file(cwd + '/dispatch.x')
    if default is not None:
        input = f90nml.read(default + '/params.nml')
    else:
        input = f90nml.Namelist()

    input_file = input_name + '.nml'
    if not os.path.isfile(input_file) or overwrite:
        if nml_settings is not None:
            input.patch(nml_settings)

        input.write(input_file, force=True)

        if log_file is None:
            command = './dispatch.x {:}'.format(input_file)
        else:
            command = './dispatch.x {:} > {:}'.format(input_file, log_file)

        print(input_file, ' launched')
        os.environ['OMP_NUM_THREADS'] = str(threads)
        iostat = os.system(command)
        if iostat!=0 and log_file is not None: 
            with open(log_file) as f:
                for line in f:
                    pass
                last_line = line
                print(last_line)
    else:
        print(input_file, ' already exists')
        input = f90nml.read(input_file)

    datadir = 'data/'
    if 'io_params' in input:
        if 'datadir' in input['io_params']:
            datadir = input['io_params']['datadir']

    return os.path.join(datadir, input_name)

#=====================================================================
class read(object):
    #=================================================================
    def __init__(self, folder, verbose=0, **kwargs):
        #-------------------------------------------------------------
        folder = folder + os.sep
        sfolder = folder + 'save' + os.sep
        self.sfolder = sfolder
        self.folder = folder
        self.nml = f90nml.read(folder + 'params.nml')   
        df = pd.read_csv(folder + 'abund.txt', header=None, sep='\s+', na_filter=False)
        self.abundances = {k.lower().capitalize():v for k, v in df.values}
        
        maxiter = self.nml['m3d_params']['maxiter']
        depart = self.nml['atom_params']['precomputed_depart'].strip()
        if (maxiter > 0):
            self.mode = 'NLTE'
        elif(depart==''):
            self.mode = 'LTE'
        else:
            self.mode = 'NLTE'
        
        atom_file = folder + 'atom.txt'
        if os.path.isfile(atom_file):
            self.qnorm = self.nml['atom_params']['qnorm']
            self.atom = atom_class(atom_file, qnorm=self.qnorm, verbose=verbose, **kwargs)
            loggf_cut = self.nml['atom_params']['logf_cut']
            red_cut   = self.nml['atom_params']['red_cut']
            blue_cut  = self.nml['atom_params']['blue_cut']
            in_abnd   = self.nml['atom_params']['abundance']
            
            if self.atom.id in self.abundances:
                self.atom.abnd = float(self.abundances[self.atom.id])
            elif in_abnd > -99:
                self.atom.abnd = in_abnd

            self.atom.title = 'A({}) = {:.2f}'.format(self.atom.id, self.atom.abnd)
            
            for trans in self.atom.trans:
                trans.run = self
        
            lines = self.atom.trans[:self.atom.nline]
            if 'line_mask' in self.nml.keys():
                if 'in_lines' in self.nml['line_mask'].keys():
                    in_lines = self.nml['line_mask']['in_lines']
                    if type(in_lines)==int: in_lines=[in_lines]
                    lines = [self.atom.trans[i-1] for i in in_lines]
                
            lines = [line for line in lines if np.log10(line.fosc*line.lower.g) > loggf_cut]
            lines = [line for line in lines if line.lam_red > blue_cut and line.lam_blue < red_cut]
            
            self.nline = len(lines)
            self.line = lines

            reports = glob.glob(folder + 'convergence_report_*.txt')
            if len(reports)>0:
                self.lvlerr = np.zeros(self.atom.nk)
                for report in reports:
                    df = pd.read_csv(report, sep='\s+')
                    self.lvlerr = np.maximum(self.lvlerr, df.maxconv.values)
        
        if os.path.isfile(sfolder + 'out_lam.bin'):
            if (verbose > 0): print('Reading wavelengths')
            self.lam = np.fromfile(sfolder + 'out_lam.bin', dtype=np.float64)
            self.lam = vacuum2obs(self.lam[::-1])
            self.dlam = np.diff(self.lam)

        if os.path.isfile(sfolder + 'bb_lam.bin'):
            self.bb_lam = vacuum2obs(np.fromfile(sfolder + 'bb_lam.bin', dtype=np.float64))[::-1]
        if os.path.isfile(sfolder + 'bf_lam.bin'):
            self.bf_lam = vacuum2obs(np.fromfile(sfolder + 'bf_lam.bin', dtype=np.float64))[::-1]
        if os.path.isfile(sfolder + 'all_lam.bin'):
            self.all_lam = vacuum2obs(np.fromfile(sfolder + 'all_lam.bin', dtype=np.float64))[::-1]
        
        file = folder + 'angle_quadrature.txt'
        if os.path.isfile(file):
            df = pd.read_csv(file, header=None, skiprows=1, sep='\s+', names=['ang', 'wmu', 'mux', 'muy', 'muz'])
            self.vec = df.values[:,2:]
            self.wts = df.values[:,1]
            self.n_ang = self.wts.size
            self.mus, self.imus = np.unique(1 - self.vec[:,-1], return_inverse=True)
            self.mus = 1 - self.mus
            self.nmu = self.mus.size

        if hasattr(self, 'lam'):
            if (verbose > 0): print('Reading Intensities')
            ie = self.read_patch_save('ie', fdim=0)[0]

            self.ie     = np.array(ie[:,0,0,:,0])
            self.ie_lte = np.array(ie[:,0,0,:,1])
            self.cnt    = np.array(ie[:,0,0,:,2])

            self.flux     = np.sum(self.wts * self.ie    , axis=-1)
            self.flux_lte = np.sum(self.wts * self.ie_lte, axis=-1)
            self.flux_cnt = np.sum(self.wts * self.cnt   , axis=-1)

            if self.nml['m3d_params']['save_resolved']:
                i3 = self.read_patch_save('i3', fdim=0)[0]
                self.i3     = i3[...,0]
                self.i3_lte = i3[...,1]
                self.c3     = i3[...,2]

        mesh = open(folder + 'atmos_mesh.txt', 'r')
        nx, ny, nz = [int(i) for i in mesh.readline().split()]
        dx, dy, dz = [float(i) for i in mesh.readline().split()]
        mesh.close()

        if nx*ny > 1:
            self.dim = '3D'
        else:
            self.dim = '1D'

        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.xx = np.linspace(0, (nx-1)*dx, nx)
        self.yy = np.linspace(0, (ny-1)*dy, ny)
        self.zz = np.linspace(0, (nz-1)*dz, nz)

        if os.path.isfile(sfolder + 'atmos_meta.txt'):
            mem = self.read_patch_save('atmos')[0]
            mem = mem[..., 0].squeeze()
            self.ne   = mem[...,0].compute()
            self.rho  = mem[...,1].compute()
            self.temp = mem[...,2].compute()
            self.vx   = mem[...,3].compute() * CC
            self.vy   = mem[...,4].compute() * CC
            self.vz   = mem[...,5].compute() * CC
            self.chi  = mem[...,6].compute()
            self.tau  = compute_tau_scale(self.chi, self.dz*1e8, axis=-1)
            self.ltau = np.log10(self.tau)
            
        if os.path.isfile(folder + 'atmos_data.bin'):
            dat = FortranFile(folder + 'atmos_data.bin')
            self.temp = dat.read_reals(np.float32)
            self.rho  = dat.read_reals(np.float32)
            self.px   = dat.read_reals(np.float32)
            self.py   = dat.read_reals(np.float32)
            self.pz   = dat.read_reals(np.float32)
            self.pe   = dat.read_reals(np.float32)
            self.tau  = dat.read_reals(np.float32)
            self.chi  = dat.read_reals(np.float32)
            dat.close()

            if nx*ny > 1:
                self.temp = self.temp.reshape((nx,ny,nz), order='F')
                self.rho  = self.rho .reshape((nx,ny,nz), order='F')
                self.px   = self.px  .reshape((nx,ny,nz), order='F')
                self.py   = self.py  .reshape((nx,ny,nz), order='F')
                self.pz   = self.pz  .reshape((nx,ny,nz), order='F')
                self.pe   = self.pe  .reshape((nx,ny,nz), order='F')
                self.tau  = self.tau .reshape((nx,ny,nz), order='F')
                self.chi  = self.chi .reshape((nx,ny,nz), order='F')

            self.ltau = np.log10(self.tau)

    #=====================================================================
    def get_xx(self, lam, dlam=None, xmin=None, xmax=None):
        #-----------------------------------------------------------------
        ib = None
        ir = None
        if dlam is not None:
            if isinstance(dlam, tuple):
                dl = dlam[:]

            elif not isinstance(dlam, (bool)):
                dl = (-dlam, dlam)

            ib = np.searchsorted(lam, dl[0])
            ir = np.searchsorted(lam, dl[1])

        if xmin is not None: ib = np.searchsorted(lam, xmin)
        if xmax is not None: ir = np.searchsorted(lam, xmax)

        mask = slice(ib, ir)
        return lam[mask], mask

    #=====================================================================
    def set_xaxis(self, xx, dlam=None, set_labels=True):
        #-----------------------------------------------------------------
        xmin = min(xx)
        xmax = max(xx)
        dxx = xmax - xmin
        if dlam is not None:
            plt.xlim(xmin, xmax)
            if set_labels: plt.xlabel('$\Delta\lambda$ [Å]')
        else:
            diff = 0.05*dxx
            plt.xlim(xmin-diff, xmax+diff)
            if set_labels: plt.xlabel('Wavelength [Å]')

        plt.ticklabel_format(useOffset=False)

    #=====================================================================
    def set_yaxis(self, ang, imu, norm, per_aa=False, renorm=False, set_labels=True):
        #-----------------------------------------------------------------
        if ang is not None:
            ylabel = 'Intensity µ={:.2f}'.format(self.vec[ang,-1])
            
        elif imu is not None:
            if imu < 0: imu = self.imus.max() + imu + 1
            ylabel = 'Intensity µ={:.2f}'.format(self.mus[imu])
            
        else:
            ylabel = 'Flux'

        if per_aa: norm=False

        if set_labels:
            if norm or renorm:
                ylabel = 'Normalised ' + ylabel
            else:
                units = 'erg/s/cm$^2$'

                if ylabel!='Flux':
                    units = units + '/sr'

                if per_aa:
                    units = units + '/Å'
                else:
                    units = units + '/Hz'

                ylabel = ylabel + ' [' + units + ']'

            plt.ylabel(ylabel)

    #=====================================================================
    def get_yy(self, ang=None, imu=None, norm=True, cont=False, LTE=False, mask=slice(None, None)):
        #-----------------------------------------------------------------
        if imu is not None:
            if imu >= self.nmu:
                raise ValueError('imu is out of range. Nmu = {:}'.format(self.nmu))
            elif imu < 0:
                imu = self.imus.max() + imu + 1

        if self.n_ang == 1: ang=0

        if ang is None and imu is None:
            cnt = self.flux_cnt
            if cont:
                yy = self.flux_cnt
            elif LTE:
                yy = self.flux_lte
            else:   
                yy = self.flux

            cnt = cnt[mask]
            yy  = yy[mask]

        else:
            cnt = self.cnt
            if cont:
                yy = self.cnt
            elif LTE:
                yy = self.ie_lte
            else:
                yy = self.ie

            if ang is not None:
                cnt = cnt[mask,ang]
                yy = yy[mask,ang]
            
            elif imu is not None:
                mu_mask = self.imus==imu
                cnt = np.mean(cnt[mask, mu_mask], axis=-1)
                yy  = np.mean(yy [mask, mu_mask], axis=-1)
                
        
        if norm:
            yy = yy / cnt

            
        return np.array(yy), np.array(cnt)
    
    #=====================================================================
    def crop(self, lam0=None, dlam=None, xmin=None, xmax=None, vdop=None, 
             vshift=None, convolve=None, per_aa=False, **kwargs):
        #-----------------------------------------------------------------
        xx = self.lam

        if lam0 is not None:
            xx = xx - lam0
        else:
            dlam = None
            lam0 = 0

        if per_aa:
            kwargs['norm'] = False

        xx, mask = self.get_xx(xx, dlam, xmin, xmax)
        yy, cnt = self.get_yy(mask=mask, **kwargs)
        xx = xx + lam0

        if vdop is not None or vshift is not None or convolve is not None:
            sxx, syy = convol(xx, yy, vdop=vdop, vshift=vshift, custom=convolve)
            yy = np.interp(xx, sxx, syy)

        if per_aa:
            yy = yy * CC / (xx**2) * 1e8

        return xx, yy
    
    #=====================================================================
    def interp(self, newx, vdop=None, vshift=None, convolve=None, **kwargs):
        #-----------------------------------------------------------------
        ib = np.searchsorted(self.lam, newx[0])
        ir = np.searchsorted(self.lam, newx[-1])
        mask = slice(ib, ir)

        xx = self.lam[mask]
        yy, cnt = self.get_yy(mask=mask, **kwargs)

        if vdop is not None or vshift is not None or convolve is not None:
            xx, yy = convol(xx, yy, vdop=vdop, vshift=vshift, custom=convolve)

        yy = np.interp(newx, xx, yy)
        return yy

    #=====================================================================
    def plot_spectrum(self, ang=None, imu=None, norm=True, LTE=False,
                      vdop=None, vshift=None, lam0=None, dlam=10, cont=False,
                      show_lines=False, xmin=None, xmax=None,  per_aa=False, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(lam0=lam0, dlam=dlam, ang=ang, imu=imu, norm=norm,
                           xmin=xmin, xmax=xmax, LTE=LTE, vdop=vdop, vshift=vshift,
                           per_aa=per_aa, cont=cont)
        
        self.set_xaxis(xx, dlam)
        self.set_yaxis(ang, imu, norm, per_aa=per_aa)

        plt.xlabel('Wavelength [Å]')

        if show_lines:
            self.plot_line_labels()

        return plt.plot(xx, yy, **kwargs)
    
    #=====================================================================
    def plot_line_labels(self, lam0=None, nbins=4):
        xmin, xmax, ymin, ymax = plt.axis()
        x = np.array([line.lam0 for line in self.line])
        if lam0 is not None: x = x - lam0
        mask = np.where((x > xmin) * (x < xmax))[0]
        mask = mask[np.argsort(x[mask])[::-1]]
        x = x[mask]
        n = x.size

        edges = np.linspace(xmin, xmax, nbins+1)
        bins = np.zeros_like(edges, dtype=int)
        ibin = np.digitize(x, edges)
        occu = np.zeros_like(x)

        for i, idx in enumerate(ibin):
            bins[idx] += 1
            occu[i] = bins[idx]

        occu = 1 + occu*(ymax-ymin)/20

        segs = np.zeros((n, 2, 2))
        segs[:, 0, 1] = 1.0
        segs[:, 1, 1] = occu
        segs[:, :, 0] = x[:, None]

        line_segments = LineCollection(segs, color='lightgrey')
        plt.gca().add_collection(line_segments)

        for i, idx in enumerate(mask):
            line = self.line[idx]
            x = line.lam0
            if lam0 is not None: x = x - lam0
            plt.text(x=x, y=occu[i], s=idx, rotation=45, color='grey', fontsize=8)

        plt.ylim(ymin, occu.max() + (ymax-ymin)/5)

    #=====================================================================
    def plot_surface(self, i, cmap='hot', **kwargs):
        i3 = self.i3[0]

        if not 'extent' in kwargs.keys():
            if self.xx.size > 1 and self.yy.size > 1:
                extent = [0, self.xx[-1], 0, self.yy[-1]]
                plt.xlabel('Mm')
                plt.ylabel('Mm')

        plt.imshow(i3[:,:,i], cmap=cmap, extent=extent, **kwargs)

    #=====================================================================
    def get_depart(self, ilvl=None):
        #-----------------------------------------------------------------
        mem = self.read_patch_save('atom_patch')[0]
        if mem.ndim==6: mem = mem[...,0]
        lte = mem[...,0]
        dep = mem[...,1]

        if ilvl is not None:
            lte = lte[...,ilvl]
            dep = dep[...,ilvl]

        lte = np.array(lte)
        dep = np.array(dep)
        if self.dim=='1D':
            lte = lte.squeeze()
            dep = dep.squeeze()
        return lte, dep
    
    #=====================================================================
    def plot_depart(self, ilvl, mean='optical', tau_in=None, **kwargs):
        #-----------------------------------------------------------------
        _, dep = self.get_depart(ilvl)
        plt.title('Departure coefficient of level {:.2f}[eV]'.format(self.atom.level[ilvl].ev))
        plt.ylabel(r'$\rm n\,/\,n_{\rm LTE}$')
        if mean=='optical':
            if self.dim=='3D':
                if tau_in is None:
                    ltau = np.mean(self.ltau, axis=(0, 1))[::-1]
                else:
                    ltau = tau_in

                new_dep = interp_tau(ltau.copy(), self.ltau[...,::-1].copy(), dep[...,::-1].copy())
                dep = np.mean(new_dep, axis=(0,1))
            else:
                ltau = self.ltau
                dep = self.get_depart(ilvl)

            plt.xlabel(r'log$\tau_{500}$')
            return plt.plot(ltau, dep, **kwargs)
        
        elif mean=='geometrical':
            if self.dim=='3D':
                dep = np.mean(dep, axis=(0,1))

            plt.xlabel('Height [Mm]')
            return plt.plot(self.zz, dep, **kwargs)
    
    
    #=====================================================================
    def get_cog_corr(self, element, LTE_runs, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(norm=True, **kwargs)
        weq = np.trapz(1-yy, x=xx)*1e3
        abnd = self.abundances[element]

        lte_abnd = np.zeros(len(LTE_runs))
        lte_weqs = np.zeros(len(LTE_runs))
        for i, run in enumerate(LTE_runs):
            lte_abnd[i] = run.abundances[element]
            xx, yy = run.crop(norm=True, LTE=True, **kwargs)
            lte_weqs [i] = np.trapz(1-yy, x=xx)*1e3

        order = np.argsort(lte_weqs)
        lte_weqs = lte_weqs[order]
        lte_abnd = lte_abnd[order]

        a0 = np.interp(weq, lte_weqs, lte_abnd)
        corr = abnd - a0

        return corr, lte_weqs, lte_abnd, weq, abnd

    #=====================================================================
    def plot_cog_corr(self, element, **kwargs):
        #-----------------------------------------------------------------
        corr, lte_weqs, lte_abnd, weq, abnd = self.get_cog_corr(element=element,**kwargs)
        plt.scatter(lte_weqs, lte_abnd, label='LTE', color='cornflowerblue')
        label = 'NLTE\n∆A={:.2f}'.format(corr)
        plt.scatter(weq, abnd, label=label, color='k')
        plt.xlabel(r'W$_{eq}$')
        plt.ylabel('A({})'.format(element))

        a0 = abnd - corr
        plt.scatter(weq, a0, facecolor='none', color='k')
        plt.annotate('', xy=(weq, a0), xytext=(weq, abnd),
                    arrowprops=dict(arrowstyle='<-', color='k'))
        plt.legend(loc='lower right')

    #=====================================================================
    def calc_teff(self):
        #-----------------------------------------------------------------
        if self.n_ang == 1:
            raise ValueError('n_ang = 1! Cannot return effective Temperature from disk center intensity only!')

        xx, yy = self.crop(per_aa=True)
        bol = np.trapz(yy , x=xx)
        teff = (bol / SBOL)**0.25
        return teff

    #=====================================================================
    def read_patch_save(self, kind, concat=True, fdim=3):
        #-----------------------------------------------------------------
        try:
            meta = open(self.sfolder + kind + '_meta.txt')
        except FileNotFoundError:
            print('{:} not found!'.format(self.sfolder + kind + '_meta.txt'))
            if (not kind in self.nml['m3d_params']['save_patch_kind'].split()):
                print('The patch was not saved. Check the save_patch_kind parameter in the m3d_params namelist.')
            raise FileNotFoundError()

        dct = {}
        dct_out = {}
        for line in meta:
            if not line.strip().startswith('*'):
                info = line.split()
                file = info[0]
                prec = int(info[1])
                dtype = '<f{:d}'.format(prec//8)
                dims = tuple([int(d) for d in info[2:]])

                key = file.split('.')[0].split('_')[-1]
                key = np.array([int(k) for k in [key[0:3], key[3:6]]])
                dct[tuple(key)] = np.memmap(self.sfolder + file, dtype=dtype, mode='r', shape=dims, order='F')

        meta.close()

        if concat:
            return self.concat_patches(dct, kind, fdim)
        else:
            return dct

    #=====================================================================
    def concat_patches(self, dct, kind, fdim):
        #-----------------------------------------------------------------
        all_iz  = np.unique([key[0] for key in dct.keys()])
        all_inu = np.unique([key[1] for key in dct.keys()])
        f_indep = all(all_inu == 0)

        n_nu = self.nml['m3d_params']['n_nu']

        def merge_set(inus):
            merge = []
            for iz in all_iz:
                merge_nu = []

                for inu in inus:
                    key = (iz, inu)
                    if key in dct:
                        merge_nu.append(dct[key])
                
                if len(merge_nu) > 0:
                    mem = da.concatenate(merge_nu, axis=fdim)
                    if not f_indep: mem = np.flip(mem, axis=fdim)
                    merge.append(mem)
            
            if len(merge) > 1:
                mem = da.concatenate(merge, axis=2)

            return mem

        if f_indep:
            return [merge_set([0])]
        else:
            set1 = [inu for inu in all_inu if inu > 0 and inu <= n_nu]
            set2 = [inu for inu in all_inu if inu > 0 and inu > n_nu]
            set3 = [inu for inu in all_inu if inu < 0][::-1]
            sets = [s for s in [set1, set2, set3] if len(s)>0]
        
            return [merge_set(s) for s in sets]

#=====================================================================
def weq_diff(run1, run2, marker='.', linestyle='none', **kwargs):
    #-----------------------------------------------------------------
    weq1 = [line.calc_weq() / line.calc_weq(LTE=True) for line in run1.line]
    weq2 = [line.calc_weq() / line.calc_weq(LTE=True) for line in run2.line]
    maxw = np.max([weq1, weq2])
    minw = np.min([weq1, weq2])
    plt.plot([minw, maxw], [minw, maxw], zorder=-10, color='lightgrey', linestyle='--')
    plt.plot(weq1, weq2, marker='.', linestyle='none', **kwargs)

    def compare_dict(x, y):
        diff_x = {k: x[k] for k in x if k in y and x[k] != y[k]}
        diff_y = {k: y[k] for k in y if k in y and x[k] != y[k]}
        return diff_x, diff_y

    def cstrip(val):
        if type(val)==str:
            val = val.strip()
        return val

    x = run1.nml
    y = run2.nml

    diff_dict_x = {}
    diff_dict_y = {}

    for key, val in x.items():
        diff_x, diff_y = compare_dict(x[key], y[key])
        if len(diff_x) > 0:
            diff_dict_x = dict(diff_dict_x, **diff_x)
            diff_dict_y = dict(diff_dict_y, **diff_y)

    xlabel = '\n'.join(['{:}={:}'.format(key, cstrip(val)) for key, val in diff_dict_x.items()])
    ylabel = '\n'.join(['{:}={:}'.format(key, cstrip(val)) for key, val in diff_dict_y.items()])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)