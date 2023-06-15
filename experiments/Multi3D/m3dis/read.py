from .atom import atom_class, CC
from .trans import trans_class
from .auxfuncs import vacuum2obs, continuum, freq_quad, convol
import numpy as np
import pandas as pd
import dask.array as da
import matplotlib.pyplot as plt
import glob
import f90nml
from scipy.signal import argrelextrema
from scipy.io import FortranFile
from matplotlib.collections import LineCollection
import os

cc = 2.99792458e10 # Speed of light [cm/s]

#=====================================================================
class read(object):
    #=================================================================
    def __init__(self, folder, verbose=0):
        #-------------------------------------------------------------
        folder = folder + os.sep
        sfolder = folder + 'save' + os.sep
        self.sfolder = sfolder
        self.nml = f90nml.read(folder + 'params.nml')
        
        df = pd.read_csv(folder + 'abund.txt', header=None, sep='\s+', na_filter=False)
        self.abundances = {k.lower().capitalize():v for k, v in df.values}
        
        if os.path.isfile(folder + 'atom.txt'):
            self.qnorm = self.nml['atom_params']['qnorm']
            self.atom = atom_class(folder + 'atom.txt', qnorm=self.qnorm, verbose=verbose)
            loggf_cut = self.nml['atom_params']['logf_cut']
            red_cut   = self.nml['atom_params']['red_cut']
            blue_cut  = self.nml['atom_params']['blue_cut']
            in_abnd   = self.nml['atom_params']['abundance']
            
            if self.atom.id in self.abundances:
                self.atom.abnd = self.abundances[self.atom.id]
            elif in_abnd > -99:
                self.atom.abnd = in_abnd

            self.atom.title = 'A({}) = {:.2f}'.format(self.atom.id, self.atom.abnd)
        
            if 'line_mask' in self.nml.keys():
                in_lines = self.nml['line_mask']['in_lines']
                lines = [self.atom.trans[i-1] for i in in_lines]
            else:
                lines = self.atom.trans[:self.atom.nline]
                
            lines = [line for line in lines if np.log10(line.fosc*line.lower.g) > loggf_cut]
            lines = [line for line in lines if line.lam_red > blue_cut and line.lam_blue < red_cut]
            
            self.nline = len(lines)
            self.line = lines
        
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
        
        if (verbose > 0): print('Reading Intensities')
        file = folder + 'angle_quadrature.txt'
        df = pd.read_csv(file, header=None, skiprows=1, sep='\s+', names=['ang', 'wmu', 'mux', 'muy', 'muz'])
        self.vec = df.values[:,2:]
        self.wts = df.values[:,1]
        self.n_ang = self.wts.size
        self.mus, self.imus = np.unique(1 - self.vec[:,-1], return_inverse=True)
        self.mus = 1 - self.mus
        self.nmu = self.mus.size
        
        self.i3 = [None]*self.n_ang
        self.c3 = [None]*self.n_ang
        self.ie  = [None]*self.n_ang
        self.cnt = [None]*self.n_ang
        self.i3_lte = [None]*self.n_ang
        self.ie_lte = [None]*self.n_ang

        out_rad = self.read_patch_save('out_rad', concat=False)
        ie_dict = self.read_patch_save('ie',      concat=False)
        ipos = list(out_rad.keys())[0][0]
        nmu  = np.min([key[2] for key in out_rad.keys()])

        self.flux = np.zeros_like(self.lam)
        self.flux_lte = np.zeros_like(self.lam)
        self.flux_cnt = np.zeros_like(self.lam)

        out_rad_keys = list(out_rad.keys())
        min_ang = min([key[-1] for key in out_rad_keys])

        for i in range(self.n_ang):
            key = out_rad_keys[i]
            ipos, inu, ang = key
            mem = out_rad[key]
            ang = ang - min_ang
            self.i3[ang]     = mem[...,0,:,0]
            self.i3_lte[ang] = mem[...,0,:,1]
            self.c3[ang]     = mem[...,0,:,2]

            mem = ie_dict[key]
            self.ie[ang]     = mem[0,0,0,:,0]
            self.ie_lte[ang] = mem[0,0,0,:,1]
            self.cnt[ang]    = mem[0,0,0,:,2]

            self.flux     += self.wts[ang] * self.ie[ang]
            self.flux_lte += self.wts[ang] * self.ie_lte[ang]
            self.flux_cnt += self.wts[ang] * self.cnt[ang]
        
        dim = self.i3[0].shape
            
        if (np.product(dim[:2]) > 1):
            self.dim = '3D'
        else:
            self.dim = '1D'
            
        if hasattr(self, 'atom'):
            #print('Prepairing transitions')
            for line in self.line:
                line.run = self
                line.prepared = False
            
        if os.path.isfile(folder + 'atmos_mesh.txt'):
            mesh = open(folder + 'atmos_mesh.txt', 'r')
            nx, ny, nz = [int(i) for i in mesh.readline().split()]
            dx, dy, dz = [float(i) for i in mesh.readline().split()]
            mesh.close()
            
            self.xx = np.linspace(0, (nx-1)*dx, nx)
            self.yy = np.linspace(0, (ny-1)*dy, ny)
            self.zz = np.linspace(0, (nz-1)*dz, nz)
            
            dat = FortranFile(folder + 'atmos_data.bin')
            self.temp = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.rho  = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.px   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.py   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.pz   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.pe   = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.tau  = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.chi  = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            self.eta  = dat.read_reals(np.float32).reshape((nx,ny,nz), order='F')
            dat.close()
        
        maxiter = self.nml['m3d_params']['maxiter']
        depart = self.nml['atom_params']['precomputed_depart'].strip()
        if (maxiter > 0):
            self.mode = 'NLTE'
        elif(depart==''):
            self.mode = 'LTE'
        else:
            self.mode = 'NLTE'

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
    def set_xaxis(self, xx, dlam=None):
        #-----------------------------------------------------------------
        if dlam is not None:
            plt.xlim(min(xx), max(xx))
            plt.xlabel('$\Delta\lambda$ [Å]')
        else:
            diff = 0.05*(xx[0] - xx[-1])
            plt.xlim(xx[0]-diff, xx[-1]+diff)
            plt.xlabel('Wavelength [Å]')

        plt.ticklabel_format(useOffset=False)

    #=====================================================================
    def set_yaxis(self, ang, imu, norm, renorm=False):
        #-----------------------------------------------------------------
        if ang is not None:
            ylabel = 'Intensity µ={:.2f}'.format(self.vec[ang,-1])
            
        elif imu is not None:
            mus, idx = np.unique(1-self.vec[:,-1], return_inverse=True)
            mu = 1-mus[imu]
            ylabel = 'Intensity µ={:.2f}'.format(mu)
            
        else:
            ylabel = 'Flux'

        if norm or renorm:
            plt.ylabel('Normalised ' + ylabel)
        else:
            plt.ylabel(ylabel)

    #=====================================================================
    def get_yy(self, ang=None, imu=None, norm=True, LTE=False, mask=slice(None, None)):
        #-----------------------------------------------------------------
        if imu is not None:
            if imu >= self.nmu: raise ValueError('imu is out of range. Nmu = {:}'.format(self.nmu))
        if self.n_ang == 1: ang=0

        if ang is None and imu is None:
            cnt = self.flux_cnt
            if LTE:
                yy = self.flux_lte
            else:
                yy = self.flux

            cnt = cnt[mask]
            yy  = yy[mask]

        else:
            cnt = self.cnt
            if LTE:
                yy = self.ie_lte
            else:
                yy = self.ie

            if ang is not None:
                cnt = cnt[ang][mask]
                yy = yy[ang][mask]
            
            elif imu is not None:
                cnt = [c[mask] for i, c in enumerate(cnt) if self.imus[i]==imu]
                yy  = [y[mask] for i, y in enumerate(yy)  if self.imus[i]==imu]
                cnt = np.mean(cnt, axis=0)
                yy  = np.mean(yy,  axis=0)
        
        if norm:
            yy = yy / cnt
            
        return yy, cnt
    
    #=====================================================================
    def crop(self, lam0=None, dlam=None, xmin=None, xmax=None, vdop=None, 
             vshift=None, convolve=None, **kwargs):
        #-----------------------------------------------------------------
        xx = self.lam

        if lam0 is not None:
            xx = xx - lam0
        else:
            dlam = None
            lam0 = 0

        xx, mask = self.get_xx(xx, dlam, xmin, xmax)
        yy, cnt = self.get_yy(mask=mask, **kwargs)
        xx = xx + lam0

        if vdop is not None or vshift is not None or convolve is not None:
            sxx, syy = convol(xx, yy, vdop=vdop, vshift=vshift, custom=convolve)
            yy = np.interp(xx, sxx, syy)

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
                      vdop=None, vshift=None, lam0=None, dlam=10, 
                      show_lines=False, xmin=None, xmax=None,  **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.crop(lam0=lam0, dlam=dlam, ang=ang, imu=imu, norm=norm,
                           xmin=xmin, xmax=xmax, LTE=LTE, vdop=vdop, vshift=vshift)
        
        self.set_xaxis(xx, dlam)
        self.set_yaxis(ang, imu, norm)

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
    def get_populations(self):
        mem = self.read_patch_save('atom_patch')[0]
        lte  = mem[...,0]
        nlte = mem[...,1]
        return lte, nlte

    #=====================================================================
    def get_depart(self, ilvl, mean=True, xxin=None):
        atom_patch = self.read_patch_save('atom_patch', concat=False)
        key = [1,0,0]
        xx = []
        yy = []
        bot = 0
        bots = []
        while tuple(key) in atom_patch:
            val = atom_patch[tuple(key)]
            nx, ny, nz = val.shape[:3]
            level = nx/5
            xx = np.concatenate([xx, bot + np.arange(nz)/level])
            bot = bot + nz/level
            nlte = val[:,:,:,ilvl,1]
            lte  = val[:,:,:,ilvl,0]

            if mean:
                nlte = np.mean(nlte, axis=(0,1))
                lte  = np.mean(lte,  axis=(0,1))

            if yy == []:
                yy = nlte/lte
            else:
                yy = np.concatenate([yy, nlte/lte], axis=-1)

            key[0] = key[0] + 1
            bots.append(bot)

        xx = xx/xx[-1]

        if xxin is not None:
            xx = xxin

        return xx, yy, bots

    #=====================================================================
    def plot_depart(self, ilvl, show_splits=False, xxin=None, **kwargs):
        xx, yy, bots = self.get_depart(ilvl, mean=True, xxin=xxin)

        plt.plot(xx, yy, **kwargs)

        if show_splits:
            bots = np.array(bots) / xx[-1]
            for bot in bots:
                plt.axvline(bot, linewidth=0.2)

        return xx, yy
    
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
    def read_patch_save(self, kind, concat=True):
        #-----------------------------------------------------------------
        meta = open(self.sfolder + kind + '_meta.txt')

        dct = {}
        dct_out = {}
        for line in meta:
            if not line.strip().startswith('*'):
                info = line.split()
                file = info[0]
                prec = int(info[1])
                dims = [int(d) for d in info[2:]]

                if prec == 32: dtype=np.float32
                if prec == 64: dtype=np.float64

                key = file.split('.')[0].split('_')[-1]
                key = [int(k) for k in [key[0:3], key[3:6], key[6:9]]]

                #if '_out' in file: key = key + ['out']

                dct[tuple(key)] = np.memmap(self.sfolder + file, dtype=dtype, mode='r').reshape(dims, order='F')

        meta.close()

        all_inu = np.unique([key[1] for key in dct.keys()])
        f_indep = all_inu.sum() == 0
        if kind in ['out_rad', 'ie']: f_indep = False
        if kind in ['atom_rate']: f_indep = True

        if concat:
            return self.concat_patches(dct, kind, f_indep)

        elif not f_indep:
            for key, val in dct.items():
                dct[key] = val[..., ::-1, :]
            
            return dct

        else:
            return dct

    #=====================================================================
    def concat_patches(self, dct, kind, f_indep):
        #-----------------------------------------------------------------
        all_iz  = np.unique([key[0] for key in dct.keys()])
        all_inu = np.unique([key[1] for key in dct.keys()])
        all_imu = np.unique([key[2] for key in dct.keys()])

        n_nu = self.nml['m3d_params']['n_nu']

        internal = {}
        output = {}

        for imu in all_imu:
            merge_in = []
            merge_out = []
            for iz in all_iz:
                if f_indep:
                    key = (iz, all_inu[0], imu)
                    if key in dct:
                        merge_out.append(dct[key])

                else:
                    merge_nu = []
                    for inu in range(1, n_nu+1):
                        key = (iz, inu, imu)
                        if key in dct:
                            merge_nu.append(dct[key])

                    if len(merge_nu) > 0:
                        internal[(iz,imu)] = da.concatenate(merge_nu, axis=-2)
                        merge_in.append(internal[(iz,imu)])
                    
                    merge_nu = []
                    for inu in range(n_nu+1, 2*n_nu+1):
                        key = (iz, inu, imu)
                        if key in dct:
                            merge_nu.append(dct[key])

                    if len(merge_nu) > 0:
                        output[(iz,imu)] = da.concatenate(merge_nu, axis=-2)
                        merge_out.append(output[(iz,imu)])

            if len(merge_in) > 0:
                internal[(imu)] = da.concatenate(merge_in, axis=2)

            if len(merge_out) > 0:
                output[(imu)] = da.concatenate(merge_out, axis=2)

        if f_indep:
            return output

        else:
            for d in [internal, output]:
                keys = [key for key in d.keys() if type(key) is tuple]
                for key in keys: del d[key]

                for key, val in d.items():
                    d[key] = val[..., ::-1, :]

            return internal, output
