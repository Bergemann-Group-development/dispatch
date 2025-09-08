from .atom import atom_class
from .trans import trans_class
from .auxfuncs import *
from .plotting_scripts import *
from .fast_atom import fast_atom, pandas_atom
from .fast_trans import fast_trans, llist_trans
from .bokeh_scripts import bkh
from .linelist import LineList
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mpl
import f90nml
from scipy.signal import argrelextrema
from scipy.io import FortranFile
from matplotlib.collections import LineCollection
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
import os
import copy
import textwrap
import pathlib

def option(nml, key):
    if key in nml.keys():
        return nml[key]
    else:
        return False

#=====================================================================
class tumult(object):
    @classmethod
    #=================================================================
    def read(cls, folder, verbose=0, read_atmos=True, read_atom=True, fast_read=False, read_llist=False, **kwargs):
        self = cls() # create a new instance of the class
        folder = pathlib.Path(folder)
        sfolder = folder / 'save'
        self.sfolder = sfolder
        self.folder = folder
        self.nml = f90nml.read(folder / 'params.nml')
        self.datadir = pathlib.Path(self.nml['io_params']['datadir'].strip())
        
        df = pd.read_csv(folder / 'abund.txt', header=None, sep=r'\s+', na_filter=False)
        self.abundances = {k.lower().capitalize():v for k, v in df.values}
        
        if os.path.isfile(sfolder / 'out_lam.bin'):
            if (verbose > 0): print('Reading wavelengths')
            self.lam = np.fromfile(sfolder / 'out_lam.bin', dtype=np.float64)
            self.lam = vacuum2obs(self.lam[::-1])
            self.dlam = np.diff(self.lam)

        if os.path.isfile(sfolder / 'bb_lam.bin'):
            self.bb_lam = vacuum2obs(np.fromfile(sfolder / 'bb_lam.bin', dtype=np.float64))[::-1]
        if os.path.isfile(sfolder / 'bf_lam.bin'):
            self.bf_lam = vacuum2obs(np.fromfile(sfolder / 'bf_lam.bin', dtype=np.float64))[::-1]
        if os.path.isfile(sfolder / 'all_lam.bin'):
            self.all_lam = vacuum2obs(np.fromfile(sfolder / 'all_lam.bin', dtype=np.float64))[::-1]
        
        self.has_spectrum = hasattr(self, 'lam')
        if self.has_spectrum:
            if len(self.lam) == 0: self.has_spectrum = False
            if self.nml['m3d_params']['long_scheme'].strip() == 'none': self.has_spectrum = False

        self.read_atom(read_atom, fast_read, verbose, **kwargs)

        file = folder / 'angle_quadrature.txt'
        if os.path.isfile(file):
            df = pd.read_csv(file, header=None, skiprows=1, sep=r'\s+', names=['ang', 'wmu', 'mux', 'muy', 'muz'])
            self.vec = df.values[:,2:]
            self.wts = df.values[:,1] * 2*np.pi
            self.n_ang = self.wts.size
            self.mus, self.imus = np.unique(1 - self.vec[:,-1], return_inverse=True)
            self.mus = 1 - self.mus
            self.nmu = self.mus.size

        make_eos = 'make_eos' in self.nml['m3d_params']
        if make_eos: make_eos = self.nml['m3d_params']['make_eos']
        
        if make_eos:
            mem = self.read_patch_save('eos', fdim=None, lazy=False)[0]
            self.temp = mem[0]
            self.rho = mem[1]
            self.ne = mem[2]
            self.pg = mem[3]
            self.E = mem[4]
            self.chi = mem[5]

        else:
            if self.has_spectrum:
                if (verbose > 0): print('Reading Intensities')
                ie = self.read_patch_save('ie', lazy=False)[0]
                saved_LTE = ie.shape[-1] == 3
                self.mode = 'NLTE' if saved_LTE else 'LTE'

                self.ie     = ie[:,0,0,:,0]
                self.cnt    = ie[:,0,0,:,1]

                self.flux     = np.sum(self.wts * self.ie    , axis=-1)
                self.flux_cnt = np.sum(self.wts * self.cnt   , axis=-1)

                if saved_LTE:
                    self.ie_lte = ie[:,0,0,:,2]
                    self.flux_lte = np.sum(self.wts * self.ie_lte, axis=-1)

                if self.nml['m3d_params']['save_resolved']:
                    i3 = self.read_patch_save('i3')[0]
                    self.i3     = i3[...,0]
                    self.c3     = i3[...,1]
                    if saved_LTE:
                        self.i3_lte = i3[...,2]

            if os.path.isfile(folder / 'linelist_out.txt') and read_llist:
                self.read_llist()
                
            if read_atmos:
                self.read_atmos()
                
        return self

    #=====================================================================
    def read_atom(self, read_atom=True, fast_read=True, verbose=True, **kwargs):
        self.atoms = []
        atom_files = self.folder.glob('atoms/*')
        moldict = {'CH':'C', 'OH':'O', 'NH':'N', 'CN':'N', 'CO':'C',
                    'NO':'N', 'H2O':'O', 'C2':'C', 'N2':'N', 'O2':'O'}

        if 'atom_params' in self.nml:
            self.qnorm = self.nml['atom_params']['qnorm']
            nq_fac = self.nml['atom_params']['nq_fac']
            if self.nml['m3d_params']['m1d_legacy_mode']: nq_fac *= 2
            in_abnd = self.nml['atom_params']['abundance']
            #loggf_cut  = self.nml['atom_params']['logf_cut']
            #red_cut    = self.nml['atom_params']['red_cut']
            #blue_cut   = self.nml['atom_params']['blue_cut']
        else:
            depart = ''
            in_abnd = None

        for i, atom_file in enumerate(atom_files):
            if atom_file.is_file():
                if not read_atom:
                    with open(atom_file, 'r') as f:
                        line = f.readline()
                        while line.startswith('*'):
                            line = f.readline()

                        atom = empty()
                        atom.id = line.split()[0].capitalize()
                        atom.molec = atom.id in moldict.keys()

                else:
                    pandas_format = not atom_file.name.startswith('atom.')

                    #if 'line_mask' in self.nml.keys():
                    #    if 'in_lines' in self.nml['line_mask'].keys():
                    #        in_lines = self.nml['line_mask']['in_lines']
                    #        if type(in_lines)==int: in_lines=[in_lines]
                    #        self.in_lines = np.array(in_lines)-1

                    if pandas_format or fast_read:
                        if pandas_format:
                            atom = pandas_atom.read(atom_file)
                        else:
                            atom = fast_atom(atom_file, **kwargs)

                        atom.run = self
                        if not hasattr(atom, 'in_lines'): atom.in_lines = slice(None, atom.nline)

                        if self.has_spectrum:
                            atom.line = atom.trans[atom.in_lines]

                            if option(self.nml['atom_params'], 'overlapping_bb_in_output'):
                                atom.load_sp(self.lam, qnorm=self.qnorm)
                            else:
                                nnu = self.lam.size
                                ib, ir = atom.get_non_overlapping_indices(atom.line, nnu, nq_fac)
                                atom.bbt['ir'] = 0
                                atom.bbt['ib'] = 0
                                atom.bbt.loc[atom.in_lines, 'ir'] = ir
                                atom.bbt.loc[atom.in_lines, 'ib'] = ib

                            #mask = bbt.loggf > loggf_cut
                            #mask = mask & (bbt.lam_vac_red > blue_cut) 
                            #mask = mask & (bbt.lam_vac_blue < red_cut) 
                            #lines = lines_class(self.atom, bbt[mask], self)
                            #self.nline = len(lines)
                            #self.line = lines

                    else:
                        atom = atom_class(atom_file, qnorm=self.qnorm, verbose=verbose, **kwargs)
                        atom.run = self

                        if not hasattr(atom, 'in_lines'): atom.in_lines = slice(None, atom.nline)
                        atom.line = np.array(atom.trans)[atom.in_lines]
                        atom.nline = len(atom.line)
                            
                        #lines = [line for line in lines if np.log10(line.fosc*line.lower.g) > loggf_cut]
                        #lines = [line for line in lines if line.lam_red > blue_cut and line.lam_blue < red_cut]

                        if self.has_spectrum:
                            if not option(self.nml['atom_params'], 'overlapping_bb_in_output'):
                                nnu = self.lam.size
                                ib, ir = self.atom.get_non_overlapping_indices(atom.line, nnu, nq_fac)
                                atom.ib = ib.astype(int)
                                atom.ir = ir.astype(int)

                if atom.id in self.abundances:
                    atom.abnd = float(self.abundances[atom.id])
                elif atom.molec:
                    atom.abnd = self.abundances[moldict[atom.id]]
                elif in_abnd and in_abnd[i] > -99:
                    atom.abnd = in_abnd[i]

                self.atoms.append(atom)

        self.atoms = np.array(self.atoms)[np.argsort([atom.id for atom in self.atoms])]
        reports = self.folder.glob('convergence_report_*.txt')
        self.conv_report = []
        for report in reports:
            df = pd.read_csv(report, sep=r'\s+')
            self.conv_report.append(df)

        self.has_atom = read_atom and len(self.atoms) > 0
    
    #=====================================================================
    def read_atmos(self):
        atmos_params = self.nml['atmos_params']

        if atmos_params['atmos_format'].strip().lower() == 'save':
            atmos_folder = pathlib.Path(atmos_params['atmos_file'].strip())
            mesh = open(atmos_folder / 'atmos_mesh.txt', 'r')
            afolder = atmos_folder / 'save'
        else:
            mesh = open(self.folder / 'atmos_mesh.txt', 'r')
            afolder = self.sfolder

        amr = 'amr' in atmos_params
        if amr: amr = atmos_params['amr']

        dims = atmos_params['dims']
        self.atmos_cuts = []
        if amr:
            nz = atmos_params['nz']
            self.zz = np.zeros(nz)
            j = 0
            off = 0
            for k in range(dims):
                nx, ny, nnz = [int(i) for i in mesh.readline().split()]
                dx, dy, dz = [float(i) for i in mesh.readline().split()]
                self.zz[j:j+nnz] = np.arange(1,nnz+1)*dz + off
                j += nnz
                off = self.zz[j-1]
                self.atmos_cuts.append(nnz)

            self.atmos_cuts = np.cumsum(self.atmos_cuts)

        else:
            nx, ny, nz = [int(i) for i in mesh.readline().split()]
            dx, dy, dz = [float(i) for i in mesh.readline().split()]
            self.dz = dz
            self.zz = np.arange(1, nz+1)*dz
            self.atmos_cuts = np.cumsum(near_equal_split(nz, dims))

        self.dx = dx
        self.dy = dy
        self.xx = np.arange(nx)*dx
        self.yy = np.arange(ny)*dy
        mesh.close()

        if nx*ny > 1:
            self.dim = '3D'
        else:
            self.dim = '1D'

        if os.path.isfile(afolder / 'atmos_meta.txt'):
                mem = self.read_patch_save('atmos', folder=afolder, zdim=2)[0]
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

    #=====================================================================
    def read_llist(self):
        self.llist = LineList.from_dispatch_format(self.folder / 'linelist_out.txt')
        init = np.frompyfunc(llist_trans, 2, 1)
        self.llist_line = init(np.arange(len(self.llist)), self)

        T = 4000
        elements = [ion.split()[0] for ion in self.llist.ion]
        abnd = np.array([self.abundances[el] if el in self.abundances.keys() else -np.inf for el in elements])
        #strength = np.log10(self.llist.lam0*1e-8 * np.exp(-self.llist.evi*EK/T)) + self.llist.loggf + abnd
        #self.llist['strength2'] = strength

    #=====================================================================
    def get_line(self, lam0, tol=0.05):
        if not hasattr(self.atoms[0], 'line'): self.read_atom()
        all_lines = np.concatenate([[line for line in atom.line] for atom in self.atoms])
        all_lam0 = np.array([line.lam0 for line in all_lines])
        i = np.argmin(abs(all_lam0-lam0))
        line = all_lines[i]
        if abs(line.lam0 - lam0) > tol:
            print(f'WARNING: Line not found! Closest line to {lam0:.2f}Å found at {line.lam0:.2f}Å')

        return line

    #=====================================================================
    def plot_spectrum(self, set_labels=True, show_lines=False, lam_offset=False,
                      label=None, bokeh=False, **kwargs):
        #-----------------------------------------------------------------
        if bokeh:
            pkwargs = bkh.plot_kwargs
        else:
            pkwargs = plt_kwargs

        crop_kwargs = {k:v for k, v in kwargs.items() if k not in pkwargs}
        plot_kwargs = {k:v for k, v in kwargs.items() if k not in crop_kwargs}

        if bokeh:
            fig = bkh.figure(tools='xpan, xwheel_zoom, box_zoom, reset')
            xcrop, _ = self.get_xx(lam_offset=lam_offset, **crop_kwargs)
            self.set_xaxis(xcrop, lam_offset=lam_offset, set_labels=set_labels,
                           bokeh=bokeh, **crop_kwargs)
            crop_kwargs['dlam'] = None
            xx, yy = self.crop(lam_offset=lam_offset, **crop_kwargs)
        else:
            xx, yy = self.crop(lam_offset=lam_offset, **crop_kwargs)
            
            self.set_xaxis(xx, lam_offset=lam_offset, set_labels=set_labels,
                           bokeh=bokeh, **crop_kwargs)
            
        self.set_yaxis(bokeh=bokeh, **crop_kwargs)

        if show_lines:
            self.plot_lline_labels()

        if type(label) is str:
            if '{WEQ}' in label:
                weq = np.trapz(1-yy, x=xx)*1e3
                label = label.replace('{WEQ}', '{:.2f}'.format(weq))

        if bokeh:
            p = bkh.plot(xx, yy, label=label, **plot_kwargs)
        else:
            p = plt.plot(xx, yy, label=label, **plot_kwargs)
            
        return xyplot(p, xx, yy)


    #=====================================================================
    def set_xaxis(self, xx, dlam=None, velocity_units=False,
                  set_labels=True, lam_offset=True, bokeh=False, **kwargs):
        #-----------------------------------------------------------------
        xmin = min(xx)
        xmax = max(xx)
        dxx = xmax - xmin
        offset = lam_offset and dlam

        if velocity_units: 
            xlabel = r'$v_{\rm dop}$ [km/s]'
        elif offset:
            xlabel = r'∆λ [Å]'
        else:
            xlabel = 'Wavelength [Å]'

        if bokeh:
            if set_labels:
                bkh.xlabel(xlabel)

            bkh.xlim(xmin, xmax)
        else:
            if set_labels:
                ax = plt.gca()
                sub = ax.get_subplotspec()
                grid = sub.get_gridspec()
                if sub.rowspan.stop == grid._nrows:
                    plt.xlabel(xlabel)

            plt.xlim(xmin, xmax)
            plt.ticklabel_format(useOffset=False)
            adjust_max_xticks()

    #=====================================================================
    def set_yaxis(self, ang=None, imu=None, norm=True, per_aa=False, renorm=None,
                  set_labels=True, bokeh=False, **kwargs):
        #-----------------------------------------------------------------
        if ang:
            ylabel = 'Intensity µ={:.2f}'.format(self.vec[ang,-1])
        elif imu:
            if imu < 0: imu = self.imus.max() + imu + 1
            ylabel = 'Intensity µ={:.2f}'.format(self.mus[imu])
        else:
            ylabel = 'Flux'

        if per_aa: norm=False

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

        if set_labels:
            if bokeh:
                bkh.ylabel(ylabel)
            else:
                ax = plt.gca()
                sub = ax.get_subplotspec()
                if sub.colspan.start == 0:
                    plt.ylabel(ylabel)

    #=====================================================================
    def get_xx(self, slc=None, lam0=None, dlam=None, xmin=None, xmax=None,
               vshift=None, velocity_units=False, lam_offset=True, **kwargs):
        #-----------------------------------------------------------------
        lam = self.lam

        if vshift:
            vshift = vshift / CC_KMS # vshift has to be given in km/s
            lshift = 1 + vshift
            lam = lam * lshift
    
        if lam0:
            lam = lam - lam0
            if xmin: xmin = xmin - lam0
            if xmax: xmax = xmax - lam0
        else:
            dlam = None
            lam0 = 0

        if slc:
            ib = slc.start
            ir = slc.stop
        else:
            ib = None
            ir = None

        if dlam:
            if not isinstance(dlam, (bool)):
                if hasattr(dlam, '__iter__'):
                    dl = dlam[:]
                else:
                    dl = (dlam, dlam)

                ib = np.searchsorted(lam, -dl[0])
                ir = np.searchsorted(lam, dl[1])

        if xmin:
            ib = max(np.searchsorted(lam, xmin)-1, 0)
            if lam[ib] < xmin-10: ib += 1

        if xmax:
            ir = min(np.searchsorted(lam, xmax)+1, lam.size-1)
            if lam[ir-1] > xmax+10: ir -= 1

        mask = slice(ib, ir)
        xx = lam[mask]

        if lam0:
            xx = xx + lam0

            if velocity_units:
                xx = xx / lam0
                if dlam: xx = xx - 1
                xx = xx * CC_KMS
                
            elif lam_offset and dlam:
                xx = xx - lam0

        return xx, mask

    #=====================================================================
    def get_yy(self, ang=None, imu=None, norm=True, cont=False, LTE=False,
               resolved=False, mask=slice(None, None)):
        #-----------------------------------------------------------------
        if not LTE and self.mode=='LTE':
            print('Cannot get NLTE spectrum, because run was in LTE!')
            LTE = True

        LTE = LTE and hasattr(self, 'ie_lte')

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
            if resolved:
                if cont:
                    yy = self.c3
                elif LTE:
                    yy = self.i3_lte
                else:
                    yy = self.i3

                if ang is None:
                    raise ValueError('ang is required for resolved spectrum')
                else:
                    cnt = np.array(cnt[mask, ..., ang]).T
                    yy = np.array(yy[mask, ..., ang]).T

            else:
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
            
        return yy
    
    #=====================================================================
    def crop(self, slc=None, lam0=None, dlam=None, xmin=None, xmax=None, vdop=None,
             vrot=None, vshift=None, convolve=None, zeta_rt=None, lam_offset=True,
             per_aa=False, renorm=None, fac=None, tilt=None, velocity_units=False, 
             ang=None, imu=None, norm=True, cont=False, LTE=False, line=None, Res=None,  
             resolved=False, **kwargs):
        #-----------------------------------------------------------------
        if line is not None:
            line.prepare()
            slc = slice(line.ib, line.ir)
            lam0 = line.lam0

        if per_aa:
            norm = False

        xx, mask1 = self.get_xx(slc=slc, lam0=lam0, dlam=dlam, xmin=xmin, xmax=xmax,
                               vshift=vshift, velocity_units=velocity_units, lam_offset=lam_offset)
        lam = self.lam[mask1]

        if Res:
            Res = CC_KMS / Res
            if vdop is None:
                vdop = Res
            else:
                vdop = np.sqrt(vdop**2 + Res**2)

        vels = np.array([vdop, vrot, zeta_rt])
        vmask = vels != None
        if any(vmask): # we need to add padding to avoid convolution artefacts at the edges
            uconv = 5*lam[lam.size//2]/CC_KMS
            pad = uconv * np.sum(vels[vmask])

            if dlam:
                dlam = list(dlam + np.array([pad, pad]))

            if xmin: xmin = xmin - pad
            if xmax: xmax = xmax + pad

        _, mask2 = self.get_xx(slc=slc, lam0=lam0, dlam=dlam, xmin=xmin, xmax=xmax,
                               vshift=vshift, velocity_units=velocity_units, lam_offset=lam_offset)
        padded_lam = self.lam[mask2]

        get_yy_dict = {'LTE':LTE, 'norm':norm, 'cont':cont, 'resolved':resolved, 'mask':mask2}
        if hasattr(ang, '__iter__'):
            yy = [self.get_yy(ang=i, **get_yy_dict) for i in ang]
        elif hasattr(imu, '__iter__'):
            yy = [self.get_yy(imu=i, **get_yy_dict) for i in imu]
        else:
            yy = [self.get_yy(ang=ang, imu=imu, **get_yy_dict)]

        yy = np.array(yy)

        if renorm:
            if renorm == 'max':
                yy = yy / np.max(yy)
            elif renorm == 'edge':
                yy = yy / continuum(xx, yy)
            elif renorm > 0:
                _, mnorm = self.get_xx(lam0=lam0, dlam=renorm)
                ynorm = self.get_yy(ang=ang, imu=imu, norm=norm, cont=cont, LTE=LTE, mask=mnorm)
                yy = yy / np.max(ynorm)

        if set([vdop, vrot, zeta_rt, convolve]) != {None}:
            newy = []
            for y in yy:
                sxx, syy = convol(padded_lam, y, vdop=vdop, vrot=vrot, zeta_rt=zeta_rt, custom=convolve)
                newy.append(np.interp(lam, sxx, syy))

            yy = np.array(newy)

        if per_aa:
            yy = yy * CC / (lam**2) * 1e8

        if fac:
            yy = yy * fac
        
        if tilt:
            yy = yy * (1 + tilt*(xx - (xx[0]+xx[-1])/2))

        if len(yy)==1:
            yy = yy[0]

        return xx, yy
    
    #=====================================================================
    def interp(self, newx, **kwargs):
        #-----------------------------------------------------------------
        xmin = newx.min()
        xmax = newx.max()

        xx, yy = self.crop(xmin=xmin, xmax=xmax, **kwargs)

        if yy.ndim==1:
            yy = np.interp(newx, xx, yy)
        else:
            yy = np.array([np.interp(newx, xx, y) for y in yy])

        return yy
    
    #=====================================================================
    def calc_weq(self, **kwargs):
        #-----------------------------------------------------------------
        return calc_weq(self, **kwargs)
            

    #=====================================================================
    def plot_resolved(self, ang=0, f=0, LTE=False, cont=False, **kwargs):
        #-----------------------------------------------------------------
        if LTE:
            i3 = self.i3_lte
        elif cont:
            i3 = self.c3
        else:
            i3 = self.i3

        Z = i3[f, :, :, ang].T
        self.plot_surface(Z, ang=ang, **kwargs)

    #=====================================================================
    def plot_surface(self, Z, ang=None, rasterized=True, shading='nearest', 
                     show_labels=True, cmap='hot', rot90=0, **kwargs):
        #-----------------------------------------------------------------
        xx = self.xx + self.dx/2
        yy = self.yy + self.dy/2

        if ang is None:
            mux, muy, muz = (0.0, 0.0, 1.0)
        else:
            mux, muy, muz = self.vec[ang]

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
            plt.xlabel('X [Mm]')
            plt.ylabel('Y [Mm]')
            plt.title(f'({mux:.2f},{muy:.2f},{muz:.2f})')
        else:
            plt.xticks([])
            plt.yticks([])

    #=====================================================================
    def get_toth(self):
        #-----------------------------------------------------------------
        xmy = np.array([10**(self.abundances[key]-12)*atomic_weights[key]
                        for key in self.abundances.keys()])
        grph = np.sum(xmy) * UU
        toth = self.rho / grph
        return toth

    #=====================================================================
    def get_pop(self, atom=0, ilvl=None):
        #-----------------------------------------------------------------
        if type(atom) is int: atom = self.atoms[atom]
        if ilvl is not None:
            if ilvl < 0: ilvl = atom.nk + ilvl
        mem = self.read_patch_save(f'pop_{atom.id}', zdim=2)[0]
        if mem.ndim==6: mem = mem[...,0]
        lte = mem[...,0]
        dep = mem[...,1]

        if ilvl is not None:
            lte = lte[...,ilvl]
            dep = dep[...,ilvl]

        lte = np.array(lte).squeeze()
        dep = np.array(dep).squeeze()
        nlte = dep*lte

        return lte, nlte

    #=====================================================================
    def get_depart(self, atom=0, ilvl=None, mean=False):
        #-----------------------------------------------------------------
        if type(atom) is int:
            atomid = self.atoms[atom].id
        else:
            atomid = atom

        if hasattr(self, 'runs'):
            xy = np.array([run.get_depart(ilvl, mean=mean) for run in self.runs])

            if mean:
                ltau = np.mean(xy[:,0], axis=0)
                yy1d = [np.interp(ltau, xy[i,0], xy[i,1]) for i in range(len(xy))]
                yy1d = np.mean(yy1d, axis=0)
                return ltau, yy1d
            else:
                return xy[:,0], xy[:,1]
            
        mem = self.read_patch_save(f'pop_{atomid}', zdim=2)[0]

        if mem.ndim==6: mem = mem[...,0]
        dep = mem[...,1]

        if ilvl is not None:
            dep = dep[...,ilvl]

        dep = np.array(dep)[:,:,::-1]
        ltau = self.ltau[...,::-1]

        if self.dim=='1D':
            dep = dep.squeeze()
        elif mean:
            ltau, dep = self.get_ltau_mean(dep, zaxis=2)

        return ltau, dep
    
    #=====================================================================
    def plot_depart(self, atom=0, ilvl=0, mean=True, **kwargs):
        #-----------------------------------------------------------------
        if type(atom) is int:
            atomid = self.atoms[atom].id
        else:
            atomid = atom

        xx, yy = self.get_depart(atomid, ilvl, mean=mean)
        plt.title(f'Departure coefficient of {atomid:} level {ilvl:}')
        plt.ylabel(r'$\rm n\,/\,n_{\rm LTE}$')
        plt.xlabel(r'log$\tau_{500}$')
        return plt.plot(xx.ravel(), yy.ravel(), **kwargs)
    
    
    #=====================================================================
    def get_cog_corr(self, element, LTE_runs, **kwargs):
        #-----------------------------------------------------------------
        weq = self.calc_weq(norm=True, **kwargs)
        abnd = self.abundances[element]

        lte_abnd = np.zeros(len(LTE_runs))
        lte_weqs = np.zeros(len(LTE_runs))
        for i, run in enumerate(LTE_runs):
            lte_abnd[i] = run.abundances[element]
            lte_weqs[i] = run.calc_weq(norm=True, LTE=True, **kwargs)

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
    def calc_teff(self, **kwargs):
        #-----------------------------------------------------------------
        if self.n_ang == 1:
            raise ValueError('n_ang = 1! Cannot return effective Temperature from disk center intensity only!')

        xx, yy = self.crop(norm=False, per_aa=True, **kwargs)
        bol = np.trapz(yy , x=xx)
        teff = (bol / SBOL)**0.25
        return teff

    #=====================================================================
    def read_patch_save(self, kind, folder=None, concat=True, fdim=0, zdim=None, lazy=True):
        #-----------------------------------------------------------------
        if folder: 
            folder = pathlib.Path(folder)
        else:
            folder = self.sfolder

        meta_file = folder / f'{kind}_meta.txt'

        try:
            meta = open(meta_file)
        except FileNotFoundError:
            print(f'{meta_file} not found!')
            if (not kind in self.nml['m3d_params']['save_patch_kind'].split()):
                print('The patch was not saved. Check the save_patch_kind parameter in the m3d_params namelist.')
            raise FileNotFoundError()

        dct = {}
        for line in meta:
            if not line.strip().startswith('*'):
                info = line.split()
                file = info[0]
                prec = int(info[1])
                dtype = '<f{:d}'.format(prec//8)
                dims = tuple([int(d) for d in info[2:]])

                key = file.split('.')[0].split('_')[-1]
                key = np.array([int(k) for k in textwrap.wrap(key, 3)])
                dct[tuple(key)] = np.memmap(folder / file, dtype=dtype, mode='r', shape=dims, order='F')

        meta.close()
        if zdim is None: zdim = len(dims)-1

        if concat:
            return self.concat_patches(dct, kind, fdim, zdim, lazy=lazy)
        else:
            return dct

    #=====================================================================
    def concat_patches(self, dct, kind, fdim, zdim, lazy):
        #-----------------------------------------------------------------
        all_iz  = np.unique([key[0] for key in dct.keys()])
        all_inu = np.unique([key[1] for key in dct.keys()])
        f_indep = all(all_inu == 0) or fdim is None

        n_nu = self.nml['m3d_params']['n_nu']

        if lazy:
            concat = memmap_list
        else:
            concat = np.concatenate

        def merge_set(inus):
            merge = []
            for iz in all_iz:
                merge_nu = []

                for inu in inus:
                    key = (iz, inu)
                    if key in dct:
                        merge_nu.append(dct[key])
                
                if len(merge_nu) > 0:
                    if len(merge_nu) == 1:
                        mem = merge_nu[0]
                    else:
                        mem = concat(merge_nu, axis=fdim)
                        
                    if not f_indep: mem = np.flip(mem, axis=fdim)
                    merge.append(mem)
            
            if len(merge) > 1:
                mem = concat(merge, axis=zdim)

            return mem

        if f_indep:
            return [merge_set(all_inu)]
        else:
            set1 = [inu for inu in all_inu if inu > 0 and inu <= n_nu]
            set2 = [inu for inu in all_inu if inu > 0 and inu > n_nu]
            set3 = [inu for inu in all_inu if inu < 0][::-1]
            sets = [s for s in [set1, set2, set3] if len(s)>0]
        
            return [merge_set(s) for s in sets]
        

    #=====================================================================
    def interp_ltau_surf(self, cube, newx=None, log=False, zaxis=-1):
        #-----------------------------------------------------------------
        if zaxis != -1: cube = np.moveaxis(cube, zaxis, -1)
        ltau = self.ltau
        step = -1 if ltau[...,0].ravel()[0] > ltau[...,-1].ravel()[0] else 1
        xx = ltau[...,::step].copy()
        yy = cube[...,::step].copy()
        if log: yy = np.log(yy)
        if newx is None:
            axis = tuple(i for i in range(xx.ndim-1))
            newx = np.mean(xx, axis=axis)
        newy = interp_tau(newx, xx, yy)
        if log: newy = np.exp(newy)
        if zaxis != -1: newy = np.moveaxis(newy, -1, zaxis)
        return newx, newy
    
    #=====================================================================
    def plot_geo_cube(self, cube, marker=',', linestyle='none', rasterized=True,
                       xmin=None, xmax=None, ymin=None, ymax=None, **kwargs):
        #-----------------------------------------------------------------
        xx = (self.height * np.ones_like(self.ltau)).ravel()
        yy = cube.ravel()
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
    def plot_ltau_cube(self, cube, marker=',', linestyle='none', rasterized=True,
                       xmin=None, xmax=None, ymin=None, ymax=None, density=False,
                       cmap='GnBu', **kwargs):
        #-----------------------------------------------------------------
        xx = self.ltau.ravel()
        yy = cube.ravel()
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        plt.xlim(np.nanmin(xx), np.nanmax(xx))
        plt.ylim(np.nanmin(yy), np.nanmax(yy))

        if density:
            return scatter_density(xx, yy, cmap=cmap, **kwargs)
        else:
            return plt.plot(xx, yy, marker=marker, linestyle=linestyle,
                            rasterized=rasterized, **kwargs)

    #=====================================================================
    def get_ltau_mean(self, cube, log=False, newx=None, **kwargs):
        #-----------------------------------------------------------------
        if type(cube) is str:
            cube = getattr(self, cube)

        if self.dim == '3D':
            newx, yy3d = self.interp_ltau_surf(cube, log=log, newx=newx, **kwargs)
            if log: yy3d = np.log(yy3d)
            axis = tuple(i for i in range(yy3d.ndim-1))
            yy1d = np.mean(yy3d, axis=axis)
            if log: yy1d = np.exp(yy1d)
        else:
            if newx is None: newx = self.ltau[::-1]
            step = -1 if self.ltau[0] > self.ltau[-1] else 1
            yy1d = np.interp(newx, self.ltau[::step], cube[::step])

        return newx, yy1d
    
    #=====================================================================
    def plot_ltau_mean(self, yy, log=False, newx=None, xmin=None, xmax=None,
                       ymin=None, ymax=None, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.get_ltau_mean(yy, log=log, newx=newx)
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        return plt.plot(xx, yy, **kwargs)
    
    #=====================================================================
    def get_geo_mean(self, yy, newx=None, log=False):
        if hasattr(self, 'runs'):
            xx = np.mean([run.height for run in self.runs], axis=0)
            yy = np.mean([run.get_geo_mean(yy, newx=xx) for run in self.runs], axis=0)
            return xx, yy
        
        if type(yy) is str:
            yy = getattr(self, yy)

        if newx is None:
            xx = self.height

        if self.dim == '3D':
            if log: yy = np.log(yy)
            yy = np.mean(yy, axis=(0,1))
            if log: yy = np.exp(yy)
            
        if xx is not self.height:
            yy = np.interp(xx, self.height, yy)

        return xx, yy
    
    #=====================================================================
    def plot_geo_mean(self, yy, log=False, xmin=None, xmax=None,
                       ymin=None, ymax=None, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.get_geo_mean(yy, log=log)
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        return plt.plot(xx, yy, **kwargs)
    
    #=================================================================
    def get_ltau_surface(self, ltau, cube, log=False, zaxis=-1):
        if hasattr(self, 'runs'):
            yy = np.array([run.get_ltau_surface(ltau, cube, log=log, zaxis=zaxis) for run in self.runs])
        else:
            if type(cube) is str:
                cube = getattr(self, cube)

            if self.dim == '3D':
                newx = np.array([ltau])
                _, yy = self.interp_ltau_surf(cube, newx=newx, log=log, zaxis=zaxis)
            else:
                ltau1d = self.ltau
                step = -1 if ltau1d[0] > ltau1d[-1] else 1
                yy = np.interp(ltau, ltau1d[::step], cube[::step])
                
        return yy

    #=================================================================
    def plot_ltau_hist(self, ltau, cube, log=False, zaxis=-1, show_mean=False, color=None, **kwargs):
        yy = self.get_ltau_surface(ltau, cube, log=log, zaxis=zaxis)

        if show_mean:
            if log:
                mean = 10**np.mean(np.log10(yy))
            else:
                mean = np.mean(yy)
                
            plt.axvline(mean, ls='--', color=color)

        return plt.hist(yy.ravel(), color=color, **kwargs)
    
    #=================================================================
    def get_clv(self, center_diff=False, **kwargs):
        xx = self.mus
        yy = self.calc_weq(imu=np.arange(len(xx)), **kwargs)
        if center_diff: yy = yy - yy[0]
        return xx, yy
    
    #=====================================================================
    def plot_line_labels(self, offset=None, tol=0.02, show_atom=True,
                         show_llist=True, dry=False, xmin=None, xmax=None, **kwargs):
        #-----------------------------------------------------------------
        if xmin is None or xmax is None:
            ax = plt.gca()
            xmin, xmax = ax.get_xlim()
            if offset is not None:
                xmin = xmin - offset
                xmax = xmax - offset

        if hasattr(self, 'llist') and show_llist:
            if tol is None:
                llist = self.llist
                mask = (llist.lam0 > xmin) & (llist.lam0 < xmax)
                llist = llist[mask]
            else:
                df = self.llist[['ion', 'lam0', 'strength', 'loggf', 'ga']]
                mask = (df.lam0 > xmin) & (df.lam0 < xmax)
                df = df[mask]
                df = df.sort_values(by=['ion', 'lam0']).reset_index(drop=True)
                species = np.unique(df['ion'])

                llist = pd.DataFrame()
                for s in species:
                    sub = df[df['ion'] == s]
                    if len(sub) == 0:
                        continue
                    elif len(sub) == 1:
                        llist = pd.concat([llist, sub])
                    else:
                        groups = assign_group_indices(sub['lam0'].values, tol) - 1
                        for i in range(groups.max()+1):
                            group = sub[groups == i]
                            if len(group) == 1:
                                llist = pd.concat([llist, group])
                            else:
                                strength = group.strength.sum()
                                wts = group.strength / strength
                                lam0 = np.sum(group.lam0 * wts)
                                loggf = np.log10(np.sum(10**group.loggf))
                                new_row = pd.DataFrame({'ion': [s], 'lam0': [lam0], 'strength': [strength], 'loggf': [loggf]})
                                llist = pd.concat([llist, new_row])

            max_strength = llist.strength.max()
        else:
            show_atom = True
            max_strength = 0.1
            llist = pd.DataFrame()

        if self.has_atom and show_atom:
            bbt = self.atom.bbt[['lam0', 'loggf', 'ga']].copy().reset_index(drop=True)
            mask = (bbt.lam0 > xmin) & (bbt.lam0 < xmax)
            bbt = bbt[mask]
            bbt['ion'] = self.atom.id
            bbt['strength'] = max_strength
            if len(bbt) == 1:
                llist = pd.concat([llist, bbt])
            elif len(bbt) > 1 and tol is not None:
                groups = assign_group_indices(bbt['lam0'].values, tol) - 1
                for i in range(groups.max()+1):
                    group = bbt[groups == i]
                    if len(group) == 1:
                        llist = pd.concat([llist, group])
                    else:
                        strength = group.strength.sum()
                        wts = group.strength / strength
                        lam0 = np.sum(group.lam0 * wts)
                        loggf = np.log10(np.sum(10**group.loggf))
                        new_row = pd.DataFrame({'ion': [self.atom.id], 'lam0': [lam0], 'strength': [strength], 'loggf': [loggf]})
                        llist = pd.concat([llist, new_row])

        if len(llist) > 0:    
            llist = llist.sort_values(by='lam0').reset_index(drop=True)
            if offset is not None: llist['lam0'] = llist['lam0'] + offset
            if not dry:
                plot_line_annotations(llist.lam0.values, llist.ion.values, llist.strength.values, **kwargs)
        
        return llist
    
    #=====================================================================
    def plot_line_stem(self, bottom=1.0, scale=None, offset=None, xmin=None, xmax=None, **kwargs):
        def cog_func(xx):
            mask = xx < np.exp(0.5)
            yy = np.zeros_like(xx)
            yy[mask] = xx[mask] / (2*np.exp(1))**0.5
            yy[~mask] = np.log(xx[~mask])**0.5
            return yy
        
        if scale is None:
            ymin, ymax = plt.gca().get_ylim()
            scale = - 0.2 * (ymax - ymin)
        
        llist = self.plot_line_labels(offset=offset, xmin=xmin, xmax=xmax, show_atom=False, dry=True)
        xx = llist.lam0
        yy = cog_func(llist.strength.values)
        yy = yy - np.min(yy)
        yy = yy / np.max(yy)
        yy = bottom + yy * scale
        plt.plot([xx, xx], [bottom*np.ones_like(xx), yy], **kwargs)
    
    #=================================================================
    def fit_observed(self, oxx, oyy, fit=None, bounds=None, p0=None, sigma=None,
                     syst=None, method='curve_fit', tol=None, Res='auto', **kwargs):
        #-----------------------------------------------------------------
        fit = fit or []
        bounds = bounds or {}
        p0 = p0 or {}
        
        if len(fit)==0: return {}, None
        lbounds = {'vshift':-3, 'vdop':1e-3, 'vrot':1e-3, 'zeta_rt':1e-3, 'fac':0.5, 'tilt':-0.2}
        ubounds = {'vshift':+3, 'vdop':1e1,  'vrot':1e1,  'zeta_rt':1e1,  'fac':2.0, 'tilt':+0.2}
        p0s     = {'vshift':0,  'vdop':1,    'vrot':1,    'zeta_rt':1,    'fac':1.0, 'tilt':0.0}
        
        p0s.update(p0)
        for key in bounds.keys():
            lbounds[key], ubounds[key] = bounds[key]

        p0 = [p0s[f] for f in fit]
        bounds = ([lbounds[f] for f in fit], [ubounds[f] for f in fit])

        if method != 'curve_fit':
            bounds = np.vstack(bounds).T

        if Res == 'auto':
            Res = np.mean(oxx) / np.mean(np.diff(oxx))

        def opt_func(newx, *args):
            ka = {f:args[i] for i, f in enumerate(fit)}
            return self.interp(newx, **ka, Res=Res, **kwargs)
        
        def fit_params(x, y):
            if method == 'curve_fit':
                pfit, _ = optimize.curve_fit(opt_func, x, y, p0=p0, bounds=bounds, sigma=sigma)

            else:
                def opt_func2(params, newx):
                    res = y - opt_func(newx, *params)
                    if sigma is not None: res = res/sigma
                    return np.sum(res**2)
                
                pfit = optimize.minimize(opt_func2, p0, args=(x), bounds=bounds, method=method, tol=tol).x

            return pfit

        pfit = fit_params(oxx, oyy)
        fitted = opt_func(oxx, *pfit)
        best = {f:pfit[i] for i, f in enumerate(fit)}

        if syst is None:
            return best, fitted
        else:
            pfits = np.zeros((len(pfit), 3))
            pfits[:, 0] = fit_params(oxx, oyy*(1-syst))
            pfits[:, 1] = pfit
            pfits[:, 2] = fit_params(oxx, oyy*(1+syst))
            perr = [pchip_sigma([-syst, 0, syst], p, 0, syst) for p in pfits]
            perr = {f:perr[i] for i, f in enumerate(fit)}
            return best, fitted, perr



    
    #=================================================================
    def plot_fit(self, oxx, oyy, fit=[], bounds=None, p0=None, sigma=None, method='curve_fit', label=None, **kwargs):
        crop_kwargs = {k:v for k, v in kwargs.items() if k not in plt_kwargs}
        plot_kwargs = {k:v for k, v in kwargs.items() if k not in crop_kwargs}
        p, fitted = self.fit_observed(oxx, oyy, fit=fit, bounds=bounds, p0=p0, sigma=sigma, method=method, **crop_kwargs)
        if label is None:
            label = ', '.join([f'{k}={v:.2f}' for k, v in p.items()])
        else:
            label = ', '.join([label,] + [f'{k}={v:.2f}' for k, v in p.items()])
        return plt.plot(oxx, fitted, label=label, **plot_kwargs)
    
    #=================================================================
    def to_multi1d_atmos(self, dir, name, vz=0.0, vmic=1.0):
        df = pd.DataFrame()
        df['Height'] = self.zz*1e3
        df['Temp'] = self.temp
        df['Ne'] = self.ne
        df['V'] = vz
        df['Vturb'] = vmic

        f = open(dir + f'atmos.{name:}', 'w+')
        f.write(f'{name:}\n')
        f.write('HEIGHT\n')
        f.write('4.4\n')
        f.write(f'{len(df):}\n')
        text = df[::-1].to_string(index=False, float_format='%.4e')
        text = '*' + text[1:]
        f.write(text)
        f.close()

        f = open(dir + f'dscale.{name:}', 'w+')
        f.write(f'{name:}\n')
        f.write('HEIGHT\n')
        f.write(f'{len(df):} -1.0000e+01\n')
        text = df.Height[::-1].to_string(index=False, float_format='%.4e')
        f.write(text)
        f.close()

    #=================================================================
    def to_dispatch1d_atmos(self, file, vmic=1.0, FeH=''):
        ltau, temp = self.get_ltau_mean('temp', log=True)
        _, ne      = self.get_ltau_mean('ne',   log=True)
        _, rho     = self.get_ltau_mean('rho',  log=True)

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
    def plot_slice(self, var, ix=None, iy=None, iz=None, shading='nearest',
                cmap='hot', rasterized=True, logNorm=False, contrast=False,
                show_cuts=False, title=None, tau_contours=False):
    
        if type(var) is str:
            var = getattr(self, var)

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
    
    #=================================================================
    def get_duration(self, string=True):
        from datetime import timedelta

        file = self.folder / 'rank_00000.disp'
        with open(file, 'rb') as f:
            lines = f.readlines()
            seconds = float(lines[-1].strip().split()[0])
            time_delta = timedelta(seconds=seconds)
            if string: time_delta = str(time_delta)
            return time_delta

    def __str__(self):
        if hasattr(self, 'label'):
            return self.label
        else:
            return self.folder
        
    def rename_folder(self, new_folder):
        new_folder = pathlib.Path(new_folder)
        self.folder = self.folder.rename(new_folder)
        self.sfolder = new_folder / 'save'
    

    @classmethod
    #=================================================================
    def combine(self, runs, weights=None):
        #-------------------------------------------------------------
        self = tumult() # create a new instance of the class
        self.run0 = runs[0]
        self.folder = '\n'.join([str(run.folder) for run in runs])

        def check_and_combine(key):
            if hasattr(self.run0, key):
                atts = [getattr(run, key) for run in runs]
                same = all([np.array_equal(a, atts[0]) for a in atts])
                if not same:
                    print('Warning: {:} is not the same for all runs!'.format(key))
                else:
                    setattr(self, key, atts[0])

        self.runs = runs
        self.weights = weights
        self.nml = self.run0.nml
        check_and_combine('mode')
        check_and_combine('mus')
        check_and_combine('imus')
        check_and_combine('nmu')
        check_and_combine('n_ang')
        check_and_combine('vec')
        check_and_combine('wts')
        check_and_combine('nline')
        check_and_combine('dim')

        if hasattr(runs[0], 'ltau'):
            if weights is not None:
                taumin = np.min([run.ltau for run in runs])
                taumax = np.max([run.ltau for run in runs])
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
                check_and_combine('abundances')
                self.ne = np.array([run.ne for run in runs])
                self.rho = np.array([run.rho for run in runs])
                self.temp = np.array([run.temp for run in runs])
                self.vx = np.array([run.vx for run in runs])
                self.vy = np.array([run.vy for run in runs])
                self.vz = np.array([run.vz for run in runs])
                self.chi = np.array([run.chi for run in runs])
                self.tau = np.array([run.tau for run in runs])
                self.ltau = np.array([run.ltau for run in runs])
                weights = np.ones(len(runs))/len(runs)

        if self.run0.has_atom:
            self.atom = copy.copy(self.run0.atom)
            self.atom.run = self

            if type(self.atom.trans[0]) is fast_trans:
                init = np.frompyfunc(fast_trans, 2, 1)
                self.atom.trans = init(np.arange(self.atom.nline), self.atom)
            else:
                self.atom.trans = [copy.copy(trans) for trans in self.run0.atom.trans]

            self.in_lines = self.run0.in_lines
            self.line = np.array(self.atom.trans)[self.in_lines]

        if self.run0.has_spectrum:
            self.lam    = self.run0.lam
            self.ie     = np.sum([run.ie*w for run, w in zip(runs, weights)], axis=0)
            self.cnt    = np.sum([run.cnt*w for run, w in zip(runs, weights)], axis=0)

            self.flux     = np.sum([run.flux*w for run, w in zip(runs, weights)], axis=0)
            self.flux_cnt = np.sum([run.flux_cnt*w for run, w in zip(runs, weights)], axis=0)

            if hasattr(self.run0, 'ie_lte'):
                self.ie_lte = np.sum([run.ie_lte*w for run, w in zip(runs, weights)], axis=0)
                self.flux_lte = np.sum([run.flux_lte*w for run, w in zip(runs, weights)], axis=0)

        return self

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