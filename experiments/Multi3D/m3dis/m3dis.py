from .atom import atom_class
from .trans import trans_class
from .auxfuncs import *
from .fast_atom import fast_atom, pandas_atom
from .fast_trans import fast_trans, llist_trans
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mpl
import glob
import f90nml
from scipy.signal import argrelextrema
from scipy.io import FortranFile
from matplotlib.collections import LineCollection
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
import os
import copy

def option(nml, key):
    if key in nml.keys():
        return nml[key]
    else:
        return False

def launch(input_name, nml_settings=None, default=None, overwrite=False, threads=32, log_file=None, slurm_script=None, verbose=True, dry=False):
    cwd = os.getcwd()
    check_file(cwd + '/dispatch.x')
    if default is not None:
        if type(default) is f90nml.namelist.Namelist:
            input = copy.deepcopy(default)
        elif default.endswith('.nml'):
            input = f90nml.read(default)
        else:
            input = f90nml.read(default + '/params.nml')
    else:
        input = f90nml.Namelist()

    input_file = input_name + '.nml'
    if not os.path.isfile(input_file) or overwrite:
        if nml_settings is not None:
            input.patch(nml_settings)

        input.write(input_file, force=True)

        if not dry:
            if slurm_script is None:
                if log_file is None:
                    command = './dispatch.x {:}'.format(input_file)
                else:
                    command = './dispatch.x {:} > {:}'.format(input_file, log_file)
            else:
                command = 'sbatch {:} {:}'.format(slurm_script, input_file)

            if verbose:
                print(input_file, ' launched')
            else:
                command = command + ' 2> /dev/null'

            os.environ['OMP_NUM_THREADS'] = str(threads)
            iostat = os.system(command)
            if iostat!=0 and log_file is not None: 
                with open(log_file) as f:
                    for line in f:
                        pass
                    last_line = line
                    print(last_line)
        else:
            print(input_file, ' written')

    else:
        print(input_file, ' already exists')
        input = f90nml.read(input_file)

    datadir = 'data/'
    if 'io_params' in input:
        if 'datadir' in input['io_params']:
            datadir = input['io_params']['datadir']

    return os.path.join(datadir.strip(), input_name)

#=====================================================================
class m3dis(object):
    @classmethod
    #=================================================================
    def read(self, folder, verbose=0, read_atmos=True, read_atom=True, fast_read=False, read_llist=False, **kwargs):
        #-------------------------------------------------------------
        self = m3dis()
        folder = folder + os.sep
        sfolder = folder + 'save' + os.sep
        self.sfolder = sfolder
        self.folder = folder
        self.nml = f90nml.read(folder + 'params.nml')   
        df = pd.read_csv(folder + 'abund.txt', header=None, sep=r'\s+', na_filter=False)
        self.abundances = {k.lower().capitalize():v for k, v in df.values}
        
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
        
        depart = ''
        atom_file = folder + 'atom.txt'
        if read_atom and os.path.isfile(atom_file):
            self.qnorm = self.nml['atom_params']['qnorm']
            loggf_cut  = self.nml['atom_params']['logf_cut']
            red_cut    = self.nml['atom_params']['red_cut']
            blue_cut   = self.nml['atom_params']['blue_cut']
            in_abnd    = self.nml['atom_params']['abundance']
            depart     = self.nml['atom_params']['precomputed_depart'].strip()

            pandas_format = 'pandas_format' in self.nml['atom_params']
            if pandas_format: pandas_format = pandas_format & self.nml['atom_params']['pandas_format']

            nq_fac = self.nml['atom_params']['nq_fac']
            if self.nml['m3d_params']['m1d_legacy_mode']: nq_fac *= 2

            if 'line_mask' in self.nml.keys():
                if 'in_lines' in self.nml['line_mask'].keys():
                    in_lines = self.nml['line_mask']['in_lines']
                    if type(in_lines)==int: in_lines=[in_lines]
                    self.in_lines = np.array(in_lines)-1

            if pandas_format or fast_read:
                if pandas_format:
                    self.atom = pandas_atom.read(atom_file)
                else:
                    self.atom = fast_atom(atom_file, **kwargs)

                self.atom.run = self
                if not hasattr(self, 'in_lines'): self.in_lines = slice(None, self.atom.nline)

                if hasattr(self, 'lam'):
                    self.line = self.atom.trans[self.in_lines]

                    if option(self.nml['atom_params'], 'overlapping_bb_in_output'):
                        self.atom.load_sp(self.lam, qnorm=self.qnorm)
                    else:
                        nnu = self.lam.size
                        ib, ir = self.atom.get_non_overlapping_indices(self.line, nnu, nq_fac)
                        self.atom.bbt['ir'] = 0
                        self.atom.bbt['ib'] = 0
                        self.atom.bbt.loc[self.in_lines, 'ir'] = ir
                        self.atom.bbt.loc[self.in_lines, 'ib'] = ib

                    #mask = bbt.loggf > loggf_cut
                    #mask = mask & (bbt.lam_vac_red > blue_cut) 
                    #mask = mask & (bbt.lam_vac_blue < red_cut) 
                    #lines = lines_class(self.atom, bbt[mask], self)
                    #self.nline = len(lines)
                    #self.line = lines

            else:
                self.atom = atom_class(atom_file, qnorm=self.qnorm, verbose=verbose, **kwargs)
                self.atom.run = self

                if not hasattr(self, 'in_lines'): self.in_lines = slice(None, self.atom.nline)
                lines = np.array(self.atom.trans)[self.in_lines]
                    
                #lines = [line for line in lines if np.log10(line.fosc*line.lower.g) > loggf_cut]
                #lines = [line for line in lines if line.lam_red > blue_cut and line.lam_blue < red_cut]
                self.nline = len(lines)
                self.line = lines

                if hasattr(self, 'lam'):
                    if not option(self.nml['atom_params'], 'overlapping_bb_in_output'):
                        nnu = self.lam.size
                        ib, ir = self.atom.get_non_overlapping_indices(self.line, nnu, nq_fac)
                        self.atom.ib = ib.astype(int)
                        self.atom.ir = ir.astype(int)


            if self.atom.id in self.abundances:
                self.atom.abnd = float(self.abundances[self.atom.id])
            elif in_abnd > -99:
                self.atom.abnd = in_abnd
            elif self.atom.molec:
                moldict = {'CH':'C', 'OH':'O', 'NH':'N', 'CN':'N', 'CO':'C',
                           'NO':'N', 'H2O':'O', 'C2':'C', 'N2':'N', 'O2':'O'}
                self.atom.abnd = self.abundances[moldict[self.atom.id]]

            #self.atom.title = 'A({}) = {:.2f}'.format(self.atom.id.strip(), self.atom.abnd)

            reports = glob.glob(folder + 'convergence_report_*.txt')
            if len(reports)>0:
                self.lvlerr = np.zeros(self.atom.nk)
                for report in reports:
                    df = pd.read_csv(report, sep=r'\s+')
                    self.lvlerr = np.maximum(self.lvlerr, df.maxconv.values)
            
        if (self.nml['m3d_params']['maxiter'] > 0):
            self.mode = 'NLTE'
        elif(depart==''):
            self.mode = 'LTE'
        else:
            self.mode = 'NLTE'

        file = folder + 'angle_quadrature.txt'
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
            if hasattr(self, 'lam'):
                if (verbose > 0): print('Reading Intensities')
                ie = self.read_patch_save('ie', lazy=False)[0]

                self.ie     = ie[:,0,0,:,0]
                self.ie_lte = ie[:,0,0,:,1]
                self.cnt    = ie[:,0,0,:,2]

                self.flux     = np.sum(self.wts * self.ie    , axis=-1)
                self.flux_lte = np.sum(self.wts * self.ie_lte, axis=-1)
                self.flux_cnt = np.sum(self.wts * self.cnt   , axis=-1)

                if self.nml['m3d_params']['save_resolved']:
                    i3 = self.read_patch_save('i3')[0]
                    self.i3     = i3[...,0]
                    self.i3_lte = i3[...,1]
                    self.c3     = i3[...,2]

            if os.path.isfile(folder + 'linelist_out.txt') and read_llist:
                self.llist = pd.read_csv(folder + 'linelist_out.txt', sep=r'\s+', quotechar="'")
                self.llist['lam0'] = vacuum2obs(self.llist.lam_vac)
                self.llist = self.llist.sort_values(by='lam0')

                init = np.frompyfunc(llist_trans, 2, 1)
                self.llist_line = init(np.arange(len(self.llist)), self)

                T = 4000
                elements = [ion.split()[0] for ion in self.llist.ion]
                abnd = np.array([self.abundances[el] if el in self.abundances.keys() else np.nan for el in elements])
                strength = np.log10(self.llist.lam_vac*1e-8 * np.exp(-self.llist.ev_i*EK/T)) + self.llist.loggf + abnd
                self.llist['strength'] = strength

            atmos_params = self.nml['atmos_params']

            if atmos_params['atmos_format'].strip().lower() == 'save':
                mesh = open(atmos_params['atmos_file'].strip() + '/atmos_mesh.txt', 'r')
                afolder = atmos_params['atmos_file'].strip() + '/save/'
            else:
                mesh = open(folder + 'atmos_mesh.txt', 'r')
                afolder = sfolder

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

            if os.path.isfile(afolder + 'atmos_meta.txt') and read_atmos:
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
                self.tau  = compute_tau_scale(self.chi, self.zz*1e8, axis=-1)
                self.ltau = np.log10(self.tau)

                zz = np.tile(self.zz, (nx, ny, 1)).squeeze()
                self.z0 = np.mean(self.get_ltau_surface(0, zz))
                self.height = self.zz - self.z0
                
            if os.path.isfile(folder + 'atmos_data.bin') and read_atmos:
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

        return self

    #=====================================================================
    def get_xx(self, lam, slc=None, dlam=None, xmin=None, xmax=None):
        #-----------------------------------------------------------------
        if slc is None:
            ib = None
            ir = None
        else:
            ib = slc.start
            ir = slc.stop

        if dlam is not None:
            if not isinstance(dlam, (bool)):
                if isinstance(dlam, tuple):
                    dl = dlam[:]
                else:
                    dl = (-dlam, dlam)

                ib = np.searchsorted(lam, dl[0])
                ir = np.searchsorted(lam, dl[1])

        if xmin is not None:
            ib = max(np.searchsorted(lam, xmin)-1, 0)
            if lam[ib] < xmin-10: ib += 1

        if xmax is not None:
            ir = min(np.searchsorted(lam, xmax)+1, lam.size-1)
            if lam[ir-1] > xmax+10: ir -= 1

        mask = slice(ib, ir)
        xx = lam[mask]
        return xx, mask

    #=====================================================================
    def set_xaxis(self, xx, dlam=None, velocity_units=False,
                  set_labels=True, lam_offset=True):
        #-----------------------------------------------------------------
        xmin = min(xx)
        xmax = max(xx)
        dxx = xmax - xmin
        offset = lam_offset and dlam is not None and dlam is not False

        if set_labels:
            ax = plt.gca()
            sub = ax.get_subplotspec()
            grid = sub.get_gridspec()
            if sub.rowspan.stop == grid._nrows:
                if velocity_units: 
                    plt.xlabel(r'$v_{\rm dop}$ [km/s]')
                elif offset:
                    plt.xlabel(r'∆λ [Å]')
                else:
                    plt.xlabel('Wavelength [Å]')

        plt.xlim(xmin, xmax)
        plt.ticklabel_format(useOffset=False)

    #=====================================================================
    def set_yaxis(self, ang, imu, norm, per_aa=False, renorm=None, set_labels=True):
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
            ax = plt.gca()
            sub = ax.get_subplotspec()
            if sub.colspan.start == 0:
                if norm or renorm is not None:
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
    def get_yy(self, ang=None, imu=None, norm=True, cont=False, LTE=False,
               resolved=False, mask=slice(None, None)):
        #-----------------------------------------------------------------
        if not LTE and self.mode=='LTE':
            print('Cannot get NLTE spectrum, because run was in LTE!')
            LTE = True

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
             ang=None, imu=None, norm=True, cont=False, LTE=False, kr=None, 
             resolved=False, **kwargs):
        #-----------------------------------------------------------------
        if kr is not None:
            line = self.atom.trans[kr]
            line.prepare()
            slc = slice(line.ib, line.ir)
            lam0 = line.lam0

        if per_aa:
            norm = False

        oxx = self.lam

        if vshift is not None:
            vshift = vshift / CC_KMS # vshift has to be given in km/s
            lshift = oxx * vshift / (1+vshift)
            oxx = oxx + lshift
    
        if lam0 is not None:
            oxx = oxx - lam0
            if xmin is not None: xmin = xmin - lam0
            if xmax is not None: xmax = xmax - lam0
        else:
            dlam = None
            lam0 = 0

        newx, mask = self.get_xx(oxx, slc=slc, dlam=dlam, xmin=xmin, xmax=xmax)

        vels = np.array([vdop, vrot, zeta_rt])
        vmask = vels != None
        if any(vmask): # we need to add padding to avoid convolution artefacts at the edges
            uconv = 5*(newx[newx.size//2] + lam0)/CC_KMS
            pad = uconv * np.sum(vels[vmask])

            if dlam is not None: dlam = dlam + pad
            if xmin is not None: xmin = xmin - pad
            if xmax is not None: xmax = xmax + pad

        xx, mask = self.get_xx(oxx, slc=slc, dlam=dlam, xmin=xmin, xmax=xmax)

        get_yy_dict = {'LTE':LTE, 'norm':norm, 'cont':cont, 'resolved':resolved, 'mask':mask}
        if hasattr(ang, '__iter__'):
            yy = [self.get_yy(ang=i, **get_yy_dict) for i in ang]
        elif hasattr(imu, '__iter__'):
            yy = [self.get_yy(imu=i, **get_yy_dict) for i in imu]
        else:
            yy = [self.get_yy(ang=ang, imu=imu, **get_yy_dict)]

        yy = np.array(yy)

        if renorm is not None:
            if renorm == 0:
                yy = yy / np.max(yy)
            elif renorm > 0:
                yy = yy / continuum(xx, yy)
            elif renorm < 0:
                xnorm, mnorm = self.get_xx(oxx, dlam=-renorm)
                ynorm = self.get_yy(ang=ang, imu=imu, norm=norm, cont=cont, LTE=LTE, mask=mnorm)
                yy = yy / np.max(ynorm)

        xx = xx + lam0
        newx = newx + lam0

        if set([vdop, vrot, zeta_rt, convolve]) != {None}:
            newy = []
            for y in yy:
                sxx, syy = convol(xx, y, vdop=vdop, vrot=vrot, zeta_rt=zeta_rt, custom=convolve)
                newy.append(np.interp(newx, sxx, syy))

            yy = np.array(newy)

        xx = newx

        if lam0 is not None:
            offset = lam_offset and dlam is not None and dlam is not False
            if offset:
                xx = xx - lam0

            if velocity_units:
                xx = xx / lam0
                if dlam is None or dlam is False:
                    xx = xx - 1
                xx = xx * CC_KMS

        if per_aa:
            yy = yy * CC / (xx**2) * 1e8

        if fac is not None:
            yy = yy * fac
        
        if tilt is not None:
            yy = yy * (1 + tilt*(xx - (xx[0]+xx[-1])/2))

        if len(yy)==1:
            yy = yy[0]

        return xx, yy
    
    #=====================================================================
    def interp(self, newx, **kwargs):
        #-----------------------------------------------------------------
        xmin = newx.min()
        xmax = newx.max()

        #if 'vdop' in kwargs.keys():
        #    if kwargs['vdop'] is not None:
        #        xmin = xmin-kwargs['vdop']
        #        xmax = xmax+kwargs['vdop']
        #if 'vrot' in kwargs.keys():
        #    if kwargs['vrot'] is not None:
        #        xmin = xmin-kwargs['vrot']
        #        xmax = xmax+kwargs['vrot']

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
    def plot_spectrum(self, velocity_units=None, set_labels=True,
                      lam_offset=None, ang=None, imu=None, norm=True,
                      dlam=None, renorm=None, per_aa=False, show_lines=False,
                      label=None, **kwargs):
        #-----------------------------------------------------------------
        crop_kwargs = {k:v for k, v in kwargs.items() if k not in plt_kwargs}
        plot_kwargs = {k:v for k, v in kwargs.items() if k not in crop_kwargs}
        
        xx, yy = self.crop(dlam=dlam, ang=ang, imu=imu, norm=norm, renorm=renorm,
                           per_aa=per_aa, velocity_units=velocity_units, lam_offset=lam_offset,
                           **crop_kwargs)
        
        self.set_xaxis(xx, dlam, velocity_units=velocity_units, set_labels=set_labels,
                       lam_offset=lam_offset)
        
        self.set_yaxis(ang, imu, norm, renorm=renorm, per_aa=per_aa)

        if show_lines:
            self.plot_lline_labels()

        if type(label) is str:
            if '{WEQ}' in label:
                weq = np.trapz(1-yy, x=xx)*1e3
                label = label.replace('{WEQ}', '{:.2f}'.format(weq))

        return plt.plot(xx, yy, label=label, **plot_kwargs)

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
    def plot_surface(self, Z, ang=0, rasterized=True, shading='nearest', 
                     show_labels=True, cmap='hot', rot90=0, **kwargs):
        #-----------------------------------------------------------------
        xx = self.xx + self.dx/2
        yy = self.yy + self.dy/2
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
        xmy = np.array([10**(self.abundances[key]-12)*atomic_weigths[key]
                        for key in self.abundances.keys()])
        grph = np.sum(xmy) * UU
        toth = self.rho / grph
        return toth

    #=====================================================================
    def get_pop(self, ilvl=None):
        #-----------------------------------------------------------------
        if ilvl is not None:
            if ilvl < 0: ilvl = self.atom.nk + ilvl
        mem = self.read_patch_save('atom_patch', zdim=2)[0]
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
    def get_depart(self, ilvl=None, mean=False):
        #-----------------------------------------------------------------
        if ilvl is not None:
            if ilvl < 0: ilvl = self.atom.nk + ilvl

        if hasattr(self, 'runs'):
            xy = np.array([run.get_depart(ilvl, mean=mean) for run in self.runs])

            if mean:
                ltau = np.mean(xy[:,0], axis=0)
                yy1d = [np.interp(ltau, xy[i,0], xy[i,1]) for i in range(len(xy))]
                yy1d = np.mean(yy1d, axis=0)
                return ltau, yy1d
            else:
                return xy[:,0], xy[:,1]

        mem = self.read_patch_save('atom_patch', zdim=2)[0]
        if mem.ndim==6: mem = mem[...,0]
        dep = mem[...,1]

        if ilvl is not None:
            dep = dep[...,ilvl]

        dep = np.array(dep)
        ltau = self.ltau

        if self.dim=='1D':
            dep = dep.squeeze()
        elif mean:
            ltau, dep = self.get_ltau_mean(dep, zaxis=2)

        return ltau, dep
    
    #=====================================================================
    def plot_depart(self, ilvl, mean=True, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.get_depart(ilvl, mean=mean)
        plt.title('Departure coefficient of level {:}'.format(ilvl))
        plt.ylabel(r'$\rm n\,/\,n_{\rm LTE}$')
        plt.xlabel(r'log$\tau_{500}$')
        return plt.plot(xx, yy, **kwargs)
    
    
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
        if folder is None: folder = self.sfolder

        try:
            meta = open(folder + kind + '_meta.txt')
        except FileNotFoundError:
            print('{:} not found!'.format(folder + kind + '_meta.txt'))
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
                dct[tuple(key)] = np.memmap(folder + file, dtype=dtype, mode='r', shape=dims, order='F')

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
        step = -1 if ltau[0,0,0] > ltau[0,0,-1] else 1
        xx = ltau[...,::step].copy()
        yy = cube[...,::step].copy()
        if log: yy = np.log(yy)
        if newx is None: newx = np.mean(xx, axis=(0,1))
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
    def get_ltau_mean(self, cube, log=False, **kwargs):
        #-----------------------------------------------------------------
        if hasattr(self, 'runs'):
            newx = np.mean([run.ltau for run in self.runs], axis=(0,1,2))
            kwargs['newx'] = newx
            yy1d = [run.get_ltau_mean(cube, log=log, **kwargs)[1] for run in self.runs]
            yy1d = np.mean(yy1d, axis=0)
            return newx, yy1d

        if type(cube) is str:
            cube = getattr(self, cube)

        if self.dim == '3D':
            ltau1d, yy3d = self.interp_ltau_surf(cube, log=log, **kwargs)
            if log: yy3d = np.log(yy3d)
            yy1d = np.mean(yy3d, axis=(0,1))
            if log: yy1d = np.exp(yy1d)
        else:
            ltau1d = self.ltau
            yy1d = cube

        return ltau1d, yy1d
    
    #=====================================================================
    def plot_ltau_mean(self, yy, log=False, xmin=None, xmax=None,
                       ymin=None, ymax=None, **kwargs):
        #-----------------------------------------------------------------
        xx, yy = self.get_ltau_mean(yy, log=log)
        xx, yy = crop_inputs(xx, yy, xmin, xmax, ymin, ymax)
        return plt.plot(xx, yy, **kwargs)
    
    #=====================================================================
    def get_geo_mean(self, yy, log=False):
        if hasattr(self, 'runs'):
            yy1d = [run.get_geo_mean(yy, log=log)[1] for run in self.runs]
            yy1d = np.mean(yy1d, axis=0)
            newx = self.runs[0].zz - np.mean([run.z0 for run in self.runs])
            return newx, yy1d
        
        if type(yy) is str:
            yy = getattr(self, yy)

        if log: yy = np.log(yy)
        yy = np.mean(yy, axis=(0,1))
        if log: yy = np.exp(yy)
        xx = self.height
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
    def plot_line_labels(self, xx=None, labels=None, weights=None, offset=None, **kwargs):
        #-----------------------------------------------------------------
        if xx is None: xx = np.array([line.lam0 for line in self.line])
        if offset is not None: xx = xx + offset
        if labels is None: labels = np.array([self.atom.id for line in self.line])
        if weights is None:
            T = 4000
            loggf = np.array([line.loggf for line in self.line])
            lam0s = np.array([line.lam0 for line in self.line])*1e-8
            evis = np.array([line.evi for line in self.line])
            weights = np.log10(lam0s * np.exp(-evis*EK/T)) + loggf

        return self.plot_line_annotations(xx, labels, weights=weights, **kwargs)
    
    #=====================================================================
    def plot_lline_labels(self, xx=None, labels=None, weights=None, offset=None, **kwargs):
        #-----------------------------------------------------------------
        if xx is None: xx = self.llist.lam0.values
        if offset is not None: xx = xx + offset
        if labels is None: labels = self.llist.ion.values
        if weights is None: weights = self.llist.strength.values
        return self.plot_line_annotations(xx, labels, weights=weights, **kwargs)
    
    #=====================================================================
    def plot_line_annotations(self, xx, labels, weights=None, ha='center', va='bottom',
                         rotation=90, fontfamily='serif', color='k', 
                         linewidth=0.8, ntop=-1, **kwargs):
        #-----------------------------------------------------------------
        ax = plt.gca()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        h = 0.1 * (ymax - ymin)

        mask = (xx > xmin) & (xx < xmax)
        xx = xx[mask]
        labels = np.array(labels)[mask]
        if weights is not None: weights = weights[mask]

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
            mask = np.argsort(weights)[::-1][:ntop]
            mask = np.isin(np.arange(xx.size), mask)
            #mask = llist.nlargest(ntop, 'strength').sort_values('lam0').index
            #llist = llist.loc[mask]
            labels = labels[mask]
            xx = xx[mask]

        #xx = llist.lam0.values
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
    
    #=================================================================
    def fit_observed(self, oxx, oyy, fit=[], bounds=None, p0=None, sigma=None, method='curve_fit', tol=None, **kwargs):
        lbounds = {'vshift':-3, 'vdop':1e-3, 'vrot':1e-3, 'zeta_rt':1e-3, 'fac':0.5, 'tilt':-0.2}
        ubounds = {'vshift':+3, 'vdop':1e1,  'vrot':1e1,  'zeta_rt':1e1,  'fac':2.0, 'tilt':+0.2}
        p0s     = {'vshift':0,  'vdop':1,    'vrot':1,    'zeta_rt':1,    'fac':1.0, 'tilt':0.0}
        if bounds is None: bounds = ([lbounds[f] for f in fit], [ubounds[f] for f in fit])
        if p0 is None: p0 = [p0s[f] for f in fit]

        def opt_func(newx, *args):
            ka = {f:args[i] for i, f in enumerate(fit)}
            return self.interp(newx, **ka, **kwargs)
        
        if method == 'curve_fit':
            pfit, _ = optimize.curve_fit(opt_func, oxx, oyy, p0=p0, bounds=bounds, sigma=sigma)

        else:
            bounds = np.vstack(bounds).T
            def opt_func2(params, newx):
                res = oyy - opt_func(newx, *params)
                if sigma is not None: res = res/sigma
                return np.sum(res**2)
            pfit = optimize.minimize(opt_func2, p0, args=(oxx), bounds=bounds, method=method, tol=tol).x

        fitted = opt_func(oxx, *pfit)
        return {f:pfit[i] for i, f in enumerate(fit)}, fitted
    
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
    def to_dispatch1d_atmos(self, dir, name, vmic=1.0):
        ltau, chi = self.get_ltau_mean('chi',  log=True)
        _, temp   = self.get_ltau_mean('temp', log=True)
        _, ne     = self.get_ltau_mean('ne',   log=True)
        _, rho    = self.get_ltau_mean('rho',  log=True)
        tau = 10**(ltau)

        # Flip the arrays to get them in increasing tau order
        tau = tau[::-1]
        chi = chi[::-1]

        # Calculate dz = - d_tau / chi at the midpoint between each pair of tau values
        dz = - np.diff(tau) * 2 / (chi[:-1] + chi[1:])

        # Compute z by cumulatively summing dz
        zz = np.concatenate(([0], np.cumsum(dz)))  # Start from z0 = 0 or any other initial condition
        zz = - zz[::-1]
        zz = zz - np.interp(1, tau, zz)

        df = pd.DataFrame()
        df['Depth'] = zz
        df['Temp'] = temp
        df['Pe'] = ne*temp*BK
        df['Rho'] = rho
        df['Vturb'] = vmic

        f = open(dir + f'{name:}.txt', 'w+')
        f.write(f'{name:}\n')
        f.write(f'{len(df):}\n')
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

    @classmethod
    #=================================================================
    def combine(self, runs):
        #-------------------------------------------------------------
        self = m3dis()
        self.run0 = runs[0]
        self.folder = '\n'.join([run.folder for run in runs])

        def check_and_combine(key):
            if hasattr(self.run0, key):
                atts = [getattr(run, key) for run in runs]
                same = all([np.array_equal(a, atts[0]) for a in atts])
                if not same:
                    print('Warning: {:} is not the same for all runs!'.format(key))
                else:
                    setattr(self, key, atts[0])

        self.runs = runs
        self.nml = self.run0.nml
        check_and_combine('mode')
        check_and_combine('mus')
        check_and_combine('imus')
        check_and_combine('nmu')
        check_and_combine('n_ang')
        check_and_combine('vec')
        check_and_combine('wts')
        check_and_combine('nline')

        self.ne = np.array([run.ne for run in runs])
        self.rho = np.array([run.rho for run in runs])
        self.temp = np.array([run.temp for run in runs])
        self.vx = np.array([run.vx for run in runs])
        self.vy = np.array([run.vy for run in runs])
        self.vz = np.array([run.vz for run in runs])
        self.chi = np.array([run.chi for run in runs])
        self.tau = np.array([run.tau for run in runs])
        self.ltau = np.array([run.ltau for run in runs])
        
        if hasattr(self.run0, 'atom'):
            self.atom = copy.copy(self.run0.atom)
            self.atom.run = self
            self.in_lines = self.run0.in_lines

            if type(self.atom.trans[0]) is fast_trans:
                init = np.frompyfunc(fast_trans, 2, 1)
                self.atom.trans = init(np.arange(self.atom.nline), self.atom)
            else:
                self.atom.trans = [copy.copy(trans) for trans in self.run0.atom.trans]

            self.line = np.array(self.atom.trans)[self.in_lines]

        if hasattr(self.run0, 'lam'):
            self.lam    = self.run0.lam
            self.ie     = np.mean([run.ie     for run in runs], axis=0)
            self.ie_lte = np.mean([run.ie_lte for run in runs], axis=0)
            self.cnt    = np.mean([run.cnt    for run in runs], axis=0)

            self.flux     = np.mean([run.flux     for run in runs], axis=0)
            self.flux_lte = np.mean([run.flux_lte for run in runs], axis=0)
            self.flux_cnt = np.mean([run.flux_cnt for run in runs], axis=0)

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