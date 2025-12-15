from .atom import AtomClass
from .auxfuncs import *
from .decorators import broadcast_arguments
from .plotting_scripts import *
from .fast_atom import FastAtom, PandasAtom
from .linelist import LineList
from .NameList import NameList
from .spectrum import SpectrumClass
from .spectrum_grid import SpectrumGrid
from .memmap_list import MemmapList
from .atmos import AtmosClass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
import textwrap
import pathlib
import functools
from tqdm import tqdm
import multiprocessing as mp

def option(nml, key):
    if key in nml.keys():
        return nml[key]
    else:
        return False

#=====================================================================
def read_run(folder):
  return tumult.read(folder, read_atom=False)

#=====================================================================
def read_list(folders):
  with mp.Pool() as pool:
    run_list = list(tqdm(pool.imap(read_run, folders), total=len(folders)))
    return np.array(run_list)

#=====================================================================
class grid(SpectrumGrid):
  def __init__(self, runs, *args, **kwargs):
    if type(runs[0]) == str:
      runs = read_list(runs)
    
    return super().__init__(runs, *args, **kwargs)

#=====================================================================
class tumult(object):
    @property
    def flux(self):
        return self.spectrum[[('lam', '', ''), ('NLTE', 'Flux', '')]]

    @property
    def atmos(self):
        if not hasattr(self, '_atmos'):
            self._atmos = AtmosClass.from_run(self)
        return self._atmos
    
    @classmethod
    #=================================================================
    def read(cls, folder, verbose=0, read_atom=True, fast_read=False, read_llist=False, **kwargs):
        self = cls() # create a new instance of the class
        folder = pathlib.Path(folder)
        sfolder = folder / 'save'
        self.sfolder = sfolder
        self.folder = folder
        self.nml = NameList.from_file(folder / 'params.nml')
        self.datadir = pathlib.Path(self.nml['io_params']['datadir'].strip())
        
        df = pd.read_csv(folder / 'abund.txt', header=None, sep=r'\s+', na_filter=False)
        self.abundances = {k.lower().capitalize():v for k, v in df.values}

        def read_wavelength_file(file):
            if (self.sfolder / file).exists():
                return vac2air(np.fromfile(sfolder / file, dtype=np.float64))[::-1]
        
        self.lam = read_wavelength_file('out_lam.bin')
        self.bb_lam = read_wavelength_file('bb_lam.bin')
        self.bf_lam = read_wavelength_file('bf_lam.bin')
        self.all_lam = read_wavelength_file('all_lam.bin')
        self.has_spectrum = self.lam is not None and self.lam.size > 0
        self.legacy = not 'read_departures' in self.nml['atmos_params']
        
        if self.has_spectrum:
            self.has_spectrum = self.nml['m3d_params']['long_scheme'].strip() != 'none'

        self.read_atom(read_atom, fast_read, verbose, **kwargs)

        file = folder / 'angle_quadrature.txt'
        if file.exists():
            df = pd.read_csv(file, header=None, skiprows=1, sep=r'\s+', names=['ang', 'wmu', 'mux', 'muy', 'muz'])
            self.vec = df.values[:,2:]
            self.wts = df.values[:,1] * 2*np.pi
            self.n_ang = self.wts.size
            self.mus, self.imus = np.unique(-self.vec[:,-1], return_inverse=True)
            self.mus = -self.mus
            self.nmu = self.mus.size
        
        if self.nml.get_item('m3d_params', 'make_eos', False):
            mem = self.read_patch_save('eos', fdim=None, lazy=False)
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
                
                if self.lam.shape[0] != ie.shape[0]:
                    print(f'Wavelength and Intensity dimensions do not match!\nRun {folder} might have never completed!')

                if self.legacy:
                    self.ie     = ie[:,0,0,:,0]
                    self.ie_lte = ie[:,0,0,:,1]
                    self.cnt    = ie[:,0,0,:,2]
                    self.mode = 'LTE' if np.all(self.ie_lte == self.ie) else 'NLTE'
                else:
                    self.ie     = ie[:,0,0,:,0]
                    self.cnt    = ie[:,0,0,:,1]
                    if saved_LTE:
                        self.ie_lte = ie[:,0,0,:,2]
                    else:
                        self.ie_lte = self.ie

                    self.mode = 'NLTE' if saved_LTE else 'LTE'

                mus = [f'µ={mu:.2f}' for mu in self.vec[:,-1]]

                flux = {}
                flux[('NLTE','Flux', '')] = np.sum(self.wts * self.ie    , axis=-1)
                flux[('LTE' ,'Flux', '')] = np.sum(self.wts * self.ie_lte, axis=-1)
                flux[('Cont','Flux', '')] = np.sum(self.wts * self.cnt   , axis=-1)
                ints = {}
                ints.update({('NLTE', mu, ang): self.ie[:,ang] for ang, mu in enumerate(mus)})
                ints.update({('LTE' , mu, ang): self.ie_lte[:,ang] for ang, mu in enumerate(mus)})
                ints.update({('Cont', mu, ang): self.cnt[:,ang] for ang, mu in enumerate(mus)})

                #columns = pd.MultiIndex.from_tuples(data.keys())
                self.spectrum = SpectrumClass.from_arrays(wave=self.lam, flux=flux, ints=ints,
                                                   wave_unit=u.angstrom, flux_unit=u.erg/(u.s*u.cm**2*u.Hz))
                #self.spectrum = self.spectrum.sort_index(axis=1, ascending=[False, False, True])

                if self.nml['m3d_params']['save_resolved']:
                    i3 = self.read_patch_save('i3')[0]
                    self.i3     = i3[...,0]
                    self.c3     = i3[...,1]
                    if saved_LTE:
                        self.i3_lte = i3[...,2]

            if os.path.isfile(folder / 'linelist_out.txt') and read_llist:
                self.read_llist()
                
        return self

    #=====================================================================
    def read_atom(self, read_atom=True, fast_read=True, verbose=True, **kwargs):
        self.atoms = []
        atom_files = list(self.folder.glob('atoms/*'))
        atom_files += list(self.folder.glob('atom.txt'))
        moldict = {'CH':'C', 'OH':'O', 'NH':'N', 'CN':'N', 'CO':'C',
                    'NO':'N', 'H2O':'O', 'C2':'C', 'N2':'N', 'O2':'O'}

        if 'atom_params' in self.nml:
            self.qnorm = self.nml['atom_params']['qnorm']
            nq_fac = self.nml['atom_params']['nq_fac']
            if self.nml['m3d_params']['m1d_legacy_mode']: nq_fac *= 2
            in_abnd = self.nml['atom_params']['abundance']
            in_abnd = [in_abnd] if not hasattr(in_abnd, '__iter__') else in_abnd
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
                    pf = self.nml['atom_params'].get('pandas_format', None)
                    pandas_format = not atom_file.name.startswith('atom.') if pf is None else pf

                    #if 'line_mask' in self.nml.keys():
                    #    if 'in_lines' in self.nml['line_mask'].keys():
                    #        in_lines = self.nml['line_mask']['in_lines']
                    #        if type(in_lines)==int: in_lines=[in_lines]
                    #        self.in_lines = np.array(in_lines)-1

                    if pandas_format or fast_read:
                        if pandas_format:
                            atom = PandasAtom.read(atom_file)
                        else:
                            atom = FastAtom(atom_file, **kwargs)

                        atom.run = self
                        if not hasattr(atom, 'in_lines'): atom.in_lines = slice(None, None)

                        if self.has_spectrum:
                            atom.line = atom.trans.get_lines(atom.in_lines)

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
                        atom = AtomClass(atom_file, qnorm=self.qnorm, verbose=verbose, **kwargs)
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
    def read_llist(self):
        self.llist = LineList.from_dispatch_format(self.folder / 'linelist_out.txt')

    #=====================================================================
    def get_line(self, lam0, tol=0.05):
        #-----------------------------------------------------------------
        if not hasattr(self.atoms[0], 'line'): self.read_atom()
        all_lines = np.concatenate([[line for line in atom.line] for atom in self.atoms])
        all_lam0 = np.array([line.lam0 for line in all_lines])
        i = np.argmin(abs(all_lam0-lam0))
        line = all_lines[i]
        if abs(line.lam0 - lam0) > tol:
            print(f'WARNING: Line not found! Closest line to {lam0:.2f}Å found at {line.lam0:.2f}Å')

        return line
    
    #=====================================================================
    @broadcast_arguments(numeric=False)
    def crop(self, slc=None, lam0=None, ang=None, imu=None, norm=True, 
             cont=False, LTE=False, line=None, resolved=False, **kwargs):
        #-----------------------------------------------------------------
        if line is not None:
            line.prepare()
            slc = slice(line.ib, line.ir)
            lam0 = line.lam0

        if not LTE and self.mode=='LTE':
            print('Cannot get NLTE spectrum, because run was in LTE!')
            LTE = True

        if imu is not None:
            if imu >= self.nmu:
                raise ValueError('imu is out of range. Nmu = {:}'.format(self.nmu))
            elif imu < 0:
                imu = self.imus.max() + imu + 1

        if self.n_ang == 1: ang=0

        if cont: key = ('Cont',)
        elif LTE: key = ('LTE',)
        else: key = ('NLTE',)

        if ang is None and imu is None:
            key += ('Flux', '')
            cnt = self.spectrum[('Cont', 'Flux', '')]
            yy = self.spectrum[key]

        else:
            #if resolved:
            #    if cont:
            #        yy = self.c3
            #    elif LTE:
            #        yy = self.i3_lte
            #    else:
            #        yy = self.i3

            #    if ang is None:
            #        raise ValueError('ang is required for resolved spectrum')
            #    else:
            #        cnt = np.array(cnt[..., ang]).T
            #        yy = np.array(yy[..., ang]).T

            if ang is not None:
                mu = f'µ={self.vec[ang,-1]:.2f}'
                key += (mu, ang)
                yy = self.spectrum[key]
                cnt = self.spectrum[('Cont', mu, ang)]
            elif imu is not None:
                mu = f'µ={self.mus[imu]:.2f}'
                cnt = self.get_imu_mean(('Cont', mu))
                yy = self.get_imu_mean(key + (mu,))
                key += (mu, '')
            else:
                cnt = self.spectrum[('Cont', 'Flux', '')]
                
        if norm:
            yy = yy / cnt
            yy.data.columns = [u.dimensionless_unscaled]

        spec = yy.crop(slc=slc, lam0=lam0, **kwargs)
        return spec
    
    #=====================================================================
    #@functools.cache
    def get_imu_mean(self, key):
        #-----------------------------------------------------------------
        spec = self.spectrum[key]
        columns = spec.data.columns
        unit = columns[0][-1]
        spec.data = pd.DataFrame({unit:spec.data.mean(axis=1)})
        return spec

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
        self.plot_cube_slice(Z, ang=ang, **kwargs)

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
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            if ang is not None:
                plt.title(f'({mux:.2f},{muy:.2f},{muz:.2f})')
        else:
            plt.xticks([])
            plt.yticks([])

    #=====================================================================
    def get_pop(self, atom=0, ilvl=None):
        #-----------------------------------------------------------------
        if type(atom) is int: atom = self.atoms[atom]
        if ilvl is not None:
            if ilvl < 0: ilvl = atom.nk + ilvl
        mem = self.read_patch_save(f'pop_{atom.id}', zdim=2)
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
            
        mem = self.read_patch_save(f'pop_{atomid}', zdim=2)

        if mem.ndim==6: mem = mem[...,0]
        dep = mem[...,1]

        if ilvl is not None:
            dep = dep[...,ilvl]

        dep = np.array(dep)

        if self.dim=='3D' and mean:
            ltau, dep = self.atmos.get_ltau_mean(dep, zaxis=2)
        else:
            ltau = self.atmos.ltau[...,::-1].squeeze()
            dep = dep[:,:,::-1].squeeze()

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
    def read_patch_save(self, kind, folder=None, concat=True, fdim=0, zdim=None, lazy=True):
        #-----------------------------------------------------------------
        folder = pathlib.Path(folder) if folder else self.sfolder
        resolved_folder = str(folder.resolve()) + os.sep
        meta_file = resolved_folder + f'{kind.replace(' ', '_')}_meta.txt'
        dct = {}

        try:
            with open(meta_file) as meta:
                for line in meta:
                    if line.strip().startswith('*'): continue
                    info = line.split()
                    file = info[0]
                    dtype = f'<f{int(info[1]) // 8}'
                    dims  = tuple(map(int, info[2:]))
                    key = file.split('.')[0].split('_')[-1]
                    key = tuple(map(int, textwrap.wrap(key, 3)))
                    dct[key] = np.memmap(resolved_folder + file, dtype=dtype, mode='r', shape=dims, order='F')

        except FileNotFoundError:
            if (not kind in self.nml['m3d_params']['save_patch_kind'].split()):
                print('The patch was not saved. Check the save_patch_kind parameter in the m3d_params namelist.')
            raise FileNotFoundError(f'{meta_file} not found!')

        if concat:
            zdim = zdim or len(dims)-1
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
        concat = MemmapList if lazy else np.concatenate

        def merge_set(inus):
            merged = []
            for iz in all_iz:
                parts = [dct[(iz, inu)] for inu in inus if (iz, inu) in dct]
                if not parts: continue
                mem = parts[0] if len(parts) == 1 else concat(parts, axis=fdim)
                if not f_indep: mem = np.flip(mem, axis=fdim)
                merged.append(mem)
            return merged[0] if len(merged) == 1 else concat(merged, axis=zdim)

        if f_indep:
            return merge_set(all_inu)
        else:
            set1 = [inu for inu in all_inu if inu > 0 and inu <= n_nu]
            set2 = [inu for inu in all_inu if inu > 0 and inu > n_nu]
            set3 = [inu for inu in all_inu if inu < 0][::-1]
            sets = [s for s in [set1, set2, set3] if len(s)>0]
        
            return [merge_set(s) for s in sets]
    
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
                     syst=None, method='L-BFGS-B', tol=None, res='auto', **kwargs):
        #-----------------------------------------------------------------
        if len(oxx) == 0:
            raise ValueError('empty oxx array in fit_observed!')
        if len(oxx) != len(oyy):
            raise ValueError('oxx and oyy must have the same length!')

        fit = fit or []
        bounds = bounds or {}
        p0 = p0 or {}
        
        if len(fit)==0: return {}, None
        lbounds = {'vshift':-3, 'vdop':1e-3, 'vrot':1e-3, 'zeta_rt':1e-3, 'continuum':0.5, 'tilt':-0.2}
        ubounds = {'vshift':+3, 'vdop':1e1,  'vrot':1e1,  'zeta_rt':1e1,  'continuum':2.0, 'tilt':+0.2}
        p0s     = {'vshift':0,  'vdop':1,    'vrot':1,    'zeta_rt':1,    'continuum':1.0, 'tilt':0.0}
        
        p0s.update(p0)
        for key in bounds.keys():
            lbounds[key], ubounds[key] = bounds[key]

        p0 = [p0s[f] for f in fit]
        bounds = ([lbounds[f] for f in fit], [ubounds[f] for f in fit])
        bounds = np.vstack(bounds).T

        if res == 'auto':
            res = np.mean(oxx) / np.mean(np.diff(oxx))

        def evaluate(newx, *args):
            ka = {f:args[i] for i, f in enumerate(fit)}
            newy = self.crop(newx=newx, **ka, res=res, **kwargs).ydata
            if len(newy) == 0:
                raise ValueError('Observations outside of model wavelength range!')
            return newy
        
        def fit_params(x, y):
            def residuals(params, newx):
                r = y - evaluate(newx, *params)
                if sigma is not None: r = r/sigma
                return np.sum(r**2)
            
            pfit = optimize.minimize(residuals, p0, args=(x), bounds=bounds, method=method, tol=tol).x
            return pfit

        pfit = fit_params(oxx, oyy)
        fitted = evaluate(oxx, *pfit)
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
        
    #=====================================================================
    def get_slice(self, line=None, dlam=None, lam0=None, xmin=None, xmax=None):
        #-----------------------------------------------------------------
        slc = slice(None, None)

        if line is not None:
            line.prepare()
            slc = slice(line.ib, line.ir)
            lam0 = line.lam0

        if dlam is not None and lam0 is not None:
            xmin = lam0 - dlam
            xmax = lam0 + dlam

        if xmin is not None:
            ib = np.searchsorted(self.lam, xmin)
            slc = slice(ib, slc.stop)

        if xmax is not None:
            ir = np.searchsorted(self.lam, xmax)
            slc = slice(slc.start, ir)

        return slc

    #=====================================================================
    def get_cntrbf(self, ang=None, imu=None, fmean=False, gmean=False,
                   slc=None, norm=False, **kwargs):
        #-----------------------------------------------------------------
        mem  = self.read_patch_save('cntrbf', fdim=0, zdim=3, lazy=True)[-1]
        slc = slc if slc is not None else self.get_slice(**kwargs)

        lam = self.lam[slc]
        cntrbf = mem[slc]

        #if LTE:
        #    cntrbf = mem[..., 2:4, :]
        #else:
        #    cntrbf = mem[..., 0:2, :]
        
        #nf, nx, ny, nz, nset, nang = cntrbf.shape
        ltau = self.atmos.ltau
        get_mean = gmean and self.atmos.dim == '3D'

        if get_mean:
            new_ltau = np.mean(ltau, axis=(0,1))[::-1]

        all_angs = np.arange(self.n_ang)
        if ang is None and imu is None:
            angs = all_angs
        elif imu is not None:
            if imu < 0: imu = self.imus.max()+imu+1
            angs = all_angs[self.imus==imu]
        elif ang is not None:
            angs = np.array(ang).squeeze()
            if angs.ndim == 0: angs = np.array([angs])

        wts = self.wts[angs]
        wts = wts / np.sum(wts)
            
        conversion = np.log(10) * self.atmos.tau / self.atmos.chi

        first_time = True
        for ang, w in zip(angs, wts):
            cang = np.array(cntrbf[..., ang])
            cf, net_chi = np.moveaxis(cang, -1, 0)
            vec = self.vec[ang]
            zz = self.atmos.zz / vec[-1] * 1e8
            tau = compute_tau_scale(net_chi, zz, axis=-1)
            cf = conversion * cf * np.exp(-tau)
            if self.atmos.dim == '3D':
                vec = vec * np.array([-1, -1, 1])
                cf = rotate_array(self.atmos.xx, self.atmos.yy, self.atmos.zz, cf, vec)

            if get_mean:
                yy = np.array(cf[...,::-1]).copy()
                xx = ltau[...,::-1].copy()
                newcf = interp_tau(new_ltau, xx, yy)
                cf = np.mean(newcf, axis=(-3,-2))

            if first_time:
                cfunc = np.zeros_like(cf)
                first_time = False

            cfunc = cfunc + w*cf

        if get_mean: ltau = new_ltau
        cfunc = cfunc.squeeze()
        if cfunc.ndim>1 and fmean: cfunc = np.trapz(cfunc, x=lam, axis=0)
        if norm:
            normalisation = abs(np.trapz(cfunc, x=ltau))
            cfunc = cfunc / normalisation
        return ltau, cfunc.squeeze()

    #=====================================================================
    def plot_cntrbf(self, line=None, lam0=None, dlam=None, xmin=None, xmax=None,
                    slc=None, imu=None, ang=None, title=None, norm=False, **kwargs):
        #-----------------------------------------------------------------
        if slc is None:
            slc = self.get_slice(line=line, dlam=dlam, lam0=lam0, xmin=xmin, xmax=xmax)

        lam = self.lam[slc]
        ltau, cf = self.get_cntrbf(slc=slc, imu=imu, ang=ang, gmean=True, fmean=True, norm=norm)

        if type(slc) == int:
            wave = rf'λ={lam:.3f}[Å]'
        else:
            wave = rf'λ={lam[0]:.2f}-{lam[-1]:.2f}[Å]'
            
        if ang is None and imu is None:
            angle = 'Flux'
        else:
            if ang is not None:
                imu = self.imus[ang]

            angle = r'μ={:.2f}'.format(self.mus[imu])

        plt.xlabel(r'log$\tau_{500}$')
        plt.ylabel('Contribution function')
        set_title(title, default=wave + ', ' + angle)
        return plt.plot(ltau, cf, **kwargs)
    

    #=====================================================================
    def hist_cntrbf(self, imu=None, ang=None, LTE=False, dlam=0.5, title=None, cmap='binary', **kwargs):
        #-----------------------------------------------------------------
        xx, slc = self.get_xx(dlam=dlam)
        slc = slice(slc.start-self.ib, slc.stop-self.ib)

        ltau, cf = self.get_cntrbf(imu=imu, ang=ang, LTE=LTE, mean=True, fslice=slc)
            
        self.set_xaxis(xx, dlam=dlam)
        plt.ylabel(r'log$\tau_{500}$')
        set_title(title, default='Contribution function of ' + self.title)

        X, Y, Z = *np.meshgrid(xx, ltau), cf.T
        return plt.pcolormesh(X, Y, Z, cmap=cmap, edgecolors='face', **kwargs)

    
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
    def combine(cls, runs, weights=None):
        #-------------------------------------------------------------
        self = cls() # create a new instance of the class
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
        check_and_combine('abundances')

        if hasattr(runs[0], 'atmos'):
            self.atmos = AtmosClass.combine([run.atmos for run in runs], weights=weights)

        weights = weights or np.ones(len(runs))/len(runs)

        if self.run0.has_spectrum:
            self.spectrum = SpectrumClass.combine([run.spectrum for run in runs], weights=weights)

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