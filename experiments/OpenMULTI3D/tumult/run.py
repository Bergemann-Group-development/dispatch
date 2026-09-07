"""The dispatch run class (`tumult`) and its helpers.

Moved from tumult/tumult.py in the Phase 6 decomposition (that module remains as the
package's import shim). The raw patch-save reading lives in :mod:`tumult.io_binary`
and the intensity-layout / flux physics in :mod:`tumult.flux`; this class is the
orchestration and analysis facade over them (lazy atmos/linelist/atom properties,
spectra, departures, contribution functions, plotting).
"""
from .atom import AtomClass
from .auxfuncs import *
from .decorators import broadcast_arguments
from .plotting_scripts import *
from .fast_atom import FastAtom, PandasAtom
from .flux import integrate_flux, read_angle_quadrature, split_ie
from .io_binary import read_patch_save as _read_patch_save
from .linelist import LineList
from .NameList import NameList
from .spectrum import Spectrum
from .spectrum_grid import SpectrumGrid
from .atmos import AtmosClass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
import pathlib
import functools
from tqdm import tqdm
import multiprocessing as mp

trapz = getattr(np, 'trapz', np.trapezoid)

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
    if isinstance(runs[0], (str, pathlib.Path)):
      runs = read_list(runs)

    return super().__init__(runs, *args, **kwargs)

#=====================================================================
class tumult(object):
    @property
    def atmos(self):
        if not hasattr(self, '_atmos'):
            if hasattr(self, 'runs'):
                self._atmos = AtmosClass.combine([run.atmos for run in self.runs], weights=self.weights)
            else:
                self._atmos = AtmosClass.from_run(self)
        return self._atmos

    @atmos.setter
    def atmos(self, value):
        self._atmos = value

    @property
    def llist(self):
        if not hasattr(self, '_llist'):
            run = self.runs[0] if hasattr(self, 'runs') else self
            self._llist = LineList.from_tumult_format(run.folder / 'linelist_out.txt')
        return self._llist

    @llist.setter
    def llist(self, value):
        self._llist = value

    @property
    def atoms(self):
        if not hasattr(self, '_atoms'):
            try:
              run = self.runs[0] if hasattr(self, 'runs') else self
              self._atoms = run.read_atom(fast_read=True)
            except Exception as e:
              print(f'Error reading atoms: {e}\nTrying again with fast_read=False...')
              self._atoms = self.read_atom(fast_read=False)

        return self._atoms

    @atoms.setter
    def atoms(self, value):
        self._atoms = value

    @classmethod
    #=================================================================
    def read(cls, folder, verbose=0, **kwargs):
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
                return vac2air(np.fromfile(sfolder / file, dtype=np.float64))

        self.lam = read_wavelength_file('out_lam.bin')
        self.bb_lam = read_wavelength_file('bb_lam.bin')
        self.bf_lam = read_wavelength_file('bf_lam.bin')
        self.all_lam = read_wavelength_file('all_lam.bin')
        self.has_spectrum = self.lam is not None and self.lam.size > 0
        self.legacy = not 'read_departures' in self.nml['atmos_params']

        if self.has_spectrum:
            self.has_spectrum = self.nml['rt_output_params']['angle_quad'].strip() != 'none'

        quad = read_angle_quadrature(folder / 'angle_quadrature.txt')
        if quad is not None:
            self.vec = quad["vec"]
            self.wts = quad["wts"]
            self.n_ang = quad["n_ang"]
            self.mus = quad["mus"]
            self.imus = quad["imus"]
            self.nmu = quad["nmu"]

        if self.nml['run_params'].get('make_eos', False):
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

                if self.lam.shape[0] != ie.shape[0]:
                    print(f'Wavelength and Intensity dimensions do not match!\nRun {folder} might have never completed!')

                # explicit layout dispatch (legacy [ie, lte, cnt] vs modern
                # [ie, cnt(, lte)]); fails loudly on an unrecognised component count
                self.ie, self.ie_lte, self.ie_cnt, self.mode = split_ie(ie, self.legacy)

                try:
                    # the Fortran side now integrates the mu-quadrature itself and
                    # saves the result directly (same [ie, lte, cnt] component
                    # layout as 'ie', just with the ray axis already collapsed)
                    flux = self.read_patch_save('flux', lazy=False)[0]
                    flux_i, flux_lte, flux_cnt, _ = split_ie(flux, self.legacy)
                    self.flux     = flux_i[:, 0]
                    self.flux_lte = flux_lte[:, 0]
                    self.flux_cnt = flux_cnt[:, 0]
                except FileNotFoundError:
                    # older runs, from before the Fortran side saved 'flux'
                    # directly: fall back to integrating 'ie' here
                    self.flux     = integrate_flux(self.ie,     self.wts)
                    self.flux_lte = integrate_flux(self.ie_lte, self.wts)
                    self.flux_cnt = integrate_flux(self.ie_cnt, self.wts)

                if self.nml['output_params']['save_resolved']:
                    i3 = self.read_patch_save('i3')[0]
                    self.i3     = i3[...,0]
                    self.c3     = i3[...,1]
                    if self.mode == 'NLTE' and i3.shape[-1] > 2:
                        self.i3_lte = i3[...,2]

        return self

    #=====================================================================
    def read_atom(self, read_atom=True, fast_read=True, verbose=True, **kwargs):
        atoms = []
        atom_files = list(self.folder.glob('atoms/*'))
        atom_files += list(self.folder.glob('atom.txt'))
        moldict = {'CH':'C', 'OH':'O', 'NH':'N', 'CN':'N', 'CO':'C',
                    'NO':'N', 'H2O':'O', 'C2':'C', 'N2':'N', 'O2':'O'}

        if 'atom_params' in self.nml:
            self.qnorm = self.nml['atom_params']['qnorm']
            nq_fac = self.nml['atom_params']['nq_fac']
            if self.nml['atom_params']['legacy_freql']: nq_fac *= 2
            in_abnd = self.nml['atom_params']['abundance']
            in_abnd = [in_abnd] if not hasattr(in_abnd, '__iter__') else in_abnd
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

                    if pandas_format or fast_read:
                        if pandas_format:
                            atom = PandasAtom.read(atom_file)
                        else:
                            atom = FastAtom(atom_file, **kwargs)

                        atom.run = self
                        if not hasattr(atom, 'in_lines'): atom.in_lines = slice(None, None)

                        if self.has_spectrum:
                            atom.line = atom.trans.get_lines(atom.in_lines)

                            if option(self.nml['output_params'], 'overlapping_bb_in_output'):
                                atom.load_sp(self.lam, qnorm=self.qnorm)
                            else:
                                nnu = self.lam.size
                                ib, ir = atom.get_non_overlapping_indices(atom.trans, nnu, nq_fac)
                                atom.bbt['ir'] = 0
                                atom.bbt['ib'] = 0
                                atom.bbt.loc[atom.in_lines, 'ir'] = ir[:atom.nline][atom.in_lines]
                                atom.bbt.loc[atom.in_lines, 'ib'] = ib[:atom.nline][atom.in_lines]

                    else:
                        atom = AtomClass(atom_file, qnorm=self.qnorm, verbose=verbose, **kwargs)
                        atom.run = self

                        if not hasattr(atom, 'in_lines'): atom.in_lines = slice(None, atom.nline)
                        atom.line = np.array(atom.trans)[atom.in_lines]
                        atom.nline = len(atom.line)

                        if self.has_spectrum:
                            if not option(self.nml['output_params'], 'overlapping_bb_in_output'):
                                nnu = self.lam.size
                                ib, ir = atom.get_non_overlapping_indices(atom.trans, nnu, nq_fac)
                                atom.ib = ir.astype(int)
                                atom.ir = ib.astype(int)

                if atom.id in self.abundances:
                    atom.abnd = float(self.abundances[atom.id])
                elif atom.molec:
                    atom.abnd = self.abundances[moldict[atom.id]]
                elif in_abnd and in_abnd[i] > -99:
                    atom.abnd = in_abnd[i]

                atoms.append(atom)

        atoms = np.array(atoms)[np.argsort([atom.id for atom in atoms])]
        reports = self.folder.glob('convergence_report_*.txt')
        self.conv_report = []
        for report in reports:
            df = pd.read_csv(report, sep=r'\s+')
            self.conv_report.append(df)

        return atoms

    #=====================================================================
    @functools.cache
    def get_line(self, lam0, tol=0.05, verbose=True):
        #-----------------------------------------------------------------
        if not hasattr(self.atoms[0], 'line'): self.read_atom()
        all_lines = np.concatenate([[line for line in atom.line] for atom in self.atoms])
        all_lam0 = np.array([line.lam0 for line in all_lines])
        i = np.argmin(abs(all_lam0-lam0))
        line = all_lines[i]
        if verbose and abs(line.lam0 - lam0) > tol:
            print(f'WARNING: Line not found! Closest line to {lam0:.2f}Å found at {line.lam0:.2f}Å')

        return line

    #=====================================================================
    def get_spectrum(self, ang=None, imu=None, LTE=False, cont=False, resolved=False):
        """
        Returns a Spectrum object for the given angle (ang) or mu-index (imu).
        If ang and imu are both None, returns the flux spectrum.

        Parameters:
        -----------
        ang: int or None
            Index of the angle in the angle quadrature. Default is None.
        imu: int or None
            Index of the mu value in the angle quadrature. Default is None.
        LTE: bool
            If True, returns the LTE spectrum. If False, returns the NLTE spectrum. Default is False.
        cont: bool
            If True, returns the continuum spectrum. Default is False.
        resolved: bool
            If True, returns the resolved spectrum. Default is False.
        """
        if not LTE and self.mode=='LTE':
            print('Cannot get NLTE spectrum, because run was in LTE!')
            LTE = True

        if imu is not None:
            if imu >= self.nmu:
                raise ValueError('imu is out of range. Nmu = {:}'.format(self.nmu))
            elif imu < 0:
                imu = self.imus.max() + imu + 1

        if self.n_ang == 1: ang=0

        flux = ang is None and imu is None
        yunit = u.erg/(u.s*u.cm**2*u.angstrom)
        if not flux: yunit = yunit / u.sr

        if cont:
            if flux:
                yy = self.flux_cnt
            elif ang is not None:
                yy = self.ie_cnt[:, ang]
            elif imu is not None:
                imu_mask = self.imus == imu
                yy = np.mean(self.ie_cnt[:, imu_mask], axis=1)
        elif LTE:
            if flux:
                yy = self.flux_lte
            elif ang is not None:
                yy = self.ie_lte[:, ang]
            elif imu is not None:
                imu_mask = self.imus == imu
                yy = np.mean(self.ie_lte[:, imu_mask], axis=1)
        else:
            if flux:
                yy = self.flux
            elif ang is not None:
                yy = self.ie[:, ang]
            elif imu is not None:
                imu_mask = self.imus == imu
                yy = np.mean(self.ie[:, imu_mask], axis=1)

        return Spectrum(self.lam, yy, xunit=u.angstrom, yunit=yunit)

    #=====================================================================
    def crop(self, slc=None, x0=None, ang=None, imu=None, norm=True,
             cont=False, LTE=False, line=None, **kwargs):
        #-----------------------------------------------------------------
        if line is not None:
            line.prepare()
            slc = slice(line.ib, line.ir)
            x0 = line.lam0

        spec = self.get_spectrum(ang=ang, imu=imu, LTE=LTE, cont=cont)
        spec = spec.crop(slc=slc, x0=x0, **kwargs)

        if norm:
          cntm = self.get_spectrum(ang=ang, imu=imu, LTE=LTE, cont=True)
          cntm = cntm.crop(slc=slc, x0=x0, **kwargs)
          spec = spec / cntm

        return spec

    #=================================================================
    def fit(self, spec,  **kwargs):
        return Spectrum.fit(self, spec, **kwargs)

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
        self.atmos.plot_cube_slice(Z, ang=ang, **kwargs)

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
    def get_depart(self, ilvl=None, iatom=0, mean=False):
        #-----------------------------------------------------------------
        atomid = self.atoms[iatom].id if type(iatom) is int else iatom

        if hasattr(self, 'runs'):
            xy = np.array([run.get_depart(ilvl, atomid, mean=mean) for run in self.runs])

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

        if mean and self.atmos.dim == '3D':
            ltau, dep = self.atmos.get_surf_mean('ltau', dep, axis=2)
        else:
            ltau = self.atmos.ltau.squeeze()
            dep = dep.squeeze()

        return ltau, dep

    #=====================================================================
    def plot_depart(self, ilvl=0, iatom=0, mean=True, **kwargs):
        #-----------------------------------------------------------------
        atomid = self.atoms[iatom].id if type(iatom) is int else iatom
        xx, yy = self.get_depart(ilvl, atomid, mean=mean)
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
        """Load ``kind`` from this run's save folder (see :mod:`tumult.io_binary`)."""
        folder = pathlib.Path(folder) if folder else self.sfolder
        return _read_patch_save(folder, kind, n_nu=self.nml['run_params']['n_nu'],
                                concat=concat, fdim=fdim, zdim=zdim, lazy=lazy)

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
                llist = self.llist.reset_index()
                mask = (llist.lam0 > xmin) & (llist.lam0 < xmax)
                llist = llist[mask]
            else:
                df = self.llist.reset_index()[['ion', 'lam0', 'strength', 'loggf', 'ga']]
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

        if show_atom:
            for atom in self.atoms:
                bbt = atom.bbt[['lam0', 'loggf', 'ga']].copy().reset_index(drop=True)
                mask = (bbt.lam0 > xmin) & (bbt.lam0 < xmax)
                bbt = bbt[mask]
                bbt['ion'] = atom.id
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
                            new_row = pd.DataFrame({'ion': [atom.id], 'lam0': [lam0], 'strength': [strength], 'loggf': [loggf]})
                            llist = pd.concat([llist, new_row])

        if len(llist) > 0:
            llist = llist.sort_values(by='lam0').reset_index(drop=True)
            llist['ion'] = llist['ion'].str.strip()
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
        if cfunc.ndim>1 and fmean: cfunc = trapz(cfunc, x=lam, axis=0)
        if norm:
            normalisation = abs(trapz(cfunc, x=ltau))
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
        self.folder = [run.folder for run in runs]

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

        weights = weights or np.ones(len(runs))/len(runs)

        if self.run0.has_spectrum:
            def weight_runs(attr):
                wsum = np.sum([getattr(r, attr) * w for r, w in zip(runs, weights)], axis=0)
                setattr(self, attr, wsum)

            weight_runs('ie')
            weight_runs('ie_lte')
            weight_runs('ie_cnt')
            weight_runs('flux')
            weight_runs('flux_lte')
            weight_runs('flux_cnt')
            self.lam = self.run0.lam

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
