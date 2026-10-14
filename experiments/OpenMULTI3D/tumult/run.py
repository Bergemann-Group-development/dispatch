"""The dispatch run class (`Run`) and its helpers.

Moved from tumult/tumult.py in the Phase 6 decomposition (that module remains as the
package's import shim). The raw patch-save reading lives in :mod:`tumult.io_binary`
and the intensity-layout / flux physics in :mod:`tumult.flux`; this class is the
orchestration and analysis facade over them (lazy atmos/linelist/atom properties,
spectra, departures, contribution functions, plotting).
"""
from .atom import Atom, Level
from .auxfuncs import (assign_group_indices, interp_tau, optical_depth_from_top,
                       rotate_array)
from ._abundance_numerics import cog_invert
from .const import CC
from .physics import vac2air
from .decorators import broadcast_arguments
from .flux import (base_set_indices, integrate_flux, ladder_sets,
                   read_angle_quadrature, read_segments, read_set_manifest,
                   split_ie)
from .io_binary import read_patch_save as _read_patch_save
from .linelist import LineList
from .namelist import NameList
from .spectrum import Spectrum
from .spectrum_grid import SpectrumGrid
from .atmos import AtmosClass
import numpy as np
import pandas as pd
from astropy import units as u
import pathlib
import functools
from tqdm import tqdm
import multiprocessing as mp

trapz = getattr(np, 'trapz', np.trapezoid)

# kinds written by the output RT patches, whose i_nu axis is the output split rather
# than the internal n_windows
OUTPUT_PATCH_KINDS = ('flux', 'ie')

def option(nml, key):
    if key in nml.keys():
        return nml[key]
    else:
        return False

#=====================================================================
def read_run(folder):
  return Run.read(folder, read_atom=False)

#=====================================================================
def read_list(folders):
  with mp.Pool() as pool:
    run_list = list(tqdm(pool.imap(read_run, folders), total=len(folders)))
    return np.array(run_list)

#=====================================================================
class RunGrid(SpectrumGrid):
  def __init__(self, runs, *args, **kwargs):
    if isinstance(runs[0], (str, pathlib.Path)):
      runs = read_list(runs)

    return super().__init__(runs, *args, **kwargs)

#=====================================================================
class Run(object):
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
            run = self.runs[0] if hasattr(self, 'runs') else self
            self._atoms = run.read_atom()

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

        # out_lam.bin holds the vacuum grid the run computed on; keep it beside the
        # air grid so vacuum quantities need no air2vac round trip
        lam_file = sfolder / 'out_lam.bin'
        self.lam_vac = np.fromfile(lam_file, dtype=np.float64) if lam_file.exists() else None
        self.lam = vac2air(self.lam_vac) if self.lam_vac is not None else None
        self.bb_lam = read_wavelength_file('bb_lam.bin')
        self.bf_lam = read_wavelength_file('bf_lam.bin')
        self.all_lam = read_wavelength_file('all_lam.bin')
        self.has_spectrum = self.lam is not None and self.lam.size > 0
        self.legacy = not 'read_departures' in self.nml['atmos_params']
        # the run's own description of its set axis, where it wrote one; None
        # falls split_ie back on the historical count-based layout guess
        self.sets = read_set_manifest(folder)
        # the windows a multi-window run synthesised; None for an ordinary single
        # window, whose grid has no gaps to respect
        self.segments = read_segments(folder)
        self.ladder = ladder_sets(self.sets)

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
                self._ie_sets = ie[:, 0, 0, :, :]

                if self.lam.shape[0] != ie.shape[0]:
                    print(f'Wavelength and Intensity dimensions do not match!\nRun {folder} might have never completed!')

                # explicit layout dispatch (legacy [ie, lte, cnt] vs modern
                # [ie, cnt(, lte)]); fails loudly on an unrecognised component count
                self.ie, self.ie_lte, self.ie_cnt, self.mode = split_ie(
                    ie, self.legacy, self.sets)

                try:
                    # the Fortran side now integrates the mu-quadrature itself and
                    # saves the result directly (same [ie, lte, cnt] component
                    # layout as 'ie', just with the ray axis already collapsed)
                    flux = self.read_patch_save('flux', lazy=False)[0]
                    self._flux_sets = flux[:, 0, 0, 0, :]
                    flux_i, flux_lte, flux_cnt, _ = split_ie(flux, self.legacy, self.sets)
                    self.flux     = flux_i[:, 0]
                    self.flux_lte = flux_lte[:, 0]
                    self.flux_cnt = flux_cnt[:, 0]
                except FileNotFoundError:
                    # older runs, from before the Fortran side saved 'flux'
                    # directly: fall back to integrating 'ie' here
                    self.flux     = integrate_flux(self.ie,     self.wts)
                    self.flux_lte = integrate_flux(self.ie_lte, self.wts)
                    self.flux_cnt = integrate_flux(self.ie_cnt, self.wts)
                    self._flux_sets = np.stack(
                        [integrate_flux(self._ie_sets[..., i], self.wts)
                         for i in range(self._ie_sets.shape[-1])], axis=-1)

                if self.nml['output_params']['save_resolved']:
                    i3 = self.read_patch_save('i3')[0]
                    # 'i3' carries the same set axis as 'ie', so it takes the same
                    # layout: a manifest lookup where the run wrote one, and the
                    # historical [intensity, continuum, lte] order otherwise
                    if self.sets:
                        j, j_lte, j_cnt, _ = base_set_indices(self.sets, i3.shape[-1])
                    else:
                        j, j_cnt, j_lte = 0, 1, 2
                    self.i3     = i3[...,j]
                    self.c3     = i3[...,j_cnt]
                    if self.mode == 'NLTE' and i3.shape[-1] > 2:
                        self.i3_lte = i3[...,j_lte]

        return self

    #=====================================================================
    def read_atom(self):
        """The model atoms in the run's ``atoms/`` folder, sorted by id. A run keeps
        each atom in the table format; a classic MULTI copy (``atoms/atom.*``, from
        runs made before dispatch.x converted its atoms) is converted on reading."""
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
            in_abnd = None

        for i, atom_file in enumerate(atom_files):
            if atom_file.is_file():
                atom = Atom.read(atom_file)
                atom.run = self

                if self.has_spectrum:
                    atom.load_sp(self.lam, qnorm=self.qnorm)
                    if not option(self.nml['output_params'], 'overlapping_bb_in_output'):
                        atom.load_non_overlapping(self.lam.size, nq_fac)

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
        """The bound-bound transition of the run's model atoms nearest to ``lam0`` [Å,
        air], a :class:`tumult.trans.Transition`; ``verbose`` warns when it is further
        than ``tol`` [Å] away."""
        trans = [t for atom in self.atoms for t in atom.bbtrans]
        all_lam0 = np.array([t.lam0 for t in trans])
        line = trans[np.argmin(abs(all_lam0-lam0))]
        if verbose and abs(line.lam0 - lam0) > tol:
            print(f'WARNING: Line not found! Closest line to {lam0:.2f}Å found at {line.lam0:.2f}Å')

        return line

    #=====================================================================
    def get_spectrum(self, ang=None, imu=None, LTE=False, cont=False, resolved=False,
                     index=None):
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
        index: int or None
            Zero-based component of the saved set axis (see :attr:`ladder`) instead of the
            run's own spectrum. A set is one specific component, so ``LTE`` does not apply.
        """
        yy, yunit = self._component(ang, imu, LTE, cont, index)
        return Spectrum(self.lam, yy, xunit=u.angstrom, yunit=yunit,
                        segments=self.segments)

    def _component(self, ang, imu, LTE, cont, index):
        """The raw array and unit :meth:`get_spectrum` wraps."""
        if index is not None:
            if not hasattr(self, '_ie_sets'):
                raise RuntimeError('run has no saved spectrum sets')
            index = int(index)
            if index < 0 or index >= self._ie_sets.shape[-1]:
                raise IndexError(f'set index {index} is out of range')
        elif not LTE and not cont and self.mode == 'LTE':
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

        if index is not None and not cont:
            fluxes, intensities = self._flux_sets[:, index], self._ie_sets[..., index]
        elif cont:
            fluxes, intensities = self.flux_cnt, self.ie_cnt
        elif LTE:
            fluxes, intensities = self.flux_lte, self.ie_lte
        else:
            fluxes, intensities = self.flux, self.ie

        if flux:
            yy = fluxes
        elif ang is not None:
            yy = intensities[:, ang]
        else:
            yy = np.mean(intensities[:, self.imus == imu], axis=1)
        return yy, yunit

    #=====================================================================
    def crop(self, slc=None, x0=None, ang=None, imu=None, norm=True,
             cont=False, LTE=False, line=None, index=None, xmin=None, xmax=None,
             window=None, **kwargs):
        """The spectrum (``index``: that set component, see :meth:`get_spectrum`) cropped
        by :meth:`Spectrum.crop` and, with ``norm``, divided by the continuum.

        A plain ``xmin``/``xmax`` or ``window`` crop slices the stored arrays before
        building the spectrum, so it never materialises the whole run.
        """
        #-----------------------------------------------------------------
        if line is not None:
            slc = slice(line.ib, line.ir)
            x0 = line.lam0

        if window is not None:
            xmin, xmax = window
        if slc is None and x0 is None and not kwargs and xmin is not None \
                and xmax is not None:
            lo, hi = (v.to_value(u.angstrom, equivalencies=u.spectral())
                      if isinstance(v, u.Quantity) else float(v) for v in (xmin, xmax))
            # ascending views, as Spectrum.crop searches and returns them
            order = slice(None) if self.lam[-1] > self.lam[0] else slice(None, None, -1)
            lam = self.lam[order]
            part = slice(np.searchsorted(lam, lo, side='left'),
                         np.searchsorted(lam, hi, side='right'))
            yy, yunit = self._component(ang, imu, LTE, cont, index)
            yy = np.asarray(yy[order][part])
            if norm:
                cntm = self._component(ang, imu, LTE, True, None)[0]
                yy = yy / np.asarray(cntm[order][part])
                yunit = u.dimensionless_unscaled
            return Spectrum(np.asarray(lam[part]), yy, xunit=u.angstrom, yunit=yunit)

        crop = dict(slc=slc, x0=x0, xmin=xmin, xmax=xmax, **kwargs)
        spec = self.get_spectrum(ang=ang, imu=imu, LTE=LTE, cont=cont, index=index)
        spec = spec.crop(**crop)

        if norm:
          cntm = self.get_spectrum(ang=ang, imu=imu, LTE=LTE, cont=True)
          cntm = cntm.crop(**crop)
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
        """The departure coefficients of level ``ilvl`` (all levels when None) against
        log tau_500, as ``(ltau, dep)``. ``ilvl`` is an index from 0 into the levels of
        atom ``iatom`` (its position in :attr:`atoms` or its id), or a
        :class:`tumult.atom.Level`, which names its atom itself (``line.lower``)."""
        #-----------------------------------------------------------------
        if isinstance(ilvl, Level):
            ilvl, iatom = ilvl.i, ilvl.atom.id
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
        """Draw the departure coefficient of one level (see :meth:`get_depart`)."""
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        if isinstance(ilvl, Level):
            ilvl, iatom = ilvl.i, ilvl.atom.id
        atomid = self.atoms[iatom].id if type(iatom) is int else iatom
        xx, yy = self.get_depart(ilvl, atomid, mean=mean)
        plt.title(f'Departure coefficient of {atomid:} level {ilvl:}')
        plt.ylabel(r'$\rm n\,/\,n_{\rm LTE}$')
        plt.xlabel(r'$\log_{10}\tau_{500}$')
        return plt.plot(xx.ravel(), yy.ravel(), **kwargs)


    #=====================================================================
    def get_cog_corr(self, element, LTE_runs, **kwargs):
        """NLTE abundance correction of this run's spectrum on the LTE curve of growth
        of ``LTE_runs``: the LTE abundance that reproduces this run's equivalent width,
        subtracted from this run's abundance of ``element``.

        ``kwargs`` select the line as in :meth:`crop` (``window=``, ``line=``, ...); the
        equivalent widths [mÅ] are trapezoidal integrals of the normalised depression
        over that range. ``LTE_runs`` may be NLTE runs, whose LTE sets are used.
        """
        #-----------------------------------------------------------------
        weq = self.crop(**kwargs).measure(profile="trapz").maa
        abnd = self.abundances[element]
        lte_abnd = np.array([r.abundances[element] for r in LTE_runs])
        lte_weqs = np.array([r.crop(LTE=True, **kwargs).measure(profile="trapz").maa
                             for r in LTE_runs])
        corr = abnd - cog_invert(lte_abnd, lte_weqs, weq)["a"]

        return corr, lte_weqs, lte_abnd, weq, abnd

    #=====================================================================
    def plot_cog_corr(self, element, **kwargs):
        """Draw the NLTE correction of :meth:`get_cog_corr`: the LTE curve of growth,
        this run's NLTE equivalent width at its abundance, and the LTE abundance that
        reproduces that width, NLTE in C0 and LTE in C1."""
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        corr, lte_weqs, lte_abnd, weq, abnd = self.get_cog_corr(element=element, **kwargs)
        order = np.argsort(lte_abnd)
        plt.plot(np.asarray(lte_abnd)[order], np.asarray(lte_weqs)[order], 's--',
                 color='C1', ms=4, label='LTE curve of growth')
        plt.plot(abnd, weq, 'o', color='C0', label='NLTE at A({})'.format(element))
        a0 = abnd - corr
        plt.plot(a0, weq, 'o', mfc='none', color='C1', label='LTE, same W')
        plt.annotate('', xy=(abnd, weq), xytext=(a0, weq),
                     arrowprops=dict(arrowstyle='->', color='k'))
        plt.annotate('$\\Delta A={:.2f}$'.format(corr).replace('-', '\u2212'),
                     xy=(0.5*(a0 + abnd), weq), xytext=(0, -14), textcoords='offset points',
                     ha='center', va='top')
        plt.xlabel('A({})'.format(element))
        plt.ylabel(r'$W$ [mÅ]')
        plt.legend(loc='lower right')

    #=====================================================================
    def read_patch_save(self, kind, folder=None, concat=True, fdim=0, zdim=None, lazy=True):
        """Load ``kind`` from this run's save folder (see :mod:`tumult.io_binary`)."""
        folder = pathlib.Path(folder) if folder else self.sfolder
        n_windows = self.nml['run_params']['n_windows']
        if self.segments and kind in OUTPUT_PATCH_KINDS:
            # A multi-window run splits its output spectrum per window, into as many
            # patches as the Fortran side chose (composition_pre_init), and saves each
            # under abs(i_nu). Every one of those ids belongs to the output spectrum, so
            # the split point must lie above all of them; the grid size is such a bound,
            # since no patch is empty. The internal kinds keep n_windows.
            n_windows = self.segments[-1]['last']
        return _read_patch_save(folder, kind, n_windows=n_windows,
                                concat=concat, fdim=fdim, zdim=zdim, lazy=lazy)

    #=================================================================
    def get_clv(self, center_diff=False, **kwargs):
        """Equivalent width [mÅ] of the spectrum :meth:`crop` selects with ``kwargs``,
        at each μ of the run; ``center_diff`` subtracts the width at the first μ."""
        xx = self.mus
        yy = np.array([self.crop(imu=i, **kwargs).measure(profile="trapz").maa
                       for i in range(len(xx))])
        if center_diff: yy = yy - yy[0]
        return xx, yy

    #=====================================================================
    def plot_line_labels(self, offset=None, tol=0.02, show_atom=True,
                         show_llist=True, dry=False, xmin=None, xmax=None, **kwargs):
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        from .plotting_scripts import plot_line_annotations
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
                df = self.llist.reset_index()[['ion', 'lam0', 'strength', 'loggf', 'grad']]
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
                bbt = atom.bbt[['lam0', 'loggf', 'grad']].copy().reset_index(drop=True)
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
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
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
        """Index range of ``self.lam`` selected by a line, a window or a centre and
        half-width; a centre ``lam0`` alone selects the single nearest point (an int)."""
        slc = slice(None, None)

        if line is not None:
            slc = slice(line.ib, line.ir)
            lam0 = line.lam0
        elif lam0 is not None and dlam is None and xmin is None and xmax is None:
            return int(np.argmin(np.abs(self.lam - lam0)))

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
        """Contribution function per unit log tau_500 of the main spectrum set.

        Returns ``(ltau, cf)`` with ``cf`` shaped (wavelength, depth) in the units of
        ``run.flux`` (no ray selection) or ``run.ie`` (``imu``/``ang``, averaged over the
        selected rays). Integrated over ``ltau`` it gives the line depression
        F_cnt - F (or I_cnt - I) under the default ``depression_cntrbf=T``, and the
        emergent flux (intensity) otherwise.

        dispatch saves, per ray and point, the source term of that quantity per unit
        path length and per Hz together with the total opacity. Here the term is
        attenuated by exp(-tau_nu) along the ray, mapped from path length to
        log tau_500 (ds = ln10 tau_500 / (chi_500 mu) dlog tau_500) and converted to
        per Angstrom. ``fmean`` integrates over the selected wavelengths, ``gmean``
        averages a 3D run horizontally on surfaces of constant tau_500, and ``norm``
        scales the result to unit area in ``ltau``.
        """
        #-----------------------------------------------------------------
        mem  = self.read_patch_save('cntrbf', fdim=0, zdim=3, lazy=True)[-1]
        slc = slc if slc is not None else self.get_slice(**kwargs)

        lam = self.lam[slc]
        # an int selects one wavelength; read it as a length-one range
        cntrbf = mem[slice(slc, slc+1) if isinstance(slc, (int, np.integer)) else slc]
        per_aa = (CC*1e8 / np.atleast_1d(self.lam_vac[slc])**2).reshape(
            (-1,) + (1,)*(cntrbf.ndim - 3))

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

        # the flux quadrature (run.flux = run.ie @ run.wts), or the mean of the rays
        wts = self.wts[angs]
        if ang is not None or imu is not None:
            wts = wts / np.sum(wts)

        conversion = np.log(10) * self.atmos.tau / self.atmos.chi

        first_time = True
        for ang, w in zip(angs, wts):
            cang = np.array(cntrbf[..., ang])
            cf, net_chi = np.moveaxis(cang, -1, 0)
            vec = self.vec[ang]
            path = self.atmos.zz / vec[-1] * 1e8
            tau = optical_depth_from_top(net_chi, path, axis=-1)
            cf = per_aa * conversion / vec[-1] * cf * np.exp(-tau)
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
        ltau = ltau.squeeze()
        cfunc = cfunc.squeeze()
        if cfunc.ndim>1 and fmean: cfunc = trapz(cfunc, x=lam, axis=0)
        if norm:
            normalisation = abs(trapz(cfunc, x=ltau, axis=-1))
            cfunc = cfunc / np.expand_dims(normalisation, -1)
        return ltau, cfunc

    #=====================================================================
    def _cntrbf_label(self):
        """Axis label for the kind of contribution function the run saved."""
        if self.nml['output_params'].get('depression_cntrbf', True):
            return 'Depression contribution function'
        return 'Intensity contribution function'

    #=====================================================================
    def plot_cntrbf(self, line=None, lam0=None, dlam=None, xmin=None, xmax=None,
                    slc=None, imu=None, ang=None, title=None, norm=False, **kwargs):
        """Plot the contribution function against log tau_500 at one wavelength
        (``lam0`` alone), or integrated over a range (``line``, ``lam0`` and ``dlam``,
        ``xmin``/``xmax``); see :meth:`get_cntrbf`."""
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        from .plotting_scripts import set_title
        if slc is None:
            slc = self.get_slice(line=line, dlam=dlam, lam0=lam0, xmin=xmin, xmax=xmax)

        lam = self.lam[slc]
        ltau, cf = self.get_cntrbf(slc=slc, imu=imu, ang=ang, gmean=True, fmean=True, norm=norm)

        if isinstance(slc, (int, np.integer)):
            wave = rf'λ={lam:.3f}[Å]'
        else:
            wave = rf'λ={lam[0]:.2f}-{lam[-1]:.2f}[Å]'

        if ang is None and imu is None:
            angle = 'Flux'
        else:
            if ang is not None:
                imu = self.imus[ang]

            angle = r'μ={:.2f}'.format(self.mus[imu])

        plt.xlabel(r'$\log_{10}\tau_{500}$')
        plt.ylabel(self._cntrbf_label())
        set_title(title, default=wave + ', ' + angle)
        return plt.plot(ltau, cf, **kwargs)


    #=====================================================================
    def plot_cntrbf_map(self, line=None, lam0=None, dlam=None, xmin=None, xmax=None,
                        imu=None, ang=None, title=None, cmap='binary', **kwargs):
        """Map the contribution function over wavelength and log tau_500, divided by
        the continuum at each wavelength, so that a column integrates to the depth
        1 - F/F_cnt of the normalised spectrum (with ``depression_cntrbf=T``)."""
        #-----------------------------------------------------------------
        import matplotlib.pyplot as plt  # deferred: only plotting pays matplotlib's import cost
        from .plotting_scripts import set_title
        slc = self.get_slice(line=line, dlam=dlam, lam0=lam0, xmin=xmin, xmax=xmax)
        lam = self.lam[slc]
        ltau, cf = self.get_cntrbf(slc=slc, imu=imu, ang=ang, gmean=True)
        cnt = np.asarray(self._component(ang, imu, False, True, None)[0])[slc]

        plt.xlabel(r'$\lambda$ [Å]')
        plt.ylabel(r'$\log_{10}\tau_{500}$')
        set_title(title, default=self._cntrbf_label())
        return plt.pcolormesh(lam, ltau, (cf / cnt[:, None]).T, cmap=cmap,
                              shading='nearest', **kwargs)


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
        check_and_combine('segments')

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
            self.lam_vac = self.run0.lam_vac

        return self
