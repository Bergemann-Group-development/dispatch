"""In-process interface to the HOP EOS/opacity engine (libhop.so).

Wraps the standalone HOP library (built with ``make -C HOP/standalone``) via
ctypes, exposing the same composition_type solves that dispatch.x uses:

- :class:`HopEOS` -- allocate/initialize a composition instance, solve the EOS
  on arbitrary (T, rho) or (T, ne) points, compute continuum absorption and
  scattering opacities on a wavelength grid.
- :func:`make_eos_table` -- replicate the Fortran ``make_eos`` table
  (``init_eos_comp``): the 6 planes (temp, rho, ne, pg, E, chi500) on a given
  set of (T, rho) points, so EOS tables can be produced without a dispatch.x
  run.  :func:`logstep` builds the ln(T) x ln(rho) node grid dispatch uses.
- :func:`make_opac_table` -- the other half of ``make_eos``: the 'opa' and
  'scat' tables of ``save_opac_table``, on the same points, optionally with LTE
  line-list opacity, split over the wavelength axis and computed in a pool of
  worker processes.

Limitations: one instance per process (HOP keeps configuration in module-level
globals); hard failures inside the Fortran (missing data files, solver
non-convergence) terminate the process, so validate inputs beforehand.  Thread
count is controlled with OMP_NUM_THREADS (latched at the first parallel region).
If numpy links a different OpenMP runtime (MKL/iomp), import this module before
heavy numpy use or set MKL_THREADING_LAYER=GNU in case of runtime clashes.
"""
import ctypes
import multiprocessing
import os
from pathlib import Path

import numpy as np

from .io_binary import write_patch_bin, write_patch_meta, write_patch_save

__all__ = ['HopEOS', 'HopError', 'make_eos_table', 'make_opac_table',
           'logstep', 'wavelength_grid', 'freq_split', 'find_library']

_EXP_DIR = Path(__file__).resolve().parents[1]
_DEFAULT_LIB = _EXP_DIR / 'HOP' / 'standalone' / 'build' / 'libhop.so'

_F32_HUGE = np.float32(np.finfo(np.float32).max)


class HopError(RuntimeError):
    """Nonzero status returned by the HOP C API."""


def find_library():
    """Locate libhop.so: $HOP_LIB if set, else HOP/standalone/build/libhop.so."""
    lib = os.environ.get('HOP_LIB', '')
    path = Path(lib) if lib else _DEFAULT_LIB
    if not path.is_file():
        raise HopError(f"libhop.so not found at {path}; build it with "
                       f"'make -C {_EXP_DIR / 'HOP' / 'standalone'}' "
                       f"or point $HOP_LIB at it")
    return path


def _c_str(s):
    return str(s).encode() + b'\0'


def _verbosity(verbose):
    """HOP's integer verbosity from a bool: False = -1 (nothing at all), True = 1.

    Any other value passes through, so the dispatch levels (0 = the linelist
    line count only, 1 = progress messages, 2 = per-solve diagnostics) remain
    available as ints.
    """
    return (1 if verbose else -1) if isinstance(verbose, bool) else int(verbose)


_f64 = np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS')
_f32 = np.ctypeslib.ndpointer(np.float32, flags='C_CONTIGUOUS')

_c = ctypes.c_char_p
_i = ctypes.c_int
_r4 = ctypes.c_float
_r8 = ctypes.c_double

#: The libhop entry points, mirroring the BIND(c, name=...) declarations in
#: HOP/standalone/hop_c_api.f90 -- name: (restype, argtypes).  tests/test_hop.py
#: checks this against the Fortran source, so keep it in sync when the API grows.
_SIGNATURES = {
    'hop_init':             (_i, [_c, _c, _c, _c, _r8, _r4, _r4, _r4]),
    'hop_free':             (_i, []),
    'hop_set_option_str':   (_i, [_c, _c]),
    'hop_set_option_real':  (_i, [_c, _r8]),
    'hop_set_option_int':   (_i, [_c, _i]),
    'hop_load_state':       (_i, [_i, _f64, _f64, _f64, _r4, _i, _i]),
    'hop_get_n':            (_i, []),
    'hop_get_field':        (_i, [_c, _i, _f64]),
    'hop_get_grph':         (_r8, []),
    'hop_keep_species':     (_i, [_c]),
    'hop_get_species':      (_i, [_c, _i, _f64]),
    'hop_absorption':       (_i, [_i, _f64, _f32, _i]),
    'hop_line_absorption':  (_i, [_i, _f64, _f32, _r8, _i]),
    'hop_init_linelist':    (_i, [_c, _r8, _i, _f64]),
    'hop_scattering':       (_i, [_i, _f64, _f32]),
}


class HopEOS:
    """The HOP composition_type, driven from Python.

    Parameters mirror the &composition_params namelist; data-file paths default
    to the standard files under <experiment>/input/.  ``options`` is a dict of
    HOP module-level settings applied *before* the Fortran init (required for
    init-time options such as ``abundances`` or ``exclude_elements``); the
    recognized names are the ``case`` labels of ``hop_set_option_str/real/int``
    in ``HOP/standalone/hop_c_api.f90``.

    The library is silent by default (``verbose=False``, i.e. HOP verbosity -1);
    ``verbose=0`` restores what a dispatch run prints at its own default -- the
    'N of M lines read from linelist' line -- ``verbose=True``/1 adds HOP's
    progress messages (which abundances were overridden, ...) and ``verbose=2``
    the per-solve diagnostics.

    One live instance per process -- constructing a second without closing the
    first raises.  Use as a context manager or call close().
    """

    _open = False  # class-level: tracks the single Fortran-side instance

    def __init__(self, absdat_file=None, abund_file=None, feh=0.0,
                 isotope_file=None, molec_file=None, tmolim=-99.0,
                 alpha_enhancement=0.0, r_process_enhancement=0.0,
                 verbose=False, **options):
        self._lib = None
        self._owns = False
        if HopEOS._open:
            raise HopError('another HopEOS instance is already open in this '
                           'process (HOP holds global state); close() it first')

        absdat_file = Path(absdat_file or _EXP_DIR / 'input' / 'TS_absdat.dat')
        abund_file = Path(abund_file or _EXP_DIR / 'input' / 'abund' / 'abund_magg22')
        isotope_file = Path(isotope_file) if isotope_file else _EXP_DIR / 'input' / 'isotopes_blp'
        molec_file = Path(molec_file) if molec_file else _EXP_DIR / 'input' / 'molec_data_BC.txt'
        for p in (absdat_file, abund_file, isotope_file, molec_file):
            if not p.is_file():
                raise HopError(f'data file not found: {p}')

        self._lib = self._load_library()

        options = dict(options or {})
        options['verbose'] = _verbosity(verbose)
        # the VCS hydrogen-profile tables sit next to the other input data, but
        # HOP's default path is cwd-relative -- make it absolute unless overridden
        options.setdefault('vcs_data_dir', str(_EXP_DIR / 'input' / 'vcs'))
        for name, value in options.items():
            self.set_option(name, value)

        status = self._lib.hop_init(
            _c_str(absdat_file), _c_str(abund_file),
            _c_str(isotope_file), _c_str(molec_file),
            float(feh), float(tmolim),
            float(alpha_enhancement), float(r_process_enhancement))
        if status != 0:
            raise HopError(f'hop_init failed with status {status}')
        HopEOS._open = True
        self._owns = True

    @staticmethod
    def _load_library():
        lib = ctypes.CDLL(str(find_library()), mode=ctypes.RTLD_GLOBAL)
        for name, (restype, argtypes) in _SIGNATURES.items():
            func = getattr(lib, name)
            func.restype = restype
            func.argtypes = argtypes
        return lib

    # -- configuration ------------------------------------------------------

    def set_option(self, name, value):
        """Set a HOP module-level option (namelist-style global).

        Which of the three libhop setters is used follows the Python type of
        `value`: str -> hop_set_option_str, bool/int -> ..._int (retried as
        ..._real, so an integer-valued `ref_lam=5000` works), everything else
        -> ..._real.  libhop itself is the authority on option names; an
        unrecognized one comes back as status 1 and is raised here.
        """
        key = _c_str(name.lower())
        if isinstance(value, str):
            attempts = [('hop_set_option_str', _c_str)]
        elif isinstance(value, (bool, int, np.integer)):
            # ints are valid input for a real-valued option too (ref_lam=5000)
            attempts = [('hop_set_option_int', int), ('hop_set_option_real', float)]
        else:
            attempts = [('hop_set_option_real', float)]

        for setter, cast in attempts:
            status = getattr(self._lib, setter)(key, cast(value))
            if status == 0:
                return
            if status != 1:
                break  # name was recognized; retrying the other bucket would mask why
        if status == 2:
            raise HopError(f'option {name} needs an initialized HOP instance '
                           f'(set it after hop_init, not before)')
        raise HopError(f'option {name} rejected by libhop '
                       f'(unknown name or wrong type bucket, status {status})')

    def init_linelist(self, files, dlam=100.0, lam=None):
        """Read line-list file(s) into this instance (the &linelist_params path).

        files: a file, a list of files, or a directory (read in sorted order --
        dispatch reads linelist_folder in the order os%get_filenames returns, so
        ties between identical lines may end up ordered differently).  lam, when
        given, restricts the lines kept to [lam[0]-dlam, lam[-1]+dlam], exactly as
        composition_pre_init does with the output spectrum's range -- pass the
        wavelengths you are going to compute opacities on.

        Must be called before load_state: each solve works on a copy of the
        configured instance, so lines read afterwards would not be seen.
        """
        if isinstance(files, (str, Path)):
            p = Path(files)
            files = sorted(f for f in p.iterdir() if f.is_file()) if p.is_dir() else [p]
        if lam is None:
            lam_arr, nlam = np.zeros(1), 0
        else:
            lam_arr = np.ascontiguousarray(np.atleast_1d(lam), np.float64)
            nlam = lam_arr.size
        for f in files:
            f = Path(f)
            if not f.is_file():
                raise HopError(f'line-list file not found: {f}')
            status = self._lib.hop_init_linelist(_c_str(f), float(dlam), nlam, lam_arr)
            if status == 2:
                raise HopError('hop_init_linelist needs an initialized instance')
            if status != 0:
                raise HopError(f'hop_init_linelist({f}) failed with status {status}')
        return self

    # -- EOS ----------------------------------------------------------------

    def keep_species(self, ids):
        """Mark species so their number densities survive subsequent load_state
        calls (retrievable with species_density).  ids: a species id or list of
        ids -- 'Fe I', 'H -', 'CH', ... as named in the atomic/molecular data.
        """
        if isinstance(ids, str):
            ids = [ids]
        for sid in ids:
            status = self._lib.hop_keep_species(_c_str(sid))
            if status != 0:
                raise HopError(f'unknown species id {sid!r} (use e.g. "Fe I" '
                               f'for ions -- a bare element id is not a species '
                               f'-- or a molecule formula like "CH")')

    def species_density(self, species_id):
        """Number density [1/cm3] of one species ('Fe I', 'CH', ...) in the
        loaded state.  The species must have been kept before the load_state
        call (species= argument or keep_species()), or be one HOP retains
        anyway (continuum-opacity / line-list species)."""
        n = self.n
        if n == 0:
            raise HopError('no state loaded; call load_state first')
        out = np.empty(n, np.float64)
        status = self._lib.hop_get_species(_c_str(species_id), n, out)
        if status == 4:
            raise HopError(f'density of {species_id!r} was not retained; pass '
                           f'species=[...] to load_state (or call keep_species) '
                           f'before solving')
        if status != 0:
            raise HopError(f'unknown species id {species_id!r} (use e.g. "Fe I" '
                           f'or a molecule formula like "CH")')
        return out

    def load_state(self, temp, rho=None, ne=None, vmic=-1.0, species=None):
        """Solve the EOS on a set of cells.

        temp [K] plus either rho [g/cm3] (electron density solved by Newton
        iteration; pass ne as well to provide a starting guess) or ne [1/cm3]
        (rho follows from particle conservation).  vmic [km/s] only matters for
        later line-opacity calls; < 0 means unset (the make_eos default).
        species: id(s) to keep for species_density() ('Fe I', 'CH', ...).
        Returns self; results are exposed as the array properties temp, ne,
        rho, pe, pg, E, eh, chi500.
        """
        if species is not None:
            self.keep_species(species)
        temp = np.ascontiguousarray(np.atleast_1d(temp), np.float64)
        n = temp.size
        use_rho, use_ne = rho is not None, ne is not None
        if not (use_rho or use_ne):
            raise HopError('load_state needs rho and/or ne')
        rho_in = np.ascontiguousarray(np.atleast_1d(rho if use_rho else np.zeros(n)), np.float64)
        ne_io = np.ascontiguousarray(np.atleast_1d(ne if use_ne else np.zeros(n)), np.float64)
        if rho_in.size != n or ne_io.size != n:
            raise HopError('temp, rho, ne must have equal size')
        if not np.all(np.isfinite(temp)) or np.any(temp <= 0):
            raise HopError('temp must be finite and positive')

        status = self._lib.hop_load_state(n, temp, ne_io, rho_in,
                                          float(vmic), int(use_rho), int(use_ne))
        if status != 0:
            raise HopError(f'hop_load_state failed with status {status}')
        return self

    @property
    def n(self):
        """Number of cells in the currently loaded state (0 if none)."""
        return self._lib.hop_get_n()

    def _get_field(self, name):
        n = self.n
        if n == 0:
            raise HopError('no state loaded; call load_state first')
        out = np.empty(n, np.float64)
        status = self._lib.hop_get_field(_c_str(name), n, out)
        if status != 0:
            raise HopError(f'hop_get_field({name}) failed with status {status}')
        return out

    temp = property(lambda self: self._get_field('temp'), doc='temperature [K]')
    ne = property(lambda self: self._get_field('ne'), doc='electron density [1/cm3]')
    rho = property(lambda self: self._get_field('rho'), doc='mass density [g/cm3]')
    pe = property(lambda self: self._get_field('pe'), doc='electron pressure [dyn/cm2]')
    pg = property(lambda self: self._get_field('pg'), doc='gas pressure [dyn/cm2]')
    E = property(lambda self: self._get_field('e'), doc='internal energy [erg/g]')
    eh = property(lambda self: self._get_field('eh'), doc='hydrogen ionization/dissociation energy term')
    chi500 = property(lambda self: self._get_field('chi500'),
                      doc='extinction at the reference wavelength (air 5000 A) [1/cm]')

    @property
    def grph(self):
        """Grams of matter per hydrogen nucleus."""
        return self._lib.hop_get_grph()

    # -- opacities ----------------------------------------------------------

    def _check_lam(self, lam_vac):
        lam = np.ascontiguousarray(np.atleast_1d(lam_vac), np.float64)
        if lam.ndim != 1 or lam.size < 1:
            raise HopError('lam_vac must be a 1D array')
        if lam.size > 1 and not np.all(np.diff(lam) > 0):
            raise HopError('lam_vac must be strictly ascending (vacuum Angstrom)')
        n = self.n
        if n == 0:
            raise HopError('no state loaded; call load_state first')
        return lam, n

    def _opacity(self, cfunc, lam_vac, *extra):
        lam, n = self._check_lam(lam_vac)
        out = np.zeros((n, lam.size), np.float32)  # C order == Fortran (nlam, n)
        status = cfunc(lam.size, lam, out, *extra)
        if status != 0:
            raise HopError(f'{cfunc.__name__} failed with status {status}')
        return out

    def absorption(self, lam_vac, include_absmet=False, include_linelist=False,
                   dlam=100.0):
        """Absorption extinction chi [1/cm] on vacuum wavelengths.

        Returns float32 (n_cells, n_lam).  include_absmet adds the sampled
        metal-line background (osmet) contribution.  include_linelist adds the
        LTE line-list opacity (init_linelist must have been called), reproducing
        what rt_patch_mod%generate_opacity_table bakes into the make_eos table;
        dlam [AA] is the &linelist_params line-wing search radius.
        """
        if include_linelist:
            return self._opacity(self._lib.hop_line_absorption, lam_vac,
                                 float(dlam), int(bool(include_absmet)))
        return self._opacity(self._lib.hop_absorption, lam_vac,
                             int(bool(include_absmet)))

    def scattering(self, lam_vac):
        """Electron (Thomson) + Rayleigh scattering extinction [1/cm], float32 (n_cells, n_lam)."""
        return self._opacity(self._lib.hop_scattering, lam_vac)

    # -- lifetime -----------------------------------------------------------

    def close(self):
        # only the instance that actually owns the Fortran state may free it:
        # a closed instance being garbage-collected must not tear down the
        # successor that has meanwhile taken over the process-wide globals
        if self._owns:
            self._owns = False
            HopEOS._open = False
            self._lib.hop_free()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

def logstep(start, stop, step):
    """Like np.logspace but with a ln-space step size rather than a count."""
    lmin = np.log(start)
    lmax = np.log(stop)
    n = int(np.ceil((lmax - lmin) / step))
    return np.exp(lmin + step * np.arange(n))


def _points(temp, rho):
    """The (T, rho) table points as float64 arrays, which must share their shape."""
    temp, rho = np.asarray(temp, np.float64), np.asarray(rho, np.float64)
    if temp.shape != rho.shape:
        raise HopError(f'temp and rho must have the same shape, got '
                       f'{temp.shape} and {rho.shape}')
    return temp, rho


def _flat_cells(temp, rho):
    """The (temp, rho) points as flat cell lists, in Fortran order -- so on the
    usual meshgrid(tab_temp, tab_rho, indexing='ij') grid the temperature index
    runs fastest, as it does in the cells of a dispatch make_eos run."""
    temp, rho = _points(temp, rho)
    return temp.ravel(order='F'), rho.ravel(order='F')


def _unflatten(arr, shape):
    """Fold a (nlam, ncell) array back onto (nlam,) + shape, cells in Fortran order."""
    nd = len(shape)
    arr = arr.reshape((arr.shape[0],) + tuple(shape[::-1]))
    return np.ascontiguousarray(arr.transpose((0,) + tuple(range(nd, 0, -1))))


def freq_split(nnu, i, n_nu):
    """Wavelength range [start, stop) of split i (1-based) out of n_nu.

    Python (0-based, half-open) port of HOP's unequal_intervals, which is what
    ray_patch_mod%set_freq_indices uses to hand each rt_patch its slice -- so
    split i here holds exactly the wavelengths of the patch saved with |i_nu| = i.
    """
    nv, md = divmod(nnu, n_nu)
    j = (abs(i) - 1) % n_nu + 1
    return (j - 1) * nv + min(md, j - 1), j * nv + min(md, j)


def wavelength_grid(aa_blue, aa_red, R=0, daa=0, nnu=0, in_air=True):
    """The &spectrum_params output wavelength grid, in vacuum Angstrom.

    Mirrors spectrum_type%init: with R, a geometric (constant-R) grid built
    between the air->vacuum converted end points; otherwise an arithmetic grid
    from daa or nnu, converted to vacuum afterwards.
    """
    from .physics import air2vac
    if R:
        blue, red = (air2vac(aa_blue), air2vac(aa_red)) if in_air else (aa_blue, aa_red)
        fac = 1.0 + 1.0 / R
        n = int(np.ceil(np.log(red / blue) / np.log(fac))) + 1
        # accumulate rather than blue*fac**k: the Fortran recurrence is
        # lam(i) = lam(i-1)*fac, and the two differ in the last bits
        lam = np.full(n, fac)
        lam[0] = blue
        return np.multiply.accumulate(lam)
    if daa:
        nnu = int(round((aa_red - aa_blue) / daa)) + 1
    elif nnu:
        daa = (aa_red - aa_blue) / (nnu - 1)
    else:
        raise HopError('give R, daa or nnu together with aa_blue/aa_red')
    lam = aa_blue + daa * np.arange(nnu)
    return air2vac(lam) if in_air else lam


def make_eos_table(hop, temp, rho, vmic=None):
    """Replicate the Fortran make_eos table (init_eos_comp) in Python.

    temp/rho: the temperature [K] and density [g/cm^3] of the table points, two
    arrays of the same shape -- normally the 2D (ntemp, nrho) node grid of
    ``np.meshgrid(logstep(...), logstep(...), indexing='ij')``, which is the grid dispatch
    itself tabulates on.  vmic defaults to the Fortran opac_table_vmic default
    (-1 = unset); note that a dispatch run does not necessarily use that value
    -- read_external overwrites opac_table_vmic with the mean of the
    atmosphere's own vmic whenever any of it is positive.  It does not enter the
    EOS, but it does broaden the lines of an opacity table, so see
    make_opac_table.

    Returns the float32 (6,) + temp.shape table with the planes
    (temp, rho, ne, pg, E, chi500) exactly as dispatch.x writes them.
    """
    tt, rr = _flat_cells(temp, rho)
    hop.load_state(tt, rho=rr, vmic=-1.0 if vmic is None else vmic)

    shape = np.shape(temp)
    eos = np.empty((6,) + shape, np.float32)
    for k, field in enumerate(('temp', 'rho', 'ne', 'pg', 'e', 'chi500')):
        vals = hop._get_field(field)
        if field != 'chi500':  # chi500 is float32 on the Fortran side already
            vals = np.clip(vals, -_F32_HUGE, _F32_HUGE)
        eos[k] = vals.reshape(shape, order='F').astype(np.float32)
    return eos


def _opac_chunk(job, hop=None):
    """One wavelength chunk of make_opac_table.

    Module level and argument-only so it survives pickling to a spawned worker,
    which calls it with hop=None and so opens (and closes) its own instance.  In
    the sequential path the caller passes an already configured instance --
    including its line list, read once over the full range rather than per chunk.
    """
    os.environ.setdefault('OMP_NUM_THREADS', '1')  # before the library is loaded
    lam = job['lam']
    shape = job['temp'].shape
    tt, rr = _flat_cells(job['temp'], job['rho'])

    own = hop is None
    hop = HopEOS(**job['hop_kwargs']) if own else hop
    try:
        if own and job['linelist'] is not None:
            hop.init_linelist(job['linelist'], dlam=job['dlam'], lam=lam)
        hop.load_state(tt, rho=rr, vmic=job['vmic'])
        chi = hop.absorption(lam, include_absmet=job['include_absmet'],
                             include_linelist=job['linelist'] is not None,
                             dlam=job['dlam'])
        sig = hop.scattering(lam)
    finally:
        if own:
            hop.close()

    # (ncell, nlam) [1/cm] -> (nlam, ntemp, nrho) [cm2/g], as save_opac_table
    # writes it: cells run temperature-fastest, hence the Fortran-order unfold
    over_rho = (1.0 / rr).astype(np.float32)
    out = [_unflatten((arr * over_rho[:, None]).T, shape) for arr in (chi, sig)]
    if job['folder'] is None:
        return tuple(out)
    return tuple(write_patch_bin(job['folder'], kind, job['key'], arr)
                 for kind, arr in zip(('opa', 'scat'), out))


def make_opac_table(temp, rho, lam, linelist=None, dlam=100.0,
                    include_absmet=False, vmic=None, nchunks=None, nprocs=1,
                    folder=None, iz=1, hop=None, write_eos=True, verbose=False,
                    **hop_kwargs):
    """Build the make_eos opacity tables ('opa', 'scat') without running dispatch.x.

    Reproduces rt_patch_mod%generate_opacity_table + save_opac_table: 'opa' is
    the total absorption chi/rho [cm2/g] (continuum, plus the LTE line list when
    one is given), 'scat' is (electron + Rayleigh scattering)/rho, both float32
    (nlam,) + temp.shape.

    temp/rho: the table points, two arrays of the same shape -- normally the 2D
    (ntemp, nrho) node grid of ``np.meshgrid(logstep(...), logstep(...), indexing='ij')``,
    which is what dispatch tabulates on and what the save layout below expects.
    lam: vacuum wavelengths, e.g. from :func:`wavelength_grid`.  linelist: file,
    list of files or folder for the line-list opacity (None = continuum only).
    hop_kwargs go to HopEOS (feh, abund_file, module options, ...) and are used
    for the workers as well.

    verbose: the HOP verbosity of the instances made here, so False (the
    default) keeps the Fortran side quiet -- with a line list that is one 'N of
    M lines read from linelist' line per chunk -- and True (or an int, see
    :class:`HopEOS`) lets it print.  An instance passed as ``hop`` keeps the
    verbosity it was constructed with.

    vmic [km/s] broadens the lines and must match what dispatch would use, or
    the line opacity comes out wrong (an order of magnitude in the cool, dense
    corner): the &opac_table_params default is -1 (thermal only), but
    read_external_mod overrides it with the mean of the atmosphere's own vmic
    whenever any of that is positive -- so pass that mean here for an
    atmosphere that carries microturbulence.

    The wavelength axis is split into ``nchunks`` pieces (default: ``nprocs``)
    the same way dispatch splits it over its n_nu rt_patches, and the pieces run
    in ``nprocs`` processes -- HOP keeps process-wide global state, so each
    worker holds its own instance and reads only the lines its chunk needs.
    With ``nprocs=1`` everything runs in this process, on ``hop`` if you pass an
    open instance (which is then the one the line list is read into) and on a
    temporary one otherwise.  Because the workers are spawned, a script calling
    this with nprocs > 1 must guard its top level with ``if __name__ ==
    '__main__':`` -- without it each worker re-runs the script on import and the
    pool spawns forever (notebooks and pytest are unaffected).

    With ``folder``, each chunk is written straight to disk in the dispatch save
    layout (``opa_<iz><chunk>.bin`` + ``opa_meta.txt``, likewise 'scat', plus the
    'eos' table when ``write_eos``), so ``tumult.run``'s read_patch_save reads
    the result unchanged and the parent never holds the whole table.

    Returns a dict with lam/temp/rho, eos (unless write_eos is False), folder,
    and opa/scat -- the stitched arrays when writing to memory, None when they
    went to ``folder``.
    """
    temp, rho = _points(temp, rho)  # checked here, before anything is computed
    hop_kwargs = dict(hop_kwargs, verbose=verbose)
    lam = np.ascontiguousarray(np.atleast_1d(lam), np.float64)
    nchunks = min(nchunks or nprocs, lam.size)
    if linelist is not None and not isinstance(linelist, (str, Path)):
        linelist = [str(f) for f in linelist]
    elif linelist is not None:
        linelist = str(linelist)

    jobs = []
    for i in range(1, nchunks + 1):
        fs, fe = freq_split(lam.size, i, nchunks)
        jobs.append(dict(lam=lam[fs:fe], temp=temp, rho=rho,
                         linelist=linelist, dlam=dlam,
                         vmic=-1.0 if vmic is None else vmic,
                         include_absmet=include_absmet, folder=folder,
                         key=(iz, i), hop_kwargs=hop_kwargs))

    own_hop = False
    try:
        if nprocs > 1:
            # spawn, not fork: a forked child would inherit the parent's Fortran
            # globals (and its live instance, which HopEOS allows only one of)
            ctx = multiprocessing.get_context('spawn')
            pool = ctx.Pool(min(nprocs, nchunks))
            try:
                results = pool.map(_opac_chunk, jobs, chunksize=1)
                # close/join, not the with-block's terminate(): the workers have
                # to exit normally for the Fortran side to flush what it wrote
                pool.close()
                pool.join()
            finally:
                pool.terminate()
        else:
            own_hop = hop is None
            if own_hop:
                hop = HopEOS(**hop_kwargs)
            if linelist is not None:
                hop.init_linelist(linelist, dlam=dlam, lam=lam)
            results = [_opac_chunk(job, hop) for job in jobs]

        out = dict(lam=lam, temp=temp, rho=rho, folder=folder,
                   opa=None, scat=None)
        if folder is None:
            out['opa'] = np.concatenate([r[0] for r in results], axis=0)
            out['scat'] = np.concatenate([r[1] for r in results], axis=0)
        else:
            for k, kind in enumerate(('opa', 'scat')):
                write_patch_meta(folder, kind, [r[k] for r in results])

        if write_eos:
            # after the pool, so no worker can inherit the parent's own instance
            own_hop = own_hop or hop is None
            if hop is None:
                hop = HopEOS(**hop_kwargs)
            out['eos'] = make_eos_table(hop, temp, rho, vmic=vmic)
            if folder is not None:
                write_patch_save(folder, 'eos', {(iz, 0): out['eos']})
    finally:
        if own_hop and hop is not None:
            hop.close()
    return out
