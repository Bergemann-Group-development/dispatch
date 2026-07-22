"""In-process interface to the HOP EOS/opacity engine (libhop.so).

Wraps the standalone HOP library (built with ``make -C HOP/standalone``) via
ctypes, exposing the same composition_type solves that dispatch.x uses:

- :class:`HopEOS` -- allocate/initialize a composition instance, solve the EOS
  on arbitrary (T, rho) or (T, ne) points, compute continuum absorption and
  scattering opacities on a wavelength grid.
- :func:`make_eos_table` -- replicate the Fortran ``make_eos`` table
  (``init_eos_comp``): a regular ln(T) x ln(rho) grid with the 6 planes
  (temp, rho, ne, pg, E, chi500), so EOS tables can be produced without a
  dispatch.x run.

Limitations: one instance per process (HOP keeps configuration in module-level
globals); hard failures inside the Fortran (missing data files, solver
non-convergence) terminate the process, so validate inputs beforehand.  Thread
count is controlled with OMP_NUM_THREADS (latched at the first parallel region).
If numpy links a different OpenMP runtime (MKL/iomp), import this module before
heavy numpy use or set MKL_THREADING_LAYER=GNU in case of runtime clashes.
"""
import ctypes
import os
from pathlib import Path

import numpy as np

__all__ = ['HopEOS', 'HopError', 'make_eos_table', 'find_library']

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

    Initialization is silent by default; ``verbose=1`` turns on HOP's own
    progress messages (which abundances were overridden, ...), ``verbose=2``
    adds the per-solve diagnostics.

    One live instance per process -- constructing a second without closing the
    first raises.  Use as a context manager or call close().
    """

    _open = False  # class-level: tracks the single Fortran-side instance

    def __init__(self, absdat_file=None, abund_file=None, feh=0.0,
                 isotope_file=None, molec_file=None, tmolim=-99.0,
                 alpha_enhancement=0.0, r_process_enhancement=0.0,
                 **options):
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

    def absorption(self, lam_vac, include_absmet=False):
        """Continuum absorption extinction chi [1/cm] on vacuum wavelengths.

        Returns float32 (n_cells, n_lam).  include_absmet adds the sampled
        metal-line background (osmet) contribution.
        """
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


def make_eos_table(hop, temp, rho, ldtemp=0.05, ldrho=0.25, vmic=None):
    """Replicate the Fortran make_eos table (init_eos_comp) in Python.

    temp/rho: atmosphere temperature/density values (any shape) or (min, max)
    pairs -- only their extrema matter.  ldtemp/ldrho are the ln-space grid
    steps from &composition_params/&opac_table_params.  vmic defaults to the
    Fortran opac_table_vmic default (-1 = unset), which is what make_eos uses.

    Returns (eos, tab_temp, tab_rho): eos is float32 (6, ntemp, nrho) with the
    planes (temp, rho, ne, pg, E, chi500) exactly as dispatch.x writes them;
    tab_temp/tab_rho are the 1D grid node values.
    """
    lnt0, lnt1 = np.log(np.min(temp)), np.log(np.max(temp))
    lnr0, lnr1 = np.log(np.min(rho)), np.log(np.max(rho))
    ntemp = int(np.ceil((lnt1 - lnt0) / ldtemp))
    nrho = int(np.ceil((lnr1 - lnr0) / ldrho))
    tab_temp = np.exp(lnt0 + ldtemp * np.arange(ntemp))
    tab_rho = np.exp(lnr0 + ldrho * np.arange(nrho))

    # flatten Fortran-style: temperature index fastest, matching (ntemp, nrho, 1)
    tt = np.tile(tab_temp, nrho)
    rr = np.repeat(tab_rho, ntemp)
    hop.load_state(tt, rho=rr, vmic=-1.0 if vmic is None else vmic)

    eos = np.empty((6, ntemp, nrho), np.float32)
    for k, field in enumerate(('temp', 'rho', 'ne', 'pg', 'e', 'chi500')):
        vals = hop._get_field(field)
        if field != 'chi500':  # chi500 is float32 on the Fortran side already
            vals = np.clip(vals, -_F32_HUGE, _F32_HUGE)
        eos[k] = vals.reshape(ntemp, nrho, order='F').astype(np.float32)
    return eos, tab_temp, tab_rho
