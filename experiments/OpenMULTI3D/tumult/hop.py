"""In-process interface to the HOP EOS/opacity engine (libhop.so).

Wraps the standalone HOP library (built with ``make -C HOP/standalone``) via
cffi in ABI mode, exposing the same composition_type solves that dispatch.x uses:

- :class:`HopEOS` -- allocate/initialize a composition instance, solve the EOS on
  arbitrary (T, rho) or (T, ne) points (``load_state``, after which the state is
  readable as the array properties temp/rho/ne/pe/pg/energy/chi500), and compute
  opacities on a wavelength grid (``opacity``, the one entry point for continuum,
  sampled metal-line and line-list opacity).
- :func:`logstep` / :func:`wavelength_grid` build the ln(T) x ln(rho) node grid
  and the &spectrum_params wavelength grid a dispatch run would tabulate on, so a
  dispatch EOS/opacity table can be reproduced without running dispatch.x.

Limitations: one instance per process (HOP keeps configuration in module-level
globals); hard failures inside the Fortran (missing data files, solver
non-convergence) terminate the process, so validate inputs beforehand.  Thread
count is controlled with OMP_NUM_THREADS (latched at the first parallel region).
If numpy links a different OpenMP runtime (MKL/iomp), import this module before
heavy numpy use or set MKL_THREADING_LAYER=GNU in case of runtime clashes.
"""
import multiprocessing
import os
from pathlib import Path

import cffi
import numpy as np

__all__ = ['HopEOS', 'HopError', 'logstep', 'wavelength_grid', 'freq_split',
           'find_library', 'find_header', 'open_library']

_EXP_DIR = Path(__file__).resolve().parents[1]
_STANDALONE = _EXP_DIR / 'HOP' / 'standalone'
_DEFAULT_LIB = _STANDALONE / 'build' / 'libhop.so'
_DEFAULT_HEADER = _STANDALONE / 'libhop.h'


class HopError(RuntimeError):
    """Nonzero status returned by the HOP C API."""


def find_library():
    """Locate libhop.so: $HOP_LIB if set, else HOP/standalone/build/libhop.so."""
    lib = os.environ.get('HOP_LIB', '')
    path = Path(lib) if lib else _DEFAULT_LIB
    if not path.is_file():
        raise HopError(f"libhop.so not found at {path}; build it with "
                       f"'make -C {_STANDALONE}' or point $HOP_LIB at it")
    return path


def find_header():
    """Locate libhop.h: $HOP_HEADER if set, else HOP/standalone/libhop.h.

    The header is generated from the BIND(C) declarations by the Makefile's ``$(HDR)``
    rule and committed, so it is the single source of truth for the ABI -- nothing here
    hand-writes a signature.
    """
    hdr = os.environ.get('HOP_HEADER', '')
    path = Path(hdr) if hdr else _DEFAULT_HEADER
    if not path.is_file():
        raise HopError(f"libhop.h not found at {path}; generate it with "
                       f"'make -C {_STANDALONE} bindings'")
    return path


def _c_str(s):
    """A NUL-terminated C string; cffi passes bytes straight through to ``char *``."""
    return str(s).encode()


def _opt_value(value):
    """Format an option value for hop_set_option, which takes everything as a string.

    Formatted so the Fortran ``read`` on the other side round-trips: ``repr`` for floats
    (shortest round-trippable form) and the lowercase words the Fortran accepts for bools.
    """
    if isinstance(value, (bool, np.bool_)):
        return 'true' if value else 'false'
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, str):
        return value
    return repr(float(value))


def _last_error(lib):
    """The library's message for the most recent nonzero status ('' if none)."""
    buf = _ffi.new('char[1024]')
    lib.hop_last_error(buf, 1024)
    return _ffi.string(buf).decode()


def _check(lib, status, what):
    """Raise with the library's own message, falling back to the bare status."""
    if status != 0:
        raise HopError(_last_error(lib) or f'{what} failed with status {status}')


def _verbosity(verbose):
    """HOP's integer verbosity from a bool: False = -1 (nothing at all), True = 1.

    Any other value passes through, so the dispatch levels (0 = the linelist
    line count only, 1 = progress messages, 2 = per-solve diagnostics) remain
    available as ints.
    """
    return (1 if verbose else -1) if isinstance(verbose, bool) else int(verbose)


#: cffi's FFI, built once per process by open_library() from the generated header.
_ffi = None


def open_library():
    """Open libhop.so against the generated header; returns ``(ffi, lib)``.

    The one place that knows how to load the library, so anything reaching into libhop
    directly (the hydrogen-Stark study scripts, say) binds against the same generated
    header as HopEOS rather than hand-writing a second set of signatures.

    RTLD_GLOBAL so the OpenMP runtime the library pulls in stays visible to it.
    """
    global _ffi
    if _ffi is None:
        ffi = cffi.FFI()
        ffi.cdef(find_header().read_text())
        _ffi = ffi
    return _ffi, _ffi.dlopen(str(find_library()), _ffi.RTLD_GLOBAL)


def _buf(arr, dtype, ctype):
    """A C pointer to a numpy array's data, dtype-checked.

    ``ffi.from_buffer`` reinterprets raw bytes and checks neither the dtype nor the
    layout, so the check lives here: every caller goes through ``np.ascontiguousarray``
    first, and this is what keeps a wrong-width array from being reinterpreted silently.
    """
    if arr.dtype != dtype or not arr.flags['C_CONTIGUOUS']:
        raise HopError(f'internal: expected a C-contiguous {np.dtype(dtype)} array, '
                       f'got {arr.dtype} (contiguous={arr.flags["C_CONTIGUOUS"]})')
    return _ffi.from_buffer(ctype, arr)


class HopEOS:
    """The HOP composition_type, driven from Python.

    Parameters mirror the &composition_params namelist; data-file paths default
    to the standard files under <experiment>/input/.  ``options`` is a dict of
    HOP module-level settings applied *before* the Fortran init (required for
    init-time options such as ``abundances`` or ``exclude_elements``); the
    recognized names are exactly what ``option_names()`` reports; the value is
    parsed by libhop against the option's own type.

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

        options = dict(options or {})
        # the VCS hydrogen-profile tables sit next to the other input data, but
        # HOP's default path is cwd-relative -- make it absolute unless overridden
        options.setdefault('vcs_data_dir', str(_EXP_DIR / 'input' / 'vcs'))

        # everything needed to build an identical instance elsewhere: opacity()
        # hands this to its spawned workers, and uses it to rebuild this instance
        # when a line list has to be replaced (see _load_linelist)
        self._recipe = dict(options,
                            absdat_file=str(absdat_file), abund_file=str(abund_file),
                            isotope_file=str(isotope_file), molec_file=str(molec_file),
                            feh=float(feh), tmolim=float(tmolim),
                            alpha_enhancement=float(alpha_enhancement),
                            r_process_enhancement=float(r_process_enhancement),
                            verbose=_verbosity(verbose))
        self._linelist = None   # (files, dlam, lam_min, lam_max) read into the instance
        self._state = None      # the load_state arguments, replayed after a rebuild
        self._shape = ()        # shape of the state, given to the field properties
        self._species = []      # keep_species ids, replayed after a rebuild

        self._lib = self._load_library()
        self._init_fortran()

    @staticmethod
    def _load_library():
        return open_library()[1]

    def _init_fortran(self):
        """Apply the recipe's options and run hop_init.  The Fortran side must be
        free (fresh instance, or right after hop_free in _rebuild)."""
        recipe = dict(self._recipe)
        files = {k: recipe.pop(k) for k in ('absdat_file', 'abund_file',
                                            'isotope_file', 'molec_file')}
        scalars = {k: recipe.pop(k) for k in ('feh', 'tmolim', 'alpha_enhancement',
                                              'r_process_enhancement')}
        for name, value in recipe.items():   # verbose included, before hop_init
            self.set_option(name, value)

        status = self._lib.hop_init(
            _c_str(files['absdat_file']), _c_str(files['abund_file']),
            _c_str(files['isotope_file']), _c_str(files['molec_file']),
            scalars['feh'], scalars['tmolim'],
            scalars['alpha_enhancement'], scalars['r_process_enhancement'])
        _check(self._lib, status, 'hop_init')
        HopEOS._open = True
        self._owns = True

    def _rebuild(self):
        """Start the Fortran side over from the recipe, keeping the kept-species marks.

        Line lists accumulate onto the never-solved master instance and cannot be
        removed, so replacing one means rebuilding -- see _load_linelist.

        The old line list is deliberately *not* replayed: _load_linelist is the only
        caller, and it reads the replacement (and re-solves the state) itself right
        after.  Replaying here would put the old list underneath the replacement and
        double-count every line the two windows share.  _state is left in place so
        _load_linelist can re-solve it once, after the new lines are in.
        """
        species = list(self._species)
        self._lib.hop_free()
        HopEOS._open = False
        self._owns = False
        self._linelist, self._species = None, []
        self._init_fortran()
        if species:
            self.keep_species(species)

    # -- configuration ------------------------------------------------------

    def set_option(self, name, value):
        """Set a HOP module-level option (namelist-style global).

        The value crosses as a string and libhop parses it against the option's own
        type, so there is nothing to guess here -- an unknown name and a value of the
        wrong type come back distinguishable, each with a message.  ``option_names()``
        lists what libhop accepts.
        """
        _check(self._lib,
               self._lib.hop_set_option(_c_str(name.lower()), _c_str(_opt_value(value))),
               f'set_option({name!r})')

    def option_names(self):
        """The option names libhop accepts, straight from the library."""
        buf = _ffi.new('char[4096]')
        self._lib.hop_option_names(buf, 4096)
        return _ffi.string(buf).decode().split()

    @staticmethod
    def _linelist_files(files):
        """The line-list argument as a tuple of file paths: a file, a list of
        files, or a directory (read in sorted order -- dispatch reads
        linelist_folder in the order os%get_filenames returns, so ties between
        identical lines may end up ordered differently)."""
        if isinstance(files, (str, Path)):
            p = Path(files)
            files = sorted(f for f in p.iterdir() if f.is_file()) if p.is_dir() else [p]
        out = []
        for f in files:
            f = Path(f)
            if not f.is_file():
                raise HopError(f'line-list file not found: {f}')
            out.append(str(f))
        return tuple(out)

    def _load_linelist(self, files, dlam, lam):
        """Make sure the line list (files, dlam) covering lam is the one read in.

        Only lines inside [lam[0]-dlam, lam[-1]+dlam] are kept, exactly as
        composition_pre_init does with the output spectrum's range.  A line list
        can only be *added* to the never-solved master instance -- and the solve
        works on a copy of it -- so a different list, or one read over too narrow
        a window, means rebuilding the instance and re-solving the state.
        """
        files = self._linelist_files(files)
        lam_min, lam_max = float(np.min(lam)), float(np.max(lam))
        if self._linelist is not None:
            have_files, have_dlam, have_min, have_max = self._linelist
            if (files, dlam) == (have_files, have_dlam) and \
                    have_min <= lam_min and lam_max <= have_max:
                return                       # already covered
            self._rebuild()                  # lines cannot be unread
        lam_arr = np.ascontiguousarray([lam_min, lam_max], np.float64)
        for f in files:
            status = self._lib.hop_init_linelist(
                _c_str(f), float(dlam), 2, _buf(lam_arr, np.float64, 'double *'))
            _check(self._lib, status, f'hop_init_linelist({f})')
        self._linelist = (files, dlam, lam_min, lam_max)
        if self._state is not None:
            self.load_state(**self._state)   # the solve must see the new lines

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
            _check(self._lib, status, f'keep_species({sid!r})')
            if sid not in self._species:
                self._species.append(sid)

    def species_density(self, species_id):
        """Number density [1/cm3] of one species ('Fe I', 'CH', ...) in the
        loaded state.  The species must have been kept before the load_state
        call (species= argument or keep_species()), or be one HOP retains
        anyway (continuum-opacity / line-list species)."""
        n = self.n
        if n == 0:
            raise HopError('no state loaded; call load_state first')
        out = np.empty(n, np.float64)
        status = self._lib.hop_get_species(_c_str(species_id), n,
                                           _buf(out, np.float64, 'double *'))
        _check(self._lib, status, f'species_density({species_id!r})')
        return out.reshape(self._shape)

    def load_state(self, temp, rho=None, ne=None, vmic=-1.0, species=None):
        """Solve the EOS on a set of cells.

        temp [K] plus either rho [g/cm3] (electron density solved by Newton
        iteration; pass ne as well to provide a starting guess) or ne [1/cm3]
        (rho follows from particle conservation).  The arrays may have any shape,
        as long as they share it -- the cells are solved independently, and the
        result properties come back in that same shape.  vmic [km/s] only matters
        for later line-opacity calls; < 0 means unset (thermal broadening only),
        and opacity(vmic=...) overrides it without re-solving.  species: id(s) to
        keep for species_density() ('Fe I', 'CH', ...).

        Returns self; results are exposed as the array properties temp, ne, rho,
        pe, pg, energy, eh, chi500.
        """
        if species is not None:
            self.keep_species(species)
        temp_a = np.ascontiguousarray(np.atleast_1d(temp), np.float64)
        shape = temp_a.shape
        n = temp_a.size
        use_rho, use_ne = rho is not None, ne is not None
        if not (use_rho or use_ne):
            raise HopError('load_state needs rho and/or ne')
        rho_in = np.ascontiguousarray(np.atleast_1d(rho if use_rho else np.zeros(shape)),
                                      np.float64)
        ne_io = np.ascontiguousarray(np.atleast_1d(ne if use_ne else np.zeros(shape)),
                                     np.float64)
        if rho_in.shape != shape or ne_io.shape != shape:
            raise HopError('temp, rho, ne must have the same shape')
        if not np.all(np.isfinite(temp_a)) or np.any(temp_a <= 0):
            raise HopError('temp must be finite and positive')

        # the Fortran side takes a flat cell list; C order throughout, so the
        # unflattening in _get_field is its exact inverse
        status = self._lib.hop_load_state(
            n, _buf(temp_a, np.float64, 'double *'), _buf(ne_io, np.float64, 'double *'),
            _buf(rho_in, np.float64, 'double *'), float(vmic),
            int(use_rho), int(use_ne))
        _check(self._lib, status, 'load_state')
        self._shape = shape
        self._state = dict(temp=temp_a, rho=rho_in if use_rho else None,
                           ne=ne_io if use_ne else None, vmic=vmic)
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
        status = self._lib.hop_get_field(_c_str(name), n,
                                         _buf(out, np.float64, 'double *'))
        _check(self._lib, status, f'hop_get_field({name})')
        return out.reshape(self._shape)

    temp = property(lambda self: self._get_field('temp'), doc='temperature [K]')
    ne = property(lambda self: self._get_field('ne'), doc='electron density [1/cm3]')
    rho = property(lambda self: self._get_field('rho'), doc='mass density [g/cm3]')
    pe = property(lambda self: self._get_field('pe'), doc='electron pressure [dyn/cm2]')
    pg = property(lambda self: self._get_field('pg'), doc='gas pressure [dyn/cm2]')
    energy = property(lambda self: self._get_field('e'), doc='internal energy [erg/g]')
    eh = property(lambda self: self._get_field('eh'),
                  doc='hydrogen ionization/dissociation energy term')
    chi500 = property(lambda self: self._get_field('chi500'),
                      doc='extinction at the reference wavelength (air 5000 A) [1/cm]')

    @property
    def grph(self):
        """Grams of matter per hydrogen nucleus."""
        return self._lib.hop_get_grph()

    # -- opacity ------------------------------------------------------------

    def _set_vmic(self, vmic):
        status = self._lib.hop_set_vmic(self.n, float(vmic))
        _check(self._lib, status, 'hop_set_vmic')
        if self._state is not None:
            self._state['vmic'] = float(vmic)

    def _opacity_call(self, lam, include_absmet, include_linelist, dlam,
                      line_scattering):
        """One hop_opacity call over the whole of lam, as flat (ncell, nlam)."""
        n = self.n
        chi = np.zeros((n, lam.size), np.float32)   # C order == Fortran (nlam, n)
        scat = np.zeros_like(chi)
        status = self._lib.hop_opacity(lam.size, _buf(lam, np.float64, 'double *'),
                                       _buf(chi, np.float32, 'float *'),
                                       _buf(scat, np.float32, 'float *'),
                                       int(bool(include_absmet)),
                                       int(bool(include_linelist)), float(dlam),
                                       int(bool(line_scattering)))
        _check(self._lib, status, 'hop_opacity')
        return chi, scat

    def scattering(self, lam):
        """Continuum scattering coefficient [1/cm]: Thomson + Rayleigh only.

        Shape ``<state shape> + (nlam,)``, float32, on vacuum wavelengths.  This is the
        scattering part of the continuum on its own -- ``opacity()``'s ``scat`` return
        adds the (1-eps)-weighted line terms on top when ``line_scattering`` is set.
        """
        lam = self._check_lam(lam)
        n = self.n
        if n == 0:
            raise HopError('no state loaded; call load_state first')
        sig = np.zeros((n, lam.size), np.float32)   # C order == Fortran (nlam, n)
        status = self._lib.hop_scattering(lam.size, _buf(lam, np.float64, 'double *'),
                                          _buf(sig, np.float32, 'float *'))
        _check(self._lib, status, 'hop_scattering')
        return sig.reshape(self._shape + (lam.size,))

    def opacity(self, lam, linelist=None, dlam=100.0, include_absmet=False,
                line_scattering=False, vmic=None, nchunks=None, nprocs=1):
        """Extinction chi and its scattering part [1/cm] on vacuum wavelengths.

        Both come back float32 with shape ``<state shape> + (nlam,)``, i.e. per
        volume -- divide by ``hop.rho[..., None]`` for the per-gram opacity a
        dispatch 'opa'/'scat' table holds.

        chi is the total: continuum (with Thomson + Rayleigh scattering folded in,
        as HOP's scattering_absorption does), the sampled metal-line background
        with ``include_absmet``, and the LTE line list with ``linelist``.  scat is
        the scattering part of it: Thomson + Rayleigh, plus -- with
        ``line_scattering`` -- the (1-eps)-weighted absmet and line-list opacity,
        the van Regemorter split composition_type%opacity applies.

        linelist: a file, a list of files or a folder (the &linelist_params
        linelist_file / linelist_folder), with ``dlam`` [AA] the line-wing search
        radius; only lines within dlam of the requested wavelengths are read.
        Switching to a different list rebuilds the instance and re-solves the
        state, since HOP can only ever add lines.

        vmic [km/s] broadens the lines and must match what dispatch would use, or
        the line opacity comes out wrong (an order of magnitude in the cool, dense
        corner): the &opac_table_params default is -1 (thermal only), but
        read_external_mod overrides it with the mean of the atmosphere's own vmic
        whenever any of that is positive -- so pass that mean here for an
        atmosphere that carries microturbulence.  It affects neither the EOS nor
        the continuum, and is applied without re-solving the state.

        The wavelength axis is split into ``nchunks`` pieces (default: ``nprocs``)
        the same way dispatch splits it over its n_nu rt_patches, and the pieces
        run in ``nprocs`` processes -- HOP keeps process-wide global state, so each
        worker holds its own instance and reads only the lines its chunk needs.
        With ``nprocs=1`` the chunks are computed here, in one instance holding the
        whole line list.  Because the workers are spawned, a script calling this
        with nprocs > 1 must guard its top level with ``if __name__ ==
        '__main__':`` -- without it each worker re-runs the script on import and
        the pool spawns forever (notebooks and pytest are unaffected).
        """
        lam = self._check_lam(lam)
        if self._state is None:
            raise HopError('no state loaded; call load_state first')
        if linelist is not None:
            linelist = self._linelist_files(linelist)
        nchunks = max(1, min(nchunks or nprocs, lam.size))
        nprocs = max(1, min(nprocs, nchunks))

        jobs = []
        for i in range(1, nchunks + 1):
            fs, fe = freq_split(lam.size, i, nchunks)
            jobs.append(dict(lam=lam[fs:fe], recipe=self._recipe, state=self._state,
                             species=list(self._species), linelist=linelist, dlam=dlam,
                             vmic=vmic, include_absmet=include_absmet,
                             line_scattering=line_scattering))

        if nprocs > 1:
            # spawn, not fork: a forked child would inherit the parent's Fortran
            # globals (and its live instance, which HopEOS allows only one of)
            ctx = multiprocessing.get_context('spawn')
            pool = ctx.Pool(nprocs)
            try:
                results = pool.map(_opacity_chunk, jobs, chunksize=1)
                # close/join, not the with-block's terminate(): the workers have
                # to exit normally for the Fortran side to flush what it wrote
                pool.close()
                pool.join()
            finally:
                pool.terminate()
        else:
            if linelist is not None:
                self._load_linelist(linelist, dlam, lam)
            if vmic is not None:
                self._set_vmic(vmic)
            results = [_opacity_chunk(job, self) for job in jobs]

        shape = self._shape + (lam.size,)
        return tuple(np.concatenate([r[k] for r in results], axis=1).reshape(shape)
                     for k in (0, 1))

    def _check_lam(self, lam_vac):
        lam = np.ascontiguousarray(np.atleast_1d(lam_vac), np.float64)
        if lam.ndim != 1 or lam.size < 1:
            raise HopError('lam_vac must be a 1D array')
        if lam.size > 1 and not np.all(np.diff(lam) > 0):
            raise HopError('lam_vac must be strictly ascending (vacuum Angstrom)')
        return lam

    # -- lifetime -----------------------------------------------------------

    def close(self):
        # only the instance that actually owns the Fortran state may free it:
        # a closed instance being garbage-collected must not tear down the
        # successor that has meanwhile taken over the process-wide globals
        if self._owns:
            self._owns = False
            HopEOS._open = False
            self._linelist, self._state, self._shape = None, None, ()
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


def _opacity_chunk(job, hop=None):
    """One wavelength chunk of HopEOS.opacity.

    Module level and argument-only so it survives pickling to a spawned worker,
    which calls it with hop=None and so opens (and closes) its own instance, built
    from the parent's recipe and given the parent's state.  In the sequential path
    the caller passes its own instance -- including its line list, read once over
    the full wavelength range rather than per chunk.
    """
    os.environ.setdefault('OMP_NUM_THREADS', '1')  # before the library is loaded
    own = hop is None
    hop = HopEOS(**job['recipe']) if own else hop
    try:
        if own:
            if job['linelist'] is not None:
                hop._load_linelist(job['linelist'], job['dlam'], job['lam'])
            if job['species']:
                hop.keep_species(job['species'])
            hop.load_state(**job['state'])
            if job['vmic'] is not None:
                hop._set_vmic(job['vmic'])
        return hop._opacity_call(job['lam'], job['include_absmet'],
                                 job['linelist'] is not None, job['dlam'],
                                 job['line_scattering'])
    finally:
        if own:
            hop.close()


def logstep(start, stop, step):
    """Like np.logspace but with a ln-space step size rather than a count."""
    lmin = np.log(start)
    lmax = np.log(stop)
    n = int(np.ceil((lmax - lmin) / step))
    return np.exp(lmin + step * np.arange(n))


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
