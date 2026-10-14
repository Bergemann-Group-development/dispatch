"""Drive MULTI 2.3 (1D) from tumult and read its results as tumult spectra.

`Multi1D.from_tumult_run(run)` sets up and runs MULTI for the atom and atmosphere of a 1D
OpenMULTI3D run, so the two codes can be compared without hand-matching inputs; the
`multi1d` preset makes OpenMULTI3D reproduce MULTI's frequency quadrature, angles and
composition. `Multi1D.read(folder)` reads such a result again without rerunning MULTI.

Point `MULTI1D_PATH` at the MULTI tree, either by setting the environment variable or by
assigning the module constant. The executable is `source_dp/mul23lus.x` and the default
`ABSDAT`, `ABSMET` and `ABUND` files are in `input/` there.

MULTI writes one wavelength grid per line, so each line is its own `Spectrum` on its own
grid (`Multi1D.lines`, `Multi1D.get_line`). The windows of a multiplet overlap, which is
why they are not concatenated into one segmented spectrum. MULTI works in frequency on
vacuum wavelengths; the spectra here are on the air axis of an OpenMULTI3D run, and
absolute fluxes are converted from per Hz to per Å. MULTI defines the flux as
2π Σ w_μ μ I(μ), as `SRC/rt_patch_mod.f90` does, so absolute fluxes compare directly. A
normalised spectrum is divided by MULTI's own continuum, the background-only flux at line
centre (`FLUX(0)`), which is also what MULTI's equivalent widths use.

MULTI computes the LTE profiles of an NLTE atom only in a separate run with `ITMAX=0`;
`from_tumult_run` makes both, in `NLTE/` and `LTE/`, and `crop(LTE=True)` reads the
second.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from scipy.integrate import trapezoid

from .auxfuncs import marcs_name_to_params
from .const import CC
from .physics import vac2air
from .spectrum import Spectrum

MULTI1D_PATH = os.environ.get('MULTI1D_PATH')

# MULTI's INPUT for a comparison run; from_tumult_run overrides the iteration and
# convergence keys from the run's namelist
DEFAULT_INPUT = {
    'DIFF': 5.0, 'ELIM1': 0.01, 'ELIM2': 0.001, 'QNORM': 1.00, 'THIN': 0.1,
    'IATOM2': 0, 'IHSE': 0, 'IOPAC': 1, 'ISTART': 1, 'ISUM': 0, 'ITRAN': 0, 'NMU': 3,
    'ICONV': 1, 'ILAMBD': 0, 'ITMAX': 0,
    'IWABND': 0, 'IWATMS': 0, 'IWATOM': 0, 'IWCHAN': 0, 'IWDAMP': 0, 'IWEMAX': 1,
    'IWEQW': 0, 'IWEVEC': 0, 'IWHEAD': 0, 'IWHSE': 0, 'IWLGMX': 1, 'IWLINE': 0,
    'IWLTE': 0, 'IWN': 0, 'IWNIIT': 0, 'IWOPAC': 0, 'IWRAD': 0, 'IWRATE': 0, 'IWSTRT': 0,
    'IWTAUQ': 0, 'IWTEST': 0, 'IWWMAT': 0, 'IWJFIX': 0, 'IWARN': 2, 'IOPACL': 0,
    'ISCAT': 0, 'INCRAD': 0, 'INGACC': 1, 'ICRSW': 0, 'IOSMET': 0, 'EOSMET': 0,
    'IDL1': 1, 'IDLNY': 0, 'IDLCNT': 0, 'IDLOPC': 0,
}

# scratch files of a finished run, about 2 GB for a large atom; the reader needs none
SCRATCH = ['OPC', 'PHI', 'JNY', 'DUMS', 'DUMS5', 'DUMI', 'SWIDE', 'XWIDE']

FLUX_UNIT = u.erg / (u.s * u.cm**2 * u.AA)


#=====================================================================
class _Records:
    """Sequential reader of a little-endian Fortran unformatted file (4-byte markers)."""
    def __init__(self, path):
        self.file = open(path, 'rb')

    def size(self):
        """Byte length of the next record, without consuming it."""
        pos = self.file.tell()
        size = int(np.frombuffer(self.file.read(4), '<i4')[0])
        self.file.seek(pos)
        return size

    def read(self, dtype):
        size = int(np.frombuffer(self.file.read(4), '<i4')[0])
        data = np.frombuffer(self.file.read(size), dtype)
        self.file.read(4)
        return data

    def ints(self):
        return self.read('<i4')

    def reals(self):
        return self.read('<f4').astype(np.float64)

    def text(self):
        return self.read('u1').tobytes().decode('ascii', errors='replace')

    def skip(self, n=1):
        for _ in range(n):
            size = int(np.frombuffer(self.file.read(4), '<i4')[0])
            self.file.seek(size + 4, 1)

    def close(self):
        self.file.close()


#=====================================================================
def read_idl1(path):
    """The quantities of a MULTI ``IDL1`` file (``mul23_writeidl.f``, ``IFMT=0``) that a
    comparison needs, as a dict of arrays; per-level arrays are ``(nk, ndep)`` and line
    arrays ``(..., nrad)``, in Fortran order.

    Two writers are in use. One writes ``ITMAX`` after ``NSTAR`` and the NLTE
    populations only when ``ITMAX > 0``; an older one writes no ``ITMAX`` and always
    writes the populations. The record after ``NSTAR`` tells them apart: an ``ITMAX``
    record holds one integer. Without it, the run counts as LTE when its populations
    equal the LTE ones.
    """
    rec = _Records(path)
    try:
        ndep, nk, nline, nwide, nrad, nrfix, nmu, mq = rec.ints()
        if ndep == 0:
            raise ValueError(f'{path} is a compact IDL1 file (IFMT > 0), which carries '
                             'no profiles of unselected lines; rerun MULTI with IDL1=1')
        out = dict(ndep=ndep, nk=nk, nline=nline, nrad=nrad, nmu=nmu, mq=mq)
        out['nq'] = rec.ints() if nrad > 0 else np.zeros(0, int)
        out['qnorm'] = rec.reals()[0]
        out['abnd'], out['awgt'] = rec.reals()
        out['ev'] = rec.reals()
        out['g'] = rec.reals()
        out['ion'] = rec.ints()
        rec.skip()                                           # HN3C2
        if nrad > 0:
            rec.skip()                                       # KTRANS
            out['jrad'] = rec.ints()
            out['irad'] = rec.ints()
            out['f'] = rec.reals()
            rec.skip()                                       # IWIDE
            out['ga'] = rec.reals()
            out['gw'] = rec.reals()
            out['gq'] = rec.reals()
        rec.skip(2)                                          # KRAD, Z
        if nwide > 0:
            rec.skip()                                       # ALFAC
        rec.skip()                                           # HNY4P
        if nrad > 0:
            out['alamb'] = rec.reals()
        if nline > 0:
            rec.skip()                                       # A
        rec.skip(2)                                          # B, TOTN
        if nrad > 0:
            rec.skip()                                       # BP
        out['nstar'] = rec.reals().reshape((nk, ndep), order='F')
        if rec.size() == 4:
            itmax = int(rec.ints()[0])
            if itmax > 0:
                out['n'] = rec.reals().reshape((nk, ndep), order='F')
                rec.skip()                                   # C
            else:
                out['n'] = out['nstar']
        else:
            out['n'] = rec.reals().reshape((nk, ndep), order='F')
            rec.skip()                                       # C
            itmax = 0 if np.array_equal(out['n'], out['nstar']) else None
        out['mode'] = 'LTE' if itmax == 0 else 'NLTE'
        if nrfix > 0:
            rec.skip(6)                                      # JFX ... ITRAD
        rec.skip()                                           # DNYD
        if nline > 0:
            rec.skip()                                       # ADAMP
        label = rec.text()
        width = len(label) // nk
        out['label'] = [label[i*width:(i+1)*width].strip() for i in range(nk)]
        out['atomid'] = rec.text().strip()
        rec.skip()                                           # CROUT
        out['grav'] = rec.reals()[0]
        rec.skip()                                           # CMASS
        out['temp'] = rec.reals()
        out['ne'] = rec.reals()
        rec.skip()                                           # VEL
        out['tau'] = rec.reals()
        rec.skip(2)                                          # XNORM, HEIGHT
        out['atmosid'] = rec.text().split()[0]
        out['vturb'] = rec.reals()
        rec.skip(2)                                          # BH, NH
        out['rho'] = rec.reals()
        if nrad > 0:
            out['qmax'] = rec.reals()
            out['q0'] = rec.reals()
            out['ind'] = rec.ints()
            rec.skip()                                       # DIFF
            out['q'] = rec.reals().reshape((mq, nrad), order='F')
            rec.skip()                                       # WQ
        rec.skip()                                           # WQMU
        if nwide > 0:
            rec.skip()                                       # FRQ
        if nrad > 0:
            rec.skip(2)                                      # WPHI, SL
        if nline > 0:
            out['weqlte'] = rec.reals()
            out['weq'] = rec.reals()
        if nrad > 0:
            rec.skip(2)                                      # RIJ, RJI
            out['flux'] = rec.reals().reshape((mq + 1, nrad), order='F')
            out['outint'] = rec.reals().reshape((mq + 1, nmu, nrad), order='F')
            rec.skip()                                       # COOL
        out['xmu'] = rec.reals()
        out['wmu'] = rec.reals()
    finally:
        rec.close()
    return out


#=====================================================================
class Multi1DLine:
    """Bound-bound transition ``kr`` (from 0) of a MULTI run, on its own wavelength grid.

    ``i`` and ``j`` are the lower and upper level indices from 0, ``lam0`` the air and
    ``lam0_vac`` the vacuum wavelength [Å]. ``lam`` is the ascending air grid; ``flux``
    [erg/s/cm²/Å], ``ie`` [erg/s/cm²/Å/sr, (nlam, nmu)] and their continua ``flux_cnt``
    and ``ie_cnt`` are on it.
    """
    def __init__(self, parent, kr, data):
        self.parent = parent
        self.kr = kr
        self.i = int(data['irad'][kr]) - 1
        self.j = int(data['jrad'][kr]) - 1
        self.lam0_vac = float(data['alamb'][kr])
        self.lam0 = float(vac2air(np.array([self.lam0_vac]))[0])
        self.loggf = float(np.log10(data['f'][kr] * data['g'][self.i]))
        self.ga, self.gw, self.gq = (float(data[k][kr]) for k in ('ga', 'gw', 'gq'))
        self.qmax, self.q0 = float(data['qmax'][kr]), float(data['q0'][kr])

        nq = int(data['nq'][kr])
        qq = data['q'][:nq, kr]
        flux = data['flux'][1:nq + 1, kr]
        ie = data['outint'][1:nq + 1, :, kr]
        if data['ind'][kr] != 2:
            # a symmetric profile stores one half, from line centre outwards
            qq = np.concatenate([-qq[:0:-1], qq])
            flux = np.concatenate([flux[:0:-1], flux])
            ie = np.concatenate([ie[:0:-1], ie])

        # q is the frequency offset in units of qnorm: nu = nu0 (1 + q qnorm/c)
        qn = data['qnorm'] * 1e5 / CC
        lam_vac = self.lam0_vac / (1 + qq * qn)
        order = np.argsort(lam_vac)
        lam_vac = lam_vac[order]
        self.lam_vac = lam_vac
        self.lam = vac2air(lam_vac)

        # F_lambda = F_nu c / lambda^2, with lambda in vacuum
        jac = CC * 1e8 / lam_vac**2
        self.flux = flux[order] * jac
        self.ie = ie[order] * jac[:, None]
        self.flux_cnt = data['flux'][0, kr] * jac
        self.ie_cnt = data['outint'][0, :, kr][None, :] * jac[:, None]
        if 'weq' in data:
            # MULTI's own equivalent widths [Å]
            self.ew = float(data['weq'][kr])
            self.ew_lte = float(data['weqlte'][kr])

    @property
    def title(self):
        return 'Line = {:.2f} [Å]'.format(self.lam0)

    #=====================================================================
    def get_spectrum(self, imu=None, cont=False):
        """Flux, or the intensity at ``xmu[imu]``, of this run (not its LTE twin)."""
        if imu is None:
            yy = self.flux_cnt if cont else self.flux
            yunit = FLUX_UNIT
        else:
            yy = (self.ie_cnt if cont else self.ie)[:, imu]
            yunit = FLUX_UNIT / u.sr
        return Spectrum(self.lam, yy, xunit=u.AA, yunit=yunit)

    #=====================================================================
    def crop(self, imu=None, norm=True, cont=False, LTE=False, **kwargs):
        """The line's spectrum cropped by :meth:`Spectrum.crop`, centred on ``lam0`` as
        :meth:`tumult.run.Run.crop` centres a line; ``norm`` divides by the continuum.
        ``LTE`` takes the line from the LTE run of an NLTE comparison."""
        line = self
        if LTE and self.parent.mode != 'LTE':
            if self.parent.lte is None:
                raise ValueError(f'{self.parent.folder} has no LTE run')
            line = self.parent.lte.lines[self.kr]
        kwargs.setdefault('x0', line.lam0)
        spec = line.get_spectrum(imu=imu, cont=cont).crop(**kwargs)
        if norm:
            spec = spec / line.get_spectrum(imu=imu, cont=True).crop(**kwargs)
        return spec


#=====================================================================
class Multi1D:
    """A MULTI 2.3 run (``IDL1``), with its LTE twin in ``lte`` when it is NLTE."""
    lte = None

    @classmethod
    #=================================================================
    def read(cls, folder):
        """A comparison folder written by :meth:`from_tumult_run` (``NLTE/`` and
        ``LTE/``, or ``LTE/`` alone), or a single MULTI run directory."""
        folder = Path(folder)
        if (folder / 'IDL1').exists():
            return cls._read_run(folder)
        if (folder / 'NLTE' / 'IDL1').exists():
            self = cls._read_run(folder / 'NLTE')
            if (folder / 'LTE' / 'IDL1').exists():
                self.lte = cls._read_run(folder / 'LTE')
            self.folder = folder
            return self
        if (folder / 'LTE' / 'IDL1').exists():
            self = cls._read_run(folder / 'LTE')
            self.folder = folder
            return self
        raise FileNotFoundError(f'no IDL1 in {folder}, {folder}/NLTE or {folder}/LTE')

    @classmethod
    def _read_run(cls, folder):
        self = cls()
        self.folder = folder
        data = read_idl1(folder / 'IDL1')
        self.data = data
        self.mode = data['mode']
        self.atomid = data['atomid']
        self.abnd = data['abnd']
        self.mus = data['xmu']
        self.ltau = np.log10(data['tau'])
        self.lines = [Multi1DLine(self, kr, data) for kr in range(data['nline'])]
        return self

    #=====================================================================
    def get_line(self, lam0, tol=0.05, verbose=True):
        """The line nearest to ``lam0`` [Å, air]; ``verbose`` warns when it is further
        than ``tol`` [Å] away."""
        lam0s = np.array([l.lam0 for l in self.lines])
        line = self.lines[int(np.argmin(np.abs(lam0s - lam0)))]
        if verbose and abs(line.lam0 - lam0) > tol:
            print(f'WARNING: Line not found! Closest line to {lam0:.2f}Å found at {line.lam0:.2f}Å')
        return line

    #=====================================================================
    def get_depart(self, ilvl=None):
        """Departure coefficients ``n/n*`` against log tau_500, as ``(ltau, dep)``:
        ``dep`` is ``(ndep, nk)``, or ``(ndep,)`` for level ``ilvl`` (from 0), the layout
        of :meth:`tumult.run.Run.get_depart` for a 1D run."""
        dep = (self.data['n'] / self.data['nstar']).T
        return self.ltau, dep if ilvl is None else dep[:, ilvl]

    #=====================================================================
    def compare(self, run):
        """Line by line against an OpenMULTI3D ``run`` of the same atom, as a DataFrame.

        Both normalised flux profiles are taken on MULTI's grid, OpenMULTI3D's
        interpolated onto it, over the part of MULTI's window that the run's line window
        covers: ``ew_*`` [mÅ] integrate 1 - F there and ``max_dflux`` is the largest
        difference of the normalised fluxes. A run with a line list has blends in its
        windows that MULTI does not compute; compare such runs line by line with care.
        """
        atom = self._run_atom(run)
        rows = []
        for line in self.lines:
            tran = atom.bbtrans[line.kr]
            if abs(tran.lam0 - line.lam0) > 0.05:
                raise ValueError(f'line {line.kr} is at {line.lam0:.3f} Å in MULTI and '
                                 f'{tran.lam0:.3f} Å in the run: not the same atom')
            row = {'kr': line.kr, 'lam0': line.lam0, 'i': line.i, 'j': line.j}
            modes = [True] if run.mode == 'LTE' else [False, True]
            for lte in modes:
                tag = '_lte' if lte else ''
                syn = line.crop(LTE=lte)
                om3d = tran.crop(LTE=lte)
                xx, yy = syn.xdata.value, syn.ydata.value
                ox, oy = om3d.xdata.value, om3d.ydata.value
                mask = (xx >= ox.min()) & (xx <= ox.max())
                oy = np.interp(xx[mask], ox, oy)
                row[f'ew_m1d{tag}'] = trapezoid(1 - yy[mask], xx[mask]) * 1e3
                row[f'ew_om3d{tag}'] = trapezoid(1 - oy, xx[mask]) * 1e3
                row[f'max_dflux{tag}'] = np.max(np.abs(oy - yy[mask]))
            rows.append(row)
        return pd.DataFrame(rows)

    @staticmethod
    def _run_atom(run):
        if len(run.atoms) != 1:
            raise ValueError(f'MULTI computes one atom; the run has {len(run.atoms)}')
        return run.atoms[0]

    @classmethod
    #=================================================================
    def from_tumult_run(cls, run, folder=None, atom_file=None, logg=None, settings=None,
                        overwrite=False, keep_scratch=False, verbose=False):
        """Run MULTI for the atom, atmosphere and abundances of the 1D OpenMULTI3D
        ``run`` and read the result.

        The atmosphere is the run's own (``Atmos.to_multi1d_atmos``) and the abundances
        are the run's, written into ``ATMOS``, ``DSCALE`` and ``ABUND`` and into the
        abundance line of the copied ``ATOM``. ``atom_file`` is the classic MULTI atom;
        by default the run's namelist entry, relative to the working directory. ``logg``
        defaults to the value in a MARCS file name. ``settings`` overrides keys of the
        MULTI ``INPUT`` (:data:`DEFAULT_INPUT`). An NLTE run gets an NLTE and an LTE
        MULTI run in ``folder/NLTE`` and ``folder/LTE``, an LTE run only the second;
        ``folder`` defaults to ``<run folder>.m1d`` beside the run.
        """
        if run.atmos.dim != '1D':
            raise ValueError('MULTI is a 1D code; the run is on a 3D atmosphere')
        atom = cls._run_atom(run)
        nml = run.nml
        folder = Path(folder) if folder else run.folder.with_name(run.folder.name + '.m1d')
        if folder.exists():
            if not overwrite:
                raise FileExistsError(f'{folder} exists; read it with Multi1D.read or '
                                      'pass overwrite=True')
            shutil.rmtree(folder)

        atom_file = Path(atom_file or nml['atom_params']['atom_file'].strip())
        if not atom_file.is_file():
            raise FileNotFoundError(f'MULTI atom {atom_file} not found; pass atom_file=')
        atmos_file = Path(nml['atmos_params']['atmos_file'].strip())
        if logg is None:
            try:
                logg = marcs_name_to_params(atmos_file.name)[1]
            except ValueError:
                raise ValueError(f'cannot read logg from {atmos_file.name}; pass logg=') from None

        cfg = dict(DEFAULT_INPUT)
        cfg['ITMAX'] = nml['rt_iteration_params']['maxiter']
        cfg['ILAMBD'] = nml['rt_output_params']['maxiter']
        cfg['ELIM1'] = cfg['ELIM2'] = nml['nlte_params']['convlim']
        cfg['QNORM'] = nml['atom_params']['qnorm']
        absmet_file = nml['opac_table_params']['absmet_file'].strip()
        cfg['IOSMET'] = 1 if absmet_file else 0
        cfg.update(settings or {})

        runs = {'LTE': {**cfg, 'ITMAX': 0}}
        if run.mode == 'NLTE':
            runs = {'NLTE': cfg, **runs}
        for name, run_cfg in runs.items():
            run_dir = folder / name
            cls._prepare(run_dir, run, atom, atom_file, logg, run_cfg,
                         Path(absmet_file) if absmet_file else None)
            cls.launch(run_dir, keep_scratch=keep_scratch, verbose=verbose)
        return cls.read(folder)

    @staticmethod
    def _prepare(run_dir, run, atom, atom_file, logg, cfg, absmet_file):
        top = _multi1d_path()
        run_dir.mkdir(parents=True)
        run.atmos.to_multi1d_atmos(run_dir, 'model', logg=logg, abund=run.abundances)
        for src, dst in (('atmos', 'ATMOS'), ('dscale', 'DSCALE'), ('abund', 'ABUND')):
            (run_dir / f'{src}.model').rename(run_dir / dst)

        # the abundance is the first number of the atom file, after its id line
        text = atom_file.read_text()
        text = re.sub(r'(^\s*)([0-9]*\.[0-9]*)', rf'\g<1>{run.abundances[atom.id]}',
                      text, count=1, flags=re.MULTILINE)
        (run_dir / 'ATOM').write_text(text)

        absmet_file = absmet_file.resolve() if absmet_file else top / 'input' / 'absmet'
        (run_dir / 'ABSMET').symlink_to(absmet_file)
        (run_dir / 'ABSDAT').symlink_to(top / 'input' / 'absdat')
        (run_dir / 'INPUT').write_text(',\n'.join(f'{k}={v:G}' if isinstance(v, float)
                                                  else f'{k}={v}' for k, v in cfg.items()))

    @staticmethod
    #=================================================================
    def launch(run_dir, keep_scratch=False, verbose=False):
        """Run MULTI in a prepared ``run_dir`` (``ATOM``, ``ATMOS``, ``DSCALE``,
        ``ABUND``, ``ABSMET``, ``ABSDAT``, ``INPUT``), logging to ``JOBLOG``. Scratch
        files are removed afterwards unless ``keep_scratch``."""
        run_dir = Path(run_dir)
        exe = _multi1d_path() / 'source_dp' / 'mul23lus.x'
        (run_dir / 'IDL1').unlink(missing_ok=True)
        log_file = run_dir / 'JOBLOG'
        with open(log_file, 'w') as log:
            proc = subprocess.run([str(exe)], stdout=log, stderr=subprocess.STDOUT,
                                  cwd=run_dir)
        log = log_file.read_text()
        if verbose:
            print(log)
        idl1 = run_dir / 'IDL1'
        if proc.returncode != 0 or not idl1.exists() or idl1.stat().st_size == 0:
            tail = '\n'.join(log.splitlines()[-15:])
            raise RuntimeError(f'MULTI failed in {run_dir} (exit {proc.returncode}). '
                               f'Last JOBLOG lines:\n{tail}')
        if not keep_scratch:
            for name in SCRATCH:
                (run_dir / name).unlink(missing_ok=True)


def _multi1d_path():
    if not MULTI1D_PATH:
        raise RuntimeError('set MULTI1D_PATH (environment or tumult.multi1d.MULTI1D_PATH) '
                           'to the MULTI 2.3 tree')
    return Path(MULTI1D_PATH).expanduser()
