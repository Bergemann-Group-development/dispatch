"""Drive Korg.jl from tumult, with the interface `turbospectrum.py` uses.

`Korg.from_tumult_run(run)` takes an OpenMULTI3D run and synthesises the same window
with the same transitions, so the two codes can be compared without hand-matching
inputs. `launch` accepts a keyword dict the way the Turbospectrum wrapper accepts a
babsma/bsyn input block; the keys here are Korg's, and `korg_driver.jl` translates them
into a `Korg.synthesize` call.

Korg is a Julia package, so a run means starting julia on the driver script. Point
`KORG_PATH` at the Korg source tree and `JULIA_BIN` at the julia binary, either by
setting the environment variables or by assigning the module constants. Starting julia
costs a few seconds of package load and compilation per process, which dwarfs a small
synthesis; prefer one wide window over many narrow ones, and `nthreads` only for
windows large enough to repay it.

Wavelengths are air on the tumult side and vacuum inside Korg. The driver converts the
window bounds and every line once, and writes the spectrum back on the air axis, so
`lam` here lines up sample-for-sample with an OpenMULTI3D run over the same window.
"""
from .spectrum import Spectrum
from .auxfuncs import (atomic_weights, get_particle_weight, molecule_dict,
                       near_equal_split, tmp_txtfile)
from .const import BK
from .atmos import AtmosClass
from astropy import units as u
import numpy as np
import pandas as pd
import subprocess
import os
import math
from pathlib import Path
import multiprocessing as mp
from contextlib import ExitStack
from functools import partial
from io import StringIO
import pickle

KORG_PATH = os.environ.get('KORG_PATH', '~/Documents/Synthesis Codes/Korg.jl-main')
JULIA_BIN = os.environ.get('JULIA_BIN', '~/.juliaup/bin/julia')
DRIVER = Path(__file__).with_name('korg_driver.jl')

# HOP/const_mod.f90, so the broadening recipes below use HOP's own constants
EM, UU, MH, PI = 9.109534e-28, 1.6605655e-24, 1.6726e-24, math.pi
C19 = math.log10(39.0 * (8.0 / (PI * EM * BK ** 5)) ** (1.0 / 6.0))   # Gray p.237, ~19.4

# the atmosphere dump the driver writes to stdout, one row per layer
ATMOS_KEY = 'tau,T,ne,nt,z'


#=====================================================================
def vbar2(temp, awgt):
    """HOP ``ion_mod.f90 compute_vbar2`` -- mean *square* relative speed, cm^2/s^2."""
    return (8.0 * BK / (PI * UU)) * temp * (UU / MH + 1.0 / awgt)


#=====================================================================
def convert_broadening(ion, gvdw, grad, gstk):
    """One line list record's ``(gvdw, grad, gstk)`` as Korg's ``(gamma_rad, gamma_stark, vdW)``.

    HOP and Korg both descend from Turbospectrum and share the ABO packing, the
    T^(1/6) Stark scaling and the T^0.3 van der Waals scaling about 10,000 K, so most
    codes pass through untouched. The arithmetic below mirrors ``prepare_stark`` and
    ``prepare_vanderwaals`` in ``HOP/linelist_mod.f90``. NaN asks Korg to apply its own
    approximation, which is what HOP does for the corresponding zero code.

    One code does not map exactly. ``-20 < gvdw < 0`` is Turbospectrum's 'W' recipe, which
    HOP makes temperature independent by dividing out the ``vbar2**0.3`` it applies
    later, while Korg scales it as T^0.3; the value returned matches HOP at 10,000 K and
    drifts from it elsewhere. It is rare: in a 200 A extract of `master_line.txt` it
    covers 1 of 6168 lines.
    """
    is_mol = ion in molecule_dict
    awgt = get_particle_weight(ion if is_mol else ion.split()[0])

    if gstk < -10.0:                       # recipe c4 SIU: gstk is log10(C4)
        gamma_stark = BK * 10.0 ** (C19 + (2.0 / 3.0) * gstk) * 1e4 ** (1.0 / 6.0)
    elif gstk < 0.0:                       # log10(gamma_stark) per electron at 1e4 K
        gamma_stark = 10.0 ** gstk
    elif gstk == 0.0:
        gamma_stark = float('nan')       # HOP uses its Cowley-like branch, Korg its own
    else:
        gamma_stark = gstk

    if gvdw >= 20.0:                       # ABO, packed identically in both codes
        vdw = gvdw
    elif is_mol:                         # HOP applies no Unsold to molecules
        vdw = 0.0
    elif gvdw < -20.0:                     # recipe c6 SIU: gvdw is log10(C6)
        vdw = math.log10(17.0) + 0.4 * gvdw + 0.3 * math.log10(vbar2(1e4, awgt))
    elif gvdw < 0.0:                       # Turbospectrum 'W'; matched at 1e4 K only
        vdw = gvdw - 0.3 * math.log10(vbar2(1e4, awgt))
    elif gvdw == 0.0:
        vdw = 0.0
    else:                                # 0 < gvdw < 20: Unsold fudge factor, verbatim
        vdw = gvdw

    # grad passes through as HOP reads it: prepare_radiative treats a non-zero value as a
    # rate in s^-1 and falls back to the classical width when it is zero.
    return (grad if grad != 0.0 else float('nan')), gamma_stark, vdw


#=====================================================================
def write_linelist(llist, filename):
    """Write a tumult `LineList` as the CSV `korg_driver.jl` reads.

    Columns are ``lam0, ion, evi, loggf, gamma_rad, gamma_stark, vdw``, sorted by
    wavelength because Korg requires a sorted list. Wavelengths stay in air; the driver
    converts them.
    """
    df = llist.copy()
    if isinstance(df.index, pd.MultiIndex) or 'ion' not in df.columns:
        df = df.reset_index()
    for col in ('gvdw', 'grad', 'gstk'):
        if col not in df.columns:
            df[col] = 0.0

    broadening = [convert_broadening(row.ion, row.gvdw, row.grad, row.gstk)
                  for row in df.itertuples()]
    out = pd.DataFrame({'lam0': df.lam0.values,
                        'ion': df.ion.values,
                        'evi': df.evi.values,
                        'loggf': df.loggf.values})
    out[['gamma_rad', 'gamma_stark', 'vdw']] = np.asarray(broadening, dtype=float)
    out = out.sort_values('lam0')
    out.to_csv(filename, index=False)
    return out


#=====================================================================
from ._spectrum_table import SpectrumTable


class Korg(SpectrumTable):
    @property
    #=================================================================
    def atmos(self):
        if not hasattr(self, '_atmos'):
            self._atmos = self.read_atmos()
        return self._atmos

    @classmethod
    #=================================================================
    def read(cls, filename):
        self = cls()
        self.filename = filename

        with open(filename, 'r') as f:
            line = f.readline()
            self.mus = [float(s) for s in line.strip('# mu-points').split()]
            df = pd.read_csv(f, sep=r'\s+', header=None)
            self.df = df.apply(pd.to_numeric, errors='coerce')
            self.lam = self.df.iloc[:, 0].values

        return self

    @classmethod
    #=================================================================
    def launch(cls, input, llist=None, cwd=None, settings=None, verbose=False):
        cwd = Path(cwd or KORG_PATH).expanduser()

        with ExitStack() as stack:
            if settings is not None:
                input.update(settings)

            # an already-set LINELIST (e.g. from from_tumult_run, which converts the
            # list once up front) is used as-is; only convert here when a llist object
            # is passed and no linelist file has been prepared yet
            if llist is not None and 'LINELIST' not in input:
                input['LINELIST'] = cls._write_linelist_file(
                        stack, llist, input['LAMBDA_MIN'], input['LAMBDA_MAX'])

            # julia runs with cwd=KORG_PATH, so anything relative in the input would
            # resolve against the Korg tree rather than the caller's directory. Make
            # the inputs absolute here (relative to where python is running) and, like
            # the Turbospectrum wrapper, place a relative result inside cwd.
            for key in ('ATMOSPHERE', 'LINELIST'):
                if key in input:
                    input[key] = str(Path(input[key]).expanduser().resolve())
            result_path = Path(input['RESULTFILE']).expanduser()
            if not result_path.is_absolute():
                result_path = cwd / result_path
            input['RESULTFILE'] = str(result_path)

            # Remove a previous result *before* running, so that a silent failure can
            # never be mistaken for a fresh spectrum -- cls.read below would otherwise
            # happily pick up the stale file.
            result_path.unlink(missing_ok=True)

            input_path = tmp_txtfile(stack)
            cls._write_input_file(input, input_path)
            stdout, stderr = run_julia(DRIVER, input_path, cwd=cwd, verbose=verbose)

            # julia exits 0 on some errors it catches and reports, so check that the
            # spectrum is actually there and carries data before trusting it
            if not result_path.exists() or result_path.stat().st_size == 0:
                tail = '\n'.join(stdout.splitlines()[-15:])
                raise RuntimeError(
                        f"Korg produced no usable spectrum ({result_path}): julia exited "
                        f"0 but wrote no data. Last output:\n{tail}\n"
                        f"--- stderr ---\n{stderr}")

            self = cls.read(result_path)
            self.input = input
            self.llist = llist
            self.stdout = stdout
            self.stderr = stderr
            return self

    @staticmethod
    #=================================================================
    def _write_input_file(input, filename):
        # One "KEY value" per line. A list value repeats its key, which is how the
        # driver receives MU_VALUES and the per-element abundances without needing a
        # JSON package that the Korg environment does not carry.
        with open(filename, 'w') as f:
            for key, value in input.items():
                if isinstance(value, dict):
                    for name, item in value.items():
                        f.write(f'{key} {name} {item}\n')
                elif isinstance(value, (list, tuple, np.ndarray)):
                    for item in value:
                        f.write(f'{key} {item}\n')
                elif isinstance(value, bool):
                    f.write(f'{key} {str(value).lower()}\n')
                else:
                    f.write(f'{key} {value}\n')

    @staticmethod
    #=================================================================
    def _write_linelist_file(stack, llist, lam_min, lam_max):
        # Convert a tumult LineList (masked to the synthesis window) into the CSV the
        # driver reads, and return its path. The temp file lives as long as `stack`.
        llist_file = tmp_txtfile(stack)
        mask = (llist.lam0 >= lam_min) & (llist.lam0 <= lam_max)
        write_linelist(llist[mask], llist_file)
        return str(llist_file)

    @classmethod
    #=================================================================
    def concatenate(cls, runs):
        self = cls()
        self.runs = runs
        self.mus = runs[0].mus
        self.df = pd.concat([run.df for run in runs], axis=0)
        self.lam = self.df.iloc[:, 0].values
        return self

    @classmethod
    #=================================================================
    def from_tumult_run(cls, run, nthreads=0, **kwargs):
        nml = run.nml

        # Korg keys abundances by element symbol and reads them as A(X) when the driver
        # passes solar_relative=false, so the run's own dict goes through unchanged.
        # H is excluded because A(H) = 12 by definition and Korg rejects setting it.
        abnds = {k: v for k, v in run.abundances.items()
                 if k in atomic_weights and k != 'H'}

        vmic = nml['atmos_params']['vmic']
        vmic = vmic if vmic > 0 else np.mean(run.atmos.vmic)

        llist = getattr(run, 'llist', None)

        # Take the grid from the run's own wavelengths rather than from daa, which is
        # zero whenever spectrum_params set a resolution instead. Korg needs a linearly
        # spaced grid, so an R-based (constant-velocity) run is matched on sample count
        # and end points; across a narrow window the two grids differ by a fraction of
        # a step in the middle, and a comparison sensitive to that should interpolate.
        lam = np.asarray(run.lam)
        lmin, lmax, npoints = float(lam[0]), float(lam[-1]), lam.size
        lpad = nml['linelist_params']['dlam']

        input = {
                'ATMOSPHERE': nml['atmos_params']['atmos_file'],
                'RESULTFILE': 'tmp.spec',
                'LAMBDA_MIN': lmin,
                'LAMBDA_MAX': lmax,
                'NPOINTS': npoints,
                'VMIC': f'{vmic:.3f}',
                # HOP has already applied its own selection window when writing
                # linelist_out.txt, so keep every line the list carries
                'LINE_BUFFER': lpad,
        }
        input['ABUNDANCE'] = abnds

        # bind kwargs to cls.launch method
        launch = partial(cls.launch, **kwargs)

        with ExitStack() as stack:
            # Convert the run's line list once, up front, rather than re-converting it
            # inside every (possibly parallel) launch. Masked to the padded window; the
            # temp file lives until this ExitStack closes, i.e. past all launches.
            if llist is not None:
                input['LINELIST'] = cls._write_linelist_file(
                        stack, llist, lmin - lpad, lmax + lpad)

            if nthreads == 0:
                return launch(input)

            else:
                # Split on sample index so the chunks tile the run's grid exactly, with
                # no shared sample and none dropped between them. LINE_BUFFER then
                # applies to each chunk's own bounds rather than the whole window, so a
                # split spectrum differs from an unsplit one at the ~1e-3 level in
                # normalised flux: both are the same finite-buffer approximation, taken
                # around different edges. Synthesise in one piece when that matters.
                edges = np.concatenate([[0], np.cumsum(near_equal_split(npoints, nthreads))])

                args = []
                for i, (first, last) in enumerate(zip(edges[:-1], edges[1:] - 1)):
                    input_chunk = input.copy()
                    input_chunk['LAMBDA_MIN'] = float(lam[first])
                    input_chunk['LAMBDA_MAX'] = float(lam[last])
                    input_chunk['NPOINTS'] = int(last - first + 1)
                    input_chunk['RESULTFILE'] = f'tmp_{i:03d}.spec'
                    args.append(input_chunk)

                with mp.Pool(nthreads) as pool:
                    results = pool.starmap(launch, [(arg,) for arg in args])

                return cls.concatenate(results)

    #=================================================================
    #=================================================================
    def read_atmos(self):
        lines = self.stdout.split('\n')
        data = [l[len(ATMOS_KEY):] for l in lines if l.startswith(ATMOS_KEY)]
        string = '\n'.join(data)
        df = pd.read_csv(StringIO(string), sep=r'\s+', header=None,
                         names=['tau', 'temp', 'ne', 'nt', 'z'])

        atmos = AtmosClass()
        atmos.dim = '1D'
        for col in df.columns:
            setattr(atmos, col, df[col].values)

        atmos.ltau = np.log10(atmos.tau)
        atmos.pe = atmos.ne * BK * atmos.temp
        atmos.pg = atmos.nt * BK * atmos.temp
        return atmos

    #=================================================================
    def save(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    #=================================================================
    def load(cls, filename):
        with open(filename, 'rb') as f:
            self = pickle.load(f)
        return self


#=====================================================================
def run_julia(script, input_file, cwd=None, verbose=True):
    cwd = Path(cwd or KORG_PATH).expanduser()
    julia = Path(JULIA_BIN).expanduser()
    command = [str(julia), f'--project={cwd}', str(script), str(input_file)]

    pr = subprocess.run(command, cwd=str(cwd), capture_output=True, text=True)
    if verbose:
        print(pr.stdout)
        print(pr.stderr)

    if pr.returncode != 0:
        raise RuntimeError(
                f"julia exited with code {pr.returncode} (cwd={cwd}).\n"
                f"--- stdout ---\n{pr.stdout}\n--- stderr ---\n{pr.stderr}")

    return pr.stdout, pr.stderr
