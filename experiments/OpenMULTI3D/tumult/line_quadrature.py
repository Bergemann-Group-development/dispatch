"""Bound-bound wavelength grids for NLTE runs, generated per atom and atmosphere.

The frequency quadrature of an NLTE run is normally assembled from the ``nq``, ``qmax``
and ``q0`` columns of the model atom, which put the same points (in Doppler units) on
every line regardless of where and how it forms. This module builds the grid from the
lines' actual opacity profiles on the target atmosphere instead: HOP evaluates the
atom's transitions as an LTE line list on a coarse sampling of the model, a depth
envelope of the line-to-continuum opacity ratio is formed, and a fine grid is thinned
wherever the log of that envelope is linear between kept points. That is the criterion
matched to the trapezoidal weights ``set_weights`` gives each line on the shared grid.

:func:`generate_line_quadrature` writes two files a run then uses through
``&atom_params line_grid_file=... atom_file=...``: the grid (vacuum Angstrom) and a copy
of the atom in which each line's ``nq qmax q0`` are rewritten so that its
``lam_min``/``lam_max`` span the interval of the grid its rates are integrated over. The
copy differs from the original only in those three columns.

Only the classic (``atom.*``) model-atom format is handled, and only single-component
lines (``iw <= 1``).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .const import CC_KMS, HCE
from .physics import vac2air

CM_TO_EV = HCE * 1e-8            # energy in cm^-1 -> eV
QNORM_KMS = 10.0                 # the Fortran &atom_params qnorm default the atom copy assumes
RATIO_EPS = 1e-6                 # chi_line/chi_cont below this is float32 rounding noise
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X']


# --- the model atom -------------------------------------------------------------------
@dataclass
class ClassicAtom:
    """The parts of a classic model atom the grid generator needs.

    ``lines`` holds one row per bound-bound transition, in file order; ``records`` the
    line index of every text line of the file that is a bb record, so the file can be
    rewritten with those records replaced and everything else verbatim.
    """
    element: str
    ev_cm: np.ndarray            # level energies [cm^-1]
    g: np.ndarray
    ion: np.ndarray              # ionisation stage per level (1 = neutral)
    lines: np.ndarray            # structured: i, j (0-based), f, nq, qmax, q0, iw, ga, gv, gq
    text: list[str]
    records: dict[int, int]      # text line index -> bb line index

    @property
    def lam0_vac(self):
        """Ritz vacuum wavelengths [AA] of the lines."""
        return 1e8 / (self.ev_cm[self.lines['j']] - self.ev_cm[self.lines['i']])


def _is_entry(text):
    t = text.lstrip()
    return len(t.strip()) > 0 and not t.startswith('*')


def read_classic_atom(path):
    """Parse the header, levels and bound-bound records of a classic model atom.

    Follows atom_mod's read_classic: entries are the non-blank lines not starting with
    ``*``; hyperfine components (``iw > 1``) would follow their line record without
    comment skipping and are not supported here.
    """
    text = Path(path).read_text().splitlines()
    entries = [k for k, t in enumerate(text) if _is_entry(t)]
    pos = 0

    def entry():
        nonlocal pos
        k = entries[pos]
        pos += 1
        return k, text[k].split()

    _, (element, *_rest) = entry()
    _, _abund = entry()
    _, counts = entry()
    nk, nline = int(counts[0]), int(counts[1])

    ev_cm, g, ion = np.empty(nk), np.empty(nk), np.empty(nk, int)
    for n in range(nk):
        _, w = entry()
        ev_cm[n], g[n] = float(w[0]), float(w[1])
        ion[n] = int(w[-1])

    dtype = [('i', int), ('j', int), ('f', float), ('nq', int), ('qmax', float),
             ('q0', float), ('iw', int), ('ga', float), ('gv', float), ('gq', float)]
    lines = np.empty(nline, dtype=dtype)
    records = {}
    for n in range(nline):
        k, w = entry()
        jj, ii = int(w[0]), int(w[1])
        lines[n] = (min(ii, jj) - 1, max(ii, jj) - 1, float(w[2]), int(w[3]), float(w[4]),
                    float(w[5]), int(w[6]), float(w[7]), float(w[8]), float(w[9]))
        if lines[n]['iw'] > 1:
            raise NotImplementedError(f'line {n + 1} of {path} has {lines[n]["iw"]} '
                                      'hyperfine components; not supported')
        records[k] = n
    return ClassicAtom(element.strip(), ev_cm, g, ion, lines, text, records)


def write_atom_linelist(atom, path):
    """Write the atom's lines as an input line list, one block per ionisation stage.

    The record layout is the one linelist_out.txt uses (``lam_air evi loggf gv ga gq i
    j``): the lower-level energy is relative to the ground level of the line's ionisation
    stage, and the broadening columns are the atom's own, since atom_mod and the line-list
    reader feed them to the same prepare_* routines.
    """
    lam_air = vac2air(atom.lam0_vac)
    i, j = atom.lines['i'], atom.lines['j']
    loggf = np.log10(atom.g[i] * atom.lines['f'])
    stage = atom.ion[i]
    ground = {s: np.flatnonzero(atom.ion == s)[0] for s in np.unique(stage)}
    evi = np.array([(atom.ev_cm[ii] - atom.ev_cm[ground[s]]) * CM_TO_EV
                    for ii, s in zip(i, stage)])
    with Path(path).open('w') as f:
        for s in sorted(ground):
            sel = np.flatnonzero(stage == s)
            sel = sel[np.argsort(lam_air[sel])]
            f.write(f"'{atom.element} {ROMAN[s - 1]}'\n'0' {sel.size}\n")
            for n in sel:
                f.write(f"{lam_air[n]:12.4f}{evi[n]:10.4f}{loggf[n]:9.3f}"
                        f"{atom.lines['gv'][n]:14.6f}{atom.lines['ga'][n]:12.3E}"
                        f"{atom.lines['gq'][n]:12.3E}{i[n] + 1:6d}{j[n] + 1:6d}\n")
    return path


def write_atom_copy(atom, path, half_width_kms, qnorm):
    """Copy the atom with each line's ``nq qmax q0`` set to a 3-point linear window.

    With ``q0 = qmax`` the Fortran builds a linear grid, so ``qmax`` in units of ``qnorm``
    [km/s] places ``lam_min``/``lam_max`` exactly ``half_width_kms`` from the line centre;
    ``nq = 3`` keeps the atom's own grid to the centre and the two edges, which the run
    ignores anyway when it takes the grid from a file. Every other line of the file, and
    the records of lines whose half width is NaN, are written unchanged.
    """
    q = np.asarray(half_width_kms, float) / qnorm
    out = []
    for k, t in enumerate(atom.text):
        n = atom.records.get(k)
        if n is None or not np.isfinite(q[n]):
            out.append(t)
            continue
        w = t.split()
        w[3], w[4], w[5] = '3', f'{q[n]:.4f}', f'{q[n]:.4f}'
        out.append('  ' + '  '.join(w))
    Path(path).write_text('\n'.join(out) + '\n')
    return path


# --- the opacity envelope --------------------------------------------------------------
def merge_windows(lo, hi):
    """Sorted, merged ``(lo, hi)`` intervals from per-line bounds."""
    order = np.argsort(lo)
    merged = []
    for a, b in zip(lo[order], hi[order]):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def fine_grid(windows, lam0, R):
    """Constant-``R`` vacuum grid over each window, with every line centre inserted."""
    parts = []
    for lo, hi in windows:
        n = int(np.ceil(np.log(hi / lo) * R)) + 1
        parts.append(lo * np.exp(np.arange(n) / R))
        parts.append(lam0[(lam0 >= lo) & (lam0 <= hi)])
    return np.unique(np.concatenate(parts))


def line_envelope(hop, lam, linelist, dlam, windows, vmic, chunk=500_000):
    """``max`` over the loaded depths of ``chi_line/chi_cont`` on ``lam``.

    The depth maximum keeps the narrowest core and the widest wing the line shows
    anywhere in the structure, which is what one grid serving every depth has to
    resolve. A product (or sum of logs) over depth would instead be dominated by the
    steepest layer and compress the wings away.
    """
    env = np.empty(lam.size)
    for a in range(0, lam.size, chunk):
        b = min(a + chunk, lam.size)
        cont, _ = hop.opacity(lam[a:b])
        tot, _ = hop.opacity(lam[a:b], linelist=linelist, dlam=dlam, vmic=vmic,
                             windows=windows)
        ratio = np.asarray(tot, float) / np.asarray(cont, float) - 1.0
        env[a:b] = np.maximum(ratio, 0.0).max(axis=0)
    return env


# --- thinning ----------------------------------------------------------------------------
def thin_linear(x, g, tol):
    """Indices of the points of ``(x, g)`` kept so that linear interpolation between
    consecutive kept points reproduces ``g`` within ``tol`` (Douglas-Peucker).

    The end points are always kept. Piecewise linearity of the integrand between nodes
    is what makes trapezoidal weights exact, so this is the thinning matched to
    ``set_weights``.
    """
    keep = np.zeros(x.size, bool)
    keep[[0, -1]] = True
    stack = [(0, x.size - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        xs, gs = x[a + 1:b], g[a + 1:b]
        dev = np.abs(gs - (g[a] + (xs - x[a]) * (g[b] - g[a]) / (x[b] - x[a])))
        k = int(np.argmax(dev))
        if dev[k] <= tol:
            continue
        m = a + 1 + k
        keep[m] = True
        stack.append((a, m))
        stack.append((m, b))
    return np.flatnonzero(keep)


def line_extents(lam, env, floor, lam0, cap_lo, cap_hi, rise_dex=0.05):
    """Per-line ``(lo, hi)`` [AA]: from each centre outwards to where the envelope falls
    below the floor, rises again by ``rise_dex`` above its running minimum (a neighbour
    taking over; the line's own wing beyond that point is at or below the minimum), or
    reaches the generous window.

    Stopping at the minimum between blended lines is what keeps a line's rate integral
    to its own profile: in a line forest the total envelope never drops to the floor,
    and walking on would give every line the cores of all its neighbours.
    """
    rise = 10.0 ** rise_dex
    lo = np.empty(lam0.size)
    hi = np.empty(lam0.size)
    ic = np.searchsorted(lam, lam0)

    def walk(seg, fl):
        """Index within seg (0 = centre) of the last point of the line's own wing."""
        if seg.size == 0 or seg[0] <= fl[0]:
            return 0
        below = np.flatnonzero(seg <= fl)
        end = below[0] if below.size else seg.size
        seg = seg[:end]
        runmin = np.minimum.accumulate(seg)
        up = np.flatnonzero(seg > runmin * rise)
        if up.size:
            end = int(np.argmin(seg[:up[0]]))
            return end
        return end - 1

    for n, k in enumerate(ic):
        a = np.searchsorted(lam, cap_lo[n])
        b = np.searchsorted(lam, cap_hi[n], side='right')
        r = walk(env[k:b], floor[k:b])
        l = walk(env[k:a - 1 if a > 0 else None:-1], floor[k:a - 1 if a > 0 else None:-1])
        lo[n] = lam[k - l]
        hi[n] = lam[k + r]
    return lo, hi


# --- driver -----------------------------------------------------------------------------
def generate_line_quadrature(atom_file, atmos_file, out_dir, *, qnorm=QNORM_KMS, vmic=1.0,
                             feh=0.0, abundances='', relative_abundances=True,
                             half_window_kms=300.0, R=1e6, tol_dex=0.05,
                             rel_floor=1e-3, abs_floor=1e-2, lam_range=None,
                             verbose=True):
    """Write ``<out_dir>/line_grid.txt`` and ``<out_dir>/atom.<El>`` for a run.

    Steps: the atom's lines become a HOP line list; their opacity is evaluated on
    :data:`linestrength.NDEPTH` depths of ``atmos_file`` (a MARCS model) over
    ``+-half_window_kms`` around every line at resolving power ``R``; the depth
    envelope of ``chi_line/chi_cont`` is floored per line at ``min(rel_floor * peak,
    abs_floor)``, never below :data:`RATIO_EPS`, and thinned in log10 to ``tol_dex``;
    each line's window is where its envelope stays above that floor. ``qnorm`` must match the run's ``&atom_params``
    value, since the atom copy expresses the windows in those units. ``lam_range``
    (vacuum AA) restricts the lines considered, e.g. to a run's ``blue_cut``/``red_cut``.

    The relative floor sets how far down its own profile a line is resolved, which is
    what its (normalised) rate integral needs; the absolute one keeps the wings of
    lines that stay above ``abs_floor`` of the continuum farther out than that, for the
    spectrum they produce.
    """
    from .hop import HopEOS
    from .linestrength import read_marcs
    from .synthesis.runtime import HOP_LOCK

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    atom = read_classic_atom(atom_file)
    lam0_all = atom.lam0_vac
    use = np.ones(lam0_all.size, bool)
    if lam_range is not None:
        use &= (lam0_all >= lam_range[0]) & (lam0_all <= lam_range[1])
    lam0 = lam0_all[use]
    if verbose:
        print(f'{atom.element}: {lam0.size} of {lam0_all.size} lines, '
              f'{lam0.min():.1f}-{lam0.max():.1f} AA')

    fac = half_window_kms / CC_KMS
    cap_lo, cap_hi = lam0 * (1 - fac), lam0 * (1 + fac)
    windows = merge_windows(cap_lo, cap_hi)
    lam = fine_grid(windows, lam0, R)
    dlam = float(max(b - a for a, b in windows))
    if verbose:
        print(f'{len(windows)} windows, {lam.size} fine-grid points')

    linelist = out_dir / f'{atom.element}_lines.txt'
    write_atom_linelist(atom, linelist)
    T, rho = read_marcs(atmos_file)
    options = {}
    if abundances:
        options = {'abundances': str(abundances),
                   'relative_abundances': bool(relative_abundances)}
    with HOP_LOCK:
        hop = HopEOS(feh=float(feh), **options)
        try:
            hop.load_state(T, rho=rho, vmic=float(vmic))
            env = line_envelope(hop, lam, linelist, dlam, windows, vmic)
        finally:
            hop.close()

    # each point's floor comes from the peak of the nearest line centre
    order = np.argsort(lam0)
    centres = lam0[order]
    ic = np.searchsorted(lam, centres)               # centres are fine-grid points
    peak = env[ic]
    right = np.clip(np.searchsorted(centres, lam), 0, centres.size - 1)
    left = np.clip(right - 1, 0, centres.size - 1)
    near = np.where(np.abs(lam - centres[left]) < np.abs(lam - centres[right]), left, right)
    # HOP returns float32 opacities, so the ratio is quantised at ~1.2e-7 of the
    # continuum: below RATIO_EPS the envelope is rounding noise (0 or one quantum),
    # and a line that weak gets the minimal window rather than a grid fitted to noise
    floor = np.clip(rel_floor * peak[near], RATIO_EPS, abs_floor)

    lo, hi = line_extents(lam, env, floor, centres, cap_lo[order], cap_hi[order])
    g = np.log10(np.maximum(env, floor))
    keep = []
    for a, b in windows:
        sel = np.flatnonzero((lam >= a) & (lam <= b))
        keep.append(sel[thin_linear(lam[sel], g[sel], tol_dex)])
    grid = np.unique(np.concatenate([lam[np.concatenate(keep)], centres]))
    # snap every line's window outwards onto grid points (at least one beyond the centre
    # on each side), so that lam_min/lam_max land on them; the atom copy is symmetric,
    # so the wider side wins
    ig = np.searchsorted(grid, centres)
    lo = grid[np.minimum(np.searchsorted(grid, lo, side='right') - 1, ig - 1)]
    hi = grid[np.minimum(np.maximum(np.searchsorted(grid, hi), ig + 1), grid.size - 1)]
    half = np.maximum(centres - lo, hi - centres) / centres * CC_KMS     # [km/s]
    half_all = np.full(lam0_all.size, np.nan)
    half_all[np.flatnonzero(use)[order]] = half

    grid_file = out_dir / 'line_grid.txt'
    with grid_file.open('w') as f:
        f.write(f'# bound-bound grid for {Path(atom_file).name} on {Path(atmos_file).name}'
                f' [vacuum AA]\n')
        f.write(f'# R={R:g} tol_dex={tol_dex} rel_floor={rel_floor} abs_floor={abs_floor}'
                f' ratio_eps={RATIO_EPS:g}'
                f' half_window_kms={half_window_kms} vmic={vmic} qnorm={qnorm}\n')
        np.savetxt(f, grid, fmt='%.6f')
    atom_copy = out_dir / f'atom.{atom.element.lower()}'
    write_atom_copy(atom, atom_copy, half_all, qnorm)
    if verbose:
        lo, hi = centres * (1 - half / CC_KMS), centres * (1 + half / CC_KMS)
        per_line = np.searchsorted(grid, hi * (1 + 1e-9), side='right') - \
            np.searchsorted(grid, lo * (1 - 1e-9))
        print(f'{grid.size} grid points; per line mean {per_line.mean():.1f}, median '
              f'{np.median(per_line):.0f}, min {per_line.min()}, max {per_line.max()}')
        print(f'wrote {grid_file} and {atom_copy}')
    return grid_file, atom_copy


def line_windows(atom_copy, qnorm=QNORM_KMS):
    """``{(i, j): (lam0, lo, hi)}`` in vacuum AA from an atom copy's ``qmax`` columns, keyed
    by the 1-based lower/upper level indices linelist_out.txt records carry as ``i j``."""
    atom = read_classic_atom(atom_copy)
    lam0 = atom.lam0_vac
    fac = atom.lines['qmax'] * qnorm / CC_KMS
    return {(int(i) + 1, int(j) + 1): (float(l0), float(l0 * (1 - f)), float(l0 * (1 + f)))
            for i, j, l0, f in zip(atom.lines['i'], atom.lines['j'], lam0, fac)}


@dataclass(frozen=True)
class PreparedLineQuadrature:
    """The cached output of :func:`prepare_line_quadrature` for one atom on one model."""
    folder: str
    grid_file: str
    atom_file: str          # the copy carrying the line windows
    source_atom: str
    atmos_file: str
    settings: dict

    def to_dict(self):
        return dict(folder=self.folder, grid_file=self.grid_file, atom_file=self.atom_file,
                    source_atom=self.source_atom, atmos_file=self.atmos_file,
                    settings=dict(self.settings))

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: d[k] for k in ('folder', 'grid_file', 'atom_file', 'source_atom',
                                        'atmos_file', 'settings')})


def prepare_line_quadrature(output_root, atom_file, atmos_file, *, vmic=1.0, feh=0.0,
                            abundances='', relative_abundances=True, **params):
    """Generate the grid and atom copy once per (atom, model, settings) under ``output_root``.

    The digest covers the identities (path, size, mtime) of the atom and the model, and
    every generator setting, so a changed input never reuses a stale grid. A complete
    preparation is marked by its ``quadrature.json``; anything else is regenerated.
    """
    from .linestrength import _file_identity

    atom_id = _file_identity(atom_file)
    atmos_id = _file_identity(atmos_file)
    settings = dict(vmic=float(vmic), feh=float(feh), abundances=str(abundances),
                    relative_abundances=bool(relative_abundances), qnorm=QNORM_KMS,
                    ratio_eps=RATIO_EPS,
                    **{k: (float(v) if isinstance(v, (int, float)) else v)
                       for k, v in params.items()})
    identity = {'atom': atom_id, 'atmos': atmos_id, 'settings': settings}
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True, default=str)
                            .encode()).hexdigest()[:16]
    folder = Path(output_root).expanduser().resolve() / f'quadrature-{digest}'
    manifest = folder / 'quadrature.json'
    if manifest.is_file():
        prepared = PreparedLineQuadrature.from_dict(json.loads(manifest.read_text()))
        if Path(prepared.grid_file).is_file() and Path(prepared.atom_file).is_file():
            return prepared
    grid_file, atom_copy = generate_line_quadrature(
        atom_id['path'], atmos_id['path'], folder, vmic=vmic, feh=feh, abundances=abundances,
        relative_abundances=relative_abundances, verbose=False, **params)
    prepared = PreparedLineQuadrature(str(folder), str(grid_file), str(atom_copy),
                                      atom_id['path'], atmos_id['path'], settings)
    tmp = manifest.with_name(f'{manifest.name}.{os.getpid()}.partial')
    tmp.write_text(json.dumps(prepared.to_dict(), indent=2, sort_keys=True) + '\n')
    tmp.replace(manifest)
    return prepared


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('atom_file')
    p.add_argument('atmos_file', help='MARCS model')
    p.add_argument('out_dir')
    p.add_argument('--qnorm', type=float, default=QNORM_KMS, help="the run's &atom_params qnorm")
    p.add_argument('--vmic', type=float, default=1.0)
    p.add_argument('--feh', type=float, default=0.0)
    p.add_argument('--abundances', default='', help='HOP abundance overrides, e.g. Fe=+0.4')
    p.add_argument('--half-window-kms', type=float, default=300.0)
    p.add_argument('--R', type=float, default=1e6)
    p.add_argument('--tol-dex', type=float, default=0.05)
    p.add_argument('--rel-floor', type=float, default=1e-3)
    p.add_argument('--abs-floor', type=float, default=1e-2)
    p.add_argument('--lam-range', type=float, nargs=2, metavar=('BLUE', 'RED'),
                   help='vacuum AA')
    a = p.parse_args(argv)
    generate_line_quadrature(a.atom_file, a.atmos_file, a.out_dir, qnorm=a.qnorm,
                             vmic=a.vmic, feh=a.feh, abundances=a.abundances,
                             half_window_kms=a.half_window_kms, R=a.R, tol_dex=a.tol_dex,
                             rel_floor=a.rel_floor, abs_floor=a.abs_floor,
                             lam_range=a.lam_range)


if __name__ == '__main__':
    main()
