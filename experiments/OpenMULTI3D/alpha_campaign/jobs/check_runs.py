"""Sanity checks of finished campaign runs, read through tumult.read().crop().

    python check_runs.py RUN_DIR [RUN_DIR ...]           # one line per run, exit 1 on a failure
    python check_runs.py --gate RUN_DIR                    # the 3D Fe gate on a 3D Mg run

A run fails when its log does not end in "OpenMULTI3D terminated", when a normalised
spectrum (any set component) holds a non-finite value or lies outside [0, 1.05], or, for
an NLTE run, when its departure coefficients are not finite and positive.  NLTE runs also
report whether the iterations met convlim and the largest |F_NLTE - F_LTE| of the
normalised spectra (both from the same run, so with the same atomic data).  The gate
additionally requires a run that met convlim or ended below a change of GATE_CHANGE, so
that a slow Mg run using up maxiter close to its convlim does not block Fe, and an
NLTE-LTE difference between 1e-4 and 0.5: zero would mean the departures were not used,
and more than half the continuum points to broken populations.  Logs are found as
logs/<set>_<model>_<run>.log (common.sh's run()).
"""
import re
import sys
from pathlib import Path

import numpy as np

import tumult

GATE_CHANGE = 1e-3     # the gate accepts a run that used up maxiter below this change
CAMPAIGN = Path('/mnt/beegfs/gemini/groups/bergemann/users/hoppe/alpha_campaign')


def log_of(run):
    rel = run.resolve().relative_to(CAMPAIGN / 'runs')
    return CAMPAIGN / 'logs' / ('_'.join(rel.parts) + '.log')


def check(run, gate=False):
    problems, notes = [], []
    log = log_of(run)
    text = log.read_text(errors='replace') if log.exists() else ''
    if 'OpenMULTI3D terminated' not in text:
        problems.append(f'no termination line in {log.name}')
    converged = 'terminated successfully' in text
    last = re.search(r'largest relative change was ([\d.]+) %', text)
    change = 0.0 if converged or not last else float(last[1]) / 100
    r = tumult.read(str(run))
    nlte = r.mode == 'NLTE'
    lo, hi = (r.lam.min(), r.lam.max())
    dmax, fmin, nsets = 0.0, 1.0, len(r.sets or [])
    for blue, red in [(lo, hi)]:            # the whole saved range, gaps included
        variants = [dict(LTE=lte) for lte in ([False, True] if nlte else [True])]
        variants += [dict(index=i) for i in range(len(r.sets or []))
                     if r.sets[i]['kind'] != 'continuum']
        for kw in variants:
            f = np.asarray(r.crop(xmin=blue, xmax=red, **kw).ydata.value)
            what = f"{blue:.0f}-{red:.0f} {kw}"
            if f.size == 0 or not np.isfinite(f).all():
                problems.append(f'empty or non-finite flux in {what}')
            elif f.min() < 0 or f.max() > 1.05:
                problems.append(f'flux range {f.min():.3f}-{f.max():.3f} in {what}')
            else:
                fmin = min(fmin, float(f.min()))
        if nlte:
            dn = np.asarray(r.crop(xmin=blue, xmax=red).ydata.value)
            dl = np.asarray(r.crop(xmin=blue, xmax=red, LTE=True).ydata.value)
            dmax = max(dmax, float(np.nanmax(np.abs(dn - dl))))
    notes.append(f'sets={nsets} minF={fmin:.3f}')
    if nlte:
        for atom in r.atoms:
            _, b = r.get_depart(iatom=atom.id)
            b = np.asarray(b)
            if not np.isfinite(b).all() or (b <= 0).any():
                problems.append(f'{atom.id} departures not finite and positive')
            else:
                notes.append(f'{atom.id} b {b.min():.2e}..{b.max():.2e}')
        notes.append(f'converged={converged} last change={change:.1e} max|NLTE-LTE|={dmax:.4f}')
        if gate and change >= GATE_CHANGE:
            problems.append(f'stopped at a change of {change:.1e}')
        if gate and not 1e-4 < dmax < 0.5:
            problems.append(f'NLTE-LTE difference {dmax:.2e} outside 1e-4..0.5')
    return problems, notes


def main(args):
    gate = '--gate' in args
    bad = 0
    for run in (Path(a) for a in args if a != '--gate'):
        try:
            problems, notes = check(run, gate)
        except Exception as exc:                      # unreadable output is a failure
            problems, notes = [f'{type(exc).__name__}: {exc}'], []
        bad += bool(problems)
        print(f"{'FAIL' if problems else 'ok  '} {run}: {'; '.join(problems + notes)}", flush=True)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
