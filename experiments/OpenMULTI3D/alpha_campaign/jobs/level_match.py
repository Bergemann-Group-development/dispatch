"""How the line-list lines of a species map onto a model atom's levels in a spectrum on
precomputed departures (stage 2 of the NLTE route).

    python level_match.py LINELIST ATOM SPECIES [SPECIES ...] [--windows FILE]

HOP gives each line of the atom's element the departure coefficients of the levels
nearest in energy (HOP/linelist_mod.f90, find_closest over the departure set's energies,
tagged 100 x stage): the lower level at E_low, the upper at E_low + h nu.  There is no
tolerance, so a line whose levels the atom lacks silently takes its neighbours'.  This
script reports, per species, the share of lines (all lines, and the lines with
log gf - E_low[eV] > -3 as a crude strength cut) whose lower and upper energies both lie
within 0.01, 0.05 and 0.2 eV of an atom level, and prints the matched levels of the
species' strongest lines.  Energies are referred to each stage's ground level as in
atom_t%ev_keys.
"""
import sys
from pathlib import Path

import numpy as np

from tumult.atom import AtomClass

HC = 12398.419843320026        # eV A


def windows(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith('#')]
    return [(float(r[0]), float(r[1])) for r in rows]


def block(lines, species, win):
    i = lines.index(f"'{species}'")
    n = int(lines[i + 1].split()[1])
    d = np.array([[float(x) for x in l.split()[:3]] for l in lines[i + 2:i + 2 + n]])
    keep = np.zeros(len(d), bool)
    for lo, hi in win:
        keep |= (d[:, 0] >= lo) & (d[:, 0] <= hi)
    return d[keep]


def main(args):
    win = [(0, 1e9)]
    if '--windows' in args:
        k = args.index('--windows')
        win = windows(args[k + 1])
        args = args[:k] + args[k + 2:]
    linelist, atom_file, species = args[0], args[1], args[2:]
    a = AtomClass(atom_file)
    ev = np.array([l.ev for l in a.level])
    ion = np.array([l.ion for l in a.level])
    label = [str(l.label).strip("' ") for l in a.level]
    keys = ev.copy()
    ground = ev[0]
    for k in range(len(ev)):
        if k and ion[k] > ion[k - 1]:
            ground = ev[k]
        keys[k] = 100 * ion[k] + ev[k] - ground
    lines = Path(linelist).read_text().splitlines()
    print(f'{Path(atom_file).name}: {len(ev)} levels')
    for sp in species:
        stage = sp.count('I')
        d = block(lines, sp, win)
        if not len(d):
            print(f'  {sp}: no lines in the windows')
            continue
        lam, elo, lgf = d.T
        klo = 100 * stage + elo
        kup = klo + HC / lam
        il = np.abs(keys[None, :] - klo[:, None]).argmin(1)
        iu = np.abs(keys[None, :] - kup[:, None]).argmin(1)
        dl, du = np.abs(keys[il] - klo), np.abs(keys[iu] - kup)
        strong = lgf - elo > -3
        nlev = int((ion == stage).sum())
        for name, sel in (('all', np.ones(len(d), bool)), ('strong', strong)):
            n = int(sel.sum())
            frac = [((dl[sel] < t) & (du[sel] < t)).mean() if n else np.nan for t in (0.01, 0.05, 0.2)]
            print(f'  {sp:6s} {name:6s} n={n:5d} (atom levels of this stage: {nlev:3d})  '
                  f'both within 0.01/0.05/0.2 eV: {frac[0]:.2f} / {frac[1]:.2f} / {frac[2]:.2f}')
        for j in np.argsort(lgf - elo)[::-1][:4]:
            print(f'    {lam[j]:9.3f} loggf {lgf[j]:+.2f} Elow {elo[j]:.3f} -> {label[il[j]]} '
                  f'(dE {dl[j]:.3f}) .. {label[iu[j]]} (dE {du[j]:.3f})')


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1:])
