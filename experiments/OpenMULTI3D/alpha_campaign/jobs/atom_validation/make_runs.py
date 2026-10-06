#!/usr/bin/env python3
"""Namelists of the atom validation: 1D NLTE syntheses of the benchmark stars and the
solar centre-to-limb variation, for the original and the repaired atom of each element.

One run per (star, element, atom, atmosphere) solves the statistical equilibrium at the
star's reference abundance A_ref and synthesises the element's windows
(``select_lines.py``) with the full FGK line list. Its output sets (``sets.txt``) are a
curve-of-growth ladder on those departure coefficients: NLTE and LTE at
A_ref + ``LADDER``, the unperturbed continuum, and the element removed (``da`` =
``BLEND_FREE``, the blends alone). A star whose result lands far from A_ref gets a second
pass solved at the result (``--a-ref FILE``).

A_ref is the solar abundance of the synthesis (abund_magg22, the code's default) plus
the benchmark [X/H]: [Fe/H] for Fe and Al, the GBS v2.1 [X/H] (Jofre et al. 2015, CDS
III/281) for Mg, Si, Ca and Mn. The background composition is scaled by [Fe/H] with the
star's [alpha/Fe]; vmic is the benchmark value in both atmospheres.

The solar CLV runs (batch ``clv``) give intensities at the IAG atlas's mu values and two
more limb points for Takeda's data; their flux is not used (the custom quadrature is not
a flux quadrature).

    python make_runs.py REF_DIR LINES_DIR ATMOS_DIR OUT_DIR [--a-ref FILE --batch NAME]

OUT_DIR mirrors $C/runs/atom_validation: nml/<batch>/*.nml, sets/, windows/, atmos/,
runs_<batch>.csv and tasks_<batch>.txt (one namelist per line, longest first).
"""
import shutil
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
import make_nml  # noqa: E402
from select_lines import ATOMS, ATOM_FILES  # noqa: E402

make_nml.DATADIR = 'data'          # outputs under <run dir>/data/nml/<batch>/<name>/
SOLAR = {'Fe': 7.51, 'Mg': 7.58, 'Ca': 6.32, 'Si': 7.567, 'Al': 6.48, 'Mn': 5.53}   # abund_magg22
GBS_ORDER = ['Mg', 'Si', 'Ca', 'Ti', 'Sc', 'V', 'Cr', 'Mn', 'Co', 'Ni']
LADDER = (-0.6, -0.4, -0.2, 0.2, 0.4, 0.6)
BLEND_FREE = -5.0
CLV_MU = (1.0, 0.99, 0.98, 0.97, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1)
ATMOS = ('av3d', 'marcs')
R_SAMPLING = 2e5
TIO_BELOW = 4200                   # Teff [K] below which the TiO bands are added
# relative cost per run (fe607a ~300 s, the others ~90 s at 2 threads in the audit)
COST = {'Fe': 4, 'Mg': 1, 'Ca': 1, 'Si': 1, 'Al': 1, 'Mn': 1}


def reference_abundances(ref, stars):
    """A_ref per star and element, with the GBS value and its uncertainty where it exists."""
    rows = {}
    for line in (Path(ref) / 'III_281' / 'gbs.dat').read_text().splitlines():
        v = line[91:].split()
        xh = v[26::2]
        exh = v[27::2]
        rows[line[:10].replace(' ', '')] = {el: (float(xh[i]), float(exh[i])) for i, el in enumerate(GBS_ORDER)}
    out = []
    for s in stars.itertuples():
        for el in ATOMS:
            xh, err = rows.get(s.star, {}).get(el, (None, None))
            if el in ('Fe', 'Al') or xh is None:
                xh, err = s.feh, None
            if s.star == 'Sun':
                xh = 0.0
            out.append(dict(star=s.star, element=el, xh_ref=xh, e_xh=err, a_ref=round(SOLAR[el] + xh, 3)))
    return pd.DataFrame(out)


def set_file(path):
    """The ladder manifest every run reads (one element per file)."""
    rows = ['# continuum, the NLTE run itself, LTE, the ladder, then the element removed',
            '1 continuum - -', '2 nlte - -', '3 lte - -']
    n = 3
    for kind in ('nlte', 'lte'):
        for da in LADDER:
            n += 1
            rows.append(f'{n} {kind} {{el}} {da:+.2f}')
    rows.append(f'{n + 1} lte {{el}} {BLEND_FREE:+.1f}')
    for el in ATOMS:
        (path / f'ladder_{el.lower()}.txt').write_text('\n'.join(rows).format(el=el) + '\n')


def linelists(teff):
    files = list(make_nml.FGK_LINELISTS)
    if teff < TIO_BELOW:
        files += [f'{m}TiO.txt' for m in (46, 47, 48, 49, 50)]
    return [f'{make_nml.MASTER}/{f}' for f in files]


def one(name, star, element, atom, atmos, a_ref, rt_output=None):
    atom_file = ATOM_FILES.get(atom, atom)
    atom_path = f'atoms/{atom_file}' if not atom_file.startswith('public/') else atom_file
    return make_nml.namelist(
        dict(atmos_format='text', atmos_file=f'atmos/{name}_{atmos}.atmos', use_ne=False, amr='T',
             vmic=float(star.vmic)),
        window=f'windows/windows_{element.lower()}.txt', R=R_SAMPLING, linelist=linelists(star.teff),
        atom=dict(file=atom_path, abundance=float(a_ref), nlte=dict(convlim=1e-5), maxiter=400),
        feh=float(star.feh), output=dict(set_file=f'sets/ladder_{element.lower()}.txt'),
        rt_output=rt_output, composition=dict(alpha_enhancement=float(star.alpha)),
        absmet=make_nml.absmet_for(star.feh))


def main(ref, lines_dir, atmos_dir, out, a_ref_file=None, batch='b'):
    out = Path(out)
    for sub in ('sets', 'windows', 'atmos'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    set_file(out / 'sets')
    for f in Path(lines_dir).glob('windows_*.txt'):
        shutil.copy(f, out / 'windows' / f.name)
    for f in Path(atmos_dir).glob('*.atmos'):
        shutil.copy(f, out / 'atmos' / f.name)
    shutil.copy(Path(atmos_dir) / 'atmospheres.log', out / 'atmos')
    stars = pd.read_csv(HERE / 'stars.csv', comment='#')
    refs = reference_abundances(ref, stars).set_index(['star', 'element'])
    if a_ref_file:          # second pass: A_ref from a table of (star, element, atom, atmos, a_ref)
        override = pd.read_csv(a_ref_file).set_index(['star', 'element', 'atom', 'atmos']).a_ref
    runs = []
    for s in stars.itertuples():
        for el, atoms in ATOMS.items():
            for atom in atoms:
                for atmos in ATMOS:
                    a = refs.loc[(s.star, el), 'a_ref']
                    if a_ref_file:
                        key = (s.star, el, atom, atmos)
                        if key not in override.index:
                            continue
                        a = override[key]
                    runs.append(dict(batch=batch, star=s.star, element=el, atom=atom, atmos=atmos, a_ref=a,
                                     kind='flux', name=f'{s.star}_{atom}_{atmos}'))
                    if s.star == 'Sun' and not a_ref_file:
                        runs.append(dict(batch='clv', star=s.star, element=el, atom=atom, atmos=atmos, a_ref=a,
                                         kind='clv', name=f'Sun_{atom}_{atmos}'))
    runs = pd.DataFrame(runs)
    star_by = stars.set_index('star')
    clv = dict(angle_quad='custom', custom_mu=' '.join(f'{m:g}' for m in CLV_MU), nphi=1)
    for r in runs.itertuples():
        path = out / 'nml' / r.batch / f'{r.name}.nml'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(one(r.star, star_by.loc[r.star], r.element, r.atom, r.atmos, r.a_ref,
                            rt_output=clv if r.kind == 'clv' else None))
    runs['nml'] = 'nml/' + runs.batch + '/' + runs.name + '.nml'
    runs['cost'] = runs.element.map(COST)
    for b, sel in runs.groupby('batch'):
        sel.to_csv(out / f'runs_{b}.csv', index=False)
        (out / f'tasks_{b}.txt').write_text('\n'.join(sel.sort_values('cost', ascending=False,
                                                                       kind='stable').nml) + '\n')
        print(f'{b}: {len(sel)} runs')
    refs.reset_index().to_csv(out / 'reference_abundances.csv', index=False)


if __name__ == '__main__':
    args = sys.argv[1:]
    opts = {}
    for key in ('--a-ref', '--batch'):
        if key in args:
            k = args.index(key)
            opts[key] = args[k + 1]
            del args[k:k + 2]
    main(*args, a_ref_file=opts.get('--a-ref'), batch=opts.get('--batch', 'b'))
