#!/usr/bin/env python3
"""Candidate lines of the atom validation and the synthesis windows that cover them.

Sources, all restricted to the elements of the validated atoms (Fe, Mg, Ca, Si, Al, Mn):

- ``gbs``: the lines the Gaia benchmark analyses used, per star: Fe from Jofre et al.
  (2014, CDS J/A+A/564/A133, table6 usage flags), the others from Jofre et al. (2015,
  J/A+A/582/A81, the per-star final line tables). Their EWs (several methods per line)
  are kept in ``gbs_ews.csv`` as an independent check of our own measurements.
- ``diag``: the campaign's diagnostic Mg/Fe set (``../windows/diagnostic_mgfe.csv``).
- ``takeda``: the lines of Takeda & UeNo's solar centre-to-limb study (``tableE.dat``).
- ``al``: Al I lines in the benchmark range (the GBS analyses have no Al).

Every candidate is identified with its entry in the campaign's VALD list (the list the
syntheses read), which gives the wavelength, lower level energy and log gf used
throughout. Candidates are only candidates: which lines count as clean is decided per
star after the syntheses (``analyse.py``).

The windows (``windows_<el>.txt``) are each line +- ``half_width`` A, merged where they
overlap. ``half_width`` is 0.8 A, or 2.5 A for lines strong enough to have damping wings
in the cool giants (excitation-corrected strength log gf - 0.87 chi + A(X) - 7.5 > -3).

The NLTE syntheses give each line-list line the departure coefficients of the atom levels
nearest in energy to its lower and upper level (find_closest in HOP/linelist_mod.f90, no
tolerance). ``dE_lo_<atom>`` / ``dE_up_<atom>`` are those energy offsets for every
validated atom, and ``mismatch_<atom>`` flags a line with either above 0.1 eV (audit 2,
item 4): its NLTE result is reported separately.

    python select_lines.py REF_DIR VALD_LIST TAKEDA_TABLE ATOM_DIR OUT_DIR
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
ELEMENTS = ('Fe', 'Mg', 'Ca', 'Si', 'Al', 'Mn')
ION = {'I': 1, 'II': 2}
ROMAN = {1: 'I', 2: 'II'}
ATOM_FILES = {'mg86d': 'public/atom.mg86d'}       # relative to the atom folder, if not the name
Z = {'Fe': 26, 'Mg': 12, 'Ca': 20, 'Si': 14, 'Al': 13, 'Mn': 25}
A_SUN = {'Fe': 7.46, 'Mg': 7.55, 'Ca': 6.30, 'Si': 7.51, 'Al': 6.43, 'Mn': 5.42}
LAM_RANGE = (4775.0, 6825.0)          # covered by every benchmark spectrum
CLV_RANGE = (4200.0, 8000.0)          # the IAG centre-to-limb atlas
AL_LINES = (5557.06, 6696.02, 6698.67, 7835.31, 7836.13)
MATCH_AA = 0.04
MISMATCH_EV = 0.1
# validated atoms per element: the original (public, as its pandas conversion, which
# agrees with the classic file within 0.003 dex) and the repaired one; for Mg the public
# mg86d, the mg86e it was updated to, and the repaired mg86e
ATOMS = {'Fe': ('fe607a_pd', 'fe607a_fix'), 'Mg': ('mg86d', 'mg86e_pd', 'mg86e_fixed'),
         'Ca': ('ca105b_pd', 'ca105b_fix'), 'Si': ('si340_pd', 'si340_ce'),
         'Al': ('al_qmh_pd', 'al_qmh_ci'), 'Mn': ('mn281kbc_pd', 'mn281kbc_dedup')}
HW, HW_STRONG = 0.8, 2.5

# Jofre et al. (2014) table6 star columns, in order
GBS_FE_STARS = ['18Sco', '61CygA', '61CygB', 'alfCenA', 'alfCenB', 'alfCet', 'alfTau', 'Arcturus',
                'betAra', 'betGem', 'betHyi', 'betVir', 'delEri', 'epsEri', 'epsFor', 'epsVir', 'etaBoo',
                'gamSge', 'Gmb1830', 'HD107328', 'HD122563', 'HD140283', 'HD220009', 'HD22879', 'HD49933',
                'HD84937', 'ksiHya', 'muAra', 'muCas', 'muLeo', 'Procyon', 'psiPhe', 'Sun', 'tauCet']


def star_key(name):
    """Our star key for a CDS star name ('alf Cen A' -> 'alfCenA', 'mu CasA' -> 'muCas')."""
    key = name.replace(' ', '')
    return {'muCasA': 'muCas', 'Gmb1830': 'Gmb1830'}.get(key, key)


def read_vald(path):
    """The VALD list as a table: element, ion, lam (air), evi [eV], loggf."""
    rows, species = [], None
    with open(path) as f:
        for line in f:
            if line.startswith("'"):
                text = line.strip().strip("'")
                parts = text.split()
                if len(parts) == 2 and parts[1] in ION:
                    species = (parts[0], ION[parts[1]])
                elif text.startswith('0'):
                    continue
                else:
                    species = None
                continue
            if species and species[0] in ELEMENTS:
                v = line.split()
                rows.append((species[0], species[1], float(v[0]), float(v[1]), float(v[2])))
    return pd.DataFrame(rows, columns=['element', 'ion', 'lam', 'evi', 'loggf'])


def match_vald(vald, element, ion, lam, evi=None):
    """The VALD entry of a line quoted at ``lam`` (and ``evi``): the nearest wavelength
    within MATCH_AA, among entries within 0.06 eV of ``evi`` when it is given."""
    sel = vald[(vald.element == element) & (vald.ion == ion) & ((vald.lam - lam).abs() < MATCH_AA)]
    if evi is not None and len(sel):
        near = sel[(sel.evi - evi).abs() < 0.06]
        sel = near if len(near) else sel
    if not len(sel):
        return None
    return sel.loc[(sel.lam - lam).abs().idxmin()]


def gbs_fe(ref):
    """{(element, ion, lam): set of stars} from Jofre (2014), and the EW table."""
    ews = []
    for line in (ref / 'J_A_A_564_A133' / 'ew.dat').read_text().splitlines():
        vals = [float(x) for x in line[37:].split()]
        good = [v for v in vals if v > 0]
        ews.append(dict(star=star_key(line[:9].strip()), element='Fe',
                        ion=1 if line[10:14].strip() == 'FeI' else 2,
                        lam=float(line[17:24]), evi=float(line[29:34]),
                        ew_gbs=np.median(good) if good else np.nan,
                        ew_gbs_spread=np.std(good) if len(good) > 1 else np.nan,
                        n_methods=len(good)))
    ews = pd.DataFrame(ews)
    # table6 quotes 0.1 A; the EW table gives the line's 0.01 A wavelength and EP
    known = ews.groupby(['ion', 'lam']).evi.first().reset_index()
    usage = {}
    for line in (ref / 'J_A_A_564_A133' / 'table6.dat').read_text().splitlines():
        ion = 1 if line[:4].strip() == 'FeI' else 2
        lam = float(line[5:11])
        flags = line[14:].split()
        hit = known[(known.ion == ion) & ((known.lam - lam).abs() < 0.051)]
        if len(hit) != 1:
            print(f'Jofre (2014) Fe {ROMAN[ion]} {lam}: {len(hit)} EW-table entries, skipped')
            continue
        usage[('Fe', ion, float(hit.lam.iloc[0]), float(hit.evi.iloc[0]))] = \
            {GBS_FE_STARS[i] for i, f in enumerate(flags) if f == '1'}
    return usage, ews


def gbs_others(ref):
    """{(element, ion, lam): set of stars} from Jofre (2015), the line data and the EWs."""
    info, element = [], None
    for line in (ref / 'J_A_A_582_A81' / 'infoline.dat').read_text().splitlines():
        if line.startswith('# element:'):
            element = line.split()[-1]
        elif not line.startswith('#') and line.strip():
            v = line.split()
            info.append(dict(element=element, ion=int(v[1]), lam=float(v[0]),
                             loggf=float(v[2]), evi=float(v[3])))
    info = pd.DataFrame(info)
    stars = {}
    for line in (ref / 'J_A_A_582_A81' / 'stars.dat').read_text().splitlines():
        stars[line.split()[-1]] = star_key(line[24:36].strip())
    zel = {v: k for k, v in Z.items()}
    usage, ews = {}, []
    for fname, star in stars.items():
        for line in (ref / 'J_A_A_582_A81' / fname).read_text().splitlines():
            if line.startswith('#') or not line.strip():
                continue
            left, mid, _ = line.split('|')
            v = left.split()
            element = zel.get(int(v[0]))
            if element is None:
                continue
            lam = float(v[1])
            cand = info[(info.element == element) & ((info.lam - lam).abs() < 0.02)]
            ion = int(cand.ion.iloc[0]) if len(cand) else 1
            usage.setdefault((element, ion, lam), set()).add(star)
            good = [float(x) for x in mid.split() if float(x) > 0]
            ews.append(dict(star=star, element=element, ion=ion, lam=lam,
                            evi=float(cand.evi.iloc[0]) if len(cand) else np.nan,
                            ew_gbs=np.median(good) if good else np.nan,
                            ew_gbs_spread=np.std(good) if len(good) > 1 else np.nan,
                            n_methods=len(good)))
    return usage, info, pd.DataFrame(ews)


def takeda_lines(table):
    rows = []
    for line in Path(table).read_text().splitlines():
        code = float(line[21:26])
        z, ion = int(code), int(round(100 * (code - int(code)))) + 1
        element = {v: k for k, v in Z.items()}.get(z)
        if element is None:
            continue
        rows.append(dict(element=element, ion=ion, lam=float(line[27:35]),
                         evi=float(line[35:41]), code=line[:12].strip()))
    return pd.DataFrame(rows)


def level_offsets(lines, atom_dir):
    """Energy offsets of the departure levels each line takes, per atom (see above)."""
    from tumult.atom_audit import level_keys, load_atom
    for element, names in ATOMS.items():
        sel = lines.element == element
        sub = lines[sel]
        klo = 100 * sub.ion + sub.evi
        kup = klo + 12398.42 / (sub.lam * 1.000277)
        for name in names:
            keys = level_keys(load_atom(Path(atom_dir) / ATOM_FILES.get(name, name)))
            dlo = np.abs(keys[None, :] - klo.to_numpy()[:, None]).min(1)
            dup = np.abs(keys[None, :] - kup.to_numpy()[:, None]).min(1)
            lines.loc[sel, 'dE_lo_' + name] = dlo
            lines.loc[sel, 'dE_up_' + name] = dup
            lines.loc[sel, 'mismatch_' + name] = (dlo > MISMATCH_EV) | (dup > MISMATCH_EV)
    return lines


def windows(lines):
    """Merged (blue, red) windows over every line's +- half_width."""
    spans = sorted((r.lam - r.half_width, r.lam + r.half_width) for r in lines.itertuples())
    merged = [list(spans[0])]
    for lo, hi in spans[1:]:
        if lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return merged


def main(ref, vald_path, takeda_table, atom_dir, out):
    ref, out = Path(ref), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    vald = read_vald(vald_path)
    fe_usage, fe_ews = gbs_fe(ref)
    oth_usage, info, oth_ews = gbs_others(ref)
    ours = set(pd.read_csv(HERE / 'stars.csv', comment='#').star)
    cands = {}

    def add(element, ion, lam, evi, source, tag=''):
        hit = match_vald(vald, element, ion, lam, evi)
        if hit is None:
            print(f'no VALD match: {element} {ROMAN[ion]} {lam:.3f} ({source})')
            return None
        key = (element, ion, round(float(hit.lam), 3))
        row = cands.setdefault(key, dict(element=element, ion=ion, lam=float(hit.lam), evi=float(hit.evi),
                                         loggf=float(hit.loggf), gbs='', diag='', takeda=''))
        if tag:
            row[source] = ';'.join(sorted(set(filter(None, row[source].split(';') + tag.split(';')))))
        return key

    for (element, ion, lam, evi), stars in fe_usage.items():
        used = sorted(stars & ours)
        if used:
            add(element, ion, lam, evi, 'gbs', ';'.join(used))
    for (element, ion, lam), stars in oth_usage.items():
        used = sorted(stars & ours)
        if element in ELEMENTS and used:
            src = info[(info.element == element) & ((info.lam - lam).abs() < 0.02)]
            add(element, ion, lam, float(src.evi.iloc[0]) if len(src) else None, 'gbs', ';'.join(used))
    diag = pd.read_csv(HERE.parent / 'windows' / 'diagnostic_mgfe.csv')
    for r in diag.itertuples():
        if LAM_RANGE[0] < r.lam_air < CLV_RANGE[1]:
            add(r.element, int(r.ion), r.lam_air, r.evi, 'diag', r.surveys)
    for r in takeda_lines(takeda_table).itertuples():
        add(r.element, r.ion, r.lam, r.evi, 'takeda', r.code)
    for lam in AL_LINES:
        add('Al', 1, lam, None, 'diag', 'al')

    lines = pd.DataFrame(cands.values()).sort_values(['element', 'lam']).reset_index(drop=True)
    lines['eup'] = lines.evi + 12398.42 / (lines.lam * 1.000277)
    strength = lines.loggf - 0.87 * lines.evi + lines.element.map(A_SUN) - 7.5
    lines['half_width'] = np.where(strength > -3.0, HW_STRONG, HW)
    lines['benchmark'] = (lines.lam > LAM_RANGE[0]) & (lines.lam < LAM_RANGE[1])
    lines['clv'] = (lines.lam > CLV_RANGE[0]) & (lines.lam < CLV_RANGE[1])
    lines = level_offsets(lines, atom_dir)
    lines.to_csv(out / 'lines.csv', index=False, float_format='%.4f')
    pd.concat([fe_ews, oth_ews[oth_ews.element.isin(ELEMENTS)]]).to_csv(out / 'gbs_ews.csv', index=False,
                                                                         float_format='%.3f')
    for element in ELEMENTS:
        sel = lines[(lines.element == element) & (lines.benchmark | lines.clv)]
        spans = windows(sel)
        (out / f'windows_{element.lower()}.txt').write_text(
            '# blue_AA  red_AA (air)\n' + ''.join(f'{lo:.2f} {hi:.2f}\n' for lo, hi in spans))
        print(f'{element}: {len(sel)} lines ({sel.benchmark.sum()} in the benchmark range, '
              f'{(sel.gbs != "").sum()} gbs, {(sel.diag != "").sum()} diag, {(sel.takeda != "").sum()} takeda), '
              f'{len(spans)} windows, {sum(h - l for l, h in spans):.0f} A; '
              + ', '.join(f'{n} mismatched {int(sel["mismatch_" + n].sum())}' for n in ATOMS[element]))


if __name__ == '__main__':
    main(*sys.argv[1:6])
