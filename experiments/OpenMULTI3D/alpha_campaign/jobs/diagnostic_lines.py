"""Diagnostic line set for diagnostic-aware atom pruning (Mg and Fe).

Writes ``windows/diagnostic_mgfe.csv`` with one row per line: element, ion, lam_air [AA],
evi [eV above the stage ground], loggf, half_width [AA] of its synthesis window, and
``surveys``, the ';'-joined sources that use or would use the line, so that the set can
be filtered. The pruning (tumult.atom_rates.diagnostic_relevance) takes the levels of
these lines; the 1D validation synthesises the windows. Other targets reuse the
machinery with their own csv in the same columns.

Sources:
  audit    the atom audit's windows (alpha_campaign/reports/2026-10-06_atom_audit.md):
           in each window every Mg/Fe line within 1.5 dex of the strongest
  sed      LINES_SED of delta_maps.py, the campaign's SED lines
  4most    4MOST HRS arms (windows/hrs.txt)
  galah    GALAH/HERMES arms (Buder et al. 2021)
  ges      Gaia-ESO UVES 580 nm setup and GIRAFFE HR21, from the Gaia-ESO line list
           (Heiter et al. 2021; input/linelists/nlte_ges_linelist_*)
  apogee   APOGEE H band, the three detector ranges, from the APOGEE DR17 atomic list
           (Smith et al. 2021) as Korg.jl distributes it
The surveys publish ranges and line lists, not one list of diagnostic lines, so each
survey's lines follow one rule within its ranges: Fe I with strength s in (-1.3, -0.2)
and no Fe I line within 0.6 AA (scaled with wavelength) and 1 dex of it, Fe II with s in
(-3.0, -1.0), Mg I with s > -3, where s = loggf + log10(lam) - (5040/5777) evi is a solar
weak-line strength. Entries of different lists closer than 0.03 AA and 0.03 eV are one
line.

    python diagnostic_lines.py <master_line.txt> <ges list> <apogee TS atoms list> [out.csv]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from tumult.linelist import LineList
from delta_maps import LINES_SED

HERE = Path(__file__).parent
THETA = 5040 / 5777

AUDIT = {
    'Mg': [(2848.0, 2856.0), (3827.0, 3841.0), (4166.3, 4168.3), (4570.4, 4571.8), (4702.0, 4704.0),
           (4729.3, 4730.8), (5164.0, 5189.0), (5527.4, 5529.4), (5710.1, 5712.1), (6318.0, 6319.6),
           (7386.9, 7388.5), (7690.6, 7692.6), (8735.0, 8737.0), (8804.0, 8810.0), (8922.6, 8924.6),
           (11825.0, 11831.5), (15038.5, 15042.0)],
    'Fe': [(l - 0.6, l + 0.6) for l in
           (4994.13, 5225.53, 5247.05, 5250.21, 5383.37, 5434.52, 5560.21, 5661.35, 6082.71, 6151.62,
            6173.34, 6200.31, 6240.65, 6252.56, 6265.13, 6336.82, 6703.57, 6750.15, 7090.38, 8688.62,
            5197.58, 5234.62, 5264.81, 5414.07, 6084.11, 6149.26, 6247.56, 6432.68, 6456.38, 6516.08)],
}
# survey: (line list, ranges [AA, air])
SURVEYS = {
    'galah': ('master', [(4713, 4903), (5648, 5873), (6478, 6737), (7585, 7887)]),
    'ges': ('ges', [(4770, 5770), (5830, 6830), (8484, 9001)]),
    'apogee': ('apogee', [(15140, 15810), (15850, 16435), (16470, 16960)]),
}
# window half width at 6000 AA; it scales with wavelength
HALF_WIDTH = 0.4
# log of the singly ionised over the neutral stage of Mg and Fe in the solar line-forming
# layers, which puts strengths of the two stages on one scale
ION_OFFSET = 1.5


def strength(tab):
    return tab.loggf + np.log10(tab.lam0) - THETA * tab.evi


def half_width(lam):
    return np.round(np.maximum(HALF_WIDTH, HALF_WIDTH * np.asarray(lam) / 6000), 2)


def hrs_arms():
    rows = [r.split() for r in (HERE / 'windows/hrs.txt').read_text().splitlines()
            if r.strip() and not r.startswith('#')]
    return [(float(r[0]), float(r[1])) for r in rows]


def in_ranges(lam, ranges):
    return np.any([(lam >= b) & (lam <= r) for b, r in ranges], axis=0)


def rows(tab, element, ion, width, source):
    return pd.DataFrame(dict(element=element, ion=ion, lam_air=tab.lam0.round(3), evi=tab.evi.round(4),
                             loggf=tab.loggf, half_width=width, surveys=source))


def isolated(tab, lam_sep=0.6, dex=1.0):
    """Lines without a line of the same table within lam_sep [AA, at 6000 AA] and dex of
    their strength."""
    lam, s = tab.lam0.to_numpy(), strength(tab).to_numpy()
    order = np.argsort(lam)
    lam, s = lam[order], s[order]
    keep = np.ones(len(lam), bool)
    for k in range(len(lam)):
        sep = lam_sep * max(1.0, lam[k] / 6000)
        lo, hi = np.searchsorted(lam, [lam[k] - sep, lam[k] + sep])
        near = np.r_[lo:k, k + 1:hi]
        keep[k] = not np.any(s[near] > s[k] - dex)
    out = np.zeros(len(lam), bool)
    out[order] = keep
    return out


def survey_lines(lists, ranges, source):
    """The Mg I, Fe I and Fe II lines a survey with these ranges would use."""
    parts = []
    fe1 = lists['Fe'][1][in_ranges(lists['Fe'][1].lam0, ranges)]
    fe1 = fe1[isolated(fe1) & strength(fe1).between(-1.3, -0.2)]
    fe2 = lists['Fe'][2][in_ranges(lists['Fe'][2].lam0, ranges)]
    fe2 = fe2[strength(fe2).between(-3.0, -1.0)]
    mg1 = lists['Mg'][1][in_ranges(lists['Mg'][1].lam0, ranges)]
    mg1 = mg1[strength(mg1) > -3]
    for tab, el, ion in ((fe1, 'Fe', 1), (fe2, 'Fe', 2), (mg1, 'Mg', 1)):
        parts.append(rows(tab, el, ion, half_width(tab.lam0), source))
    return parts


def by_stage(llist):
    """{element: {stage: table}} for Mg and Fe, whatever the case of the species label."""
    ion = llist.ion.str.upper().str.split()
    out = {}
    for el in ('Mg', 'Fe'):
        out[el] = {st: llist[(ion.str[0] == el.upper()) & (ion.str[1] == 'I' * st)].reset_index(drop=True)
                   for st in (1, 2)}
    return out


def merge(parts):
    """One row per line: entries of the same element and stage within 0.03 AA and 0.03 eV
    are one line (lists differ by up to 0.02 eV in the energy of one level), which keeps the first entry's data, the widest window and all sources."""
    tab = pd.concat(parts, ignore_index=True).sort_values(['element', 'ion', 'lam_air'], kind='stable')
    out = []
    for _, grp in tab.groupby(['element', 'ion'], sort=False):
        lam, ev = grp.lam_air.to_numpy(), grp.evi.to_numpy()
        new = np.r_[True, (np.diff(lam) > 0.03)]
        line = np.cumsum(new)
        # a second level within a wavelength group splits it
        for k in range(1, len(lam)):
            if not new[k] and abs(ev[k] - ev[k - 1]) > 0.03:
                line[k:] += 1
        for _, one in grp.groupby(line):
            first = one.iloc[0].copy()
            first['half_width'] = one.half_width.max()
            first['surveys'] = ';'.join(dict.fromkeys(';'.join(one.surveys).split(';')))
            out.append(first)
    return pd.DataFrame(out).sort_values(['element', 'lam_air'], ignore_index=True)


def build(master_line, ges_line, apogee_line):
    lists = {'master': by_stage(LineList.from_tumult_format(master_line).reset_index()),
             'ges': by_stage(LineList.from_tumult_format(ges_line).reset_index()),
             'apogee': by_stage(LineList.from_ts_format(apogee_line).reset_index())}
    master = lists['master']
    parts = []
    for el, wins in AUDIT.items():
        for b, r in wins:
            # every line of the element in the window within 1.5 dex of the strongest,
            # the singly ionised stage's strengths raised by ION_OFFSET
            win = pd.concat([tab[(tab.lam0 >= b) & (tab.lam0 <= r)].assign(stage=ion)
                             for ion, tab in master[el].items()])
            if not len(win):
                continue
            s = strength(win) + ION_OFFSET * (win.stage - 1)
            for ion, sub in win[s >= s.max() - 1.5].groupby('stage'):
                parts.append(rows(sub, el, ion, round((r - b) / 2, 3), 'audit'))
    for el, lam, hw in LINES_SED:
        ion = 2 if el.endswith('II') else 1
        tab = master[el[:2]][ion]
        near = tab[np.abs(tab.lam0 - lam) < 0.05]
        parts.append(rows(near.iloc[[np.argmax(strength(near))]], el[:2], ion, hw, 'sed'))
    parts += survey_lines(master, hrs_arms(), '4most')
    for survey, (source, ranges) in SURVEYS.items():
        parts += survey_lines(lists[source], ranges, survey)
    return merge(parts)


def windows(diag, element):
    """Merged synthesis windows [AA] of an element's diagnostic lines."""
    sub = diag[diag.element == element].sort_values('lam_air')
    out = []
    for lam, hw in zip(sub.lam_air, sub.half_width):
        b, r = lam - hw, lam + hw
        if out and b <= out[-1][1]:
            out[-1][1] = max(out[-1][1], r)
        else:
            out.append([b, r])
    return out


if __name__ == '__main__':
    out = Path(sys.argv[4]) if len(sys.argv) > 4 else HERE / 'windows/diagnostic_mgfe.csv'
    diag = build(*(Path(a) for a in sys.argv[1:4]))
    diag.to_csv(out, index=False)
    for survey in ('audit', 'sed', '4most', 'galah', 'ges', 'apogee'):
        sub = diag[diag.surveys.str.split(';').apply(lambda s: survey in s)]
        print(survey, sub.groupby(['element', 'ion']).size().to_dict())
    print('union', diag.groupby(['element', 'ion']).size().to_dict())
