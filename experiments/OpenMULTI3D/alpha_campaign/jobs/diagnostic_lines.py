"""Diagnostic line set for diagnostic-aware atom pruning (Mg and Fe).

Writes ``windows/diagnostic_mgfe.csv`` with one row per line: element, ion, lam_air [AA],
evi [eV above the stage ground], loggf, half_width [AA] of its window and source. The
pruning (tumult.atom_rates.diagnostic_relevance) takes the levels of these lines; the 1D
validation synthesises the windows. Other targets reuse the machinery with their own
csv in the same columns.

Sources:
  audit    the atom audit's windows (alpha_campaign/reports/2026-10-06_atom_audit.md):
           in each window every line of the element within 1.5 dex of its strongest
  sed      LINES_SED of delta_maps.py
  hrs      4MOST HRS arms (windows/hrs.txt): Fe I with strength s in (-1.3, -0.2) and no
           Fe I line within 0.6 AA and 1 dex of it, Fe II with s in (-3.0, -1.0), Mg I
           with s > -3
with s = loggf + log10(lam) - (5040/5777) evi, a solar weak-line strength.

    python diagnostic_lines.py <master_line.txt> [out.csv]
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
HRS_HALF_WIDTH = 0.4
# log of the singly ionised over the neutral stage of Mg and Fe in the solar line-forming
# layers, which puts strengths of the two stages on one scale
ION_OFFSET = 1.5


def strength(tab):
    return tab.loggf + np.log10(tab.lam0) - THETA * tab.evi


def hrs_arms():
    rows = [r.split() for r in (HERE / 'windows/hrs.txt').read_text().splitlines()
            if r.strip() and not r.startswith('#')]
    return [(float(r[0]), float(r[1])) for r in rows]


def in_windows(lam, windows):
    return np.any([(lam >= b) & (lam <= r) for b, r in windows], axis=0)


def rows(tab, element, ion, half_width, source):
    return pd.DataFrame(dict(element=element, ion=ion, lam_air=tab.lam0.round(3), evi=tab.evi.round(4),
                             loggf=tab.loggf, half_width=round(half_width, 3), source=source))


def isolated(tab, lam_sep=0.6, dex=1.0):
    """Lines without a line of the same table within lam_sep [AA] and dex of their strength."""
    lam, s = tab.lam0.to_numpy(), strength(tab).to_numpy()
    order = np.argsort(lam)
    lam, s = lam[order], s[order]
    keep = np.ones(len(lam), bool)
    for k in range(len(lam)):
        lo, hi = np.searchsorted(lam, [lam[k] - lam_sep, lam[k] + lam_sep])
        near = np.r_[lo:k, k + 1:hi]
        keep[k] = not np.any(s[near] > s[k] - dex)
    out = np.zeros(len(lam), bool)
    out[order] = keep
    return out


def build(master_line):
    llist = LineList.from_tumult_format(master_line).reset_index()
    lists = {el: {ion: llist[llist.ion == f'{el} {"I" * ion}'].reset_index(drop=True) for ion in (1, 2)}
             for el in ('Mg', 'Fe')}
    parts = []
    for el, wins in AUDIT.items():
        for b, r in wins:
            # every line of the element in the window within 1.5 dex of the strongest,
            # the singly ionised stage's strengths raised by ION_OFFSET
            win = pd.concat([tab[(tab.lam0 >= b) & (tab.lam0 <= r)].assign(stage=ion)
                             for ion, tab in lists[el].items()])
            s = strength(win) + ION_OFFSET * (win.stage - 1)
            for ion, sub in win[s >= s.max() - 1.5].groupby('stage'):
                parts.append(rows(sub, el, ion, (r - b) / 2, 'audit'))
    for el, lam, hw in LINES_SED:
        ion = 2 if el.endswith('II') else 1
        tab = lists[el[:2]][ion]
        near = tab[np.abs(tab.lam0 - lam) < 0.05]
        best = near.iloc[[np.argmax(strength(near))]]
        parts.append(rows(best, el[:2], ion, hw, 'sed'))
    arms = hrs_arms()
    fe1, fe2, mg1 = lists['Fe'][1], lists['Fe'][2], lists['Mg'][1]
    arm = fe1[in_windows(fe1.lam0, arms)]
    sel = arm[isolated(arm) & strength(arm).between(-1.3, -0.2)]
    parts.append(rows(sel, 'Fe', 1, HRS_HALF_WIDTH, 'hrs'))
    parts.append(rows(fe2[in_windows(fe2.lam0, arms) & strength(fe2).between(-3.0, -1.0)], 'Fe', 2,
                      HRS_HALF_WIDTH, 'hrs'))
    parts.append(rows(mg1[in_windows(mg1.lam0, arms) & (strength(mg1) > -3)], 'Mg', 1, HRS_HALF_WIDTH, 'hrs'))
    out = pd.concat(parts, ignore_index=True)
    # a line selected twice keeps its first source and widest window
    out['half_width'] = out.groupby(['element', 'ion', 'lam_air']).half_width.transform('max')
    return out.drop_duplicates(['element', 'ion', 'lam_air']).sort_values(['element', 'lam_air'],
                                                                          ignore_index=True)


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
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / 'windows/diagnostic_mgfe.csv'
    diag = build(Path(sys.argv[1]))
    diag.to_csv(out, index=False)
    print(diag.groupby(['element', 'ion', 'source']).size().to_string())
