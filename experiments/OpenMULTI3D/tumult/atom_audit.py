"""Sanity checks for model atoms.

The checks work on a :class:`tumult.fast_atom.PandasAtom` (``load_atom`` converts the
classic MULTI format on the fly) and reproduce what dispatch.x does with the data:
level energies define line and edge frequencies only through ``atom_mod%read_classic``,
while ``read_pandas`` takes each line's frequency from its tabulated air wavelength; bound-
free cross sections are interpolated with the monotone cubic Hermite scheme of
``trans_mod%hermite_interp`` and set to zero outside the tabulated range; collision
coefficients are interpolated in temperature with the same scheme (on the logarithm for
the keys listed in ``LOG_KEYS``) and held constant outside the tabulated range.

Each ``check_*`` function returns a table with one row per level, transition or
continuum, and a list of :class:`Finding`. ``audit`` runs all of them; the command line

    python -m tumult.atom_audit input/atoms/atom.mg86g --nist DIR --linelist FILE

prints the findings. Reference data (NIST ASD levels, line lists, original cross
sections) are read by :mod:`tumult.atomic_data`.
"""
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import exp1

from .const import BK, CC, EE, EM, HH, PI, RYDV
from .physics import air2vac

CM_TO_EV = HH * CC / EE          # 1 cm^-1 in eV
HC_AA_EV = HH * CC / EE * 1e8    # E [eV] = HC_AA_EV / lam [AA]
H_MINUS_EV = 0.754               # electron affinity of hydrogen [eV]
A0 = 0.529177e-8                 # Bohr radius [cm]
KRAMERS = 2.815e29               # hydrogenic photoionisation constant [cm^2 Hz^3]

# Keys whose tabulated values dispatch interpolates in log (read_classic); read_pandas
# takes the log of every key.
LOG_KEYS = ('OHM', 'OHMEGA', 'CE', 'CI', 'CH0')
SEVERITIES = ('serious', 'moderate', 'minor', 'info')


@dataclass
class Finding:
    """One result of a check: what was tested, how bad, and the evidence."""
    check: str
    severity: str
    message: str
    count: int = 0
    items: list = field(default_factory=list)

    def __str__(self):
        text = f'[{self.severity}] {self.check}: {self.message}'
        if self.items:
            text += '\n      e.g. ' + '; '.join(str(i) for i in self.items[:5])
        return text


# ---------------------------------------------------------------------------------------
# dispatch.x numerics
# ---------------------------------------------------------------------------------------
def hermite_slopes(xx, yy):
    """Node slopes of ``interpol_mod%pchip_slopes8``: weighted harmonic mean of the
    neighbouring secants (zero across a sign change or a flat segment), and quadratic
    end slopes."""
    xx = np.asarray(xx, float)
    yy = np.asarray(yy, float)
    n = len(xx)
    dx = np.diff(xx)
    dy = np.diff(yy)
    slopes = np.zeros(n)
    if n == 2:
        slopes[:] = dy[0] / dx[0]
        return slopes
    d1 = dy[:-1] / dx[:-1]
    d2 = dy[1:] / dx[1:]
    alpha = (1 + dx[:-1] / (dx[1:] + dx[:-1])) / 3
    same = (dy[:-1] != 0) & (dy[1:] != 0) & (np.sign(d1) == np.sign(d2))
    with np.errstate(divide='ignore', invalid='ignore'):
        harm = 1 / ((1 - alpha) / d1 + alpha / d2)
    slopes[1:-1] = np.where(same, harm, 0)
    slopes[0] = 2 * dy[0] / dx[0] - slopes[1]
    slopes[-1] = 2 * dy[-1] / dx[-1] - slopes[-2]
    return slopes


def hermite_eval(xx, yy, slopes, xnew):
    """``interpol_mod%hermite_interpolation``: cubic Hermite on the nodes, held at the
    end values outside ``[xx[0], xx[-1]]``."""
    xx = np.asarray(xx, float)
    yy = np.asarray(yy, float)
    xnew = np.asarray(xnew, float)
    i = np.clip(np.searchsorted(xx, xnew, side='right') - 1, 0, len(xx) - 2)
    h = xx[i + 1] - xx[i]
    t = (xnew - xx[i]) / h
    t2 = t * t
    t3 = t2 * t
    p = 3 * t2 - 2 * t3
    val = (1 - p) * yy[i] + p * yy[i + 1] + h * ((t3 - 2 * t2 + t) * slopes[i] + (t3 - t2) * slopes[i + 1])
    val = np.where(xnew <= xx[0], yy[0], val)
    return np.where(xnew >= xx[-1], yy[-1], val)


def bf_interpolate(lam_table, alpha_table, lam):
    """Cross section that dispatch.x puts on the wavelengths ``lam`` [AA] for a tabulated
    continuum (``trans_mod%hermite_interp``): Hermite in frequency, zero outside."""
    lam_table = np.asarray(lam_table, float)
    order = np.argsort(-lam_table)                  # increasing frequency
    nu_tab = CC * 1e8 / lam_table[order]
    alpha_tab = np.asarray(alpha_table, float)[order]
    nu = CC * 1e8 / np.asarray(lam, float)
    if len(nu_tab) == 1:
        # a one-point table is nonzero only at its own frequency
        return np.where(nu == nu_tab[0], alpha_tab[0], 0.0)
    alpha = hermite_eval(nu_tab, alpha_tab, hermite_slopes(nu_tab, alpha_tab), nu)
    return np.where((nu < nu_tab[0]) | (nu > nu_tab[-1]), 0.0, alpha)


def collision_interpolate(temps, values, temp, key):
    """Tabulated collision value at ``temp`` as gencol forms it (log interpolation for
    LOG_KEYS, constant outside the table)."""
    temps = np.asarray(temps, float)
    values = np.asarray(values, float)
    if key in LOG_KEYS:
        yy = np.log(values)
        return np.exp(hermite_eval(temps, yy, hermite_slopes(temps, yy), temp))
    return hermite_eval(temps, values, hermite_slopes(temps, values), temp)


def planck_photon_weight(nu, temp):
    """Photon-number weight B_nu/(h nu) up to a constant, the integrand weight of a
    photoionisation rate in a Planck field (and of the recombination rate in LTE)."""
    xx = HH * nu / (BK * temp)
    with np.errstate(over='ignore'):
        return nu**2 / np.expm1(xx)


def kramers(lam, binding_ev, charge):
    """Hydrogenic cross section [cm^2] with Gaunt factor 1 for a level bound by
    ``binding_ev`` in an ion of core charge ``charge`` (1 for a neutral atom)."""
    nu = CC * 1e8 / np.asarray(lam, float)
    neff = charge * np.sqrt(RYDV / binding_ev)
    return KRAMERS * charge**4 / (neff**5 * nu**3)


# ---------------------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------------------
def is_pandas_format(path):
    """True for an atom written by PandasAtom.write (a header row naming 'ev')."""
    with open(path) as f:
        rows = [r for r in f if r.strip() and not r.lstrip().startswith('*')][:4]
    return len(rows) == 4 and 'ev' in rows[3].split()


def load_atom(path):
    """A PandasAtom from either atom format."""
    from .fast_atom import PandasAtom
    path = Path(path)
    if is_pandas_format(path):
        return PandasAtom.read(path)
    return PandasAtom.from_old(path)


def ion_ground(atom):
    """Index of the first level of each ion stage, keyed by stage."""
    ions = atom.lvl.ion.to_numpy()
    return {int(s): int(np.argmax(ions == s)) for s in np.unique(ions)}


def continuum_level(atom):
    """For each level, the level that dispatch ionises it to (first level of the next
    stage), or -1."""
    grounds = ion_ground(atom)
    return np.array([grounds.get(int(s) + 1, -1) for s in atom.lvl.ion])


# ---------------------------------------------------------------------------------------
# levels
# ---------------------------------------------------------------------------------------
def check_levels(atom, energy_tolerance=0.01):
    """Level structure: ordering, duplicates, and connections of every level.

    ``energy_tolerance`` [cm^-1] decides when two levels of an ion count as duplicates.
    """
    lvl = atom.lvl.reset_index(drop=True).copy()
    nk = len(lvl)
    lvl['ecm'] = lvl.ev / CM_TO_EV
    bbt = getattr(atom, 'bbt', None)
    bft = getattr(atom, 'bft', None)
    col = getattr(atom, 'col', None)
    lvl['n_bb'] = 0
    if bbt is not None:
        lvl['n_bb'] = np.bincount(bbt.i, minlength=nk) + np.bincount(bbt.j, minlength=nk)
    lvl['has_bf'] = False
    if bft is not None:
        lvl.loc[np.unique(bft.i), 'has_bf'] = True
    lvl['n_col'] = 0
    if col is not None:
        lvl['n_col'] = np.bincount(col.i, minlength=nk) + np.bincount(col.j, minlength=nk)
    cont = continuum_level(atom)
    lvl['continuum'] = cont
    lvl['binding_ev'] = np.where(cont >= 0, lvl.ev.to_numpy()[np.maximum(cont, 0)] - lvl.ev, np.nan)

    findings = []
    bad = []
    for stage, grp in lvl.groupby('ion'):
        dif = np.diff(grp.ev.to_numpy())
        bad += list(grp.index[1:][dif < 0])
    if bad:
        findings.append(Finding('level order', 'serious',
                                f'{len(bad)} levels lower than their predecessor in the same ion',
                                len(bad), bad))
    dup = []
    for stage, grp in lvl.groupby('ion'):
        ecm = grp.ecm.to_numpy()
        gg = grp.g.to_numpy()
        idx = grp.index.to_numpy()
        for k in np.where(np.abs(np.diff(ecm)) < energy_tolerance)[0]:
            dup.append((int(idx[k]), int(idx[k + 1]), float(ecm[k]), float(gg[k]), float(gg[k + 1])))
    if dup:
        findings.append(Finding('duplicate levels', 'serious',
                                f'{len(dup)} pairs of levels within {energy_tolerance} cm^-1 in one ion',
                                len(dup), dup))
    badg = lvl.index[(lvl.g < 1) | (np.abs(lvl.g - np.round(lvl.g)) > 1e-6)].tolist()
    if badg:
        findings.append(Finding('statistical weights', 'serious',
                                f'{len(badg)} levels with g < 1 or non-integer g', len(badg), badg))

    top = lvl.ion == lvl.ion.max()
    no_bb = lvl.index[(lvl.n_bb == 0) & ~top].tolist()
    if no_bb:
        findings.append(Finding('levels without lines', 'minor',
                                f'{len(no_bb)} levels below the top ion have no bound-bound transition',
                                len(no_bb), no_bb))
    no_bf = lvl.index[~lvl.has_bf & (cont >= 0)].tolist()
    if no_bf:
        findings.append(Finding('levels without photoionisation', 'moderate',
                                f'{len(no_bf)} ionisable levels have no bound-free table '
                                '(dispatch deactivates them unless kramer_bf)', len(no_bf), no_bf))
    isolated = lvl.index[(lvl.n_bb == 0) & ~lvl.has_bf & (lvl.n_col == 0) & ~top].tolist()
    if isolated:
        findings.append(Finding('isolated levels', 'serious',
                                f'{len(isolated)} levels have no radiative or collisional connection',
                                len(isolated), isolated))
    above = lvl.index[(lvl.binding_ev <= 0)].tolist()
    if above:
        findings.append(Finding('levels above the limit', 'serious',
                                f'{len(above)} levels at or above the ionisation limit of their ion',
                                len(above), above))
    stages = lvl.groupby('ion').size()
    findings.append(Finding('ion stages', 'info',
                            ', '.join(f'ion {s}: {n} levels' for s, n in stages.items())))
    return lvl, findings


def check_limits(atom, nist_limits, hop_limits=None, tolerance_ev=0.005):
    """Ionisation energies of the atom against NIST (and the EOS of HOP).

    ``nist_limits`` and ``hop_limits`` map an ion stage to its ionisation energy [eV].
    An offset changes the LTE ratio of the stages by exp(offset/kT): 0.01 eV is 2% at
    5800 K. dispatch only warns when the EOS value differs by more than 10%.
    """
    grounds = ion_ground(atom)
    rows = []
    findings = []
    for stage in sorted(grounds)[:-1]:
        lim = atom.lvl.ev.iloc[grounds[stage + 1]] - atom.lvl.ev.iloc[grounds[stage]]
        row = dict(ion=stage, atom_ev=lim, nist_ev=nist_limits.get(stage, np.nan))
        if hop_limits is not None:
            row['hop_ev'] = hop_limits.get(stage, np.nan)
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab['offset_ev'] = tab.atom_ev - tab.nist_ev
    fractions = stage_fractions(atom, {r.ion: r.nist_ev for r in tab.itertuples()})
    tab['upper_fraction'] = [fractions.get(r.ion + 1, np.nan) for r in tab.itertuples()]
    for r in tab.itertuples():
        if abs(r.offset_ev) > tolerance_ev:
            sev = 'serious' if abs(r.offset_ev) > 0.05 else 'moderate'
            # the offset matters through the transitions into the upper stage, and those
            # only as far as that stage holds the element
            if r.upper_fraction < 1e-3:
                sev = 'minor'
            findings.append(Finding('ionisation energy', sev,
                                    f'ion {r.ion}: {r.atom_ev:.4f} eV in the atom, {r.nist_ev:.4f} eV '
                                    f'in NIST ({r.offset_ev:+.4f} eV, LTE stage ratio x'
                                    f'{np.exp(-r.offset_ev / (BK * 5800 / EE)):.3f} at 5800 K); ion {r.ion + 1} '
                                    f'holds {r.upper_fraction:.2g} of the element at {STAGE_TEMP:.0f} K, '
                                    f'Pe = {STAGE_PE:g} dyn/cm2'))
        if hop_limits is not None and abs(r.atom_ev - r.hop_ev) > tolerance_ev:
            findings.append(Finding('ionisation energy vs EOS', 'moderate',
                                    f'ion {r.ion}: atom {r.atom_ev:.4f} eV, HOP absdat {r.hop_ev:.4f} eV'))
    return tab, findings


STAGE_TEMP = 6500.0    # [K] and
STAGE_PE = 10.0        # [dyn/cm2]: a warm dwarf photosphere, to weigh ion stages


def stage_fractions(atom, limits, temp=STAGE_TEMP, pe=STAGE_PE):
    """LTE (Saha) fraction of the element in each ion stage of the atom, with the
    atom's own partition functions and the ionisation energies ``limits`` (stage ->
    eV) at ``temp`` and electron pressure ``pe``."""
    zz = partition_functions(atom.lvl, [temp]).iloc[0]
    ne = pe / (BK * temp)
    saha = 2 * (2 * PI * EM * BK * temp / HH**2)**1.5 / ne
    stages = sorted(zz.index)
    rel = [1.0]
    for s in stages[:-1]:
        rel.append(rel[-1] * zz[s + 1] / zz[s] * saha * np.exp(-limits.get(s, np.inf) * EE / (BK * temp)))
    rel = np.array(rel) / np.sum(rel)
    return dict(zip(stages, rel))


def partition_functions(lvl, temps, reference=None):
    """Partition function of each ion stage over ``temps`` from a level table with
    columns ev, g and ion. Energies count from ``reference`` (ion stage -> energy of
    its ground state on the table's scale), by default from each stage's lowest level;
    pass the reference of the unmerged atom when the ground term is a superlevel."""
    temps = np.atleast_1d(np.asarray(temps, float))
    out = {}
    for stage, grp in lvl.groupby('ion'):
        ev = grp.ev.to_numpy() - (grp.ev.min() if reference is None else reference[stage])
        out[int(stage)] = (grp.g.to_numpy()[:, None] * np.exp(-ev[:, None] * EE / (BK * temps[None, :]))).sum(0)
    return pd.DataFrame(out, index=pd.Index(temps, name='temp'))


def check_completeness(atom, nist_levels, temps=(4000., 5000., 6000., 8000., 10000.),
                       top_margin_ev=0.0):
    """Levels and partition function of each ion against a NIST level table.

    ``nist_levels`` has columns ion, ev (relative to the stage's ground), g and term, with
    one row per fine-structure level below the limit. The atom's partition function of
    each stage, relative to the NIST sum, is the fraction of the stage's population that
    the atom represents in LTE. NIST is incomplete at high excitation for complex ions,
    so the ratio can exceed 1 where the atom carries predicted levels.
    """
    zatom = partition_functions(atom.lvl, temps)
    nl = nist_levels[nist_levels.ion.isin(zatom.columns)]
    znist = partition_functions(nl, temps)
    ratio = zatom / znist
    findings = []
    for stage in ratio.columns:
        worst = ratio[stage].min()
        top = ratio[stage].max()
        sev = 'info'
        if worst < 0.9 or top > 1.1:
            sev = 'moderate'
        if worst < 0.5:
            sev = 'serious'
        n_atom = int((atom.lvl.ion == stage).sum())
        g_atom = atom.lvl.g[atom.lvl.ion == stage].sum()
        g_nist = nl.g[nl.ion == stage].sum()
        findings.append(Finding('partition function', sev,
                                f'ion {stage}: U_atom/U_NIST = {ratio[stage].iloc[0]:.3f} .. '
                                f'{ratio[stage].iloc[-1]:.3f} for T = {temps[0]:.0f} .. {temps[-1]:.0f} K; '
                                f'{n_atom} levels with sum g = {g_atom:.0f}, NIST sum g = {g_nist:.0f}'))
    return ratio, findings


def match_levels(atom, nist_levels, tolerance_cm=50.0):
    """Assign each NIST fine-structure level to the atom level of the same stage that
    holds it: the nearest level in energy within ``tolerance_cm`` (a merged term or a
    superlevel holds several). Returns the NIST table with the assigned index (-1 if
    none) and the energy offset."""
    nl = nist_levels.copy()
    nl['level'] = -1
    nl['offset_cm'] = np.nan
    for stage, grp in nl.groupby('ion'):
        sel = atom.lvl.index[atom.lvl.ion == stage]
        if len(sel) == 0:
            continue
        ground = atom.lvl.ev[sel].min()
        ecm = (atom.lvl.ev[sel].to_numpy() - ground) / CM_TO_EV
        nist_cm = grp.ev.to_numpy() / CM_TO_EV
        k = np.abs(nist_cm[:, None] - ecm[None, :]).argmin(1)
        off = nist_cm - ecm[k]
        ok = np.abs(off) <= tolerance_cm
        nl.loc[grp.index, 'level'] = np.where(ok, sel.to_numpy()[k], -1)
        nl.loc[grp.index, 'offset_cm'] = off
    return nl


def check_nist_levels(atom, nist_levels, tolerance_cm=50.0, below_limit_ev=None):
    """NIST levels that no atom level holds, and atom levels whose statistical weight
    differs from the NIST levels assigned to them.

    ``below_limit_ev`` restricts the missing-level search to levels bound by more than
    this (the atom is not expected to hold the levels just below the limit).
    """
    nl = match_levels(atom, nist_levels, tolerance_cm)
    lim = {s: g.ev.max() for s, g in nist_levels.groupby('ion')}
    findings = []
    if below_limit_ev is not None and 'limit_ev' in nl:
        nl['bound_ev'] = nl.limit_ev - nl.ev
        cand = nl[(nl.level < 0) & (nl.bound_ev > below_limit_ev)]
    else:
        cand = nl[nl.level < 0]
    for stage, grp in cand.groupby('ion'):
        if stage not in set(atom.lvl.ion):
            continue
        top = atom.lvl.ev[atom.lvl.ion == stage].max() - atom.lvl.ev[atom.lvl.ion == stage].min()
        inside = grp[grp.ev < top]
        if len(inside):
            findings.append(Finding('missing NIST levels', 'moderate',
                                    f'ion {stage}: {len(inside)} NIST levels below the atom\'s highest level '
                                    f'({top:.3f} eV) are not within {tolerance_cm:.0f} cm^-1 of an atom level',
                                    len(inside),
                                    [f"{r.term} {r.ev:.3f} eV" for r in inside.itertuples()]))
    gsum = nl[nl.level >= 0].groupby('level').g.sum()
    gat = atom.lvl.g.loc[gsum.index]
    over = gsum.index[gsum > gat * 1.0001]
    if len(over):
        findings.append(Finding('level weights vs NIST', 'minor',
                                f'{len(over)} atom levels are nearer to NIST levels of larger summed g '
                                'than their own (merged FS split across atom levels or misassignment)',
                                len(over),
                                [f"{i}: g={atom.lvl.g[i]:.0f} vs NIST {gsum[i]:.0f}" for i in over]))
    return nl, findings


# ---------------------------------------------------------------------------------------
# bound-bound transitions
# ---------------------------------------------------------------------------------------
def line_table(atom):
    """Bound-bound transitions with the quantities dispatch.x derives from them."""
    bbt = atom.bbt.copy()
    lvl = atom.lvl
    ev = lvl.ev.to_numpy()
    gg = lvl.g.to_numpy()
    bbt['ev_i'] = ev[bbt.i]
    bbt['ev_j'] = ev[bbt.j]
    bbt['g_i'] = gg[bbt.i]
    bbt['g_j'] = gg[bbt.j]
    bbt['ion_i'] = lvl.ion.to_numpy()[bbt.i]
    bbt['ion_j'] = lvl.ion.to_numpy()[bbt.j]
    bbt['lam_ritz_vac'] = HC_AA_EV / (bbt.ev_j - bbt.ev_i)
    # read_pandas: lam0 = air2vac(lam_air) and the upper energy follows from it
    bbt['lam_vac'] = air2vac(bbt.lam0.to_numpy())
    bbt['f'] = 10**bbt.loggf / bbt.g_i
    bbt['A'] = 6.671e15 * bbt.f * bbt.g_i / bbt.g_j / bbt.lam_vac**2
    return bbt


def check_lines(atom, lam_tolerance=1e-4, evi_tolerance=1e-3, classic=False):
    """Bound-bound transitions: level indices, frequency consistency, f-values,
    duplicates.

    For a ``classic`` atom read_classic derives line frequencies from the level energies,
    so the wavelength checks only matter once the atom is converted to the pandas format
    (whose read_pandas converts tabulated air wavelengths); they are then minor."""
    pandas_only = ' (classic format: dispatch uses the level energies; this applies after conversion to pandas)'
    pandas_only = pandas_only if classic else ''
    bbt = line_table(atom)
    findings = []
    nk = len(atom.lvl)
    bad = bbt.index[(bbt.i < 0) | (bbt.j >= nk) | (bbt.i >= bbt.j)].tolist()
    if bad:
        findings.append(Finding('line levels', 'serious',
                                f'{len(bad)} transitions with invalid or reversed level indices', len(bad), bad))
    cross = bbt.index[bbt.ion_i != bbt.ion_j].tolist()
    if cross:
        findings.append(Finding('line levels', 'serious',
                                f'{len(cross)} bound-bound transitions connect two ion stages',
                                len(cross), cross))
    rel = bbt.lam_vac / bbt.lam_ritz_vac - 1
    bbt['lam_offset'] = rel
    mis = bbt[np.abs(rel) > lam_tolerance]
    if len(mis):
        # an energy offset dE breaks detailed balance of the radiative rates in LTE by
        # exp(dE/kT): the line frequency and the level energies disagree
        dev = (HC_AA_EV / mis.lam_vac - HC_AA_EV / mis.lam_ritz_vac).abs()
        sev = 'serious' if (dev > 0.01).any() else 'moderate'
        sev = 'minor' if classic else sev
        uv = mis[mis.lam_ritz_vac < 2000]
        findings.append(Finding('line wavelength vs levels', sev,
                                f'{len(mis)} lines whose wavelength (as dispatch reads it) differs from the '
                                f'level energies by more than {lam_tolerance:.0e} (max {dev.max():.3f} eV; '
                                f'{len(uv)} of them below 2000 AA){pandas_only}', len(mis),
                                [f"{r.Index}: {r.lam_vac:.3f} vs Ritz {r.lam_ritz_vac:.3f}" for r in
                                 mis.reindex(dev.sort_values(ascending=False).index).itertuples()]))
    if 'evi' in bbt:
        off = (bbt.evi - bbt.ev_i).abs()
        badevi = bbt.index[off > evi_tolerance].tolist()
        if badevi:
            findings.append(Finding('line lower energy', 'moderate',
                                    f'{len(badevi)} lines whose evi column differs from the lower level '
                                    f'energy by more than {evi_tolerance} eV', len(badevi), badevi))
    nonpos = bbt.index[~np.isfinite(bbt.f) | (bbt.f <= 0)].tolist()
    if nonpos:
        findings.append(Finding('f-values', 'serious', f'{len(nonpos)} lines with f <= 0 or not finite',
                                len(nonpos), nonpos))
    big = bbt.index[bbt.f > 3].tolist()
    if big:
        findings.append(Finding('f-values', 'moderate', f'{len(big)} lines with f > 3',
                                len(big), [f"{i}: f={bbt.f[i]:.2f}" for i in big]))
    dup = bbt[bbt.duplicated(['i', 'j'], keep=False)]
    if len(dup):
        findings.append(Finding('duplicate lines', 'serious',
                                f'{len(dup)} rows share a level pair; dispatch adds both to the rates',
                                len(dup), sorted(set(zip(dup.i, dup.j)))))
    pole = bbt[((bbt.lam_ritz_vac > 1590) & (bbt.lam_ritz_vac < 1620))
               | ((bbt.lam_ritz_vac > 850) & (bbt.lam_ritz_vac < 905))]
    if len(pole):
        findings.append(Finding('air wavelengths near a pole', 'minor' if classic else 'serious',
                                f'{len(pole)} lines near 877 or 1603 AA, where the air/vacuum formula of '
                                f'physics.vac2air and handyfuncs_mod%air2vac diverges{pandas_only}', len(pole),
                                [f"{r.Index}: {r.lam_ritz_vac:.2f} AA" for r in pole.itertuples()]))
    nq = bbt.nq.clip(lower=0).sum()
    findings.append(Finding('line frequencies', 'info',
                            f'{len(bbt)} lines, sum of nq = {nq} (frequency-transition pairs before nq_fac)'))
    return bbt, findings


def nist_to_atom(atom, nist_levels, stage, term_window_ev=0.15, energy_tolerance_ev=0.003):
    """For each NIST level of one stage, the atom level that represents it: the level
    of the same LS symmetry nearest in energy within ``term_window_ev`` (a term or
    superlevel holds fine-structure levels some way from its mean energy), else the
    nearest level within ``energy_tolerance_ev``, else -1."""
    lvl = atom.lvl
    sel = lvl.index[lvl.ion == stage].to_numpy()
    ev = lvl.ev[sel].to_numpy() - lvl.ev[sel].min()
    terms = [label_term(str(lvl.label[i])) for i in sel]
    nist = nist_levels[nist_levels.ion == stage]
    out = np.full(len(nist), -1)
    for k, (e, t) in enumerate(zip(nist.ev.to_numpy(), nist.term.astype(str))):
        dist = np.abs(ev - e)
        term = label_term(t)
        same = np.array([term is not None and a == term for a in terms])
        cand = np.where(same & (dist <= term_window_ev))[0]
        if len(cand):
            out[k] = sel[cand[dist[cand].argmin()]]
        elif dist.min() <= energy_tolerance_ev:
            out[k] = sel[dist.argmin()]
    return nist.ev.to_numpy(), out


def compare_linelist(atom, llist, stage, energy_tolerance_ev=0.003, strong_loggf=None,
                     strength_temperature=5000.0, nist_levels=None):
    """Atom lines of one ion stage against a line list (columns lam0 [AA] in air above
    2000 AA and vacuum below, evi [eV], loggf).

    Each line-list transition is mapped to atom levels by its lower energy and the upper
    energy E_low + hc/lambda, within ``energy_tolerance_ev`` of an atom level of that
    stage (or inside the energy span of the fine-structure levels a merged level stands
    for, when the atom gives it as a term). Mapped transitions are summed per atom level
    pair (gf adds over components) and compared with the atom's gf. A list transition
    with both levels in the atom but no atom transition is missing; its LTE strength
    proxy, gf exp(-E_low/kT), ranks it against the atom's own lines.

    With ``nist_levels`` the line-list energies are first identified with NIST levels
    (within ``energy_tolerance_ev``) and these with atom levels by :func:`nist_to_atom`,
    so that the components of a term-merged atom level find it.
    """
    lvl = atom.lvl
    sel = lvl.index[lvl.ion == stage].to_numpy()
    ground = lvl.ev[sel].min()
    ev = lvl.ev[sel].to_numpy() - ground
    from .atomic_data import ts_vacuum_wavelength
    lam_vac = ts_vacuum_wavelength(llist.lam0.to_numpy())
    lo = llist.evi.to_numpy()
    up = lo + HC_AA_EV / lam_vac

    def assign(energy):
        k = np.abs(energy[:, None] - ev[None, :]).argmin(1)
        ok = np.abs(energy - ev[k]) <= energy_tolerance_ev
        out = np.where(ok, sel[k], -1)
        if nist_levels is not None and (nist_levels.ion == stage).any():
            nist_ev, nist_atom = nist_to_atom(atom, nist_levels, stage, energy_tolerance_ev=energy_tolerance_ev)
            kn = np.abs(energy[:, None] - nist_ev[None, :]).argmin(1)
            hit = np.abs(energy - nist_ev[kn]) <= energy_tolerance_ev
            out = np.where(hit & (nist_atom[kn] >= 0), nist_atom[kn], out)
        return out

    tab = pd.DataFrame(dict(lam0=llist.lam0.to_numpy(), evi=lo, loggf=llist.loggf.to_numpy(),
                            i=assign(lo), j=assign(up)))
    tab['gf'] = 10**tab.loggf
    tab['strength'] = tab.gf * np.exp(-tab.evi * EE / (BK * strength_temperature))
    mapped = tab[(tab.i >= 0) & (tab.j >= 0) & (tab.i != tab.j)].copy()
    lo_hi = np.sort(mapped[['i', 'j']].to_numpy(), axis=1)
    mapped['i'] = lo_hi[:, 0]
    mapped['j'] = lo_hi[:, 1]
    pairs = mapped.groupby(['i', 'j']).agg(gf_list=('gf', 'sum'), strength=('strength', 'sum'),
                                          n_components=('gf', 'size'), lam0=('lam0', 'mean'))
    bbt = line_table(atom)
    bbt = bbt[bbt.ion_i == stage]
    atom_pairs = bbt.groupby(['i', 'j']).agg(gf_atom=('loggf', lambda s: (10**s).sum()))
    comp = pairs.join(atom_pairs, how='outer')
    comp['dloggf'] = np.log10(comp.gf_atom) - np.log10(comp.gf_list)
    findings = []
    both = comp.dropna(subset=['gf_atom', 'gf_list'])
    off = both[both.dloggf.abs() > 0.3]
    # a systematic offset (gf stored as f, log gf of another source) rather than scatter
    sev = 'serious' if len(both) >= 5 and abs(both.dloggf.median()) > 0.3 else ('moderate' if len(off) else 'info')
    findings.append(Finding('f-values vs line list', sev,
                            f'ion {stage}: {len(both)} atom transitions found in the line list; '
                            f'{len(off)} differ by more than 0.3 dex in gf (median {both.dloggf.median():+.3f} dex)',
                            len(off), [f"{i}-{j} {r.lam0:.2f} AA: {r.dloggf:+.2f}" for (i, j), r in
                                       off.reindex(off.dloggf.abs().sort_values(ascending=False).index).iterrows()]))
    missing = comp[comp.gf_atom.isna()].sort_values('strength', ascending=False)
    if strong_loggf is None and len(bbt):
        ref = (10**bbt.loggf * np.exp(-bbt.ev_i.sub(ground) * EE / (BK * strength_temperature)))
        strong = ref.quantile(0.5)
    else:
        strong = 10**(strong_loggf or 0)
    big = missing[missing.strength > strong]
    sev = 'serious' if len(big) > 10 else ('moderate' if len(big) else 'info')
    findings.append(Finding('missing strong lines', sev,
                            f'ion {stage}: {len(missing)} line-list level pairs with both levels in the atom '
                            f'have no atom transition; {len(big)} of them are stronger than the median atom line',
                            len(big), [f"{i}-{j} {r.lam0:.2f} AA gf={r.gf_list:.3g}" for (i, j), r in big.iterrows()]))
    unmatched = int(((tab.i < 0) | (tab.j < 0)).sum())
    findings.append(Finding('line list mapping', 'info',
                            f'ion {stage}: {len(tab)} line-list transitions, {unmatched} with a level '
                            f'outside the atom (within {energy_tolerance_ev} eV)'))
    return comp, findings


def level_keys(atom):
    """dispatch's departure keys (atom_mod%ev_keys): 100 * stage plus the energy above
    the ground level of the level's stage [eV]."""
    lvl = atom.lvl
    grounds = ion_ground(atom)
    ev0 = np.array([lvl.ev[grounds[int(s)]] for s in lvl.ion])
    return 100 * lvl.ion.to_numpy(float) + lvl.ev.to_numpy(float) - ev0


def check_level_matching(atom, nist_levels, llist, stage, strong_loggf=-1.0, tolerance_ev=0.003):
    """Which atom levels the line list's lines take their departure coefficients from.

    An NLTE synthesis gives each line-list line the departure coefficients of the atom
    levels nearest in energy to its lower and upper level (``find_closest`` in
    HOP/linelist_mod.f90, without a tolerance). For the lines with loggf above
    ``strong_loggf`` this finds the NIST level of each end (within ``tolerance_ev``) and
    reports ends whose matched atom level has a different LS symmetry
    (:func:`label_term`); a level of the other multiplicity within 0.002 eV counts as a
    match, since singlets and triplets of high l are degenerate and often merged."""
    from .atomic_data import ts_vacuum_wavelength
    keys = level_keys(atom)
    nist = nist_levels[nist_levels.ion == stage]
    lines = llist[llist.loggf > strong_loggf]
    ev_lo = lines.evi.to_numpy(float)
    ev_up = ev_lo + HC_AA_EV / ts_vacuum_wavelength(lines.lam0)
    terms = [label_term(t) for t in nist.term]
    rows = []
    for end, evs in (('lower', ev_lo), ('upper', ev_up)):
        near = np.abs(keys[None, :] - (100 * stage + evs)[:, None]).argmin(1)
        dn = np.abs(nist.ev.to_numpy()[None, :] - evs[:, None])
        kn = dn.argmin(1)
        for k, (lam, ev, i, n) in enumerate(zip(lines.lam0, evs, near, kn)):
            if dn[k, n] > tolerance_ev or int(atom.lvl.ion[i]) != stage:
                continue
            atom_term = label_term(str(atom.lvl.label[i]))
            rows.append(dict(lam0=lam, loggf=lines.loggf.iloc[k], end=end, ev=ev, level=i,
                             level_label=str(atom.lvl.label[i]).strip(), offset_ev=keys[i] - 100 * stage - ev,
                             nist_term=nist.term.iloc[n],
                             # singlet and triplet of a high-l configuration are degenerate
                             same_term=atom_term == terms[n] or (
                                 atom_term is not None and terms[n] is not None
                                 and atom_term[1:] == terms[n][1:] and abs(keys[i] - 100 * stage - ev) < 0.002)))
    tab = pd.DataFrame(rows)
    findings = []
    # without labels: how far the level whose departure an end takes lies from it
    sel = atom.lvl.ion.to_numpy() == stage
    top = keys[sel].max() - 100 * stage
    offs = []
    for evs in (ev_lo, ev_up):
        near = keys[sel][np.abs(keys[sel][None, :] - (100 * stage + evs)[:, None]).argmin(1)]
        offs.append(near - 100 * stage - evs)
    offs = np.abs(np.concatenate(offs))
    above = int((np.concatenate([ev_lo, ev_up]) > top + 0.01).sum())
    if len(offs):
        frac10 = (offs > 0.1).mean()
        sev = 'serious' if frac10 > 0.05 else ('moderate' if (offs > 0.01).mean() > 0.05 else 'info')
        findings.append(Finding('line-list level offsets', sev,
                                f'ion {stage}: of {len(offs)} ends of the {len(lines)} lines with loggf > '
                                f'{strong_loggf}, {(offs > 0.01).mean():.1%} take the departure of a level more '
                                f'than 0.01 eV away and {frac10:.1%} more than 0.1 eV; {above} ends lie above '
                                f'the atom\'s highest level of the stage ({top:.3f} eV)', int((offs > 0.1).sum())))
    if len(tab):
        bad = tab[~tab.same_term & tab.level_label.map(lambda l: label_term(l) is not None)
                  & tab.nist_term.map(lambda t: label_term(str(t)) is not None)]
        sev = 'serious' if len(bad) else 'info'
        findings.append(Finding('line-list level matching', sev,
                                f'ion {stage}: {len(tab)} ends of {len(lines)} lines with loggf > {strong_loggf} '
                                f'have a NIST level; {len(bad)} take the departure of a level of another LS term '
                                f'(median |offset| of those {bad.offset_ev.abs().median() if len(bad) else 0:.3f} eV)',
                                len(bad), [f'{r.lam0:.2f} AA {r.end} {r.nist_term} {r.ev:.3f} eV -> '
                                           f'{r.level_label} ({r.offset_ev:+.3f} eV)'
                                           for r in bad.sort_values('loggf', ascending=False).itertuples()]))
    return tab, findings


def check_comments(atom_file, lam_tolerance=1e-3, loggf_tolerance=0.05):
    """Lines of a classic MULTI atom whose comment (the source line list, as most atom
    builders write it) disagrees with the data dispatch.x uses: the Ritz wavelength of
    the level energies and the f-value. A disagreement marks a misassigned level, a
    merged multiplet whose level energies or f-value were not updated, or a vacuum
    wavelength in the comment."""
    import contextlib
    import io
    from .atom import AtomClass
    with contextlib.redirect_stdout(io.StringIO()):
        atom = AtomClass(atom_file, qnorm=10., read_bf=False, read_coll=False, sanity_check=True)
    rows = []
    for t in atom.trans[:atom.nline]:
        rows.append(dict(kr=t.kr, i=t.i, j=t.j, lam_ritz=t.lam0, lam_comment=t.c_lam0,
                         loggf=np.log10(t.fosc * t.lower.g), loggf_comment=t.c_loggf))
    tab = pd.DataFrame(rows)
    tab['lam_offset'] = tab.lam_comment / tab.lam_ritz - 1
    tab['dloggf'] = tab.loggf - tab.loggf_comment
    findings = []
    lam = tab[tab.lam_offset.abs() > lam_tolerance]
    if len(lam):
        uv = lam[lam.lam_ritz < 2000]
        findings.append(Finding('comment wavelength vs levels', 'moderate',
                                f'{len(lam)} lines whose comment wavelength differs from the Ritz wavelength by '
                                f'more than {lam_tolerance:.0e} ({len(uv)} below 2000 AA, where comments may be '
                                'vacuum)', len(lam),
                                [f"{r.kr} ({r.i}-{r.j}): Ritz {r.lam_ritz:.2f} vs {r.lam_comment:.2f} AA" for r in
                                 lam.reindex(lam.lam_offset.abs().sort_values(ascending=False).index).itertuples()]))
    gf = tab[tab.dloggf.abs() > loggf_tolerance]
    if len(gf):
        findings.append(Finding('comment loggf vs f-value', 'moderate',
                                f'{len(gf)} lines whose f-value differs from the comment loggf by more than '
                                f'{loggf_tolerance} dex', len(gf),
                                [f"{r.kr} ({r.i}-{r.j}) {r.lam_ritz:.2f} AA: {r.dloggf:+.2f}" for r in
                                 gf.reindex(gf.dloggf.abs().sort_values(ascending=False).index).itertuples()]))
    return tab, findings


# ---------------------------------------------------------------------------------------
# bound-free transitions
# ---------------------------------------------------------------------------------------
def bf_tables(atom):
    """Tabulated continua as {lower level: (lam [AA] increasing, alpha [cm^2])}."""
    out = {}
    if not hasattr(atom, 'bft'):
        return out
    for i, grp in atom.bft.groupby('i'):
        lam = grp.lam_air.to_numpy(float)
        alpha = grp.alpha.to_numpy(float)
        order = np.argsort(lam, kind='stable')
        lam = lam[order]
        # read_pandas and read_classic separate repeated wavelengths by 1e-6 AA
        for k in range(1, len(lam)):
            if lam[k] <= lam[k - 1]:
                lam[k:] += lam[k - 1] - lam[k] + 1e-6
        out[int(i)] = (lam, alpha[order])
    return out


def bf_rate_integral(lam_nodes, alpha_nodes, temp, lam_eval):
    """Photoionisation rate in a Planck field up to a constant, with the cross section
    as dispatch interpolates it onto ``lam_eval`` and integrated in frequency."""
    alpha = bf_interpolate(lam_nodes, alpha_nodes, lam_eval)
    nu = CC * 1e8 / lam_eval
    order = np.argsort(nu)
    return np.trapezoid((alpha * planck_photon_weight(nu, temp))[order], nu[order])


def dense_grid(lam_nodes, per_interval=20):
    """Wavelengths that subdivide every tabulation interval, for evaluating an
    interpolant between its nodes."""
    lam = np.sort(np.asarray(lam_nodes, float))
    frac = np.linspace(0, 1, per_interval, endpoint=False)
    inner = (lam[:-1, None] + np.diff(lam)[:, None] * frac[None, :]).ravel()
    return np.append(inner, lam[-1])


def leave_one_out_error(lam_nodes, alpha_nodes, temp):
    """Interpolation error of a table at its own resolution: each interior node is
    predicted from the table without it. Returns the largest relative error at a node
    and the relative change of the Planck-weighted rate when every second node is
    dropped (a table that resolves its shape changes little)."""
    lam = np.asarray(lam_nodes, float)
    alpha = np.asarray(alpha_nodes, float)
    n = len(lam)
    if n < 4:
        return np.nan, np.nan
    errs = []
    for k in range(1, n - 1):
        keep = np.arange(n) != k
        pred = bf_interpolate(lam[keep], alpha[keep], lam[k:k + 1])[0]
        errs.append(abs(pred - alpha[k]) / max(alpha[k], 1e-40))
    grid = dense_grid(lam)
    full = bf_rate_integral(lam, alpha, temp, grid)
    half = np.r_[np.arange(0, n - 1, 2), n - 1]
    thin = bf_rate_integral(lam[half], alpha[half], temp, grid)
    return float(np.max(errs)), float(thin / full - 1) if full > 0 else np.nan


def check_bf(atom, temp=5800.0, threshold_tolerance=0.02):
    """Bound-free tables: threshold position and value, coverage towards the blue,
    resonance structure and how well dispatch's interpolation reproduces the table.

    The threshold value is compared with the hydrogenic cross section of the level's
    binding energy (Gaunt factor 1); orders of magnitude apart flag a unit error or a
    misassigned table. The blue coverage is the fraction of the Planck-weighted
    hydrogenic rate (at ``temp``) that lies beyond the table's blue end, where dispatch
    sets the cross section to zero.
    """
    lvl = atom.lvl
    cont = continuum_level(atom)
    tables = bf_tables(atom)
    rows = []
    for i, (lam, alpha) in tables.items():
        stage = int(lvl.ion[i])
        jc = cont[i]
        bind = lvl.ev[jc] - lvl.ev[i] if jc >= 0 else np.nan
        lam_edge = HC_AA_EV / bind if jc >= 0 else np.nan
        row = dict(i=i, ion=stage, label=lvl.label[i], npts=len(lam), lam_edge=lam_edge,
                   lam_red=lam[-1], lam_blue=lam[0], alpha_red=alpha[-1], alpha_max=alpha.max())
        row['edge_offset'] = lam[-1] / lam_edge - 1
        row['kramers_ratio'] = alpha[-1] / kramers(lam_edge, bind, stage) if jc >= 0 else np.nan
        # Planck-weighted hydrogenic rate missing beyond the blue end of the table
        nu0 = CC * 1e8 / lam_edge
        nu_blue = CC * 1e8 / lam[0]
        nu = np.geomspace(nu0, nu0 * 50, 4000)
        weight = planck_photon_weight(nu, temp) / nu**3
        total = np.trapezoid(weight, nu)
        beyond = np.trapezoid(np.where(nu > nu_blue, weight, 0), nu)
        row['blue_missing'] = beyond / total if total > 0 else np.nan
        nnu = CC * 1e8 / lam
        row['n_maxima'] = int(((alpha[1:-1] > alpha[:-2]) & (alpha[1:-1] > alpha[2:])).sum())
        row['max_jump'] = float(np.max(np.abs(np.diff(np.log(np.maximum(alpha, 1e-40)))))) if len(alpha) > 1 else 0
        row['dnu_min_rel'] = float(np.min(np.abs(np.diff(nnu))) / nnu.max()) if len(lam) > 1 else np.nan
        loo, thin = leave_one_out_error(lam, alpha, temp)
        row['loo_max'] = loo
        row['rate_halved'] = thin
        row['nonpositive'] = int((alpha <= 0).sum())
        rows.append(row)
    tab = pd.DataFrame(rows)
    findings = []
    if tab.empty:
        findings.append(Finding('bound-free tables', 'serious', 'the atom has no bound-free tables'))
        return tab, findings
    mis = tab[tab.edge_offset.abs() > threshold_tolerance]
    if len(mis):
        red = mis[mis.edge_offset < 0]
        findings.append(Finding('bf threshold vs levels', 'serious' if len(red) else 'moderate',
                                f'{len(mis)} tables end more than {threshold_tolerance:.0%} from the edge '
                                f'given by the level energies; {len(red)} stop short of the edge, so the '
                                'region next to the threshold has no cross section', len(mis),
                                [f"{r.i} {r.label.strip()}: table {r.lam_red:.1f} vs edge {r.lam_edge:.1f} AA"
                                 for r in mis.reindex(mis.edge_offset.abs().sort_values(ascending=False).index).itertuples()]))
    kr = tab[(tab.kramers_ratio < 1e-2) | (tab.kramers_ratio > 1e2)]
    if len(kr):
        findings.append(Finding('bf threshold value', 'moderate',
                                f'{len(kr)} threshold cross sections differ from hydrogenic by more than 100x',
                                len(kr), [f"{r.i} {r.label.strip()}: {r.kramers_ratio:.2g}" for r in kr.itertuples()]))
    blue = tab[tab.blue_missing > 0.05]
    if len(blue):
        findings.append(Finding('bf blue coverage', 'serious' if (tab.blue_missing > 0.2).any() else 'moderate',
                                f'{len(blue)} tables stop so close to the edge that more than 5% of a hydrogenic '
                                f'rate at {temp:.0f} K lies beyond them (max {tab.blue_missing.max():.0%})', len(blue),
                                [f"{r.i} {r.label.strip()}: {r.blue_missing:.0%}" for r in
                                 blue.sort_values('blue_missing', ascending=False).itertuples()]))
    loo = tab[tab.rate_halved.abs() > 0.05]
    if len(loo):
        findings.append(Finding('bf sampling', 'moderate',
                                f'{len(loo)} tables change their Planck-weighted rate by more than 5% when every '
                                'second point is dropped: the tabulation does not resolve the shape it encodes',
                                len(loo), [f"{r.i} {r.label.strip()}: {r.rate_halved:+.1%}" for r in
                                           loo.reindex(loo.rate_halved.abs().sort_values(ascending=False).index).itertuples()]))
    nonpos = tab[tab.nonpositive > 0]
    if len(nonpos):
        findings.append(Finding('bf values', 'moderate', f'{len(nonpos)} tables contain zero or negative values',
                                len(nonpos), nonpos.i.tolist()))
    findings.append(Finding('bound-free tables', 'info',
                            f'{len(tab)} tables, {int(tab.npts.sum())} points (median {int(tab.npts.median())})'))
    return tab, findings


def compare_bf_reference(lam_table, alpha_table, lam_ref, alpha_ref, temps=(4000., 6000., 10000.)):
    """A tabulated cross section against the original data it was made from.

    The original is taken as given (its own nodes, linearly between them); the table as
    dispatch interpolates it. Returns the relative difference of the Planck-weighted
    rate at each temperature, integrated over the overlap of the two, and the fraction
    of the original's rate outside the table's range."""
    lam_ref = np.asarray(lam_ref, float)
    order = np.argsort(lam_ref)
    lam_ref = lam_ref[order]
    alpha_ref = np.asarray(alpha_ref, float)[order]
    grid = np.union1d(lam_ref, dense_grid(lam_table))
    nu = CC * 1e8 / grid
    ref = np.interp(grid, lam_ref, alpha_ref, left=0, right=0)
    tab = bf_interpolate(lam_table, alpha_table, grid)
    inside = (grid >= np.min(lam_table)) & (grid <= np.max(lam_table))
    out = {}
    order = np.argsort(nu)
    for t in temps:
        ww = planck_photon_weight(nu, t)
        r_ref = np.trapezoid((ref * ww)[order], nu[order])
        r_in = np.trapezoid((np.where(inside, ref, 0) * ww)[order], nu[order])
        r_tab = np.trapezoid((tab * ww)[order], nu[order])
        out[f'rate_{t:.0f}'] = r_tab / r_ref - 1
        out[f'outside_{t:.0f}'] = 1 - r_in / r_ref
    return out


L_LETTERS = 'SPDFGHIKLMNOQ'


def label_term(label):
    """(2S+1, L, parity) of the last LS term in a level label such as '3s3p 3P*' or
    'Level 43 = 7d1D', or None."""
    import re
    terms = re.findall(r'(\d)([SPDFGHIKLMNOQ])\d*(\*?)', label)
    if not terms:
        return None
    mult, ell, odd = terms[-1]
    return int(mult), L_LETTERS.index(ell), int(odd == '*')


def term_ranks(atom, stage, fine_structure_ev=0.02):
    """Rank of each level among the terms of its LS symmetry in one ion stage: levels
    of one term (same symmetry, within ``fine_structure_ev`` of the previous one) share
    a rank, and ranks count terms upwards from 1 as TOPbase's ilv does."""
    lvl = atom.lvl[atom.lvl.ion == stage]
    ranks = {}
    last = {}
    for i, r in lvl.sort_values('ev').iterrows():
        term = label_term(r.label)
        if term is None:
            continue
        rank, ev = last.get(term, (0, -np.inf))
        if r.ev - ev > fine_structure_ev:
            rank += 1
        last[term] = (rank, r.ev)
        ranks[i] = (term, rank)
    return ranks


def compare_bf_topbase(atom, topbase, stage, temps=(4000., 6000., 10000.), threshold_tolerance=0.25):
    """Each tabulated continuum of an ion stage against the TOPbase cross section of the
    same LS term (read by :func:`tumult.atomic_data.read_topbase_photo`).

    Atom levels are matched to TOPbase by symmetry and rank (:func:`term_ranks`), which
    assumes that the atom holds every lower term of the same symmetry; matches whose
    thresholds differ by more than ``threshold_tolerance`` (relative) are dropped. The
    TOPbase cross section is shifted so that its threshold falls on the atom's edge (OP
    energies are theoretical), and the table as dispatch interpolates it is compared
    with it by :func:`compare_bf_reference`."""
    from .atomic_data import RYD_AA, RYD_EV, MB
    lvl = atom.lvl
    cont = continuum_level(atom)
    tables = bf_tables(atom)
    ranks = term_ranks(atom, stage)
    rows = []
    for i, (lam, alpha) in tables.items():
        if i not in ranks or cont[i] < 0:
            continue
        term, rank = ranks[i]
        key = (*term, rank)
        if key not in topbase:
            continue
        energy, sigma = topbase[key]
        bind_ry = (lvl.ev[cont[i]] - lvl.ev[i]) / RYD_EV
        off = energy[0] / bind_ry - 1
        if abs(off) > threshold_tolerance:
            continue
        lam_ref = RYD_AA / (energy - energy[0] + bind_ry)
        row = dict(i=i, label=lvl.label[i].strip(), topbase=key, threshold_offset=off,
                   npts=len(lam), npts_ref=len(energy))
        row.update(compare_bf_reference(lam, alpha, lam_ref, sigma * MB, temps))
        rows.append(row)
    tab = pd.DataFrame(rows)
    findings = []
    if tab.empty:
        return tab, findings
    key = f'rate_{temps[1]:.0f}'
    bad = tab[tab[key].abs() > 0.2]
    sev = 'serious' if len(bad) > 0.2 * len(tab) else ('moderate' if len(bad) else 'info')
    findings.append(Finding('bf vs TOPbase', sev,
                            f'ion {stage}: {len(tab)} tables matched to TOPbase terms; the Planck-weighted rate at '
                            f'{temps[1]:.0f} K of the table as dispatch interpolates it differs by more than 20% '
                            f'for {len(bad)} (median {tab[key].median():+.0%})', len(bad),
                            [f"{r.i} {r.label}: {getattr(r, key):+.0%}" for r in
                             bad.reindex(bad[key].abs().sort_values(ascending=False).index).itertuples()]))
    return tab, findings


# ---------------------------------------------------------------------------------------
# collisions
# ---------------------------------------------------------------------------------------
def collision_temperatures(atom):
    return np.array([float(c) for c in atom.col.columns[3:]])


def exothermic_coefficient(key, values, temps, dev, g_lower, g_upper):
    """The tabulated quantity of a gencol key turned into a quantity that should vary
    slowly with temperature: the rate coefficient of the exothermic direction, up to
    temperature-independent factors (Saha factors of ionising processes cancel or are
    removed). A strong power-law trend points to an inverted table or a unit error.

    ``dev`` is the energy gap [eV]; for CH0 the H- affinity is subtracted.
    """
    tt = np.asarray(temps, float)
    vv = np.asarray(values, float)
    boltz = np.exp(dev * EE / (BK * tt))
    if key == 'CE':
        return g_lower / g_upper * vv * np.sqrt(tt)
    if key in ('OHM', 'OHMEGA'):
        return vv / np.sqrt(tt)
    if key == 'CH':
        return vv
    if key == 'CI':
        # cup = ct sqrt(T) exp(-dE/kT); three-body recombination ~ cup exp(dE/kT) T^-3/2
        return vv / tt
    if key == 'CH0':
        de = dev - H_MINUS_EV
        return vv * np.exp(max(de, 0) * EE / (BK * tt)) if de > 0 else vv
    if key == 'CHI':
        return vv * boltz / tt**1.5
    return vv


def van_regemorter(f, dev, temps):
    """Upward excitation coefficient per electron [cm^3/s] of van Regemorter (1962) in
    the form of atom_mod%classic_col."""
    tt = np.asarray(temps, float)
    c0 = PI * A0**2 * np.sqrt(8 * BK / (PI * EM))
    ryd_k = RYDV * EE / BK
    u0 = dev * EE / (BK * tt)
    pp = np.maximum(0.2 * np.exp(-u0), np.sqrt(3) / (2 * PI) * exp1(u0))
    return 14.5 * c0 * ryd_k / np.sqrt(tt) * RYDV / dev * f * pp


def check_ch_multiplicity(atom, col, vals, rtol=1e-6):
    """Hydrogen collision records that repeat one term-to-term rate for every
    fine-structure pair. gencol takes a CH value as the downward rate of that pair, so a
    copied term rate multiplies the total rate out of each upper level by the number of
    lower components; the term rate should be shared by statistical weight."""
    from .atom_build import fine_structure_groups
    group = np.empty(len(atom.lvl), int)
    for k, members in enumerate(fine_structure_groups(atom)):
        group[members] = k
    ch = col[col.key == 'CH'].copy()
    if not len(ch):
        return []
    ch['glo'] = group[ch.lo]
    ch['ghi'] = group[ch.hi]
    copied = []
    factors = []
    for (glo, ghi), grp in ch.groupby(['glo', 'ghi']):
        if len(grp) < 2 or glo == ghi:
            continue
        rows = vals[col.index.get_indexer(grp.index)]
        per_upper = grp.groupby('hi').size().mean()
        if per_upper > 1 and np.allclose(rows, rows[0], rtol=rtol, atol=0):
            copied.append((int(grp.lo.min()), int(grp.hi.min()), len(grp)))
            factors.append(per_upper)
    if not copied:
        return []
    return [Finding('hydrogen collision multiplicity', 'moderate',
                    f'{len(copied)} term pairs repeat one CH value for each of their fine-structure pairs; '
                    f'the downward rate out of an upper level is then {np.mean(factors):.1f}x the term rate '
                    f'on average (max {np.max(factors):.0f}x)', len(copied),
                    [f"{a}-{b} ({n} records)" for a, b, n in copied])]


def check_collisions(atom, slope_limit=2.5):
    """Collision data: coverage of the radiative transitions, duplicates, positivity,
    plausibility of the temperature dependence and agreement of CE with van Regemorter
    for allowed transitions.

    The slope is d ln q / d ln T of :func:`exothermic_coefficient`, fitted over the
    table; ``slope_limit`` bounds what physical exothermic coefficients do.
    """
    findings = []
    if not hasattr(atom, 'col'):
        findings.append(Finding('collisions', 'serious', 'the atom has no collision data'))
        return pd.DataFrame(), findings
    col = atom.col.copy()
    temps = collision_temperatures(atom)
    vals = col.iloc[:, 3:].to_numpy(float)
    lvl = atom.lvl
    ev = lvl.ev.to_numpy()
    gg = lvl.g.to_numpy()
    lo = np.minimum(col.i, col.j).to_numpy()
    hi = np.maximum(col.i, col.j).to_numpy()
    col['lo'] = lo
    col['hi'] = hi
    col['dev'] = ev[hi] - ev[lo]
    nonpos = col.index[(vals <= 0).any(1) | ~np.isfinite(vals).all(1)]
    if len(nonpos):
        findings.append(Finding('collision values', 'serious',
                                f'{len(nonpos)} records with zero, negative or non-finite values; read_pandas takes '
                                'their log (-inf), read_classic aborts', len(nonpos),
                                [f"{col.key[k]} {col.i[k]}-{col.j[k]}" for k in nonpos]))
    dup = col[col.duplicated(['key', 'lo', 'hi'], keep=False)]
    if len(dup):
        findings.append(Finding('duplicate collisions', 'serious',
                                f'{len(dup)} records repeat a key and level pair; gencol adds them, so these rates '
                                'are counted twice or more', len(dup),
                                sorted(set(zip(dup.key, dup.lo, dup.hi)))[:20]))
    ee = col[col.key.isin(['CE', 'OHM', 'OHMEGA'])]
    both = ee[ee.duplicated(['lo', 'hi'], keep=False) & ~ee.duplicated(['key', 'lo', 'hi'], keep=False)]
    if len(both):
        findings.append(Finding('duplicate collisions', 'serious',
                                f'{len(both)} level pairs carry both CE and OHM records (both are added)',
                                len(both)))
    cross = col[col.key.isin(['CE', 'OHM', 'OHMEGA', 'CH']) & (lvl.ion.to_numpy()[lo] != lvl.ion.to_numpy()[hi])]
    if len(cross):
        findings.append(Finding('collision levels', 'serious',
                                f'{len(cross)} excitation records connect two ion stages', len(cross),
                                [f"{r.key} {r.lo}-{r.hi}" for r in cross.itertuples()]))
    # temperature behaviour: the slope d ln q / d ln T over the whole table, and the
    # activation energy E_act of q ~ T^b exp(-E_act/kT) fitted over the photospheric
    # range. An exothermic coefficient has |E_act| of a few tenths of an eV at most;
    # E_act = -dE means the table carries a spurious exp(+dE/kT) (an upward rate stored
    # as a downward one, fe607a), E_act >> 0 a rate suppressed at low temperature.
    lnt = np.log(temps)
    win = (temps >= 2000) & (temps <= 20000)
    if win.sum() < 4:
        win = np.ones_like(temps, bool)
    design = np.vstack([np.ones(win.sum()), lnt[win], -EE / (BK * temps[win])]).T
    logq = np.full(vals.shape, np.nan)
    for k, r in enumerate(col.itertuples()):
        if np.all(vals[k] > 0):
            logq[k] = np.log(exothermic_coefficient(r.key, vals[k], temps, r.dev, gg[r.lo], gg[r.hi]))
    ok = np.isfinite(logq).all(1)
    col['slope'] = np.nan
    col['e_act'] = np.nan
    if ok.any():
        col.loc[ok, 'slope'] = np.polyfit(lnt, logq[ok].T, 1)[0]
        if len(temps) >= 4:
            col.loc[ok, 'e_act'] = np.linalg.lstsq(design, logq[ok][:, win].T, rcond=None)[0][2]
    findings.append(Finding('collision temperatures', 'info',
                            f'{len(temps)} temperatures {temps.min():.0f} .. {temps.max():.0f} K; '
                            'outside this range gencol holds the tabulated value constant'))
    for key, grp in col.groupby('key'):
        steep = grp[grp.slope.abs() > slope_limit]
        sev = 'info'
        if len(steep):
            sev = 'serious' if len(steep) > 0.1 * len(grp) else 'moderate'
        findings.append(Finding(f'collision T-dependence {key}', sev,
                                f'{len(grp)} records, exothermic-coefficient slope median {grp.slope.median():+.2f}, '
                                f'{len(steep)} beyond +-{slope_limit}', len(steep),
                                [f"{r.lo}-{r.hi} dE={r.dev:.2f} eV slope {r.slope:+.1f}" for r in
                                 steep.reindex(steep.slope.abs().sort_values(ascending=False).index).itertuples()]))
        act = grp[grp.e_act.notna()]
        upward = act[(act.dev > 0.3) & ((act.e_act + act.dev).abs() < 0.3 * act.dev)]
        if key in ('CE', 'OHM', 'OHMEGA', 'CH') and len(upward):
            findings.append(Finding(f'collision Boltzmann factor {key}', 'serious',
                                    f'{len(upward)} of {len(grp)} records (dE > 0.3 eV) have an activation '
                                    'energy within 30% of -dE: the table carries an extra exp(+dE/kT), as an '
                                    'upward rate stored where gencol expects a downward one would', len(upward),
                                    [f"{r.lo}-{r.hi} dE={r.dev:.2f} eV E_act={r.e_act:+.2f} eV" for r in
                                     upward.sort_values('dev', ascending=False).itertuples()]))
        hot = act[act.e_act > 0.5]
        if len(hot):
            dex = (hot.e_act * EE / (BK * 5000.0)).median() / np.log(10)
            findings.append(Finding(f'collision activation energy {key}', 'serious' if len(hot) > 0.1 * len(act)
                                    else 'moderate',
                                    f'{len(hot)} of {len(act)} records have an exothermic coefficient that '
                                    f'falls as exp(-E_act/kT) with E_act > 0.5 eV (median {hot.e_act.median():.2f} eV, '
                                    f'i.e. {dex:.1f} dex suppression at 5000 K relative to a flat rate)',
                                    len(hot), [f"{r.lo}-{r.hi} dE={r.dev:.2f} eV E_act={r.e_act:+.2f} eV" for r in
                                               hot.sort_values('e_act', ascending=False).itertuples()]))
    findings += check_ch_multiplicity(atom, col, vals)
    # coverage of radiative transitions
    if hasattr(atom, 'bbt'):
        pairs = set(zip(atom.bbt.i, atom.bbt.j))
        for keys, name in ((('CE', 'OHM', 'OHMEGA'), 'electron'), (('CH',), 'hydrogen')):
            have = set(zip(col.lo[col.key.isin(keys)], col.hi[col.key.isin(keys)]))
            miss = pairs - have
            sev = 'info' if not miss else ('moderate' if len(miss) > 0.1 * len(pairs) else 'minor')
            if pairs and len(miss) == len(pairs):
                sev = 'serious'
            findings.append(Finding(f'{name} excitation coverage', sev,
                                    f'{len(pairs) - len(miss)} of {len(pairs)} radiative pairs have {"/".join(keys)} '
                                    'data (the rest get none unless use_classic_col)', len(miss)))
    cont = continuum_level(atom)
    ionisable = [i for i in range(len(lvl)) if cont[i] >= 0]
    for keys, name in ((('CI',), 'electron ionisation'), (('CH0', 'CHI'), 'hydrogen ionisation/charge transfer')):
        have = set(col.lo[col.key.isin(keys)])
        miss = [i for i in ionisable if i not in have]
        sev = 'info' if not miss else ('moderate' if name.startswith('electron') else 'minor')
        if ionisable and len(miss) == len(ionisable) and name.startswith('electron'):
            sev = 'serious'
        findings.append(Finding(f'{name} coverage', sev,
                                f'{len(ionisable) - len(miss)} of {len(ionisable)} ionisable levels have '
                                f'{"/".join(keys)} data', len(miss), miss))
    # CE against van Regemorter for allowed transitions
    if hasattr(atom, 'bbt'):
        bbt = line_table(atom)
        fval = {(r.i, r.j): r.f for r in bbt.itertuples()}
        ce = col[(col.key == 'CE')]
        ratios = []
        tmid = np.array([6000.0])
        for k, r in zip(ce.index, ce.itertuples()):
            f = fval.get((r.lo, r.hi))
            if f is None or r.dev <= 0:
                continue
            ct = collision_interpolate(temps, vals[k], tmid, 'CE')
            cdn = gg[r.lo] / gg[r.hi] * ct * np.sqrt(tmid)
            cup = cdn * gg[r.hi] / gg[r.lo] * np.exp(-r.dev * EE / (BK * tmid))
            ratios.append((r.lo, r.hi, float(np.log10(cup / van_regemorter(f, r.dev, tmid))[0])))
        if ratios:
            rr = np.array([x[2] for x in ratios])
            far = [x for x in ratios if abs(x[2]) > 2]
            findings.append(Finding('CE vs van Regemorter', 'moderate' if len(far) > 0.1 * len(rr) else 'info',
                                    f'{len(rr)} allowed transitions: log(CE/vR) at 6000 K median {np.median(rr):+.2f}, '
                                    f'16-84% {np.percentile(rr, 16):+.2f} .. {np.percentile(rr, 84):+.2f}; '
                                    f'{len(far)} beyond 2 dex', len(far),
                                    [f"{a}-{b}: {c:+.1f}" for a, b, c in sorted(far, key=lambda x: -abs(x[2]))]))
    return col, findings


# ---------------------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------------------
def audit(atom, nist_levels=None, nist_limits=None, linelists=None, temp=5800.0, classic=False):
    """Run every check that the given data allow. ``linelists`` maps an ion stage to a
    line-list table for that stage; ``classic`` says the atom was read from the classic
    format. Returns {name: table} and the findings."""
    tables = {}
    findings = []
    tables['levels'], fnd = check_levels(atom)
    findings += fnd
    tables['lines'], fnd = check_lines(atom, classic=classic)
    findings += fnd
    tab, fnd = check_bf(atom, temp=temp)
    tables['bf'] = tab
    findings += fnd
    tab, fnd = check_collisions(atom)
    tables['collisions'] = tab
    findings += fnd
    if nist_limits is not None:
        tab, fnd = check_limits(atom, nist_limits)
        tables['limits'] = tab
        findings += fnd
    if nist_levels is not None:
        tab, fnd = check_completeness(atom, nist_levels)
        tables['partition'] = tab
        findings += fnd
        tab, fnd = check_nist_levels(atom, nist_levels)
        tables['nist_match'] = tab
        findings += fnd
    for stage, llist in (linelists or {}).items():
        tab, fnd = compare_linelist(atom, llist, stage, nist_levels=nist_levels)
        tables[f'linelist_{stage}'] = tab
        findings += fnd
        if nist_levels is not None:
            tab, fnd = check_level_matching(atom, nist_levels, llist, stage)
            tables[f'level_matching_{stage}'] = tab
            findings += fnd
    findings.sort(key=lambda f: SEVERITIES.index(f.severity))
    return tables, findings


def main(argv=None):
    import argparse
    from . import atomic_data
    parser = argparse.ArgumentParser(description='Sanity checks for a model atom')
    parser.add_argument('atom_file', type=Path)
    parser.add_argument('--nist', type=Path, help='folder with NIST ASD level tables, e.g. MgI_levels.csv (atomic_data.read_nist_levels)')
    parser.add_argument('--linelist', type=Path, nargs='+',
                        help='line lists in Turbospectrum format (one file or one per species)')
    parser.add_argument('--temp', type=float, default=5800.0)
    args = parser.parse_args(argv)
    atom = load_atom(args.atom_file)
    nist_levels = nist_limits = None
    stages = sorted(set(atom.lvl.ion))
    if args.nist:
        nist_levels, nist_limits = atomic_data.nist_levels_for_atom(args.nist, atom.id, stages)
    linelists = None
    if args.linelist:
        linelists = {}
        for path in args.linelist:
            linelists.update(atomic_data.linelist_by_stage(path, atom.id, stages[:-1]))
    classic = not is_pandas_format(args.atom_file)
    _, findings = audit(atom, nist_levels, nist_limits, linelists, temp=args.temp, classic=classic)
    if classic:
        findings += check_comments(args.atom_file)[1]
    print(f'{args.atom_file}: {len(atom.lvl)} levels, {len(getattr(atom, "bbt", []))} lines')
    for f in findings:
        print(f)


if __name__ == '__main__':
    main()
