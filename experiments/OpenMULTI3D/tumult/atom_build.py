"""Tools that change or assemble model atoms in the pandas format.

Bound-free tables are written for dispatch's interpolation: :func:`resample_for_hermite`
picks the nodes from which ``trans_mod%hermite_interp`` (monotone cubic Hermite in
frequency, zero outside the table) reproduces a given cross section to a tolerance.
Cross sections with resonances are first smoothed (:func:`smooth_cross_section`), as
Lind & Amarsi (2024) and FORMATO3 do, because the positions of computed resonances are
uncertain and a table cannot carry them at an affordable number of frequencies.
"""
import numpy as np
import pandas as pd

from .fast_atom import PandasAtom
from .atom_audit import (CC, bf_interpolate, bf_tables, compare_bf_reference, continuum_level, ion_ground,
                         term_ranks)


def smooth_cross_section(lam, sigma, width=0.03):
    """Gaussian smoothing of a cross section in ln(frequency) with standard deviation
    ``width``, on the original points. The edge (largest wavelength) stays the edge:
    points near it are averaged only over the part of the kernel inside the table."""
    lam = np.asarray(lam, float)
    order = np.argsort(lam)[::-1]                  # increasing frequency
    xx = np.log(CC * 1e8 / lam[order])
    yy = np.asarray(sigma, float)[order]
    # trapezoid weights of the nodes, so that unevenly sampled resonances count by extent
    dx = np.gradient(xx)
    out = np.empty_like(yy)
    for k, x0 in enumerate(xx):
        sel = np.abs(xx - x0) < 4 * width
        ww = np.exp(-0.5 * ((xx[sel] - x0) / width)**2) * dx[sel]
        out[k] = np.sum(ww * yy[sel]) / np.sum(ww)
    result = np.empty_like(out)
    result[order] = out
    return result


def resample_for_hermite(lam, sigma, rtol=0.03, floor=1e-3, max_points=200):
    """Nodes for a bound-free table: the subset of (lam, sigma) from which dispatch's
    interpolation reproduces sigma at every input point within ``rtol`` relative to
    max(sigma, floor * max(sigma)). Starts from the two ends and inserts the worst point
    until the tolerance or ``max_points`` is reached. Returns lam increasing."""
    lam = np.asarray(lam, float)
    order = np.argsort(lam)
    lam = lam[order]
    sigma = np.asarray(sigma, float)[order]
    scale = np.maximum(sigma, floor * sigma.max())
    keep = {0, len(lam) - 1}
    while len(keep) < min(max_points, len(lam)):
        idx = np.array(sorted(keep))
        pred = bf_interpolate(lam[idx], sigma[idx], lam)
        err = np.abs(pred - sigma) / scale
        err[idx] = 0
        worst = int(np.argmax(err))
        if err[worst] <= rtol:
            break
        keep.add(worst)
    idx = np.array(sorted(keep))
    return lam[idx], sigma[idx]


def topbase_continua(atom, topbase, stage, width=0.03, rtol=0.03, threshold_tolerance=0.25):
    """Bound-free tables for the levels of one ion stage from TOPbase cross sections.

    Levels are matched by LS symmetry and rank (:func:`tumult.atom_audit.term_ranks`),
    the TOPbase threshold is moved to the atom's edge, the cross section smoothed and
    resampled for dispatch's interpolation. Returns {level: (lam, sigma)} for matched
    levels and a table that compares each new and old table with the smoothed and with
    the unsmoothed reference."""
    from .atomic_data import MB, RYD_AA, RYD_EV
    lvl = atom.lvl
    cont = continuum_level(atom)
    old = bf_tables(atom)
    ranks = term_ranks(atom, stage)
    new = {}
    rows = []
    for i, (term, rank) in ranks.items():
        if cont[i] < 0 or (*term, rank) not in topbase:
            continue
        energy, sigma = topbase[(*term, rank)]
        bind_ry = (lvl.ev[cont[i]] - lvl.ev[i]) / RYD_EV
        off = energy[0] / bind_ry - 1
        if abs(off) > threshold_tolerance:
            continue
        lam_ref = RYD_AA / (energy - energy[0] + bind_ry)
        sig_ref = sigma * MB
        smooth = smooth_cross_section(lam_ref, sig_ref, width)
        lam_new, sig_new = resample_for_hermite(lam_ref, smooth, rtol)
        new[i] = (lam_new, sig_new)
        row = dict(i=i, label=lvl.label[i].strip(), threshold_offset=off, npts_old=len(old.get(i, ((), ()))[0]),
                   npts_new=len(lam_new), npts_ref=len(lam_ref))
        for name, table in (('new', (lam_new, sig_new)), ('old', old.get(i))):
            if table is None:
                continue
            cmp = compare_bf_reference(*table, lam_ref, sig_ref, temps=(6000.,))
            row[f'{name}_vs_raw'] = cmp['rate_6000']
        row['smooth_vs_raw'] = compare_bf_reference(lam_ref, smooth, lam_ref, sig_ref, temps=(6000.,))['rate_6000']
        rows.append(row)
    return new, pd.DataFrame(rows)


def replace_continua(atom, tables):
    """A copy of the atom's bound-free table with the continua of the given levels
    replaced ({level: (lam, sigma)})."""
    bft = atom.bft[~atom.bft.i.isin(list(tables))]
    add = [pd.DataFrame(dict(i=i, lam_air=lam, alpha=sig)) for i, (lam, sig) in tables.items()]
    return pd.concat([bft] + add, ignore_index=True)


def drop_duplicate_collisions(col):
    """Collision records with each (key, lower, upper) once, keeping the first: gencol
    adds repeated records."""
    lo = np.minimum(col.i, col.j)
    hi = np.maximum(col.i, col.j)
    return col[~pd.DataFrame(dict(k=col.key, lo=lo, hi=hi)).duplicated()].reset_index(drop=True)



def prune(atom, keep_lines, keep_levels):
    """The atom with only the bound-bound transitions in ``keep_lines`` (mask over
    atom.bbt) and the levels in ``keep_levels`` (mask over atom.lvl). Every stage's ground
    level is kept, since the continua and the collisional ionisation end there. A removed
    level takes its lines, bound-free table and collision records with it; the remaining
    levels are renumbered."""
    keep_levels = np.asarray(keep_levels, bool).copy()
    keep_levels[list(ion_ground(atom).values())] = True
    gone = np.flatnonzero(~keep_levels)
    out = PandasAtom.from_pandas(atom.id, atom.abnd, atom.ai, atom.lvl.copy(),
                                 bbt=atom.bbt[np.asarray(keep_lines, bool)
                                              & keep_levels[atom.bbt.i] & keep_levels[atom.bbt.j]].copy(),
                                 bft=atom.bft[keep_levels[atom.bft.i]].copy(),
                                 col=atom.col[keep_levels[atom.col.i] & keep_levels[atom.col.j]].copy())
    out.remove_levels(gone)
    out.nk = len(out.lvl)
    out.nline = len(out.bbt)
    out.ncol = len(out.col)
    return out


# ---------------------------------------------------------------------------------------
# superlevels
# ---------------------------------------------------------------------------------------
OHM_C0 = 8.629e-6      # gencol OHM constant: q_down = C0 / g_u * Omega / sqrt(T)
H_MINUS_EV = 0.754


def fine_structure_groups(atom, max_spread_ev=0.5):
    """Groups of the levels of one LS term: same ion and label up to the J value
    (e.g. 'a5D4', 'a5D3', ... -> 'a5D'), each group's levels within ``max_spread_ev``."""
    import re
    lvl = atom.lvl
    keys = [re.sub(r'([SPDFGHIKLMNOQ])\d+(/\d)?(\*?)\s*$', r'\1\3', str(l).split('=')[-1].strip())
            for l in lvl.label]
    groups = {}
    for i, (k, s, e) in enumerate(zip(keys, lvl.ion, lvl.ev)):
        groups.setdefault((s, k), []).append(i)
    out = []
    for members in groups.values():
        members = sorted(members, key=lambda m: lvl.ev[m])
        start = 0
        for k in range(1, len(members) + 1):
            if k == len(members) or lvl.ev[members[k]] - lvl.ev[members[start]] > max_spread_ev:
                out.append(members[start:k])
                start = k
    return sorted(out, key=lambda g: lvl.ev[g[0]])


def energy_bin_groups(atom, groups, binding_ev, width_ev, by_parity=True):
    """Merge the given groups further: groups of an ion bound by less than
    ``binding_ev`` are collected into energy bins of ``width_ev`` (separately for each
    parity when ``by_parity``)."""
    from .atom_audit import continuum_level
    lvl = atom.lvl
    cont = continuum_level(atom)
    out = []
    bins = {}
    for grp in groups:
        i = grp[0]
        if cont[i] < 0:
            out.append(grp)
            continue
        ev = np.average(lvl.ev[grp], weights=lvl.g[grp])
        bind = lvl.ev[cont[i]] - ev
        if bind >= binding_ev:
            out.append(grp)
            continue
        parity = int('*' in str(lvl.label[i])) if by_parity else 0
        key = (int(lvl.ion[i]), parity, int(np.floor(bind / width_ev)))
        bins.setdefault(key, []).extend(grp)
    out += [sorted(g) for g in bins.values()]
    return sorted(out, key=lambda g: lvl.ev[g[0]])


def _boltzmann(gg, ev, temp):
    ww = gg * np.exp(-(ev - ev.min()) * 11604.518 / temp)
    return ww / ww.sum()


def merge_levels(atom, groups, temp_lines=5500.0, rtol_bf=0.03):
    """A PandasAtom whose levels are the given groups (lists of level indices covering
    every level once).

    Following Lind & Amarsi (2024), a superlevel has g = sum g_i. Its energy keeps the
    partition function at ``temp_lines``, g_I exp(-E_I/kT) = sum g_i exp(-E_i/kT); the
    g-weighted mean energy of Lind & Amarsi loses 6% of U(Fe II) at 4000 K when the
    ground term a6D is merged. Lines between two superlevels merge into one, at the Ritz wavelength of the
    superlevels, with f_IJ = sum_i w_i sum_j f_ij and Boltzmann weights w_i at
    ``temp_lines``; it keeps the broadening, nq, qmax and q0 of the strongest
    component. Lines inside a superlevel are dropped. Continua merge with the same
    weights, each component shifted to the superlevel's edge, and are resampled for
    dispatch's interpolation. Collision records turn into upward rate coefficients per
    collider at each tabulated temperature, are summed over the upper and Boltzmann-
    weighted over the lower components at that temperature, and turn back into a CE, CH,
    CI, CH0 or CHI record (OHM becomes CE), so the merged upward rate is exact when the
    components are in relative LTE.
    """
    from .fast_atom import PandasAtom
    from .physics import vac2air
    from .atom_audit import HC_AA_EV, bf_interpolate, bf_tables
    lvl = atom.lvl.reset_index(drop=True)
    nk = len(lvl)
    owner = np.full(nk, -1)
    for k, grp in enumerate(groups):
        owner[grp] = k
    if (owner < 0).any():
        raise ValueError(f'levels {np.where(owner < 0)[0]} are in no group')
    gg = lvl.g.to_numpy(float)
    ev = lvl.ev.to_numpy(float)
    new_lvl = []
    for grp in groups:
        g_tot = gg[grp].sum()
        label = lvl.label[grp[0]] if len(grp) == 1 else f'{str(lvl.label[grp[0]]).strip()}+{len(grp) - 1}'
        # the energy that keeps the superlevel's Boltzmann factor at temp_lines
        boltz = np.sum(gg[grp] * np.exp(-(ev[grp] - ev[grp].min()) * 11604.518 / temp_lines)) / g_tot
        e_sup = ev[grp].min() - np.log(boltz) * temp_lines / 11604.518
        new_lvl.append(dict(ev=e_sup, g=g_tot, ion=int(lvl.ion[grp[0]]), label=label))
    new_lvl = pd.DataFrame(new_lvl)
    weights_a = np.zeros(nk)
    for grp in groups:
        weights_a[grp] = _boltzmann(gg[grp], ev[grp], temp_lines)

    # lines
    bbt = atom.bbt.copy()
    bbt['I'] = owner[bbt.i]
    bbt['J'] = owner[bbt.j]
    bbt = bbt[bbt.I != bbt.J]
    bbt['wf'] = weights_a[bbt.i] * 10**bbt.loggf / gg[bbt.i]
    bbt['gf'] = 10**bbt.loggf
    rows = []
    for (I, J), comp in bbt.groupby(['I', 'J']):
        lo, up = (I, J) if new_lvl.ev[I] < new_lvl.ev[J] else (J, I)
        strongest = comp.loc[comp.gf.idxmax()]
        f_ij = comp.wf.sum()
        lam_vac = HC_AA_EV / (new_lvl.ev[up] - new_lvl.ev[lo])
        row = strongest.copy()
        row['i'], row['j'] = lo, up
        # read_pandas applies air2vac at every wavelength
        row['lam0'] = float(vac2air(lam_vac))
        row['evi'] = new_lvl.ev[lo]
        row['loggf'] = np.log10(f_ij * new_lvl.g[lo]) if lo == I else np.log10(f_ij * new_lvl.g[I])
        rows.append(row)
    new_bbt = pd.DataFrame(rows)[list(atom.bbt.columns)].reset_index(drop=True)
    new_bbt[['i', 'j', 'nq']] = new_bbt[['i', 'j', 'nq']].astype(int)

    # continua
    tables = bf_tables(atom)
    new_bf = []
    cont_ev = {}
    for i in range(nk):
        # continuum of each level: first level of the next ion
        nxt = np.where(lvl.ion.to_numpy() == lvl.ion[i] + 1)[0]
        cont_ev[i] = ev[nxt[0]] if len(nxt) else np.nan
    for k, grp in enumerate(groups):
        comps = [i for i in grp if i in tables]
        if not comps:
            continue
        edge_I = HC_AA_EV / (cont_ev[grp[0]] - new_lvl.ev[k])
        shifted = []
        for i in comps:
            lam, alpha = tables[i]
            # photon energy above the component's own edge, put on the superlevel's
            de = HC_AA_EV / lam - (cont_ev[i] - ev[i])
            shifted.append((HC_AA_EV / (de + HC_AA_EV / edge_I), alpha))
        grid = np.unique(np.concatenate([s[0] for s in shifted]))
        wsum = weights_a[comps].sum()
        sigma = sum(weights_a[i] / wsum * bf_interpolate(lam_s, alpha, grid) for i, (lam_s, alpha) in zip(comps, shifted))
        good = sigma > 0
        if good.sum() < 2:
            continue
        lam_n, sig_n = resample_for_hermite(grid[good], sigma[good], rtol_bf)
        new_bf.append(pd.DataFrame(dict(i=k, lam_air=lam_n, alpha=sig_n)))
    new_bft = pd.concat(new_bf, ignore_index=True)

    # collisions
    col = atom.col
    temps = np.array([float(t) for t in col.columns[3:]])
    vals = col.iloc[:, 3:].to_numpy(float)
    kt = 11604.518 / temps                      # 1/kT [1/eV]
    sqt = np.sqrt(temps)
    wt = np.zeros((nk, len(temps)))
    for grp in groups:
        for m, t in enumerate(temps):
            wt[grp, m] = _boltzmann(gg[grp], ev[grp], t)
    acc = {}
    for r, v in zip(col.itertuples(index=False), vals):
        key, i, j = r[0], int(r[1]), int(r[2])
        lo, up = min(i, j), max(i, j)
        if owner[lo] == owner[up]:
            continue
        de = ev[up] - ev[lo]
        if key in ('OHM', 'OHMEGA'):
            cup = OHM_C0 * v / (gg[lo] * sqt) * np.exp(-de * kt)
            kind = 'CE'
        elif key == 'CE':
            cup = v * sqt * np.exp(-de * kt)
            kind = 'CE'
        elif key == 'CH':
            cup = v * gg[up] / gg[lo] * np.exp(-de * kt)
            kind = 'CH'
        elif key == 'CI':
            cup = v * sqt * np.exp(-de * kt)
            kind = 'CI'
        elif key in ('CH0', 'CHI'):
            cup = v
            kind = key
        else:
            raise ValueError(f'collision key {key} cannot be merged')
        index = (kind, owner[lo], owner[up])
        acc[index] = acc.get(index, 0) + wt[lo] * cup
    rows = []
    for (kind, lo, up), cup in acc.items():
        de = new_lvl.ev[up] - new_lvl.ev[lo]
        if kind == 'CE':
            ct = cup * np.exp(de * kt) / sqt
        elif kind == 'CH':
            ct = cup * new_lvl.g[lo] / new_lvl.g[up] * np.exp(de * kt)
        elif kind == 'CI':
            ct = cup * np.exp(de * kt) / sqt
        else:
            ct = cup
        rows.append([kind, lo, up] + list(ct))
    new_col = pd.DataFrame(rows, columns=list(col.columns))
    new_col = new_col[(new_col.iloc[:, 3:] > 0).all(axis=1)].reset_index(drop=True)
    return PandasAtom.from_pandas(atom.id, atom.abnd, atom.ai, new_lvl, new_bbt, new_bft, new_col)


def ce_from_upsilon(upsilon, g_lower, temps):
    """gencol CE value for an effective collision strength: the upward rate per electron
    8.629e-6 Y / (g_l sqrt(T)) exp(-dE/kT) is ct sqrt(T) exp(-dE/kT) in gencol's form."""
    return OHM_C0 * np.asarray(upsilon, float) / (g_lower * np.asarray(temps, float))


def ce_van_regemorter(f, de, temps):
    """gencol CE value of the van Regemorter (1962) rate of an allowed transition."""
    from .atom_audit import van_regemorter
    temps = np.asarray(temps, float)
    return van_regemorter(f, de, temps) * np.exp(de * 11604.518 / temps) / np.sqrt(temps)


def term_members(atom):
    """{(ion, term label without J): [levels]} from the level labels."""
    import re
    out = {}
    for i, l in enumerate(atom.lvl.label):
        key = re.sub(r'([SPDFGHIKLMNOQ])\d+(/\d)?(\*?)\s*$', r'\1\3', str(l).split('=')[-1].strip())
        out.setdefault((int(atom.lvl.ion[i]), key), []).append(i)
    return out


# ---------------------------------------------------------------------------------------
# assembling an atom
# ---------------------------------------------------------------------------------------
def match_lines(atom, lines, stage, tolerance_ev=0.005):
    """Bound-bound table rows for the lines of one ion stage from a line list
    (:func:`tumult.atomic_data.linelist_by_stage`: lam0 in TS convention, evi above the
    stage's ground [eV], loggf, gvdw, grad, gstk).

    The lower level is the level of that stage nearest to evi and the upper level the
    one nearest to evi + h nu, both within ``tolerance_ev``. The list's upper_g column
    is 1 throughout the Turbospectrum lists and is not used. Components between the same levels
    (isotopes, hyperfine structure) add their gf. The wavelength is the Ritz one of the
    two levels, converted to air at every wavelength because read_pandas converts
    lam_air back with air2vac throughout; that conversion fails near its poles at 877 and
    1603 AA (:func:`tumult.atom_audit.check_lines`). Returns the rows and the lines left
    unmatched."""
    from .atom_audit import HC_AA_EV
    from .atomic_data import ts_vacuum_wavelength
    from .physics import vac2air
    lvl = atom.lvl
    sel = np.where(lvl.ion.to_numpy() == stage)[0]
    ev0 = lvl.ev[sel[0]]
    ev_rel = lvl.ev.to_numpy(float)[sel] - ev0
    lam_vac = ts_vacuum_wavelength(lines.lam0)
    ev_up = lines.evi.to_numpy(float) + HC_AA_EV / lam_vac
    lo = np.abs(ev_rel[None, :] - lines.evi.to_numpy(float)[:, None]).argmin(1)
    dup = np.abs(ev_rel[None, :] - ev_up[:, None])
    up = dup.argmin(1)
    good = ((np.abs(ev_rel[lo] - lines.evi) < tolerance_ev)
            & (dup[np.arange(len(up)), up] < tolerance_ev) & (up != lo))
    tab = lines[good].copy()
    tab['i'], tab['j'] = sel[lo[good]], sel[up[good]]
    tab['gf'] = 10**tab.loggf
    rows = []
    for (i, j), comp in tab.groupby(['i', 'j']):
        strongest = comp.loc[comp.gf.idxmax()]
        lam = HC_AA_EV / (lvl.ev[j] - lvl.ev[i])
        rows.append(dict(i=i, j=j, lam0=float(vac2air(lam)), evi=lvl.ev[i],
                         loggf=np.log10(comp.gf.sum()), grad=strongest.grad, gvdw=strongest.gvdw,
                         gstk=strongest.gstk))
    return pd.DataFrame(rows), lines[~good]


def seaton_ci(atom, temps, gaunt=(0.1, 0.2, 0.3)):
    """gencol CI records for every level with a continuum from the approximation of
    Seaton (1962) as used by Allen (1973): the ionisation rate per electron is
    1.55e13 g sigma_0 exp(-u) / (u sqrt(T)) cm^3/s with u = E_bind / kT, sigma_0 the
    threshold photoionisation cross section [cm^2] of the atom's own table and g the
    mean Gaunt factor of the stage (``gaunt``: neutral, singly, doubly ionised)."""
    from .atom_audit import bf_tables
    temps = np.asarray(temps, float)
    cont = continuum_level(atom)
    tables = bf_tables(atom)
    lvl = atom.lvl
    rows = []
    for i, (lam, alpha) in tables.items():
        if cont[i] < 0:
            continue
        sigma0 = alpha[np.argmax(lam)]
        u = (lvl.ev[cont[i]] - lvl.ev[i]) * 11604.518 / temps
        # cup = ct sqrt(T) exp(-dE/kT) in gencol's form
        ct = 1.55e13 * gaunt[int(lvl.ion[i]) - 1] * sigma0 / (u * temps)
        rows.append(['CI', i, cont[i]] + list(ct))
    return pd.DataFrame(rows, columns=['key', 'i', 'j'] + [f'{t:.1f}' for t in temps])
