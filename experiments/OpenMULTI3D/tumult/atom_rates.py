"""Per-transition rates of an NLTE run and the reduction criteria built on them.

dispatch.x writes the rate files with ``&output_params save_net_rates=T`` (see
``atom_mod%write_rates``): per patch, ``save/rates_<atom>_<zs>.bin`` with the populations
and the radiative rate coefficients of the last update, and ``save/colrates_<atom>_<zs>.bin``
with the upward collision coefficients. :func:`read_rates` joins the patches of a 1D run.

Statistical equilibrium makes the net rates into and out of each level cancel. A
transition whose net rate is a negligible part of the net flow through both its levels
does not shape the departure coefficients: removing it changes the equilibrium only
through the opacity it adds to other transitions (Lind et al. 2017, Mallinson et al.
2022). :func:`transition_relevance` measures that share per transition.
"""
from pathlib import Path

import numpy as np
import pandas as pd


def _read_patch(path, col_path):
    with open(path, 'rb') as f:
        nk, nrad, nx, ny, nz, it = np.fromfile(f, np.int32, 6)
        irad = np.fromfile(f, np.int32, nrad) - 1
        jrad = np.fromfile(f, np.int32, nrad) - 1
        ncell = nx * ny * nz
        pop = np.fromfile(f, np.float64, ncell * nk).reshape((nk, nz, ny, nx))
        nstar = np.fromfile(f, np.float64, ncell * nk).reshape((nk, nz, ny, nx))
        rij = np.fromfile(f, np.float32, ncell * nrad).reshape((nrad, nz, ny, nx))
        rji = np.fromfile(f, np.float32, ncell * nrad).reshape((nrad, nz, ny, nx))
    with open(col_path, 'rb') as f:
        nkc, cx, cy, cz = np.fromfile(f, np.int32, 4)
        cup = np.fromfile(f, np.float32, cx * cy * cz * nk * (nk - 1) // 2).reshape((-1, cz, cy, cx))
    # one column: (level or transition, z)
    first = (slice(None), slice(None), 0, 0)
    return dict(iter=int(it), irad=irad, jrad=jrad, pop=pop[first], nstar=nstar[first],
                rij=rij[first].astype(float), rji=rji[first].astype(float), cup=cup[first].astype(float))


def read_rates(run_dir, atom_id):
    """Rates of the first column of a run, patches joined in z (index 0 at the bottom).
    Returns a dict with pop, nstar (nk, nz), rij, rji (nrad, nz), cup (npair, nz), irad,
    jrad and the iteration."""
    save = Path(run_dir) / 'save'
    files = sorted(save.glob(f'rates_{atom_id.strip()}_*.bin'))
    if not files:
        raise FileNotFoundError(f'no rates_{atom_id}_*.bin in {save}')
    parts = [_read_patch(f, f.with_name('col' + f.name)) for f in files]
    out = dict(iter=parts[0]['iter'], irad=parts[0]['irad'], jrad=parts[0]['jrad'])
    for key in ('pop', 'nstar', 'rij', 'rji', 'cup'):
        out[key] = np.concatenate([p[key] for p in parts], axis=1)
    return out


def pair_index(i, j):
    """Index of the level pair i < j in the j-major order of the collision file."""
    return j * (j - 1) // 2 + i


def net_rates(rates):
    """Net radiative rate of every transition, n_i R_ij - n_j R_ji [cm^-3 s^-1], and the
    net collisional rate of every level pair, (nrad, nz) and (npair, nz)."""
    pop, nstar = rates['pop'], rates['nstar']
    irad, jrad = rates['irad'], rates['jrad']
    rad = pop[irad] * rates['rij'] - pop[jrad] * rates['rji']
    nk = pop.shape[0]
    jj, ii = np.array([(j, i) for j in range(1, nk) for i in range(j)]).T
    cup = rates['cup']
    cdn = cup * nstar[ii] / np.maximum(nstar[jj], 1e-300)
    col = pop[ii] * cup - pop[jj] * cdn
    return rad, col, ii, jj


def transition_relevance(rates, departure_floor=1e-3):
    """Shares of the net flow that each radiative transition carries, maximised over
    depth.

    For level k the net flow is N_k = 1/2 sum |net| over all radiative and collisional
    transitions of k. ``share`` is the local measure, the larger of |net_t| / N_i and
    |net_t| / N_j: a transition that dominates the balance of a level it touches has a
    share near one, however weakly populated that level is. ``share_atom`` divides by
    the largest N_k at that depth instead, which ranks transitions by their part in the
    flow through the whole atom, so that the many transitions among sparsely populated
    high levels rank low. Depths where no level departs from LTE by more than
    ``departure_floor`` (|ln b| <= floor) are left out, since every net rate there is
    rounding noise. The table also holds the depth where ``share`` peaks, the net
    radiative and collisional rates of the pair there, and whether collisions dominate
    the pair's gross rate."""
    rad, col, ii, jj = net_rates(rates)
    nk, nz = rates['pop'].shape
    irad, jrad = rates['irad'], rates['jrad']
    flow = np.zeros((nk, nz))
    np.add.at(flow, irad, np.abs(rad))
    np.add.at(flow, jrad, np.abs(rad))
    np.add.at(flow, ii, np.abs(col))
    np.add.at(flow, jj, np.abs(col))
    flow *= 0.5
    with np.errstate(divide='ignore', invalid='ignore'):
        lnb = np.abs(np.log(rates['pop'] / rates['nstar']))
        active = np.nanmax(lnb, axis=0) > departure_floor
        share = np.maximum(np.abs(rad) / flow[irad], np.abs(rad) / flow[jrad])
        share_atom = np.abs(rad) / flow.max(0)
    share = np.where(active[None, :] & np.isfinite(share), share, 0)
    share_atom = np.where(active[None, :] & np.isfinite(share_atom), share_atom, 0)
    kz = share.argmax(1)
    kk = np.arange(len(irad))
    pidx = pair_index(irad, jrad)
    colnet = col[pidx, kz]
    gross_rad = rates['pop'][irad, kz] * rates['rij'][kk, kz] + rates['pop'][jrad, kz] * rates['rji'][kk, kz]
    gross_col = rates['pop'][irad, kz] * rates['cup'][pidx, kz] * 2
    return pd.DataFrame(dict(i=irad, j=jrad, share=share[kk, kz], share_atom=share_atom.max(1),
                             z=kz, net_rad=rad[kk, kz], net_col=colnet,
                             collisional=gross_col > gross_rad))


def departures(rates):
    """Departure coefficients b = n/n* (nk, nz)."""
    return rates['pop'] / rates['nstar']



def diagnostic_relevance(rates, sources, departure_floor=1e-3):
    """Influence of every level and radiative transition on the departure coefficients
    of the diagnostic levels, maximised over depth.

    ``sources`` (nk, nz) marks, per depth, the levels whose b the diagnostic lines take
    there (both ends of each line, over the depths where it forms). Statistical
    equilibrium splits the net flow F_k = 1/2 sum |net| through level k over its level
    pairs; c(a -> b) = (|net_rad(a,b)| + |net_col(a,b)|) / F_a is the part of a's flow
    that passes through b. A level's influence I_k is the largest product of c along a
    path from a source level to k (1 for the sources), found as a shortest path in -ln c:
    a level whose flow reaches the sources only through branches that each carry a small
    part of the flow there cannot move their b much. A radiative transition (i, j)
    carries I_i |net| / F_i through its lower and I_j |net| / F_j through its upper level,
    its share of the flow of that level weighted by the level's influence; for a
    transition between two sources the larger of the two is the local share of
    :func:`transition_relevance`. The two are kept apart so that a pruning can add up
    what it removes from each level. Depths without sources, or where no level departs
    from LTE by more than ``departure_floor`` (|ln b|), do not count.

    Returns (influence (nk,), relevance through the lower and the upper level, each
    maximised over depth (nrad, 2))."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    rad, col, ii, jj = net_rates(rates)
    nk, nz = rates['pop'].shape
    # a continuum without a final level (the top stage's ground) has jrad = -1 and no rates
    valid = rates['jrad'] >= 0
    irad, jrad = rates['irad'][valid], rates['jrad'][valid]
    rad = rad[valid]
    with np.errstate(divide='ignore', invalid='ignore'):
        lnb = np.abs(np.log(rates['pop'] / rates['nstar']))
    a = np.concatenate([irad, ii])
    b = np.concatenate([jrad, jj])
    influence = np.zeros(nk)
    relevance = np.zeros((len(irad), 2))
    for z in range(nz):
        src = np.flatnonzero(sources[:, z])
        if len(src) == 0 or not np.nanmax(lnb[:, z]) > departure_floor:
            continue
        w = np.concatenate([np.abs(rad[:, z]), np.abs(col[:, z])])
        pair = coo_matrix((w, (a, b)), shape=(nk, nk)).tocsr()
        pair = (pair + pair.T).tocoo()
        flow = np.zeros(nk)
        np.add.at(flow, pair.row, 0.5 * pair.data)
        frac = pair.data / np.maximum(flow[pair.row], 1e-300)
        ok = frac > 0
        # -ln(frac) >= 0 since |net(a,b)| <= 2 F_a; the offset keeps zero-cost edges,
        # which a sparse graph would drop, in the graph
        cost = np.maximum(-np.log(np.minimum(frac[ok], 1.0)), 0) + 1e-12
        graph = coo_matrix((cost, (pair.row[ok], pair.col[ok])), shape=(nk, nk)).tocsr()
        infl = np.exp(-dijkstra(graph, directed=True, indices=src, min_only=True))
        infl[src] = 1.0
        influence = np.maximum(influence, infl)
        with np.errstate(divide='ignore', invalid='ignore'):
            net = np.abs(rad[:, z])
            rel = np.stack([np.nan_to_num(infl[irad] * net / flow[irad]),
                            np.nan_to_num(infl[jrad] * net / flow[jrad])], axis=1)
        relevance = np.maximum(relevance, rel)
    out = np.zeros((len(valid), 2))
    out[valid] = relevance
    return influence, out


def level_index(keys, stage, ev):
    """Atom levels that line-list levels take their departure coefficients from: the
    nearest of the departure ``keys`` (100 * stage + energy above the stage ground [eV],
    :func:`tumult.atom_audit.level_keys`) to 100 * stage + ev, as ``find_closest`` in
    HOP/linelist_mod.f90 matches them."""
    target = 100 * np.asarray(stage, float) + np.asarray(ev, float)
    return np.abs(keys[None, :] - target[:, None]).argmin(1)


def removable(ends, lower, upper, nk, budget):
    """Transitions that can go together: in order of rising relevance, each transition is
    removed while the relevance removed through each of its levels (``ends`` (n, 2)
    through the lower and the upper level, as :func:`diagnostic_relevance` gives them)
    adds up to at most ``budget`` per level. Summing bounds the joint effect of many small
    removals on one level, which a cut on each transition's relevance alone does not.
    Returns a mask of the removable transitions."""
    removed = np.zeros(nk)
    out = np.zeros(len(ends), bool)
    for t in np.argsort(ends.max(1), kind='stable'):
        i, j = lower[t], upper[t]
        if removed[i] + ends[t, 0] <= budget and removed[j] + ends[t, 1] <= budget:
            removed[i] += ends[t, 0]
            removed[j] += ends[t, 1]
            out[t] = True
    return out
