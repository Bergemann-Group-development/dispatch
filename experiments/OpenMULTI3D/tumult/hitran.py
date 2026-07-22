"""Assemble master-format molecular line lists from HITRAN/HITEMP data.

HITRAN distributes no radiative damping constant in the form the master line
lists use.  What it does carry is the Einstein A coefficient of every tabulated
transition, and a level identification (the global/local quanta strings) that is
shared between the lines that start and end on that level.  Summing A over all
transitions out of a level therefore gives the inverse radiative lifetime of
that level, restricted to the transitions the database happens to contain:

    gamma_rad(u->l) = sum_k A_uk + sum_k A_lk                    [1/s]

Both sums are lower bounds -- transitions outside the database's wavelength
coverage or below its intensity cutoff are missing, and for a level that only
ever appears as a lower state the sum is empty.

HOW MUCH THIS MATTERS: for the infrared rovibrational lines these databases
contain, essentially not at all.  At 3.3 um and 4000 K the Doppler width is
~0.22 AA, while gamma_rad ~ 15 1/s contributes a Lorentz width of 4e-10 AA --
nine orders of magnitude smaller -- and van der Waals broadening at
photospheric densities (~1e9 1/s) is another seven orders above that.  Even a
gamma wrong by 1e4 could not move the profile.  The classical fallback that
ga = 0 selects, f_conv * gf * nu^2, is likewise tiny for a weak line (it scales
with gf: ~0.06 1/s at gf = 1e-8), so it is not a hazard either.

The reason to do this properly anyway is electronic transitions, where
gamma_rad ~ 1e8 1/s genuinely shapes the wings.  For a HITEMP rovibrational
list, treat the ga column as cosmetic and spend the effort on --strength-cut,
which controls the line count and therefore the actual cost.

Usage
-----
From a HITEMP archive as downloaded from hitran.org -- the .par.bz2 is read
directly, so hitemp_bz2_extract.py is not needed:

    python -m tumult.hitran --par HITEMP/06_HITEMP2020.par.bz2 --id CH4 \\
        --iso 1 --strength-cut -10 --out input/linelists/hitemp/CH4.txt

Straight from the HITRAN API via hapi2 (needs HAPI2_API_KEY in the
environment, or --api-key).  Note this gets HITRAN, not HITEMP, which is
room-temperature-complete only -- for CO2/CH4/NO/H2O at photospheric
temperatures use the HITEMP archives instead:

    python -m tumult.hitran --fetch NH3 --id NH3 --out .../NH3.txt

Two builders share the same physics:

  build_linelist()            one DataFrame in memory; fine to a few million
                              lines, and the reference implementation.
  build_linelist_streaming()  two passes, memory-bounded, for the 1e7-1e8 line
                              HITEMP files.  Selected automatically by file
                              size; --stream/--no-stream overrides.

Both agree bit-for-bit on lam0/evi/loggf (ga differs at ~1e-15 relative, from
summation order in the level accumulation).

--iso selects one HITRAN isotopologue and the result is written with iso tag
'0', which linelist_mod reads as iso_frac = 1, i.e. the full species density.
Minor isotopologues therefore must NOT simply be merged into such a file: they
would each be given the full number density of the parent species.
"""

import argparse
import bz2
import os
import sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# physical constants (cgs) and unit conversions
# ---------------------------------------------------------------------------

# g_lo * f_lu = GF_CONST * lam_vac[AA]^2 * g_up * A[1/s]   (= m_e c / 8 pi^2 e^2)
GF_CONST = 1.499193e-16

# E[eV] = CM1_TO_EV * E[cm^-1]
CM1_TO_EV = 1.2398419843320026e-4

# linelist_mod.f90 reads the ga column as: 0 -> classical estimate f_conv*gf*nu^2,
# nonzero -> the radiative damping rate in 1/s, used verbatim.  So write gamma
# straight through.  (An earlier revision of the reader rescaled 0 < ga < 100 by
# 1e8; that branch has been removed, since no real line list used it and it
# clashed with molecular gamma values that legitimately fall in that range.)
GA_ENCODE = 1.0


def airvac_corr(w):
    """Refractive index of air, identical to HOP/handyfuncs_mod.f90."""
    s2 = (1e4 / w) ** 2
    return (1.00008336624212083
            + 0.02408926869968 / (130.1065924522 - s2)
            + 0.0001599740894897 / (38.92568793293 - s2))


def vac2air(vac):
    """Vacuum -> air wavelength, the exact inverse of the reader's air2vac.

    Applied at every wavelength, with no 2000 AA cutoff, because
    linelist_mod.f90 calls air2vac unconditionally on the lam_air column.
    Matching it is what makes the round trip exact; the convention question
    is moot since no molecular band of interest sits in the vacuum UV.
    """
    return vac / airvac_corr(vac)


# ---------------------------------------------------------------------------
# input: native HITRAN/HITEMP .par records
# ---------------------------------------------------------------------------

# (name, start, stop, dtype) over the 160-character HITRAN2004 record
PAR_FIELDS = [
    ('molec_id',            0,   2, int),
    ('local_iso_id',        2,   3, str),   # '0' means isotopologue 10
    ('nu',                  3,  15, float),
    ('sw',                 15,  25, float),
    ('a',                  25,  35, float),
    ('gamma_air',          35,  40, float),
    ('gamma_self',         40,  45, float),
    ('elower',             45,  55, float),
    ('n_air',              55,  59, float),
    ('delta_air',          59,  67, float),
    ('global_upper_quanta', 67,  82, str),
    ('global_lower_quanta', 82,  97, str),
    ('local_upper_quanta',  97, 112, str),
    ('local_lower_quanta', 112, 127, str),
    ('gp',                146, 153, float),
    ('gpp',               153, 160, float),
]


"""Numeric-only field subset: everything the streaming path needs.  Skipping
the four quanta strings is most of the parsing cost and nearly all of the
memory, and the streaming path keys levels by energy anyway."""
PAR_NUMERIC = ['local_iso_id', 'nu', 'a', 'elower', 'gp', 'gpp']


def iter_par(filename, numin=None, numax=None, iso=None, fields=None,
             chunksize=32_000_000):
    """Stream a .par or .par.bz2 file as DataFrames, one chunk at a time.

    Records are nu-sorted, so numax permits an early break.  Two full passes
    over a .bz2 cost two decompressions, which is still cheaper than writing
    and re-reading a multi-GB uncompressed intermediate.
    """
    wanted = [f for f in PAR_FIELDS if fields is None or f[0] in fields]
    opener = bz2.open if filename.endswith('.bz2') else open
    with opener(filename, 'rt') as fh:
        while True:
            block = fh.readlines(chunksize)
            if not block:
                return
            block = [r for r in block if len(r.rstrip('\r\n')) >= 160]
            if not block:
                continue
            nu = np.array([float(r[3:15]) for r in block])
            stop = numax is not None and nu[-1] > numax
            keep = np.ones(len(block), dtype=bool)
            if numin is not None:
                keep &= nu >= numin
            if numax is not None:
                keep &= nu <= numax
            if keep.any():
                block = [r for r, k in zip(block, keep) if k]
                cols = {}
                for name, i0, i1, dtype in wanted:
                    raw = [r[i0:i1] for r in block]
                    cols[name] = (raw if dtype is str else
                                  pd.to_numeric(pd.Series(raw).str.strip(),
                                                errors='coerce'))
                df = pd.DataFrame(cols)
                if 'local_iso_id' in df:
                    df['local_iso_id'] = _decode_iso(df['local_iso_id'])
                    if iso is not None:
                        df = df[df['local_iso_id'] == iso]
                df = df.dropna(subset=[c for c in ('nu', 'a', 'elower', 'gp')
                                       if c in df])
                if len(df):
                    yield df.reset_index(drop=True)
            if stop:
                return


def _decode_iso(col):
    """HITRAN writes isotopologue 10 as '0' and 11+ as 'A', 'B', ..."""
    return [10 if c == '0' else (int(c) if c.isdigit() else ord(c) - ord('A') + 11)
            for c in col.astype(str).str.strip()]


def read_par(filename, numin=None, numax=None, chunksize=32_000_000):
    """Parse a HITRAN/HITEMP .par file into a DataFrame.

    Accepts the .par.bz2 archives from hitran.org directly, so there is no need
    to run hitemp_bz2_extract.py first -- that helper writes an uncompressed
    intermediate, which for CH4 or CO2 is several GB of scratch we never need.
    numin/numax do the same wavenumber windowing it offers, including the early
    break, which is valid because HITEMP records are sorted by wavenumber.

    Read in chunks so that the multi-million-line files (CO2 is ~11e6 records,
    CH4 ~31e6) do not exist as Python strings all at once.
    """
    opener = bz2.open if filename.endswith('.bz2') else open
    frames = []
    with opener(filename, 'rt') as fh:
        done = False
        while not done:
            block = fh.readlines(chunksize)
            if not block:
                break
            block = [r for r in block if len(r.rstrip('\r\n')) >= 160]
            if numin is not None or numax is not None:
                nu = np.array([float(r[3:15]) for r in block])
                keep = np.ones(len(block), dtype=bool)
                if numin is not None:
                    keep &= nu >= numin
                if numax is not None:
                    keep &= nu <= numax
                    # records are nu-sorted, so once past numax we are finished
                    done = bool(len(nu)) and nu[-1] > numax
                block = [r for r, k in zip(block, keep) if k]
            if not block:
                continue
            cols = {}
            for name, i0, i1, dtype in PAR_FIELDS:
                raw = [r[i0:i1] for r in block]
                if dtype is str:
                    cols[name] = raw
                else:
                    cols[name] = pd.to_numeric(
                        pd.Series(raw).str.strip(), errors='coerce')
            frames.append(pd.DataFrame(cols))
    df = pd.concat(frames, ignore_index=True)
    # HITRAN writes isotopologue 10 as '0' and 11+ as 'A', 'B', ...
    iso = df['local_iso_id'].str.strip()
    df['local_iso_id'] = [
        10 if c == '0' else (int(c) if c.isdigit() else ord(c) - ord('A') + 11)
        for c in iso
    ]
    return df.dropna(subset=['nu', 'a', 'elower', 'gp', 'gpp'])


def fetch_hitran(molecule, lam_min_aa=900.0, lam_max_aa=3e5, api_key=None):
    """Pull all transitions of every isotopologue of `molecule` via hapi2."""
    from hapi2 import SETTINGS, Molecule, fetch_molecules, fetch_isotopologues, \
        fetch_transitions

    key = api_key or os.environ.get('HAPI2_API_KEY')
    if not key:
        raise SystemExit(
            'no HITRAN API key: pass --api-key or set HAPI2_API_KEY')
    SETTINGS['api_key'] = key

    numin = int(1e8 / lam_max_aa)
    numax = int(1e8 / lam_min_aa) + 1

    fetch_molecules()
    fetch_isotopologues([Molecule(molecule)])
    fetch_transitions(Molecule(molecule).isotopologues, numin=numin,
                      numax=numax, llst_name=molecule + '_llist')
    return pd.DataFrame([t.dump() for t in Molecule(molecule).transitions])


# ---------------------------------------------------------------------------
# radiative damping from Einstein A coefficients
# ---------------------------------------------------------------------------

def _quanta_keys(df):
    """Level identifiers from HITRAN's own quanta strings.

    The same physical level carries byte-identical (global, local) quanta
    wherever it appears -- as the upper state of one transition and the lower
    state of another.  That exact-match property is what lets the lower-level
    sum be assembled at all.
    """
    def norm(col):
        return df[col].astype(str).str.split().str.join(' ')

    iso = df['local_iso_id'].astype(str)
    up = iso + '|' + norm('global_upper_quanta') + '|' + norm('local_upper_quanta')
    lo = iso + '|' + norm('global_lower_quanta') + '|' + norm('local_lower_quanta')
    return up, lo


def _energy_keys(df, tol=1e-3):
    """Level identifiers from term energy, for lists with unusable quanta.

    Clusters all level energies (upper and lower alike) by single-linkage on a
    sorted array rather than rounding to a grid: rounding would split a level
    whose two occurrences straddle a bin edge, which is the failure mode that
    matters here, since a split level loses part of its A sum.
    """
    e_up = df['elower'].values + df['nu'].values
    e_lo = df['elower'].values
    iso = df['local_iso_id'].values

    keys = []
    for e in (e_up, e_lo):
        keys.append(np.stack([iso, e], axis=1))
    allk = np.concatenate(keys)

    order = np.lexsort((allk[:, 1], allk[:, 0]))
    s = allk[order]
    new = np.ones(len(s), dtype=bool)
    new[1:] = (s[1:, 0] != s[:-1, 0]) | (np.diff(s[:, 1]) > tol)
    cluster = np.empty(len(s), dtype=np.int64)
    cluster[order] = np.cumsum(new) - 1

    n = len(df)
    return pd.Series(cluster[:n]), pd.Series(cluster[n:])


def _level_keys(df, method='auto', verbose=True):
    """Pick a level-identification scheme; see the two helpers above.

    Quanta are preferred when they work, being exact.  But bulk HITEMP files
    are not uniform in this: some carry blank or truncated quanta fields, and a
    blank field silently collapses every level of the molecule into one key,
    which would hand back a single enormous gamma for the whole list.  Detect
    that by checking whether lines sharing a key agree on upper-state energy,
    and fall back on energy clustering if they do not.
    """
    if method == 'energy':
        return _energy_keys(df)

    up, lo = _quanta_keys(df)
    e_up = df['elower'].values + df['nu'].values
    spread = pd.Series(e_up).groupby(up.values).agg(lambda x: x.max() - x.min())
    frac_bad = float((spread > 1e-3).mean()) if len(spread) else 1.0

    if method == 'quanta':
        return up, lo
    if frac_bad > 0.01:
        if verbose:
            print(f"  quanta strings resolve levels poorly "
                  f"({100*frac_bad:.1f}% of keys span >1e-3 cm^-1); "
                  f"falling back to energy clustering")
        return _energy_keys(df)
    return up, lo


def radiative_damping(df, method='auto', verbose=True):
    """gamma_rad [1/s] per line, as sum_k A_uk + sum_k A_lk over the list.

    Returns the damping array plus a dict of diagnostics describing how much of
    each sum is actually populated, since that is the only handle on how much
    of a lower bound the answer is.
    """
    up, lo = _level_keys(df, method=method, verbose=verbose)

    # inverse lifetime of every level that appears as an upper state
    a_out = df.groupby(up.values)['a'].sum()

    gamma_up = up.map(a_out).fillna(0.0).values
    gamma_lo = lo.map(a_out).fillna(0.0).values
    gamma = gamma_up + gamma_lo

    # Residual check on whatever keying survived.  Under energy clustering a
    # few keys legitimately exceed the tolerance -- single-linkage chains
    # through near-degenerate levels -- so this is a scale indicator, not an
    # error: a couple of percent is normal, a large fraction means the levels
    # are not being resolved and the sums are being over-collected.
    e_up = df['elower'].values + df['nu'].values
    spread = pd.Series(e_up).groupby(up.values).agg(lambda x: x.max() - x.min())
    bad = (spread > 1e-3).sum()

    diag = {
        'n_lines': len(df),
        'n_levels': len(a_out),
        'lower_populated': float(np.mean(gamma_lo > 0)),
        'inconsistent_levels': int(bad),
        'gamma_median': float(np.median(gamma[gamma > 0])) if (gamma > 0).any() else 0.0,
        'gamma_max': float(gamma.max()),
    }
    if verbose:
        print(f"  levels identified      : {diag['n_levels']}")
        print(f"  lines with lower-state A: {100*diag['lower_populated']:.1f}%")
        print(f"  gamma_rad median / max : {diag['gamma_median']:.3e} / "
              f"{diag['gamma_max']:.3e} 1/s")
        if bad:
            print(f"  {bad} of {diag['n_levels']} level keys "
                  f"({100*bad/max(diag['n_levels'],1):.1f}%) span >1e-3 cm^-1 "
                  f"in upper-state energy")
    return gamma, diag


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

def build_linelist(df, tref=4000.0, strength_cut=None, lam_min=None,
                   lam_max=None, verbose=True):
    """HITRAN transitions -> master-format columns, sorted by air wavelength."""
    df = df.reset_index(drop=True)

    gamma, diag = radiative_damping(df, verbose=verbose)

    # air-pressure shift is a laboratory quantity and does not belong in a
    # stellar list; use the unshifted line centre.
    lam_vac = 1e8 / df['nu'].values                       # AA
    lam_air = vac2air(lam_vac)
    evi = CM1_TO_EV * df['elower'].values

    # g_lo f_lu from the Einstein A of the transition.  Note this uses gp, the
    # UPPER state degeneracy -- the gf/A relation is g_lo f_lu = k lam^2 g_up A.
    gf = GF_CONST * lam_vac**2 * df['gp'].values * df['a'].values
    with np.errstate(divide='ignore'):
        loggf = np.log10(gf)

    out = pd.DataFrame({
        'lam0': lam_air,
        'evi': evi,
        'loggf': loggf,
        'gv': 0.0,
        'ga': np.where(gamma > 0, gamma * GA_ENCODE, 0.0),
        'gq': 0.0,
    })
    out = out[np.isfinite(out['loggf'])]

    if lam_min is not None:
        out = out[out['lam0'] >= lam_min]
    if lam_max is not None:
        out = out[out['lam0'] <= lam_max]

    if strength_cut is not None:
        # Boltzmann-weighted oscillator strength at tref, relative to the
        # strongest line in the list.  This is the meaningful way to prune a
        # HITEMP file: a raw loggf cut discards strong hot-band lines whose
        # small A is compensated by a highly populated lower state, and keeps
        # cold lines that never contribute at photospheric temperatures.
        logs = out['loggf'] - out['evi'] * (5040.0 / tref)
        keep = logs >= logs.max() + strength_cut
        if verbose:
            print(f"  strength cut {strength_cut:+.1f} dex at {tref:.0f} K: "
                  f"{len(out)} -> {int(keep.sum())} lines")
        out = out[keep]

    return out.sort_values('lam0').reset_index(drop=True), diag


def _cluster(energies, tol=1e-3):
    """Single-linkage cluster of a 1-D energy array; returns (sorted_uniq, ids).

    Membership is looked up later with searchsorted, which is exact because
    every query energy is itself one of the values that went into the cluster.
    """
    uniq = np.unique(energies)
    ids = np.cumsum(np.r_[True, np.diff(uniq) > tol]) - 1
    return uniq, ids


def build_linelist_streaming(filename, iso=None, tref=4000.0, strength_cut=None,
                             lam_min=900.0, lam_max=200000.0, tol=1e-3,
                             chunksize=32_000_000, verbose=True):
    """Two-pass, memory-bounded version of build_linelist for the big files.

    CH4 HITEMP2020 is ~31e6 records and CO2 HITEMP2024 more again; holding
    those as one DataFrame (with the quanta strings) is gigabytes for no gain.
    Instead:

      pass 1  accumulate sum(A) per exact upper-state energy, collect the set
              of level energies, and track the peak line strength;
      then    cluster the energies globally -- this is why it cannot be done in
              one pass, since a level's A sum is only complete once the whole
              file has been seen;
      pass 2  re-stream, look each line's two levels up in the cluster table,
              and keep only lines surviving the strength cut.

    Both passes see identical E = elower + nu values for a given record, so the
    searchsorted lookup in pass 2 always hits.
    """
    numin, numax = 1e8 / lam_max, 1e8 / lam_min
    kw = dict(numin=numin, numax=numax, iso=iso, fields=PAR_NUMERIC,
              chunksize=chunksize)

    # ---- pass 1 -----------------------------------------------------------
    e_parts, a_parts, elo_parts = [], [], []
    n_in, max_logs, next_report = 0, -np.inf, 0
    for df in iter_par(filename, **kw):
        n_in += len(df)
        e_up = df['elower'].values + df['nu'].values
        # collapse within the chunk so what we retain is levels, not lines
        u, inv = np.unique(e_up, return_inverse=True)
        e_parts.append(u)
        a_parts.append(np.bincount(inv, weights=df['a'].values, minlength=len(u)))
        elo_parts.append(np.unique(df['elower'].values))

        lam_vac = 1e8 / df['nu'].values
        gf = GF_CONST * lam_vac**2 * df['gp'].values * df['a'].values
        with np.errstate(divide='ignore', invalid='ignore'):
            logs = np.log10(gf) - CM1_TO_EV * df['elower'].values * (5040.0 / tref)
        if np.isfinite(logs).any():
            max_logs = max(max_logs, np.nanmax(logs[np.isfinite(logs)]))
        if verbose and n_in // 5_000_000 > next_report:
            next_report = n_in // 5_000_000
            print(f"    pass 1: {n_in/1e6:.1f}M lines read", flush=True)

    if not e_parts:
        raise SystemExit('no transitions in range')

    e_bucket = np.concatenate(e_parts)
    a_bucket = np.concatenate(a_parts)
    # a level can appear in several chunks; merge those partial sums
    e_bucket, inv = np.unique(e_bucket, return_inverse=True)
    a_bucket = np.bincount(inv, weights=a_bucket, minlength=len(e_bucket))

    all_e = np.concatenate([e_bucket, np.concatenate(elo_parts)])
    uniq, cid = _cluster(all_e, tol=tol)
    a_level = np.bincount(cid[np.searchsorted(uniq, e_bucket)],
                          weights=a_bucket, minlength=cid[-1] + 1)

    thresh = max_logs + strength_cut if strength_cut is not None else -np.inf
    if verbose:
        print(f"  pass 1: {n_in} transitions, {len(a_level)} levels, "
              f"peak strength {max_logs:.2f}")
        if strength_cut is not None:
            print(f"  keeping lines with strength >= {thresh:.2f} "
                  f"({strength_cut:+g} dex at {tref:.0f} K)")

    # ---- pass 2 -----------------------------------------------------------
    keep_cols, n_kept, gam_lo_hits, n_gam = [], 0, 0, 0
    for df in iter_par(filename, **kw):
        e_up = df['elower'].values + df['nu'].values
        e_lo = df['elower'].values
        g_up = a_level[cid[np.searchsorted(uniq, e_up)]]
        g_lo = a_level[cid[np.searchsorted(uniq, e_lo)]]
        gamma = g_up + g_lo
        gam_lo_hits += int((g_lo > 0).sum()); n_gam += len(gamma)

        lam_vac = 1e8 / df['nu'].values
        evi = CM1_TO_EV * e_lo
        gf = GF_CONST * lam_vac**2 * df['gp'].values * df['a'].values
        with np.errstate(divide='ignore', invalid='ignore'):
            loggf = np.log10(gf)
        keep = np.isfinite(loggf) & ((loggf - evi * (5040.0 / tref)) >= thresh)
        if keep.any():
            keep_cols.append((vac2air(lam_vac[keep]), evi[keep], loggf[keep],
                              gamma[keep]))
            n_kept += int(keep.sum())

    lam, evi, loggf, gamma = (np.concatenate(c) for c in zip(*keep_cols))
    out = pd.DataFrame({
        'lam0': lam, 'evi': evi, 'loggf': loggf, 'gv': 0.0,
        'ga': np.where(gamma > 0, gamma * GA_ENCODE, 0.0), 'gq': 0.0,
    })
    # input is nu-ascending, so this is wavelength-descending; reversing is
    # enough, but verify rather than assume -- the reader bisects on order.
    out = out[::-1].reset_index(drop=True)
    if not out['lam0'].is_monotonic_increasing:
        out = out.sort_values('lam0').reset_index(drop=True)

    diag = {'n_lines': n_in, 'n_levels': len(a_level), 'n_kept': n_kept,
            'lower_populated': gam_lo_hits / max(n_gam, 1)}
    if verbose:
        nz = gamma[gamma > 0]
        print(f"  pass 2: kept {n_kept} of {n_in} lines")
        print(f"  lines with lower-state A: {100*diag['lower_populated']:.1f}%")
        if len(nz):
            print(f"  gamma_rad median / max : {np.median(nz):.3e} / "
                  f"{nz.max():.3e} 1/s")
    return out, diag


def write_linelist(out, filename, species, comment, iso='0'):
    """Write one block in the format init_linelist2 / LineList parse."""
    with open(filename, 'w') as fh:
        fh.write(f"'{species}' '{comment}'\n")
        fh.write(f"'{iso}' {len(out)}\n")
        for r in out.itertuples(index=False):
            fh.write(f"{r.lam0:12.3f} {r.evi:7.3f} {r.loggf:8.3f} "
                     f"{r.gv:.0f} {r.ga:9.3e} {r.gq:.0f} 0 0\n")


# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--par', help='local HITRAN/HITEMP .par file')
    src.add_argument('--fetch', help='HITRAN molecule to download via hapi2')
    p.add_argument('--id', required=True,
                   help='species id for the header, as in molec_data_BC.txt')
    p.add_argument('--out', required=True, help='output line list')
    p.add_argument('--iso', type=int, default=None,
                   help='keep only this HITRAN local_iso_id (default: all)')
    p.add_argument('--comment', default=None, help='header provenance string')
    p.add_argument('--lam-min', type=float, default=900.,
                   help='AA, air (default 900)')
    p.add_argument('--lam-max', type=float, default=200000.,
                   help='AA, air (default 200000, as the existing master lists)')
    p.add_argument('--tref', type=float, default=4000.0,
                   help='reference temperature for the strength cut [K]')
    p.add_argument('--strength-cut', type=float, default=None,
                   help='keep lines within this many dex of the strongest '
                        'Boltzmann-weighted line at --tref, e.g. -8')
    p.add_argument('--stream', dest='stream', action='store_true', default=None,
                   help='force the two-pass memory-bounded path')
    p.add_argument('--no-stream', dest='stream', action='store_false',
                   help='force the single-pass in-memory path')
    p.add_argument('--api-key', default=None)
    args = p.parse_args(argv)

    common = dict(tref=args.tref, strength_cut=args.strength_cut,
                  lam_min=args.lam_min, lam_max=args.lam_max)

    if args.par:
        comment = args.comment or f'HITEMP {os.path.basename(args.par)}'
        if args.iso is not None:
            comment += f' (iso {args.iso})'
        # anything past a few million records goes down the two-pass path;
        # .bz2 size is a rough proxy, HITEMP compresses ~8x
        big = os.path.getsize(args.par) > 50e6
        stream = big if args.stream is None else args.stream
        print(f"reading {args.par} ({'streaming' if stream else 'in-memory'})")
        if stream:
            out, _ = build_linelist_streaming(args.par, iso=args.iso, **common)
        else:
            # window in wavenumber while reading, so the out-of-range far-IR
            # tail never becomes rows; 1e8/lam[AA] = nu[cm^-1]
            df = read_par(args.par, numin=1e8 / args.lam_max,
                          numax=1e8 / args.lam_min)
            if args.iso is not None:
                df = df[df['local_iso_id'] == args.iso]
            print(f"  {len(df)} transitions")
            out, _ = build_linelist(df, **common)
    else:
        print(f"fetching {args.fetch} from HITRAN")
        df = fetch_hitran(args.fetch, api_key=args.api_key)
        comment = args.comment or f'HITRAN online, {args.fetch}'
        if args.iso is not None:
            df = df[df['local_iso_id'] == args.iso]
            comment += f' (iso {args.iso})'
        print(f"  {len(df)} transitions")
        out, _ = build_linelist(df, **common)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    write_linelist(out, args.out, args.id, comment)
    print(f"wrote {len(out)} lines to {args.out} "
          f"({out['lam0'].min():.1f} - {out['lam0'].max():.1f} AA)")


if __name__ == '__main__':
    sys.exit(main())
