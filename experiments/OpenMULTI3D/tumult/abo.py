"""Python replication of Turbospectrum's ABO collisional broadening (anstee.f).

Turbospectrum resolves van der Waals broadening for atomic lines in `anstee.f`,
which is only reachable when the line list carries fdamp < 20. It has two paths:

  1. Bilinear interpolation of the Anstee/Barklem/O'Mara tables (recipes A/B/C)
     for neutral s-p, p-s, p-d, d-p, d-f and f-d transitions, keyed on the
     effective principal quantum numbers of the two levels.
  2. A hardcoded per-line override (recipe S) held in three tables inside
     DATA/Anstee-*.dat -- 5024 "special" lines, 24188 Fe II and 13165 Cr II --
     matched on wavelength and lower-level energy alone. This takes precedence.

Both paths are reproduced here so the coefficients can be baked into a line list
as packed `sigma.alpha` (fdamp >= 20), which makes Turbospectrum skip anstee.f
entirely and lets other codes see the same numbers.

Reference: Anstee & O'Mara 1995; Barklem & O'Mara 1997, 1998; Barklem, O'Mara &
Ross 1998; Barklem, Piskunov & O'Mara 2000; Barklem & Aspelund-Johansson 2005.
"""
import numpy as np
import pandas as pd
import re

TS_ROOT = '/mnt/beegfs/gemini/groups/bergemann/users/hoppe/TSFitPy/turbospectrum'
ANSTEE_DAT = f'{TS_ROOT}/DATA/Anstee-1802-2014.dat'
IRWIN_DAT = f'{TS_ROOT}/DATA/IRWIN_atoms_v07.3.dat'

# anstee.f matching tolerances and its cm^-1 -> eV conversion
DLAM, DEV, CM_PER_EV = 0.05, 0.02, 8067.0
RYD_H = 13.595                 # H I ionisation energy, homogenised in anstee.f
HC_EV_A = 12399.8              # anstee.f line 79: eup = 12399.8/xlb + elo

# Fixed-width layout of the per-line blocks: format(i2,3x,i2,f10.3,f11.3,59x,f8.2,f8.3)
_SPECIAL_COLS = [(0, 2), (5, 7), (7, 17), (17, 28), (87, 95), (95, 103)]
_SPECIAL_NAMES = ['iel', 'ion', 'xlb', 'elo', 'sigma', 'alpha']
# (first data line, number of rows) -- verified against the nspecial/nfe2/ncr2 parameters
_BLOCKS = {'special': (226, 5024), 'fe2': (5255, 24188), 'cr2': (29449, 13165)}

# Grid tables: (name, row index range, column index range). Indices are 10*n*.
# table1/2 s-p, table3/4 p-d, table5/6 d-f; sigma tables are odd, alpha even.
_GRIDS = {
    'sp': dict(sig=1, alp=2, i=(10, 30), j=(13, 30)),
    'pd': dict(sig=3, alp=4, i=(13, 30), j=(23, 40)),
    'df': dict(sig=5, alp=6, i=(23, 40), j=(33, 50)),
}

# Transition -> (grid, which of the pair is "nstar1"). anstee.f always passes the
# lower-l orbital as nstar1 regardless of which end of the transition it is on.
_TRANS = {
    ('s', 'p'): ('sp', 'lo'), ('p', 's'): ('sp', 'up'),
    ('p', 'd'): ('pd', 'lo'), ('d', 'p'): ('pd', 'up'),
    ('d', 'f'): ('df', 'lo'), ('f', 'd'): ('df', 'up'),
}
_RECIPE = {'sp': 'A', 'pd': 'B', 'df': 'C'}
# valid n* windows, mirroring the guards around each interpolation call
_RANGE = {'sp': ((1.0, 3.0), (1.3, 3.0)), 'pd': ((1.3, 3.0), (2.3, 4.0)),
          'df': ((2.3, 4.0), (3.3, 5.0))}


def _read_grids(path=ANSTEE_DAT):
    """Parse the six n*-indexed tables from the head of the Anstee data file."""
    rows = open(path).readlines()
    nums = re.compile(r'^\s*[-+]?[\d.]')   # alpha tables start with a bare ".246"
    grids, k = {}, 0
    for t in range(1, 7):
        g = _GRIDS['sp' if t <= 2 else 'pd' if t <= 4 else 'df']
        i0, i1 = g['i']
        j0, j1 = g['j']
        while not nums.match(rows[k]):     # skip the commentary before each table
            k += 1
        block = []
        for _ in range(i1 - i0 + 1):
            block.append([float(x) for x in rows[k].split()])
            k += 1
        arr = np.array(block)
        assert arr.shape == (i1 - i0 + 1, j1 - j0 + 1), (t, arr.shape)
        grids[t] = arr
    return grids


def _read_special(path=ANSTEE_DAT):
    tabs = {}
    for name, (start, n) in _BLOCKS.items():
        df = pd.read_fwf(path, colspecs=_SPECIAL_COLS, names=_SPECIAL_NAMES,
                         skiprows=start - 1, nrows=n)
        assert len(df) == n and not df.isna().any().any(), name
        df['elo_ev'] = df.elo / CM_PER_EV
        tabs[name] = df.reset_index(drop=True)
    return tabs


def _read_ionpot(path=IRWIN_DAT):
    """Neutral ionisation energies, eV. partf.f stores the I->II potential on the
    singly-ionised record (spec = Z.01), so that is what we pick up."""
    # species names are quoted and may contain blanks, so match around them
    rec = re.compile(r"^\s*\d+\s+(\d+\.\d+)\s+'[^']*'\s+([-+]?[\d.eE]+)")
    ip = {}
    for ln in open(path):
        m = rec.match(ln)
        if not m:
            continue
        spec, pot = float(m.group(1)), float(m.group(2))
        z = int(spec)
        if abs(spec - z - 0.01) < 1e-6:
            ip[z] = pot
    return ip


class ABO:
    def __init__(self):
        self.grids = _read_grids()
        self.tabs = _read_special()
        self.ionpot = _read_ionpot()

    # ---- path 1: table interpolation (recipes A/B/C) -------------------------
    def _interp(self, key, n1, n2):
        """Bilinear interpolation on the 0.1-spaced n* grid, as anstee.f does."""
        g = _GRIDS[key]
        a, b = self.grids[g['sig']], self.grids[g['alp']]
        i0, j0 = g['i'][0], g['j'][0]

        i1 = (10.0 * n1).astype(int)
        i2 = (10.0 * n2).astype(int)
        f1 = 10.0 * n1 - i1
        f2 = 10.0 * n2 - i2
        r, c = i1 - i0, i2 - j0
        # guard the +1 neighbour at the top edge; fract is 0 there so it is unused
        rc = np.clip(r, 0, a.shape[0] - 2)
        cc = np.clip(c, 0, a.shape[1] - 2)

        def bilin(t):
            x = t[rc, cc] * (1 - f2) + t[rc, cc + 1] * f2
            y = t[rc + 1, cc] * (1 - f2) + t[rc + 1, cc + 1] * f2
            return x * (1 - f1) + y * f1

        return bilin(a), bilin(b)

    def interpolated(self, iel, ion, lam, elo, levlo, levup):
        """Vectorised recipes A/B/C. Returns (sigma, alpha, ok) with ok=False
        where the transition type or n* range is not covered."""
        n = len(lam)
        sigma = np.full(n, np.nan)
        alpha = np.full(n, np.nan)
        recipe = np.full(n, 'U', dtype='<U1')

        eion = np.array([self.ionpot.get(z, np.nan) for z in iel])
        eup = HC_EV_A / lam + elo
        neutral = (ion == 1) & np.isfinite(eion)

        lo = np.char.lower(levlo.astype(str))
        up = np.char.lower(levup.astype(str))

        with np.errstate(invalid='ignore', divide='ignore'):
            nlo = 1.0 / np.sqrt((eion - elo) / RYD_H)
            nup = 1.0 / np.sqrt((eion - eup) / RYD_H)

        for (a, b), (key, which) in _TRANS.items():
            m = neutral & (lo == a) & (up == b) & np.isfinite(nlo) & np.isfinite(nup)
            if not m.any():
                continue
            # nstar1 is the lower-l orbital of the pair
            n1 = np.where(which == 'lo', nlo, nup)[m]
            n2 = np.where(which == 'lo', nup, nlo)[m]
            (l1, h1), (l2, h2) = _RANGE[key]
            good = (n1 >= l1) & (n1 <= h1) & (n2 >= l2) & (n2 <= h2)
            if not good.any():
                continue
            s, v = self._interp(key, n1[good], n2[good])
            idx = np.nonzero(m)[0][good]
            sigma[idx], alpha[idx], recipe[idx] = s, v, _RECIPE[key]

        # last resort in anstee.f: d-p lines with nstard < 2.0 use the upper p
        # level alone (Barklem, priv. comm. 1997), piecewise-linear in nstarp
        m = neutral & (recipe == 'U') & (lo == 'd') & (up == 'p') & (nlo < 2.0)
        if m.any():
            p = nup[m]
            s = np.full(p.shape, np.nan)
            v = np.full(p.shape, np.nan)
            for lo_e, hi_e, s0, ds, v0, dv in [
                    (2.6, 2.7, 762., 960.0, 0.270, -0.04),
                    (2.7, 2.8, 858., 1110., 0.266, +0.02),
                    (2.8, 2.9, 969., 610.0, 0.268, -0.10),
                    (2.9, 3.0, 1030., 440.0, 0.258, -0.04)]:
                # first bin is closed at both ends, the rest are half-open
                w = (p >= lo_e) & (p <= hi_e) if lo_e == 2.6 else (p > lo_e) & (p <= hi_e)
                s[w] = s0 + ds * (p[w] - lo_e)
                v[w] = v0 + dv * (p[w] - lo_e)
            idx = np.nonzero(m)[0]
            ok = np.isfinite(s)
            sigma[idx[ok]], alpha[idx[ok]], recipe[idx[ok]] = s[ok], v[ok], 'S'

        return sigma, alpha, recipe

    # ---- path 2: hardcoded per-line override (recipe S) ----------------------
    def special(self, iel, ion, lam, elo):
        """Per-line lookup, |dlam| <= 0.05 A and |dE| < 0.02 eV, species-branched.
        First match in original table order wins, mirroring anstee.f's goto."""
        sigma = np.full(len(lam), np.nan)
        alpha = np.full(len(lam), np.nan)

        for (z, i), idx in pd.DataFrame({'z': iel, 'i': ion}).groupby(['z', 'i']).groups.items():
            idx = np.asarray(idx)
            if (z, i) == (26, 2):
                tab = self.tabs['fe2']
            elif (z, i) == (24, 2):
                tab = self.tabs['cr2']
            else:
                tab = self.tabs['special']
                tab = tab[(tab.iel == z) & (tab.ion == i)]
                if not len(tab):
                    continue
            order = np.argsort(tab.xlb.values, kind='stable')
            xlb = tab.xlb.values[order]
            eev = tab.elo_ev.values[order]
            sg, al = tab.sigma.values[order], tab.alpha.values[order]
            rank = tab.index.values[order]          # original table position

            l, e = lam[idx], elo[idx]
            klo = np.searchsorted(xlb, l - DLAM, 'left')
            khi = np.searchsorted(xlb, l + DLAM, 'right')
            for k in np.nonzero(khi > klo)[0]:
                sl = slice(klo[k], khi[k])
                hit = np.nonzero(np.abs(e[k] - eev[sl]) < DEV)[0]
                if not len(hit):
                    continue
                j = hit + klo[k]
                best = j[np.argmin(rank[j])]
                sigma[idx[k]], alpha[idx[k]] = sg[best], al[best]
        return sigma, alpha

    # ---- combined ------------------------------------------------------------
    def damping(self, iel, ion, lam, elo, levlo, levup):
        """Full anstee.f result. Returns (sigma, alpha, recipe); recipe 'U' means
        no ABO data and Turbospectrum would fall back to Unsoeld."""
        sigma, alpha, recipe = self.interpolated(iel, ion, lam, elo, levlo, levup)
        s2, a2 = self.special(iel, ion, lam, elo)
        ok = np.isfinite(s2)                        # the override wins where present
        sigma[ok], alpha[ok], recipe[ok] = s2[ok], a2[ok], 'S'
        return sigma, alpha, recipe


def packed(sigma, alpha):
    """Pack into the fdamp >= 20 form Turbospectrum reads as sigma.alpha."""
    return np.where(np.isfinite(sigma), np.floor(sigma) + alpha, np.nan)
