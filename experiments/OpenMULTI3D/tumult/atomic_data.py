"""Readers for the reference data that model atoms are checked against and built from.

- NIST ASD level tables, saved as CSV from the levels form
  (https://physics.nist.gov/PhysRefData/ASD/levels_form.html, format "CSV (to file)"),
  one file per spectrum named like ``MgI_levels.csv``.
- Turbospectrum-format line lists (:meth:`tumult.linelist.LineList.from_ts_format`).
- TOPbase photoionisation tables (``p`` files from
  https://cdsweb.u-strasbg.fr/topbase/topbase.html) and the Fe I cross sections of
  Bautista et al. (2017, ApJ 838, 152; CDS J/ApJ/838/152).
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .const import CC, EE, HH
from .linelist import LineList

CM_TO_EV = HH * CC / EE
RYD_EV = 13.605693
RYD_AA = 911.267                 # vacuum wavelength of 1 Ry above the threshold [AA]
MB = 1e-18                       # megabarn [cm^2]
ROMAN = ('I', 'II', 'III', 'IV', 'V', 'VI')
_NUMBER = re.compile(r'[-+]?\d+\.?\d*')


def element_symbol(atom_id):
    """'MG', 'Mg' or ' mg ' -> 'Mg'."""
    return atom_id.strip().capitalize()


def _unwrap(cell):
    """NIST CSV cells are written as ="..." so that spreadsheets keep them as text."""
    if not isinstance(cell, str):
        return cell
    cell = cell.strip()
    if cell.startswith('="') and cell.endswith('"'):
        cell = cell[2:-1]
    return cell.strip()


def read_nist_levels(path):
    """Fine-structure levels of one spectrum from a NIST ASD CSV.

    Returns the levels (columns config, term, J, g, ecm, ev) sorted by energy, and the
    first ionisation limit [eV]. Levels with an uncertain or relative energy (brackets,
    '+x') keep their numerical value; rows without one are dropped.
    """
    raw = pd.read_csv(path, dtype=str).map(_unwrap)
    limit = raw[raw.Term == 'Limit']
    limit_cm = float(_NUMBER.search(limit['Level (cm-1)'].iloc[0]).group()) if len(limit) else np.nan
    lvl = raw[(raw.Term != 'Limit') & raw.g.notna()].copy()
    ecm = lvl['Level (cm-1)'].map(lambda s: float(m.group()) if isinstance(s, str) and (m := _NUMBER.search(s)) else np.nan)
    lvl = pd.DataFrame(dict(config=lvl.Configuration.values, term=lvl.Term.values, J=lvl.J.values,
                            g=pd.to_numeric(lvl.g, errors='coerce').values, ecm=ecm.values))
    lvl = lvl.dropna(subset=['g', 'ecm'])
    lvl = lvl[lvl.ecm < limit_cm] if np.isfinite(limit_cm) else lvl
    lvl['ev'] = lvl.ecm * CM_TO_EV
    return lvl.sort_values('ecm', ignore_index=True), limit_cm * CM_TO_EV


def nist_levels_for_atom(folder, atom_id, stages):
    """NIST levels of the given ion stages (1 = neutral) of an element as one table with
    an ``ion`` column, and the ionisation energy of each stage [eV]."""
    sym = element_symbol(atom_id)
    tabs = []
    limits = {}
    for stage in stages:
        path = Path(folder) / f'{sym}{ROMAN[stage - 1]}_levels.csv'
        if not path.exists():
            continue
        lvl, lim = read_nist_levels(path)
        lvl['ion'] = stage
        lvl['limit_ev'] = lim
        tabs.append(lvl)
        limits[stage] = lim
    return pd.concat(tabs, ignore_index=True), limits


def _vald_designation(row):
    """'  LS        3d6.(5D).4s\\ (4D).5s g5D' -> ('LS', '3d6.(5D).4s\\ (4D).5s', 'g5D').
    Extracts returned by email quote the row; the quotes are dropped."""
    part = row.strip().strip("'").split(None, 1)
    if len(part) < 2:
        return (part[0] if part else ''), '', ''
    cfg, _, term = part[1].rpartition(' ')
    return part[0], cfg.strip(), term


def read_vald_long(path, species=None):
    """VALD3 extract in long format (air wavelengths, energies in eV; ``.gz`` read as
    such): one row per transition with element, ion, lam [A], loggf, elow, eup [eV],
    j_lo, j_up, lgrad, stark, vdw (VALD's damping columns), coupling, config and term of
    both levels, and the reference row. species: e.g. {'Fe 1', 'Fe 2'} to keep only those.
    """
    import gzip
    path = Path(path)
    text = (gzip.open(path, 'rt') if path.suffix == '.gz' else open(path)).read().splitlines()
    rows = []
    for k, r in enumerate(text):
        if not re.match(r"'[A-Z][a-z]? +\d+'", r):
            continue
        p = [x.strip() for x in r.split(',')]
        name = p[0].strip("'")
        if species and name not in species:
            continue
        el, ion = name.split()
        lo, up = _vald_designation(text[k + 1]), _vald_designation(text[k + 2])
        rows.append(dict(element=el, ion=int(ion), lam=float(p[1]), loggf=float(p[2]), elow=float(p[3]),
                         j_lo=float(p[4]), eup=float(p[5]), j_up=float(p[6]), lgrad=float(p[10]),
                         stark=float(p[11]), vdw=float(p[12]), coupling_lo=lo[0], config_lo=lo[1], term_lo=lo[2],
                         coupling_up=up[0], config_up=up[1], term_up=up[2],
                         ref=' '.join(text[k + 3].strip("'_ ").split())))
    return pd.DataFrame(rows)


def linelist_by_stage(path, atom_id, stages):
    """Lines of an element from a Turbospectrum line list, one table per ion stage with
    columns lam0 [AA], evi [eV], loggf. Isotopic components stay separate rows."""
    llist = LineList.from_ts_format(path)
    sym = element_symbol(atom_id)
    out = {}
    for stage in stages:
        label = f'{sym} {ROMAN[stage - 1]}'
        if label in llist.index.get_level_values(0):
            out[stage] = llist.loc[label].reset_index(drop=True)
    return out


def ts_vacuum_wavelength(lam):
    """Vacuum wavelength of Turbospectrum/VALD line-list entries: air above 2000 AA,
    vacuum below."""
    from .physics import air2vac
    lam = np.asarray(lam, float)
    return np.where(lam < 2000, lam, air2vac(np.maximum(lam, 2000)))


def read_topbase_photo(path):
    """TOPbase photoionisation file (``p`` table).

    Returns {(2S+1, L, parity, ilv): (e_photon, sigma)}: the cross section sigma [Mb]
    at photon energies [Ry], starting at the threshold of the Opacity Project's own
    level energy. OP energies are not experimental; shift a table by the difference of
    the thresholds before using it for a level with known energy.
    """
    rows = [r.split() for r in Path(path).read_text().splitlines() if r.strip()]
    out = {}
    k = 1                                    # first row: Z, number of electrons, 'P'
    while k + 2 < len(rows):
        mult, ell, par, ilv = map(int, rows[k])
        npts = int(rows[k + 1][1])
        data = np.array(rows[k + 3:k + 3 + npts], dtype=float)
        out[(mult, ell, par, ilv)] = (data[:, 0], data[:, 1])
        k += 3 + npts
    return out


BAUTISTA_TEMPS = np.array([3000., 5000., 8000., 10000., 15000., 20000.])


def read_bautista_ecs(folder):
    """Fe I effective collision strengths of Bautista et al. (2017) between LS terms:
    the term table (index -> term label such as 'a5D' or 'z7D*') and a table of
    (lower, upper, upsilon at BAUTISTA_TEMPS)."""
    terms = {}
    for row in (Path(folder) / 'table4.dat').read_text().splitlines():
        term = row[27:36].replace('^', '').replace(' ', '')
        terms[int(row[0:2])] = term[:-1] + '*' if term.endswith('o') else term
    ups = np.loadtxt(Path(folder) / 'table5.dat')
    tab = pd.DataFrame(ups[:, 2:], columns=BAUTISTA_TEMPS)
    tab.insert(0, 'upper', [terms[int(k)] for k in ups[:, 1]])
    tab.insert(0, 'lower', [terms[int(k)] for k in ups[:, 0]])
    return tab
