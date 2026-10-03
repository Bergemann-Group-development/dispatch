"""Gaia XP stellar parameters for the 4MOST HRS stars, from the S4 target tables.

    python s4_labels.py PHOT_DIR S4_DIR DISK_CSV CRUDE_LABELS OUT.csv

PHOT_DIR holds the <set>_photometry.csv files of the HRS spectra (one row per star:
OBJECT, RA, DECL, the spectrum ids, and for 0301/0401/0408 the Gaia DR3 designation).
Each star gets a Gaia DR3 source_id from DR3Name, or else by position (1 arcsec) in the
S4 target cross-match XMatch_S4_GXP.fits (S4_DIR).  Labels, by source_id:

  disk:  the S4 4MIDABLE-HR disk-sample table (DISK_CSV, March 2026): Teff, log g, [M/H]
         of the XGBoost model of Andrae, Rix & Chandra (2023) on Gaia XP spectra, and
         Gaia GSP-Phot Teff and log g (its mg_gspphot is the absolute G magnitude, not a
         metallicity);
  s4xm:  the XGBoost labels of the S4 target cross-match (XMatch_S4_GXP.fits, 2024);
  andrae: the full XGBoost catalogue (table_1_catwise.fits), for stars outside the S4
          target list.

None of the tables quotes per-star uncertainties.  CRUDE_LABELS is place_stars.py's
labels.csv; OUT gets one row per star with the crude fit of its highest-S/N spectrum
(teff_crude, logg_crude, feh_crude: log g held at Gaia where possible) for comparison.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

SETS = ('0201', '0301', '0401', '0408')


def stars(phot_dir):
    rows = []
    for s in SETS:
        d = pd.read_csv(Path(phot_dir) / f'{s}_photometry.csv')
        sid = (d['DR3Name'].str.extract(r'(\d+)$')[0].astype('Int64') if 'DR3Name' in d
               else pd.Series(pd.NA, index=d.index, dtype='Int64'))
        rows.append(pd.DataFrame(dict(set=s, OBJECT=d['OBJECT'], ra=d['RA'], dec=d['DECL'],
                                      source_id=sid, how=np.where(sid.notna(), 'DR3Name', ''))))
    return pd.concat(rows, ignore_index=True).drop_duplicates('OBJECT')


def columns(path, names):
    with fits.open(path, memmap=True) as h:
        return {n: np.asarray(h[1].data[n]) for n in names}


def native(a):
    """FITS columns are big-endian; pandas needs native byte order."""
    return a.astype(a.dtype.newbyteorder('=')) if a.dtype.kind in 'iuf' else a


def by_position(st, s4_dir):
    """source_id of stars without a DR3 designation, by position in the S4 cross-match."""
    todo = st['source_id'].isna()
    if not todo.any():
        return st
    from scipy.spatial import cKDTree
    c = columns(Path(s4_dir) / 'XMatch_S4_GXP.fits', ['RA', 'DEC', 'source_id_1'])

    def xyz(ra, dec):
        ra, dec = np.radians(ra), np.radians(dec)
        return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    tree = cKDTree(xyz(c['RA'], c['DEC']))
    d, i = tree.query(xyz(st.loc[todo, 'ra'].to_numpy(), st.loc[todo, 'dec'].to_numpy()))
    ok = d < np.radians(1 / 3600)
    idx = st.index[todo][ok]
    st.loc[idx, 'source_id'] = c['source_id_1'][i[ok]].astype(np.int64)
    st.loc[idx, 'how'] = 'position'
    return st


def labels(st, s4_dir, disk_csv):
    ids = st['source_id'].dropna().astype(np.int64).unique()
    disk = pd.read_csv(disk_csv, usecols=['source_id', 'SUBSURVEY', 'teff_xgboost', 'logg_xgboost',
                                          'mh_xgboost', 'teff_gspphot', 'logg_gspphot', 'mg_gspphot'])
    disk = disk[disk['source_id'].isin(ids)].drop_duplicates('source_id')
    out = st.merge(disk.rename(columns={'teff_xgboost': 'teff_disk', 'logg_xgboost': 'logg_disk',
                                        'mh_xgboost': 'mh_disk', 'mg_gspphot': 'MG_gspphot'}),
                   on='source_id', how='left')
    c = columns(Path(s4_dir) / 'XMatch_S4_GXP.fits',
                ['source_id_1', 'SUBSURVEY', 'teff_xgboost', 'logg_xgboost', 'mh_xgboost'])
    keep = np.isin(c['source_id_1'], ids)
    xm = pd.DataFrame({k: native(v[keep]) for k, v in c.items()}).drop_duplicates('source_id_1')
    xm['SUBSURVEY'] = xm['SUBSURVEY'].astype(str).str.strip()
    xm = xm.rename(columns={'source_id_1': 'source_id', 'SUBSURVEY': 'subsurvey_s4',
                            'teff_xgboost': 'teff_s4xm', 'logg_xgboost': 'logg_s4xm',
                            'mh_xgboost': 'mh_s4xm'})
    out = out.merge(xm, on='source_id', how='left')
    rest = np.setdiff1d(ids, xm['source_id'].to_numpy())
    c = columns(Path(s4_dir) / 'table_1_catwise.fits', ['source_id'])
    keep = np.isin(c['source_id'], rest)
    del c
    a = columns(Path(s4_dir) / 'table_1_catwise.fits',
                ['source_id', 'teff_xgboost', 'logg_xgboost', 'mh_xgboost'])
    an = pd.DataFrame({k: native(v[keep]) for k, v in a.items()}).rename(
        columns={'teff_xgboost': 'teff_andrae', 'logg_xgboost': 'logg_andrae', 'mh_xgboost': 'mh_andrae'})
    out = out.merge(an, on='source_id', how='left')
    # adopted label: the newest table that has the star
    out['label_src'] = ''
    for k in ('andrae', 's4xm', 'disk'):
        has = out[f'teff_{k}'].notna()
        for q in ('teff', 'logg', 'mh'):
            out.loc[has, q] = out.loc[has, f'{q}_{k}']
        out.loc[has, 'label_src'] = k
    return out


def crude(out, crude_csv):
    d = pd.read_csv(crude_csv, dtype={'specid': str, 'set': str})
    d['set'] = d['set'].str.zfill(4)
    n = d.groupby('OBJECT').size().rename('n_spectra')
    d = d.sort_values('snr').drop_duplicates('OBJECT', keep='last')
    held = d['teff_c'].notna()
    d['teff_crude'] = np.where(held, d['teff_c'], d['teff'])
    d['logg_crude'] = np.where(held, d['logg_c'], d['logg'])
    d['feh_crude'] = np.where(held, d['feh_c'], d['feh'])
    d['logg_held'] = held
    cols = ['OBJECT', 'snr', 'teff_crude', 'logg_crude', 'feh_crude', 'logg_held']
    return out.merge(d[cols], on='OBJECT', how='left').merge(n, on='OBJECT', how='left')


def main(phot_dir, s4_dir, disk_csv, crude_csv, out):
    st = by_position(stars(phot_dir), s4_dir)
    df = crude(labels(st, s4_dir, disk_csv), crude_csv)
    df.to_csv(out, index=False, float_format='%.5g')
    df['n_spectra'] = df['n_spectra'].fillna(0)
    lab = df['label_src'] != ''
    print(f'{len(df)} stars, {int(df["n_spectra"].sum())} spectra')
    print(df.groupby('set').apply(lambda g: pd.Series(dict(
        stars=len(g), with_source_id=g['source_id'].notna().sum(),
        labelled=(g['label_src'] != '').sum(),
        spectra=int(g['n_spectra'].sum()),
        spectra_labelled=int(g.loc[g['label_src'] != '', 'n_spectra'].sum()),
        disk=(g['label_src'] == 'disk').sum(), s4xm=(g['label_src'] == 's4xm').sum(),
        andrae=(g['label_src'] == 'andrae').sum(),
        gspphot=g['teff_gspphot'].notna().sum()))).to_string())
    print('source_id by', df['how'].value_counts().to_dict())
    both = df['teff_disk'].notna() & df['teff_s4xm'].notna()
    for q in ('teff', 'logg', 'mh'):
        dd = df.loc[both, f'{q}_disk'] - df.loc[both, f'{q}_s4xm']
        print(f'disk - s4xm {q}: n={both.sum()} median {dd.median():+.3f} sd {dd.std():.3f} '
              f'max|d| {dd.abs().max():.3f}')
    g = lab & df['teff_gspphot'].notna()
    for q, gq in (('teff', 'teff_gspphot'), ('logg', 'logg_gspphot')):
        dd = df.loc[g, gq] - df.loc[g, q]
        print(f'GSP-Phot - XGBoost {q}: n={g.sum()} median {dd.median():+.3f} '
              f'robust sd {1.4826 * (dd - dd.median()).abs().median():.3f}')
    c = lab & df['teff_crude'].notna()
    for s in ('all',) + SETS:
        sel = c & ((df['set'] == s) if s != 'all' else True)
        if not sel.any():
            continue
        parts = []
        for q, cq in (('teff', 'teff_crude'), ('logg', 'logg_crude'), ('mh', 'feh_crude')):
            dd = df.loc[sel, cq] - df.loc[sel, q]
            parts.append(f'{q} {dd.median():+.2f} (robust sd {1.4826 * (dd - dd.median()).abs().median():.2f})')
        print(f'crude - XGBoost, {s} (n={sel.sum()}): ' + ', '.join(parts))


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
