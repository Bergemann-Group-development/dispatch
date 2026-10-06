"""Effective temperatures of the <3D> averaging study's flux runs.

    python avg3d_study_teff.py LAM_FILE RUN [RUN ...] > teff.txt

LAM_FILE is the vacuum wavelength grid the runs were computed on (flx_wavelengths.vac);
RUN are run folders (data/avg3d_study/{1d,3d}/<model>/<name>).  Teff = (F/sigma)^(1/4)
with F the trapezoid integral of the LTE physical flux over the vacuum grid (130 nm to
20 micron), sorted as the runs sort it.  The flux beyond 20 micron is left out; for a
4500-5800 K black body it is 1-2e-4 of the total, 0.1-0.3 K in Teff.
"""
import sys
from pathlib import Path

import astropy.units as u
import numpy as np

import tumult
from tumult.utils import calc_teff


def teff(run_folder, lam_vac):
    run = tumult.read(str(run_folder), verbose=0)
    flux = np.asarray(run.flux_lte, float)
    if flux.size != lam_vac.size:
        raise ValueError(f'{run_folder}: {flux.size} fluxes for {lam_vac.size} wavelengths')
    return calc_teff(lam_vac * u.AA, flux * u.erg / u.s / u.cm**2 / u.AA).to_value(u.K)


def main(lam_file, runs):
    lam_vac = np.sort(np.loadtxt(lam_file))
    print('# model run Teff[K]')
    for r in map(Path, runs):
        print(f'{r.parent.name} {r.name} {teff(r, lam_vac):.2f}', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
