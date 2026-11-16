# OpenMULTI3D

OpenMULTI3D synthesises stellar spectra from 1D and 3D model atmospheres, in LTE and NLTE. It solves the statistical equilibrium of model atoms and the radiative transfer on hydrostatic and hydrodynamic model atmospheres (MARCS, Stagger, M3DIS and others), with the task-based parallelisation of the DISPATCH framework.

This repository holds DISPATCH and, in [`experiments/OpenMULTI3D/`](experiments/OpenMULTI3D/), the OpenMULTI3D application:

- **`dispatch.x`**, the Fortran executable built from `SRC/` and `HOP/`, which runs one synthesis or NLTE solve from a namelist.
- **HOP**, the equation of state, continuous and line opacity library, also packaged on its own as `hopacity`.
- **`tumult`**, the Python package that prepares and launches runs, reads their output and provides analysis workflows such as abundance fits and curves of growth.
- **TumultStudio**, a browser GUI (Bokeh) built on `tumult` for interactive spectroscopic analysis.

## Range of stellar parameters

OpenMULTI3D is designed for diagnostic spectroscopy of late-B and AFGKM-type stars: main-sequence stars, subgiants, red giants, AGB stars and red supergiants. It has been tested with 1D (MARCS, MAFAGS-OS, M1DIS) and 3D (M3DIS) model atmospheres in this range:

| Parameter | Range |
|---|---|
| Effective temperature | 2500 K ≲ T<sub>eff</sub> ≲ 8000 K |
| Surface gravity | −1.0 ≲ log g ≲ 5.5 |
| Metallicity | −10 ≲ [Fe/H] ≲ +1 |

The limits in temperature and metallicity follow mainly from the microphysics: the partition functions, the ions and molecules covered by the line lists and the continuum cross-sections. Towards low log g the limit is the plane-parallel geometry of the radiative transfer, which neglects the sphericity of extended giant and supergiant atmospheres. At [Fe/H] = −10 the atmospheres are effectively metal-free apart from Li, so the spectra show only H I and Li I features.

## Citing OpenMULTI3D

If you use OpenMULTI3D, HOP or `hopacity` in published work, please cite the code release paper, [Hoppe et al. (2026), arXiv:2611.09370](https://arxiv.org/abs/2611.09370), accepted for publication in A&A:

```bibtex
@ARTICLE{HoppeBergemann2026,
       author = {{Hoppe}, Richard and {Bergemann}, Maria and {Eitner}, Philipp and {Nordlund}, {\r{A}}ke and {Leenaarts}, Jorrit and {Popovas}, Andrius and {Storm}, Nicholas},
        title = "{The OpenMULTI3D stellar analysis package. A 3D non-LTE spectrum synthesis code for AFGKM-type stars within the DISPATCH framework}",
      journal = {arXiv e-prints},
         year = 2026,
        month = nov,
          eid = {arXiv:2611.09370},
        pages = {arXiv:2611.09370},
          doi = {10.48550/arXiv.2611.09370},
archivePrefix = {arXiv},
       eprint = {2611.09370},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026arXiv261109370H}
}
```

## Getting started

OpenMULTI3D builds with DISPATCH's shared configuration, so clone the whole repository and work in the experiment directory:

```bash
git clone https://github.com/Bergemann-Group-development/dispatch.git
cd dispatch/experiments/OpenMULTI3D
```

Then follow these guides:

- [README_INSTALL.md](experiments/OpenMULTI3D/README_INSTALL.md): requirements, the Python environment, the build, input data downloads and a verification run.
- [QUICKSTART.md](experiments/OpenMULTI3D/QUICKSTART.md): a first LTE and NLTE synthesis, read and plotted in Python and opened in the GUI.
- [docs/manual/](experiments/OpenMULTI3D/docs/manual/): the LaTeX manual, covering the calculations, TumultStudio, the Python workflows, running on a cluster, the physics and the namelist reference.

## DISPATCH

DISPATCH is a task-based framework for astrophysical simulations, primarily 3D radiation-hydrodynamics, by Åke Nordlund, Jon Ramsey, Andrius Popovas, Troels Haugbølle and Michael Kuffmeier. Its sources sit at the top level of this repository and its documentation in [`docs/`](docs/).

## License

OpenMULTI3D is distributed under the GNU General Public License version 3 only ([`experiments/OpenMULTI3D/LICENSE`](experiments/OpenMULTI3D/LICENSE), with a copy in [`LICENSE-GPL`](LICENSE-GPL)); its authors and the licenses of the code and data it builds on are listed in [`experiments/OpenMULTI3D/AUTHORS`](experiments/OpenMULTI3D/AUTHORS). DISPATCH is distributed under the BSD-3-Clause license in [`LICENSE`](LICENSE).
