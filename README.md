# OpenMULTI3D

OpenMULTI3D synthesises stellar spectra from 1D and 3D model atmospheres, in LTE and NLTE. It solves the statistical equilibrium of model atoms and the radiative transfer on hydrostatic and hydrodynamic model atmospheres (MARCS, Stagger, M3DIS and others), with the task-based parallelisation of the DISPATCH framework.

This repository holds DISPATCH and, in [`experiments/OpenMULTI3D/`](experiments/OpenMULTI3D/), the OpenMULTI3D application:

- **`dispatch.x`**, the Fortran executable built from `SRC/` and `HOP/`, which runs one synthesis or NLTE solve from a namelist.
- **HOP**, the equation of state, continuous and line opacity library, also packaged on its own as `hopacity`.
- **`tumult`**, the Python package that prepares and launches runs, reads their output and provides analysis workflows such as abundance fits and curves of growth.
- **TumultStudio**, a browser GUI (Bokeh) built on `tumult` for interactive spectroscopic analysis.

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

OpenMULTI3D is distributed under the GNU General Public License version 3 only ([`experiments/OpenMULTI3D/LICENSE`](experiments/OpenMULTI3D/LICENSE)); its authors and the licenses of the code and data it builds on are listed in [`experiments/OpenMULTI3D/AUTHORS`](experiments/OpenMULTI3D/AUTHORS). DISPATCH is distributed under the BSD-3-Clause license in [`LICENSE`](LICENSE).
