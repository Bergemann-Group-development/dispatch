# OpenMULTI3D source reference

OpenMULTI3D synthesises stellar spectra from 1D and 3D model atmospheres, in LTE and NLTE. This reference is generated from the doc comments of its Fortran and C sources in `SRC/` and `HOP/`. Installation, running and the Python interface (`tumult` and the TumultStudio GUI) are described in the [installation guide](https://github.com/Bergemann-Group-development/dispatch/blob/main/experiments/OpenMULTI3D/README.md) and the [manual](https://github.com/Bergemann-Group-development/dispatch/tree/main/experiments/OpenMULTI3D/docs/manual) of the [repository](https://github.com/Bergemann-Group-development/dispatch/tree/main/experiments/OpenMULTI3D).

The **Modules** and **Data Types** menus list every module and derived type; **Files** follows the source tree.

## SRC: the application

`SRC/` builds `dispatch.x` on the DISPATCH framework, whose task list schedules the patches below. The main program is \ref dispatch.f90.

- \ref m3d_patch_mod holds `m3d_patch_t`, the base type of all OpenMULTI3D tasks, and the run parameters.
- \ref read_external_mod and \ref read_atmos_mod read the model atmosphere and build the task list; \ref experiment_mod holds the depth patches.
- \ref composition_patch_mod evaluates the equation of state and the opacity tables of a depth slab.
- \ref nlte_patch_mod solves the statistical equilibrium of the model atoms of \ref atom_mod, accelerated by \ref ng_mod.
- \ref rt_patch_mod computes the opacities, formal solution, mean intensity and rates of one frequency window, direction and depth slab; \ref ray_patch_mod, \ref ray_mod, \ref rt_column_mod and \ref rt_params_mod supply windows, quadratures, sweeps and settings.

## HOP: equation of state and opacities

`HOP/` computes the equation of state and the continuous and line opacities. It also builds on its own as `libhop` and the Python package [hopacity](https://github.com/Bergemann-Group-development/hopacity), whose C interface is declared in `HOP/hopacity/libhop.h`.

- Equation of state: \ref composition_mod, \ref element_mod, \ref ion_mod, \ref molecule_mod.
- Continuous opacity: \ref bftrans_mod, \ref osmet_mod, \ref optab_mod.
- Line opacity: \ref linelist_mod, \ref line_types_mod, \ref linelist_reader_mod, the Voigt approximations of \ref qlg, the hydrogen profiles of \ref hprofl_mod and \ref hvcs_mod, and the hydrogen-collision cross-sections of \ref anstee_mod.
- Model atoms and atmospheres: \ref trans_mod, \ref bbtrans_mod, \ref level_mod, \ref atom_convert_mod, \ref atmos_reader_mod.

## Not covered

The formal solvers (short and long characteristics; linear, quadratic and Hermite) are implemented in the submodule `rt_solvers_smod` of `SRC/rt_solvers_mod.f90`, which Doxygen cannot parse; their interfaces are listed under \ref rt_patch_mod. The DISPATCH framework itself (`gpatch_t`, the task list, the MPI and OpenMP scheduling) lives outside `experiments/OpenMULTI3D` and is not part of this reference.
