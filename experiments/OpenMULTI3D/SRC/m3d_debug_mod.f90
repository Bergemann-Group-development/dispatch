!===============================================================================
!> Debug probe indices and diagnostic switches, read from &m3d_debug_params by
!> m3d_patch_mod.
!>
!> pcheck, fcheck, linecheck, lvlcheck and atomcheck pick the cell, frequency,
!> transition, level and atom whose values the patch holding that cell prints.
!> check_NaN adds NaN checks to the opacity, population and intensity arrays;
!> diag_pops, diag_db and diag_no_scat_loop turn on single diagnostics of the
!> statistical-equilibrium and scattering iterations.
!===============================================================================
module m3d_debug_mod
  implicit none

  real(8), public, dimension(6) :: memcheck !< values at the probe, gathered for printing
  integer, public :: pcheck(3)=[1,1,-1]     !< probe cell (x, y, z)
  integer, public :: fcheck=1               !< probe frequency index
  integer, public :: linecheck=1            !< probe radiative transition index
  integer, public :: lvlcheck=1             !< probe level index
  integer, public :: atomcheck=1            !< probe atom index
  logical, public :: check_NaN=.false.      !< check opacities, populations and intensities for NaN
  !> diagnose the statistical-equilibrium linear solve: per patch and iteration,
  !> report how many level corrections hit the exponential damping floor and how
  !> large the residual of the LU solve is at the worst cell
  logical, public :: diag_pops=.false.
  !> one-shot detailed-balance dump of the radiative rates at pcheck
  logical, public :: diag_db=.false.
  !> Cut the continuum-scattering Lambda-loop, without touching any opacity.
  !> eta = chi_tot*Bnu + sigma_scat*(J_old - Bnu) becomes eta = chi_tot*Bnu, i.e.
  !> scattering is treated as pure absorption.  chi (which already contains
  !> sigma_scat) is untouched, so every optical depth, every Psi and every formal
  !> solution keeps its geometry; the only thing removed is the dependence of
  !> the source function on the previous iteration's mean intensity.  A model
  !> change, not a fix: use it to decide whether an instability is fed by the
  !> unaccelerated scattering iteration.
  logical, public :: diag_no_scat_loop=.false.
  !> global z of the patch currently being updated, so the diag_pops prints and the
  !> statistical-equilibrium error (which happen down in atom_mod, where only
  !> patch-local indices exist) can be placed in the atmosphere; set by
  !> nlte_patch_t%update before update_populations
  integer, public :: diag_zs=0
  !$omp threadprivate(memcheck, diag_zs)

contains

end module m3d_debug_mod
