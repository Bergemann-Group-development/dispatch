module m3d_debug_mod
  implicit none

  real(8), public, dimension(6) :: memcheck
  integer, public :: pcheck(3)=[1,1,-1]
  integer, public :: fcheck=1, linecheck=1, lvlcheck=1, atomcheck=1
  logical, public :: check_NaN=.false.
  ! diagnose the statistical-equilibrium linear solve: per patch and iteration,
  ! report how many level corrections hit the exponential damping floor and how
  ! large the residual of the LU solve is at the worst cell
  logical, public :: diag_pops=.false.
  ! one-shot detailed-balance dump of the radiative rates at pcheck
  logical, public :: diag_db=.false.
  ! Cut the continuum-scattering Lambda-loop, without touching any opacity.
  ! eta = chi_tot*Bnu + sigma_scat*(J_old - Bnu) becomes eta = chi_tot*Bnu, i.e.
  ! scattering is treated as pure absorption.  chi (which already contains
  ! sigma_scat) is untouched, so every optical depth, every Psi and every formal
  ! solution keeps its geometry -- the ONLY thing removed is the dependence of
  ! the source function on the previous iteration's mean intensity.  A model
  ! change, not a fix: use it to decide whether an instability is fed by the
  ! unaccelerated scattering iteration.
  logical, public :: diag_no_scat_loop=.false.
  ! global z of the patch currently being updated, so the diag_pops prints (which
  ! happen down in atom_mod, where only patch-local indices exist) can be placed
  ! in the atmosphere; set by nlte_patch_t%update before update_populations
  integer, public :: diag_zs=0
  !$omp threadprivate(memcheck, diag_zs)

contains

end module m3d_debug_mod
