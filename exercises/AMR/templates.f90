!> =============================================================================
!> The patch list module deals only with patch, which can be updated, and does not
!> need to know any internal details about patches, which just need to respond to
!> a call patch%update.  The refine%needed(patch) returns .true. if the patch needs
!> refinement, in which case the patch_list updater will assume that new patches
!> are injected on the patch list
!> =============================================================================
MODULE patch_list_mod
  USE timestep_mod
  USE refine_mod
  type:: patch_list_t
  end type
CONTAINS
SUBROUTINE update (self)
  type (patch_list):: self
  type (timestep_t), pointer:: patch
  ...
  patch => link%patch
  if (patch%refine(patch)) then
    call split (link)
  else
    call patch%update
  end
  ...
END
END MODULE patch_list_mod

!> =============================================================================
!> adams_bashforth.f90
!> The timestep module is an (optional) extension to the MHD module, which can
!> take higher order time steps, e.g. with the Adams-Bashforth methods.  In the
!> case where a Riemann solver is used, the entire thing may be a module that
!> defines a timestep_t object, which responds to a call patch%update.  However,
!> it must always save a number of previous time slices, even if they might not
!> be used in the update.  It may thus be sensible to put the time slice rotation
!> in an underlying service layer; i.e., in the generic patch_t, which is overlaid
!> by the mhd_t and timestep_t layers.
!> =============================================================================
MODULE timestep_mod
  type, extends(mhd_t):: timestep_t
  ...
  end type
CONTAINS
!> =============================================================================
!> The time step update procedure should also take care of the periodic wrapping
!> of time indices to memory slots 1:ntime.  The 'it' slot should always be the
!> most 'current' complete time slot, while the 'new' time slot is the one that
!> receives the updated values.  After a time step is complete it becomes the new
!> 'it' time slot, in the periodic mapping.
!> =============================================================================
SUBROUTINE update (self)
  type(timestep_t):: self
  ...
  call self%pde
  ...
END SUBROUTINE
END MODULE timestep_mod

!> =============================================================================
!> file: primitive_mod.f90
!> The MHD module should already contain an update procedure, which just takes a
!> 1st order time step.  The MHD patch should also contain a Courant procedure,
!> which uses detailed knowledge of the MHD solver to evaluate the allowed time
!> step size.
!> =============================================================================
MODULE mhd_mod
  type, extends(patch_t):: mhd_t
    type(scalar_t), pointer:: Rho, Pg, Pm, ...
    type(vector_t), pointer:: U, B, ...
  contains
    procedure:: init
    procedure:: pde
  end type
CONTAINS
!> =============================================================================
!> The PDE init routine maps the variable components to memory slot, so they
!> can be enumerated in the timestop routine, while referred to by variably
!> names in the PDE routine
!> =============================================================================
SUBROUTINE init (self)
  p%x => patch%mem(:,:,:,4,patch%it,1)
  B%x => patch%mem(:,:,:,7,patch%it,1)
  ...
END
!> =============================================================================
!> The PDE procedure should use the MHD variables, in scalar and vector form, to
!> compute time derivative value (mapped into patch%mem by the mhd_t%init procedure)
!> =============================================================================
SUBROUTINE pde (self)
  type(mhd_t):: self
  ...
  dBdt = -curl(E)
END
MODULE mhd_mod
