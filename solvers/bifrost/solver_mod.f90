!> =============================================================================
!> =============================================================================
MODULE solver_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE download_mod
  USE mpi_mesg_mod
  USE bifrost_mod
  USE mesh_interface_mod
  USE bifrost_rt_mod
  implicit none
  private
  type, extends(patch_t), public:: solver_t
    type(bifrost_rt_t):: bifrost_rt
    type(bifrost_mesh_t):: bifrost_mesh
  contains
    procedure:: init
    procedure:: get
    procedure:: put
    procedure:: dnload
    procedure:: update
    procedure:: pack
    procedure:: unpack
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: gas_pressure
    procedure:: gas_velocity_vector
    procedure:: gas_velocity_scalar
    generic:: gas_velocity => gas_velocity_vector, gas_velocity_scalar
  end type
  type(solver_t), public:: solver
  character(len=80), save:: id = &
  '$Id$'
CONTAINS

!> =============================================================================
!> Initialize the Bifrost patch, including reading in the initial state, which
!> is then copied to self%mem, self%mesh, and self%time
!> =============================================================================
SUBROUTINE init (self)
  class(solver_t):: self
  ! --------------------------------------------------------------------
  integer :: n(3)
  call trace%begin ('solver_t%init')
  call bifrost%init (n,self%bifrost_mesh)                                                 ! initialize the Bifrost state
  self%nv = 8                                   ! set number of variables
  self%nt = 7                                   ! set number of time slots
  self%ng = 10                                  ! set number of guard zones
  self%nw = 1                                   ! time derivatives allocated
  self%n = n
  call self%patch_t%init                        ! allocate memory and set size
  associate (l => self%li, u => self%ui)
  call bifrost%get (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,self%it,1), self%bifrost_mesh, self%time, self%dtime)
  call self%bifrost_rt%init
  call self%bifrost_rt%get
  end associate
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Get mesh details from Bifrost
!> =============================================================================
SUBROUTINE get (self)
  class(solver_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%get')
  associate (l=>self%li, u=> self%ui)
  call bifrost%get (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,self%new,1), self%bifrost_mesh, self%time, self%dtime)
  call self%bifrost_rt%get
  end associate
  call trace%end
END SUBROUTINE get

!> =============================================================================
!> Get mesh details from Bifrost
!> =============================================================================
SUBROUTINE put (self)
  class(solver_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%get')
  associate (l=>self%li, u=> self%ui)
  call bifrost%put (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,self%it,1), self%bifrost_mesh, self%time)
  call self%bifrost_rt%put
  end associate
  call trace%end
END SUBROUTINE put

!> =============================================================================
!> Download guard zone values
!> =============================================================================
SUBROUTINE dnload (self, only)
  class(solver_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%dnload')
  call download%download_link (self%link)
  call trace%end
END SUBROUTINE dnload

!> =============================================================================
!> Take one time step, by copying relevant data to BIFROST, calling its update
!> procedure, and then copying the patch data back from Bifrost.
!> =============================================================================
SUBROUTINE update (self)
  class(solver_t):: self
  !.............................................................................
  call trace%begin ('solver_t%update')
  call self%put
  call bifrost%update (self%time, self%dtime)
  call self%get
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Pack information about a patch into a buffer.  Let the patch_t%pack take
!> care of most of it, and then append whatever more is needed
!> =============================================================================
SUBROUTINE pack (self, mesg)
  class(solver_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%pack')
  call self%patch_t%pack (mesg)
  ! ...
  call trace%end
END SUBROUTINE pack

!> =============================================================================
!> Unpack an MPI buffer.  The patch_t%unpack takes care of locating the patch,
!> and unpacking the patch related header, after which remains only to unpack
!> whatever else we appended
!> =============================================================================
SUBROUTINE unpack (self, mesg)
  class(solver_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%unpack')
  call self%patch_t%unpack (mesg)
  ! ...
  call trace%end
END SUBROUTINE unpack

!> =============================================================================
SUBROUTINE void (self, out)
  class (solver_t):: self
  real, dimension(:,:,:):: out
END SUBROUTINE

FUNCTION gas_pressure (self) RESULT (pg)
  class(solver_t):: self
  real, dimension(:,:,:), pointer:: d, s
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  integer, save:: nprint=3
  pg =0
END FUNCTION gas_pressure

!> =============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer :: d, lnd
  real, dimension(:,:,:,:), pointer:: p, ld, dd
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
  v=0
END FUNCTION gas_velocity_vector

!> =============================================================================
FUNCTION gas_velocity_scalar (self, idir) RESULT (v)
  class (solver_t):: self
  integer:: idir
  real, dimension(:,:,:), pointer :: d, lnd, ld, dd
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: v
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  v=0
  nullify (d)
END FUNCTION gas_velocity_scalar

END MODULE solver_mod
