!===============================================================================
!===============================================================================
MODULE solver_mod
  !Dispatch modules
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE download_mod
  USE mpi_mesg_mod
  USE bifrost_mod
  USE mesh_interface_mod
  USE bifrost_rt_mod
  !Bifrost modules
  USE mhd_mod
  USE extras_mod
  USE radiation_mod
  USE eos_mod
  USE index_mod
  implicit none
  private
  type, extends(patch_t), public:: solver_t
    type(bifrost_rt_t):: bifrost_rt
    type(bifrost_mesh_t):: bifrost_mesh
  contains
    procedure:: init
    procedure:: dnload
    procedure:: get
    procedure:: init
    procedure:: pack
    procedure:: put
    procedure:: unpack
    procedure:: update
    ! ----------------------------
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

!=======================================================================
!> Organize calls to the extras and the mhd solver init
!=======================================================================
SUBROUTINE init (self)
  class(solver_t) :: self
  !.....................................................................
  call self%mhd_t%init
  call self%extras_t%init
END SUBROUTINE init

!===============================================================================
!> Initialize the Bifrost patch, including reading in the initial state, which
!> is then copied to self%mem, self%mesh, and self%time
!===============================================================================
SUBROUTINE init (self)
  class(solver_t):: self
  ! --------------------------------------------------------------------
  integer:: nv, nt, ng(3), nw, n(3), l(3), u(3),i
  call trace%begin ('solver_t%init')
  self%nv = 8                                                           ! set number of variables
  self%nt = 5                                                           ! set number of time slots
  self%ng = 5                                                           ! set number of guard zones
  self%nw = 1                                                           ! time derivatives allocated
  ng=5
  call self%idx%init(ienergy=5,self%mhd)                                ! setup index_t indices
  call self%patch_t%setup(nv=8, nt=5, ng=ng,nw=1)                       ! allocate memory and set size            !FIXME: hardcoded defaults for now
  call bifrost%init                                                     ! initialize the Bifrost state
  call self%bifrost_mesh%init(self%mesh)                                ! initialize the Bifrost mesh
  call self%bifrost_mesh%put
  !call self%bifrost_mesh%check
  ! Initialize Bifrost-local modules
  call init_mhd
  call init_eos
  call init_radiation
  call init_snap
  call init_extras
  call init_boundaries
  ! the next one might not be needed in DISPATCH
  call update_mhd
  self%bifrost_mesh%myid = self%id
  do i=1,3
    l(i)=self%mesh(i)%lb
    u(i)=self%mesh(i)%ub
  end do
  call bifrost%get (self%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,self%it,1), self%mesh, self%time, self%dtime)
  call self%bifrost_rt%init
  call self%bifrost_rt%get
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Get mesh details from Bifrost
!===============================================================================
SUBROUTINE get (self)
  class(solver_t):: self
  !---------------------------------------------------------------------
  call trace%begin ('solver_t%get')
  call bifrost%get (self%mem(:,:,:,:,self%new,1), self%mesh, self%time, self%dtime)
  call self%bifrost_rt%get
  call trace%end()
END SUBROUTINE get

!===============================================================================
!> Get mesh details from Bifrost
!===============================================================================
SUBROUTINE put (self)
  class(solver_t):: self
  !---------------------------------------------------------------------
  call trace%begin ('solver_t%put')
  call self%bifrost_mesh%put
  call bifrost%put (self%mem(:,:,:,:,self%it,1), self%mesh, self%time)
  call self%bifrost_rt%put
  call trace%end()
END SUBROUTINE put

!===============================================================================
!> Download guard zone values
!===============================================================================
SUBROUTINE dnload (self, only)
  class(solver_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%dnload')
  call download%download_link (self%link)
  call trace%end()
END SUBROUTINE dnload

!===============================================================================
!> Organize calls to the extras and the hydro solver.  To be compatible with
!> Riemann solvers, which initially copy the %it slot to %new, and expect
!> pre_update() methods to refer to %new, we do the same copy operation here.
!===============================================================================
SUBROUTINE update (self)
  class(solver_t):: self
  !.............................................................................
  call trace%begin ('solver_t%update')
  call self%mem_lock(self%new)%set ('update')
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  call self%mem_lock(self%new)%unset ('update')
  call self%pre_update
  call self%put
  call bifrost%alloc_scratches
  call bifrost%update (self%time, self%dtime)
  call self%get
  call bifrost%dealloc_scratches
  call self%post_update
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Pack information about a patch into a buffer.  Let the patch_t%pack take
!> care of most of it, and then append whatever more is needed
!===============================================================================
SUBROUTINE pack (self, mesg)
  class(solver_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%pack')
  call self%patch_t%pack (mesg)
  ! ...
  call trace%end()
END SUBROUTINE pack

!===============================================================================
!> Unpack an MPI buffer.  The patch_t%unpack takes care of locating the patch,
!> and unpacking the patch related header, after which remains only to unpack
!> whatever else we appended
!===============================================================================
SUBROUTINE unpack (self, mesg)
  class(solver_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('solver_t%unpack')
  call self%patch_t%unpack (mesg)
  ! ...
  call trace%end()
END SUBROUTINE unpack

!===============================================================================
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

!===============================================================================
FUNCTION gas_velocity_vector (self) RESULT (v)
  class (solver_t):: self
  real, dimension(:,:,:), pointer :: d, lnd
  real, dimension(:,:,:,:), pointer:: p, ld, dd
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3):: v
  !-----------------------------------------------------------------------------
  v=0
END FUNCTION gas_velocity_vector

!===============================================================================
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
