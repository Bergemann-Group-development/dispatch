!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE patch_mod
  USE units_mod
  ! ----------------------------------------------------------------------------
  ! Experiment - specific USE modules
  ! ----------------------------------------------------------------------------
  !USE lower_damp_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    logical :: lower=.false.
    logical :: upper=.false.
    procedure(mirror), pointer     :: boundary
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: pre_init
    procedure:: init
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: mirror
    procedure:: no_boundary
  end type
  integer, save:: verbose=1
CONTAINS

!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(boundaries_t) :: self
  class(patch_t), target :: patch
  ! ----------------------------------------------------------------------------
  call trace%begin('boundaries_t%pre_init')
  call self%standard_boundaries_t%pre_init(patch)

  self%lower = self%patch%origin(3) <= 0
  self%upper = self%patch%origin(3)+self%patch%size(3) >= self%patch%box(3)

  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the boundary data type.
!===============================================================================
SUBROUTINE init (self)
  class(boundaries_t) :: self
  ! ----------------------------------------------------------------------------
  class(patch_t), pointer:: patch
  integer:: iostat, condition=-1
  logical, save:: first_time=.true.
  namelist /boundary_params/ verbose, condition
  !-----------------------------------------------------------------------------
  call trace%begin('boundaries_t%init')
  !.............................................................................
  call io%lock%set('boundaries_t%init')
  patch => self%patch
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, boundary_params)
    end if
    call io%namelist_warning ('boundary_params', iostat)
    write (io%namelist, boundary_params)
    first_time = .false.
  end if
  call io%lock%unset('boundaries_t%init')
  if (.not. io%master) verbose=0

  select case(patch%kind)
    case ('hlld_mhd_patch')
      self%boundary => mirror
    case default
      self%boundary => no_boundary
  end select

  call self%standard_boundaries_t%init ()
  call trace%end ()
END SUBROUTINE init

!===============================================================================
SUBROUTINE mirror(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer :: i, b

    if (self%upper) then
        b = ubound(mem,3)
        mem(:,:,b,:,:) = mem(:,:,b-1,:,:)
    endif

    if (self%lower) then
        b = lbound(mem,3)
        mem(:,:,b,:,:) = mem(:,:,b+1,:,:)
    endif
END SUBROUTINE

!===============================================================================
SUBROUTINE no_boundary(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  ! do nothing
END SUBROUTINE

END MODULE boundaries_mod
