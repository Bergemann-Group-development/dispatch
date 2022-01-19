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
    logical, dimension(:), pointer, contiguous :: topbot
    procedure(mirror), pointer     :: boundary
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: mirror
    procedure:: mirror_bot_zero_top
    procedure:: radiation
  end type
  integer, save:: verbose=1
CONTAINS

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

  allocate(self%topbot(2))
  self%topbot = .false. ! default value

  select case (condition)
  case(1)
    select case(patch%kind)
      case ('angle_rad')
        self%boundary => mirror_bot_zero_top
      case default
        self%boundary => mirror 
    end select
  case(-1)
    self%boundary => mirror
  end select

  call self%standard_boundaries_t%init ()
END SUBROUTINE init

!===============================================================================
SUBROUTINE mirror(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer :: i, b

    if (self%topbot(1)) then
        b = ubound(mem,3)
        mem(:,:,b,:,:) = mem(:,:,b-1,:,:)
    endif

    if (self%topbot(2)) then
        b = lbound(mem,3)
        mem(:,:,b,:,:) = mem(:,:,b+1,:,:)
    endif
END SUBROUTINE

!===============================================================================
SUBROUTINE mirror_bot_zero_top(self, mem)
  class(boundaries_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer :: i, b

    if (self%topbot(1)) then
        b = ubound(mem,3)
        mem(:,:,b,:,:) = 0
    endif

    if (self%topbot(2)) then
        b = lbound(mem,3)
        mem(:,:,b,:,:) = mem(:,:,b+1,:,:)
    endif
END SUBROUTINE


SUBROUTINE radiation (self, mem, src, rk, tt, dir, mu)
  class(boundaries_t) :: self
  real, dimension(:,:,:,:), pointer:: mem, src, rk
  real, dimension(:,:,:), pointer:: tt,mem2
  real   :: mu(3)
  integer:: dir
  !.....................................................................
  ! keep this just for rt_mod to compile
END SUBROUTINE radiation

END MODULE boundaries_mod
