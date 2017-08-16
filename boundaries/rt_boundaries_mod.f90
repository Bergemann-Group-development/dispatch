!> =============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!> =============================================================================
MODULE rt_boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mpi_mod
  USE bits_mod
  USE mesh_mod
  USE eos_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  !-----------------------------------------------------------------------
  type, public :: rt_outer_bndry_t
    integer:: id = 0
    integer:: bits = 0
    logical:: unset = .true.
    real, dimension(:,:,:), pointer :: xl => null()
    real, dimension(:,:,:), pointer :: xu => null()
    real, dimension(:,:,:), pointer :: yl => null()
    real, dimension(:,:,:), pointer :: yu => null()
    real, dimension(:,:,:), pointer :: zl => null()
    real, dimension(:,:,:), pointer :: zu => null()
  contains
    procedure:: init => init_outer_bndry
    procedure:: set
    procedure:: is_set
    procedure:: is_clear
    procedure:: condition                                               ! modify outer BCs if necessary
  end type
!-----------------------------------------------------------------------

  type :: rt_inner_bndry_t
    real(8) :: position(3)
    real(8) :: velocity(3)
  end type
!-----------------------------------------------------------------------
  type, extends(rt_inner_bndry_t) :: rt_inner_sphere_t
    real(8) :: radius
  end type
!-----------------------------------------------------------------------
  type, extends(rt_inner_bndry_t) :: rt_inner_cuboid_t
    real(8) :: llc(3), urc(3)
  end type
!-----------------------------------------------------------------------
  type, public:: rt_boundaries_t
    class(rt_inner_bndry_t), dimension(:), pointer :: inner_bndry => null()
    real, dimension(:), pointer    :: xl, xu, yl, yu, zl, zu
    real(8) :: llc(3), urc(3)
  contains
    procedure:: init
    procedure:: update
  end type
  type(rt_boundaries_t), public :: rt_boundaries
CONTAINS

!> =====================================================================
!> Initialize the boundary data type
!> =====================================================================
SUBROUTINE init (self, llc, urc)
  class(rt_boundaries_t)              :: self
  real(8):: llc(3), urc(3)
  !.....................................................................
  real, save :: hot_temp = 1000, cool_temp = 10
  logical, save:: first_time=.true.
  real, dimension(:,:,:), pointer :: lnd, ee
  real, dimension(:,:,:,:), pointer :: src
  namelist /rt_bc_params/ hot_temp, cool_temp
  !---------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_bc_params)
    if (io%master) write (*, rt_bc_params)
    call scaling%init()
  allocate(lnd(1,1,2))
  allocate(ee(1,1,2))
  allocate(src(1,1,2, eos%n_lambda))
  allocate(self%xl(eos%n_lambda))
  allocate(self%xu(eos%n_lambda))
  allocate(self%yl(eos%n_lambda))
  allocate(self%yu(eos%n_lambda))
  allocate(self%zl(eos%n_lambda))
  allocate(self%zu(eos%n_lambda))
  lnd = 0.
  ee(:,:,1) = (cool_temp/scaling%temp)/(eos%gamma - 1.)
  ee(:,:,2) = (hot_temp/scaling%temp)/(eos%gamma - 1.)
  call eos%lookup (shape(ee), x=ee, lny=lnd, &                          ! call the table lookup, ...
                    src=src)  
  self%xl = src(1,1,2,:) - src(1,1,1,:)
  self%xu = 0.0
  self%yl = 0.0
  self%yu = 0.0
  self%zl = 0.0
  self%zu = 0.0
  deallocate(lnd,ee,src)
  self%llc = llc
  self%urc = urc
  end if
  !$omp end critical (input_cr)
  print *, 'cool', self%xu
  print *, 'hot', self%xl
  print *, 'rt boundary llc, urc', self%llc, self%urc
  call trace_end
END SUBROUTINE init

SUBROUTINE update(self)
  class(rt_boundaries_t) :: self
  !---------------------------------------------------------------------
END SUBROUTINE update

!> =====================================================================
!> Initialize the outer boundary shells for rt patches
!> =====================================================================
SUBROUTINE init_outer_bndry (self, n)
  class(rt_outer_bndry_t) :: self
  integer :: n(3)
  !.....................................................................
  integer :: i
  !---------------------------------------------------------------------
  if (self%unset) then
    if (self%is_set(bits%xl)) then
      allocate(self%xl(n(2),n(3), eos%n_lambda))
      do i=1, eos%n_lambda
        self%xl(:,:,i) = rt_boundaries%xl(i)
      end do
    end if
    if (self%is_set(bits%xu)) then
      allocate(self%xu(n(2),n(3), eos%n_lambda))
      do i=1, eos%n_lambda
        self%xu(:,:,i) = rt_boundaries%xu(i)
      end do
    end if
    if (self%is_set(bits%yl)) then
      allocate(self%yl(n(3),n(1), eos%n_lambda))
      do i=1, eos%n_lambda
        self%yl(:,:,i) = rt_boundaries%yl(i)
      end do
    end if
    if (self%is_set(bits%yu)) then
      allocate(self%yu(n(3),n(1), eos%n_lambda))
      do i=1, eos%n_lambda
        self%yu(:,:,i) = rt_boundaries%yu(i)
      end do
    end if
    if (self%is_set(bits%zl)) then
      allocate(self%zl(n(1),n(2), eos%n_lambda))
      do i=1, eos%n_lambda
        self%zl(:,:,i) = rt_boundaries%zl(i)
      end do
    end if
    if (self%is_set(bits%zu)) then
      allocate(self%zu(n(1),n(2), eos%n_lambda))
      do i=1, eos%n_lambda
        self%zu(:,:,i) = rt_boundaries%zu(i)
      end do
    end if
    self%unset = .false.
  end if
END SUBROUTINE init_outer_bndry

!> =============================================================================
!> Set bits in the status word
!> =============================================================================
SUBROUTINE set (self, bits)
  class(rt_outer_bndry_t) :: self
  integer                :: bits
  !-----------------------------------------------------------------------------
  self%bits = ior (self%bits, bits)
END SUBROUTINE set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
LOGICAL FUNCTION is_set (self, bits)
  class(rt_outer_bndry_t) :: self
  integer               :: bits
  !-----------------------------------------------------------------------------
  is_set = (iand (self%bits, bits) /=0)
END FUNCTION is_set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bits)
  class(rt_outer_bndry_t) :: self
  integer               :: bits
  !-----------------------------------------------------------------------------
  is_clear = (iand (self%bits, bits) ==0)
END FUNCTION is_clear

!> =============================================================================
!> Apply a boundary condition, using the procedure pointer "operator"
!> =============================================================================
SUBROUTINE condition (self, time)
  class(rt_outer_bndry_t) :: self
  real(8) :: time
  !.............................................................................

END SUBROUTINE condition

END MODULE rt_boundaries_mod
