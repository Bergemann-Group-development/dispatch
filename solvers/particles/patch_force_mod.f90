!===============================================================================
!> Data type for storing the time slices of forces from a patch onto a sink 
!> particle
!===============================================================================
MODULE patch_force_mod
  USE dll_mod
  USE trace_mod
  implicit none
  private
  type, public, extends(dll_node_t):: patch_force_t
    real, allocatable:: f(:,:)
    real(8), allocatable:: t(:)
    integer:: it=0, nt=5
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: update
    procedure:: force
  end type
CONTAINS

!===============================================================================
!> Initialize storage for one  patch
!===============================================================================
SUBROUTINE init (self)
  class(patch_force_t):: self
  !.............................................................................
  call trace%begin ('patch_force_t%init')
  allocate (self%f(self%nt,3), self%t(self%nt))
  self%it = 1
  self%t = 0d0
  self%f = 0.0
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Remove storage for one  patch
!===============================================================================
SUBROUTINE dealloc (self)
  class(patch_force_t):: self
  !.............................................................................
  deallocate (self%f, self%t)
END SUBROUTINE dealloc

!===============================================================================
!> Update storage for one  patch
!===============================================================================
SUBROUTINE update (self, force, time)
  class(patch_force_t):: self
  real:: force(3)
  real(8):: time
  !.............................................................................
  self%it = mod(self%it,self%nt) + 1
  self%f(self%it,:) = force
  self%t(self%it) = time
END SUBROUTINE update

!===============================================================================
!> Interpolate stored forces to the given time
!===============================================================================
FUNCTION force (self, time)
  class(patch_force_t):: self
  real:: force(3)
  real(8):: time
  !.............................................................................
  real(8):: t(self%nt), p
  integer:: i, j(self%nt), k
  !-----------------------------------------------------------------------------
  ! Construct non-decreasing sequence of stored times
  !-----------------------------------------------------------------------------
  do i=1,self%nt
    j(i) = mod(self%it+i-1,self%nt) + 1
  end do
  t = self%t(j)
  !-----------------------------------------------------------------------------
  ! Search for interpolation interval (k,k+1). with k in 1...nt-1
  !-----------------------------------------------------------------------------
  k = 1
  do i=1,self%nt-1
    if (t(i) < time) then
      k = i
    else
      exit
    end if
  end do
  p = (time-t(k))/max(t(k+1)-t(k),1d-30)
  p = max(min(p,1d0),0d0)
  force = (1d0-p)*self%f(j(k),:) + p*self%f(j(k+1),:)
END FUNCTION force

END MODULE patch_force_mod
