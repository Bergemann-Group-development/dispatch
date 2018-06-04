!> =============================================================================
!> Data type for computin global average, while using only nearest node nbor
!> communication
!> =============================================================================
MODULE global_mod
  USE trace_mod
  USE lagrange_mod
  USE mpi_mod
  implicit none
  private
  type, public:: global_t
    integer:: nMPI, nt
    real(8):: cadence
    real(8), dimension(:,:), allocatable:: v, w, t
  contains
    procedure:: init
    procedure:: update
    procedure:: merge
    procedure:: average
  end type
CONTAINS

!> =============================================================================
!> Initialize -- one instance is needed for each global average
!> =============================================================================
SUBROUTINE init (self, nMPI, nt, cadence)
  class(global_t):: self
  integer:: nMPI, nt
  real(8):: cadence
  !.............................................................................
  call trace%begin ('global_t%init')
  self%cadence = cadence
  self%nMPI = nMPI
  self%nt = nt
  allocate (self%v(nt,nMPI))
  allocate (self%w(nt,nMPI))
  allocate (self%t(nt,nMPI))
  self%v = 0d0
  self%w = 0d0
  self%t = 0d0
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Keep updating the first time slot until it ís more than "cadence" ahead, 
!> then rotate.
!> =============================================================================
SUBROUTINE update (self, v, w, t)
  class(global_t):: self
  real(8):: v, w, t
  !.............................................................................
  call trace%begin ('global_t%update')
  if (t > self%t(2,mpi%rank)+self%cadence) then
    self%v(2:self%nt,mpi%rank) = self%v(1:self%nt-1,mpi%rank)
    self%w(2:self%nt,mpi%rank) = self%w(1:self%nt-1,mpi%rank)
    self%t(2:self%nt,mpi%rank) = self%t(1:self%nt-1,mpi%rank)
  end if
  self%v(self%nt,mpi%rank) = v
  self%w(self%nt,mpi%rank) = w
  self%t(self%nt,mpi%rank) = t
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> Merge information with that from a node neighbor
!> =============================================================================
SUBROUTINE merge (self, nbor)
  class(global_t):: self, nbor
  !.............................................................................
  integer:: i, j
  !-----------------------------------------------------------------------------
  call trace%begin ('global_t%merge')
  do j=1,self%nMPI
  do i=1,self%nt
    if (nbor%t(i,j) > self%t(i,j)) then
      self%v(i,j) = nbor%v(i,j)
      self%w(i,j) = nbor%w(i,j)
      self%t(i,j) = nbor%t(i,j)
    end if
  end do
  end do
  call trace%end
END SUBROUTINE merge

!> =============================================================================
!> Return a global average
!> =============================================================================
FUNCTION average (self, t)
  class(global_t):: self
  real(8):: t, average
  !.............................................................................
  real(8):: v, w, weight
  integer:: i, j
  !-----------------------------------------------------------------------------
  call trace%begin ('global_t%average')
  weight = 0d0
  average = 0d0
  do j=1,self%nMPI
    v = lagrange%interpolation(t, self%t(:,i), self%v(:,i))
    w = lagrange%interpolation(t, self%t(:,i), self%w(:,i))
    weight = weight + w
    average = average + v*w
  end do
  average = average
  call trace%end
END FUNCTION average

END MODULE global_mod
