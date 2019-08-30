!===============================================================================
!> Data type for computing global average, while using only nearest node nbor
!> communication
!===============================================================================
MODULE global_mod
  USE io_mod
  USE trace_mod
  USE lagrange_mod
  USE mpi_mod
  implicit none
  private
  type, public:: global_t
    integer:: nMPI, nt
    integer:: order=2
    integer:: istep=0
    real(8):: cadence
    real(8), dimension(:,:), allocatable:: v, w, t
    character(len=80):: label
  contains
    procedure:: init
    procedure:: update
    procedure:: merge
    procedure:: average
  end type
CONTAINS

!===============================================================================
!> Initialize -- one instance is needed for each global average
!===============================================================================
SUBROUTINE init (self, label, nMPI, nt, cadence)
  class(global_t):: self
  character(len=*):: label
  integer:: nMPI, nt
  real:: cadence
  !.............................................................................
  call trace%begin ('global_t%init')
  self%label = label
  self%cadence = cadence
  self%nMPI = nMPI
  self%nt = nt
  allocate (self%v(nt,nMPI))
  allocate (self%w(nt,nMPI))
  allocate (self%t(nt,nMPI))
  self%v = 0d0
  self%w = 0d0
  self%t = 0d0
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Keep updating the first time slot until it ís more than "cadence" ahead, 
!> then rotate.
!===============================================================================
SUBROUTINE update (self, v, w, t)
  class(global_t):: self
  real(8):: v, w, t, t1
  !.............................................................................
  call trace%begin ('global_t%update')
  if (t > self%t(self%nt-1,mpi%rank)+self%cadence) then
    self%v(1:self%nt-1,mpi%rank) = self%v(2:self%nt,mpi%rank)
    self%w(1:self%nt-1,mpi%rank) = self%w(2:self%nt,mpi%rank)
    self%t(1:self%nt-1,mpi%rank) = self%t(2:self%nt,mpi%rank)
    self%v(self%nt,mpi%rank) = v
    self%w(self%nt,mpi%rank) = w
    self%t(self%nt,mpi%rank) = t
    t1 = (t+self%t(self%nt-1,mpi%rank))/2d0
    self%istep = self%istep + 1
    if (io%master) then
      print '(a,9f12.6)', 'new times:', self%t(:,mpi%rank), t1
      print '(a,9f12.6)', ' averages:', self%v(:,mpi%rank), self%average (t1)
    end if
  else
    self%v(self%nt,mpi%rank) = v
    self%w(self%nt,mpi%rank) = w
    self%t(self%nt,mpi%rank) = t
  end if
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Merge information with that from a node neighbor
!===============================================================================
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
  call trace%end()
END SUBROUTINE merge

!===============================================================================
!> Return a global average
!===============================================================================
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
    v = lagrange%sequence(t, self%t(:,i), self%v(:,i), min(self%order,self%istep))
    w = lagrange%sequence(t, self%t(:,i), self%w(:,i), min(self%order,self%istep))
    weight = weight + w
    average = average + v*w
  end do
  average = average
  call trace%end()
END FUNCTION average

END MODULE global_mod
