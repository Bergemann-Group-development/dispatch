!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Generic module for MHD task.  To be overloaded by an MHD layer and a timestep
!> layer.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE trace_mod
  USE io_mod
  implicit none
  integer, parameter:: ntime=7
  private
  type, public:: task_t
    real, dimension(3):: position, size, ds
    integer:: nt, nv
    integer:: n(3), ng(3)=2, gn(3), li(3), ui(3)
    integer:: it=1, new=2
    integer, pointer:: iit(:)
    real, pointer:: t(:)
    real, pointer:: dt(:)
    real:: time, dtime, end_time=1.0
    real:: out_time=0.1, out_next=0.0           ! controls output
    real:: courant=0.2                          ! courant constant
    integer:: istep=0                           ! time step counter
    integer:: iout=0                            ! output counter
    real, dimension(:,:,:,:,:,:), pointer:: mem
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  logical:: is_allocated = .false.
CONTAINS
!> =============================================================================
!> Initialize the task memory and task mesh.  For this we need to be told the
!> task dimensions (which include the number of variable and the number of time
!> slices), the task position, and the task size
!> =============================================================================
SUBROUTINE init (self, dims, position, size)
  class(task_t):: self
  integer:: dims(5)
  real:: position(3), size(3)
  !.............................................................................
  self%position = position
  if (is_allocated) return                      ! check
  is_allocated = .true.                         ! mark
  call trace_begin ('task_t%init')
  self%size = size                              ! patch size
  self%n  = dims(1:3)                           ! patch dimensions
  self%gn = self%n + 2*self%ng                  ! including guard cells
  self%li = 1       + self%ng                   ! lower interior index
  self%ui = self%gn - self%ng                   ! upper interior index
  self%nt = dims(4)                             ! number of time slices
  self%nv = dims(5)                             ! number of variables
  allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
  associate (nx=>self%gn(1), ny=>self%gn(2), nz=>self%gn(3))
    allocate (self%mem(nx,ny,nz,self%nv,self%nt,2))
  end associate
  call trace_end
END SUBROUTINE init
!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  !.............................................................................
  call trace_begin ('task_t%update')
  if (self%time >= self%out_time) then
    call self%output
    self%out_next = self%out_next + self%out_time
  end if
  self%istep = self%istep + 1
  self%dt(self%it) = self%dtime
  self%t(self%new) = self%time
  self%it = self%new
  self%new = mod(self%new,ntime)+1
  self%t(self%new) = self%time
  call trace_end
END SUBROUTINE update
!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self)
  class(task_t):: self
  !.............................................................................
  if (io%master) print*,'task_t%output: should do output at t =', self%time
END SUBROUTINE output
!> =============================================================================
END MODULE
