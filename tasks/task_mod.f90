!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE mpi_mod
  USE bits_mod
  USE trace_mod
  USE scene_mod
  USE mpi_mesg_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type, public:: task_t
    integer:: id = 0
    integer:: status = 0
    integer:: rank = 0
    integer:: n_check = 0
    integer:: n_nbors = 0
    integer:: istep = 0
    integer:: iout = 0
    integer:: level
    integer:: nt=3
    integer:: it=1
    integer:: new=2
    integer:: verbose=1
    integer:: nq=999
    integer:: n_behind=0
    real(8):: time = 0d0
    real(8):: end_time = 1d0
    real(8):: out_next = 0d0
    real(8):: out_time = 1d0
    real(8):: print_next = 0d0
    real(8):: print_time = 1d0
    real(8):: dtime = 0d0
    real(4):: grace = 0.0
    real(8):: position(3) = 0d0
    real(8):: velocity(3) = 0d0
    real(8):: size(3) = 0.9999d0
    real(8):: gsize(3)                  ! convenience
    real(8):: box(3)=1d0                ! the periodic wrapping size
    real(8):: wallclock
    integer, dimension(:), pointer:: iit
    real(8), dimension(:), pointer:: t
    real(8), dimension(:), pointer:: dt
    character(len=32), pointer:: name
    class(task_t), pointer:: parent
    type(lock_t):: lock                 ! OMP lock
    character(len=4):: kind='task'
    character(len=12):: type='task_t'
    real:: latency=0.0                  ! MPI message latency
  contains
    procedure:: init
    procedure:: dnload => void
    procedure:: update => void
    procedure:: rotate
    procedure:: output
    procedure:: info
    procedure:: is_ahead_of
    procedure:: set
    procedure:: clear
    procedure:: print
    procedure:: is_set
    procedure:: is_clear
    procedure:: allocate => void
    procedure:: deallocate => void
    procedure:: unique_id
    procedure:: load_balance
    procedure:: has_finished
    procedure:: overlaps
    procedure:: distance                ! return the distance btw two patches
    procedure:: pack
    procedure:: unpack
  end type
  integer, save:: id=0
CONTAINS

!> =============================================================================
SUBROUTINE void (self)
  class (task_t):: self
  print*, 'ERROR: task_t%void called'
END SUBROUTINE

!> =============================================================================
!> Define unique patch ids.  A critical region is OK; this is only done once.
!> =============================================================================
SUBROUTINE init (self)
  class (task_t):: self
  logical:: first_time=.true.
  !.............................................................................
  call trace_begin ('task_t%init')
  !$omp critical
  if (first_time) then
    first_time = .false.
    if (self%id==0) then
      id = mpi%rank+1
      self%id = id
    end if
  else if (self%id==0) then
    id = id+mpi%size
    self%id = id
  end if
  id = self%id
  scene%n_patch = self%id
  !$omp end critical
  if (io%verbose>2) &
    write (io%log_unit,*) 'patch id', self%id
  call self%lock%init (self%kind)
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE rotate (self)
  class(task_t):: self
  integer:: iv, nv, new, i
  real(8):: time
  integer, save:: itimer=0
  real(8), pointer:: rptr
  integer, pointer:: iptr
  !.............................................................................
  call trace_begin ('task_mod::rotate',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Call timestep update procedure to fill the next time slot with values
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,'(a,7i10)')   'task_mod::rotate iit =', self%iit
    write (io_unit%log,'(a,7f10.6)') 'task_mod::rotate   t =', self%t(self%iit)
  end if
  !-----------------------------------------------------------------------------
  ! Time slot rotation.  Note that self%time should be the LAST of these
  ! updates, since we do not want to enforce a critical region, and the
  ! new self%time is what a task is judged on!  Also, note that the only time
  ! slot information accessed from other tasks is self%it, self%time, self%iit,
  ! self%t(iit(1:nt-1)), and self%dt(iit(1:nt-1)).  The iit((nt) entries are
  ! not accessed.
  !-----------------------------------------------------------------------------
  !!omp critical (time_cr)
  time = self%time + self%dtime
  rptr => self%dt(self%it)
  !$omp atomic write
  rptr = self%dtime                                     ! just updated
  self%dt(self%new)= self%dtime                         ! next estimate
  self%t(self%new) = time                               ! initial time
  new = self%new
  !$omp atomic write
  self%it = new                                         ! update time slot, new
  self%new = mod(self%new,self%nt)+1                    ! increment / rotate
  do i=1,self%nt-1
    iptr => self%iit(i)
    !$omp atomic write
    iptr = self%iit(i+1)
  end do
  iptr => self%iit(self%nt)
  !$omp atomic write
  iptr = self%new                                       ! new right-most slot
  !$omp atomic update
  self%istep = self%istep + 1
  !$omp atomic write
  self%time = time
  if (io%verbose>2) &
   write (io_unit%log,*) 'task id:',self%id,'  it:',self%it,'  new:',self%new
  if (io%verbose>3) &
   write (io_unit%log,*) 'task iit:', self%iit
  if (io%verbose>3) &
   write (io_unit%log,*) 'task t(iit):', self%t(self%iit)
  !!omp end critical (time_cr)
  call trace_end (itimer)
END SUBROUTINE rotate

!> =============================================================================
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class (task_t):: self
  character(len=64), optional:: experiment_name
  !.............................................................................
  call trace_begin ('task_t%output')
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Print info to stdout
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(task_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
END SUBROUTINE info

!> =============================================================================
!> Check if self (which is a nbor task) is ahead of task (which is the one to
!> possibly move to the ready queue), using self%dtime*task%grace as the grace
!> period, since we want to limit the extrapolation in the nbor task to at most
!> task%grace*self%dtime.
!> =============================================================================
LOGICAL FUNCTION is_ahead_of (self, task)
  class (task_t):: self, task
  real(8):: nbtime, nbdtime
  !.............................................................................
  !$omp atomic read
  nbtime = self%time
  !$omp atomic read
  nbdtime = self%dtime
  is_ahead_of = nbtime + nbdtime*task%grace >= task%time
  if (task%id==io%id_debug.or.io_unit%verbose>4) &
    print'(i6,i4,2x,a,i6,3f9.5,l3)', self%id, omp_mythread, 'mk is_ahead_of: ', &
    task%id, self%time, self%dtime*task%grace, task%time, is_ahead_of
END FUNCTION

SUBROUTINE print (self)
  class (task_t):: self
END SUBROUTINE

!> =============================================================================
!> Set status bits
!> =============================================================================
SUBROUTINE set (self, bits)
  class (task_t):: self
  integer:: bits
  !.............................................................................
!  call trace_begin ('task_t%set')
  !$omp atomic
  self%status = ior(self%status, bits)
!  call trace_end
END SUBROUTINE

!> =============================================================================
!> Clear status bits
!> =============================================================================
SUBROUTINE clear (self, bits)
  class (task_t):: self
  integer:: bits
  !.............................................................................
!  call trace_begin ('task_t%clear')
  !$omp atomic
  self%status = iand(self%status, not(bits))
!  call trace_end
END SUBROUTINE

!> =============================================================================
!> Check if ANY of the bits is set (argument can be sum of bits)
!> =============================================================================
LOGICAL FUNCTION is_set (self, bits)
  class (task_t):: self
  integer:: bits, status
  !.............................................................................
!  call trace_begin ('task_t%is_set')
  !$omp atomic read
  status = self%status
  is_set = (iand(status, bits) /= 0)
!  call trace_end
END FUNCTION

!> =============================================================================
!> Check that ALL bits are clear (argument can be sum of bits)
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bits)
  class (task_t):: self
  integer:: bits, status
  !.............................................................................
!  call trace_begin ('task_t%is_clear')
  !$omp atomic read
  status = self%status
  is_clear = (iand(status, bits) == 0)
!  call trace_end
END FUNCTION

!> =============================================================================
LOGICAL FUNCTION has_finished (self)
  class (task_t):: self
  integer:: bit
  !.............................................................................
  has_finished = (self%time + 1e-6*self%dtime > self%end_time) &
    .or. self%is_set(bits%virtual)
  if (has_finished) call self%set (bits%frozen)
END FUNCTION

!> =============================================================================
FUNCTION unique_id (self)
  class (task_t):: self
  integer:: unique_id
  !$omp atomic
  id = id+1             ! FIXME
  unique_id = id
END FUNCTION unique_id

!> =============================================================================
!> The load balancing procedure is part of the negotiation between a patch and
!> its virtual neighbors.  We use the pro-active approach, where the decision
!> to migrate to a different rank is taken by the task itself.  It can then be
!> implemented with a minimum of complication, in that one of the last things
!> the task does before broadcasting itself to its neighbors is to change its
!> rank.  The new rank should check if an incoming virtual patch in fact has
!> opted to become a boundart patch instead.  It comes with a complete list of
!> neighbors already.
!>
!> To decide whether to migrate to another rank, the following parameters should
!> be considered, for both tasks:
!>
!> task%load    : the load of the task, in milliseconds per update
!> task%dt      : the time step per update
!> task%time    : the current time of the task
!> rank%time    : the current time of the neighbor rank
!> rank%surplus : the extra capacity available on the nbor rank
!>
!> If the time difference (rank%time-task%time)/task%dt is positive it means
!> that even the most delayed task on the other rank is ahead of the current
!> task, which speaks in favor of migrating.
!>
!> If the nbor rank has extra capacity, this also speaks in favor of migration
!>
!> A good measure of the load of a rank is the number of wall clock seconds per
!> simulation time unit.  This can be computed by running through the list of
!> tasks on the rank, adding for each task its cost in ms/time_unit, and finally
!> dividing by the number of (real) cores available to work on the tasks.
!> =============================================================================
SUBROUTINE load_balance (self)
class(task_t):: self

END SUBROUTINE load_balance

!> =============================================================================
!> Check if there is overlap, using the full size (including guard zones), since
!> the guard zones may be used in interpolations for qualified time slices. The
!> definition of overlap also MUST be symmetric, since it is being used to define
!> neighbor relations.
!> =============================================================================
FUNCTION overlaps (self, task)
  class(task_t):: self
  class(task_t):: task
  logical:: overlaps
  !.............................................................................
  overlaps = all(abs(distance(self,task)) <= 0.5000000001_8*(self%gsize+task%gsize))
END FUNCTION overlaps

!> =============================================================================
!> Signed distance between the centers of two tasks in the periodic box
!> =============================================================================
FUNCTION distance (self, task) RESULT (out)
  class(task_t):: self
  class(task_t):: task
  real(8):: out(3)
  !.............................................................................
  out = self%position-task%position
  out = modulo(out+0.5_8*self%box,max(self%box,1d-30))-0.5_8*self%box
END FUNCTION distance


SUBROUTINE pack (self, mesg)
  class(task_t):: self
  class(mesg_t), pointer:: mesg
  call mpi%abort ('task_t%pack called')
END SUBROUTINE

SUBROUTINE unpack (self, buffer)
  class(task_t):: self
  integer, dimension(:), pointer:: buffer
  call mpi%abort ('task_t%unpack called')
END SUBROUTINE

END MODULE task_mod
