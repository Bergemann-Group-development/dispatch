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
    real(8):: wallclock
    character(len=32), pointer:: name
    class(task_t), pointer:: parent
    class(*), pointer:: link            ! this could be useful
    type(lock_t):: lock                 ! OMP lock
    character(len=4):: kind='task'
  contains
    procedure:: init
    procedure:: update
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
  end type
  !-----------------------------------------------------------------------------
  ! Time slot info and procedures added in tasks_t
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: tasks_t
    real(8), dimension(:), pointer:: t
    real(8), dimension(:), pointer:: dt
    integer, dimension(:), pointer:: iit
    integer:: nt=3
    integer:: it=1
    integer:: new=2
    integer:: verbose=1
  end type
  integer, save:: id=0
  integer, save, public:: id_debug=0
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
    id = mpi%rank+1
    self%id = id
  else if (self%id==0) then
    id = id+mpi%size
    self%id = id
  end if
  scene%n_patch = self%id
  !$omp end critical
  write (io%log_unit,*) 'obtained patch id', self%id
  call self%lock%init (self%kind)
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class (task_t):: self
  !.............................................................................
  call trace_begin ('task_t%update')
  self%time = self%time + self%dtime
  if (io%master .and. io%verbose > 1) &
    print *, 'task_t%update: id, time =', self%id, self%time, self%name
  call trace_end
END SUBROUTINE

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
!> =============================================================================
LOGICAL FUNCTION is_ahead_of (self, task)
  class (task_t):: self, task
  !.............................................................................
  !$omp flush
  is_ahead_of = self%time + max(self%dtime,task%dtime)*task%grace >= task%time
  if (task%id==id_debug.or.io%verbose>4) &
    print'(i6,i4,2x,a,i6,3f9.5,l3)', self%id, omp_mythread(), 'mk is_ahead_of: ', &
    task%id, self%time, self%dtime*task%grace, task%time, is_ahead_of
END FUNCTION

SUBROUTINE print (self)
  class (task_t):: self
END SUBROUTINE

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%set')
  !$omp atomic
  self%status = ior(self%status, bit)
!  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE clear (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%clear')
  !$omp atomic
  self%status = iand(self%status, not(bit))
!  call trace_end
END SUBROUTINE
!> =============================================================================
LOGICAL FUNCTION is_set (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%is_set')
  is_set = (iand(self%status, bit) /= 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%is_clear')
  is_clear = (iand(self%status, bit) == 0)
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

END MODULE task_mod
