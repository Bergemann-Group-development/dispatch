!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE mpi_mod
  USE bits_mod
  USE trace_mod
  USE scene_mod
  USE mpi_mesg_mod
  USE random_mod
  USE dll_mod
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
    integer:: level = 0
    integer:: nt=0
    integer:: it=1
    integer:: new=2
    integer:: verbose=1
    integer:: nq=999
    integer:: n_behind=0
    integer:: restart=-9
    integer:: state=0
    integer:: parentid=0
    integer:: thread=-1
    integer:: mem_thread=-1
    integer:: n_dump=0                  ! number of dump files
    integer:: ip                        ! rank local task id
    integer:: n_needed=1                ! number of tasks that needs it
    integer:: seq=0                     ! incremental sequence number
    real(8):: time = 0d0
    real(8):: atime = 0d0               ! for active queue
    real(8):: out_next = 0d0
    real(8):: print_next = 0d0
    real(8):: dtime = 0d0
    real(4):: grace = 0.0
    real(8):: position(3) = 0d0
    real(8):: origin(3)=0d0             ! the cartesian box origin (for BCs)
    real(8):: velocity(3) = 0d0         ! velocity of task in Cartesian coords.
    real(8):: size(3) = 0.9999d0
    real(8):: gsize(3)                  ! convenience
    real(8):: ds(3)                     ! convenience
    real(8):: box(3)=1d0                ! the periodic wrapping size
    real(8):: wallclock
    real(8):: unpack_time=0d0
    real(4):: quality = 0.0
    integer, dimension(:), pointer:: iit
    real(8), dimension(:), pointer:: t
    real(8), dimension(:), pointer:: dt
    character(len=64), pointer:: name
    class(task_t), pointer:: parent
    character(len=64):: kind='task'     ! solver name
    character(len=12):: type='task_t'   ! task type
    real:: latency=0.0                  ! MPI message latency
    real(8):: wc_last=0.0_8
    class(mesg_t), pointer:: mesg=>null()
    real(8):: sync_time=0d0
    real(8):: update_last=0.0_8
    real:: update_cadence=0.0
    logical:: syncing=.false.
    logical:: track=.false.
    logical:: periodic(3)=.true.
    logical:: rotated=.false.
    type(lock_t):: lock                 ! OMP lock
    type(random_t):: random             ! pseudo-random number generator
  contains
    procedure:: init
    procedure:: init_unique
    procedure:: init_id_reset
    procedure:: dealloc
    procedure:: dnload
    procedure:: update => void
    procedure:: rotate
    procedure:: info
    procedure:: timeslots
    procedure:: is_ahead_of
    procedure:: set
    procedure:: clear
    procedure:: print
    procedure:: is_set
    procedure:: is_clear
    procedure:: allocate => void
    procedure:: deallocate => void
    procedure:: allocate_mesg => void
    procedure:: unique_id
    procedure:: load_balance
    procedure:: has_finished
    procedure:: overlaps
    procedure:: distance                ! return the distance btw two patches
    procedure:: pack
    procedure:: unpack
    procedure:: is_relevant
    procedure:: test_update
    procedure:: solver_is
    procedure:: debug
  end type
  integer, dimension(:), pointer, save:: id => null()
  type(dll_t), pointer:: deleted_list
CONTAINS

!===============================================================================
SUBROUTINE void (self)
  class (task_t):: self
  print*, 'ERROR: task_t%void called'
END SUBROUTINE void

!===============================================================================
SUBROUTINE dnload (self, only)
  class (task_t):: self
  integer, optional:: only
  call mpi%abort('ERROR: task_t%dnload called')
END SUBROUTINE dnload

!===============================================================================
!> Define unique task id, unless already set before entry.  The task 
!> IDs must be unique across MPI ranks, and must initially satisfy 
!> modulo(self%id, mpi%size) == self%rank + 1
!===============================================================================
SUBROUTINE init (self)
  class (task_t):: self
  character(len=120):: id = &
    '$Id$ tasks/task_mod.f90'
  !-----------------------------------------------------------------------------
  call trace_begin ('task_t%init',2)
  call trace%print_id (id)
  call init_unique (self, keep=.true.)
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Define unique patch ids.  A critical region is OK; this is only done once.
!===============================================================================
SUBROUTINE init_unique (self, keep)
  class (task_t):: self
  logical, optional:: keep
  !.............................................................................
  integer:: rank
  logical, save:: first_time=.true.
  class(*), pointer:: idp
  !-----------------------------------------------------------------------------
  call trace_begin ('task_t%init_id')
  !$omp critical (task_init_cr)
  if (first_time) then
    first_time = .false.
    allocate (id(0:mpi%size-1))
    do rank=0,mpi%size-1
      id(rank) = 1+rank
    end do
    self%id = id(self%rank)
    !print *, self%rank, 'IDDBG first time, id =', id(self%rank)
    id(self%rank) = id(self%rank) + mpi%size
  !else if (deleted_list%n > 0) then
  !  idp => deleted_list%pop()
  !  select type (idp)
  !  type is (integer)
  !  id = idp
  !  end select
  !  dellocate (idp%car)
  !  dellocate (idp)
  !-----------------------------------------------------------------------------
  ! If ID already assigned and keep==.true., check consistency with MPI, and
  ! adjust the saved ID accordingly
  !-----------------------------------------------------------------------------
  else if (present(keep)) then
    if (self%id > 0) then
      if (mod(self%id-1,mpi%size) == self%rank) then
        id(self%rank) = max(id(self%rank),self%id)
      else if (io%verbose > 1) then
        print '("enforced ID =",i7," not consistent with MPI rank/size",2i6)', &
          self%id, self%rank, mpi%size
      end if
      !print *, self%rank, 'IDDBG present(keep) 1, id =', id(self%rank)
    else
      self%id = id(self%rank)
      !print *, self%rank, 'IDDBG present(keep) 2, id =', id(self%rank)
      id(self%rank) = id(self%rank) + mpi%size
    end if
  !-----------------------------------------------------------------------------
  ! Here we enforce a unique and MPI-consistent ID
  !-----------------------------------------------------------------------------
  else
    self%id = id(self%rank)
    !print *, self%rank, 'IDDBG generic, id =', id(self%rank)
    id(self%rank) = id(self%rank) + mpi%size
  end if
  if (self%id==1) self%track = .true.
  if (self%id==io%id_debug) self%track = .true.
  scene%n_patch = self%id
  !$omp end critical (task_init_cr)
  if (io%verbose>2) &
    write (io%log_unit,*) 'obtained patch id', self%id
  call self%lock%init (self%kind(1:4), id=self%id+2)
  call self%random%init
  call trace_end
END SUBROUTINE init_unique

!===============================================================================
!> Reset unique patch ids.  A critical region is OK; this is only done once.
!> This is done after, for example, `download%test`, which creates test patches.
!> It is set to one rather than zero because there must have already been a
!> patch with ID = 1 when `download%test` was called.
!===============================================================================
SUBROUTINE init_id_reset (self)
  class (task_t):: self
  logical:: first_time=.true.
  !.............................................................................
  call trace_begin ('task_t%init_reset')
  !$omp critical (task_init_cr)
  id = 1 + self%rank
  !$omp end critical (task_init_cr)
  call trace_end
END SUBROUTINE init_id_reset

!===============================================================================
!> Deallocate permanent arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class (task_t):: self
  !.............................................................................
  call trace_begin ('task_t%dealloc')
  if (associated(self%t  )) deallocate (self%t  )
  if (associated(self%dt )) deallocate (self%dt )
  if (associated(self%iit)) deallocate (self%iit)
  if (associated(self%mesg)) deallocate (self%mesg)
  call trace_end
END SUBROUTINE dealloc

!===============================================================================
!> Save the id of a deleted task, for resuse
!===============================================================================
SUBROUTINE save_id (self)
  class (task_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: node
  integer, pointer:: id
  !-----------------------------------------------------------------------------
  call trace_begin ('task_t%save_id')
  if (.not.associated(deleted_list%head)) then
    call deleted_list%init
  end if
  allocate (node, id)
  id = self%id
  node%car => id
  call deleted_list%append (node)
  call trace_end
END SUBROUTINE save_id

!===============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!===============================================================================
SUBROUTINE rotate (self)
  class(task_t):: self
  integer:: iv, nv, new, i
  real(8):: time
  integer, save:: itimer=0
  real(8), pointer:: rptr
  integer, pointer:: iptr
  !-----------------------------------------------------------------------------
  ! Implement setting of per-patch parameters here
  !-----------------------------------------------------------------------------
  if (io%processing > 0d0) then
    if (io%out_next > 0.0) self%out_next = io%out_next
    if (io%grace    > 0.0) self%grace    = io%grace
  end if
  if (self%rotated) return
  call trace_begin ('task_mod::rotate',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Call timestep update procedure to fill the next time slot with values
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,'(a,7i10)')   'task_mod::rotate iit =', self%iit
    flush (io_unit%log)
  end if
  if (omp_lock%tasks) then
    call self%lock%set ('rotate')
  end if
  !-----------------------------------------------------------------------------
  ! If this task has self%syncing set, then the time step has been set so as to
  ! arrive exactly ast sync_time.  To precent problems with round-off we set the
  ! time exactly.
  !-----------------------------------------------------------------------------
  if (self%syncing) then
    time = self%sync_time
    self%syncing = .false.
    if (io%verbose>2) &
      write (io_unit%mpi,*) 'task_t%rotate:: turning off syncing for now', self%id
  else
    time = self%time + self%dtime
  end if
  !-----------------------------------------------------------------------------
  ! Time slot rotation.  Note that self%time should be the LAST of these
  ! updates, since we do not want to enforce a critical region, and the
  ! new self%time is what a task is judged on!  Also, note that the only time
  ! slot information accessed from other tasks is self%it, self%time, self%iit,
  ! self%t(iit(1:nt-1)), and self%dt(iit(1:nt-1)).  The iit((nt) entries are
  ! not accessed.
  !-----------------------------------------------------------------------------
  rptr => self%dt(self%it)
  !$omp atomic write
  rptr = self%dtime                                     ! just updated
  self%dt(self%new)= self%dtime                         ! next estimate
  self%t(self%new) = time                               ! initial time
  new = self%new
  !$omp atomic write
  self%it = new                                         ! update time slot, new
  new = mod(new,self%nt)+1                              ! increment / rotate
  associate (snew=>self%new)
  !$omp atomic write
  snew = new
  end associate
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
  nullify (rptr, iptr)
  self%rotated = .true.
  if (omp_lock%tasks) then
    call self%lock%unset ('rotate')
  end if
  call trace_end (itimer)
END SUBROUTINE rotate

!===============================================================================
!> Print info to stdout
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(task_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
END SUBROUTINE info

!===============================================================================
!> Get timeslot info atomically.  Note that, while the thread updating the task
!> set the time and dtime of slot (nt) which, in terms of the indices retrieved
!> here is outside the range of concern, in rare cases it may be updated again,
!> in which case it will access the current slot (1) the next time, then (2)
!> and so on; therefore, it is prudent to use a task lock shared with rotate
!===============================================================================
SUBROUTINE timeslots (self, slots, times)
  class(task_t):: self
  real(8):: times(self%nt-1)
  integer:: slots(self%nt-1)
  integer:: i, j, it
  !-----------------------------------------------------------------------------
  !$omp atomic read
  it = self%it
  do j=1,self%nt-1
    i = 1 + mod(j+it,self%nt)
    if (j < self%nt-self%istep) then
      i = 1
    end if
    slots(j) = i
    !$omp atomic read
    times(j)  = self%t(i)
  end do
END SUBROUTINE timeslots

!===============================================================================
!> Check if self (which is a nbor task) is ahead of target (which is the one to
!> possibly move to the ready queue), using self%dtime*target%grace as the grace
!> period, since we want to limit the extrapolation in the nbor task to at most
!> target%grace*self%dtime.
!===============================================================================
LOGICAL FUNCTION is_ahead_of (self, target)
  class (task_t):: self, target
  real(8):: nbtime, nbdtime
  !.............................................................................
  !$omp atomic read
  nbtime = self%time
  !$omp atomic read
  nbdtime = self%dtime
  ! -- if a task is frozen (in time) assume it is forever ahead of other tasks
  if (self%is_set(bits%frozen)) then
    is_ahead_of = .true.
  ! -- a negative nbor time signals an invalid (e.g. initial RT) state
!  else if (nbtime < 0.0) then
!    is_ahead_of = .false.
  ! -- do not use a grace interval for different levels or for the first steps
  else if (self%level /= target%level .or. self%istep < 3) then
    is_ahead_of = nbtime >= target%time
  ! -- use a grace interval that is fraction of the nbor time step
  else
    is_ahead_of = nbtime + nbdtime*target%grace > target%time
  end if
  if (io_unit%verbose>1) then
    if (target%id==io%id_debug.or.io_unit%verbose>4) &
    print'(i6,i4,2x,a,i6,3f9.5,l3)', self%id, omp_mythread, 'mk is_ahead_of: ', &
    target%id, self%time, self%dtime*target%grace, target%time, is_ahead_of
  end if
END FUNCTION

SUBROUTINE print (self)
  class (task_t):: self
END SUBROUTINE

!===============================================================================
!> Set status bits
!===============================================================================
SUBROUTINE set (self, bits)
  class (task_t):: self
  integer:: bits
  !.............................................................................
!  call trace_begin ('task_t%set')
  !$omp atomic
  self%status = ior(self%status, bits)
!  call trace_end
END SUBROUTINE

!===============================================================================
!> Clear status bits
!===============================================================================
SUBROUTINE clear (self, bits)
  class (task_t):: self
  integer:: bits
  !.............................................................................
!  call trace_begin ('task_t%clear')
  !$omp atomic
  self%status = iand(self%status, not(bits))
!  call trace_end
END SUBROUTINE

!===============================================================================
!> Check if ANY of the bits is set (argument can be sum of bits)
!===============================================================================
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

!===============================================================================
!> Check that ALL bits are clear (argument can be sum of bits)
!===============================================================================
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

!===============================================================================
LOGICAL FUNCTION has_finished (self)
  class (task_t):: self
  integer:: bit
  !.............................................................................
  has_finished = (self%time > io%end_time)
END FUNCTION

!===============================================================================
FUNCTION unique_id (self)
  class (task_t):: self
  integer:: unique_id
  !$omp atomic
  id(self%rank) = id(self%rank) + mpi%size
  unique_id = id(self%rank)
END FUNCTION unique_id

!===============================================================================
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
!===============================================================================
SUBROUTINE load_balance (self)
class(task_t):: self

END SUBROUTINE load_balance

!===============================================================================
!> Check if there is overlap, using the size plus an extra cell margin. The
!> definition of overlap MUST be symmetric, since it is being used to define
!> neighbor relations.
!>
!> NOTE 1: This function was in the past NOT used by e.g. list_t%init_nbors, which
!> was instead using list_t%overlaps, based on the same criterion as below.
!>
!> NOTE 2: Since overlap and nbors can have a generalized definition it makes
!> sense that is valid for the whole task_t class.  If needed, the definition
!> can be overloaded with alternative definitions for particular tasks.
!===============================================================================
FUNCTION overlaps (self, task)
  class(task_t):: self
  class(task_t):: task
  logical:: overlaps
  !.............................................................................
  call trace%begin ('task_t%overlaps',2)
  overlaps = all(abs(distance(self,task)) <= 0.5_8*(self%size+task%size+self%ds+task%ds))
  call trace%end()
END FUNCTION overlaps

!===============================================================================
!> Signed distance between the centers of two tasks in a possibly periodic box
!>
!> NOTE: This function is normally NOT used (but remains so it can be overloaded
!> or explicitly called), since it is overloaded by the patch_t%distance, which
!> is patch geometry aware.
!===============================================================================
FUNCTION distance (self, task) RESULT (out)
  class(task_t):: self
  class(task_t):: task
  real(8):: out(3)
  !.............................................................................
  call trace%begin ('task_t%distance',2)
  out = self%position-task%position
  ! account for periodic wrapping
  where (self%periodic)
    out = modulo(out+0.5_8*self%box,max(self%box,1d-30))-0.5_8*self%box
  end where
  if (io%verbose>1 .and. (task%id==io%id_debug .or. self%id==io%id_debug)) then
    print *,'distance: self', self%position, self%id
    print *,'distance: task', task%position, task%id
    print *,'distance:  box', self%box
  end if
  call trace%end()
END FUNCTION distance

SUBROUTINE pack (self, mesg)
  class(task_t):: self
  class(mesg_t), pointer:: mesg
  call mpi%abort ('task_t%pack called')
END SUBROUTINE

SUBROUTINE unpack (self, mesg)
  class(task_t):: self
  class(mesg_t), pointer:: mesg
  call mpi%abort ('task_t%unpack called')
END SUBROUTINE

!===============================================================================
!> Stub, to be overloaded in experiment_mod to drop patches
!===============================================================================
FUNCTION is_relevant (self)
  class(task_t):: self
  logical:: is_relevant
  !-----------------------------------------------------------------------------
  is_relevant = .true.
END FUNCTION is_relevant

!===============================================================================
!> Test update procedure for dispatcher tests
!===============================================================================
SUBROUTINE test_update (self)
  class(task_t):: self
  self%dtime = self%random%ran1()
  self%dt(self%it) = self%dtime
END SUBROUTINE test_update

!===============================================================================
!> Test if solver leading characters match a given string
!===============================================================================
FUNCTION solver_is (self, kind)
  class(task_t):: self
  logical solver_is
  character(len=*):: kind
  solver_is = self%kind(1:len(kind)) == kind
END FUNCTION solver_is

!===============================================================================
!> Turn on debug output at and beyond a given level, and for id_debug task
!===============================================================================
FUNCTION debug (self, level)
  class(task_t):: self
  integer:: level
  logical debug
  !-----------------------------------------------------------------------------
  debug = (io%verbose >= level) .or. (self%id == io%id_debug)
END FUNCTION debug

END MODULE task_mod
