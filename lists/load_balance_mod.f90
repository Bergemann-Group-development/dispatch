!===============================================================================
!> Keep track of neighbor ranks and their loads, and redistribute.
!>
!> Each rank regularly updates a local queue measure -- the detailed definition 
!> of which may vary.  Each time a boundary task updates, this module is called,
!> and makes a decision whether to swap the ownership to one of the nbor ranks
!> that it normally sends a copy of the task to.
!>
!> The load measure is stored in the an mpi_window_t data type, which is also 
!> able to retrieve the measure from other ranks.  This requires only the rank
!> number, and results from an MPI window request to the rank in question.
!> There is thus no need to maintain one data type instance per rank.
!>
!> A running mean of task_list%nq -- the task list queue size -- is a useful
!> ingredient in the decision making. The running mean should be at least
!> several times the size of the task list, task_list_n, with a reasonable
!> default.
!>
!> An instantaneous queue measure is total_cost/(current_cost+average_cost),
!> where the total_load is the update cost of all tasks in the task list,
!> average cost is the total load divided by the number of tasks, and
!> current_cost is the load of the tasks in the ready queue.  The actual queue 
!> measure is a running mean of this instantaneous measure.
!>
!> That queue measure, which increases dramatically as the current_cost drops
!> towards zero, has a value that lies between essentially unity and task_list%n,
!> and there is no reason to load balance until the measure becomes rather large.
!>
!> Whether it is critical enough to warrant starting task swapping depends on
!> If the queue measure on average is small, it means that all ranks have plenty
!> of tasks to update.  Hence the need for load balancing should be based on
!> the actual value of the measure, except when it is high on all ranks the
!> critical value should be even larger, to avoid violent task swapping when all
!> ranks have short queues.
!>
!> One can turn the running mean measure back into an effective queue size, since
!> the measure may be written as N*average_cost/(current_cost+average_cost) =
!> N/(Q+1), so Q = measure/N - 1.  Then, a possible parametrization is to say
!> swapping happens when an nbor rank has a Q smaller than critical, AND the
!> local Q is larger than the nbor Q by at least a certain factor (such as 2).
!>
!> This has several benefits:
!> 1) when a rank starts getting a Q below the critical one, it will be given
!>    extra tasks first by the nbor rank that has the largest Q
!> 2) when all ranks have large Q values, it is unlikely that a single one will
!>    have an even much larger Q, so task swapping will be unlikely
!> 3) if all ranks have very low Q values, task swapping will only happen by
!>    taking tasks that consistently have at least somewhat larger values
!>
!> This strategy thus works well under both high and low rank load conditions.
!>
!> With regard to the sequence numbering of packets, the event sequence could be
!> as follows, for example:
!>
!>     rank A                                   rank B
!> event  tag bits                          event  tag bits
!> update
!> send    10 boundary+swap_request   -->   recv    10 virtual+swap_request
!>                                          update
!> recv    11 virtual                <--    send    11 boundary
!> 
!>                                          update
!> recv    12 virtual                <--    send    12 boundary
!===============================================================================
MODULE load_balance_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE list_mod
  USE link_mod
  USE bits_mod
  USE shared_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The load_balance_t data type holds the externally visible procedures
  !-----------------------------------------------------------------------------
  type, public:: load_balance_t
    logical                :: is_on, dry_run
    real                   :: queue
    type(lock_t)           :: lock
    class(list_t), pointer :: task_list
  contains
    procedure:: init
    procedure:: set
    procedure:: measure_queue
    procedure:: check_queue
    procedure:: print
  end type
  !-----------------------------------------------------------------------------
  ! Static module parameters and data types
  !-----------------------------------------------------------------------------
  integer, save :: cadence=1              ! how often load balancing is done
  integer, save :: length=5               ! length of running mean
  real, save    :: duration=0.0           ! limit the duration, if non-zero
  real, save    :: q_min=5.               ! limit below which ranks get tasks
  real, save    :: q_max=50.              ! limit above which ranks loose tasks
  integer, save :: verbose=0              ! verbosity
  type(load_balance_t), save, public  :: load_balance    ! public instance
CONTAINS

!===============================================================================
!> Initialize the load balance parameters
!===============================================================================
SUBROUTINE init (self, task_list)
  class(load_balance_t)   :: self
  class(list_t), pointer  :: task_list
  !.............................................................................
  class(link_t), pointer :: link
  logical, save :: is_on=.false.          ! task swapping actually done
  logical, save :: dry_run=.false.        ! load balancing is active
  integer       :: iostat, na
  real          :: value
  namelist /load_balance_params/ &
    is_on, dry_run, cadence, length, duration, q_min, q_max, verbose
  !-----------------------------------------------------------------------------
  if (mpi%size==1) &
    return
  call trace%begin ('load_balance_t%init')
  self%task_list => task_list
  if (verbose > 3) &
    write (io_unit%log, *) 'load_balance_t%init: associated(task_list) =', &
      associated(self%task_list)
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, load_balance_params, iostat=iostat)
  call io%namelist_warning ('load_balance_params', iostat)
  write (io%namelist, load_balance_params)
  !-----------------------------------------------------------------------------
  allocate (self%task_list%running_mean)
  na = self%task_list%na
  call self%task_list%running_mean%init (length*na, 1.0, verbose)
  call self%lock%init ('load')
  call self%task_list%mpi_window%init (real(na))
  if (verbose > 0) then
    value = self%task_list%running_mean%value(0.0)
    write (io_unit%log, *) 'load_balance_t%init: running_mean =', &
      associated(self%task_list%running_mean%ring), value
  end if
  self%is_on = is_on
  self%dry_run = dry_run
  call task_list%lock%set ('load_balance_t%init')
  link => task_list%head
  do while (associated(link))
    link%task%next_lb = link%task%nt
    link => link%next
  end do
  call task_list%lock%unset ('load_balance_t%init')
  call trace%end ()
END SUBROUTINE init
 
!===============================================================================
!> Turn load balancing on or off
!===============================================================================
SUBROUTINE set (self, on, active)
  class(load_balance_t) :: self
  logical, optional     :: on, active
  !.............................................................................
  if (present(on)) then
    !$omp atomic write
    self%is_on = on
  end if
  if (present(active)) then
    !$omp atomic write
    self%dry_run = active
  end if
END SUBROUTINE set

!===============================================================================
!> Decide whether to give up ownership of a task, to increase the load of an
!> nbor task that needs more load.
!===============================================================================
SUBROUTINE check_queue (self, link)
  class(load_balance_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  real:: my_queue, nb_queue
  logical:: printed, virtual
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%is_on) then
    call self%task_list%lock%set ('check_queue')
    call self%lock%set ('check_queue')
    call trace%begin ('load_balance%check_queue', itimer=itimer)
    task => link%task
    if (task%istep > task%next_lb) then
      task%next_lb = task%next_lb + cadence
      !-------------------------------------------------------------------------
      ! Start by computing the local queue measure, and go ahead if large enough
      !-------------------------------------------------------------------------
      my_queue = self%measure_queue ()
      if (self%queue > q_max .or. self%dry_run) then
        !-----------------------------------------------------------------------
        ! Check the queue measure on all rank nbors and swap if too low
        !-----------------------------------------------------------------------
        nbor => link%nbor
        virtual = .false.
        do while (associated(nbor))
          if (nbor%task%is_set (bits%virtual)) then
            virtual = .true.
            nbtask => nbor%task
            !--- use hash_table here, to reduce MPI window calls ---
            nb_queue = self%task_list%mpi_window%get (nbtask%rank)
            !nb_queue = nbtask%nq
            printed = .false.
            if (verbose > 1) call print
            if (nb_queue < q_min .and. task%rank_steps > task%nt) then
            !if (nb_queue < q_min) then
              if (verbose > 0) call print
              if (.not. self%dry_run) then
                write (io_unit%mpi,*) task%id, 'swap to', nbtask%rank
                call self%task_list%swap (link, .true., nbtask%rank)
              end if
              exit
            end if
          end if
          nbor => nbor%next
        end do
      end if
    end if
    call trace%end(itimer)
    call self%lock%unset ('check_queue')
    call self%task_list%lock%unset  ('check_queue')
  end if
contains
  subroutine print
    logical:: swap
    swap = (my_queue > q_max) .and. (nb_queue < q_min)
    if (swap) then
      1 format(1x,a,2(2x,2i6,f11.6,f7.1),2x,2f5.0,l4,i5)
      write (io_unit%log,1) &
       'load_balance_t%check_queue: task, my_q, nbor, nb_q, min, max =', &
        task%id, task%istep, task%time, my_queue, & 
        nbor%task%id, nbor%task%istep, nbor%task%time, nb_queue, &
        q_min, q_max, swap, task%rank_steps
    else
      write (io_unit%log,1) &
       'load_balance_t%check_queue: task, my_q, nbor, nb_q, min, max =', &
        task%id, task%istep, task%time, my_queue, & 
        nbor%task%id, nbor%task%istep, nbor%task%time, nb_queue, &
        q_min, q_max, swap, task%rank_steps
    end if
    printed = .true.
  end subroutine
END SUBROUTINE check_queue

!===============================================================================
!> Compute a measure of the task load, which for patches is the number of cells
!> divided by the time step.
!===============================================================================
FUNCTION measure_queue (self) RESULT (queue)
  class(load_balance_t)  :: self
  real                   :: queue
  !.............................................................................
  class(link_t), pointer :: link
  class(task_t), pointer :: task
  real                   :: task_cost, average_cost, measure
  real(8)                :: sum, rank_cost, queue_cost
  integer                :: n_task
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('load_balance_t%measure_queue', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Measure the cost of all tasks on the rank
  !-----------------------------------------------------------------------------
  if (verbose > 3) then
    write (io_unit%log, *) 'check_queue: associated(link) =', &
      associated(self%task_list%head)
    flush (io_unit%log)
  end if
  n_task = 0
  rank_cost = 0d0
  queue_cost = 0d0
!  call self%task_list%lock%set ('measure_queue')
  link => self%task_list%head
  do while (associated(link))
    task => link%task
    if (verbose > 3) then
      write (io_unit%log, *) link%task%id, task%is_clear (bits%virtual)
      flush (io_unit%log)
    end if
    if (task%is_clear (bits%virtual)) then
      select type (task)
      class is (experiment_t)
      task_cost = task%cost()
      end select
      n_task = n_task + 1
      rank_cost = rank_cost + task_cost
      if (task%is_set (bits%ready)) then
        queue_cost = queue_cost + task_cost
      end if
      if (verbose > 2) then
        write (io_unit%log,'(1x,a,1p,4g12.2)') &
         'load_balance_t%measure_queue: t_cost, r_cost, q_cost =', &
          task_cost, rank_cost, queue_cost
      end if
    end if
    link => link%next
  end do
!  call self%task_list%lock%unset ('measure_queue')
  !-----------------------------------------------------------------------------
  ! Define a queue measure, which goes to N when the queue size is zero, 
  ! and is N/(N+1) when the queue size is N, if task_cost is constant 
  !-----------------------------------------------------------------------------
  average_cost = rank_cost/n_task
  queue = queue_cost/average_cost
  measure = n_task/(queue+1.0)
  !-----------------------------------------------------------------------------
  ! A running average evens out fluctuations, but gives high weight to periods
  ! when the queue is short
  !-----------------------------------------------------------------------------
  measure = self%task_list%running_mean%value (measure)
  !-----------------------------------------------------------------------------
  ! Transform the measure back into an effective queue size, which equals the
  ! queue size when it is steady, but emphasizes periods of short queue
  !-----------------------------------------------------------------------------
  queue = n_task/measure - 1.0
  self%queue = queue
  !-----------------------------------------------------------------------------
  ! The current queue value is passed on via task buffers, and via mpi_window.
  ! The shared%n_queue value is used in patch_t%pack, replacing self%nq.
  !-----------------------------------------------------------------------------
  call self%task_list%mpi_window%put (mpi%rank, queue)
  !$omp atomic write
  shared%n_queue = queue
  if (verbose > 1) then
    write (io_unit%log,'(1x,a,l4,1p,4g12.2,2i5)') &
     'load_balance_t%measure_queue: a_cost, q_cost, measure, queue, nq, na =', &
      associated(self%task_list%running_mean%ring), average_cost, queue_cost, &
      measure, queue, self%task_list%nq, self%task_list%na
    flush (io_unit%log)
  end if
  call trace%end (itimer)
END FUNCTION measure_queue
 
!===============================================================================
!> Print a rank_info list
!===============================================================================
SUBROUTINE print (self)
  class(load_balance_t) :: self
  real                  :: queue
  integer               :: rank
  !.............................................................................
  do rank=0,mpi%size-1
    queue = self%task_list%mpi_window%get (rank)
    write (io_unit%log,'()') 'load_balance_t%print: rank, queue =',rank, queue
  end do
END SUBROUTINE print

END MODULE load_balance_mod
