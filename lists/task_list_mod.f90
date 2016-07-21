!> ==================-===========================================================
!> $Id$
!> A task_list_t has the distinction that it must be able to call an initial
!> condition routine, which only exists in experiment_t tasks.  Procedures in
!> the module use these foreign procedures:
!>
!>    task_t%init         : to initialize a task
!>    task_t%update       : to update a task
!>    task_t%send_to      : send data to foreign task
!>    list_t%remove_task  : remove task from list
!>    list_t%remove_head  : remove head from list
!>
!> =============================================================================
MODULE task_list_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE trace_mod
  USE refine_mod
  USE download_mod
  USE timestep_mod
  USE buffered_output_mod
  USE timer_mod
  USE rank_info_mod
  USE random_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
  contains
    procedure:: init
    procedure:: execute
    procedure:: update
    procedure:: find_task
    procedure:: unpack
    procedure:: check_mpi
    procedure:: split
    procedure:: output
    procedure:: check_flags
    procedure:: print => print_tasklist
    procedure:: print_tasks
  end type
  type, public:: load_balance_t
    logical:: on, only_initial
    real:: cadence, info_cadence, duration
    real:: next_time, next_info
    real:: threshold, grace
    integer:: excess=0
  contains
    procedure:: check => check_load
    procedure:: measure => measure_load
  end type
  type(load_balance_t):: load_balance
  public link_t
CONTAINS

!> =============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!> =============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*), optional:: name
  class(link_t), pointer:: link
  integer:: iostat
  integer:: debug=0
  real, save:: cadence=10., info_cadence=1., threshold=10., grace=0.3
  real, save:: duration=0.0
  logical:: on=.false., only_initial=.false.
  namelist /task_list_params/ debug
  namelist /load_balance_params/ on, cadence, info_cadence, threshold, grace, &
    duration, only_initial
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace_begin('task_list_t%init ')
  call self%list_t%init (name)
  call mpi_mesg%init
  rewind (io%input); read(io%input, task_list_params, iostat=iostat)
  if (iostat==0 .and. io%master) then
    self%debug = debug
    if (io%master) write (*,task_list_params)
  end if
  rewind (io%input); read(io%input, load_balance_params, iostat=iostat)
  if (io%master) then
    if (io%master) write (*,load_balance_params)
  end if
  load_balance%on = on
  load_balance%only_initial = only_initial
  load_balance%cadence = cadence
  load_balance%info_cadence = info_cadence
  load_balance%next_time = cadence
  load_balance%next_info = info_cadence
  load_balance%threshold = threshold
  load_balance%grace     = grace
  load_balance%duration  = duration
  io%do_trace = .true.
  call rank_info_list%init
  io%do_trace = .false.
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!> =============================================================================
SUBROUTINE execute (self)
  USE omp_timer_mod
  class(task_list_t):: self
  class(link_t), pointer:: link, rem
  !-----------------------------------------------------------------------------
  ! Construct neighbor lists used to perform checks for ready status
  !-----------------------------------------------------------------------------
  call self%init_all_nbors
  call self%count_status
  !-----------------------------------------------------------------------------
  ! Set virtual and boundary bits
  !-----------------------------------------------------------------------------
  call self%reset_status
  !-----------------------------------------------------------------------------
  ! Check the ready status of all tasks
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call self%check_ready (link)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Start timer, and run until the task count is zero
  !-----------------------------------------------------------------------------
  call mpi%barrier ('task_list_mod::execute')
  call tic
  call timer%print
  !$omp parallel
  do while (self%na > 0)
    call self%update
  end do
  !$omp end parallel
  call timer%print
  call toc ('rank time', n_update)
  !if (io%master) &
  !  write (io_unit%log,*) 'download_differ % =', 100.*download%n_differ/ &
  !    max(1,(download%n_differ+download%n_same))
  call buffered_output%close
END SUBROUTINE execute

!> =============================================================================
!> Update the state of the task list, taking the steps necessary to update the
!> head task on the ready_queue, and check for consequences
!> =============================================================================
SUBROUTINE update (self)
  class(task_list_t):: self
  !.............................................................................
  class(link_t), pointer:: head
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg
  integer, save:: stalled=0, max_stalled=10000
  !.............................................................................
  call trace_begin('task_list_t%update')
  !call self%check_flags
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  call self%check_mpi                                 ! check incoming MPI
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in a very brief critical region
  !-----------------------------------------------------------------------------
  !$omp critical (list_cr)
  head => self%queue                                    ! queue start
  if (associated(head)) then
    !---------------------------------------------------------------------------
    ! Set the head of the queue to the next task in time.  The ready bit needs
    ! to remain on for now, to prevent other threads from checking this task
    ! before it is updated.  After the update, the clearing of the ready bit
    ! allows the state of the current task to be evaluated by other threads, as
    ! well as by the current thread.  Locking should be used to allow only a
    ! single thread and task to check and change the ready bit status, at any
    ! one time.
    !---------------------------------------------------------------------------
    self%queue => head%next_time
    head%task%nq = self%nq                              ! for info print
    !$omp atomic
    self%nq = self%nq-1                                 ! decrement queue count
  end if
  !$omp end critical (list_cr)
  if (associated(head)) then
    task => head%task                                   ! head task
    !---------------------------------------------------------------------------
    ! If the task is frozen (e.g. because it has finished), return
    !---------------------------------------------------------------------------
    if (task%is_set (bits%frozen+bits%virtual)) then
      if (task%iout == 0) call task%output(self%name)
      call trace_end
      return
    end if
    !----------------------------------------------------------------------------
    ! Update the task, whatever that means
    !----------------------------------------------------------------------------
    call download%download_link (head, all_cells=.false.)
    call task%output (self%name)                        ! check output
    call task%update                                    ! update the task
    !----------------------------------------------------------------------------
    ! The rotate procedure in the generic patch_t module is responsible for
    ! updating the time, and rotating the memory slots where information about
    ! previous time steps are saved.
    !----------------------------------------------------------------------------
    if (io%verbose>1) &
      write (io_unit%log,'(a,i4,2x,a,i7,2x,a,i7,2x,a,i4,i4,l3)') &
      'thread', omp_mythread, 'task', task%id, &
      'step', task%istep, 'n, nq', self%n, self%nq, associated(head%nbor)
    select type (task)
    class is (patch_t)
    call task%rotate
    !$omp atomic
    n_update = n_update + product(task%n)
    end select
    call task%info (self%nq, self%na)                   ! print info on stdout
    !-----------------------------------------------------------------------------
    ! Tasks that have finished are subtracted from the task list count but are
    ! not removed, since their data may be needed by other tasks (including on
    ! other ranks).  They are, however, not added back onto the ready queue.
    !-----------------------------------------------------------------------------
    if (task%has_finished()) then  ! finished:
      load_balance%on = .false.                         ! turn of load balancing
      call task%output (self%name)                      ! check output
      !$omp atomic
      self%na = self%na-1
    end if
    !-----------------------------------------------------------------------------
    ! If the task is a boundary patch, first check if it should be given to
    ! another rank.  A patch that has been given to another rank becomes
    ! a virtual patch, and the load balance procedure has already sent those
    ! of its nbors that became new boundary patches over to relevant ranks,
    ! while the patch itself is sent here, with bits%virtual+swap_reqest set.
    !-----------------------------------------------------------------------------
    if (task%is_set(bits%boundary)) then              ! boundary patch?
      head%task%nq = self%nq                          ! make sure to pass on
      call load_balance%check (self, head)            ! sell?
      call head%send_to_vnbors                        ! send to virtual nbors
    end if
    !-----------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the
    ! queue, and its neighbors should also be checked for return to the queue.
    ! If the task was turned into a virtual task by the load balance, its virtual
    ! bit will prevent it being checked, but its nbors still need to be checked
    !-----------------------------------------------------------------------------
    call self%check_nbors (head)                      ! any nbors ready?
    !$omp atomic write
    stalled = 0
    call mpi_mesg%diagnostics (0)
  else
    !$omp critical (stalled_cr)
    stalled = stalled+1
    if (stalled == 100) then
       call self%check_all
    else if (stalled == 101) then
       call self%check_oldest
    else
      call mpi_mesg%delay (stalled)
    end if
    if (stalled == max_stalled/2) call mpi_mesg%diagnostics (1)
    if (stalled > max_stalled) call mpi%abort ('exceeded max_stalled')
    !$omp end critical (stalled_cr)
  end if
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Compute a measure of the task load, which for patches is the number of cells
!> divided by the time step.
!> =============================================================================
SUBROUTINE measure_load (self, task_list)
  class (load_balance_t):: self
  class (task_list_t):: task_list
  class (link_t), pointer:: link
  class (task_t), pointer:: task
  real:: load, sum, sum_cost, cells, sum_cells, ready
  real(8):: dt, sum_dt, wc
  integer:: n, n1
  integer, save:: itimer
  logical:: ok, active
  !-----------------------------------------------------------------------------
  call trace_begin ('task_list_mod::measure_load', itimer=itimer)
  !
  rank_info%cadence = load_balance%info_cadence
  ok = .true.
  link => task_list%head
  do while (associated(link))
    ok = merge(.false., ok, link%task%rank==mpi%rank .and. link%task%dtime==0d0)
    link => link%next
  end do
  rank_info%ok = ok                             ! active patches have non-zero dtime
  if (.not. ok) then
    rank_info%dtime = 1.0
    call trace_end (itimer)
    return
  end if
  !
  sum_dt = 0.0
  sum_cost = 0.0
  sum_cells = 0.0
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      if (task%rank==mpi%rank) then
        cells = product(task%mesh%n)
        sum_cells = sum_cells + cells
        sum_cost = sum_cost + cells/task%dtime
      end if
    end select
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! cost = number of cell updates per code time unit
  !-----------------------------------------------------------------------------
  rank_info%cost = sum_cost                     ! total rank update cost
  rank_info%rank = mpi%rank                     ! local rank
  rank_info%dtime = sum_cells/sum_cost          ! cell-weighted dtime
  !$omp critical (list_cr)
  if (associated(task_list%queue)) then         ! queue info, if available
    wc = wallclock()
    if (wc > rank_info%wall) then
      rank_info%owall = rank_info%wall          ! previous wall time
      rank_info%wall  = wc                      ! wall time
      rank_info%otime = rank_info%time          ! previsou queue time
      rank_info%time  = 0.9*rank_info%time + &  ! smoothed ..
        0.1*task_list%queue%task%time           ! .. queue time
      rank_info%nq    = task_list%nq            ! queue size
    end if
  end if
  !$omp end critical (list_cr)
  if (io%verbose > 1) &
    write (io_unit%log,'(a,3i5,1p,3e12.3)') 'measure_load: n,n1,nq,load,cost,dt =', &
      n, n1, task_list%nq, rank_info%cost, rank_info%dtime
  call trace_end (itimer)
END SUBROUTINE measure_load

!> =============================================================================
!> Decide whether to give up ownership of a patch, to increase the load of an
!> nbor task that needs more load. After a rank passes the critical threshold
!> below, there is a delay before it's boundary patches are sent to its nbor
!> ranks, are unpacked, and end up in a load_balance comparison.  This will 
!> lead to a burst of patch transfers, and oscillatory behavior. It might thus
!> be better to transfer the LB info separately, in small packages that do not
!> take as long to transfer.  This can be done in the same send_to_vnbors. 
!>
!> If one instead gives up ownership of a patch as soon the threshold is hit,
!> and one counts this immediately as an improved load balance, then one can
!> prevent the oscillatory behavior.  The 'pending' improvement would be
!> cancelled when the patch update first comes back.  This speaks in favor of
!> a measure of load balance that both includes a measure of actual load (e.g.
!> wall clock seconds per code time unit), and a measure of the ready queue.
!> A too small ready queue should result in immediate action, and in immediate
!> improvement of the measure of load balance.
!>
!> The sum of 1./dtime for all tasks in the task list is a measure of load.
!> The part of this in the ready queue is a measure of the urgency of LB. The
!> boundary patch that was just updated, and which triggered the call to LB,
!> is/was counted towards that ready-load
!> =============================================================================
SUBROUTINE check_load (self, task_list, link)
  class (load_balance_t):: self
  class (task_list_t):: task_list
  class (link_t), pointer:: link, nbor
  class (task_t), pointer:: task
  class (patch_t), pointer:: patch
  class(rank_info_t), pointer:: nbor_info
  integer,save:: delay=0
  logical:: load_condition, time_condition
  real(8):: wc
  real:: load_diff, cost, time_diff
  logical, save:: make_estimate=.true.
  !-----------------------------------------------------------------------------
  wc = wallclock()
  patch => task2patch(link%task)
  !-----------------------------------------------------------------------------
  ! Update the load balance calculation for the local rank, and enter into the
  ! rank_info_list, which is also updated via MPI
  !-----------------------------------------------------------------------------
  if (wc > load_balance%next_info) then
    load_balance%next_info = load_balance%next_info + load_balance%info_cadence
    call load_balance%measure (task_list)
    call rank_info_list%send
    rank_info%patch_cost = product(patch%mesh%gn)/patch%dtime
    call rank_info_list%print (link%task%time)
  end if
  !-----------------------------------------------------------------------------
  ! If load balancing is on, and we have information, and it's time, or if an
  ! earlier case has left an excess, go ahead
  !-----------------------------------------------------------------------------
  if (.not.load_balance%on .or. &
      .not.rank_info%ok  .or. &
      wc > load_balance%duration) then
    return
  end if
  !-----------------------------------------------------------------------------
  ! Among the virtual patches that surround the link boundary patch, find one
  ! that needs more load and give the link patch to it
  !-----------------------------------------------------------------------------
  nbor => link%nbor
  do while (associated(nbor))
    if (io%verbose>1) &
      write (io_unit%log,'(a,2i9,3i6,2f12.6,l3)') 'LB: ids,nqs,rank,times', &
        link%task%id, nbor%task%id, task_list%nq, nbor%task%nq, nbor%task%rank, &
        nbor%task%time, link%task%time, nbor%task%is_set(bits%virtual)
    if (nbor%task%is_set (bits%virtual)) then
      nbor_info => rank_info_list%find (nbor%task%rank)
      if (nbor%task%rank == mpi%rank) then
        write (io_unit%log,*) &
          'task_list_mod::check_load virtual bit on my rank', nbor%task%id
        flush (io_unit%log)
      else if (associated(nbor_info)) then
        time_diff = rank_info%time_diff (nbor_info)/load_balance%grace
        time_diff = max(time_diff, -1.)
        if (rank_info%ok.and.nbor_info%ok.and. &
            random%ran1() < 0.05*(1.-exp(-time_diff))) then
          !$omp critical (loadb_cr)
          if (io%verbose>0) then
            print 1, &
              wc,'LB: rank',mpi%rank,' gives patch',link%task%id,' to',nbor%task%rank, &
              load_balance%excess, nbor_info%n_swap, &
              nbor_info%cost, rank_info%cost
              1 format(f12.6,2x,a,i6,a,i9,a,i6,2i5,1p,2g12.3)
            write (io_unit%log,2) &
              wc,'load_balance: giving patch',link%task%id, ' to',nbor%task%rank, &
              load_balance%excess, nbor_info%n_swap, &
              nbor_info%cost, rank_info%cost
              2 format(f12.6,2x,a,i9,a,i6,2i5,1p,2g12.3)
            flush (io_unit%log)
          end if
          call link%give_to (link, nbor%task%rank)
          if (io%verbose>0) &
            write (io_unit%log,*) 'swapped boundary to virtual:', link%task%id
          !$omp end critical (loadb_cr)
          call task_list%count_status
          exit
        end if
      else
        write (io_unit%log,*) &
          'check_load  ERROR: nbor_info not associated, rank', nbor%task%rank
        nbor_info => rank_info_list%find (nbor%task%rank, debug=.true.)
        flush (io_unit%log)
      end if
    else if (nbor%task%is_set (bits%external)) then
      write (io_unit%log,*) nbor%task%id, 'LB: external status on nbor from rank', nbor%task%rank
    else if (nbor%task%rank /= mpi%rank) then
      write (io_unit%log,*) nbor%task%id, 'LB: inconsistent status on nbor from rank', nbor%task%rank
    end if
    nbor => nbor%next
  end do
END SUBROUTINE check_load

!> =============================================================================
!> Find the link with a specified task id.  If there isn't any, create one!
!> If the patch already existed, it's nbor list should be OK, but the status
!> bit of it, and all its nbors, need to updated.
!> =============================================================================
FUNCTION find_task (self, id) RESULT (link)
  class(task_list_t):: self
  integer:: id
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  class(experiment_t), pointer:: task
  !.............................................................................
  call trace_begin ('task_list_mod::find_task')
  link => self%head
  do while (associated(link))
!write (io_unit%log,*) mpi%rank,' find task: id, link%id, virtual =', id, link%task%id, link%task%is_set(bits%virtual)
    patch => task2patch(link%task)
    if (patch%id == id) then
      if (io%verbose>1) &
        write (io_unit%log,*) mpi%rank,' find task: match', id, associated(patch%mem)
      flush(io_unit%log)
      call trace_end
      return
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! No task found, which means it's a new virtual task that needs to be created
  !-----------------------------------------------------------------------------
  allocate (link)
  allocate (task)
  link%task => task
  task%id = id
  task%box = self%size
  call task%patch_t%init
  if (io%verbose>1) &
    write (io_unit%log,*) mpi%rank,' find task: new task', id, associated(task%mem)
  flush(io_unit%log)
  !-----------------------------------------------------------------------------
  ! Append to the task_list.  The calling procedure must create an nbor list,
  ! and must also update the nbor list of the nbors
  !-----------------------------------------------------------------------------
  call self%append_link (link)
  call trace_end
END FUNCTION find_task

!> =============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!> =============================================================================
SUBROUTINE unpack (self, mesg)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  if (mesg%nbuf < 40) then
    call rank_info_list%unpack (mesg%buffer)
    return
  end if
  link => self%find_task (mesg%id)
  task => link%task
  !-----------------------------------------------------------------------------
  ! Guard against lingering extra messages to a swapped patch
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%boundary)) then
    write (io_unit%log,'(f12.6,2x,a,i9)') wallclock(), &
      'task_list_mod::unpack ERROR, received mpi_mesg for boundary task:', task%id
    deallocate (mesg%buffer)
    deallocate (mesg)
    return
  end if
  if (task%is_set (bits%ready)) then
    write (io_unit%log,'(f12.6,2x,a,i9)') wallclock(), &
      'task_list_mod::unpack ERROR, received mpi_mesg for task in ready queue:', task%id
    deallocate (mesg%buffer)
    deallocate (mesg)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Elevate the type of pointer, to make sure to use the experiment_t%unpack,
  ! if it exists -- if not the call will go to the generic patch_t%unpack.
  !-----------------------------------------------------------------------------
  select type (task)
  class is (experiment_t)
    call task%unpack (mesg%buffer)
    !---------------------------------------------------------------------------
    ! If the boundary+swap bits are set, this is a task that has just changed
    ! rank, and it needs to have its nbor relations re-initialized. This includes
    ! resorting (removing + re-adding) the nbor's nbor lists in rank order.
    !---------------------------------------------------------------------------
    if (task%is_set(bits%swap_request) .and. task%is_set(bits%boundary)) then
      call self%init_nbors (link)
      call task%clear (bits%swap_request+bits%ready)
      call link%update_nbor_status (link)
      call self%count_status
      if (io%verbose>1) & 
        print '(f12.6,2x,a,i6,a,i9,a,i6)', &
        wallclock(),'LB: rank',mpi%rank,' given patch',task%id,' by',mesg%sender
      if (io%verbose>0) &
        write (io_unit%log,*) 'swapped virtual to boundary:', task%id
    !---------------------------------------------------------------------------
    ! If the link has no nbors it is a newly created virtual task. Does it need
    ! an nbor list?  At least we can use the nbor list to check that the link is
    ! in its nbors nbor lists.  A new virtual task (where no task existed) means
    ! that some nbor of it has changed from internal to boundary, which will be
    ! checked by the test_nbor_status call below, but only if an nbor list exists.
    !---------------------------------------------------------------------------
    !else if (.not.associated(link%nbor)) then
    else if (task%is_set(bits%swap_request) .and. task%is_set(bits%virtual)) then
      if (io%verbose>0) &
        write (io_unit%log,*) 'new virtal patch:', task%id
      call self%init_nbors (link)
      call task%clear (bits%swap_request+bits%ready)
      call link%update_nbor_status (link)
      call self%count_status
    end if
    deallocate (mesg%buffer)
    deallocate (mesg)
  end select
  !-----------------------------------------------------------------------------
  ! As the task has now been updated, we need to check if any of the neighbors
  ! have become ready to update. If the task has just been swapped into being
  ! a boundary patch, it could possibly also be ready to update.
  !-----------------------------------------------------------------------------
  if (io%verbose>1) &
    call link%test_nbor_status
  call self%check_nbors (link)
END SUBROUTINE unpack

!> =============================================================================
!> Check for and unpack MPI messages
!> =============================================================================
SUBROUTINE check_mpi (self)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg
  logical:: wait
  !.............................................................................
  if (mpi%size <= 1) return
  call trace_begin ('task_list_mod::check_mpi', 1)
  !-----------------------------------------------------------------------------
  ! Check for new or completed incoming and outgoing messages
  !-----------------------------------------------------------------------------
  wait = self%nq < 4
  wait = .false.
  call mpi_mesg%check (wait)
  !-----------------------------------------------------------------------------
  ! Get and unpack available message buffers
  !-----------------------------------------------------------------------------
  call mpi_mesg%get (mesg)
  do while (associated(mesg))
    call self%unpack (mesg)
    call mpi_mesg%get (mesg)
  end do
  call trace_end
END SUBROUTINE check_mpi

!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split (self, this)
  class(task_list_t):: self
  class(link_t), pointer:: this, link
  integer:: nsplit(3)
  class(patch_t), pointer:: patch, this_patch
  class(task_t), pointer:: task
  class(link_t), pointer:: tail
  integer dir, i, ix, iy, iz
  real:: pos(3)
  !.............................................................................
  if (this%task%is_set(bits%not_leaf)) return
  call trace_begin('task_list_t%split', 1)
  nsplit = msplit
  task => this%task
  select type (task)
  class is (patch_t)
  this_patch => task
  end select
  where (this_patch%size==0 .or. this_patch%n==1)
    nsplit = 1
  end where
  do ix=1,nsplit(1)
    pos(1) = ix-(nsplit(1)+1)*0.5
    do iy=1,nsplit(2)
      pos(2) = iy-(nsplit(2)+1)*0.5
      do iz=1,nsplit(3)
        pos(3) = iz-(nsplit(3)+1)*0.5
        !------------------------------------------------------------------------
        ! Allocate a new task, 1/nsplit the size of original one, and adjust its
        ! position in each direction with "pos", in units of the new size.
        !------------------------------------------------------------------------
        allocate (patch, source=this_patch)
        patch%mem_allocated = .false.
        nullify (patch%mem)
        patch%size = this_patch%size/nsplit             ! new size
        patch%position = this_patch%position + pos*patch%size
        call patch%init                                 ! patch, for a new id
        patch%t   = this_patch%t                        ! newly allocated
        patch%dt  = this_patch%dt                       ! newly allocated
        patch%iit = this_patch%iit                      ! newly allocated
        !write (io_unit%log,*) 'mk1: id,it,it,new,new',patch%id,patch%it,this_patch%it,patch%new,this_patch%new
        !write (io_unit%log,*) 'mk2; id,it,iit',patch%id,this_patch%it, this_patch%iit
        patch%iit(1:patch%nt-1) = this_patch%it         ! fake history
        patch%out_next = patch%time                     ! force output now
        !------------------------------------------------------------------------
        ! Append links to the new patch in the patch list, sorted by quality
        !------------------------------------------------------------------------
        allocate (link)
        link%task => patch
        link%parent => this                     ! link is the link to the new patch
        call self%add_by_quality (link)         ! which is added to the patch list
        call link%make_new_nbors (link, this)   ! make a new nbors list for link
        !------------------------------------------------------------------------
        ! Download values from parent and neighbors, except at time==0.0, to
        ! allow recursive refinement there
        !------------------------------------------------------------------------
        if (this_patch%time > 0.0) call download%download_link (link, all_cells=.true.)
        !------------------------------------------------------------------------
        ! Insert the new patch at the head of the ready queue
        !------------------------------------------------------------------------
        call self%queue_by_time (link)          ! add patch to queue
        if (io%verbose>0) then
          print 1,'task_list_t%split:  patch, level, iout, time, size =', &
            patch%id, patch%level, patch%iout, patch%time, patch%size
          1 format(a,3i6,f10.6,2(2x,3f8.4))
          print 1,'task_list_t%split: parent, level, iout, time, size =', &
            this_patch%id, this_patch%level, this_patch%iout, this_patch%time, this_patch%size
        end if
      end do
    end do
  end do
  call self%check_queue
  call this_patch%set(bits%not_leaf)
  call trace_end
END SUBROUTINE split

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE output (self)
  class (task_list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call trace_begin('task_list%output')
  link => self%head
  do while (associated(link))
    call link%task%output(self%name)
    link => link%next
  end do
  call trace_end
END SUBROUTINE output

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_tasks (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !if (io%verbose < 2) return
  if (present(label)) then
    write (io_unit%log,'(a,i4,":")') 'task_list: '//trim(label), self%n
  else
    write (io_unit%log,'(a,i4,":")') 'task_list', self%n
  end if
  link => self%head
  do while (associated(link))
    write (io_unit%log,'(i4,g15.6,l5)') link%task%id, link%task%time, link%task%is_set(bits%ready)
    link => link%next
  end do
  write (io_unit%log,*) ''
END SUBROUTINE print_tasks

!> =============================================================================
!> Print a table with task id and all neighbors
!> =============================================================================
SUBROUTINE print_tasklist (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%master) write (io_unit%log,*) label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
  call self%print_queue
END SUBROUTINE print_tasklist

!> =============================================================================
!> Check debug flag files
!> =============================================================================
SUBROUTINE check_flags (self)
  class (task_list_t):: self
  logical:: changed
  !.............................................................................
  if (io%verbose > 0 .and. io%do_flags) then             ! avoid speed impact
    call trace_begin ('task_list_mod::check_flags')
    call icheck ('id_debug.flag', io%id_debug, changed)
    call icheck ('verbose.flag' , io%verbose , changed)
    call lcheck ('do_trace.flag', io%do_trace, changed)
    call trace_end
  end if
contains
  !=============================================================================
  logical function exists(file)
  character(len=*):: file
  inquire (file=file, exist=exists)
  end function
  !=============================================================================
  subroutine icheck (file, value, changed)
  character(len=*):: file
  integer:: value
  logical:: changed
  if (exists(trim(file))) then
    open (io_unit%flag, file=trim(file), form='formatted', status='old')
    read (io_unit%flag,*) value
    close (io_unit%flag,status='delete')
    write (io_unit%log,*) 'flag check:', file, value
    changed = .true.
  else
    changed = .false.
  end if
  end subroutine
  !=============================================================================
  subroutine lcheck (file, value, changed)
  character(len=*):: file
  logical:: value
  logical:: changed
  if (exists(trim(file))) then
    open (io_unit%flag, file=trim(file), form='formatted', status='old')
    read (io_unit%flag,*) value
    close (io_unit%flag,status='delete')
    write (io_unit%log,*) 'flag check:', file, value
    changed = .true.
  else
    changed = .false.
  end if
  end subroutine
END SUBROUTINE check_flags

END MODULE
