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
  namelist /task_list_params/ debug
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
  !-----------------------------------------------------------------------------
  ! Set virtual and boundary bits
  !-----------------------------------------------------------------------------
  !call self%reset_status
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
  do while (self%n > 0)
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
  !.............................................................................
  call trace_begin('task_list_t%update')
  call self%check_flags
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
      write (io_unit%log,'(a,i4,2x,a,i7,2x,a,i7,2x,a,i4,i4)') &
      'thread', omp_mythread, 'task', task%id, &
      'step', task%istep, 'n, nq', self%n, self%nq
    select type (task)
    class is (patch_t)
    call task%rotate
    !$omp atomic
    n_update = n_update + product(task%n)
    end select
    if (task%is_set(bits%boundary)) then                ! boundary patch?
      !-------------------------------------------------------------------------
      ! Pack the experiment patch data into a mesg (allocates mesg and mesg%reqs),
      ! elevating the task type to make sure to use the most comprehensive pack
      ! procedure (in case there is an experiment specific one)
      !-------------------------------------------------------------------------
      select type (task)
      class is (experiment_t)
        call task%pack (mesg)                           ! allocates mesg
      end select
      call head%send_to_vnbors (mesg)                   ! send to virtual nbors
    end if
    !-----------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the queue,
    ! and its neighbors should also be checked for return to the queue.
    !-----------------------------------------------------------------------------
    call task%clear(bits%ready)                         ! needs new check
    call self%check_nbors (head)                        ! any nbor ready?
    !-----------------------------------------------------------------------------
    ! Tasks that have finished are subtracted from the task list count but are
    ! not removed, since their data may be needed by other tasks.  They are,
    ! however, not added back onto the ready queue.
    !-----------------------------------------------------------------------------
    if (task%has_finished()) then  ! finished:
      call task%output (self%name)                      ! check output
      !$omp atomic
      self%n = self%n-1
      !$omp flush
    end if
    call task%info (self%nq, self%n)                    ! print info on stdout
  end if
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Find the link with a specified task id.  If there isn't any, create one!
!> =============================================================================
FUNCTION find_task (self, id) RESULT (link)
  class(task_list_t):: self
  integer:: id
  class(link_t), pointer:: link
  class(experiment_t), pointer:: task
  !.............................................................................
  link => self%head
  do while (associated(link))
!write (io_unit%log,*) mpi%rank,' find task: id, link%id, virtual =', id, link%task%id, link%task%is_set(bits%virtual)
    if (link%task%id == id) then
!write (io_unit%log,*) mpi%rank,' find task: match', id; flush(6)
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
write (io_unit%log,*) mpi%rank,' find task: new task', id, task%box; flush(6)
  call task%init
  !-----------------------------------------------------------------------------
  ! Create a neighbor list, and update the nbor lists of the nbors, before adding
  ! the link to the task list
  !-----------------------------------------------------------------------------
  call self%init_nbors (link)
  call link%update_nbors
  call self%append_link (link)
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
  class(experiment_t), pointer:: patch
  !-----------------------------------------------------------------------------
  link => self%find_task (mesg%id)
  task => link%task
  select type (task)
  class is (experiment_t)
    patch => task
    call patch%unpack (mesg%buffer)
    deallocate (mesg%buffer)
    deallocate (mesg)
  end select
  call self%check_nbors (link)
END SUBROUTINE unpack

!> =============================================================================
!> Check for and unpack MPI messages
!> =============================================================================
SUBROUTINE check_mpi (self)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg, this
  integer, save:: itimer=0
  !.............................................................................
  if (mpi%size <= 1) return
  call trace_begin ('task_list_mod::check_mpi', 1, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check for new or completed incoming and outgoing messages
  !-----------------------------------------------------------------------------
  call mpi_mesg%check
  !-----------------------------------------------------------------------------
  ! Get and unpack available message buffers
  !-----------------------------------------------------------------------------
  call mpi_mesg%get (mesg)
  do while (associated(mesg))
    this => mesg
    if (io%verbose>0) then
      write (io_unit%log,*) 'unpack msg id, sender', mesg%id, mesg%sender
      flush (io_unit%log)
    end if
    call self%unpack (mesg)
    call mpi_mesg%get (mesg)
  end do
  call trace_end (itimer)
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
