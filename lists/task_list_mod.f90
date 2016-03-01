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
  USE omp_timer_mod
  USE mpi_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE bits_mod
  USE trace_mod
  USE refine_mod
  USE download_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
  contains
    procedure:: init
    procedure:: execute
    procedure:: update
    procedure:: check_mpi
    procedure:: split
    procedure:: output
    procedure:: print => print_tasklist
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
  call io%init
  if (present(name)) self%name = name
  call trace_begin('task_list_t%init ')
  rewind (io%input); read(io%input, task_list_params, iostat=iostat)
  if (iostat==0 .and. io%master) then
    self%debug = debug
    if (io%master) write (*,task_list_params)
  end if
  !-----------------------------------------------------------------------------
  ! Compute neighbor list used to perform checks for ready status
  !-----------------------------------------------------------------------------
  call trace ('call self%init_nbors')
  call self%init_nbors
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call self%check_ready (link)
    link => link%next
  end do
  call self%print
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Execute the task list, updating it until it is empty.  Wth !$omp parallel here,
!> everything local to self%update is thread private.
!> =============================================================================
SUBROUTINE execute (self)
  class(task_list_t):: self
  class(link_t), pointer:: link, rem
  !$omp parallel
  do while (self%n > 0)
    call self%update
  end do
  !$omp end parallel
  link => self%head
  do while (associated(link))
    rem => link
    link => link%next
    call self%remove_link (rem)
  end do
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
  integer:: n1
  !.............................................................................
  call trace_begin('task_list_t%update')
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  call self%check_mpi                                 ! check incoming MPI
  !-----------------------------------------------------------------------------
  ! Return and try again, if the queue is empty; this can happen if the one task
  ! that frees up the rest is still running, with no other task ready
  !-----------------------------------------------------------------------------
  !$omp critical (list_cr)
  !call self%lock
  call trace('list_cr 0 begin')
  if (self%nq < 1) then
    n1 = self%nq
    call self%check_all
  end if
  call trace('list_cr 0 end')
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in a critical region
  !-----------------------------------------------------------------------------
  call trace('list_cr 1 begin')
  head => self%queue                                    ! queue start
  if (associated(head)) then
    self%queue => self%queue%next_time                    ! update the link
    self%nq = self%nq-1
  end if
  call trace('list_cr 1 end')
  !call self%unlock
  !$omp end critical (list_cr)
  !-----------------------------------------------------------------------------
  ! Having secured a head, we need to check if it is at all associcated
  !-----------------------------------------------------------------------------
  if (associated(head)) then
    task => head%task                                   ! head task
    !$omp critical (list_cr)
    call trace('list_cr a begin')
    if (io%verbose>0) &
      print'(a,i4,2x,a,i7,2x,a,i7,2x,a,i4,i4)', &
      'thread', omp_mythread(), 'task', task%id, &
      'step', task%istep, 'n, nq', self%n, self%nq
    !-----------------------------------------------------------------------------
    ! Prepare it for update, by downloading data that it depends on.  During the
    ! download we must prevent threads in the nbor list from rotating mem slots
    !-----------------------------------------------------------------------------
    select type (task)
    class is (patch_t)
    if (io%verbose>2) print*,'dmin',omp_mythread(),minval(task%mem(:,:,:,1,task%it,1)),1
    call download_link (head, all_cells=.false.)        ! guard cells
    if (io%verbose>2) print*,'dmin',omp_mythread(),minval(task%mem(:,:,:,1,task%it,1)),1
    end select
    !-----------------------------------------------------------------------------
    ! Check if refinement is needed on the task; if so this will push new tasks
    ! onto the queue, with the same task time; i.e., to the head of the queue
    !-----------------------------------------------------------------------------
   !if (refine%needed (task)) then                      ! if refinement needed
   !  call self%split (head)                            ! refined tasks go back to queue
   !end if
    !-----------------------------------------------------------------------------
    ! To ensure the initial state is output, check for doing output, then update,
    ! then check again.  The check is fast, but one could move the first call
    !-----------------------------------------------------------------------------
    call task%output(self%name)                         ! check output
    call trace('list_cr a end')
    !$omp end critical (list_cr)
    call task%update                                    ! update the task
    !----------------------------------------------------------------------------
    ! The update procedure in the generic patch_t module is responsible for
    ! updating the time, and rotating the memory slots where information about
    ! previous time steps are saved.
    !----------------------------------------------------------------------------
    !$omp critical (list_cr)
    call trace('list_cr 2 begin')
    select type (task)
    class is (patch_t)
    call task%rotate
    end select
    !-----------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the queue,
    ! and its neighbors should be checked for return to the queue.
    !-----------------------------------------------------------------------------
    call task%clear(bits%ready)                         ! needs new check
    call self%check_nbors (head)                        ! any nbor ready?
    call trace('list_cr 2 end')
    call task%output(self%name)                         ! check output
    !-----------------------------------------------------------------------------
    ! Tasks that have finished are subtracted from the task list count
    !-----------------------------------------------------------------------------
    if (task%time+0.5*task%dtime > task%end_time) then   ! finished:
      call trace('list_cr 3 begin')
      self%n = self%n-1
      call trace('list_cr 3 end')
    end if
    !$omp end critical (list_cr)
  end if
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Check for MPI messages
!> =============================================================================
SUBROUTINE check_mpi (self)
  class(task_list_t):: self
  !.............................................................................
  call trace_begin ('check_mpi',2)
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
  type(patch_t), pointer:: patch, this_patch
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
        !print*,'mk1: id,it,it,new,new',patch%id,patch%it,this_patch%it,patch%new,this_patch%new
        !print*,'mk2; id,it,iit',patch%id,this_patch%it, this_patch%iit
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
        if (this_patch%time > 0.0) call download_link (link, all_cells=.true.)
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
SUBROUTINE print_tasks (self)
  class (task_list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  link => self%head
  print'(a,i4,":",$)', 'task_list', self%n
  do while (associated(link))
    print '(i4,$)', link%task%id
    link => link%next
  end do
  print*,''
END SUBROUTINE print_tasks

!> =============================================================================
!> Print a table with task id and all neighbors
!> =============================================================================
SUBROUTINE print_tasklist (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  print*, label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
  call self%print_queue
END SUBROUTINE print_tasklist

END MODULE
