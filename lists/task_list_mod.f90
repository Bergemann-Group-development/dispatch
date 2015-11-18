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
  USE list_mod
  USE task_mod
  USE tasks_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE io_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
    real(8):: size(3) = 1.0_8
    integer:: n_fail=0, n_ready=0, n_nbor=0
  contains
    procedure:: init
    procedure:: update
    procedure:: check_mpi
    procedure:: check_nbors
    procedure:: check_ready
    procedure:: check_all
    procedure:: init_nbors
    procedure:: overlaps
    procedure:: distance
    procedure:: consistency
    procedure:: qconsistency
!    procedure:: print
  end type
CONTAINS

!> =============================================================================
!> Initialize the task list, by first initializing the list tasks, and then
!> checking if they are ready to execute (in general they all are)
!> =============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*):: name
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin('task_list_t%init '//name)
  call self%init_nbors
  link => self%head
  do while (associated(link))
    call link%task%init
    link => link%next
  end do
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call self%check_ready (link)
    !print *,'q', link%task%id, self%nq
    link => link%next
  end do
  link => self%queue
  do while (associated(link))
    print *,'q:', link%task%id
    link => link%next_time
  end do

  call trace_end
END SUBROUTINE init

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
  call trace_begin('task_list_t%update '//self%name)
  if (self%nq < 1) then
    n1 = self%nq
    call self%check_all
!    print *,'check_all', n1, self%nq
  end if
  call self%qconsistency(0)
  head => self%queue                                    ! queue start
  if (associated(head)) then
    !print*,head%task%id, self%nq
    self%queue => self%queue%next_time                  ! update the link
    self%nq = self%nq-1
    call self%check_mpi                                 ! check incoming MPI
    task => head%task                                   ! head task
    call task%update                                    ! update the task
    if (task%time > task%end_time) then                 ! finished:
      call self%remove_task (task)                      ! remove from self
    else                                                ! not finished:
      call task%clear(bits%ready)                       ! needs new check
      call self%check_nbors (head)                      ! any nbor ready?
    end if
  else
    print*,'no self%queue', self%nq, self%n
  end if
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Check for MPI messages
!> =============================================================================
SUBROUTINE check_mpi (self)
  class(task_list_t):: self
  !.............................................................................
  call trace_begin ('check_mpi')
  call trace_end
END SUBROUTINE check_mpi

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  !.............................................................................
  call trace_begin('task_list_t%check_nbors '//self%name)
  task => link%task                                   ! main task
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task                               ! neighbor task
    if (nbtask%rank==mpi%rank) then                   ! local rank:
      if (nbtask%is_clear(bits%ready)) &              ! needs checking?
        call self%check_ready (nbor%link)             ! nbor ready?
    else                                              ! foreign rank:
      call task%send_to (nbtask)                      ! send task data
    end if
    nbor => nbor%nbor                                 ! next nbor
  end do
  call self%check_ready (link)                        ! finally check link task
  self%n_nbor = self%n_nbor + 1
  call trace_end
END SUBROUTINE check_nbors

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_all (self)
  class(task_list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin('task_list_t%check_all '//self%name)
  link => self%head
  do while (associated (link))                        ! keep going until end
    call self%check_ready (link)                      ! link ready?
    link => link%next
  end do
  call trace_end
END SUBROUTINE check_all

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE check_ready (self, link)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, nbtask
  class(link_t), pointer:: nbor
  !.............................................................................
  call self%qconsistency(1)
  task => link%task                                   ! already checked?
  if (task%is_set(bits%ready)) return
!  call trace_begin('task_list_t%check_ready '//self%name)
  nbor => link%nbor                                   ! start on nbor list
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task
    if (task%ahead_of(nbtask)) then
      self%n_fail = self%n_fail+1
!      print*, task%id, ' check_ready failed on task', nbtask%id
!      call trace_end
      return                                          ! give up -- do nothing
    else
!      print*,'nbor task', nbtask%id, nbtask%time
    end if
    nbor => nbor%nbor                                 ! next nbor
  end do
  call task%set(bits%ready)                           ! mark as checked
  call self%qconsistency(2)
  call self%queue_by_time (link)                      ! add task to queue
  self%n_ready = self%n_ready+1
!  print *, task%id, ' check_ready OK', self%nq
  call self%qconsistency(3)
!  call trace_end
END SUBROUTINE check_ready

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE consistency (self, link, i)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  n = 0
  task => link%task
  nbor => link%nbor                                   ! start on nbor list
  if (associated(nbor)) id1 = nbor%task%id
  do while (associated (nbor))                        ! keep going until end
    id2 = nbor%task%id
    !print*,id2
    nbor => nbor%nbor                                 ! next nbor
    n = n+1
  end do
  if (n /= link%task%nbors) then
    print*, task%id, 'ERROR: n, nbors', n, id1, id2, i
  else
!    print*, task%id, 'consistent', id1, id2, i
  end if
END SUBROUTINE consistency

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE qconsistency (self, i)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  n = 0
  link => self%queue
  do while (associated(link))
    n = n+1
    link => link%next_time
  end do
  if (n /= self%nq) then
    print *, 'ERROR: qconsistency:', n, self%nq, i
    n = 0
    link => self%queue
    do while (associated(link))
      n = n+1
      print*,n,link%task%id
      link => link%next_time
    end do
    call mpi%end ('error')
  end if
END SUBROUTINE qconsistency

!> =============================================================================
!> The very first initialization of the neighbor lists.  After this, it can be
!> maintained without a global search.
!> =============================================================================
SUBROUTINE init_nbors (self)
  class(task_list_t):: self
  type(link_t), pointer:: link, scan
  class(link_t), pointer:: head, nbor
  integer:: n_add = 0
  !.............................................................................
  call trace_begin('task_list_t%init_nbors '//self%name)
  link => self%head                                     ! start checking links
  nullify (link%nbor)                                   ! make sure!
  do while (associated(link))
    nullify (link%nbor)                                 ! initially no nbor
    head => link%nbor
    n_add = 0
    scan => self%head                                   ! scan all links (FIXME)
    do while (associated(scan))                         ! until end
      if (.not.associated(scan,link)) then
        if (self%overlaps (scan%task, link%task)) then  ! if its task overlaps
          allocate (nbor)
          nbor%task => scan%task
          nbor%link => scan
          nbor%nbor => link%nbor
          link%nbor => nbor
          n_add = n_add+1
        end if
      end if
      scan => scan%next                                 ! continue search
    end do
    link%task%nbors = n_add
    if (io%master) print*, &
      'added ', n_add, ' neighbors to', link%task%id, associated(link%nbor)
    link => link%next                                   ! continue list
  end do
  call trace_end
END SUBROUTINE init_nbors

!> =============================================================================
!> =============================================================================
FUNCTION distance (self, task1, task2)
  class(task_list_t):: self
  type(task_t), pointer:: task1, task2
  real(8):: distance(3), size(3)
  !.............................................................................
  size = self%size
  distance = task1%position - task2%position + 1.5*size
  distance = abs(modulo(distance, size) - 0.5*size)
END FUNCTION distance

!> =============================================================================
!> =============================================================================
FUNCTION overlaps (self, task1, task2)
  class(task_list_t):: self
  class(task_t), pointer:: task1, task2
  logical:: overlaps
  !.............................................................................
!  call trace_begin ('task_lis_t%overlap')
  overlaps = (all(self%distance(task1,task2) < 0.501*(task1%size+task2%size)))
!  print'(2i4,3f8.2,3x,l3)', task1%id, task2%id, &
!    task1%position/task1%size-task2%position/task2%size, overlaps
!  call trace_end
END FUNCTION overlaps

END MODULE
