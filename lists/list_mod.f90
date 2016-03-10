!> *****************************************************************************
!> Module with list handling for generic class tasks_t objects
!> *****************************************************************************
MODULE list_mod
  USE io_mod
  USE omp_mod
  USE mpi_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE patch_mod
  implicit none
  private
  type, public:: nbor_t
    class(task_t), pointer :: task => null()          ! task
    class(link_t), pointer :: nbor => null()          ! next nbor
    class(link_t), pointer :: link => null()          ! nbor's link (maybe not needed)
  end type
  type, extends(nbor_t), public:: link_t
    class(link_t),  pointer :: next => null()         ! next link in task list
    class(link_t),  pointer :: prev => null()         ! previous link in task list
    class(link_t),  pointer :: parent => null()       ! parent task link
    class(link_t),  pointer :: next_time => null()    ! next link in time order
  contains
    procedure:: download
    procedure:: add_nbor_by_quality
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    procedure:: print_nbors
  end type
  type, public:: list_t
    character(len=64):: name=''
    class(link_t), pointer:: head => null()
    class(link_t), pointer:: tail => null()
    class(link_t), pointer:: queue => null()
    integer:: n = 0
    integer:: nq = 0
    real(8):: size(3) = 1.0_8
    integer:: n_behind=0, n_check=0, n_ready=0, n_nbor=0
    integer:: debug=0
    integer:: n_tasks=1
    !---------------------------------------------------------------------------
    ! from patch_list_mod
    !---------------------------------------------------------------------------
    !integer:: dims(3)
  contains
    procedure:: init
    procedure:: update
    procedure:: reset
    procedure:: append
    procedure:: remove_link
    procedure:: remove_head
    procedure:: check_nbors
    procedure:: check_ready
    procedure:: check_all
    procedure:: check_queue
    procedure:: init_nbors
    procedure:: overlaps
    procedure:: distance
    procedure:: consistency
    procedure:: qconsistency
    procedure:: statistics
    procedure:: append_list
    procedure:: remove_task
    procedure:: add_by_quality
    procedure:: queue_by_time
    procedure:: print_queue
    procedure:: print => print_list
  end type
  integer, parameter:: mnbor=100
  public task2patch
CONTAINS

!> =============================================================================
!> Initialize list
!> =============================================================================
SUBROUTINE init (self, name)
  class(list_t):: self
  character(len=*), optional:: name
  !.............................................................................
  if (present(name)) then
    self%name = name
    print*,'init_list: ', name
  end if
END SUBROUTINE init

!> =============================================================================
!> Update simple list
!> =============================================================================
SUBROUTINE update (self)
  class(list_t):: self
END

!> =============================================================================
!> Reset list, by deallocating memory used by tasks, and deallocating tasks and
!> links
!> =============================================================================
SUBROUTINE reset (self)
  class(list_t):: self
  class(link_t), pointer:: link, old
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%reset '//self%name)
  link => self%head
  do while (associated(link))
    task => link%task
    if (associated(link)) then
      call task%deallocate
      deallocate (task)
    end if
    old => link
    link => link%next
    deallocate (old)
  end do
  nullify (self%head)
  nullify (self%tail)
  call trace_end
END SUBROUTINE reset

!> =============================================================================
!> =============================================================================
SUBROUTINE append (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin ('list_t%append '//self%name)
  !xomp critical (list_cr)
  allocate (link)
  link%task => task
  task%link => link
  if (associated(self%tail)) then
    link%prev => self%tail
    link%prev%next => link
    self%tail => link
  else
    self%head => link
    self%tail => link
  end if
  self%n = self%n+1
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE append

!> =============================================================================
!> Remove a link from the list, taking special care with the head and tail
!> =============================================================================
SUBROUTINE remove_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%remove_link '//self%name)
  !xomp critical (list_cr)
  next => link%next
  !-------------------------------------------------------------------------
  ! Check for tail task
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => link%prev                  ! point backwards, possibly null
  else
    self%tail => link%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Check for head task
  !-------------------------------------------------------------------------
  if (associated(link%prev)) then           ! not head task
    link%prev%next => next                  ! make prev point to next
  else                                      ! link was head
    self%head => next                       ! make head point to next
  end if
  task => link%task
  deallocate (link)                         ! free to remove link
  self%n = self%n-1
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Remove the head link
!> =============================================================================
SUBROUTINE remove_head (self)
  class(list_t):: self
  class(link_t), pointer:: head, link, next
  !.............................................................................
  call trace_begin ('list_t%remove_head '//self%name)
  !xomp critical (list_cr)
  head => self%head
  next => head%next
  !-------------------------------------------------------------------------
  ! Check for tail task (only one task left)
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => head%prev                  ! point backwards, possibly null
  else
    self%tail => head%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Remove the head link
  !-------------------------------------------------------------------------
  self%head => next                         ! make head point to next
  self%n = self%n-1
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE remove_head

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  integer:: n_check, n_ready, n_blocked, n_cleared, n_remain
  !.............................................................................
  call trace_begin('list_t%check_nbors '//self%name)
  !xomp critical (list_cr)
  task => link%task                                   ! main task
  n_check=0; n_ready=0; n_blocked=0; n_cleared=0; n_remain=0
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task                               ! neighbor task
    n_check = n_check+1
    if (nbtask%is_clear(bits%ready)) then             ! needs checking?
      n_blocked = n_blocked+1
      call self%check_ready (nbor%link)               ! pointer back
      if (nbtask%is_set(bits%ready)) then
        n_cleared = n_cleared+1
      else
        n_remain = n_remain+1
      end if
    else
      n_ready = n_ready+1
    end if
    nbor => nbor%next                                 ! next nbor
  end do
  call self%check_ready (link)                        ! finally check link task
  self%n_nbor = self%n_nbor + 1
  if (io%verbose>2) print '(a,i6,l3,12(2x,a,i4))', &
    'task:',task%id, task%is_set(bits%ready), 'level:',task%level, &
    'neighbors:',n_check, 'blocked:',n_blocked, 'ready:', n_ready, &
    'cleared:',n_cleared, 'remain:', n_remain
  if (io%verbose>3) call self%print
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE check_nbors

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_all (self)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin('list_t%check_all '//self%name)
  !xomp critical (list_cr)
  link => self%head
  do while (associated (link))                        ! keep going until end
    call self%check_ready (link, 1)                   ! link ready?
    link => link%next
  end do
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE check_all

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE check_ready (self, link, set_debug)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, nbtask
  class(link_t), pointer:: nbor
  integer:: set_debug, debug
  optional:: set_debug
  !.............................................................................
  call self%qconsistency(1)
  task => link%task                                     ! already checked?
  if (task%is_set(bits%ready)) return
  call trace_begin('list_t%check_ready '//self%name,2)
  nbor => link%nbor                                     ! start on nbor list
  do while (associated (nbor))                          ! keep going until end
    nbtask => nbor%task
    if (nbtask%is_ahead_of(task)) then
      nbor => nbor%next                                 ! next nbor
    else
      self%n_check = self%n_check+1
      if (self%debug==1) &
        print*, task%id, ' check_ready failed on task', nbtask%id
      call trace_end
      return                                            ! give up -- do nothing
    end if
  end do
  call task%set(bits%ready)                             ! mark as checked
  call self%qconsistency(2)
  call self%queue_by_time (link)                        ! add task to queue
  self%n_ready = self%n_ready+1
  if (self%debug>0) print *, task%id, ' found by check_ready, nq =', self%nq
  call self%qconsistency(3)
  call trace_end
END SUBROUTINE check_ready

!> =============================================================================
!> Check if the number of nbor tasks is consistent
!> =============================================================================
SUBROUTINE consistency (self, link, i)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  if (self%debug < 2) return
  call trace_begin ('list_t%consitency',1)
  n = 0
  task => link%task
  nbor => link%nbor                                   ! start on nbor list
  if (associated(nbor)) id1 = nbor%task%id
  do while (associated (nbor))                        ! keep going until end
    id2 = nbor%task%id
    !print*,id2
    nbor => nbor%next                                 ! next nbor
    n = n+1
  end do
  if (n /= link%task%n_nbors) then
    print*, task%id, 'ERROR: n, nbors', n, id1, id2, i
  else
!    print*, task%id, 'consistent', id1, id2, i
  end if
  call trace_end
END SUBROUTINE consistency

!> =============================================================================
!> Check if the number of queued tasks is consistent
!> =============================================================================
SUBROUTINE qconsistency (self, i)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  if (self%debug < 2) return
  call trace_begin ('list_t%qconsistency',0)
  n = 0
  link => self%queue
  do while (associated(link))
    n = n+1
    link => link%next_time
!    print*,n
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
  call trace_end
END SUBROUTINE qconsistency

!> =============================================================================
!> The very first initialization of the neighbor lists.  After this, it can be
!> maintained without a global search.
!> =============================================================================
SUBROUTINE init_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add = 0
  !.............................................................................
  call trace_begin('list_t%init_nbors ')
  link => self%head                                     ! start checking links
  do while (associated(link))
    nullify (link%nbor)                                 ! initially no nbor
    n_add = 0
    scan => self%head                                   ! scan all links (FIXME)
    do while (associated(scan))                         ! until end
      if (.not.associated(scan,link)) then
        if (self%overlaps (scan%task, link%task)) then  ! if its task overlaps
          allocate (nbor)
          nbor%task => scan%task
          nbor%link => scan
          nbor%next => link%nbor
          link%nbor => nbor
          n_add = n_add+1
          !---------------------------------------------------------------------
          ! Find the parent, as the patch on the next lower level that has the
          ! current path position inside it
          !---------------------------------------------------------------------
          if (nbor%task%level==link%task%level-1) then
            if (all(abs(nbor%task%position-link%task%position) < link%task%size*0.55d0)) then
              link%parent => nbor
              if (io%verbose>0) &
                print*, 'patch', link%task%id, ' parent:', nbor%task%id
            end if
          end if
        end if
      end if
      scan => scan%next                                 ! continue search
    end do
    link%task%n_nbors = n_add
    if (io%master .and. io%verbose>1) print*, &
      'added ', n_add, ' neighbors to', link%task%id, associated(link%nbor)
    link => link%next                                   ! continue list
  end do
  call trace_end
END SUBROUTINE init_nbors

!> =============================================================================
!> =============================================================================
FUNCTION distance (self, task1, task2)
  class(list_t):: self
  class(task_t), pointer:: task1, task2
  real(8):: distance(3), size(3)
  !.............................................................................
  size = self%size
  distance = task1%position - task2%position
  !distance = abs(distance)
  distance = abs(modulo(distance + 0.5*size, size) - 0.5*size)
END FUNCTION distance

!> =============================================================================
!> =============================================================================
FUNCTION overlaps (self, task1, task2, label)
  class(list_t):: self
  class(task_t), pointer:: task1, task2
  character(len=*), optional:: label
  logical:: overlaps
  real(8):: dist(3)
  !.............................................................................
  call trace_begin ('task_lis_t%overlap', 3)
  dist = self%distance(task1,task2)
  overlaps = all(self%distance(task1,task2) < 0.50000001_8*(task1%size+task2%size))
  if (self%debug>3.or.present(label)) then
    block
    real(8):: limit(3)
    limit = 0.501*(task1%size+task2%size)
    print'(a,2i4,2(3f8.2,3x),l3)', 'overlap:', task1%id, task2%id, &
     (task1%position-task2%position)/limit, dist/limit, overlaps
    end block
  end if
  call trace_end
END FUNCTION overlaps

!> =============================================================================
!> =============================================================================
SUBROUTINE statistics (self)
  class(list_t):: self
  print '(a,3(5x,a,i6))', 'output_experiment:', &
    'n_ready:', self%n_ready, &
    'n_check:', self%n_check, &
    'n_nbor:' , self%n_nbor
END SUBROUTINE statistics

!> =============================================================================
!> =============================================================================
SUBROUTINE append_list (self, task_list)
  class(list_t):: self
  class(list_t), pointer:: task_list
  class(link_t), pointer:: link
  integer:: n
  !-----------------------------------------------------------------------------
  call trace_begin ('list_t%append_list')
  !xomp critical (list_cr)
  self%tail%next => task_list%head
  self%tail      => task_list%tail
  self%n = self%n + task_list%n
  link => self%head
  n = 0
  do while (associated(link))
    link => link%next
    n = n+1
  end do
  if (n==self%n) then
    print *, trim(self%name)//': OK append of '//task_list%name
  else
    print *, trim(self%name)//': inconsistent append of '//task_list%name
  end if
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE append_list

!> =============================================================================
SUBROUTINE remove_task (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link, prev
  !.............................................................................
  call trace_begin('list_t%remove_task')
  !xomp critical (list_cr)
  link => self%head
  nullify (prev)
  do while (associated(link))
    if (link%task%id == task%id) then
      if (associated(prev)) then
        prev%next => link%next
      else
        self%head => link%next
      end if
      self%n = self%n - 1
      goto 9
    end if
    prev => link
    link => link%next
  end do
  print*,'ERROR: failed to remove from task list, task', task%id
9 continue
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE remove_task


!> =============================================================================
!> Initialize the nbors chain of a new patch (self), by
!> 1) adding the parent
!> 2) searching the parent nbors, adding patches that overlap
!> Each time it adds an nbor, it should also add itself to that nbor's nbor list
!> =============================================================================
SUBROUTINE make_new_nbors (self, selfp, parent)
  class(link_t):: self
  class(link_t), pointer:: selfp, parent, nbor
  class(patch_t), pointer:: p1, p2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbors', 1)
  !xomp critical (list_cr)
  p1 => task2patch (self%task)
  call self%make_new_nbor (selfp, parent)       ! make me a neighbor of my parent
  nbor => parent%nbor                           ! start on parent's nbors chain
  do while (associated(nbor))
    p2 => task2patch (nbor%task)
    if (.not.associated(self%task,nbor%task) &
      .and. p1%overlaps (p2)) &   ! if my task overlaps with that task
      call self%make_new_nbor (selfp, nbor%nbor)! make a neighbor (and make me one)
    nbor => nbor%next                           ! continue on parents nbors chain
  end do
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE make_new_nbors

!> =============================================================================
!> Add a link to a new neighbor to my nbors chain, and add myself into that
!> neighbor's neighbor chain.  An nbor chain link has three pointers: a pointer
!> to its task, a pointer to the next link in the chain, and a pointer back to
!> the task list link to the task.  We need to make and add two now nbor chain
!> links; one for our own nbors chain, and one for that's nbors chain.  Both
!> are added in quality order
!> =============================================================================
SUBROUTINE make_new_nbor (self, selfp, that)
  class(link_t):: self
  class(link_t), pointer:: selfp, that
  class(link_t), pointer:: link1, link2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbor', 1)
  allocate (link1)                              ! for our own nbors chain
  link1%task => that%task                       ! the link points to that task
  link1%nbor => that                            ! link to the nbor link
  link1%parent => that%parent                   ! direct link to task parent
  call self%add_nbor_by_quality (link1)         ! add to that chain in quality order
  allocate (link2)                              ! for that's nbors chain
  link2%task => self%task                       ! that's new link points to my task
  link2%nbor => selfp                           ! pointer back to self
  link2%parent => selfp%parent                  ! direct link to self parent
  call that%add_nbor_by_quality (link2)         ! add to my nbors chain in quality order
  call trace_end
END SUBROUTINE make_new_nbor

!> =============================================================================
!> Add a link node into the nbors list of self, in increasing quality order.
!> self must be a link whose 'nbors' points to a chain of nbors
!> =============================================================================
SUBROUTINE add_nbor_by_quality (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_quality')
  !xomp critical (list_cr)
  nullify (prev)
  next => self%nbor
  do while (associated(next))
    if (.not.associated(next%task,this%task)) then
      if (next%task%level >= this%task%level) then
        !-------------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-------------------------------------------------------------------------
        this%next => next
        if (associated(prev)) then
          prev%next => this
        else
          self%nbor => this
        end if
        self%task%n_nbors = self%task%n_nbors+1
        call trace_end
        goto 9
      end if
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This task has the highest quality, so append it
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => this
    nullify (this%next)
  else
    self%nbor => this
  end if
  self%task%n_nbors = self%task%n_nbors+1
9 continue
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE add_nbor_by_quality

!> =============================================================================
!> Update the nbor chains in the self list (possibly already existing) by,
!> for each link in that list, searching from the given start point for new
!> neighbors.
!> =============================================================================
SUBROUTINE search_nbors_by_overlap (self, start)
  class(link_t), target:: self
  class(link_t), pointer:: link1, link2, start, prev
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: p1, p2
  !> ===========================================================================
  !> For each task link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('list_t%search_nbors_by_overlap', 3)
  link1 => self                                 ! start on the self list
  if (io%verbose>3) print*,' id, size, position, distance:'
  do while (associated(link1))
    prev => link1%nbor                         ! existing nbor list head (empty?)
    link2 => start                              ! start searching throuch nbor chain
    do while (associated(link2))
      if (.not.associated(link1%task,link2%task)) then
        !if (io%verbose>3) &
        !  print'(2(i5,3(3f8.4,2x)),l3)', &
        !    link1%task%id, link1%task%size, link1%task%position, link1%task%distance(link2%task), &
        !    link2%task%id, link2%task%size, link2%task%position, link2%task%distance(link1%task), &
        !    link1%task%overlaps(link2%task)
        p1 => task2patch(link1%task)
        p2 => task2patch(link2%task)
        if (p1%overlaps(p2)) then
          allocate (nbor)                       ! new nbor (has no nbor)
          nbor%nbor => link2                    ! pointer back to link2
          link1%task%n_nbors = link1%task%n_nbors+1
          nbor%task => link2%task               ! but has a task
          if (associated(prev)) then            ! check if we are in list
            nbor%next => prev%next              ! could be empty
            prev%next => nbor                   ! insert the new nbor
          else                                  ! nbor list was empty
            link1%nbor => nbor                  ! start with the new nbor
          end if
          prev => nbor                          ! remember the previous nbor
        end if
      end if
      link2 => link2%next                       ! advance in the nbor chain
    end do
    link1 => link1%next                         ! advance in the self list
  end do
  call trace_end
END SUBROUTINE search_nbors_by_overlap

!> =============================================================================
!> Add a link node into the queue chain, in ascending time order.  When a link
!> to a task with time larger than this%task%time is found, or the end of the
!> list is fount, this link is inserted before.
!> =============================================================================
SUBROUTINE queue_by_time (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%queue_by_time '//self%name)
  !xomp critical (list_cr)
  call this%task%set(bits%ready)                       ! mark as checked
  if (io%verbose>2) print'(1x,a,i6,a,f10.5)', 'list_t%queue_by_time: adding task', &
    this%task%id, ' at time', this%task%time
  nullify (prev)
  next => self%queue
  do while (associated(next))
    if (associated(next%task, this%task)) then
      print*,'WARNMING: task', this%task%id, ' is already on ready queue'
    else if (next%task%time > this%task%time) then
      this%next_time => next
      if (associated(prev)) then
        prev%next_time => this
      else
        !xomp critical (queue_cr)
        self%queue => this
        !xomp end critical (queue_cr)
      end if
      self%nq = self%nq+1
      call trace_end
      goto 9
    end if
    prev => next
    next => next%next_time
  end do
  self%nq = self%nq+1
  if (associated(prev)) then
    prev%next_time => this
  else
    !xomp critical (queue_cr)
    self%queue => this
    !xomp end critical (queue_cr)
  end if
  nullify (this%next_time)
9 continue
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE queue_by_time

!> =============================================================================
!> Add a link node into the task list, in increasing quality order
!> =============================================================================
SUBROUTINE add_by_quality (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_by_quality')
  !xomp critical (list_cr)
  nullify (prev)
  next => self%head
  do while (associated(next))
    if (associated(next%task,this%task)) then
      print *,'ERROR: Trying to add task', this%task%id, ' which is already there'
      call trace_end
      goto 9
    end if
    if (next%task%level >= this%task%level) then
      !-------------------------------------------------------------------------
      ! Found a task we should be ahead of, so insert here
      !-------------------------------------------------------------------------
      this%next => next
      if (associated(prev)) then
        prev%next => this
      else
        self%head => this
      end if
      self%n = self%n+1
      call trace_end
      goto 9
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This task has the highest quality, so append it
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => this
  else
    self%head => this
  end if
  nullify (this%next)
  self%tail => this
  self%n = self%n+1
9 continue
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE add_by_quality

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE check_queue (self)
  class(list_t):: self
  class(link_t), pointer:: link
  real:: previous
  integer:: nq
  !.............................................................................
  link => self%queue
  previous = -1.0
  nq = 0
  do while (associated(link))
    if (link%task%time < previous) &
      print *,'ERROR: queue out of order at', link%task%id, link%task%time, &
      previous
    previous = link%task%time
    link => link%next_time
    nq = nq+1
  end do
  if (nq /= self%nq) print*,'ERROR: inconsistent number of tasks in queue', nq, self%nq
END SUBROUTINE check_queue

!> =============================================================================
!> Insert 'new' task before 'old' task
!> =============================================================================
SUBROUTINE insert (self, old, new)
  class(list_t):: self
  class(link_t), pointer:: old, new
  !.............................................................................
  call trace_begin ('list_t%insert '//self%name)
  !xomp critical (list_cr)
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    self%head => new
  end if
  new%next => old
  new%next%prev => new
  self%n = self%n + 1
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE insert

!> =============================================================================
FUNCTION task2patch (task) RESULT (patch)
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  select type (task)
  class is (patch_t)
  patch => task
  class default
  nullify(patch)
  stop
  end select
END FUNCTION

!> =============================================================================
!> Run through neighbor list and pick up guard cells, after which the whole
!> task is valid for download into other patches guard cells
!> =============================================================================
SUBROUTINE download (self, all_cells)
  implicit none
  class (link_t):: self
  logical, optional:: all_cells
  class (link_t),  pointer:: nbor
  class (patch_t), pointer:: patch, nbor_patch
  class (task_t), pointer:: task
  integer :: nguard, special
  real:: dmin
  !.............................................................................
  call trace_begin('link_t%download')
  patch => task2patch(self%task)
  nguard=0
  nbor => self%nbor
  do while (associated (nbor))
    task => nbor%task
    nbor_patch => task2patch(task)
    if (abs(patch%level-nbor_patch%level) < 2) then
      if (nbor_patch%level == patch%level+1) special = 6
      call nbor_patch%download_to (task, special, all_cells)
      nguard = nguard+1
    end if
    nbor => nbor%next
  end do
  if (io%verbose > 5) &
    print *, 'download: Number of neighbors with guard cells:', nguard
  call trace_end
END SUBROUTINE download

!> =============================================================================
!> Print a list of nbors
!> =============================================================================
SUBROUTINE print_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  if (present(label)) then
    if (trim(label) /= '') print'(a)', label
  end if
  print'(a,i5,2x,a,f14.8,2x,a,l3,2x,a,i3,2x,a,$)', &
    'task', self%task%id, 'time', self%task%time, &
    'ready', self%task%is_set(bits%ready), 'level', self%task%level, 'nbors:'
  if (io%verbose>4) print*,''
  link => self%nbor
  do while (associated(link))
    if (io%verbose>4) then
      print '(i7,f14.8,l3)', &
        link%task%id, link%task%time, link%task%is_set(bits%ready)
    else
      print '(i4,l2,$)', link%task%id, link%task%is_set(bits%ready)
    end if
    link => link%next
  end do
  print*,''
END SUBROUTINE print_nbors

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_queue (self)
  class (list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  link => self%queue
  print'(a,i4,":",$)', 'ready_queue', self%nq
  do while (associated(link))
    print '(i4,$)', link%task%id
    link => link%next_time
  end do
  print*,''
END SUBROUTINE print_queue

!> =============================================================================
!> =============================================================================
SUBROUTINE print_list (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t), pointer:: link, prev
  class(task_t), pointer:: task
  character(len=32):: type
  !.............................................................................
  call trace_begin ('list_t%print')
  print*,'list_t%print: ',trim(self%name), self%n
  link => self%head
  do while (associated(link))
    task => link%task
    type = ''
    !select type (task)
    !type is (task_t)
    !  type = 'task_t'
    !type is (particle_t)
    !  type = 'particle_t'
    !type is (planet_t)
    !  type = 'planet_t'
    !type is (star_t)
    !  type = 'star_t'
    !type is (pebble_t)
    !  type = 'pebble_t'
    !class is (patch_t)
      type = 'patch_t'
    !end select
    if (io%master) print '(a,i7,f10.6,3x,3f10.6,i4,i6,l5,2x,a)', 'list_t%print:', &
      task%id, task%time, task%position, task%status, task%n_nbors, associated(link%nbor), &
      type
    link => link%next
  end do
  call trace_end
END SUBROUTINE print_list

END MODULE list_mod
