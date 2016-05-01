!> *****************************************************************************
!> Module with list handling for generic class task_t objects
!> *****************************************************************************
MODULE list_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE mpi_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE download_mod
  implicit none
  private
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
    type(lock_t):: lock
    character(len=4):: kind='list'
    logical:: face_nbors=.true.
  contains
    procedure:: init
    procedure:: update
    procedure:: reset
    procedure:: append
    procedure:: append_link
    procedure:: remove_link
    procedure:: remove_head
    procedure:: check_nbors
    procedure:: check_ready
    procedure:: check_all
    procedure:: check_queue
    procedure:: init_nbors
    procedure:: init_all_nbors
    procedure:: remove_parents
    procedure:: overlaps
    procedure:: distance
    procedure:: intpos
    procedure:: consistency
    procedure:: qconsistency
    procedure:: statistics
    procedure:: append_list
    procedure:: remove_task
    procedure:: reset_status
    procedure:: add_by_quality
    procedure:: queue_by_time
    procedure:: print_queue
    procedure:: print_queue_until
    procedure:: print_queue_times
    procedure:: print => print_list
  end type
  integer, parameter:: mnbor=100
  public link_t, task2patch
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
  call self%lock%init (self%kind)
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
  allocate (link)
  link%task => task
  call self%append_link (link)
  call trace_end
END SUBROUTINE append

!> =============================================================================
!> Append to the task list, counting active tasks only
!> =============================================================================
SUBROUTINE append_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(*), pointer:: anon
  logical, save:: first_time=.true.
  !.............................................................................
  call trace_begin ('list_t%append_link '//self%name)
  !$omp critical (list_cr)
  if (associated(self%tail)) then
    link%prev => self%tail
    link%prev%next => link
    self%tail => link
  else
    self%head => link
    self%tail => link
  end if
  if (link%task%is_clear(bits%virtual)) &
    self%n = self%n+1
  if (associated(link%task)) then
    link%task%link => link
    if (first_time) then
      first_time = .false.
      anon => link%task%link
      select type (anon)
      type is (link_t)
        print*,'MK type is link_t'
      class is (link_t)
        print*,'MK class is link_t'
      class default
        print*,'MK class is unknown'
      end select
    end if
  end if
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE append_link

!> =============================================================================
!> Remove a link from the list, taking special care with the head and tail
!> =============================================================================
SUBROUTINE remove_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%remove_link '//self%name)
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
  logical:: is_clear
  !.............................................................................
  call trace_begin('list_t%check_nbors '//self%name)
  task => link%task                                   ! main task
  n_check=0; n_ready=0; n_blocked=0; n_cleared=0; n_remain=0
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task                               ! neighbor task
    n_check = n_check+1
    is_clear = nbtask%is_clear(bits%ready)
    if (is_clear) then                                ! needs checking?
      n_blocked = n_blocked+1
      call self%check_ready (nbor%link)               ! pointer back
      if (nbtask%is_set(bits%ready)) then
        n_cleared = n_cleared+1
      else
        n_remain = n_remain+1
      end if
      if (nbtask%id==id_debug) print*, id_debug, 'mk checked by:', task%id, nbtask%is_set(bits%ready)
    else
      n_ready = n_ready+1
    end if
    nbor => nbor%next                                 ! next nbor
  end do
  call self%check_ready (link)                        ! finally check link task
  if (task%id==id_debug) print*, id_debug, 'mk checked self:', task%id, task%is_set(bits%ready)
  !self%n_nbor = self%n_nbor + 1
  if (io%verbose>2) print '(a,i6,l3,12(2x,a,i4))', &
    'task:',task%id, task%is_set(bits%ready), 'level:',task%level, &
    'neighbors:',n_check, 'blocked:',n_blocked, 'ready:', n_ready, &
    'cleared:',n_cleared, 'remain:', n_remain
  if (io%verbose>3) call self%print
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
  link => self%head
  do while (associated (link))                        ! keep going until end
    call self%check_ready (link, 1)                   ! link ready?
    link => link%next
  end do
  call trace_end
END SUBROUTINE check_all

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE check_ready (self, link, set_debug)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, nbtask
  class(link_t), pointer:: nbor, scan
  integer:: set_debug, debug=0
  optional:: set_debug
  logical:: ahead, is_clear, do_download
  !.............................................................................
  call trace_begin('list_t%check_ready '//self%name,2)
  if (present(set_debug)) then
    debug = set_debug
  else
    debug = 0
  end if
  !$omp critical (list_cr)
  do_download = .false.
  call self%qconsistency(1)
  task => link%task                                     ! already checked?
  if (task%is_clear(bits%ready)) then
    nbor => link%nbor                                   ! start on nbor list
    do while (associated (nbor))                        ! keep going until end
      nbtask => nbor%task
      ahead = nbtask%is_ahead_of(task)
      if (task%id==id_debug.or.io%verbose>3) &
        print*, task%id, 'checking nbor', &
        nbtask%id, ahead, nbtask%time-task%time
      if (ahead) then
        nbor => nbor%next                               ! next nbor
      else
        goto 9
      end if
    end do
    call task%set(bits%ready)                           ! mark as checked
    call self%qconsistency(2)
    call self%queue_by_time (link)                      ! add task to queue
    call self%qconsistency(3)
    do_download = .true.
  end if
9 continue
  !$omp end critical (list_cr)
  if (do_download) then
    !-----------------------------------------------------------------------------
    ! Prepare for update, by pro-actively downloading data that it depends on.
    ! By doing it immediately, the guard zone data becomes available to neighbors
    !-----------------------------------------------------------------------------
    task => link%task
    select type (task)
    class is (patch_t)
      call task%check_density ("before download")
      call download%download_link (link, all_cells=.false.)
      call task%check_density ("after download")
    end select
  end if
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
  !$omp critical (list_cr)
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
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE consistency

!> =============================================================================
!> Check if the number of queued tasks is consistent; i.e., of the number of
!> tasks linked at self%queue is the same as nq
!> =============================================================================
SUBROUTINE qconsistency (self, i)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  real(8):: time
  !.............................................................................
  !if (self%debug < 1) return
  call trace_begin ('list_t%qconsistency',0)
  n = 0
  link => self%queue
  if (associated(link)) then
    time = link%task%time
    link => link%next_time
    n = n+1
  end if
  do while (associated(link))
    n = n+1
    if (n > self%nq) exit
    if (link%task%time < time) then
      print*,'ERROR: queue out of order'
      exit
    end if
    time = link%task%time
    link => link%next_time
  end do
  if (n /= self%nq) then
    print *, 'ERROR: qconsistency:', n, self%nq, i
    n = 0
    link => self%queue
    do while (associated(link))
      n = n+1
      print*, n, link%task%id, link%task%time
      link => link%next_time
      if (n > self%nq+1) exit
    end do
    call mpi%end ('error')
  end if
  call trace_end
END SUBROUTINE qconsistency

!> =============================================================================
!> Set non_leaf flag on non-leaf patches
!> =============================================================================
SUBROUTINE init_nonleaf (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add
  !.............................................................................
  call trace_begin('list_t%init_nonleaf ')
  !$omp critical (list_cr)
  !-----------------------------------------------------------------------------
  ! First find and mark all non-leaf patches; these are patches that the link
  ! patch is fully contained inside -- we *assume* that finer level patches
  ! cover the whole parent patch!
  !-----------------------------------------------------------------------------
  link => self%head                                     ! start checking links
  do while (associated(link))
    scan => self%head                                   ! scan all links (FIXME)
    do while (associated(scan))                         ! until end
      if (.not.associated(scan,link).and. &             ! skip same patch
        .not.scan%task%is_set(bits%not_leaf)) then      ! skip non-leaf patch
        if (self%overlaps (scan%task, link%task)) then  ! if patches overlap
          if (scan%task%level==link%task%level-1) then  ! and level is one below
            if (all(abs(scan%task%position-link%task%position) < scan%task%size*0.55d0)) then
              call scan%task%set (bits%not_leaf)        ! mark non-leaf
              if (io%verbose>0) &
                print*, 'patch', link%task%id, ' parent:', scan%task%id
            end if
          end if
        end if
      end if
      scan => scan%next                                 ! continue search
    end do
    link => link%next                                   ! continue list
  end do
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE init_nonleaf

!> =============================================================================
!> The very first initialization of the neighbor lists.  After this, it can be
!> maintained without a global search.
!> =============================================================================
SUBROUTINE init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add
  !.............................................................................
  call trace_begin('list_t%init_nbors ')
  !$omp critical (list_cr)
  nullify (link%nbor)                                 ! initially no nbor
  n_add = 0
  scan => self%head                                   ! scan all links (FIXME)
  do while (associated(scan))                         ! until end
    if (.not.associated(scan,link) .and. associated(scan%task)) then
      if (self%overlaps (scan%task, link%task)) then  ! if its task overlaps
        write (io%log_unit, '(3(i6,2(3f8.4,2x)),l4)') &
          scan%task%id, scan%task%position, scan%task%size, &
          link%task%id, link%task%position, link%task%size, &
          0, self%distance (scan%task, link%task), self%size, &
          self%overlaps (scan%task, link%task)
        allocate (nbor)
        nbor%task => scan%task
        nbor%link => scan
        call link%add_nbor_by_rank (nbor)
        n_add = n_add+1
        if (nbor%task%id==id_debug) print*, id_debug, 'mk is nbor:', link%task%id
        if (link%task%id==id_debug) print*, id_debug, 'mk has nbor:', nbor%task%id
      end if
    end if
    scan => scan%next                                 ! continue search
  end do
  link%task%n_nbors = n_add                           ! remember n_nbors
  if (io%master .and. io%verbose>1) print*, &
    'added ', n_add, ' neighbors to', link%task%id, associated(link%nbor)
  if (io%verbose>1) write (io%log_unit,*) &
    'added ', n_add, ' neighbors to', link%task%id, associated(link%nbor)
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE init_nbors

!> =============================================================================
!> The very first initialization of the neighbor lists.  After this, it can be
!> maintained without a global search.
!> =============================================================================
SUBROUTINE init_all_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin('list_t%init_all_nbors')
  call init_nonleaf (self)                              ! set non_leaf bit
  link => self%head                                     ! start checking links
  do while (associated(link))
    call self%init_nbors (link)
    link => link%next                                   ! continue list
  end do
  call trace_end
END SUBROUTINE init_all_nbors

!> =============================================================================
!> Clear and then set the virtual and boundary bits
!> =============================================================================
SUBROUTINE reset_status (self)
  class(list_t):: self
  class(link_t), pointer:: link, nbor
  integer:: n_internal, n_boundary, n_virtual
  !-----------------------------------------------------------------------------
  ! Clear virtual and boundary bits
  !-----------------------------------------------------------------------------
  write (io%log_unit,*) 'reset_status'
  n_internal = 0
  n_boundary = 0
  n_virtual = 0
  !$omp critical (list_cr)
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%virtual+bits%boundary)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Set virtual bit
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    n_internal = n_internal + 1
    if (link%task%rank /= mpi%rank) then
      call link%task%set(bits%virtual)
      n_virtual  = n_virtual  + 1
      n_internal = n_internal - 1
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Set boundary bit
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    nbor => link%nbor
    do while (associated(nbor))
      if (link%task%rank==mpi%rank .and. nbor%task%is_set(bits%virtual)) then
        call link%task%set(bits%boundary)
      end if
      write (io%log_unit,'(2(a,i7,i4,2l4,5x,a,3i4,10x))') &
        'link:', link%task%id, link%task%rank, &
                 link%task%is_set(bits%boundary), &
                 link%task%is_set(bits%virtual), &
        'pos:', self%intpos(link%task), &
        'nbor:', nbor%task%id, nbor%task%rank, &
                 nbor%task%is_set(bits%boundary), &
                 nbor%task%is_set(bits%virtual), &
        'pos:', self%intpos(nbor%task,link%task)
      nbor => nbor%next
    end do
    if (link%task%is_set(bits%boundary)) then
      n_boundary = n_boundary + 1
      n_internal = n_internal - 1
    end if
    link => link%next
  end do
  write (io%log_unit,'(4(a,i7,5x))') &
    'rank:',mpi%rank, &
    'n_internal:', n_internal, &
    'n_boundary:', n_boundary, &
    'n_virtual:' , n_virtual
  flush (io%log_unit)
  !$omp end critical (list_cr)
END SUBROUTINE reset_status

!> =============================================================================
!> Remove parent patches, which have higher level patches inside
!> =============================================================================
SUBROUTINE remove_parents (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan, parent
  !.............................................................................
  call trace_begin('list_t%remove_parents')
  !$omp critical (list_cr)
  link => self%head                                     ! start checking links
  do while (associated(link))
    scan => self%head                                   ! scan all links (FIXME)
    do while (associated(scan))                         ! until end
      !---------------------------------------------------------------------
      ! Find parents, as patches on the next lower level that have the
      ! current patch position inside
      !---------------------------------------------------------------------
      nullify (parent)
      if (scan%task%level==link%task%level-1) then
        if (all(abs(scan%task%position-link%task%position) < scan%task%size*0.55d0)) then
          parent => scan
        end if
      end if
      scan => scan%next                                 ! continue search
      if (associated(parent)) then
        print '(2(1x,a,i8,5x))', &
          'remove parent task', parent%task%id, 'level', parent%task%level
        call self%remove_link (parent)
      end if
    end do
    link => link%next                                   ! continue list
  end do
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE remove_parents

!> =============================================================================
!> =============================================================================
FUNCTION distance (self, task1, task2)
  class(list_t):: self
  class(task_t), pointer:: task1, task2
  real(8):: distance(3), size(3)
  !.............................................................................
  size = self%size
  distance = task1%position - task2%position
  distance = abs(modulo(distance + 0.5*size, size) - 0.5*size)
END FUNCTION distance

!> =============================================================================
!> =============================================================================
FUNCTION intpos (self, task1, task2) RESULT (out)
  class(list_t):: self
  class(task_t), pointer:: task1
  class(task_t), pointer, optional:: task2
  real:: dist(3)
  integer:: out(3)
  !.............................................................................
  if (present(task2)) then
    dist = task1%position - task2%position
    dist = modulo (dist + 0.5*self%size, self%size) - 0.5*self%size
    out = floor(dist/task1%size + 0.5)
  else
    out = floor(task1%position/task1%size + 0.5)
  end if
END FUNCTION intpos

!> =============================================================================
!> =============================================================================
FUNCTION overlaps (self, task1, task2, label)
  class(list_t):: self
  class(task_t), pointer:: task1, task2
  character(len=*), optional:: label
  logical:: overlaps
  real(8):: dist(3)
  !.............................................................................
  call trace_begin ('task_list_t%overlap', 3)
  if (self%face_nbors) then
    overlaps = sum(self%distance(task1,task2)/(task1%size+task2%size)) &
             < 0.50000001_8
  else
    overlaps = all(self%distance(task1,task2) < 0.50000001_8*(task1%size+task2%size))
  end if
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
  !$omp critical (list_cr)
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
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE append_list

!> =============================================================================
SUBROUTINE remove_task (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link, prev
  !.............................................................................
  call trace_begin('list_t%remove_task')
  !$omp critical (list_cr)
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
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE remove_task

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
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  ! Refuse to queue virtual tasks and finished tasks
  !-----------------------------------------------------------------------------
  if (this%task%is_set(bits%virtual) .or. this%task%has_finished()) return
  call trace_begin ('list_t%queue_by_time '//self%name)
  !-----------------------------------------------------------------------------
  ! Lock the queue while inserting the task
  !-----------------------------------------------------------------------------
  call this%task%set(bits%ready)                       ! mark as checked
  if (io%verbose>2) print'(1x,a,i6,a,f10.5)', 'list_t%queue_by_time: adding task', &
    this%task%id, ' at time', this%task%time
  nullify (prev)
  next => self%queue
  do while (associated(next))
    if (associated(next%task, this%task)) then
      print*, omp_mythread(), ' WARNING: task', this%task%id, ' is already on ready queue'
    else if (next%task%time > this%task%time) then
      this%next_time => next
      if (associated(prev)) then
        prev%next_time => this
      else
        self%queue => this
      end if
      self%nq = self%nq+1
      call trace_end
      return
    end if
    prev => next
    next => next%next_time
  end do
  self%nq = self%nq+1
  if (associated(prev)) then
    prev%next_time => this
  else
    self%queue => this
  end if
  nullify (this%next_time)
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
  !$omp critical (list_cr)
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
  !$omp end critical (list_cr)
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
  !$omp critical (list_cr)
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
  !$omp end critical (list_cr)
END SUBROUTINE check_queue

!> =============================================================================
!> Insert 'new' task before 'old' task
!> =============================================================================
SUBROUTINE insert (self, old, new)
  class(list_t):: self
  class(link_t), pointer:: old, new
  !.............................................................................
  call trace_begin ('list_t%insert '//self%name)
  !$omp critical (list_cr)
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    self%head => new
  end if
  new%next => old
  new%next%prev => new
  self%n = self%n + 1
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE insert

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_queue (self)
  class (list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  !$omp critical (list_cr)
  link => self%queue
  print'(a,2i4,":",$)', 'ready_queue', self%nq, self%n
  do while (associated(link))
    print '(i4,$)', link%task%id
    link => link%next_time
  end do
  !$omp end critical (list_cr)
  if (io%master) print*,''
END SUBROUTINE print_queue

!> =============================================================================
!> Print the queue until and including a given task
!> =============================================================================
SUBROUTINE print_queue_until (self, task)
  class(list_t):: self
  class(task_t):: task
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !$omp critical (list_cr)
  link => self%queue
  if (io%master) print*,'print_queue_until:'
  do while (associated(link))
    print '("queue: id, t0, t, t0-t =",i7,1p,3e12.3)', &
      link%task%id, task%time, link%task%time, task%time-link%task%time
    if (link%task%id==task%id) exit
    link => link%next_time
  end do
  !$omp end critical (list_cr)
END SUBROUTINE print_queue_until

!> =============================================================================
!> Print the queue
!> =============================================================================
SUBROUTINE print_queue_times (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !$omp critical (list_cr)
  link => self%queue
  if (present(label)) then
     print*,'print_queue_times: '//trim(label)
  else
     print*,'print_queue_times:'
  end if
  do while (associated(link))
    print '("queue: id, time =",i7,1p,g14.5,l5)', link%task%id, link%task%time, &
      link%task%is_set(bits%ready)
    link => link%next_time
  end do
  !$omp end critical (list_cr)
  print*,'done'
END SUBROUTINE print_queue_times

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
  if (io%master) print*,'list_t%print: ',trim(self%name), self%n
  !$omp critical (list_cr)
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
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE print_list

END MODULE list_mod
