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
  USE mpi_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE bits_mod
  USE trace_mod
  USE refine_mod
  USE experiment_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
    real(8):: size(3) = 1.0_8
    integer:: n_behind=0, n_check=0, n_ready=0, n_nbor=0
    integer:: debug=0
    !---------------------------------------------------------------------------
    ! from patch_list_mod
    !---------------------------------------------------------------------------
    integer:: dims(3)
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
    procedure:: statistics
    procedure:: append_list
    !---------------------------------------------------------------------------
    ! from patch_list_mod
    !---------------------------------------------------------------------------
    procedure:: init_cartesian
    procedure:: remove_task
    procedure:: split
    procedure:: add_by_quality
    procedure:: queue_by_time
    procedure:: output
    procedure:: check_queue
    procedure:: print_queue
    procedure:: print
  end type
  integer, parameter:: mnbor=100
  type(experiment_t), pointer, dimension(:,:,:):: patches
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
  call self%init_cartesian
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
  if (self%nq < 1) then
    n1 = self%nq
    call self%check_all
    print *,'WARNING: check_all (should not be needed), found', self%nq
  end if
  head => self%queue                                    ! queue start
  if (associated(head)) then
    self%queue => self%queue%next_time                  ! update the link
    self%nq = self%nq-1
    call self%check_mpi                                 ! check incoming MPI
    task => head%task                                   ! head task
    call head%download (all_cells=.false.)              ! guard cells
    !if (refine%needed (task)) then                      ! if refinement needed
    !  call self%split (head)                            ! refined tasks go back to queue
    !end if
    call task%output(self%name)                         ! check output
    call task%update                                    ! update the task
    call task%clear(bits%ready)                         ! needs new check
    call task%output(self%name)                         ! check output
    call self%check_nbors (head)                        ! any nbor ready?
    if (task%time+0.5*task%dtime > task%end_time) &     ! finished:
      call self%remove_task (task)                      ! remove from self
  else
    print*,'ERROR: no self%queue', self%nq, self%n
    call mpi%abort ('ERROR: ready queue empty')
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
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  integer:: nbors(mnbor), inbor
  integer:: n_check, n_ready, n_blocked, n_cleared, n_remain
  !.............................................................................
  call trace_begin('task_list_t%check_nbors '//self%name)
  task => link%task                                   ! main task
  inbor = 0
  n_check=0; n_ready=0; n_blocked=0; n_cleared=0; n_remain=0
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    nbtask => nbor%task                               ! neighbor task
    inbor = min(inbor+1,mnbor)                        ! index to diagnostic array
    n_check = n_check+1
    if (nbtask%is_clear(bits%ready)) then             ! needs checking?
      n_blocked = n_blocked+1
      nbors(inbor) = nbtask%id
      call self%check_ready (nbor%link)               ! pointer back
      if (nbtask%is_set(bits%ready)) then
        n_cleared = n_cleared+1
      else
        n_remain = n_remain+1
      end if
    else
      n_ready = n_ready+1
      nbors(inbor) = -nbtask%id
    end if
    nbor => nbor%next                                 ! next nbor
  end do
  if (self%debug>0) &
    print'(i7,a,i4,2x,26i4)', task%id, ' nbors:', inbor, nbors(1:inbor)
  call self%check_ready (link)                        ! finally check link task
  self%n_nbor = self%n_nbor + 1
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
  class(task_list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin('task_list_t%check_all '//self%name)
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
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, nbtask
  class(link_t), pointer:: nbor
  integer:: set_debug, debug, nbors(mnbor), inbor
  optional:: set_debug
  !.............................................................................
  if (present(set_debug)) then
    debug = set_debug
  else
    debug = 0
  end if
  call self%qconsistency(1)
  task => link%task                                     ! already checked?
  if (task%is_set(bits%ready)) return
  call trace_begin('task_list_t%check_ready '//self%name,2)
  nbor => link%nbor                                     ! start on nbor list
  inbor = 0
  do while (associated (nbor))                          ! keep going until end
    nbtask => nbor%task
    if (task%ahead_of(nbtask)) then
      self%n_check = self%n_check+1
      inbor = min(inbor+1,mnbor)
      nbors(inbor) = nbtask%id
      if (self%debug==1) &
        print*, task%id, ' check_ready failed on task', nbtask%id
      if (self%debug>2) &
        print '(i4,a,i3,2x,30i3)', task%id, ' failed:', inbor, nbors(1:inbor)
      call trace_end
      return                                            ! give up -- do nothing
    else
      inbor = min(inbor+1,mnbor)
      nbors(inbor) = nbtask%id
!      print*,'nbor task', nbtask%id, nbtask%time
    end if
    !nbor => nbor%nbor                                  ! next nbor
    nbor => nbor%next                                   ! next nbor
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
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  if (self%debug < 2) return
  call trace_begin ('task_list_t%consitency',1)
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
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  if (self%debug < 2) return
  call trace_begin ('task_list_t%qconsistency',0)
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
  class(task_list_t):: self
  type(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add = 0
  !.............................................................................
  call trace_begin('task_list_t%init_nbors ')
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
  class(task_list_t):: self
  class(task_t), pointer:: task1, task2
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
  real(8):: dist(3)
  !.............................................................................
  call trace_begin ('task_lis_t%overlap', 3)
  dist = self%distance(task1,task2)
  overlaps = (all(self%distance(task1,task2) < 0.501*(task1%size+task2%size)))
  if (self%debug>3) &
    print'(a,2i4,2(3f8.2,3x),l3)', 'overlap:', task1%id, task2%id, &
     task1%position/task1%size-task2%position/task2%size, dist, overlaps
  call trace_end
END FUNCTION overlaps

!> =============================================================================
!> =============================================================================
SUBROUTINE statistics (self)
  class(task_list_t):: self
  print '(a,3(5x,a,i6))', 'output_experiment:', &
    'n_ready:', self%n_ready, &
    'n_check:', self%n_check, &
    'n_nbor:' , self%n_nbor
END SUBROUTINE statistics

!> =============================================================================
!> =============================================================================
SUBROUTINE append_list (self, task_list)
  class(task_list_t):: self
  class(task_list_t), pointer:: task_list
  type(link_t), pointer:: link
  integer:: n
  !-----------------------------------------------------------------------------
  call trace_begin ('task_list_t%append_list')
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
  call trace_end
END SUBROUTINE append_list

!> =============================================================================
!> Distribute patches in a Cartesian arrangement in a box
!> =============================================================================
SUBROUTINE init_cartesian (self)
  class (task_list_t):: self
  integer:: dims(3) = [4,4,1]
  real:: origin(3) = [-0.5,-0.5,-0.5]
  real:: size(3) = [1.0,1.0,1.0]
  real:: position(3)
  type(experiment_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link
  class(task_t), pointer:: task
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  namelist /box_params/ size, dims, origin
  !.............................................................................
  call trace_begin ('task_list_t%init_cartesian')
  rewind (io%input); read (io%input, box_params)
  if (io%master) write (*, box_params)
  self%dims = dims
  self%size = size
  allocate (patches(self%dims(1),self%dims(2),self%dims(3)))
  do i=1,self%dims(1)
  do j=1,self%dims(2)
  do k=1,self%dims(3)
    !---------------------------------------------------------------------------
    ! Initialize each patch, via patch_t:init
    !---------------------------------------------------------------------------
    patch => patches(i,j,k)
    patch%box = self%size
    patch%position = ([i,j,k]-0.5_8)*self%size/self%dims + origin
    patch%size = self%size/self%dims
    call patch%init
    if (io%verbose>1) &
      print'(a,i6,2(2x,3f10.6))', 'init_cartesian: id, size, pos =', &
        patch%id, patch%size, patch%position
    patch%i = i; patch%j = j; patch%k = k
    if (io%verbose > 2) &
      print'(i5,3(3f9.4,3x))', patch%id, patch%position, patch%size
    !---------------------------------------------------------------------------
    ! Add to a list of links, with pointers from patch to patch
    !---------------------------------------------------------------------------
    task => patches(i,j,k)
    call self%append (task)
  end do
  end do
  end do
  call trace_end
END SUBROUTINE init_cartesian

!> =============================================================================
SUBROUTINE remove_task (self, task)
  class(task_list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link, prev
  !.............................................................................
  call trace_begin('task_list_t%remove_task')
  link => self%head
  nullify (prev)
  do while (associated(link))
    if (link%task%id == task%id) then
      if (associated(prev)) then
        prev%next => link%next
        print*,'removed task', task%id, '  after', prev%task%id
      else
        self%head => link%next
        print*,'removed task at head', task%id
      end if
      self%n = self%n - 1
      call trace_end
      return
    end if
    prev => link
    link => link%next
  end do
  print*,'ERROR: failed to remove from task list, task', task%id
  call trace_end
END SUBROUTINE remove_task

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
  type(link_t), pointer:: tail
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
        if (this_patch%time > 0.0) call link%download (all_cells=.true.)
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
  call trace_begin('task_list_t%make_new_nbors', 1)
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
  call trace_begin('task_list_t%make_new_nbor', 1)
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
  type(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_quality')
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
        return
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
  type(link_t), pointer:: nbor
  class(patch_t), pointer:: p1, p2
  !> ===========================================================================
  !> For each task link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('task_list_t%search_nbors_by_overlap', 3)
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
  class(task_list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%queue_by_time '//self%name)
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
  class(task_list_t):: self
  class(link_t), pointer:: this
  type(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_by_quality')
  nullify (prev)
  next => self%head
  do while (associated(next))
    if (associated(next%task,this%task)) then
      print *,'ERROR: Trying to add task', this%task%id, ' which is already there'
      call trace_end
      return
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
      return
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
  call trace_end
END SUBROUTINE add_by_quality

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE output (self)
  class (task_list_t):: self
  type(link_t) , pointer:: link, nbor
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
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE check_queue (self)
  class(task_list_t):: self
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
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_tasks (self)
  class (task_list_t):: self
  type(link_t) , pointer:: link, nbor
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
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_queue (self)
  class (task_list_t):: self
  type(link_t) , pointer:: link, nbor
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
!> Print a table with task id and all neighbors
!> =============================================================================
SUBROUTINE print (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  print*, label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
  call self%print_queue
END SUBROUTINE print

END MODULE
