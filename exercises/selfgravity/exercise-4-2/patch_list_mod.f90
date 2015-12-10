!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Handles patch_list, with links that have pointers to the next link, to a patch,
!> and to a list with neighbors to the patch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_list_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patches_mod
  USE timestep_mod
  USE refine_mod
  implicit none
  private
  !.............................................................................
  type, public:: link_t
    type(link_t), pointer:: next => null()            ! head => head%next => head%next%next
    class(link_t), pointer:: nbors => null()          ! nbors => nbors%next => ...
    class(link_t), pointer:: parent => null()
    type(link_t), pointer:: next_time => null()
    type(timestep_t), pointer:: patch => null()
    integer:: n_nbors
  contains
    procedure:: download
    procedure:: add_nbor_by_quality
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    procedure:: print_nbors
  end type
  !.............................................................................
  type, public:: patch_list_t
    type(link_t), pointer:: head=>null(), tail=>null(), queue=>null()
    real:: size(3)=1.0
    integer:: dims(3)=[4,4,1]
    character(len=120):: name = 'patch list'
    integer:: n=0, nq=0, n_nbor=0, n_check=0, n_ready=0, n_behind=0
    integer:: debug=0
  contains
    procedure:: init
    procedure:: init_cartesian
    procedure:: update
    procedure:: remove_patch
    procedure:: split
    procedure:: init_nbors_by_overlap
    procedure:: init_nbor_by_overlap
    procedure:: check_nbors
    procedure:: check_ready
    procedure:: check_all
    procedure:: add_by_quality
    procedure:: queue_by_time
    procedure:: output
    procedure:: check_queue
    procedure:: print
  end type
  type(timestep_t), pointer, dimension(:,:,:):: patches
  integer, parameter:: mnbor=200
CONTAINS
!> =============================================================================
SUBROUTINE init (self, name)
  class (patch_list_t):: self
  character(len=*):: name
  !.............................................................................
  call trace_begin('patch_list_t%init')
  call io%init
  call refine%init
  self%name = name
  call init_cartesian (self)
  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE init_cartesian (self)
  class (patch_list_t):: self
  integer:: dims(3) = [4,4,1]
  real:: origin(3) = [-0.5,-0.5,-0.5]
  real:: size(3) = [1.0,1.0,1.0]
  real:: position(3)
  type(timestep_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  namelist /box_params/ size, dims, origin
  !.............................................................................
  call trace_begin ('patch_list_t%init_cartesian')
  rewind (io%input); read (io%input, box_params)
  if (io%master) write (*, box_params)
  self%dims = dims
  self%size = size
  allocate (patches(self%dims(1),self%dims(2),self%dims(3)))
  do i=1,self%dims(1)
  do j=1,self%dims(2)
  do k=1,self%dims(3)
    !---------------------------------------------------------------------------
    ! Initialize each patch, via timestep_t:init -> pde_t:init -> patch_t:init
    !---------------------------------------------------------------------------
    patch => patches(i,j,k)
    position = ([i,j,k]-0.5_8)*self%size/self%dims + origin
    patch%box = self%size
    call patch%init (position, self%size/self%dims)
    if (io%verbose>1) &
      print*,'init_cartesian: id, size, pos =', patch%id, patch%size, patch%position
    patch%i = i; patch%j = j; patch%k = k
    if (io%verbose > 2) &
      print'(i5,3(3f9.4,3x))', patch%id, patch%position, patch%size
    !---------------------------------------------------------------------------
    ! Add to a list of links, with pointers from patch to patch
    !---------------------------------------------------------------------------
    self%n = self%n + 1
    allocate (link)
    link%patch => patches(i,j,k)
    if (associated(self%head)) then
      prev%next => link
    else
      self%head => link
    end if
    prev => link
  end do
  end do
  end do
  call self%init_nbors_by_overlap
  self%tail => link
  call trace_end
END SUBROUTINE init_cartesian

!> =============================================================================
!> Update the state of the patch list, taking the steps necessary to update the
!> head patch on the ready_queue, and check for consequences
!> =============================================================================
SUBROUTINE update (self)
  class(patch_list_t):: self
  !.............................................................................
  class(link_t), pointer:: head
  class(timestep_t), pointer:: patch
  integer:: n1
  logical:: end_game=.false.
  !.............................................................................
  call trace_begin('patch_list_t%update')
!  if (end_game) call print_lists (self)
  if (self%nq < 1) then
    n1 = self%nq
    call self%check_all
    !print *,'WARNING: check_all (should not be needed), found', self%nq, n1
  end if
  call self%check_queue
  head => self%queue                                    ! queue start
  if (associated(head)) then
    self%queue => self%queue%next_time                  ! update the link
    self%nq = self%nq-1
    call head%download                                  ! guard cells & interior
    patch => head%patch                                 ! head patch
    if (io%verbose > 1) then
      if (associated(self%queue)) then
        print*,'patch_list_t%update: removed from queue =', &
          patch%id, ' new head=', self%queue%patch%id
      else
        print*,'patch_list_t%update: removed from queue =', &
          patch%id, ' no new head'
      end if
    end if
    call patch%output                                   ! check output
    if (refine%needed (patch)) then                     ! if refinement needed
      call self%split (head)                            ! refined patchs go back to queue
      if (patch%level > refine%levelmax) then
        call self%output
        print*,'stopping after refine reached levelmax =', refine%levelmax
        stop
      end if
    end if
    call patch%update                                 ! update the patch
    call patch%output                                 ! check output
    call self%check_nbors (head)                      ! any nbor ready?
    if (patch%time > patch%end_time) then             ! finished:
      end_game = .true.
      call self%remove_patch (patch)                  ! remove from self
    else                                              ! not finished:
      call patch%clear(bits%ready)                    ! needs new check
    end if
  else
    print*,'stopping: no more tasks', self%nq, self%n
    stop
  end if
  call print_lists(self)
  call trace_end
END SUBROUTINE update

!> =============================================================================
SUBROUTINE remove_patch (self, patch)
  class(patch_list_t):: self
  class(link_t), pointer:: link, prev
  class(patch_t):: patch
  !.............................................................................
  call trace_begin('patch_list_t%remove_patch')
  link => self%head
  nullify (prev)
  do while (associated(link))
    if (link%patch%id == patch%id) then
      if (associated(prev)) then
        prev%next => link%next
        print*,'removed patch', patch%id, '  after', prev%patch%id
      else
        self%head => link%next
        print*,'removed patch at head', patch%id
      end if
      self%n = self%n - 1
      call trace_end
      return
    end if
    prev => link
    link => link%next
  end do
  print*,'ERROR: failed to remove from task list, patch', patch%id
  call trace_end
END SUBROUTINE remove_patch

!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split (self, this)
  class(patch_list_t):: self
  class(link_t), pointer:: this, link
  integer:: nsplit(3)
  type(timestep_t), pointer:: patch
  type(link_t), pointer:: tail
  integer dir, i, ix, iy, iz
  real:: pos(3)
  !.............................................................................
  call trace_begin('patch_list_t%split', 1)
  nsplit = msplit
  where (this%patch%size==0 .or. this%patch%n==1)
    nsplit = 1
  end where
  do ix=1,nsplit(1)
    pos(1) = ix-(nsplit(1)+1)*0.5
    do iy=1,nsplit(2)
      pos(2) = iy-(nsplit(2)+1)*0.5
      do iz=1,nsplit(3)
        pos(3) = iz-(nsplit(3)+1)*0.5
        !------------------------------------------------------------------------
        ! Allocate a new patch, half the size of original one, and adjust its
        ! position in each direction with "pos", in units of the new size.
        !------------------------------------------------------------------------
        allocate (patch, source=this%patch)
        patch%mem_allocated = .false.
        patch%size = this%patch%size/nsplit
        patch%position = this%patch%position + pos*patch%size
        nullify (patch%mem)
        call patch%init (patch%position, patch%size)    ! patch, for a new id
        patch%out_next = patch%time                     ! force output now
        patch%iout = this%patch%iout                    ! copy counter
        !------------------------------------------------------------------------
        ! Append links to the new patch in the patch list, sorted by quality
        !------------------------------------------------------------------------
        allocate (link)
        link%patch => patch
        link%parent => this                     ! link is the link to the new patch
        call self%add_by_quality (link)         ! which is added to the task list
        call link%make_new_nbors (link, this)   ! make a new nbors list for link
        !------------------------------------------------------------------------
        ! Download values from parent and neighbors, except at time==0.0, to
        ! allow recursive refinement there
        !------------------------------------------------------------------------
        if (this%patch%time > 0.0) call link%download
        !------------------------------------------------------------------------
        ! Insert the new patch at the head of the ready queue
        !------------------------------------------------------------------------
!        link%next_time => self%queue
!        self%queue => link
!        self%nq = self%nq+1
        if (io%verbose>0) then
          print 1,'task_list_t%split:  patch, iout, time, size, origin =', &
            patch%id, patch%iout, patch%time, patch%size, patch%origin
          1 format(a,2i6,f10.6,2(2x,3f8.4))
          print 1,'task_list_t%split: parent, iout, time, size, origin =', &
            this%patch%id, this%patch%iout, this%patch%time, this%patch%size, this%patch%origin
        end if
      end do
    end do
  end do
  call self%check_queue
  call this%patch%set(bits%not_leaf)
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
  !> ...........................................................................
  call trace_begin('patch_list_t%make_new_nbors', 1)
  call self%make_new_nbor (selfp, parent)         ! make me a neighbor of my parent
  nbor => parent%nbors                            ! start on parent's nbors chain
  do while (associated(nbor))
    if (.not.associated(self%patch,nbor%patch) &
      .and. self%patch%overlaps (nbor%patch)) &   ! if my patch overlaps with that patch
      call self%make_new_nbor (selfp, nbor%nbors) ! make a neighbor (and make me one)
    nbor => nbor%next                             ! continue on parents nbors chain
  end do
  call trace_end
END SUBROUTINE make_new_nbors

!> =============================================================================
!> Add a link to a new neighbor to my nbors chain, and add myself into that
!> neighbor's neighbor chain.  An nbor chain link has three pointers: a pointer
!> to its patch, a pointer to the next link in the chain, and a pointer back to
!> the patch list link to the patch.  We need to make and add two now nbor chain
!> links; one for our own nbors chain, and one for that's nbors chain.  Both
!> are added in quality order
!> =============================================================================
SUBROUTINE make_new_nbor (self, selfp, that)
  class(link_t):: self
  class(link_t), pointer:: selfp, that
  class(link_t), pointer:: link1, link2
  !> ...........................................................................
  call trace_begin('patch_list_t%make_new_nbor', 1)
  allocate (link1)                              ! for our own nbors chain
  link1%patch => that%patch                     ! the link points to that patch
  link1%nbors => that                           ! link to the nbor link
  link1%parent => that%parent                   ! direct link to patch parent
  call self%add_nbor_by_quality (link1)         ! add to that chain in quality order
  allocate (link2)                              ! for that's nbors chain
  link2%patch => self%patch                     ! that's new link points to my patch
  link2%nbors => selfp                          ! pointer back to self
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
  next => self%nbors
  do while (associated(next))
    if (.not.associated(next%patch,this%patch)) then
      if (next%patch%level >= this%patch%level) then
        !-------------------------------------------------------------------------
        ! Found a patch we should be ahead of, so insert here
        !-------------------------------------------------------------------------
        this%next => next
        if (associated(prev)) then
          prev%next => this
        else
          self%nbors => this
        end if
        self%n_nbors = self%n_nbors+1
        call trace_end
        return
      end if
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This patch has the highest quality, so append it
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => this
    nullify (this%next)
  else
    self%nbors => this
  end if
  self%n_nbors = self%n_nbors+1
  call trace_end
END SUBROUTINE add_nbor_by_quality

!> =============================================================================
!> Initialize nbors via a double loop over the linked list.
!> =============================================================================
SUBROUTINE init_nbors_by_overlap (self)
  class(patch_list_t):: self
  type(link_t) , pointer:: prev, link
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbors
  !> ===========================================================================
  call trace_begin('patch_list_t%init_nbors_by_overlap', 1)
  link => self%head
  do while (associated(link))
    call self%init_nbor_by_overlap (link)
    link => link%next
  end do
  call self%print (' patches in list, with nbors by overlap:')
  call trace_end
END SUBROUTINE init_nbors_by_overlap

!> =============================================================================
!> Initialize the nbor list of a single link
!> =============================================================================
SUBROUTINE init_nbor_by_overlap (self, link1)
  class(patch_list_t):: self
  type(link_t) , pointer:: prev, link1, link2, nbor
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('patch_list_t%init_nbor_by_overlap', 2)
  link1%n_nbors = 0
  nullify (link1%nbors)
  link2 => self%head
  if (io%verbose>3) print*,' ids, overlap:', &
    link2%patch%id, link1%patch%id, link1%patch%overlaps(link2%patch), &
    all(link1%patch%size==link1%patch%box)
  do while (associated(link2))
    if (link1%patch%id /= link2%patch%id) then
   !if (link1%patch%id /= link2%patch%id .or. all(link1%patch%size==link1%patch%box)) then
      if (link1%patch%overlaps(link2%patch)) then
        allocate (nbor)
        nbor%nbors => link2                   ! pointer back to link2
        link1%n_nbors = link1%n_nbors+1
        nbor%patch => link2%patch
        if (associated(link1%nbors)) then
          prev%next => nbor
        else
          link1%nbors => nbor
        end if
        prev => nbor
        if (io%verbose>=3)print*,'overlaps', link1%patch%id, link2%patch%id, ' T'
      else
        if (io%verbose>=3)print*,'overlaps', link1%patch%id, link2%patch%id, ' F'
      end if
    end if
    link2 => link2%next
  end do
  if (io%verbose > 1) print*, 'init_nbor_by_overlap, nbors found:', link1%n_nbors
  call trace_end
END SUBROUTINE init_nbor_by_overlap
!> =============================================================================
!> Update the nbor chains in the self list (possibly already existing) by,
!> for each link in that list, searching from the given start point for new
!> neighbors.
!> =============================================================================
SUBROUTINE search_nbors_by_overlap (self, start)
  class(link_t), target:: self
  class(link_t), pointer:: link1, link2, start, prev
  type(link_t), pointer:: nbor
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('patch_list_t%search_nbors_by_overlap', 3)
  link1 => self                                 ! start on the self list
  if (io%verbose>3) print*,' id, size, position, distance:'
  do while (associated(link1))
    prev => link1%nbors                         ! existing nbor list head (empty?)
    link2 => start                              ! start searching throuch nbor chain
    do while (associated(link2))
      if (.not.associated(link1%patch,link2%patch)) then
        if (io%verbose>3) &
          print'(2(i5,3(3f8.4,2x)),l3)', &
            link1%patch%id, link1%patch%size, link1%patch%position, link1%patch%distance(link2%patch), &
            link2%patch%id, link2%patch%size, link2%patch%position, link2%patch%distance(link1%patch), &
            link1%patch%overlaps(link2%patch)
        if (link1%patch%overlaps(link2%patch)) then
          allocate (nbor)                       ! new nbor (has no nbor)
          nbor%nbors => link2                   ! pointer back to link2
          link1%n_nbors = link1%n_nbors+1
          nbor%patch => link2%patch             ! but has a patch
          if (associated(prev)) then            ! check if we are in list
            nbor%next => prev%next              ! could be empty
            prev%next => nbor                   ! insert the new nbor
          else                                  ! nbor list was empty
            link1%nbors => nbor                 ! start with the new nbor
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
!> Among a patch and its neighbor patchs, move local patchs to ready_queue if they
!> are ready, and send patch data to non-locals.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(patch_list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: patch, nbpatch
  integer:: nbors(mnbor), inbor
  !.............................................................................
  call trace_begin('patch_list_t%check_nbors '//self%name)
  patch => link%patch                                 ! main patch
  nbor => link%nbors                                  ! first nbor
  inbor = 0
  self%n_check = 0
  self%n_ready = 0
  self%n_behind = 0
  do while (associated (nbor))                        ! keep going until end
    nbpatch => nbor%patch                             ! neighbor patch
    inbor = min(inbor+1,mnbor)
    if (nbpatch%is_clear(bits%ready)) then            ! needs checking?
      nbors(inbor) = nbpatch%id
      call self%check_ready (nbor%nbors)              ! pointer back
    else
      nbors(inbor) = -nbpatch%id
    end if
    nbor => nbor%next                                 ! next nbor
  end do
  if (self%debug>0) &
    print'(i7,a,i4,2x,26i4)', patch%id, ' nbors:', inbor, nbors(1:inbor)
  call self%check_ready (link)                        ! finally check link patch
  self%n_nbor = self%n_nbor + 1
  if (io%verbose>2) print '(9(1x,a,i6,1x))', &
    'patch',patch%id,'on level',patch%level,'has',link%n_nbors,'neighbors, checked', &
    self%n_check,'nbor-nbors, of which', self%n_ready,'are ready, and',self%n_behind,&
    'were behind'
  call trace_end
END SUBROUTINE check_nbors

!> =============================================================================
!> Among a patch and its neighbor patchs, move local patches to ready_queue if they
!> are ready, and send patch data to non-locals.
!> =============================================================================
SUBROUTINE check_all (self)
  class(patch_list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin('patch_list_t%check_all '//self%name)
  link => self%head
  do while (associated (link))                        ! keep going until end
    call self%check_ready (link, 1)             ! link ready?
    link => link%next
  end do
  call trace_end
END SUBROUTINE check_all

!> =============================================================================
!> Check if the link patch is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE check_ready (self, link, set_debug)
  class(patch_list_t):: self
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch, nbpatch
  class(link_t), pointer:: nbor
  integer:: set_debug, debug, nbors(mnbor), inbor
  optional:: set_debug
  !.............................................................................
  if (present(set_debug)) then
    debug = set_debug
  else
    debug = 0
  end if
  patch => link%patch                                   ! already checked?
  if (patch%is_set(bits%ready)) return                  ! skip if already checked
  if (patch%time > patch%end_time) return               ! don't bother if finished
  call trace_begin('patch_list_t%check_ready '//self%name,2)
  nbor => link%nbors                                    ! start on nbor list
  inbor = 0
  do while (associated (nbor))                          ! keep going until end
    nbpatch => nbor%patch
    inbor = min(inbor+1,mnbor)
    if (patch%ahead_of(nbpatch)) then
      self%n_check = self%n_check+1
      nbors(inbor) = nbpatch%id
      if (self%debug==1) &
        print*, patch%id, ' check_ready failed on patch', nbpatch%id
      if (self%debug>2) &
        print '(i4,a,i3,2x,30i3)', patch%id, ' failed:', inbor, nbors(1:inbor)
      call trace_end
      return                                            ! give up -- do nothing
    else
      self%n_behind = self%n_behind+1
      nbors(inbor) = nbpatch%id
!      print*,'nbor patch', nbpatch%id, nbpatch%time
    end if
    nbor => nbor%next                                   ! next nbor
  end do
  call self%queue_by_time (link)                        ! add patch to queue
  self%n_ready = self%n_ready+1
  if (self%debug>0) print *, patch%id, ' found by check_ready, nq =', self%nq
  call trace_end
END SUBROUTINE check_ready

!> =============================================================================
!> Add a link node into the queue chain, in ascending time order
!> =============================================================================
SUBROUTINE queue_by_time (self, this)
  class(patch_list_t):: self
  class(link_t), pointer:: this
  type(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%queue_by_time '//self%name)
  call this%patch%set(bits%ready)                       ! mark as checked
  if (io%verbose>2) print'(1x,a,i6,a,f10.5)', 'list_t%queue_by_time: adding patch', &
    this%patch%id, ' at time', this%patch%time
  nullify (prev)

  next => self%queue
  do while (associated(next))
    if (associated(next%patch, this%patch)) then
      print*,'WARNMING: patch', this%patch%id, ' is already on ready queue'
    else if (next%patch%time >= this%patch%time) then
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
!> Add a link node into the patch list, in increasing quality order
!> =============================================================================
SUBROUTINE add_by_quality (self, this)
  class(patch_list_t):: self
  class(link_t), pointer:: this
  type(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_by_quality')
  nullify (prev)
  next => self%head
  do while (associated(next))
    if (associated(next%patch,this%patch)) then
      print *,'ERROR: Trying to add patch', this%patch%id, ' which is already there'
      call trace_end
      return
    end if
    if (next%patch%level >= this%patch%level) then
      !-------------------------------------------------------------------------
      ! Found a patch we should be ahead of, so insert here
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
  ! This patch has the highest quality, so append it
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
!> Run through neighbor list and pick up guard cells
!> =============================================================================
SUBROUTINE download (self)
  implicit none
  class (link_t):: self
  class (link_t),  pointer:: nbor
  class (patch_t), pointer:: patch
  integer :: nguard, special
  !.............................................................................
  call trace_begin('link_t%download')
  nbor => self%nbors
  patch => self%patch
  nguard=0
  do while (associated (nbor))
    if (abs(self%patch%level-patch%level) < 2) then
      special = -1
      call nbor%patch%pde_t%download_to (patch, special)
    end if
    nbor => nbor%next
  end do
  if (io%verbose > 5) &
    print *, 'download: Number of neighbors with guard cells:', nguard
  call trace_end
END SUBROUTINE download

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE output (self)
  class (patch_list_t):: self
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  link => self%head
  do while (associated(link))
    call link%patch%output
    link => link%next
  end do
END SUBROUTINE output

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE check_queue (self)
  class(patch_list_t):: self
  class(link_t), pointer:: link
  real:: previous
  integer:: nq
  !.............................................................................
  link => self%queue
  previous = -1.0
  nq = 0
  do while (associated(link))
    if (link%patch%time < previous) &
      print *,'ERROR: queue out of order at', link%patch%id, link%patch%time, &
      previous
    previous = link%patch%time
    link => link%next_time
    nq = nq+1
  end do
  if (nq /= self%nq) print*,'ERROR: inconsistent number of patches in queue', nq, self%nq
END SUBROUTINE check_queue

!> =============================================================================
!> Print a list of nbors
!> =============================================================================
SUBROUTINE print_nbors (self, label)
  class (link_t):: self
  character(len=*):: label
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  if (label /= '') print'(a)', label
  print'(2(a,i5,2x),a,$)', 'patch', self%patch%id, 'level', self%patch%level, 'nbors:'
  link => self%nbors
  do while (associated(link))
    print '(i4,$)', link%patch%id
    link => link%next
  end do
  print*,''
END SUBROUTINE print_nbors

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_lists (self)
  class (patch_list_t):: self
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  link => self%head
  print'(a,i4,":",$)', 'task_list', self%n
  do while (associated(link))
    print '(i4,$)', link%patch%id
    link => link%next
  end do
  print*,''
  link => self%queue
  print'(a,i4,":",$)', 'ready_queue', self%nq
  do while (associated(link))
    print '(i4,$)', link%patch%id
    link => link%next_time
  end do
    print*,''
END SUBROUTINE print_lists

!> =============================================================================
!> Print a table with patch id and all neighbors
!> =============================================================================
SUBROUTINE print (self, label)
  class (patch_list_t):: self
  character(len=*):: label
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  print*, label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
END SUBROUTINE print

!> =============================================================================
END MODULE patch_list_mod
