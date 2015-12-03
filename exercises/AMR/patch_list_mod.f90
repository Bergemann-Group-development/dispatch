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
    type(link_t), pointer:: next => null()
    class(link_t), pointer:: nbors => null()
    class(link_t), pointer:: parent => null()
    type(link_t), pointer:: next_time => null()
    type(timestep_t), pointer:: patch => null()
    integer:: n_nbors
  contains
    procedure:: get_guardcells
  end type
  !.............................................................................
  type, public:: patch_list_t
    type(link_t), pointer:: head=>null(), tail=>null(), queue=>null()
    real:: size(3)=1.0
    integer:: dims(3)=[4,4,1]
    character(len=120):: name = 'patch list'
    integer:: n=0, nq=0, n_nbor=0, n_check=0, n_ready=0
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
    procedure:: queue_by_time
    procedure:: output
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
  self%name = name
  call init_cartesian (self)
  call refine%init
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
    patch%origin = origin
    call patch%init (position, self%size/self%dims)
    if (io%verbose>1) &
      print*,'init_cartesian:', patch%id, self%size, patch%size
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
  call self%print ('after init_cartesian')
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
  !.............................................................................
  call trace_begin('patch_list_t%update')
  if (self%nq < 1) then
    n1 = self%nq
    call self%check_all
    !print *,'WARNING: check_all (should not be needed), found', self%nq, n1
  end if
  head => self%queue                                    ! queue start
  if (associated(head)) then
    self%queue => self%queue%next_time                  ! update the link
    self%nq = self%nq-1
    call head%get_guardcells
    patch => head%patch                                 ! head patch
    if (refine%needed (patch)) then                     ! if refinement needed
      call self%split (head)                            ! refined patchs go back to queue
      if (patch%level > refine%levelmax) then
        call self%output
        print*,'stopping after refine'
        stop
      end if
    else
      call patch%update                                 ! update the patch
      call self%check_nbors (head)                      ! any nbor ready?
      if (patch%time > patch%end_time) then             ! finished:
        call patch%output                               ! final one
        call self%remove_patch (patch)                  ! remove from self
      else                                              ! not finished:
        call patch%clear(bits%ready)                    ! needs new check
      end if
    end if
  else
    print*,'ERROR: no self%queue', self%nq, self%n
    stop
  end if
  if (io%verbose>0) &
  print '(i8,i8,i3,f10.6,2x,3f10.6,i6,i6,2x,a)', patch%id, patch%istep, patch%it, patch%time, patch%position, &
    patch%status, self%nq
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
      else
        self%head => link%next
      end if
    end if
    prev => link
    link => link%next
  end do
  call trace_end
END SUBROUTINE remove_patch

!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split (self, this)
  class(patch_list_t):: self
  integer:: nsplit(3)
  type(timestep_t), pointer:: patch
  type(link_t), pointer:: this, link, tail, new_patches
  integer dir, i, ix, iy, iz
  real:: pos(3), position(3)
  !.............................................................................
  call trace_begin('patch_list_t%split', 1)
  new_patches => self%tail
  nsplit = msplit
  where (this%patch%size==0)
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
        ! "pos" is incremented with one unit of the new size for each new patch,
        ! and is reset back to -0.5 before continuing to the next direction
        !------------------------------------------------------------------------
        allocate (patch, source=this%patch)
        patch%mem_allocated = .false.
        nullify (patch%mem)
        position = this%patch%position + pos*patch%size
        call patch%init (position, self%size/self%dims)
        !------------------------------------------------------------------------
        ! Append links to the new patches as the list tail
        !------------------------------------------------------------------------
        allocate (link)
        link%patch => patch
        self%tail%next => link
        self%tail => link
      end do
    end do
  end do
  new_patches => new_patches%next
  this%patch%is_leaf = .false.
  call trace_end
END SUBROUTINE split

!> =============================================================================
!> Initialize nbors via a double loop over the linked list.
!> =============================================================================
SUBROUTINE init_nbors_by_overlap (self)
  class (patch_list_t):: self
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
  class (patch_list_t):: self
  type(link_t) , pointer:: prev, link1, link2, nbor
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('patch_list_t%init_nbor_by_overlap', 2)
  link1%n_nbors = 0
  nullify (link1%nbors)
  link2 => self%head
  if (io%verbose>3) print*,' id, size, position, distance:'
  do while (associated(link2))
    if (.not.associated(link1,link2)) then
!      if (io%verbose>3) &
!        print'(2(i5,3(3f8.4,2x)),l3)', &
!          link1%patch%id, link1%patch%size, link1%patch%position, link1%patch%distance(link2%patch), &
!          link2%patch%id, link2%patch%size, link2%patch%position, link2%patch%distance(link1%patch), &
!          link1%patch%overlaps(link2%patch)
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
        if (io%verbose==3)print*,'overlaps', link1%patch%id, link2%patch%id, ' T'
      else
        if (io%verbose==3)print*,'overlaps', link1%patch%id, link2%patch%id, ' F'
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
  call trace_begin('patch_list_t%search_nbor_by_overlap', 3)
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
  if (patch%is_set(bits%ready)) return
  call trace_begin('patch_list_t%check_ready '//self%name,2)
  nbor => link%nbors                                    ! start on nbor list
  inbor = 0
  do while (associated (nbor))                          ! keep going until end
    nbpatch => nbor%patch
    if (patch%ahead_of(nbpatch)) then
      self%n_check = self%n_check+1
      inbor = min(inbor+1,mnbor)
      nbors(inbor) = nbpatch%id
      if (self%debug==1) &
        print*, patch%id, ' check_ready failed on patch', nbpatch%id
      if (self%debug>2) &
        print '(i4,a,i3,2x,30i3)', patch%id, ' failed:', inbor, nbors(1:inbor)
      call trace_end
      return                                            ! give up -- do nothing
    else
      inbor = min(inbor+1,mnbor)
      nbors(inbor) = nbpatch%id
!      print*,'nbor patch', nbpatch%id, nbpatch%time
    end if
    nbor => nbor%next                                   ! next nbor
  end do
  call patch%set(bits%ready)                            ! mark as checked
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
  type(link_t), pointer:: this, next, prev
  !.............................................................................
  call trace_begin ('list_t%queue_by_time '//self%name)
  nullify (this%next_time)
  nullify (prev)
  next => self%queue
  do while (associated(next))
    this%next_time => next
    if (next%patch%time >= this%patch%time) then
      self%nq = self%nq+1
      if (associated(prev)) then
        prev%next_time => this
      else
        self%queue => this
      end if
      call trace_end
      return
    end if
    prev => next
    next => next%next_time
  end do
  self%nq = self%nq+1
  if (associated(prev)) then
    prev%next_time => this
    nullify (this%next_time)
  else
    self%queue => this
  end if
  call trace_end
END SUBROUTINE queue_by_time

!> =============================================================================
!> Run through neighbor list and pick up guard cells
!> =============================================================================
SUBROUTINE get_guardcells (self)
  implicit none
  class (link_t):: self
  class (link_t),  pointer:: nbor
  class (patch_t), pointer:: patch
  integer :: nguard
  !.............................................................................
  call trace_begin('link_t%get_guardcells')
  nbor => self%nbors
  patch => self%patch
  nguard=0
  do while (associated (nbor))
    if (nbor%patch%is_relevant_for(patch)) then
      call nbor%patch%get_guardcells_for (patch)
      nguard = nguard + 1
    endif
    nbor => nbor%next
  end do
  if (io%verbose > 5) &
    print *, 'get_guardcells: Number of neighbors with guard cells:', nguard
  call trace_end
END SUBROUTINE get_guardcells

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
SUBROUTINE print (self, label)
  class (patch_list_t):: self
  character(len=*):: label
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  link => self%head
  print*, label
  do while (associated(link))
    nbor => link%nbors
    print '(i7,1h:,3x,3f7.3,"   ",$)', link%patch%id, link%patch%position
    do while (associated(nbor))
      print '(i4,$)', nbor%patch%id
      nbor => nbor%next
    end do
    print*,''
    link => link%next
  end do
END SUBROUTINE print

!> =============================================================================
END MODULE patch_list_mod
