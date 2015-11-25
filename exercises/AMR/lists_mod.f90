!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Handles lists, with links that have pointers to the next link, to a patch,
!> and to a list with neighbors to the patch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE lists_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  implicit none
  private
  !.............................................................................
  type, public:: link_t
    type(link_t), pointer:: next => null()
    type(link_t), pointer:: nbors => null()
    type(patch_t), pointer:: patch => null()
    integer:: n_nbors = 0
  contains
    procedure:: search_nbors_by_overlap
  end type
  !.............................................................................
  type, public:: list_t
    type(link_t), pointer:: head=>null(), tail=>null()
    real:: size(3) = 1.0
    integer:: dims(3)=[4,4,1]
    character(len=120):: name = 'patch list'
    integer:: n = 0
  contains
    procedure:: init
    procedure:: init_cartesian
    procedure:: update
    procedure:: split
    procedure:: print
    procedure:: init_nbors_by_location
    procedure:: init_nbors_by_overlap
    procedure:: init_nbor_by_overlap
  end type
  type(patch_t), pointer, dimension(:,:,:):: patches
  integer:: msplit=2
CONTAINS
!> =============================================================================
SUBROUTINE init (self, name)
  class (list_t):: self
  character(len=*):: name
  !.............................................................................
  self%name = name
END SUBROUTINE init
!> =============================================================================
SUBROUTINE init_cartesian (self)
  class (list_t):: self
  integer:: dims(3) = [4,4,1]
  real:: size(3) = [1.0,1.0,1.0]
  type(patch_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link, nbor
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  namelist /box_params/ size, dims, msplit
  !.............................................................................
  rewind (io%input); read (io%input, box_params)
  if (io%master) write (*, box_params)
  self%dims = dims
  self%size = size
  allocate (patches(self%dims(1),self%dims(2),self%dims(3)))
  do i=1,self%dims(1)
  do j=1,self%dims(2)
  do k=1,self%dims(3)
    !---------------------------------------------------------------------------
    ! Initialize each patch
    !---------------------------------------------------------------------------
    patch => patches(i,j,k)
    call patch%init (self%size/self%dims)
    patch%i = i; patch%j = j; patch%k = k
    patch%mesh%p = -0.5_8*self%size + ([i,j,k]-0.5_8)*patch%mesh%s
!    print'(i5,3(3f9.4,3x))', &
!            patch%id, patch%mesh%p, patch%mesh%s, patch%mesh%d
    ! Add to a list of links, with pointers from patch to patch
    !---------------------------------------------------------------------------
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
  self%tail => link
  call self%init_nbors_by_location
  call self%init_nbors_by_overlap
END SUBROUTINE init_cartesian
!> =============================================================================
!> Initialize nbors via a double loop, using explicit knowledge about placements
!> =============================================================================
SUBROUTINE init_nbors_by_location (self)
  class (list_t):: self
  type(patch_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link, nbor
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('list_t%init_nbors_by_overlap', 1)
  link => self%head
  do i=1,self%dims(1)
  do j=1,self%dims(2)
  do k=1,self%dims(3)
    do i1=-1,1
    do j1=-1,1
    do k1=max(-1,self%dims(3)/2), min(1,self%dims(3)/2)
      i2 = modulo(i+i1-1,self%dims(1))+1
      j2 = modulo(j+j1-1,self%dims(2))+1
      k2 = modulo(k+k1-1,self%dims(3))+1
      if (i1==0.and.j1==0.and.k1==0) cycle
      allocate (nbor)
      link%n_nbors = link%n_nbors+1
      nbor%patch => patches(i2,j2,k2)
      if (associated(link%nbors)) then
        prev%next => nbor
      else
        link%nbors => nbor
      end if
      prev => nbor
    end do
    end do
    end do
    link => link%next
  end do
  end do
  end do
  call self%print (' patches in list, with nbors by location:')
  call trace_end
END SUBROUTINE init_nbors_by_location
!> =============================================================================
!> Initialize nbors via a double loop over the linked list.  This is simpler,
!> and completely general, in that no knowledge of the patch locations is needed.
!> It requires checking more patch combinations, but this is so fast, and is
!> done only at startup, that it is no real disadvantage
!> =============================================================================
SUBROUTINE init_nbors_by_overlap (self)
  class (list_t):: self
  type(link_t) , pointer:: prev, link
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbors
  !> ===========================================================================
  call trace_begin('list_t%init_nbors_by_overlap', 1)
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
  class (list_t):: self
  type(link_t) , pointer:: prev, link1, link2, nbor
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('list_t%init_nbor_by_overlap', 2)
  link1%n_nbors = 0
  nullify (link1%nbors)
  link2 => self%head
  if (io%verbose>3) print*,' id, size, position, distance:'
  do while (associated(link2))
    if (.not.associated(link1,link2)) then
      if (io%verbose>3) &
        print'(2(i5,3(3f8.4,2x)),l3)', &
          link1%patch%id, link1%patch%mesh%s, link1%patch%mesh%p, link1%patch%mesh%distance(link2%patch%mesh), &
          link2%patch%id, link2%patch%mesh%s, link2%patch%mesh%p, link2%patch%mesh%distance(link1%patch%mesh), &
          link1%patch%mesh%overlaps(link2%patch%mesh)
      if (link1%patch%mesh%overlaps(link2%patch%mesh)) then
        allocate (nbor)
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
  class (link_t), target:: self
  type(link_t) , pointer:: link1, link2, start, prev, nbor
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  call trace_begin('list_t%search_nbor_by_overlap', 3)
  link1 => self                                 ! start on the self list
  if (io%verbose>3) print*,' id, size, position, distance:'
  do while (associated(link1))
    prev => link1%nbors                         ! existing nbor list head (empty?)
    link2 => start                              ! start searching throuch nbor chain
    do while (associated(link2))
      if (.not.associated(link1%patch,link2%patch)) then
        if (io%verbose>3) &
          print'(2(i5,3(3f8.4,2x)),l3)', &
            link1%patch%id, link1%patch%mesh%s, link1%patch%mesh%p, link1%patch%mesh%distance(link2%patch%mesh), &
            link2%patch%id, link2%patch%mesh%s, link2%patch%mesh%p, link2%patch%mesh%distance(link1%patch%mesh), &
            link1%patch%mesh%overlaps(link2%patch%mesh)
        if (link1%patch%mesh%overlaps(link2%patch%mesh)) then
          allocate (nbor)                       ! new nbor (has no nbor)
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
!> Skeleton routine, for now
!> =============================================================================
SUBROUTINE update (self)
  class (list_t):: self
  !.............................................................................
  call trace_begin('list_t%update', 1)
  self%n = self%n-1
  call trace_end
END SUBROUTINE update
!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list, and initializing their nbor lists.  Also make a temporary list
!> with the new links, for use in nbor searches.  Note that the link to the
!> original patch is kept, with a "leaf" flag being cleared, indicating it is
!> not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split (self, this)
  class(list_t):: self
  integer:: nsplit(3)
  type(patch_t), pointer:: patch
  type(link_t), pointer:: this, link, nbor, tail, new_patches
  integer dir, i, ix, iy, iz
  real(8):: pos(3)
  !.............................................................................
  call trace_begin('list_t%split', 2)
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
        call patch%init (this%patch%size/nsplit)
        patch%mesh%p = this%patch%mesh%p + pos*patch%size
        !------------------------------------------------------------------------
        ! Append links to the new patches as the list tail
        !------------------------------------------------------------------------
        allocate (link)
        link%patch => patch
        self%tail%next => link
        self%tail => link
        !------------------------------------------------------------------------
        ! Prepend new nbor links to the this%nbors chain
        !------------------------------------------------------------------------
        allocate (nbor)
        nbor%patch => patch
        nbor%next => this%nbors
        this%nbors => nbor
        this%n_nbors = this%n_nbors+1
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! If all new patches have already been added to the this%nbors chain, all we
  ! need to do is to update the nbors list of the new patches by searching in
  ! the this%nbors list
  !-----------------------------------------------------------------------------
  new_patches => new_patches%next
  call new_patches%search_nbors_by_overlap (this%nbors)
  !
  this%patch%is_leaf = .false.
  call trace_end
END SUBROUTINE split
!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print (self, label)
  class (list_t):: self
  character(len=*):: label
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  link => self%head
  if (io%master) then
    print*, label
    print*,' id         position   n_nbor:  neighbors'
  end if
  do while (associated(link))
    nbor => link%nbors
    print '(i4,1h:,3f7.3, i4,":",$)', link%patch%id, link%patch%mesh%p, link%n_nbors
    do while (associated(nbor))
      print '(i4,$)', nbor%patch%id
      nbor => nbor%next
    end do
    print*,''
    link => link%next
  end do
END SUBROUTINE print
!> =============================================================================
END MODULE lists_mod
