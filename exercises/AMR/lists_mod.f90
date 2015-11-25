!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Handles lists, with links that have pointers to the next link, to a patch,
!> and to a list with neighbors to the patch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE lists_mod
  USE io_mod
  USE patches_mod
  implicit none
  private
  !.............................................................................
  type, public:: link_t
    type(link_t), pointer:: next => null()
    type(link_t), pointer:: nbor => null()
    type(patch_t), pointer:: patch => null()
  end type
  !.............................................................................
  type, public:: list_t
    type(link_t), pointer:: head=>null()
    real:: size(3) = 1.0
    integer:: dim(3)=[4,4,1]
    integer:: n
  contains
    procedure:: init
    procedure:: update
    procedure:: print
    procedure:: init_nbors_by_location
    procedure:: init_nbors_by_overlap
  end type
  type(patch_t), pointer, dimension(:,:,:):: patches
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  class (list_t):: self
  type(patch_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link, nbor
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  !.............................................................................
  self%n = 0
  allocate (patches(self%dim(1),self%dim(2),self%dim(3)))
  do i=1,self%dim(1)
  do j=1,self%dim(2)
  do k=1,self%dim(3)
    !---------------------------------------------------------------------------
    ! Initialize each patch
    !---------------------------------------------------------------------------
    patch => patches(i,j,k)
    call patch%init(self%size,self%dim,i,j,k)
    self%n = self%n+1
    !---------------------------------------------------------------------------
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
  call self%init_nbors_by_location
  call self%init_nbors_by_overlap
END SUBROUTINE
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
  link => self%head
  do i=1,self%dim(1)
  do j=1,self%dim(2)
  do k=1,self%dim(3)
    do i1=-1,1
    do j1=-1,1
    do k1=max(-1,self%dim(3)/2), min(1,self%dim(3)/2)
      i2 = modulo(i+i1-1,self%dim(1))+1
      j2 = modulo(j+j1-1,self%dim(2))+1
      k2 = modulo(k+k1-1,self%dim(3))+1
      if (i1==0.and.j1==0.and.k1==0) cycle
      allocate (nbor)
      nbor%patch => patches(i2,j2,k2)
      if (associated(link%nbor)) then
        prev%nbor => nbor
      else
        link%nbor => nbor
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
END SUBROUTINE init_nbors_by_location
!> =============================================================================
!> Initialize nbors via a double loop over the linked list.  This is simpler,
!> and completely general, in that no knowledge of the patch locations is needed.
!> It requires checking more patch combinations, but this is so fast, and is
!> done only at startup, that it is no real disadvantage
!> =============================================================================
SUBROUTINE init_nbors_by_overlap (self)
  class (list_t):: self
  type(link_t) , pointer:: prev, link1, link2, nbor
  !> ===========================================================================
  !> For each patch link, make a list of neighbors, starting with link%nbor
  !> ===========================================================================
  link1 => self%head
  do while (associated(link1))
    nullify (link1%nbor)
    link2 => self%head
    do while (associated(link2))
      if (.not.associated(link1,link2)) then
!        print'(2(i5,6f9.4,3x),l3)', &
!                link1%patch%id, link1%patch%mesh%p, link1%patch%mesh%distance(link2%patch%mesh), &
!                link2%patch%id, link2%patch%mesh%p, link2%patch%mesh%distance(link1%patch%mesh), &
!                link1%patch%mesh%overlaps(link2%patch%mesh)
        if (link1%patch%mesh%overlaps(link2%patch%mesh)) then
          allocate (nbor)
          nbor%patch => link2%patch
          if (associated(link1%nbor)) then
            prev%nbor => nbor
          else
            link1%nbor => nbor
          end if
          prev => nbor
        end if
      end if
      link2 => link2%next
    end do
    link1 => link1%next
  end do
  call self%print (' patches in list, with nbors by overlap:')
END SUBROUTINE init_nbors_by_overlap
!> =============================================================================
!> Skeleton routine, for now
!> =============================================================================
SUBROUTINE update (self)
  class (list_t):: self
  !.............................................................................
  self%n = self%n-1
END SUBROUTINE update
!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print (self, label)
  class (list_t):: self
  character(len=*):: label
  type(link_t) , pointer:: link, nbor
  !.............................................................................
  link => self%head
  print*, label
  do while (associated(link))
    nbor => link%nbor
    print '(i7,1h:,3x,3f7.3,"   ",$)', link%patch%id, link%patch%mesh%p
    do while (associated(nbor))
      print '(i4,$)', nbor%patch%id
      nbor => nbor%nbor
    end do
    print*,''
    link => link%next
  end do
END SUBROUTINE print
!> =============================================================================
END MODULE lists_mod
