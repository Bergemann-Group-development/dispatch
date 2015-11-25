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
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  class (list_t):: self
  type(patch_t), pointer, dimension(:,:,:):: patches
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
      if (i2==0.and.j2==0.and.k2==0) cycle
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
  call self%print
END SUBROUTINE init
!> =============================================================================
!> Fake routine, for now
!> =============================================================================
SUBROUTINE update (self)
  class (list_t):: self
  !.............................................................................
  self%n = self%n-1
END SUBROUTINE update
!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print (self)
  class (list_t):: self
  type(link_t) , pointer:: link, nbor
  link => self%head
  print*,' patches in list, with positions and nbors:'
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
