!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Handles patch_list, with links that have pointers to the next link, to a patch,
!> and to a list with neighbors to the patch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_list_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  USE timestep_mod
  USE refine_mod
  implicit none
  private
  !.............................................................................
  type, public:: link_t
    type(link_t), pointer:: next => null()
    type(timestep_t), pointer:: patch => null()
  end type
  !.............................................................................
  type, public:: patch_list_t
    type(link_t), pointer:: head=>null(), tail=>null()
    real:: size(3)=1.0
    integer:: dims(3)=[4,4,1]
    character(len=120):: name = 'patch list'
    integer:: n=0, nq=0
  contains
    procedure:: init
    procedure:: init_cartesian
    procedure:: update
    procedure:: remove_patch
    procedure:: split
    procedure:: print
  end type
  type(timestep_t), pointer, dimension(:,:,:):: patches
  integer:: msplit=2
CONTAINS
!> =============================================================================
SUBROUTINE init (self, name)
  class (patch_list_t):: self
  character(len=*):: name
  !.............................................................................
  call trace_begin('patch_list_t%init')
  self%name = name
  call trace_end
END SUBROUTINE init
!> =============================================================================
SUBROUTINE init_cartesian (self)
  class (patch_list_t):: self
  integer:: dims(3) = [4,4,1]
  real:: size(3) = [1.0,1.0,1.0]
  real:: position(3)
  type(timestep_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  namelist /box_params/ size, dims, msplit
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
    ! Initialize each patch
    !---------------------------------------------------------------------------
    patch => patches(i,j,k)
    position = -0.5_8*self%size + ([i,j,k]-0.5_8)*patch%size
    call patch%init_values (position, self%size/self%dims)
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
  class(link_t), pointer:: link
  class(timestep_t), pointer:: patch
  integer:: n1
  !.............................................................................
  call trace_begin('patch_list_t%update')
  link => self%head
  call trace('MK1')
  do while (associated(link))
  call trace('MK2')
    patch => link%patch                                   ! link patch
  call trace('MK3')
    if (refine%needed (patch)) then                       ! if refinement needed
  call trace('MK4')
      call self%split (link)                              ! refined patch go back to queue
  call trace('MK5')
    else
  call trace('MK6')
      call patch%update                                   ! update the patch
  call trace('MK7')
      if (patch%time > patch%end_time) then               ! finished:
  call trace('MK8')
        call self%remove_patch (patch)                    ! remove from self
  call trace('MK9')
      end if
  call trace('MK10')
    end if
  call trace('MK11')
    link => link%next
  end do
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
  !
  this%patch%is_leaf = .false.
  call trace_end
END SUBROUTINE split
!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print (self, label)
  class (patch_list_t):: self
  character(len=*):: label
  type(link_t) , pointer:: link
  !.............................................................................
  link => self%head
  if (io%master) then
    print*, label
    print*,' id         position'
  end if
  do while (associated(link))
    print '(i4,1h:,3f7.3)', link%patch%id, link%patch%position
    link => link%next
  end do
END SUBROUTINE print
!> =============================================================================
END MODULE patch_list_mod
