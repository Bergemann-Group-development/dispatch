!> =============================================================================
!> Physical boundary conditions module for a jet launching disk
!> =============================================================================
MODULE pboundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mpi_mod
  USE index_mod
  USE mesh_mod
  USE patch_mod
  implicit none
  private

  type, public:: pboundaries_t
    type(index_t):: idx
    class(mesh_t), pointer:: m(:)
  contains
    procedure:: init_bdries
    procedure:: outflow
    procedure:: axis_reflection
    procedure:: apply => apply_bcs
  end type
  type(pboundaries_t), public:: pboundaries

CONTAINS

!> =============================================================================
!> Initialise physical boundary condition data type
!> =============================================================================
SUBROUTINE init_bdries (self, mesh, position, indices)
  class(pboundaries_t)              :: self
  class(mesh_t), pointer, optional :: mesh(:)
  real(8), optional                :: position(3)
  type(index_t), optional          :: indices
  !-----------------------------------------------------------------------------
  call trace%begin ('pboundaries_t%init')

  self%m => mesh
  ! we need the variable indices to properly set the BCs
  self%idx = indices

  call trace%end
END SUBROUTINE init_bdries

!> =============================================================================
!> Apply physical boundary conditions to current patch if applicable
!> =============================================================================
SUBROUTINE apply_bcs (self, patch, cit)
  USE bits_mod
  class(pboundaries_t) :: self
  class(patch_t), target :: patch
  integer, intent(in):: cit

  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: q
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin ('pboundaries_t%apply_bcs')

  m1 => patch%mesh(1)

  if (m1%type == 'cylindrical_mesh' .or. m1%type == 'spherical_mesh') then
    if (patch%boundaries%is_set(bits%xl)) then
      ! lower x1-boundary (z or r)
      do iv=1,patch%nv
        q => patch%mem(:,:,:,iv,cit,1)
        call self%outflow(q,iv,select=1)
      end do
    else if (patch%boundaries%is_set(bits%xu)) then
      ! upper x1-boundary (z or r)
      do iv=1,patch%nv
        q => patch%mem(:,:,:,iv,cit,1)
        call self%outflow(q,iv,select=2)
      end do
    else if (patch%boundaries%is_set(bits%yl)) then
      ! lower x2-boundary (r or theta)
      do iv=1,patch%nv
        q => patch%mem(:,:,:,iv,cit,1)
        call self%axis_reflection(q,iv,select=3)
      end do
    else if (patch%boundaries%is_set(bits%yu)) then
      ! upper x2-boundary (r or theta)
      if (m1%type == 'cylindrical_mesh') then
        do iv=1,patch%nv
          q => patch%mem(:,:,:,iv,cit,1)
          call self%outflow(q,iv,select=4)
        end do
      else if (m1%type == 'spherical_mesh') then
        do iv=1,patch%nv
          q => patch%mem(:,:,:,iv,cit,1)
          call self%axis_reflection(q,iv,select=4)
        end do
      end if
    else if (patch%boundaries%is_set(bits%zl)) then
      ! should be periodic in phi
      print*,'WARNING! The problem should be periodic in phi! (zl)'
    else if (patch%boundaries%is_set(bits%zu)) then
      ! should be periodic in phi
      print*,'WARNING! The problem should be periodic in phi! (zu)'
    end if
  else
    call mpi%abort('Only cylindrical and spherical coords. are currently supported for this experiment.')
  end if

  call trace%end
END SUBROUTINE apply_bcs

!> =============================================================================
!> Standard outflow conditions
!> Apply zero gradient boundary conditions to all variables, with the exception
!> that the velocity component normal to boundary does not result in inflow.
!> =============================================================================
SUBROUTINE outflow (self, mem, iv, select)
  class(pboundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem(:,:,:)
  integer:: iv
  integer, optional:: select

  integer:: i, j, k, l_select
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin('pboundaries_t%outflow')

  m1 => self%m(1)
  m2 => self%m(2)
  m3 => self%m(3)

  l_select = 0
  if (present(select)) l_select = select

  select case (l_select)
  case (1)
    do k=m3%lb,m3%ub
      do j=m2%lb,m2%ub
        if (iv == idx%px) then
          ! prevent unphysical inflow
          do i=m1%lb,m1%lo
            mem(i,j,k) = min(mem(m1%li+(m1%lo-i),j,k), 0.0)
          end do
        else
          do i=m1%lb,m1%lo
            mem(i,j,k) = mem(m1%li+(m1%lo-i),j,k)
          end do
        end if
      end do
    end do
  case (2)
    do k=m3%lb,m3%ub
      do j=m2%lb,m2%ub
        if (iv == idx%px) then
          ! prevent unphysical inflow
          do i=m1%uo,m1%ub
            mem(i,j,k) = max(mem(m1%ui-(i-m1%uo),j,k), 0.0)
          end do
        else
          do i=m1%uo,m1%ub
            mem(i,j,k) = mem(m1%ui-(i-m1%uo),j,k)
          end do
        end if
      end do
    end do
  case(3)
    do k=m3%lb,m3%ub
      do j=m2%lb,m2%lo
        if (iv == idx%py) then
          ! prevent unphysical inflow
          mem(m1%lb:m1%ub,j,k) = min(mem(m1%lb:m1%ub,m2%li+(m2%lo-j),k), 0.0)
        else
          mem(m1%lb:m1%ub,j,k) = mem(m1%lb:m1%ub,m2%li+(m2%lo-j),k)
        end if
      end do
    end do
  case(4)
    do k=m3%lb,m3%ub
      do j=m2%uo,m2%ub
        if (iv == idx%py) then
          ! prevent unphysical inflow
          mem(m1%lb:m1%ub,j,k) = max(mem(m1%lb:m1%ub,m2%ui-(j-m2%uo),k), 0.0)
        else
          mem(m1%lb:m1%ub,j,k) = mem(m1%lb:m1%ub,m2%ui-(j-m2%uo),k)
        end if
      end do
    end do
  case default
    print*,'WARNING! No boundary condition applied!'
  end select

  nullify (m1, m2, m3)
  call trace%end
END SUBROUTINE outflow

!> =============================================================================
!> Reflection at a coordinate axis boundary conditions
!> Most variables are symmetric, but certain vectors are antisymmetric.
!> =============================================================================
SUBROUTINE axis_reflection (self, mem, iv, select)
  class(pboundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem(:,:,:)
  integer:: iv
  integer, optional:: select

  integer:: i, j, k, l_select
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin('pboundaries_t%axis_reflection')

  m1 => self%m(1)
  m2 => self%m(2)
  m3 => self%m(3)

  l_select = 0
  if (present(select)) l_select = select

  select case (l_select)
  case (3)
    do k=m3%lb,m3%ub
      if (iv == self%idx%py .or. iv == self%idx%by) then
        ! vector components which are normal to the boundary surface should
        ! be anti-symmetric across the boundary and zero at the boundary itself
        do j=m2%lb,m2%lo
          mem(m1%lb:m1%ub,j,k) = -mem(m1%lb:m1%ub,m2%li+(m2%lo-j),k)
        end do
        mem(m1%lb:m1%ub,m2%li,k) = 0.0
      else
        do j=m2%lb,m2%lo
          mem(m1%lb:m1%ub,j,k) = mem(m1%lb:m1%ub,m2%li+(m2%lo-j),k)
        end do
      end if
    end do
  case (4)
    do k=m3%lb,m3%ub
      if (iv == self%idx%py .or. iv == self%idx%by) then
        ! vector components which are normal to the boundary surface should
        ! be anti-symmetric across the boundary and zero at the boundary itself
        do j=m2%uo,m2%ub
          mem(m1%lb:m1%ub,j,k) = -mem(m1%lb:m1%ub,m2%ui-(j-m2%uo),k)
        end do
        mem(m1%lb:m1%ub,m2%li,k) = 0.0
      else
        do j=m2%uo,m2%ub
          mem(m1%lb:m1%ub,j,k) = mem(m1%lb:m1%ub,m2%ui-(j-m2%uo),k)
        end do
      end if
    end do
  case default
    print*,'WARNING! No boundary condition applied!'
  end select

  nullify (m1, m2, m3)
  call trace%end
END SUBROUTINE axis_reflection

END MODULE pboundaries_mod
