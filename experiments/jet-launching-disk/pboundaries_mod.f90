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
SUBROUTINE apply_bcs (self, patch, mhd, cit, previt)
  USE bits_mod
  class(pboundaries_t) :: self
  class(patch_t), target :: patch
  integer, intent(in):: cit
  integer, intent(in), optional:: previt
  logical, intent(in):: mhd

  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: q
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin ('pboundaries_t%apply_bcs')

  m1 => patch%mesh(1)

  if (m1%type == 'cylindrical_mesh' .or. m1%type == 'spherical_mesh') then
    if (patch%boundaries%is_set(bits%xl)) then
      ! lower x1-boundary (z or r); OUTFLOW
      ! for cylindrical, this assumes the lower z-boundary is symmetric to the upper
      ! z-boundary, ie. is not at the disk midplane.
      do iv=1,patch%nv
        q => patch%mem(:,:,:,iv,cit,1)
        call self%outflow(q,iv,select=1)
      end do
      if (mhd .and. present(previt)) call update_b(patch, cit, previt, select=1)
    end if
    if (patch%boundaries%is_set(bits%xu)) then
      ! upper x1-boundary (z or r); OUTFLOW
      do iv=1,patch%nv
        q => patch%mem(:,:,:,iv,cit,1)
        call self%outflow(q,iv,select=2)
      end do
      if (mhd .and. present(previt)) call update_b(patch, cit, previt, select=2)
    end if
    if (patch%boundaries%is_set(bits%yl)) then
      ! lower x2-boundary (r or theta); AXIS REFLECTION
      do iv=1,patch%nv
        q => patch%mem(:,:,:,iv,cit,1)
        call self%axis_reflection(q,iv,select=3)
      end do
      if (mhd .and. present(previt)) call update_b(patch, cit, previt, select=3)
    end if
    if (patch%boundaries%is_set(bits%yu)) then
      ! upper x2-boundary (r or theta); OUTFLOW (ZRP) or AXIS REFLECTION (RTP)
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
      if (mhd .and. present(previt)) call update_b(patch, cit, previt, select=4)
    end if
    if (patch%boundaries%is_set(bits%zl)) then
      ! should be periodic in phi
      print*,'WARNING! The problem should be periodic in phi! (zl)'
    end if
    if (patch%boundaries%is_set(bits%zu)) then
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
!> For the EMFs, skin values are intentionally NOT set (i.e. we keep the computed
!> values).
!> =============================================================================
SUBROUTINE outflow (self, mem, iv, select)
  class(pboundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem
  integer:: iv
  integer, optional:: select

  integer:: i, j, k, l_select
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin('pboundaries_t%outflow')

  ! Don't set boundary conditions directly on the magnetic field!
  ! Set BCs on the EMFs instead!
  if (iv >= idx%bx .and. iv <= idx%bz) return

  m1 => self%m(1)
  m2 => self%m(2)
  m3 => self%m(3)

  l_select = 0
  if (present(select)) l_select = select

  select case (l_select)
  case (1)
    do k=m3%lb,m3%ub
      if (iv == self%idx%px) then
        ! prevent unphysical inflow
        do j=m2%lb,m2%ub
          do i=m1%lb,m1%li
            mem(i,j,k) = min(mem(m1%li+1,j,k), 0.0)
          end do
        end do
      else if (iv == self%idx%ex .or. iv <= self%idx%ez) then
        do j=m2%lb,m2%gn
          do i=m1%lb,m1%lo
            mem(i,j,k) = mem(m1%li,j,k)
          end do
        end do
      else if (iv == self%idx%ey) then
        do j=m2%lb,m2%ub
          do i=m1%lb,m1%lo
            mem(i,j,k) = mem(m1%li,j,k)
          end do
        end do
      else
        do j=m2%lb,m2%ub
          do i=m1%lb,m1%lo
            mem(i,j,k) = mem(m1%li+(m1%lo-i),j,k)
          end do
        end do
      end if
    end do
  case (2)
    do k=m3%lb,m3%ub
      if (iv == self%idx%px) then
          ! prevent unphysical inflow
        do j=m2%lb,m2%ub
          do i=m1%uo,m1%gn
            mem(i,j,k) = max(mem(m1%ui,j,k), 0.0)
          end do
        end do
      else if (iv == self%idx%ex) then
        do j=m2%lb,m2%gn
          do i=m1%uo,m1%ub
            mem(i,j,k) = mem(m1%ui,j,k)
          end do
        end do
      else if (iv == self%idx%ey) then
        do j=m2%lb,m2%ub
          do i=m1%uo+1,m1%gn
            mem(i,j,k) = mem(m1%uo,j,k)
          end do
        end do
      else if (iv == self%idx%ez) then
        do j=m2%lb,m2%gn
          do i=m1%uo+1,m1%gn
            mem(i,j,k) = mem(m1%uo,j,k)
          end do
        end do
      else
        do j=m2%lb,m2%ub
          do i=m1%uo,m1%ub
            mem(i,j,k) = mem(m1%ui-(i-m1%uo),j,k)
          end do
        end do
      end if
    end do
  case(3)
    do k=m3%lb,m3%ub
      do i=m1%lb,m1%ub
        if (iv == self%idx%py) then
          ! prevent unphysical inflow
          do j=m2%lb,m2%lo
            mem(i,j,k) = min(mem(i,m2%li+1,k), 0.0)
          end do
        else if (iv >= self%idx%ex .and. iv <= self%idx%ez) then
          do j=m2%lb,m2%lo
            mem(i,j,k) = mem(i,m2%li,k)
          end do
        else
          do j=m2%lb,m2%lo
            mem(i,j,k) = mem(i,m2%li+(m2%lo-j),k)
          end do
        end if
      end do
    end do
  case(4)
    do k=m3%lb,m3%ub
      do i=m1%lb,m1%ub
        if (iv == self%idx%py) then
          ! prevent unphysical inflow
          do j=m2%uo,m2%ub
            mem(i,j,k) = max(mem(i,m2%ui,k), 0.0)
          end do
        else if (iv == self%idx%ey) then
          do j=m2%uo,m2%ub
            mem(i,j,k) = mem(i,m2%ui,k)
          end do
        else if (iv == self%idx%ex .or. iv == self%idx%ez) then
          do j=m2%uo+1,m2%gn
            mem(i,j,k) = mem(i,m2%uo,k)
          end do
        else
          do j=m2%uo,m2%ub
            mem(i,j,k) = mem(i,m2%ui-(j-m2%uo),k)
          end do
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
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem
  integer:: iv
  integer, optional:: select

  integer:: i, j, k, l_select
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin('pboundaries_t%axis_reflection')

  ! Don't set boundary conditions directly on the magnetic field!
  ! Set BCs on the EMFs instead!
  if (iv >= idx%bx .and. iv <= idx%bz) return

  m1 => self%m(1)
  m2 => self%m(2)
  m3 => self%m(3)

  l_select = 0
  if (present(select)) l_select = select

  select case (l_select)
  case (3)
    do k=m3%lb,m3%ub
      if (iv == self%idx%py) then
        ! vector components which are normal to the boundary surface should
        ! be anti-symmetric across the boundary and zero at the boundary itself
        do i=m1%lb,m1%ub
          do j=m2%lb,m2%lo
            mem(i,j,k) = -mem(i,m2%li+(m2%lo-j),k)
          end do
          mem(i,m2%li,k) = 0.0
        end do
      else if (iv == self%idx%ex) then
        do i=m1%lb,m1%ub
          do j=m2%lb,m2%lo
            mem(i,j,k) = mem(i,2*m2%li-j,k)
          end do
          mem(i,m2%li,k) = 0.0
        end do
      else if (iv == self%idx%ey) then
        do i=m1%lb,m1%gn
          do j=m2%lb,m2%lo
            mem(i,j,k) = -mem(i,m2%li+m2%lo-j,k)
          end do
        end do
      else if (iv == self%idx%ez) then
        do i=m1%lb,m1%gn
          do j=m2%lb,m2%lo
            mem(i,j,k) = -mem(i,2*m2%li-j,k)
          end do
          mem(i,m2%li,k) = 0.0
        end do
      else
        do i=m1%lb,m1%ub
          do j=m2%lb,m2%lo
            mem(i,j,k) = mem(i,m2%li+(m2%lo-j),k)
          end do
        end do
      end if
    end do
  case (4)
    do k=m3%lb,m3%ub
      if (iv == self%idx%py) then
        ! vector components which are normal to the boundary surface should
        ! be anti-symmetric across the boundary and zero at the boundary itself
        do i=m1%lb,m1%ub
          do j=m2%uo,m2%ub
            mem(i,j,k) = -mem(i,m2%ui-(j-m2%uo),k)
          end do
          mem(i,m2%uo,k) = 0.0
        end do
      else if (iv == self%idx%ex) then
        do i=m1%lb,m1%ub
          do j=m2%uo+1,m2%gn
            mem(i,j,k) = mem(i,2*m2%uo-j,k)
          end do
          mem(i,m2%uo,k) = 0.0
        end do
      else if (iv == self%idx%ey) then
        do i=m1%lb,m1%gn
          do j=m2%uo,m2%ub
            mem(i,j,k) = -mem(i,m2%ui+m2%uo-j,k)
          end do
        end do
      else if (iv == self%idx%ez) then
        do i=m1%lb,m1%gn
          do j=m2%uo+1,m2%gn
            mem(i,j,k) = -mem(i,2*m2%uo-j,k)
          end do
          mem(i,m2%uo,k) = 0.0
        end do
      else
        do i=m1%lb,m1%ub
          do j=m2%uo,m2%ub
            mem(i,j,k) = mem(i,m2%ui-(j-m2%uo),k)
          end do
        end do
      end if
    end do
  case default
    print*,'WARNING! No boundary condition applied!'
  end select

  nullify (m1, m2, m3)
  call trace%end
END SUBROUTINE axis_reflection

!> =============================================================================
SUBROUTINE update_b (patch, cit, pit, select)
  USE vector_ops, only: curl_emfs

  class(patch_t), target:: patch
  integer, intent(in):: cit, pit
  integer, optional:: select

  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: emfs, b_old, b_new
  real(kind=KindScalarVar), dimension(:,:,:,:), allocatable :: dbdt, b_copy
  integer:: l(3), u(3), u_emf(3)
  !-----------------------------------------------------------------------------
  call trace%begin('pboundaries_mod::update_b')
  if (.not. present(select)) return

  l = patch%mesh%lb
  u = patch%mesh%ub
  u_emf = patch%mesh%gn
  select case (select)
  case(1)
    u(1) = patch%mesh(1)%li
  case(2)
    l(1) = patch%mesh(1)%ui
  case(3)
    u(2) = patch%mesh(2)%li
  case(4)
    l(2) = patch%mesh(2)%uo
  case default
    print*,'WARNING! No magnetic fields will be updated!'
    return
  end select

  emfs => patch%mem(l(1):u_emf(1),l(2):u_emf(2),l(3):u_emf(3),patch%idx%ex:patch%idx%ez,cit,1)
  allocate (dbdt(l(1):u_emf(1),l(2):u_emf(2),l(3):u_emf(3),3))
  dbdt = 0.0
  dbdt = curl_emfs(patch%mesh, emfs)

  b_old => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,pit,1)
  b_new => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,cit,1)
  allocate (b_copy(size(b_new,1),size(b_new,2),size(b_new,3),3))
  b_copy = b_new

  ! reduce stencil by one zone (because of the differencing), accounting
  ! for collapsed dimensions
  !where (u > 1)
  !  u = u - 1
  !end where
  b_new(l(1):u(1),l(2):u(2),l(3):u(3),:) = b_old(l(1):u(1),l(2):u(2),l(3):u(3),:) &
                                         + patch%dtime * dbdt(l(1):u(1),l(2):u(2),l(3):u(3),:)

  deallocate(dbdt, b_copy)
  nullify(emfs, b_old, b_new)
  call trace%end
END SUBROUTINE update_b

END MODULE pboundaries_mod
