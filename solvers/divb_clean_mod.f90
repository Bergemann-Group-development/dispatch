!> =============================================================================
!> $Id$
!> =============================================================================
MODULE divb_clean_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE vector_mod
  USE scalar_mod
  USE vector_ops
  implicit none
  public
  integer, private:: verbose=0
  type, public:: divb_clean_t
    procedure(divB_clean0),   pointer    :: clean => null()
  contains
    procedure :: init => init_clean
  end type
  integer, save:: divb_cleaner=0
  integer, save:: divb_iter=4
  logical, save:: on = .false.
  type(divb_clean_t), public:: divb
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init_clean(self)
  class(divb_clean_t) :: self
  integer:: iostat
  logical, save:: first_time = .true.
  namelist /divb_clean_params/ on, divb_cleaner, divb_iter
  ! ----------------------------------------------------------------------------
  call io%lock%set ('divb_clean')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, divb_clean_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, divb_clean_params)
    end if
    call io%namelist_warning ('divb_clean_params', iostat)
    write (io%namelist, divb_clean_params,iostat=iostat)
    select case(divb_cleaner)
    case(0)
      self%clean => divb_clean0
    case(1)
      self%clean => divb_clean1
    case(2)
      self%clean => divb_clean2
    case(3)
      self%clean => divb_clean3
    case(4)
      self%clean => divb_clean3
    case default
      call io%abort('divb_clean_mod:: divb cleaner method not understood')
    end select
  end if
  call io%lock%unset ('divb_clean')
END SUBROUTINE init_clean

!===============================================================================
!> Determine the guard zone values of the magnetic field components from the
!> div(B)=0 condition.  As appreciated from the expressions below, the only
!> external knowledge that is required are the perpendicular components of the
!> magnetic field along cube edges, and values interpolated from nbors are used
!> there.
!===============================================================================
SUBROUTINE divb_clean0 (self, m, f)
  class(divb_clean_t) :: self
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer:: f
  !.............................................................................
  real, dimension(:,:,:), pointer:: b1, b2, b3
  real(8):: ds(3)
  integer:: iter, i1, i2, i3, li(3), ui(3), lb(3), ub(3), n(3)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin('mhd_t%divb_clean0', itimer=itimer)
  li = m%li
  ui = m%ui
  lb = m%lb
  ub = m%ub
  ds = m%d
  n  = m%n
  b1 => f(:,:,:,1)
  b2 => f(:,:,:,2)
  b3 => f(:,:,:,3)
  do iter=1,divb_iter
    if (n(1)>1) call iter_x
    if (n(2)>1) call iter_y
    if (n(3)>1) call iter_z
    !call iter_y
  end do
  call trace%end (itimer)
contains
  subroutine iter_x
    !---------------------------------------------------------------------------
    ! Left x-face.  Initially un-constrained: b2(z) and b3(y) along edges
    !---------------------------------------------------------------------------
    do i1=li(1),lb(1),-1
    do i3=lb(3),ub(3)-1
    do i2=lb(2),ub(2)-1
      b1(i1,i2,i3) =                 b1(i1+1,i2,i3)             &
                   + (1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3))  &
                   +  1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3)))  * ds(1)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Right x-face.  Initially unconstrained: b2(z) and b3(y) along edges
    !---------------------------------------------------------------------------
    do i1=ui(1),ub(1)-1
    do i3=lb(3),ub(3)-1
    do i2=lb(2),ub(2)-1
      b1(i1+1,i2,i3) =                           + b1(i1,i2,i3)  &
                     - (1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3)) &
                     +  1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3))) * ds(1)
    end do
    end do
    end do
  end subroutine iter_x
  subroutine iter_y
    !---------------------------------------------------------------------------
    ! Left y-face.  .  Initially unconstrained: b3(y) along edges
    !---------------------------------------------------------------------------
    do i2=li(2),lb(2),-1
    do i3=lb(3),ub(3)-1
    do i1=lb(1),ub(1)-1
      b2(i1,i2,i3) =            b2(i1,i2+1,i3)                 &
                   + (1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3)) &
                   +  1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3))) * ds(2)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Right y-face.
    !---------------------------------------------------------------------------
    do i2=ui(2),ub(2)-1
    do i3=lb(3),ub(3)-1
    do i1=lb(1),ub(1)-1
      b2(i1,i2+1,i3) =          + b2(i1,i2,i3)                   &
                     - (1./ds(3)*(b3(i1,i2,i3+1) - b3(i1,i2,i3)) &
                     +  1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3))) * ds(2)
    end do
    end do
    end do
  end subroutine iter_y
  subroutine iter_z
    !---------------------------------------------------------------------------
    ! Left z-face.  Not known: b3(:,:,l3)
    !---------------------------------------------------------------------------
    do i3=li(3),lb(3),-1
    do i2=lb(2),ub(2)-1
    do i1=lb(1),ub(1)-1
      b3(i1,i2,i3) =            b3(i1,i2,i3+1)                 &
                   + (1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3)) &
                   +  1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3))) * ds(3)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Right z-face: Not known, b3(i1,i2,i3+1)
    !---------------------------------------------------------------------------
    do i3=ui(3),ub(3)-1
    do i2=lb(2),ub(2)-1
    do i1=lb(1),ub(1)-1
      b3(i1,i2,i3+1) =          + b3(i1,i2,i3)                   &
                     - (1./ds(1)*(b1(i1+1,i2,i3) - b1(i1,i2,i3)) &
                     +  1./ds(2)*(b2(i1,i2+1,i3) - b2(i1,i2,i3))) * ds(3)
    end do
    end do
    end do
  end subroutine iter_z
END SUBROUTINE divb_clean0

!===============================================================================
!> Perform a single iteration of a div(B) clean.  The limit of stability is
!> eps=1./6.  This is intended to prevent build-up of layers of large div(B)
!> between patches. An exact solution would be to solve an internal Laplace
!> equation, with boundary conditiosn = the boundary value of div(B), and then
!> subtract the gradient of the resulting potential.
!>
!> The problem can also be greatly reduced by using a higher order time
!> interpoation scheme, thus making the electric field computed from the
!> guard zone values of B and U a better match to the neighbors E-fields.
!===============================================================================
SUBROUTINE divb_clean1(self, m, f)
  class(divb_clean_t) :: self
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer:: f
  !.............................................................................
  real, dimension(:,:,:)  , pointer:: bx, by, bz, divb
  integer:: gn(3)
  real, parameter:: eps=0.11
  !-----------------------------------------------------------------------------
  if (.not.on) return
  gn = m%gn
  bx => f(:,:,:,1)
  by => f(:,:,:,2)
  bz => f(:,:,:,3)
  call allocate_scalars (gn, divb)
  divb = ddup(bx,1) + ddup(by,2) + ddup(bz,3)
  bx = bx + eps*dddn(divb,1)
  by = by + eps*dddn(divb,2)
  bz = bz + eps*dddn(divb,3)
  call deallocate_scalars(divb)
    return
contains
  !-----------------------------------------------------------------------------
  function ddup(ff,i) result (g)
  real:: ff(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = cshift(ff,1,i) - ff
  end function
  !-----------------------------------------------------------------------------
  function dddn(ff,i) result (g)
  real:: ff(:,:,:), g(size(f,1),size(f,2),size(f,3))
  integer:: i
  !.............................................................................
  g = ff - cshift(ff,-1,i)
  end function
END SUBROUTINE divb_clean1

!===============================================================================
!> Enforce div(B)=0 in the guard zones, by succesively imposing the condition,
!> layer-by-layer, on the component perpendicular to each of the patch faces.
!===============================================================================
SUBROUTINE divb_clean2 (self, m, f)
  class(divb_clean_t) :: self
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer:: f
  !.............................................................................
  real, dimension(:,:,:), pointer:: bx, by, bz
  real, dimension(:,:,:), allocatable:: divb
  real, parameter:: eps=0.11
  integer:: l(3), u(3), ix, iy, iz, i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin('mhd_t%divb_clean2', itimer=itimer)
  bx => f(:,:,:,1)
  by => f(:,:,:,2)
  bz => f(:,:,:,3)
  l = m%lb
  u = m%ub-1
  i1 = merge(1, 0, m(1)%n > 1)
  j1 = merge(1, 0, m(2)%n > 1)
  k1 = merge(1, 0, m(3)%n > 1)
  !--- x-direction -------------------------------------------------------------
  do ix=m(1)%lo,m(1)%lb,-1
    bx(ix   ,l(2):u(2),l(3):u(3)) = &
    bx(ix+i1,l(2):u(2),l(3):u(3)) + ( &
      by(ix,l(2)+j1:u(2)+j1,l(3)   :u(3)   ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)   :u(2)   ,l(3)+k1:u(3)+k1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  do ix=m(1)%ui,m(1)%ub-1
    bx(ix+i1,l(2):u(2),l(3):u(3)) = &
    bx(ix   ,l(2):u(2),l(3):u(3)) - ( &
      by(ix,l(2)+j1:u(2)+j1,l(3)   :u(3)   ) - by(ix,l(2):u(2),l(3):u(3)) + &
      bz(ix,l(2)   :u(2)   ,l(3)+k1:u(3)+k1) - bz(ix,l(2):u(2),l(3):u(3)))
  end do
  !--- y-direction -------------------------------------------------------------
  do iy=m(2)%lo,m(2)%lb,-1
    by(l(1):u(1),iy   ,l(3):u(3)) = &
    by(l(1):u(1),iy+j1,l(3):u(3)) + ( &
      bx(l(1)+i1:u(1)+i1,iy,l(3)   :u(3)   ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)   :u(1)   ,iy,l(3)+k1:u(3)+k1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  do iy=m(2)%ui,m(2)%ub-1
    by(l(1):u(1),iy+j1,l(3):u(3)) = &
    by(l(1):u(1),iy   ,l(3):u(3)) - ( &
      bx(l(1)+i1:u(1)+i1,iy,l(3)   :u(3)   ) - bx(l(1):u(1),iy,l(3):u(3)) + &
      bz(l(1)   :u(1)   ,iy,l(3)+k1:u(3)+k1) - bz(l(1):u(1),iy,l(3):u(3)))
  end do
  !--- z-direction -------------------------------------------------------------
  do iz=m(3)%lo,m(3)%lb,-1
    bz(l(1):u(1),l(2):u(2),iz  ) = &
    bz(l(1):u(1),l(2):u(2),iz+k1) + ( &
      bx(l(1)+i1:u(1)+i1,l(2)   :u(2)   ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)   :u(1)   ,l(2)+j1:u(2)+j1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  do iz=m(3)%ui,m(3)%ub-1
    bz(l(1):u(1),l(2):u(2),iz+k1) = &
    bz(l(1):u(1),l(2):u(2),iz   ) - ( &
      bx(l(1)+i1:u(1)+i1,l(2)   :u(2)   ,iz) - bx(l(1):u(1),l(2):u(2),iz) + &
      by(l(1)   :u(1)   ,l(2)+j1:u(2)+j1,iz) - by(l(1):u(1),l(2):u(2),iz))
  end do
  call trace%end (itimer)
END SUBROUTINE divb_clean2

!===============================================================================
!> Enforce div(B)=0 in the guard zones, by succesively imposing the condition,
!> layer-by-layer, on the component perpendicular to each of the patch faces,
!> and on edges with 1, 2, or 3 guard zone surfaces.
!===============================================================================
SUBROUTINE divb_clean3 (self, m, f)
  class(divb_clean_t) :: self
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer:: f
  !.............................................................................
  real, dimension(:,:,:), pointer:: bx, by, bz
  real:: divb, divb_max
  integer:: l(3), u(3), ix, iy, iz
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin('mhd_t%divb_clean3', itimer=itimer)
  bx => f(:,:,:,1)
  by => f(:,:,:,2)
  bz => f(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Start with a cube define below.  It has div(B)=0 for all indices, except
  ! where one or more indicaes hit the upper boundary u(:). We acan extend that by
  ! fixing all of the upper boundary cases
  !-----------------------------------------------------------------------------
  l = m%li
  u = m%ui
  !-----------------------------------------------------------------------------
  ! Upper yz-surface
  !-----------------------------------------------------------------------------
  ix = u(1)
  do iz=l(3),u(3)
  do iy=l(2),u(2)
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bx(ix+1,iy,iz) = &
    bx(ix+1,iy,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,1
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Upper xz-surface
  !-----------------------------------------------------------------------------
  iy = u(2)
  do iz=l(3),u(3)
  do ix=l(1),u(1)
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    by(ix,iy+1,iz) = &
    by(ix,iy+1,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,2
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Upper xy-surface
  !-----------------------------------------------------------------------------
  iz = u(3)
  do iy=l(2),u(2)
  do ix=l(1),u(1)
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bz(ix,iy,iz+1) = &
    bz(ix,iy,iz+1) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,3
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Edge cells, Ex
  !-----------------------------------------------------------------------------
  iy = u(2)
  iz = u(3)
  do ix=m(1)%lo,m(1)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    by(ix,iy+1,iz) = &
    by(ix,iy+1,iz) - divb/2.
    bz(ix,iy,iz+1) = &
    bz(ix,iy,iz+1) - divb/2.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,4
  end do
  iy = l(2)
  iz = u(3)
  do ix=m(1)%lo,m(1)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bz(ix,iy,iz+1) = &
    bz(ix,iy,iz+1) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,5
  end do
  iy = u(2)
  iz = l(3)
  do ix=m(1)%lb,m(1)%ub-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    by(ix,iy+1,iz) = &
    by(ix,iy+1,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,6
  end do
  !-----------------------------------------------------------------------------
  ! Edge cells, By
  !-----------------------------------------------------------------------------
  ix = u(1)
  iz = u(3)
  do iy=m(2)%lo,m(2)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bx(ix+1,iy,iz) = &
    bx(ix+1,iy,iz) - divb/2.
    bz(ix,iy,iz+1) = &
    bz(ix,iy,iz+1) - divb/2.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,7
  end do
  ix = l(1)
  iz = u(3)
  do iy=m(2)%lo,m(2)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bx(ix+1,iy,iz) = &
    bx(ix+1,iy,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,8
  end do
  ix = u(1)
  iz = l(3)
  do iy=m(2)%lo,m(2)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bz(ix,iy,iz+1) = &
    bz(ix,iy,iz+1) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,9
  end do
  !-----------------------------------------------------------------------------
  ! Edge cells, Bz
  !-----------------------------------------------------------------------------
  ix = u(1)
  iy = u(2)
  do iz=m(3)%lo,m(3)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bx(ix+1,iy,iz) = &
    bx(ix+1,iy,iz) - divb/2.
    by(ix,iy+1,iz) = &
    by(ix,iy+1,iz) - divb/2.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,10
  end do
  ix = l(1)
  iy = u(2)
  do iz=m(3)%lo,m(3)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    by(ix,iy+1,iz) = &
    by(ix,iy+1,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,11
  end do
  ix = u(1)
  iy = l(2)
  do iz=m(3)%lo,m(3)%lb,-1
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    bx(ix+1,iy,iz) = &
    bx(ix+1,iy,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,12
  end do
  !-----------------------------------------------------------------------------
  ! Corner cell, all 3 upper
  !-----------------------------------------------------------------------------
  ix = u(1)
  iy = u(2)
  iz = u(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  bx(ix+1,iy,iz) = &
  bx(ix+1,iy,iz) - divb/3.
  by(ix,iy+1,iz) = &
  by(ix,iy+1,iz) - divb/3.
  bz(ix,iy,iz+1) = &
  bz(ix,iy,iz+1) - divb/3.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,13
  !-----------------------------------------------------------------------------
  ! Corner cells, 2 upper
  !-----------------------------------------------------------------------------
  ix = l(1)
  iy = u(2)
  iz = u(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  by(ix,iy+1,iz) = &
  by(ix,iy+1,iz) - divb/2.
  bz(ix,iy,iz+1) = &
  bz(ix,iy,iz+1) - divb/2.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,14
  ix = u(1)
  iy = l(2)
  iz = u(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  bx(ix+1,iy,iz) = &
  bx(ix+1,iy,iz) - divb/2.
  bz(ix,iy,iz+1) = &
  bz(ix,iy,iz+1) - divb/2.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,15
  ix = u(1)
  iy = u(2)
  iz = l(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  bx(ix+1,iy,iz) = &
  bx(ix+1,iy,iz) - divb/2.
  by(ix,iy+1,iz) = &
  by(ix,iy+1,iz) - divb/2.
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,16
  !-----------------------------------------------------------------------------
  ! Corner cells, 1 upper
  !-----------------------------------------------------------------------------
  ix = u(1)
  iy = l(2)
  iz = l(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  bx(ix+1,iy,iz) = &
  bx(ix+1,iy,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,17
  ix = l(1)
  iy = u(2)
  iz = l(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  by(ix,iy+1,iz) = &
  by(ix,iy+1,iz) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,18
  ix = l(1)
  iy = l(2)
  iz = u(3)
  divb = &
    (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
    (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
    (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
  bz(ix,iy,iz+1) = &
  bz(ix,iy,iz+1) - divb
    divb = &
      (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
      (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
      (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
    if (abs(divb) > 1e-4) print*,'ERROR',divb,19
  !-----------------------------------------------------------------------------
  ! Validate
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    divb_max = 0.0
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      divb_max = max (abs(divb), divb_max)
    end do
    end do
    end do
    print *, 'divb_max', divb_max
  end if
  !-----------------------------------------------------------------------------
  ! div(B) is now zero in all cells from l=li to u=uo, and we can extend this
  ! outwards surface by surface -- this will leave some corner value, but they
  ! are never used.  First Bx
  !-----------------------------------------------------------------------------
  if (divb_cleaner==4) then
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=u(1)+1,m(1)%ub-1
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      bx(ix+1,iy,iz) = &
      bx(ix+1,iy,iz) - divb
    end do
    end do
    end do
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1)-1,m(1)%lb,-1
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      bx(ix,iy,iz) = &
      bx(ix,iy,iz) + divb
    end do
    end do
    end do
    !-----------------------------------------------------------------------------
    ! Validate
    !-----------------------------------------------------------------------------
    if (verbose > 0) then
      divb_max = 0.0
      do iz=l(3),u(3)
      do iy=l(2),u(2)
      do ix=m(1)%lb,m(1)%ub-1
        divb = &
          (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
          (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
          (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
        divb_max = max (abs(divb), divb_max)
      end do
      end do
      end do
      print *, 'divb_max(x)', divb_max
    end if
    !-----------------------------------------------------------------------------
    ! By
    !-----------------------------------------------------------------------------
    do iz=l(3),u(3)
    do iy=u(2)+1,m(2)%ub-1
    do ix=l(1),u(1)
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      by(ix,iy+1,iz) = &
      by(ix,iy+1,iz) - divb
    end do
    end do
    end do
    do iz=l(3),u(3)
    do iy=l(2)-1,m(2)%lb,-1
    do ix=l(1),u(1)
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      by(ix,iy,iz) = &
      by(ix,iy,iz) + divb
    end do
    end do
    end do
    !-----------------------------------------------------------------------------
    ! Validate
    !-----------------------------------------------------------------------------
    if (verbose > 0) then
      divb_max = 0.0
      do iz=l(3),u(3)
      do iy=m(2)%lb,m(2)%ub-1
      do ix=l(1),u(1)
        divb = &
          (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
          (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
          (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
        divb_max = max (abs(divb), divb_max)
      end do
      end do
      end do
      print *, 'divb_max(y)', divb_max
    end if
    !-----------------------------------------------------------------------------
    ! Bz
    !-----------------------------------------------------------------------------
    do iz=u(3)+1,m(2)%ub-1
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      bz(ix,iy,iz+1) = &
      bz(ix,iy,iz+1) - divb
    end do
    end do
    end do
    do iz=l(3)-1,m(3)%lb,-1
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      divb = &
        (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
        (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
        (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
      bz(ix,iy,iz) = &
      bz(ix,iy,iz) + divb
    end do
    end do
    end do
    !-----------------------------------------------------------------------------
    ! Validate
    !-----------------------------------------------------------------------------
    if (verbose > 0) then
      divb_max = 0.0
      do iz=m(3)%lb,m(3)%ub-1
      do iy=l(2),u(2)
      do ix=l(1),u(1)
        divb = &
          (bx(ix+1,iy  ,iz  ) - bx(ix,iy,iz)) + &
          (by(ix  ,iy+1,iz  ) - by(ix,iy,iz)) + &
          (bz(ix  ,iy  ,iz+1) - bz(ix,iy,iz))
        divb_max = max (abs(divb), divb_max)
      end do
      end do
      end do
      print *, 'divb_max(z)', divb_max
    end if
  end if
   call trace%end (itimer)
END SUBROUTINE divb_clean3

END MODULE divb_clean_mod

