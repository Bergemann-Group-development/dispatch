!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE rt_boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mpi_mod
  USE bits_mod
  USE mesh_mod
  USE eos_mod
  USE units_mod
  USE scaling_mod
  USE ray_geometry_v3_mod
  USE planet_mod
  USE atmosphere_mod
  implicit none
  private
  !-----------------------------------------------------------------------
  type, public :: rt_outer_bndry_t
    integer:: id = 0
    integer:: ndir
    integer:: bits = 0
    logical:: unset = .true.
    class(mesh_t), pointer:: m(:)
    real, dimension(:), pointer :: xl => null()
    real, dimension(:), pointer :: xu => null()
    real, dimension(:), pointer :: yl => null()
    real, dimension(:), pointer :: yu => null()
    real, dimension(:), pointer :: zl => null()
    real, dimension(:), pointer :: zu => null()
  contains
    procedure:: init => init_outer_bndry
    procedure:: set
    procedure:: is_set
    procedure:: is_clear
    procedure:: condition                                               ! modify outer BCs if necessary
  end type
!-----------------------------------------------------------------------

  type :: rt_inner_bndry_t
    real(8) :: position(3)
    real(8) :: velocity(3)
  end type
!-----------------------------------------------------------------------
  type, extends(rt_inner_bndry_t) :: rt_inner_sphere_t
    real(8) :: radius
    real(kind=KindRTVar) :: src, rk
  end type
!-----------------------------------------------------------------------
  type, extends(rt_inner_bndry_t) :: rt_inner_cuboid_t
    real(8) :: llc(3), urc(3)
  end type
!-----------------------------------------------------------------------
  type, public:: rt_boundaries_t
    class(rt_inner_sphere_t), dimension(:), pointer :: inner_bndry => null()
    real(kind=KindRTVar), dimension(:), pointer :: src => null()
    real(kind=KindRTVar), dimension(:), pointer :: rk  => null()
    real(8) :: llc(3), urc(3)
    real:: diffuse_radius
  contains
    procedure:: init
    procedure:: spherical_jagged => spherical_jagged2
    procedure:: spherical_jagged0
    procedure:: update
  end type
  type(rt_boundaries_t), public :: rt_boundaries
CONTAINS

!=======================================================================
!> Initialize the boundary data type
!=======================================================================
SUBROUTINE init (self, llc, urc)
  class(rt_boundaries_t)              :: self
  real(8):: llc(3), urc(3)
  !.....................................................................
  real :: t
  logical, save:: first_time=.true.
  real, dimension(:,:,:), pointer :: d, e
  real, dimension(:,:,:,:), pointer :: src, rk
  integer :: idir
  !---------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    first_time = .false.
    allocate(d(1,1,1))
    allocate(e(1,1,1))
    allocate(src(1,1,1, eos%n_lambda))
    allocate(rk(1,1,1, eos%n_lambda))
    allocate(self%src(eos%n_lambda))
    allocate(self%rk(eos%n_lambda))
    t = atmosphere%surface%temperature
    d = atmosphere%surface%density
    e = t/(eos%gamma-1.0)
    self%llc = llc
    self%urc = urc
    call eos%lookup (shape(e), x=e, y=d, &                           ! call the table lookup, ...
                    src=src, rk=rk)
    self%rk = rk(1,1,1,:)
    self%src = src(1,1,1,:)
    print *, 'Rt boundaries: surf t, d, e:', t, d, e
    print *, 'Rt boundaries:  CGS t, d, e:', t*scaling%temp, d*scaling%d, e*scaling%e
    print *, 'Rt boundaries: self% rk, src', self%rk, self%src
    deallocate(src,d,e,rk)
  end if
  call io%lock%unset
END SUBROUTINE init

SUBROUTINE update(self)
  class(rt_boundaries_t) :: self
  !---------------------------------------------------------------------
END SUBROUTINE update

!===============================================================================
!> Impose spherical boundary condition, assuming that the sphere is jagged.
!===============================================================================
SUBROUTINE spherical_jagged (self, src, rk, m)
  class(rt_boundaries_t)                   :: self
  real(kind=KindRTVar), dimension(:,:,:,:) :: src, rk
  class(mesh_t), pointer                   :: m(:), m1, m2, m3
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x1, x2, y1, y2, z1, z2
  real(8)      :: r2
  logical      :: inside(8)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%spherical_jagged', itimer=itimer)
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  r2 = planet%radius**2
  do iz = m3%lb, m3%ub
    do iy = m2%lb, m2%ub
      do ix = m1%lb, m1%ub
        x1 = m1%p - planet%position(1) + m1%r(ix) + 0.5*m1%d
        x2 = m1%p - planet%position(1) + m1%r(ix) - 0.5*m1%d
        y1 = m2%p - planet%position(2) + m2%r(iy) + 0.5*m2%d
        y2 = m2%p - planet%position(2) + m2%r(iy) - 0.5*m2%d
        z1 = m3%p - planet%position(3) + m3%r(iz) + 0.5*m3%d
        z2 = m3%p - planet%position(3) + m3%r(iz) - 0.5*m3%d
        inside(1) = x1**2+y1**2+z1**2 < r2
        inside(2) = x2**2+y1**2+z1**2 < r2
        inside(3) = x1**2+y2**2+z1**2 < r2
        inside(4) = x2**2+y2**2+z1**2 < r2
        inside(5) = x1**2+y1**2+z2**2 < r2
        inside(6) = x2**2+y1**2+z2**2 < r2
        inside(7) = x1**2+y2**2+z2**2 < r2
        inside(8) = x2**2+y2**2+z2**2 < r2
        if (all(inside)) then
          src(ix,iy,iz,:) = self%src
          rk (ix,iy,iz,:) = self%rk
        end if
      end do
    end do
  end do
  nullify (m1, m2, m3)
  call trace%end(itimer)
END SUBROUTINE spherical_jagged

!===============================================================================
!> Make the embryo totally opaque, and make sure no light is picked up from any
!> of the optically thick cells.  This needs to be done well inside the dS/dt=0
!> boundary, to avoid the spurious cooling / heating at the discontinuity.
!===============================================================================
SUBROUTINE spherical_jagged2 (self, src,rk, m)
  class(rt_boundaries_t)                   :: self
  real(kind=KindRTVar), dimension(:,:,:,:) :: src, rk
  class(mesh_t), pointer                   :: m(:)
  !............................................................................
  real(8)                                  :: position(3), rad2, r(3), f, c
  integer                                  :: ix, iy, iz
  !-----------------------------------------------------------------------------
  call trace_begin('rt_boundaries_t%jagged2')
  rad2 = (self%diffuse_radius*planet%radius)**2
  c = 1./log(1e9)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    r(3) = m(3)%p - planet%position(3) + r3(iz)
    do iy=m(2)%lb,m(2)%ub
      r(2) = m(2)%p - planet%position(2) + r2(iy)
      do ix=m(1)%lb,m(1)%ub
        r(1) = m(1)%p - planet%position(1) + r1(ix)
        !----------------------------------------------------------------------
        ! The factor f is close to unity for r >> planet%radius, and smoothly
        ! increases the opacity and reduces the source function as r approaches
        ! planet%radius.  f is 0.61 at r=2*rp and 0.13 at r=rp. This avoids the
        ! clitches caused by a rapid change, and cuts the point source region
        ! out from the diffuse source.  f cannot become smaller than 1e-9 =
        ! exp(-log(1e9)).
        !----------------------------------------------------------------------
        f = exp(-rad2/max(c*rad2,sum(r**2)))
        rk (ix,iy,iz,:) = rk (ix,iy,iz,:)/max(f,1e-4)
        src(ix,iy,iz,:) = src(ix,iy,iz,:)*f
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE spherical_jagged2

!===============================================================================
!> Impose spherical boundary condition, assuming that the sphere is jagged.
!===============================================================================
SUBROUTINE spherical_jagged0 (self, dq, m, v)
  class(rt_boundaries_t)                   :: self
  real(kind=KindRTVar), dimension(:,:,:,:) :: dq
  real(kind=KindRTVar):: v
  class(mesh_t), pointer                   :: m(:), m1, m2, m3
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x1, x2, y1, y2, z1, z2
  real(8)      :: r2, r(3)
  logical      :: inside(8)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%spherical_jagged0', itimer=itimer)
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  r2 = planet%radius**2
  do iz = m3%lb, m3%ub
    do iy = m2%lb, m2%ub
      do ix = m1%lb, m1%ub
        x1 = m1%p - planet%position(1) + m1%r(ix) + 0.5*m1%d
        x2 = m1%p - planet%position(1) + m1%r(ix) - 0.5*m1%d
        y1 = m2%p - planet%position(2) + m2%r(iy) + 0.5*m2%d
        y2 = m2%p - planet%position(2) + m2%r(iy) - 0.5*m2%d
        z1 = m3%p - planet%position(3) + m3%r(iz) + 0.5*m3%d
        z2 = m3%p - planet%position(3) + m3%r(iz) - 0.5*m3%d
        inside(1) = x1**2+y1**2+z1**2 < r2
        inside(2) = x2**2+y1**2+z1**2 < r2
        inside(3) = x1**2+y2**2+z1**2 < r2
        inside(4) = x2**2+y2**2+z1**2 < r2
        inside(5) = x1**2+y1**2+z2**2 < r2
        inside(6) = x2**2+y1**2+z2**2 < r2
        inside(7) = x1**2+y2**2+z2**2 < r2
        inside(8) = x2**2+y2**2+z2**2 < r2
        if (all(inside)) dq(ix,iy,iz,:) = v
      end do
    end do
  end do
  nullify (m1, m2, m3)
  call trace%end(itimer)
END SUBROUTINE spherical_jagged0

!=======================================================================
!> Initialize the outer boundary shells for rt patches
!=======================================================================
SUBROUTINE init_outer_bndry (self, n, ndir, mesh)
  class(rt_outer_bndry_t) :: self
  class(mesh_t), pointer, optional :: mesh(:)
  integer :: n(3), ndir
  !.....................................................................
  integer :: i, nl, saved_status
  !---------------------------------------------------------------------
  if (self%unset) then
    if (present(mesh)) self%m => mesh
    self%unset = .false.
  end if
END SUBROUTINE init_outer_bndry

!===============================================================================
!> Set bits in the status word
!===============================================================================
SUBROUTINE set (self, bits)
  class(rt_outer_bndry_t) :: self
  integer                :: bits
  !-----------------------------------------------------------------------------
  self%bits = ior (self%bits, bits)
END SUBROUTINE set

!===============================================================================
!> Check bits in the status word
!===============================================================================
LOGICAL FUNCTION is_set (self, bits)
  class(rt_outer_bndry_t) :: self
  integer               :: bits
  !-----------------------------------------------------------------------------
  is_set = (iand (self%bits, bits) /=0)
END FUNCTION is_set

!===============================================================================
!> Check bits in the status word
!===============================================================================
LOGICAL FUNCTION is_clear (self, bits)
  class(rt_outer_bndry_t) :: self
  integer               :: bits
  !-----------------------------------------------------------------------------
  is_clear = (iand (self%bits, bits) ==0)
END FUNCTION is_clear

!=======================================================================
!> Apply a boundary condition
!=======================================================================
SUBROUTINE condition (self, f,src,rk)
  class(rt_outer_bndry_t) :: self
  real(kind=KindRTVar)           :: f(:,:,:,:), src(:,:,:,:), rk(:,:,:,:)
  !.....................................................................
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i, saved_status, il, iu, a, b, iv, nv, lam
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%condition', itimer=itimer)
  if (self%unset) then
    call mpi%abort ('rt_outer_boundaries_t%init has not been called')
  end if
  saved_status = self%bits
  nv = size(f,4)
  !> -------------------------------------------------------------------
  !> set the boundary Q values to 0.0
  !> -------------------------------------------------------------------
  if (self%is_set(bits%xl)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    do b=m3%lo,m3%uo
      do a=m2%lo,m2%uo
        f(m1%lo,a,b,:) = 0d0
      end do
    end do
  end if
  if (self%is_set(bits%xu)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    do b=m3%lo,m3%uo
      do a=m2%lo,m2%uo
        f(m1%uo,a,b,:) = 0d0
      end do
    end do
  end if
  if (self%is_set(bits%yl)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    do b=m3%lo,m3%uo
      do a=m1%lo,m1%uo
        f(a,m2%lo,b,:) = 0d0
      end do
    end do
  end if
  if (self%is_set(bits%yu)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    do b=m3%lo,m3%uo
      do a=m1%lo,m1%uo
        f(a,m2%uo,b,:) = 0d0
      end do
    end do
  end if
  if (self%is_set(bits%zl)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    do b=m2%lo,m2%uo
      do a=m1%lo,m1%uo
        do iv=1,nv,eos%n_lambda
          do lam=1,eos%n_lambda
            f(a,b,m3%lo,iv+lam-1) = -src(a,b,m3%lo,lam)
          end do
        end do
      end do
    end do
  end if
  if (self%is_set(bits%zu)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    do b=m2%lo,m2%uo
      do a=m1%lo,m1%uo
        do iv=1,nv,eos%n_lambda
          do lam=1,eos%n_lambda
            f(a,b,m3%uo,iv+lam-1) = -src(a,b,m3%uo,lam)
          end do
        end do
      end do
    end do
  end if
  self%bits = saved_status
  call trace%end(itimer)
END SUBROUTINE condition

END MODULE rt_boundaries_mod
