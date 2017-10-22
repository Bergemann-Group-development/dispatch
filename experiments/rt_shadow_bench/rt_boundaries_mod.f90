!> =============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!> =============================================================================
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
  USE rt_components_mod
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
    procedure:: condition => condition_mem                              ! modify outer BCs if necessary in no chains mode
    procedure:: condition_cell                                          ! modify outer BCs in chains mode
  end type
!-----------------------------------------------------------------------

  type :: rt_inner_bndry_t
    real(8) :: position(3)
    real(8) :: velocity(3)
  end type
!-----------------------------------------------------------------------
  type, extends(rt_inner_bndry_t) :: rt_inner_sphere_t
    real(8) :: radius
  end type
!-----------------------------------------------------------------------
  type, extends(rt_inner_bndry_t) :: rt_inner_cuboid_t
    real(8) :: llc(3), urc(3)
  end type
!-----------------------------------------------------------------------
  type, public:: rt_boundaries_t
    class(rt_inner_bndry_t), dimension(:), pointer :: inner_bndry => null()
    real, dimension(:), pointer    :: xl, xu, yl, yu, zl, zu
    real(8) :: llc(3), urc(3)
    real(8) :: position(3), radius, rk, src
    integer :: hot_dir
  contains
    procedure:: init
    procedure:: update
    procedure:: spherical_jagged
  end type
  type(rt_boundaries_t), public :: rt_boundaries
CONTAINS

!> =====================================================================
!> Initialize the boundary data type
!> =====================================================================
SUBROUTINE init (self, llc, urc)
  class(rt_boundaries_t)              :: self
  real(8):: llc(3), urc(3)
  !.....................................................................
  real, save :: hot_temp = 1000, cool_temp = 10
  real :: d0
  logical, save:: first_time=.true.
  real, dimension(:,:,:), pointer :: lnd, e
  real, dimension(:,:,:,:), pointer :: src, rk
  integer :: idir
  class(ray_dirs_t), pointer :: rd
  namelist /rt_bc_params/ hot_temp, cool_temp
  !---------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_bc_params)
    if (io%master) write (*, rt_bc_params)
    call scaling%init()
    allocate(lnd(1,1,2))
    allocate(e(1,1,2))
    allocate(src(1,1,2, eos%n_lambda))
    allocate(rk(1,1,2, eos%n_lambda))
    allocate(self%xl(eos%n_lambda))
    allocate(self%xu(eos%n_lambda))
    allocate(self%yl(eos%n_lambda))
    allocate(self%yu(eos%n_lambda))
    allocate(self%zl(eos%n_lambda))
    allocate(self%zu(eos%n_lambda))
    d0 = 1.0/scaling%d
    lnd = log(d0)
    e(:,:,1) = (cool_temp/scaling%temp)/(eos%gamma-1.0)
    e(:,:,2) = (hot_temp/scaling%temp)/(eos%gamma-1.0)
    call eos%lookup (shape(e), x=e, lny=lnd, &                           ! call the table lookup, ...
                    src=src, rk=rk)
    !src = src * rk  
    self%xl = src(1,1,2,:) - src(1,1,1,:)
    self%xu = 0.0
    self%yl = 0.0
    self%yu = 0.0
    self%zl = 0.0
    self%zu = 0.0
    deallocate(lnd,e,src, rk)
    self%llc = llc
    self%urc = urc
    ! find the hot direction ID
    rd => ray_geometry%head%ray_list
    rd_loop: do while (associated (rd))
      if (all(rd%bundle(1)%dstep == [1,0,0])) then
        self%hot_dir = rd%id
        exit rd_loop
      end if
    rd => rd%next
    end do rd_loop
  end if
  !$omp end critical (input_cr)
  print *, 'cool', self%xu
  print *, 'hot', self%xl
  print *, 'rt boundary llc, urc', self%llc, self%urc
END SUBROUTINE init

SUBROUTINE update(self)
  class(rt_boundaries_t) :: self
  !---------------------------------------------------------------------
END SUBROUTINE update

!> =====================================================================
!> Initialize the outer boundary shells for rt patches
!> =====================================================================
SUBROUTINE init_outer_bndry (self, n, ndir, mesh)
  class(rt_outer_bndry_t) :: self
  class(mesh_t), pointer, optional :: mesh(:)
  integer :: n(3), ndir
  !.....................................................................
  integer :: i, nl, saved_status
  !---------------------------------------------------------------------
  if (self%unset) then
    if (present(mesh)) self%m => mesh
    self%ndir = ndir
    nl = eos%n_lambda
    saved_status=self%bits
    if (self%is_set(bits%xl)) then
      allocate(self%xl(nl*2))
      do i=1, eos%n_lambda
        self%xl(i) = rt_boundaries%xl(i)
        self%xl(i+nl) = rt_boundaries%xu(i)
      end do
    end if
    if (self%is_set(bits%xu)) then
      allocate(self%xu(nl))
      do i=1, eos%n_lambda
        self%xu(i) = rt_boundaries%xu(i)
      end do
    end if
    if (self%is_set(bits%yl)) then
      allocate(self%yl(nl))
      do i=1, eos%n_lambda
        self%yl(i) = rt_boundaries%yl(i)
      end do
    end if
    if (self%is_set(bits%yu)) then
      allocate(self%yu(nl))
      do i=1, eos%n_lambda
        self%yu(i) = rt_boundaries%yu(i)
      end do
    end if
    if (self%is_set(bits%zl)) then
      allocate(self%zl(nl))
      do i=1, eos%n_lambda
        self%zl(i) = rt_boundaries%zl(i)
      end do
    end if
    if (self%is_set(bits%zu)) then
      allocate(self%zu(nl))
      do i=1, eos%n_lambda
        self%zu(i) = rt_boundaries%zu(i)
      end do
    end if
    self%bits=saved_status
    self%unset = .false.
  end if
END SUBROUTINE init_outer_bndry

!> =============================================================================
!> Set bits in the status word
!> =============================================================================
SUBROUTINE set (self, bits)
  class(rt_outer_bndry_t) :: self
  integer                :: bits
  !-----------------------------------------------------------------------------
  self%bits = ior (self%bits, bits)
END SUBROUTINE set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
LOGICAL FUNCTION is_set (self, bits)
  class(rt_outer_bndry_t) :: self
  integer               :: bits
  !-----------------------------------------------------------------------------
  is_set = (iand (self%bits, bits) /=0)
END FUNCTION is_set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bits)
  class(rt_outer_bndry_t) :: self
  integer               :: bits
  !-----------------------------------------------------------------------------
  is_clear = (iand (self%bits, bits) ==0)
END FUNCTION is_clear

!> =====================================================================
!> Apply a boundary condition
!> =====================================================================
SUBROUTINE condition_mem (self, f)
  class(rt_outer_bndry_t) :: self
  real(kind=KindScalarVar)           :: f(:,:,:,:)
  !.....................................................................
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i, saved_status, il, iu, a, b
  integer, save:: itimer
  !---------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%condition', itimer=itimer)
  if (self%unset) then
    call mpi%abort ('rt_outer_boundaries_t%init has not been called')
  end if
  saved_status = self%bits
  if (self%is_set(bits%xl)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    il = 1
    iu = il + eos%n_lambda-1
   ! do i=1,self%ndir
      do b=m3%lo,m3%uo
        do a=m2%lo,m2%uo
          f(m1%lo,a,b,:) = self%xl(2)
          f(m1%lo,a,b,rt_boundaries%hot_dir) = self%xl(1)
        end do
      end do
      il = il + eos%n_lambda
      iu = il + eos%n_lambda-1
  !  end do
  end if
  if (self%is_set(bits%xu)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    il = 1
    iu = il + eos%n_lambda-1
   ! do i=1,self%ndir
      do b=m3%lo,m3%uo
        do a=m2%lo,m2%uo
          f(m1%uo,a,b,:) = self%xu(1)
        end do
      end do
      il = il + eos%n_lambda
      iu = il + eos%n_lambda-1
    !end do
  end if
  if (self%is_set(bits%yl)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    il = 1
    iu = il + eos%n_lambda-1
  !  do i=1,self%ndir
      do b=m3%lo,m3%uo
        do a=m1%lo,m1%uo
          f(a,m2%lo,b,:) = self%yl(1)
        end do
      end do
      il = il + eos%n_lambda
      iu = il + eos%n_lambda-1
  !  end do
  end if
  if (self%is_set(bits%yu)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    il = 1
    iu = il + eos%n_lambda-1
  !  do i=1,self%ndir
      do b=m3%lo,m3%uo
        do a=m1%lo,m1%uo
          f(a,m2%uo,b,:) = self%yu(1)
        end do
      end do
      il = il + eos%n_lambda
      iu = il + eos%n_lambda-1
  !  end do
  end if
  if (self%is_set(bits%zl)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    il = 1
    iu = il + eos%n_lambda-1
 !   do i=1,self%ndir
      do b=m2%lo,m2%uo
        do a=m1%lo,m1%uo
          f(a,b,m3%lo,:) = self%zl(1)
        end do
      end do
      il = il + eos%n_lambda
      iu = il + eos%n_lambda-1
  !  end do
  end if
  if (self%is_set(bits%zu)) then
    m1 => self%m(1)
    m2 => self%m(2)
    m3 => self%m(3)
    il = 1
    iu = il + eos%n_lambda-1
   ! do i=1,self%ndir
      do b=m2%lo,m2%uo
        do a=m1%lo,m1%uo
          f(a,b,m3%uo,:) = self%zu(1)
        end do
      end do
      il = il + eos%n_lambda
      iu = il + eos%n_lambda-1
  !  end do
  end if
  self%bits = saved_status
  call trace%end(itimer)
END SUBROUTINE condition_mem

!> =====================================================================
!> Apply a boundary condition
!> =====================================================================
SUBROUTINE condition_cell (self, v,cell,dir)
  class(rt_outer_bndry_t) :: self
  real(kind=KindScalarVar):: v(:)
  integer:: cell(3), dir
  !.....................................................................
  integer:: i, saved_status, il, iu, a, b
  integer, save:: itimer
  !---------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%condition', itimer=itimer)
  if (self%unset) then
    call mpi%abort ('rt_outer_boundaries_t%init has not been called')
  end if
  saved_status = self%bits
  if (self%is_set(bits%xl)) then
    if (dir==rt_boundaries%hot_dir) then
      if (cell(1) < self%m(1)%li) v = self%xl(1)
    else
      if (cell(1) < self%m(1)%li) v = self%xl(2)
    end if
  end if
  if (self%is_set(bits%xu)) then
    if (cell(1) > self%m(1)%ui) v = self%xu(1)
  end if
  if (self%is_set(bits%yl)) then
    if (cell(2) < self%m(2)%li) v = self%yl(1)
  end if
  if (self%is_set(bits%yu)) then
    if (cell(2) > self%m(2)%ui) v = self%yu(1)
  end if
  if (self%is_set(bits%zl)) then
    if (cell(3) < self%m(3)%li) v = self%zl(1)
  end if
  if (self%is_set(bits%zu)) then
    if (cell(3) > self%m(3)%ui) v = self%zu(1)
  end if
  self%bits = saved_status
  call trace%end(itimer)
END SUBROUTINE condition_cell

!> =============================================================================
!> Impose spherical boundary condition, assuming that the sphere is jagged.
!> =============================================================================
SUBROUTINE spherical_jagged (self, src, rk, m)
  class(rt_boundaries_t)                   :: self
  real(kind=KindRTVar), dimension(:,:,:,:) :: src, rk
  class(mesh_t), pointer                   :: m(:), m1, m2, m3
  !............................................................................
  integer      :: ix, iy, iz
  real(8)      :: x1, x2, y1, y2, z1, z2
  real(8)      :: r2
  logical      :: inside(8)
  integer, save:: itimer
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_boundaries_t%spherical_jagged', itimer=itimer)
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  error stop "rt_boundaries_mod::spherical_jagged: this version is a dummy - self%position is not set"
  r2 = self%radius**2
  do iz = m3%lb, m3%ub
    do iy = m2%lb, m2%ub
      do ix = m1%lb, m1%ub
        x1 = m1%p - self%position(1) + m1%r(ix) + 0.5*m1%d
        x2 = m1%p - self%position(1) + m1%r(ix) - 0.5*m1%d
        y1 = m2%p - self%position(2) + m2%r(iy) + 0.5*m2%d
        y2 = m2%p - self%position(2) + m2%r(iy) - 0.5*m2%d
        z1 = m3%p - self%position(3) + m3%r(iz) + 0.5*m3%d
        z2 = m3%p - self%position(3) + m3%r(iz) - 0.5*m3%d
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
  call trace%end (itimer)
END SUBROUTINE spherical_jagged


END MODULE rt_boundaries_mod
