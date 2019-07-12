!===============================================================================
!> Add gravity to any solver
!===============================================================================
MODULE gravity_mod
  USE patch_mod
  USE units_mod
  USE scaling_mod
  USE link_mod
  USE io_mod
  USE kinds_mod
  USE trace_mod
  implicit none
  private
  type, public:: gravity_t
    procedure(const_grav), pointer :: grav
    real(kind=KindScalarVar)       :: const_grav = 0.0
    real(kind=KindScalarVar)       :: mass = -1.0
    real(8)                        :: minrad = 0.0
    real(8)                        :: position(3) = 0.0                    ! centre of gravity
    integer                        :: axis = 3
    class(patch_t), pointer:: patch
    logical:: on=.false.
  contains
    procedure:: init
    procedure:: pre_update
  end type
  type(gravity_t), public:: gravity
CONTAINS

!===============================================================================
!> Initialize gravity
!===============================================================================
SUBROUTINE init (self, link)
  class(gravity_t):: self
  class(link_t), pointer :: link
  ! ............................................................................
  class(patch_t), pointer:: patch
  integer                           :: iostat, gn(3)
  real(kind=KindScalarVar), save    :: constant = 0.0
  logical, save                     :: first_time=.true.
  logical, save :: on          = .false.
  real(kind=KindScalarVar), save:: mass = -1.0
  real(8),                  save:: minrad = 0.0
  real(8), save :: position(3) = 0.0                    ! centre of gravity
  integer, save :: axis  = 3
  namelist /gravity_params/ on, constant, axis, position, mass, minrad
  !-----------------------------------------------------------------------------
  call trace%begin ('selfgravity_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, gravity_params, iostat=iostat)
    if (io%master) write (io%output, gravity_params)
    if (constant/=0.0) constant = constant * scaling%t**2 / scaling%l
  end if
  !$omp end critical (input_cr)
  if (on) then
    self%patch => task2patch(link%task)
    self%on = on
    self%const_grav = constant
    self%axis = axis
    self%position   = position
    self%mass       = mass
    self%minrad     = minrad
    if (constant /= 0.0) then
      self%grav => const_grav
    else if (mass > -1.0) then
      self%grav => var_grav
    else
      call io%abort("gravity_mod%init:: neither constant nor mass is set.")
    end if
    if (.not.allocated(self%patch%force_per_unit_mass)) then
      gn = self%patch%gn
      allocate (self%patch%force_per_unit_mass(gn(1),gn(2),gn(3),3))
      self%patch%force_per_unit_mass = 0.0
    end if
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE pre_update(self)
  class(gravity_t) :: self
  ! ----------------------------------------------------------------------------
  call self%grav()
END SUBROUTINE pre_update
  
!===============================================================================
!===============================================================================
SUBROUTINE const_grav(self)
  class(gravity_t) :: self
  ! ----------------------------------------------------------------------------
  self%patch%force_per_unit_mass(:,:,:,self%axis) = &
  self%patch%force_per_unit_mass(:,:,:,self%axis) + self%const_grav
END SUBROUTINE const_grav

!===============================================================================
!===============================================================================
SUBROUTINE var_grav(self)
  class(gravity_t) :: self
  ! ----------------------------------------------------------------------------
  integer:: ix, iy, iz, i
  real(8):: rp
  real(8):: x, y, z, xs, ys, zs, xp, yp, zp, ds2(3)
  real(kind=KindScalarVar):: newton
  ! ............................................................................
  newton = self%mass*scaling%grav
  associate (m=>self%patch%mesh)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do i=1,3
    ds2(i) = m(i)%h(self%patch%idx%px+i-1)  *m(i)%d
  end do
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz)
    zp = z  - self%position(3)
    zs = zp + ds2(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy)
      yp = y  - self%position(2)
      ys = yp + ds2(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix)
        xp = x  - self%position(1)
        xs = xp + ds2(1)
        !-----------------------------------------------------------------------
        ! The force of gravity is included explicitly, and is oriented
        ! away from the center.
        !-----------------------------------------------------------------------
        rp = sqrt(xs**2+yp**2+zp**2)
        rp = max(rp,self%minrad)
        self%patch%force_per_unit_mass(ix,iy,iz,1) = &
        self%patch%force_per_unit_mass(ix,iy,iz,1) - newton*xs/rp**3
        rp = sqrt(xp**2+ys**2+zp**2)
        rp = max(rp,self%minrad)
        self%patch%force_per_unit_mass(ix,iy,iz,2) = &
        self%patch%force_per_unit_mass(ix,iy,iz,2) - newton*ys/rp**3
        rp = sqrt(xp**2+yp**2+zs**2)
        rp = max(rp,self%minrad)
        self%patch%force_per_unit_mass(ix,iy,iz,3) = &
        self%patch%force_per_unit_mass(ix,iy,iz,3) - newton*zs/rp**3
      end do
    end do
  end do
  end associate
  end associate
END SUBROUTINE var_grav

END MODULE gravity_mod
