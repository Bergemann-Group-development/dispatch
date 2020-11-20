!===============================================================================
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE math_mod
  USE scaling_mod
  USE index_mod
  USE vector_ops
  implicit none
  private
  type, public :: initial_t
    real :: x0 = 0.5
    real :: y0 = 0.5
    real :: z0 = 0.5
    real :: r0 = 0.25
    real :: b0 = 0.1, temp=1.0, gamma=5./3., beta
    real :: db
    integer :: expe
    logical :: mhd
    real(8) :: time
    logical :: pic_needed
    logical :: first_pic
    !procedure(initial_conditions), pointer :: condition => null()
    procedure(circular), pointer :: condition => null()
    !----Alfven wave----------------------------------------------------------
    character(len=64):: solver
    real :: alpha  = -99.0
    real :: theta  = -99.0
    real :: delta  = 0.01
    real :: b_par  = 1.0
    real :: b_perp = 0.1
    real :: vel    = 0.0
  contains
    procedure :: init
    procedure, nopass :: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, solver, gam)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gam
  real, save :: r0 = 0.25, B0 = 0.1, db=0.0, temp =1e-4, gamma=1.6667, beta=10.
  integer, save :: expe=1
  logical, save :: first_time=.true.
  logical, save :: mhd=.true.
  logical, save :: pic_needed = .false.
  logical, save :: first_pic = .true.
  namelist /IC_params/ r0, B0, db, temp, gamma, beta, expe, first_pic
  !.............................................................................
  call trace%begin ('initial_t%init')

  !$omp critical (initial_t)
  if (first_time) then
    rewind (io%input); read (io%input,initial_params)
    if (io%master) write(*,initial_params)
    first_time = .false.
  end if
  !$omp end critical (initial_t)
  self%r0 = r0/scaling%d
  self%B0 = B0/scaling%b
  self%db = db
  self%temp = temp/scaling%temp
  self%gamma = gamma
  self%mhd = mhd
  self%time = 0d0
  self%beta = beta
  self%expe = expe
  self%pic_needed = pic_needed
  self%first_pic = first_pic
  self%solver = solver
  !self%condition => initial_conditions
  self%condition => circular

  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE initial_conditions (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t), pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  !class(_t), pointer:: patch
  !.............................................................................
  integer:: ix, iy, iz                                       ! loop indices
  integer:: l(3), u(3)
  real:: x, y, z, dx, dy, dz
  real:: mu0 = 1.0, pg, d0
  real:: Flx, Fly, profile, g1, addon, pertub        ! local variables
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, Bx, By, Bz, pp
  real, dimension(:,:,:), pointer:: ex, ey, ez, jx,jy, jz
  real, dimension(:,:,:,:), pointer:: bfield, ufield, jfield
  !.............................................................................
  call trace_begin ('initial_t%init_mhd_vars')
  l = m%lb
  u = m%ub
  ! Correction due to cell size
  dx = +0.5*m(1)%d
  dy = +0.5*m(2)%d
  dz = -0.5*m(3)%d
  ! Correction due to patch size
  self%x0 = 0.5*m(1)%s
  self%y0 = 0.5*m(2)%s
  self%z0 = 0.5*m(3)%s
  d  => ff(:,:,:,idx%d)
  e  => ff(:,:,:,idx%e)
  px => ff(:,:,:,idx%px)
  py => ff(:,:,:,idx%py)
  pz => ff(:,:,:,idx%pz)
  bx => ff(:,:,:,idx%bx)
  by => ff(:,:,:,idx%by)
  bz => ff(:,:,:,idx%bz)
  ex => ff(:,:,:,idx%ex)
  ey => ff(:,:,:,idx%ey)
  ez => ff(:,:,:,idx%ez)

  ex = 0.0
  ey = 0.0
  ez = 0.0

  allocate(bfield(l(1)-2:u(1)+2,l(2)-2:u(2)+2,l(3)-2:u(3)+2,3))
  allocate(ufield(l(1)-2:u(1)+2,l(2)-2:u(2)+2,l(3)-2:u(3)+2,3))

  g1 = self%gamma-1.0

  pg = 0.5*self%b0**2*(self%beta+1.)
  d0 = pg/self%temp
  print *,' gas pressure:', pg
  print *,' mass density:', d0
  print *,' temperature :', self%temp

  select case(self%expe)

  ! The static equilibrium fluxtube experiment
  case(1)

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix)
          profile = exp(-2*(x**2+y**2)/self%r0**2)     ! Gaussian
          e(ix,iy,iz)  = self%B0**2/(2*mu0)*(1.+self%beta-profile)/g1      ! internal energy
          d(ix,iy,iz)  = e(ix,iy,iz)*g1/self%temp                          ! mass density
          px(ix,iy,iz) = 0.0                                               ! momentum
          py(ix,iy,iz) = 0.0
          pz(ix,iy,iz) = 0.0
        end do
      end do
    end do

! compute the fields
    do iz=m(3)%lb-2,m(3)%ub+2
      do iy=m(2)%lb-2,m(2)%ub+2
        if ((iy>=m(2)%lb).and.(iy<=m(2)%ub)) then
          y = m(2)%p + m(2)%r(iy)
        else
          if ((iy<m(2)%lb)) then
            y= m(2)%p + m(2)%r(m(2)%lb) - (m(2)%lb-iy)*m(2)%d
          else
            y= m(2)%p + m(2)%r(m(2)%ub) + (iy - m(2)%ub)*m(2)%d
          end if
        end if
        do ix=m(1)%lb-2,m(1)%ub+2
          if ((ix>=m(1)%lb).and.(ix<=m(1)%ub)) then
            x = m(1)%p + m(1)%r(ix)
          else
            if ((ix<m(1)%lb)) then
              x= m(1)%p + m(1)%r(m(1)%lb) - (m(1)%lb-ix)*m(1)%d
            else
              x= m(1)%p + m(1)%r(m(1)%ub) + (ix - m(1)%ub)*m(1)%d
            end if
          end if
          bfield(ix,iy,iz,1) = 0.0                                               ! Magnetic field
          bfield(ix,iy,iz,2) = 0.0
          bfield(ix,iy,iz,3) = self%B0 * exp(-((x+dx)**2+(y+dy)**2)/self%r0**2)
        end do
      end do
    end do
  ufield = 0.0
  bx = bfield(3:u(1)+2,3:u(2)+2,3:u(3)+2,1)
  by = bfield(3:u(1)+2,3:u(2)+2,3:u(3)+2,2)
  bz = bfield(3:u(1)+2,3:u(2)+2,3:u(3)+2,3)

  ! The fluxtube experiment with Gaussian sphere
  case(2)

    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(iz)
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix)
          profile = exp(-2*((x-self%x0)**2+(y-self%y0)**2)/self%r0**2)  &   ! Gaussian
               - 50 * exp(-0.5*(((x-0.9)**2+(y-0.5)**2+(z-0.5)**2)/(0.5*(self%r0-0.5)**2)))
          e(ix,iy,iz)  = self%B0**2/(2*mu0)*(1.+self%beta-profile)/g1      ! internal energy
          d(ix,iy,iz)  = e(ix,iy,iz)*g1/self%temp                          ! mass density
          px(ix,iy,iz) = 0.0                                               ! momentum
          py(ix,iy,iz) = 0.0
          pz(ix,iy,iz) = 0.0
          Bx(ix,iy,iz) = 0.0                                               ! Magnetic field
          By(ix,iy,iz) = 0.0
          Bz(ix,iy,iz) = self%B0 * exp(-((x+dx-self%x0)**2+(y+dy-self%y0)**2)/self%r0**2) 
        end do
      end do
    end do

  ! No fluxtube, to check boundaries
  case(3)


    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(iz)
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix)
          pertub = 50*exp(-0.5*(((x-0.9)**2+(y-0.5)**2+(z-0.5)**2)/ (0.5*(self%r0-0.5)**2))) 
          e(ix,iy,iz)  = self%B0**2/(2*mu0)*(1.+self%beta+pertub)/g1      ! internal energy
          d(ix,iy,iz)  = e(ix,iy,iz)*g1/self%temp                          ! mass density
          px(ix,iy,iz) = 0.0                                               ! momentum
          py(ix,iy,iz) = 0.0
          pz(ix,iy,iz) = 200.0 !-10000.0
          Bx(ix,iy,iz) = 0.0                                               ! Magnetic field
          By(ix,iy,iz) = 0.0
          Bz(ix,iy,iz) = 0.0
        end do
      end do
    end do

  end select
  

  if (io%debug(2)) then
    iz = (m(3)%lb+m(3)%ub)/2
    do ix=m(1)%li,m(1)%ui
      print *,'initial B(x):',ix,iz,bz(ix,iz,iz)
    end do
    do iy=m(1)%li,m(1)%ui
      print *,'initial B(y):',iz,iy,bz(iz,iy,iz)
    end do
    print *,'minval,maxval(d)',minval(d),maxval(d)
  end if

  deallocate(bfield, ufield)

  call trace_end
END SUBROUTINE initial_conditions

!> =============================================================================
!> Setup a patch with a circular Alfven wave 
!> can be used both with pure MHD and with MHD-PIC integration
!> =============================================================================
SUBROUTINE circular (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:), pointer :: ff, B, J
  ! ....................................................................
  integer :: ix, iy, iz
  real    ::  x,  y,  z, xs, ys, zs, xu, yu, zu
  real:: pgas, d0, theta, b0
  real:: k, d, bx, by, bz, ux, uy, uz, ex, ey, ez
  ! --------------------------------------------------------------------
  call trace%begin ('initial_t%circular')
  pgas  = 0.0002
  d0    = 1.0
  b0    = 0.1
  theta = atan(100./80.)
  k = math%pi2!sqrt((math%pi2/100.)**2+(math%pi2/80.)**2)
  !
  do iz=m(3)%lb,m(3)%ub
     z  = m(3)%p + m(3)%r(iz)
    zs  = z -0.5 * m(3)%d
    zu  = z +0.5 * m(3)%d
    do iy=m(2)%lb,m(2)%ub
      y   = m(2)%p + m(2)%r(iy)
      ys  = y -0.5 * m(2)%d
      yu  = y +0.5 * m(2)%d
      do ix=m(1)%lb,m(1)%ub
        x  = m(1)%p + m(1)%r(ix)
        xs = x -0.5 * m(1)%d
        xu = x +0.5 * m(1)%d
        !
        d  = d0 / scaling%d
        bx = b0 / scaling%b
        ux = 0.0
        if (trim(self%solver) /= 'stagger2e_pic') then
          by = b0*self%delta*sin(k*x) / scaling%b
          bz = b0*self%delta*cos(k*x) / scaling%b
          uy = self%delta*sin(k*x)    / scaling%u
          uz = self%delta*cos(k*x)    / scaling%u

          ex = uz*by - uy*bz
          ey = ux*b0*self%delta*cos(k*xs)/scaling%b - &
                     self%delta*cos(k*xs)/scaling%u * bx
          ez =       self%delta*sin(k*xs)/scaling%u * bx - &
               ux*b0*self%delta*sin(k*xs)/scaling%b
        else
          by = b0*self%delta*sin(k*xu)/ scaling%b
          bz = b0*self%delta*cos(k*xu)/ scaling%b
          uy = self%delta*sin(k*x)    / scaling%u
          uz = self%delta*cos(k*x)    / scaling%u

          ex = self%delta*cos(k*xu)/scaling%u*by - &
               self%delta*sin(k*xu)/scaling%u*bz
          ey = ux * b0*self%delta*cos(k*x)/scaling%b - &
               uz * bx
          ez = uy * bx - ux*b0*self%delta*sin(k*x)/scaling%b
        end if
        !
        if (trim(self%solver(1:6)) == 'ramses') then
          ff(ix,iy,iz,idx%e) = pgas/(self%gamma-1.0) & 
                             + 0.5*d*(ux**2 + uy**2) &
                             + 0.5 * (Bx**2 + By**2 + Bz**2) / scaling%e  
        else if (trim(self%solver) == 'stagger2_mhd_patch') then
          ff(ix,iy,iz,idx%s) = (log(pgas) - log(d) * self%gamma) &
                                * d / (self%gamma - 1.0) / scaling%e
        else if (trim(self%solver) == 'stagger2e_mhd_patch') then
          ff(ix,iy,iz,idx%e) = pgas / (self%gamma - 1.0) / scaling%e
        else if (trim(self%solver) == 'stagger2e_pic') then
          ff(ix,iy,iz,idx%e) = pgas / (self%gamma - 1.0) / scaling%e
        end if
        !
        ff(ix,iy,iz,idx%d ) = d
        ff(ix,iy,iz,idx%px) = ux * d 
        ff(ix,iy,iz,idx%py) = uy * d 
        ff(ix,iy,iz,idx%pz) = uz * d 
        ff(ix,iy,iz,idx%bx) = bx 
        ff(ix,iy,iz,idx%by) = by 
        ff(ix,iy,iz,idx%bz) = bz 
        ff(ix,iy,iz,idx%ex) = ex 
        ff(ix,iy,iz,idx%ey) = ey 
        ff(ix,iy,iz,idx%ez) = ez
      end do
    end do
  end do
  B => ff(:,:,:,idx%bx:idx%bz)
  J => ff(:,:,:,idx%jx:idx%jz)
  J = curl_dn(m%d,B)

  
  call trace%end()
END SUBROUTINE circular

END MODULE initial_mod
