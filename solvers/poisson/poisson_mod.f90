!> =============================================================================
!> Methods for solving the Poisson equation; derived from routines developed by
!> Troels Haugbølle for the course Computational Astrophysics
!> =============================================================================
MODULE poisson_mod
  USE mpi_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE patch_mod
  USE mesh_mod

  implicit none
  private
  type, public:: poisson_t
    character(len=32):: kind
    real:: omega_sor
    procedure(sor_2nd), pointer :: update => null()
    ! for use with ZEUS SOR solver
    real(8):: rjac2 ! spectral radius
    real(kind=KindScalarVar), dimension(:), allocatable:: dx1m2i, dx2m2i, dx1p2i, dx2p2i
    real(kind=KindScalarVar), dimension(:,:,:), allocatable:: eta
  contains
    procedure:: init
    procedure:: check => check_poisson
  end type

  real(8), parameter:: pi = asin(1.0) * 2.0

CONTAINS

!> =============================================================================
!> Initialize poisson_params
!> =============================================================================
SUBROUTINE init (self, patch, selfgravity)
  implicit none
  class(poisson_t)  :: self
  class(patch_t)    :: patch
  character(len=12) :: solver='sor_2nd'
  integer           :: iostat
  real,    save     :: omega_sor=0.0
  logical, save     :: first_time = .true.
  logical           :: selfgravity
  namelist /poisson_params/ solver, omega_sor
  !.............................................................................
  !$omp critical (poisson_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, poisson_params, iostat=iostat)
    if (io%master) write (io%output, poisson_params)
  end if
  !$omp end critical (poisson_cr)
  select case(trim(solver))
  case('sor_2nd')
    self%omega_sor = omega_sor
    call init_sor_2nd(self, patch)
    self%update => sor_2nd
    self%kind = solver
  case('zeus_sor')
    call init_zeus_sor(self, patch)
    self%update => zeus_sor
    self%kind = solver
  case default
    if (io%master) then
      write(*,*) 'Poisson solver '//trim(solver)//' is unknown'
      write(*,*) 'Valid solvers are: sor_2nd, zeus_sor'
    end if
    selfgravity = .false.
  end select
END SUBROUTINE init

!> =============================================================================
!> Initialize SOR factor
!> =============================================================================
SUBROUTINE init_sor_2nd (self, patch)
  implicit none
  class(poisson_t)  :: self
  class(patch_t)    :: patch
  !.............................................................................
  ! Calculate optimal guess for the SOR prefactor depending
  ! dimensionality, ds, and grid points n:
  !
  !    delta = (1/sum(1/ds^2)) sum(1/ds^2 * cos(pi/n))
  !    omega_sor = 2 / (1 + (1 - delta^2)^0.5)
  !
  ! For an even dimensionality grid
  !
  !    delta = cos(pi / n)
  !
  ! which for n > 10 is well approximated by
  !
  !    delta ~ 1 - 1/2 (pi/n)**2
  !
  ! and then
  !
  !    omega_sor = 2 / (1 + pi/n)
  !
  integer           :: i
  real              :: delta, pref
  real, parameter   :: pi=3.1415926536
  character(len=16) :: comment
  logical, save     :: first_time = .true.
  !-----------------------------------------------------------------------------
  if (self%omega_sor==0.) then
    delta = 0.0; pref = 0.0
    do i=1,3
      if (patch%n(i) > 1) then
        pref  = pref + 1.
        delta = delta + cos(pi / patch%n(i))
      endif
    enddo
    delta = delta / pref
    self%omega_sor = 2.0 / (1. + sqrt(1.0 - delta**2))
    comment =' (computed)'
  else
    comment =' (from namelist)'
  endif
  if (self%omega_sor < 0. .or. self%omega_sor>2) call mpi%abort( &
    'omega_sor in Laplace_mod has to be between 0 and 2 for stability')
  !$omp critical (poisson_cr2)
  if (first_time .and. io%master) then
    print '(a,f6.4,a)', 'poisson_mod: SOR prefactor set to: ', self%omega_sor, comment
    first_time = .false.
  endif
  !$omp end critical (poisson_cr2)
END SUBROUTINE init_sor_2nd

!> =============================================================================
!> 3pt 2nd order Red-Black Gauss-Seidel iteration with SOR acceleration
!>
!> Make one iterative update for
!>
!>   (2/dx^2 + 2/dy^2 + 2/dz^2) phi_new_ijk =
!>                                 (phi_i+1,j,k + phi_i-1,j,k)/dx^2
!>                               + (phi_i,j+1,k + phi_i,j-1,k)/dy^2
!>                               + (phi_i,j,k+1 + phi_i,j,k-1)/dz^2
!>                               - FourPiGrho
!>
!> =============================================================================
SUBROUTINE sor_2nd (self, patch, phi, resid, FourPiGrho, d, error, detailed_timer)
  class(poisson_t):: self
  class(patch_t)  :: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi, resid
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: FourPiGrho, d
  real   :: error
  logical:: detailed_timer
  integer:: i, j, k, sweep, ss, ione, jone, kone
  real   :: a, b, b1, b2, b3, ustar
  integer, save                  :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('poisson_t%sor_2nd_abs', itimer=itimer)
  associate(n=>patch%n,ds=>patch%ds,l=>patch%li,u=>patch%ui,gn=>patch%gn,omega_sor=>self%omega_sor)

  ione = 1
  jone = 1
  kone = 1
  if (n(1) <= 1) ione = 0
  if (n(2) <= 1) jone = 0
  if (n(3) <= 1) kone = 0

  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1) > 1) a = a + 2. / ds(1)**2
  if (n(2) > 1) a = a + 2. / ds(2)**2
  if (n(3) > 1) a = a + 2. / ds(3)**2
  a = 1. / a

  error = 0.0
  do sweep=0,1
    ! Add up contributions
    b1 = a / ds(1)**2 * ione
    b2 = a / ds(2)**2 * jone
    b3 = a / ds(3)**2 * kone
    b = 2.0 * (b1 + b2 + b3)
    do k=l(3),u(3)
      do j=l(2),u(2)
        ss = modulo(j + k + sweep,2)
        do i=l(1)+ss,u(1),2
          ustar = b1 * (phi(i+ione,j,k) + phi(i-ione,j,k)) &
                + b2 * (phi(i,j+jone,k) + phi(i,j-jone,k)) &
                + b3 * (phi(i,j,k+kone) + phi(i,j,k-kone)) &
                - a * FourPiGrho(i,j,k)
          phi(i,j,k)   = (1. - omega_sor) * phi(i,j,k) + omega_sor  * ustar
          resid(i,j,k) = (phi(i,j,k) - ustar) / a
          error = max(error, abs(ustar - b * phi(i,j,k)) / d(i,j,k))
        end do
      end do
    end do
  end do
  error = error / a
  end associate

  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE sor_2nd

!> =====================================================================
!> Initialise variables for the ZEUS SOR solver which are independent of
!> the MHD update and iterative calls to the SOR solver.
!> See accompanying `zeus_sor` procedure.
!> =====================================================================
SUBROUTINE init_zeus_sor (self, patch)
  class(poisson_t):: self
  class(patch_t):: patch
  class(mesh_t), pointer:: m1, m2, m3
  real(8):: q1, q2, q3, etaik, etajk, etakk, dx1a2i, dx2a2i, dx3a2i, dsi(3)
  integer:: i, j, k
  !---------------------------------------------------------------------
  call trace%begin('poisson_t%init_zeus_sor')

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  dsi = 1.0 / patch%ds

  allocate (self%dx1m2i(m1%gn), self%dx1p2i(m1%gn), self%dx2m2i(m2%gn), self%dx2p2i(m2%gn))
  allocate (self%eta(m1%gn,m2%gn,m3%gn))

  ! Compute grid-related quantities that don't change from one cycle to the next.
  ! These are required for curvilinear coords and go 1.0/ds**2 in Cartesian coords.
  ! Uniform cell sizes are assumed.
  self%dx1m2i(:) = 0.0
  self%dx1p2i(:) = 0.0
  self%dx2m2i(:) = 0.0
  self%dx2p2i(:) = 0.0
  do i=m1%li,m1%ui
    self%dx1m2i(i) = m1%h2f(i  ) * m1%h31f(i  ) / m1%h2c(i  ) / m1%h31c(i  ) * dsi(1)**2
    self%dx1p2i(i) = m1%h2f(i+1) * m1%h31f(i+1) / m1%h2c(i  ) / m1%h31c(i  ) * dsi(1)**2
  end do
  do j=m2%li,m2%ui
    self%dx2m2i(j) = m2%h32f(j  ) / m2%h32c(j) * dsi(2)**2
    self%dx2p2i(j) = m2%h32f(j+1) / m2%h32c(j) * dsi(2)**2
  end do

  !  Evaluate quantities needed by the SOR iteration loop, but are independent
  !  of phi.  These need only be calculated once for each patch. The Chebyshev
  !  acceleration factor, `omega_sor`, is also initialised.
  !
  !  Note from D. Clarke: omega_sor is started off at omega(n+1/2) in the series given
  !  by Press et al (page 859), and not omega(0) = 1.  This seems to allow
  !  for more rapid convergence.
  etakk  = 0.0
  etajk  = 0.0
  etaik  = 0.0
  q1     = 0.0
  q2     = 0.0
  q3     = 0.0
  dx1a2i = 0.0
  dx2a2i = 0.0
  dx3a2i = 0.0

  if (patch%n(3) > 1) then
    dx3a2i = dsi(3)**2
    etakk  =-2.0 * dx3a2i
    q3     = cos(pi / patch%n(3)) * dx3a2i
  end if
  if (patch%n(2) > 1) then
    dx2a2i = dsi(2)**2
    q2 = cos(pi / patch%n(2)) * dx2a2i
  end if
  if (patch%n(1) > 1) then
    dx1a2i = dsi(1)**2
    q1 = cos(pi / patch%n(1)) * dx1a2i
  end if
  self%rjac2 = ( ( q1 + q2 + q3 ) / ( dx1a2i + dx2a2i + dx3a2i ) )**2 ! Jacobi spectral radius squared
  self%omega_sor = 1.0 / ( 1.0 - 0.5 * self%rjac2 )
  do k=m3%li,m3%ui
    do j=m2%li,m2%ui
      if (patch%n(2) > 1) etajk =-self%dx2m2i(j) - self%dx2p2i(j)
      do i=m1%li,m1%ui
        if (patch%n(1) > 1) etaik =-self%dx1m2i(i) - self%dx1p2i(i)
        self%eta(i,j,k) = etaik &
                        + etajk / m1%h2c (i)**2 &
                        + etakk / (m1%h31c(i) * m2%h32c(j))**2
      end do
    end do
  end do

  call trace%end
END SUBROUTINE init_zeus_sor

!> =====================================================================
!> ZEUS-3D's version of SOR for solving Poisson's eqn.
!> Extended for curvilinear coords.
!> Originally written by Chris Howard in 2003 and modified by Alexander
!> Men'shchikov (2003), David Clarke (2004) and Nick MacDonald (2005).
!> Ported to DISPATCH by JPR (2017).
!      The SOR algorithm from Numerical Recipes (Press et al.), has been
  !  generalised for curvilinear geometries (XYZ, ZRP, RTP) and with any
  !  symmetry (1-D, 2-D, or full 3-D).
!> =====================================================================
SUBROUTINE zeus_sor (self, patch, phi, resid, FourPiGrho, d, anorm, detailed_timer)
  class(poisson_t):: self
  class(patch_t):: patch
  class(mesh_t), pointer:: m1, m2, m3
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi, resid
  real(kind=KindScalarVar), dimension(:,:,:), pointer, intent(in):: FourPiGrho, d
  real   :: anorm
  logical:: detailed_timer
  integer, save:: itimer = 0
  integer:: ione, jone, kone, i, j, k, sweep, sw1, sw2, sw3 &
          , inci, incj, inck, swi, swj, swk
  real(8):: da, dx3a2i
  !-----------------------------------------------------------------------
  if (detailed_timer) call trace%begin('poisson_mod%zeus_sor', itimer=itimer)
  associate (l => patch%li, u => patch%ui)

  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)

  ione = 1
  jone = 1
  kone = 1
  if (patch%n(1) <= 1) ione = 0
  if (patch%n(2) <= 1) jone = 0
  if (patch%n(3) <= 1) kone = 0
  ! loop strides
  inci = 1 + ione
  inck = 1 + (1-ione) * (1-jone)
  incj = 3 - inci*inck

  !  The checkerboarding described in Press et al. is performed on the
  !  innermost do-loop only, which is normally the i-loop.  However, if
  !  the problem has inherant i-symmetry, it is necessary to perform the
  !  checkerboarding on the j-loop, and if there is both i-symmetry and
  !  j-symmetry, on the k-loop.  The variables inc* and sw*, *=i,j,k, are
  !  introduced to perform the checkerboarding on the appropriate loop.
  anorm = 0.0d0
  dx3a2i = 1.0 / patch%ds(3)**2
  do sweep=1,2
    sw1 = sweep - 1
    sw2 = sw1
    swk = (inck-1) * sw1
    do k=l(3)+swk,u(3),inck
      sw3 = sw2
      swj = (incj-1) * sw2
      do j=l(2)+swj,u(2),incj
        swi = (inci-1) * sw3
        do i=l(1)+swi,u(1),inci
          da = self%eta(i,j,k) * phi(i,j,k) - FourPiGrho(i,j,k) &
             + ( phi(i+ione,j,k) * self%dx1p2i(i) + phi(i-ione,j,k) * self%dx1m2i(i) ) &
               * ione &
             + ( phi(i,j+jone,k) * self%dx2p2i(j) + phi(i,j-jone,k) * self%dx2m2i(j) ) &
               / m1%h2c(i)**2 * jone &
             + ( phi(i,j,k+kone) + phi(i,j,k-kone) ) * dx3a2i &
               / m1%h31c(i)**2 / m2%h32c(j)**2 * kone
          resid(i,j,k) = da
          phi(i,j,k)   = phi(i,j,k) - self%omega_sor / self%eta(i,j,k) * da
          anorm        = anorm + abs(da)
        end do
        sw3 = 1 - sw3
      end do
      sw2 = 1 - sw2
    end do
    !  Optimise omega using "Chebyshev acceleration".
    self%omega_sor = 1.0 / ( 1.0 - self%rjac2 * 0.25 * self%omega_sor )
  end do

  end associate
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE zeus_sor

!> =============================================================================
!> Given a grav. pot. and a density, calculate Poisson's equation and return
!> the residual.
!> =============================================================================
SUBROUTINE check_poisson (self, patch, FourPiG, resid)
  class(poisson_t):: self
  class(patch_t):: patch
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: phi
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(kind=KindScalarVar), dimension(patch%gn(1),patch%gn(2),patch%gn(3)):: resid
  real:: FourPiG
  integer:: i, j, k, ione, jone, kone
  real, dimension(3):: dsi2
  !.............................................................................
  call trace%begin('experiment_t%check_poisson')
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)

  ione = 1
  jone = 1
  kone = 1
  if (patch%n(1) <= 1) ione = 0
  if (patch%n(2) <= 1) jone = 0
  if (patch%n(3) <= 1) kone = 0
  dsi2 = 1.0 / patch%ds**2

  resid = 0.0
  associate (l => patch%li, u => patch%ui)
  do k=l(3),u(3)
    do j=l(2),u(2)
      do i=l(1),u(1)
        resid(i,j,k) = (phi(i+ione,j,k) - 2.0 * phi(i,j,k) + phi(i-ione,j,k)) * dsi2(1) * ione &
                     + (phi(i,j+jone,k) - 2.0 * phi(i,j,k) + phi(i,j-jone,k)) * dsi2(2) * jone &
                     + (phi(i,j,k+kone) - 2.0 * phi(i,j,k) + phi(i,j,k-kone)) * dsi2(3) * kone &
                     - FourPiG * d(i,j,k)
        resid(i,j,k) = abs(resid(i,j,k))
      end do
    end do
  end do
  end associate

  print*,"Max. resid.:",maxval(resid)

  call trace%end
END SUBROUTINE check_poisson

END MODULE poisson_mod
