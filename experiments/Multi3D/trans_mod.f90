! *******************************************************************
module trans_mod
  use mconst
  ! *******************************************************************
  !
  !     Contains variables connected to the atom(s) 
  !     included in the calculations
  !
  ! *******************************************************************
  implicit none

  ! *******************************************************************


  type anisotropy_type
      real(8),dimension(:,:,:),allocatable :: J00,J20
      complex(8),dimension(:,:,:),allocatable :: J2p1,J2m1,J2p2,J2m2
  end type anisotropy_type

  type transition_type
      integer :: ntrans ! transition identifier
      integer :: j   ! upper level  
      integer :: i   ! lower level
      integer :: nq  ! number of input frequencies in the line
      integer :: nnu ! number of frequencies in the line

      integer :: ired  ! start index where line is active in overall frequency grid
      integer :: iblue ! end index where line is active in overall frequency grid

      integer :: lower_ion ! ionisation stage of lower level
      real(8) :: lower_ev, upper_ev

      real(8),dimension(:),allocatable :: nu ! frequencies [Hz]
      real(8),dimension(:),allocatable :: wnu ! frequency quadrature weight [Hz]

      real(8),dimension(:,:,:),  allocatable :: rij,rji
      real(8),dimension(:,:,:),  allocatable :: qrad ! heating rate
      real(8),dimension(:,:,:),  allocatable :: gamma_ij,gamma_ji  
      logical :: shares_none  = .true.
      logical :: shared_upper = .false.
      logical :: shared_lower = .false.
      logical :: shared_cross = .false.

  contains
      procedure :: overlaps
      procedure :: name
      procedure :: check_range
  end type transition_type
  
  type, extends(transition_type) :: line_type

    character(len=72) :: profile_type ! what kind of profile to use,
                                      ! for now only 'voigt' or
                                      ! 'hlinop' for hydrogen

    real(8)  :: ga,gw,gq    ! line broadening (damping) constants
    real(8)  :: lambda0       ! line center vacuum wavelength
                            ! [cm]
    real(8)  :: nu0         ! line center frequency [Hz]
    real(8)  :: cnt_ev      ! energy of continuum above
    logical  :: hydrogen
    real(8),dimension(:),allocatable :: q ! frequencies [Doppler units]
    real(8),dimension(:),allocatable :: wq  ! frequency quadrature weight [Doppler units]

    real(4),dimension(:,:,:), pointer, contiguous :: wphi ! normalization such that
    ! \Sigma_mu=1^{nrays} \Sigma_nu=1^{nnu} phi(nu,v(mu)) = 1
    ! in case of strong velocities this does not automatically follow from
    ! the independent normalization of the angle and frequency quadrature

    real(4),dimension(:,:,:), pointer, contiguous :: adamp ! voigt damping parameter
    real(4),dimension(:,:,:), pointer, contiguous :: chi_fac, eta_fac
    real(8),dimension(:,:,:,:),allocatable :: rho ! prd rho
    real(8),dimension(:,:,:,:),allocatable :: Jrest ! Jnu in rest frame of atom
    real(8),dimension(:,:,:),  allocatable :: Qe ! prd Q elastic collision rate, used to compute coherency fraction
    real(4),dimension(:,:,:,:), pointer :: phi

    ! optional radiation anisotropy tensor
    type(anisotropy_type) :: an

    integer :: luprd=0 ! LUN for gii file
    character(len=72) :: giifile=' '

    real(8) :: Aji,Bji,Bij  ! Einstein coefficients, Bij is the
                            ! transition from the lower to the upper
                            ! level
    real(8) :: f    ! oscillator strength
    real(8) :: qmax ! width of profile in doppler units
    real(8) :: q0   ! switch from linear to logarithmic frequency grid spacing
    real(8) :: Grat ! ratio g(i)/g(j)

    ! Added to deal with multi component lines
    ! IO > 1 leads to the following being read from the atom file
    integer :: io ! number of components in term (lifted from Multi 2.3)
    real(8), dimension(:)      , allocatable :: determ
    real(8), dimension(:)      , allocatable :: wterm 
    real(8), dimension(:)      , allocatable :: fterm 
    real(8), dimension(:)      , allocatable :: gaterm
    real(8), dimension(:)      , allocatable :: gwterm
    real(8), dimension(:)      , allocatable :: gqterm
    real(8), dimension(:,:,:,:), allocatable :: adterm

    integer :: nsub  ! number of subordinate lines in case of cross redistribution
    integer, dimension(:),allocatable :: krsub ! indices of subordinate lines

    logical :: active=.true.

  contains
      procedure :: opacities => line_opacities
      procedure :: profile => line_profile
      procedure :: freql
      procedure :: LinearStarkH
      procedure :: IonicH
      procedure :: Stark
      procedure :: VanderWaals
      procedure :: damping
      procedure :: estimate_qmax
      procedure :: interpolate_phi

  end type line_type

  type, extends(transition_type) ::  cont_type

    character(len=20) :: bf_type ! either 'detail' or 'hydrogenic'

    real(8) :: nu0 ! ionization edge frequency
    real(8) :: numax ! maximum frequency 
    real(8) :: alpha0 ! cross section at ionization threshold [cm^2]
    real(4),dimension(:),allocatable :: alpha ! cross sections [cm^2]    
    real(4), dimension(:,:,:), contiguous, pointer :: ni_over_nj_lte

  contains
      procedure :: opacities => cont_opacities
      procedure :: rates => cont_rates
      procedure :: freqc
  end type cont_type

contains

  !====================================================================================
  pure subroutine cont_opacities(self, alpha, n, nstar, eta, chi)
    class(cont_type), intent(in) :: self
    real(4), dimension(:), intent(in) :: alpha
    real(4), dimension(:,:,:,:), intent(in) :: n, nstar
    real(4), dimension(:,:,:,:), intent(inout) :: eta, chi
    real(4), dimension(:,:,:), allocatable :: fakt
    integer :: i,j,k,f, nf, nx, ny, nz
    !----------------------------------------------------------------------------------
    nx = size(n, 1)
    ny = size(n, 2)
    nz = size(n, 3)
    nf = size(eta,4)

    allocate(fakt(nx,ny,nz))
    !fakt = nstar(:,:,:,self%i) / nstar(:,:,:,self%j) * n(:,:,:,self%j)
    fakt = self%ni_over_nj_lte * n(:,:,:,self%j)

    !$omp simd
    do f=1, nf
    do k=1, nz
    do j=1, ny
    do i=1, nx
      eta(i,j,k,f) = eta(i,j,k,f) + alpha(f) * fakt(i,j,k)
      chi(i,j,k,f) = chi(i,j,k,f) + alpha(f) * n(i,j,k,self%i)
    end do
    end do
    end do
    end do
  end subroutine cont_opacities
  !------------------------------------------------------------------------------------

  !====================================================================================
  subroutine line_opacities(self, nu, hnu_4pi, n, vel_dop, doppler_factor, &
                                phi, eta, chi)
    class(line_type), intent(in) :: self
    real(8), dimension(:), intent(in) :: nu, hnu_4pi
    real(4), dimension(:,:,:,:), intent(in) :: n
    real(4), dimension(:,:,:), intent(in) :: vel_dop, doppler_factor
    real(4), dimension(:,:,:,:), intent(inout) :: phi, eta, chi
    real(4), dimension(:,:,:), allocatable :: scr1, scr2
    integer :: f

    call self%profile(nu, vel_dop, doppler_factor, phi)

    do f=1, size(nu)
      phi(:,:,:,f) = phi(:,:,:,f) * hnu_4pi(f)
    end do

    scr1 = self%Aji * n(:,:,:,self%j)
    scr2 = self%Bij * n(:,:,:,self%i) - self%Bji * n(:,:,:,self%j)

    do f=1, size(nu)
      eta(:,:,:,f) = eta(:,:,:,f) + phi(:,:,:,f) * scr1
      chi(:,:,:,f) = chi(:,:,:,f) + phi(:,:,:,f) * scr2
    end do
  end subroutine line_opacities
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure subroutine cont_rates(self, walpha, fconst, nstar, blzfac, Inu, Rij, Rji)
      class(cont_type), intent(in) :: self
      real(8), dimension(:), intent(in) :: walpha, fconst
      real(4), dimension(:,:,:,:), intent(in) :: nstar
      real(8), dimension(:,:,:,:), intent(in) :: blzfac, Inu
      real(8), dimension(:,:,:), intent(out) :: Rij, Rji
      real(8), dimension(:,:,:), allocatable :: Rijf
      integer :: f

      Rij = 0
      Rji = 0

      do f=1, size(walpha)
        Rijf = walpha(f) * Inu(:,:,:,f)
        Rij = Rij + Rijf
        Rji = Rji + (fconst(f) + Rijf) * blzfac(:,:,:,f)
      end do

      Rji = Rji * nstar(:,:,:,self%i) / nstar(:,:,:,self%j)
  end subroutine cont_rates
  !------------------------------------------------------------------------------------

  !====================================================================================
  subroutine line_profile(self, nu, vel_dop, doppler_factor, phi)
      class(line_type), intent(in) :: self
      real(8), dimension(:), intent(in) :: nu
      real(4), dimension(:,:,:), intent(in) :: vel_dop, doppler_factor
      real(4), dimension(:,:,:,:), intent(inout) :: phi
      real(4), dimension(:,:,:), allocatable :: norm
      integer :: f
      !----------------------------------------------------------------------------------
      allocate(norm, mold=doppler_factor)
      norm = 1 / (self%nu0 * doppler_factor)

      do f=1, size(nu)
        phi(:,:,:,f) = nu(f)*norm - vel_dop
      end do 

      select case(trim(self%profile_type))
        ! Single-component transitions
        ! Hyperfine-Splitting not implemented yet
      case ('doppler')
        phi = exp(-phi*phi)

      case ('voigt')
        do f=1, size(nu)
            phi(:,:,:,f) = voigt_function(self%adamp, phi(:,:,:,f))
        end do
      end select

      norm = norm / sqrt(pi)

      do f=1, size(nu)
        phi(:,:,:,f) = phi(:,:,:,f)*norm! * nu_fac(f)
      end do
  end subroutine line_profile
  !------------------------------------------------------------------------------------

  !====================================================================================
  subroutine interpolate_phi(self, nu, phi)
    use op_math_mod, only: searchsorted
    class(line_type), intent(in) :: self
    real(8), dimension(:), intent(in) :: nu
    real(4), dimension(:,:,:,:), intent(out) :: phi
    integer :: i,j,k,f, i1, i2, n(4)
    real(8) :: dx, dy, norm, newq
    !----------------------------------------------------------------------------------
    n = shape(phi)
    norm = 1 / self%nu0

    do f=1, n(4)
      newq = nu(f)*norm - 1
      i1 = searchsorted(self%q, newq)

      if (i1==-1) then
        phi(:,:,:,f) = 0
        cycle 
      endif

      i2 = i1+1
      dx = (newq - self%q(i1)) / (self%q(i2)-self%q(i1))

      if (dx < 1e-6) then
        phi(:,:,:,f) = self%phi(:,:,:,i1)
        cycle
      elseif (dx > 1-1e-6) then
        phi(:,:,:,f) = self%phi(:,:,:,i2)
        cycle
      endif

      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        dy = self%phi(i,j,k,i2) - self%phi(i,j,k,i1)
        phi(i,j,k,f) = dy*dx + self%phi(i,j,k,i1)
      end do
      end do
      end do
    end do
  end subroutine interpolate_phi
  !------------------------------------------------------------------------------------

  !====================================================================================
  elemental function voigt_function(a,vs) result(h)
      !---------------------------------------------------------------------
      !
      !  computes a voigt function  h = h(a,v)
      !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
      !  traving (landolt-b\rnstein, p. 449).
      !
      !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
      !
      !---------------------------------------------------------------------
      !
      implicit none
      real(4), intent(in) :: a,vs
      real(4) :: h

      integer :: k,m,i

      real(8),dimension(5) :: a1 
      real(8) :: v,u,v2,ex,quo,h1,pqs,h1p,h2p,h3p,h4p,psi
      real(8) :: a2,u2

      real(8),dimension(19),parameter :: ak = (/ &
        -1.12470432, -0.15516677,  3.28867591, -2.34357915,             &
        0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,&
        -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,&
        -0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278 /)

      real(8),parameter :: sqp = 1.772453851
      real(8),parameter :: sq2 = 1.414213562

      v = abs(vs)
      u = a + v
      v2 = v*v

      k = 1

      if (a == 0.0) then
        h=0.
        if(v2.lt.100.)h=exp(-v2)
        return

      elseif (a > 0.2) then
        if (a > 1.4 .or. u > 3.2) then
            a2 = a*a
            u = sq2*(a2 + v2)
            u2 = 1./(u*u)
            ! a gt 1.4  or  a + v gt 3.2
            h = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) + &
                  u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2))
            return
        else
            k = 2
        endif
      elseif (v >= 5.0) then
        h = a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp)
        return
      endif

      ex=0.
      if(v2 < 100.) ex = exp(-v2)

      quo = 1.
      if (v >= 2.4) then
        quo = 1./(v2 - 1.5)
        m = 11
        !h1 = 1./(v2 - 1.5) * ( 0.554153432 + v*( 0.278711796 + v*(-0.188325687 + v*( 0.042991293 + v*-0.003278278))))
      else
        m = 6
        !h1 = -4.48480194 + v*( 9.39456063 + v*(-6.61487486 + v*( 1.98919585 + v*-0.22041650)))
        if (v < 1.3) m = 1
        !h1 = -1.12470432 + v*(-0.15516677 + v*( 3.28867591 + v*(-2.34357915 + v* 0.42139162)))
      endif

      do i=1,5
        a1(i) = ak(m)
        m = m + 1
      end do

      h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))

      if (k==1) then
        ! a le 0.2  and v lt 5.
        h = h1*a + ex*(1. + a*a*(1. - 2.*v2))
        return
      endif
      
      pqs = 2./sqp
      h1p = h1 + pqs*ex
      h2p = pqs*h1p - 2.*v2*ex
      h3p = (pqs*(1. - ex*(1. - 2.*v2)) - 2.*v2*h1p)/3. + pqs*h2p
      h4p = (2.*v2*v2*ex - pqs*h1p)/3. + pqs*h3p
      psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
      ! 0.2 < a <= 1.4  and  a + v <= 3.2
      h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
      return

  end function voigt_function
  !------------------------------------------------------------------------------------

  !====================================================================================
  function estimate_qmax(self, temp, ne, nh6, tothi, awgt) result(qmax)
    class(line_type), intent(in)  :: self
    real(4), intent(in) :: temp, ne, nh6, tothi
    real(8), intent(in) :: awgt
    real(4) :: qmax
    real(4) :: a, h0, vmax, doppler_factor, ratio
    real(8), parameter :: sqp = 1.772453851
    real(8), parameter :: sq2 = 1.414213562
    real(8), parameter :: tol = 1d-4
    !----------------------------------------------------------------------------------
    doppler_factor = sqrt(2*bk / awgt * temp) / cc
    a = self%damping(temp, ne, nh6, tothi, doppler_factor, awgt)
    h0 = voigt_function(a, 0.0)

    if (a > 1.4) then
      vmax = a * 19*sq2 / (sqp * h0 * tol)
    else
      vmax = sqrt(a / (sqp * h0 * tol))
    endif

    ratio = voigt_function(a, vmax) / h0
    qmax = vmax * doppler_factor * cc / KM_TO_CM
    !print '(I6, 6(ES17.10))', self%ntrans, self%lambda0, a, self%f, qmax, ratio
  end function estimate_qmax
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure subroutine check_range(self, fs_in, fe_in, fs_out, fe_out, skip)
      class(transition_type), intent(in) :: self
      integer, intent(in) :: fs_in, fe_in
      integer, intent(out) :: fs_out, fe_out
      logical, intent(out) :: skip

      skip = .false.

      fs_out = max(self%ired,  fs_in)
      fe_out = min(self%iblue, fe_in)

      if (fs_out > fe_out) skip = .true.
      if (fs_out > fe_in)  skip = .true.
      if (fe_out < fs_in)  skip = .true.
  end subroutine check_range
  !------------------------------------------------------------------------------------

  !====================================================================================
  logical pure function overlaps(self, other)
      class(transition_type), intent(in) :: self, other

      overlaps = .false.

      if (self%ired > other%ired .and. self%ired < other%iblue) overlaps = .true.
      if (other%ired > self%ired .and. other%ired < self%iblue) overlaps = .true.

  end function overlaps
  !------------------------------------------------------------------------------------

  !====================================================================================
  character(len=10) function name(self)
      class(transition_type), intent(in) :: self

      select type(self)
      class is (cont_type)
        write(name, '(A4, I6)') 'cont', self%ntrans
      class is (line_type)
        write(name, '(A4, I6)') 'line', self%ntrans
      end select

  end function name
  !------------------------------------------------------------------------------------


  !====================================================================================
  subroutine freql(self, qnorm, symmetric)
      !***********************************************************************
  
      !  gives  frequency quadrature points and corresponding weights.
      !  one- or two-sided cases.  the points and weights are based on
      !  the trapezoidal rule applied to a mapping function  that maps
      !  equidistant  points in  x to  equidistant points in  q  in the
      !  core (q .lt. q0), and equidistant points in log(q) in the wings
      !  (q .gt. q0).
  
      !  nq     : total number of quadrature points                (in)
      !  qmax   : maximum value of frequency (doppler units)       (in)
      !  ind    : =1 for one-sided and =2 for two-sided cases      (in)
      !  q0     : transition from core to wing (doppler units)     (in)
      !  q      : quadrature points                               (out)
      !  wq     : quadrature weights                              (out)
  
      !  coded by: ake nordlund (oct-1981).
  
      !  revised by mats carlsson (jan-1984): new options
      !  qmax.le.q0   q linear in frequency
      !  qmax or q0.lt.0  q given and only weights are calculated
  
      !:
      !: freql  88-05-03  modifications: (mats carlsson)
      !:        added options gave incorrect weights for ind=1
      !:        this affected cooling rates for lines with frequency
      !:        quadrature given with qmax=q0 or qmax.lt.0 or q0.lt.0
      !:        the corrected weights are twice the old ones
      !:
      !***********************************************************************
  
      class(line_type), intent(inout)  :: self
      real(8) :: qnorm
      real(8),dimension(:),allocatable :: q,wq
  
      real(8)  :: dq, x, x10, a, half, xmax, dx,qmax,q0, start
      integer :: nnu, j
      logical :: symmetric
  
  
      q0=self%q0
      qmax=self%qmax
      nnu=self%nnu
      allocate(q(nnu))
      allocate(wq(nnu))
  
      if (qmax.le.q0) then
        ! linear spacing
        dq = qmax / (nnu-1)
        q(1) = 0

        if (symmetric) then
          dq = 2.*dq
          q(1) = -qmax
        endif

        do j=2,nnu
            q(j) = q(j-1) + dq
        end do
  
      else if(qmax.ge.0.0 .and. q0.ge.0.0) then
        a=10.**(q0+.5)
        half=0.5
        xmax=log10(a*max(half,qmax-q0-half))
    
        dx=xmax/(nnu-1)
        start=0

        if (symmetric) then
          dx = 2.*dx
          start = -xmax
        endif

        do j=1,nnu
            x=start+(j-1)*dx
            x10=10.**x
            q(j)=x+(x10-1./x10)/a
        end do
  
      end if
  
      ! convert to [Hz]
      self%nu = self%nu0*(1.0+q*qnorm*KM_TO_CM/cc)
  
  end subroutine freql
  !------------------------------------------------------------------------------------

  !====================================================================================
  subroutine freqc(self)
      !---------------------------------------------------------------------   
      !
      !  gives quadrature points and weights for continua
      !  the points are equidistant in frequency, the weights are
      !  trapezoidal
      !  calculates crossections from nn**-3 dependence if not given
      !  explicitly
      !
      !  nnu      total number of quadrature points    
      !  nu     quadrature points (Hz)                
      !  wnu      quadrature weights                  
      !
      !---------------------------------------------------------------------

      class(cont_type), intent(inout) :: self
      
      integer :: j,nnu,ired,iblue
      real(8) :: dnu

      select case(self%bf_type)
      case('hydrogenic')
          nnu   = self%nnu
          ired  = self%ired
          iblue = self%iblue

          ! equidistant grid in nu, between nu0 and numax 
          self%nu(ired)  = self%nu0
          self%nu(iblue) = self%numax

          dnu = (self%numax-self%nu0) / (nnu-1)

          do j=ired,iblue
              self%nu(j)    = self%nu0 + (j-ired) * dnu
              self%alpha(j) = self%alpha0 * ( self%nu(j) / self%nu0 )**(-3)
          end do
      case('detail')
          ! do nothing, frequencies were computed in read_atom
      case default
          print *, 'freqc: bf_type ' // self%bf_type // ' not recognized. Stopping...'
          stop
      end select
  end subroutine freqc
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure elemental function damping (self, temp, ne, nh6, tothi, doppler_factor, awgt) result(adamp)
    class(line_type), intent(in)  :: self
    real(4), intent(in)  :: temp, ne, nh6, tothi, doppler_factor
    real(8), intent(in)  :: awgt
    real(4) :: adamp
    integer :: iterm
    !----------------------------------------------------------------------------------
    iterm = 1 ! hyperfine splitting not implemented yet
    adamp = 0

    if (self%hydrogen) then
       call self%LinearStarkH(adamp,ne)
       call self%IonicH(adamp,temp,nh6,ne)
    endif

    call self%vanderWaals(adamp,tothi,iterm,ne,temp,awgt)
    call self%Stark(adamp,iterm,ne,temp)
    
    adamp = adamp + self%ga
    adamp = adamp / (4.0*PI*self%nu0*doppler_factor)
  end function damping
  !------------------------------------------------------------------------------------

  !====================================================================================
  ! Linear Stark broadening for hydrogen, see
  ! Sutton K. (1978), JQSRT, 20, 333
  !====================================================================================
  pure elemental subroutine LinearStarkH(self, gamma, ne)
    class(line_type), intent(in)  :: self
    real(4), intent(inout) :: gamma
    real(4), intent(in)    :: ne
    real(4) :: a1
    !----------------------------------------------------------------------------------
    if(self%j-self%i.eq.1) then
      a1=0.642
    else
      a1=1.000
    endif

    gamma = gamma + a1 * 0.6 * (self%j**2-self%i**2) * ne**two_third
  end subroutine LinearStarkH
  !------------------------------------------------------------------------------------

  !====================================================================================
  ! Computes quadratic Stark broadening, see The observation and
  ! analysis of stellar photospheres, D.F. Gray, 2005, Cambridge
  ! University Press, Chapter 11, equation 11.27
  !====================================================================================
  pure elemental subroutine Stark(self, gamma, interm, ne, temp)
    class(line_type), intent(in)  :: self
    real(4),        intent(in)    :: temp,ne
    real(4),        intent(inout) :: gamma
    integer,        intent(in)    :: interm
    real(4)                       :: gq,gslg,gs
    !----------------------------------------------------------------------------------
    if (self%io .le. 1) then
      gq = self%gq
    else
      gq = self%gqterm(interm)
    endif

    if (gq .ge. 0.0) then
      gs = gq * ne
    else
      gslg = 19.4 + 2./3. * gq + log10(ne*bk) + log10(temp) / 6.
      gs = 10.0 ** gslg
    end if
    
    gamma = gamma + gs
  end subroutine Stark
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure elemental subroutine IonicH (self, gamma, temp, nh, ne) 
    class(line_type), intent(in)  :: self
    real(4),intent(in)         :: temp,nh,ne
    real(4),intent(inout)      :: gamma

    integer :: i,j,k,il,jl
    real(8) :: fn,fn1,rd,rmax,rmin,dohm,gprot,xmr,xmp,xm,rratio,xn
    real(8) :: gion
    !----------------------------------------------------------------------------------
    il = self%i
    jl = self%j

    ! ionic broadening
    if( (jl-il >= 1) .and. (jl-il <= 3) .and. (il >= 4) .and. (il <= 20) ) then

      fn1 = 2.25 * il**2 * ( 1.5 + 0.271*log(2.0/3.0*il)+  1.92/il**0.422 )
      if(jl-il.eq.1) then
        fn=fn1
      else if(jl-il.eq.2) then
        fn=fn1*(9.1333/il+5.57333)
      else if(jl-il.eq.3) then
        fn=fn1*(38.6667/il+12.9667)
      endif

      ! xm is reduced mass of radiator and perturber
      xmr=1.0
      xmp=1.0
      xm=xmp*xmr/(xmp+xmr)

      rd=6.90*sqrt( temp/nh )
      rmax=min( rd , self%lambda0 / 2. / pi )

      rmin = sqrt(2./3.) * sqrt(fn) / 2.0/pi/em*hh*sqrt(2.0*xm*uu/pi/bk/temp)
      rratio=max(rmax/rmin,one)
      dohm=2.7e8*sqrt(uu*2*xm/temp)*nh*fn*(0.5+log(rratio))
      !  dohm is half-width in angular frequency units. convert to
      !  gamma (hwhm(ohmega)=gamma/2, hwhm(nu)=gamma/4pi)
      gprot=dohm

      !  ion broadening from other ions than protons. assume an atomi! 
      !  weight of 25 and calculate n from ne-np (assuming all electrons come
      !  from singly ionized atoms
      xn=max(ne-nh,zero)
      xmp=25.0
      xm=xmp*xmr/(xmp+xmr)
      if(xn.gt.0.0) then
        rd=6.90*sqrt(temp/xn)
        rmax=min(rd,self%lambda0/2./pi)
        rmin=sqrt(2./3.)*sqrt(fn)/2./pi/em*hh*sqrt(2*xm*uu/pi/bk/temp) 
        rratio=max(rmax/rmin,one)
        dohm=2.7e8*sqrt(uu*2*xm/temp)*xn*fn*(0.5+log(rratio))
        gion=dohm
      else
        gion=0.0
      endif

      gamma = gamma + gprot + gion
    endif
  end subroutine IonicH
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure elemental subroutine VanderWaals(self,gamma,tothi,interm,ne,temp,awgt)
    class(line_type), intent(in)  :: self
    real(4), intent(in)    :: temp,ne,tothi
    real(8), intent(in)    :: awgt
    real(4), intent(inout) :: gamma
    integer, intent(in)    :: interm
    integer                :: ic
    real(4)                :: alpha,sigma,gx,gammaf,gv,c1,c625,scr
    !----------------------------------------------------------------------------------
    ! compute van der waals broadening due to hydrogen
    ! following Barklem
    if (self%io .le. 1) then
      gv = self%gw
    else
      gv = self%gwterm(interm)
    endif

    if(gv.ge.20.0) then
      sigma=int(gv)*2.80028e-21
      alpha=gv-int(gv)
      gx=2.0-alpha*0.5
      gx=gx-1.0
      gammaf=1+(-.5748646+(.9512363+(-.6998588+ &
      (.4245549-.1010678*gx)*gx)*gx)*gx)*gx

      gv = (4.0/pi)**(alpha*0.5) * gammaf * 1.e4 * sigma
      c1= 21173.*(1./1.008+1./(awgt/uu))
      scr = sqrt(c1*temp)
      scr=gv*((scr/1.e4)**(1.-alpha))
      scr=scr*tothi*1.e6*2.

    elseif(gv.ne.0.0) then
      ! Unsold approximation
      ! Jul  5 2010 jorritl: I'm sceptical about the numerical 
      ! constant  1.283984e-12. According to Gray it should be 
      ! 0.3d-30^0.4 = 6.1780059e-13 and  8.411 should be 17
      ! I haven't found the original reference to Unsold 1955 yet...
      c1 = (1./(self%cnt_ev-self%upper_ev)**2 - 1./(self%cnt_ev-self%lower_ev)**2)**.4
      c625 = 1.283984e-12 * self%lower_ion**0.8 * c1 

      ! this is v_rel^3/5
      scr = ( 8.0 * bk / pi *  (1.0/mh + 1.0/awgt) * temp )**0.3
      scr = gv * 8.411 * scr * tothi * c625
    endif

    gamma = gamma + scr
  end subroutine VanderWaals
  !------------------------------------------------------------------------------------


  end module trans_mod

! *******************************************************************
