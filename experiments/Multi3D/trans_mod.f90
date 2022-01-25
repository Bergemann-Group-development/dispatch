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
      integer :: nnu ! number of frequencies in the line

      integer :: ired  ! start index where line is active in overall frequency grid
      integer :: iblue ! end index where line is active in overall frequency grid

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
    real(8),dimension(:),allocatable :: q ! frequencies [Doppler units]
    real(8),dimension(:),allocatable :: wq  ! frequency quadrature weight [Doppler units]

    real(8),dimension(:,:,:),allocatable :: wphi ! normalization such that
    ! \Sigma_mu=1^{nrays} \Sigma_nu=1^{nnu} phi(nu,v(mu)) = 1
    ! in case of strong velocities this does not automatically follow from
    ! the independent normalization of the angle and frequency quadrature

    real(8),dimension(:,:,:),allocatable :: adamp ! voigt damping parameter
    real(8),dimension(:,:,:,:),allocatable :: rho ! prd rho
    real(8),dimension(:,:,:,:),allocatable :: Jrest ! Jnu in rest frame of atom
    real(8),dimension(:,:,:),  allocatable :: Qe ! prd Q elastic collision rate, used to compute coherency fraction

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

  contains
      procedure :: opacities => line_opacities
      procedure :: profile => line_profile

  end type line_type

  type, extends(transition_type) ::  cont_type

    character(len=20) :: bf_type ! either 'detail' or 'hydrogenic'

    real(8) :: nu0 ! ionization edge frequency
    real(8) :: numax ! maximum frequency 
    real(8) :: alpha0 ! cross section at ionization threshold [cm^2]
    real(8),dimension(:),allocatable :: alpha ! cross sections [cm^2]    

  contains
      procedure :: opacities => cont_opacities
      procedure :: rates => cont_rates
  end type cont_type

contains

  !====================================================================================
  pure subroutine cont_opacities(self, alpha, twohnu3_c2, n, nstar, blzfac, eta, chi)
    class(cont_type), intent(in) :: self
    real(8), dimension(:), intent(in) :: alpha, twohnu3_c2
    real(4), dimension(:,:,:,:), intent(in) :: n, nstar
    real(8), dimension(:,:,:,:), intent(in) :: blzfac
    real(8), dimension(:,:,:,:), intent(inout) :: eta, chi
    real(8), dimension(:,:,:), allocatable :: fakt
    integer :: f

    fakt = nstar(:,:,:,self%i) / nstar(:,:,:,self%j) * n(:,:,:,self%j)

    do f=1,size(twohnu3_c2)
      eta(:,:,:,f) = alpha(f) * fakt * blzfac(:,:,:,f)
      chi(:,:,:,f) = alpha(f) * n(:,:,:,self%i) - eta(:,:,:,f)
      eta(:,:,:,f) = twohnu3_c2(f) * eta(:,:,:,f)
    end do
  end subroutine cont_opacities
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure subroutine line_opacities(self, nu, hnu_4pi, n, vel_dop, doppler_factor, &
                                phi, eta, chi)
    class(line_type), intent(in) :: self
    real(8), dimension(:), intent(in) :: nu, hnu_4pi
    real(4), dimension(:,:,:,:), intent(in) :: n
    real(8), dimension(:,:,:), intent(in) :: vel_dop, doppler_factor
    real(8), dimension(:,:,:,:), intent(inout) :: phi, eta, chi
    real(8), dimension(:,:,:), allocatable :: scr1, scr2
    integer :: f

    call self%profile(nu, hnu_4pi, vel_dop, doppler_factor, phi)

    scr1 = self%Aji * n(:,:,:,self%j)
    scr2 = self%Bij * n(:,:,:,self%i) - self%Bji * n(:,:,:,self%j)

    do f=1, size(nu)
      eta(:,:,:,f) = phi(:,:,:,f) * scr1
      chi(:,:,:,f) = phi(:,:,:,f) * scr2
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
  pure subroutine line_rates(self, wnumu, nu, hnu_4pi, vel_dop, doppler_factor, &
                              phi, Inu, Rij, Rji)
      class(line_type), intent(in) :: self
      real(8), dimension(:), intent(in) :: wnumu, nu, hnu_4pi
      real(4), dimension(:,:,:,:), intent(in) :: Inu
      real(8), dimension(:,:,:,:), intent(inout) :: phi
      real(8), dimension(:,:,:), intent(in) :: vel_dop, doppler_factor
      real(8), dimension(:,:,:), intent(out) :: Rij, Rji
      real(8), dimension(:), allocatable :: wAji, wBij
      integer :: f

      Rij = 0
      Rji = 0

      call self%profile(nu, hnu_4pi, vel_dop, doppler_factor, phi)

      wAji = wnumu(f) * self%Aji
      wBij = wnumu(f) * self%Bij

      do f=1, size(wnumu)
        Rij = Rij + wBij(f) * Inu(:,:,:,f) * phi(:,:,:,f)
      end do

      Rji = Rij

      do f=1, size(wnumu)
        Rji = Rji + wAji(f) * phi(:,:,:,f)
      end do
  end subroutine line_rates
  !------------------------------------------------------------------------------------

  !====================================================================================
  pure subroutine line_profile(self, nu, nu_fac, vel_dop, doppler_factor, phi)
      class(line_type), intent(in) :: self
      real(8), dimension(:), intent(in) :: nu, nu_fac
      real(8), dimension(:,:,:), intent(in) :: vel_dop, doppler_factor
      real(8), dimension(:,:,:,:), intent(inout) :: phi
      real(8), dimension(:,:,:), allocatable :: norm
      integer :: f

      !if (self%io > 1) then ! io = nterm
      !   print *, 'Hyperfine-splitting not implemented yet. Stopping...'
      !   stop
      !endif

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
        phi(:,:,:,f) = phi(:,:,:,f)*norm * nu_fac(f)
      end do
  end subroutine line_profile
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
      real(8), intent(in) :: a,vs
      real(8) :: h

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
      else
        m = 6
        if (v < 1.3) m = 1
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
  pure subroutine check_range(self, fs_in, fe_in, fs_out, fe_out, skip)
      class(transition_type), intent(in) :: self
      integer, intent(in) :: fs_in, fe_in
      integer, intent(out) :: fs_out, fe_out
      logical, intent(out) :: skip

      skip = .false.

      fs_out = max(self%ired,  fs_in)
      fe_out = min(self%iblue, fe_in)

      if (fs_out >= fe_out) skip = .true.
      if (fs_out >= fe_in)  skip = .true.
      if (fe_out <= fs_in)  skip = .true.
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
  subroutine freql(line, qnorm)
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
  
      class(line_type), intent(inout)  :: line
      real(8) :: qnorm
      real(8),dimension(:),allocatable :: q,wq
  
      real(8)  :: dq, x, x10, al10, a, half, xmax, dx,qmax,q0
      integer :: nnu, j
  
  
      q0=line%q0
      qmax=line%qmax
      nnu=line%nnu
      allocate(q(nnu))
      allocate(wq(nnu))
  
      if (qmax.le.q0) then
      ! linear spacing
      dq = 2 * qmax / (nnu-1)
      q(1)=-qmax
      do j=2,nnu
          q(j) = q(j-1) + dq
      end do
  
      else if(qmax.ge.0.0 .and. q0.ge.0.0) then
  
      ! constants
      al10=log(10.)
      a=10.**(q0+.5)
      half=0.5
      xmax=log10(a*max(half,qmax-q0-half))
  
      dx=2.*xmax/(nnu-1)
      do j=1,nnu
          x=-xmax+(j-1)*dx
          x10=10.**x
          q(j)=x+(x10-1./x10)/a
      end do
  
      end if
  
      ! convert to [Hz]
      line%nu = line%nu0*(1.0+q*qnorm*KM_TO_CM/cc)
  
  end subroutine freql
  !------------------------------------------------------------------------------------

  !====================================================================================
  subroutine freqc(cnt)
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

      class(cont_type), intent(inout) :: cnt
      
      integer :: j,nnu,ired,iblue
      real(8) :: dnu

      select case(cnt%bf_type)
      case('hydrogenic')
          nnu   = cnt%nnu
          ired  = cnt%ired
          iblue = cnt%iblue

          ! equidistant grid in nu, between nu0 and numax 
          cnt%nu(ired)  = cnt%nu0
          cnt%nu(iblue) = cnt%numax

          dnu = (cnt%numax-cnt%nu0) / (nnu-1)

          do j=ired,iblue
              cnt%nu(j)    = cnt%nu0 + (j-ired) * dnu
              cnt%alpha(j) = cnt%alpha0 * ( cnt%nu(j) / cnt%nu0 )**(-3)
          end do
      case('detail')
          ! do nothing, frequencies were computed in read_atom
      case default
          print *, 'freqc: bf_type ' // cnt%bf_type // ' not recognized. Stopping...'
          stop
      end select
  end subroutine freqc
  !------------------------------------------------------------------------------------

  end module trans_mod

! *******************************************************************
