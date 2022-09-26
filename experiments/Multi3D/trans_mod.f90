! *******************************************************************
module trans_mod
  use mconst
  use bbtrans_mod
  ! *******************************************************************
  !
  !     Contains variables connected to the atom(s) 
  !     included in the calculations
  !
  ! *******************************************************************
  implicit none

  type anisotropy_type
      real(8),dimension(:,:,:),allocatable :: J00,J20
      complex(8),dimension(:,:,:),allocatable :: J2p1,J2m1,J2p2,J2m2
  end type anisotropy_type

  type, extends(bbtrans_type) :: transition_type
      integer :: ntrans ! transition identifier
      integer :: j   ! upper level  
      integer :: i   ! lower level
      integer :: nq  ! number of input frequencies in the line

      real(8) :: upper_ev

      integer :: ired  ! start index where line is active in overall frequency grid
      integer :: iblue ! end index where line is active in overall frequency grid

      real(8),dimension(:),allocatable :: nu ! frequencies [Hz]
      real(8),dimension(:),allocatable :: wnu ! frequency quadrature weight [Hz]

      real(8),dimension(:,:,:),  allocatable :: rij,rji
      real(8),dimension(:,:,:),  allocatable :: qrad ! heating rate
      real(8),dimension(:,:,:),  allocatable :: gamma_ij,gamma_ji

  contains
      procedure :: overlaps
      procedure :: check_range
  end type transition_type
  
  type, extends(transition_type) :: line_type

    character(len=72) :: profile_type ! what kind of profile to use,
                                      ! for now only 'voigt' or
                                      ! 'hlinop' for hydrogen

    real(8)  :: lambda0       ! line center vacuum wavelength
                            ! [cm]
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
      procedure :: freql
      procedure :: interpolate_phi

  end type line_type

  type, extends(transition_type) ::  cont_type

    character(len=20) :: bf_type ! either 'detail' or 'hydrogenic'

    !real(8) :: nu0 ! ionization edge frequency
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
    real(4) :: alpha2
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
      alpha2 = alpha(f)! * self%alpha0
    do k=1, nz
    do j=1, ny
    do i=1, nx
      eta(i,j,k,f) = eta(i,j,k,f) + alpha(f) * fakt(i,j,k)
      chi(i,j,k,f) = chi(i,j,k,f) + alpha2 * n(i,j,k,self%i)
    end do
    end do
    end do
    end do
  end subroutine cont_opacities
  !------------------------------------------------------------------------------------

!  !====================================================================================
!  subroutine line_opacities(self, nu, hnu_4pi, n, vel_dop, doppler_factor, &
!                                phi, eta, chi)
!    class(line_type), intent(in) :: self
!    real(8), dimension(:), intent(in) :: nu, hnu_4pi
!    real(4), dimension(:,:,:,:), intent(in) :: n
!    real(4), dimension(:,:,:), intent(in) :: vel_dop, doppler_factor
!    real(4), dimension(:,:,:,:), intent(inout) :: phi, eta, chi
!    real(4), dimension(:,:,:), allocatable :: scr1, scr2
!    integer :: f
!    !----------------------------------------------------------------------------------
!    call self%profile(nu, vel_dop, doppler_factor, self%adamp, phi)
!
!    do f=1, size(nu)
!      phi(:,:,:,f) = phi(:,:,:,f) * hnu_4pi(f)
!    end do
!
!    scr1 = self%Aji * n(:,:,:,self%j)
!    scr2 = self%Bij * n(:,:,:,self%i) - self%Bji * n(:,:,:,self%j)
!
!    do f=1, size(nu)
!      eta(:,:,:,f) = eta(:,:,:,f) + phi(:,:,:,f) * scr1
!      chi(:,:,:,f) = chi(:,:,:,f) + phi(:,:,:,f) * scr2
!    end do
!  end subroutine line_opacities
!  !------------------------------------------------------------------------------------

  !====================================================================================
  pure subroutine cont_rates(self, walpha, fconst, nstar, blzfac, Inu, Rij, Rji)
      class(cont_type), intent(in) :: self
      real(8), dimension(:), intent(in) :: walpha, fconst
      real(4), dimension(:,:,:,:), intent(in) :: nstar
      real(8), dimension(:,:,:,:), intent(in) :: blzfac, Inu
      real(8), dimension(:,:,:), intent(out) :: Rij, Rji
      real(8), dimension(:,:,:), allocatable :: Rijf
      integer :: f
      !----------------------------------------------------------------------------------
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
              ! for overflow reasons
              !self%alpha(j) = ( self%nu(j) / self%nu0 )**(-3)
          end do
      case('detail')
          ! do nothing, frequencies were computed in read_atom
      case default
          print *, 'freqc: bf_type ' // self%bf_type // ' not recognized. Stopping...'
          stop
      end select
  end subroutine freqc
  !------------------------------------------------------------------------------------

  end module trans_mod

! *******************************************************************
