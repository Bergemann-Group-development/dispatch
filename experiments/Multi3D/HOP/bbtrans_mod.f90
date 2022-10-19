MODULE bbtrans_mod
  USE prec_mod
  USE op_math_mod, only: searchsorted
  USE trace_mod
  implicit none

  type :: bbtrans_type
    real(dp) :: nu0
    integer :: nnu ! number of frequencies in the line
    real(dp) :: ga=0, gw1=0, gw2=0, gq=0 ! line broadening (damping) constants
    real(dp) :: lower_ev, upper_g
    !real(dp) :: gf
    real(dp) :: Aji, Bji, Bij ! Einstein coefficients, Bij is the
                              ! transition from the lower to the upper
                              ! level
    logical :: active=.true., unsold
    integer :: iso_idx=0 ! isotope index for lines in line list
  contains
    procedure :: profile
    procedure :: damping => line_damping
    procedure :: dampingH => line_dampingH
    procedure :: add_lte_opacities
    procedure :: LinearStarkH
    procedure :: IonicH
    procedure :: prepare_Stark
    procedure :: prepare_VanderWaals
    procedure :: estimate_qmax
    procedure :: get_range
  end type
    
contains

  !====================================================================================
  SUBROUTINE profile(self, q, vproj, inv_df, adamp, phi)
    class(bbtrans_type), intent(in) :: self
    real(8), dimension(:), intent(in) :: q
    real(4), dimension(:,:,:), intent(in) :: vproj, inv_df, adamp
    real(4), dimension(:,:,:,:), intent(inout) :: phi
    real(4), dimension(:,:,:), allocatable :: psi, damp
    integer :: f
    !----------------------------------------------------------------------------------
    allocate(psi, damp, mold=inv_df)
    damp = adamp * inv_df * over_4pi / self%nu0
    psi = psi_func(damp)

    !DIR$ FORCEINLINE
    do f=1, size(q)
      phi(:,:,:,f) = (q(f) - vproj) * inv_df
      phi(:,:,:,f) = voigt_function2(damp, phi(:,:,:,f), psi)
      phi(:,:,:,f) = phi(:,:,:,f) * inv_df
    end do 
  END SUBROUTINE profile
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE add_lte_opacities(self, nu, vproj, inv_df, EK_over_temp, inv_part, adamp, eta, chi)
    class(bbtrans_type), intent(in) :: self
    real(8), dimension(:), intent(in), contiguous :: nu
    real(4), dimension(:,:,:), intent(in), contiguous :: vproj, inv_df, inv_part, adamp
    real(8), dimension(:,:,:), intent(in), contiguous :: EK_over_temp
    real(4), dimension(:,:,:,:), intent(inout), contiguous :: eta, chi
    real(4), dimension(:,:,:), allocatable :: eta_fac, chi_fac, norm, damp, psi, phi
    real(8), dimension(:), allocatable :: q
    real(4) :: ex, ni, nj, upper_ev
    integer :: i,j,k,f, m(4)
    integer, save :: itimer=0, itimer2=0
    !----------------------------------------------------------------------------------
    call trace%begin('bbtrans_type%add_lte_opacities', itimer=itimer)
    allocate(eta_fac, chi_fac, norm, psi, damp, phi, mold=vproj)
    allocate(q, mold=nu)
    m = shape(eta)

    q = nu / self%nu0 - 1
    norm = inv_df / self%nu0 * over_4pi
    damp = adamp * norm
    psi = psi_func(damp)

    upper_ev = self%lower_ev + H_E * self%nu0

    do k=1, m(3)
    do j=1, m(2)
    do i=1, m(1)
      ex = self%lower_ev * EK_over_temp(i,j,k)
      ni = inv_part(i,j,k)
      if (ex < 80.) ni = ni * exp(-ex)
      ex = upper_ev * EK_over_temp(i,j,k)
      nj = self%upper_g * inv_part(i,j,k)
      if (ex < 80.) nj = nj * exp(-ex)

      eta_fac(i,j,k) = norm(i,j,k) *  self%Aji * nj
      chi_fac(i,j,k) = norm(i,j,k) * (self%Bij * ni - self%Bji * nj)
    end do
    end do
    end do

    call trace%end(itimer)
    call trace%begin('bbtrans_type%voigt_function', itimer=itimer2)
    
    !DIR$ FORCEINLINE
    do f=1, m(4)
      phi = (q(f) - vproj) * inv_df
      phi = voigt_function2(damp, phi, psi)
      phi = phi * nu(f)

      do k=1, m(3)
      do j=1, m(2)
      do i=1, m(1)
          eta(i,j,k,f) = eta(i,j,k,f) + eta_fac(i,j,k) * phi(i,j,k)
          chi(i,j,k,f) = chi(i,j,k,f) + chi_fac(i,j,k) * phi(i,j,k)
      end do
      end do
      end do
    end do
    call trace%end(itimer2)
  END SUBROUTINE add_lte_opacities
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION psi_func(a) result(psi)
    real(4), intent(in) :: a
    real(4) :: psi
    !----------------------------------------------------------------------------------
    psi = (0.979895023 + a*(-0.962846325 + a*(0.532770573 + a*-0.122727278))) * 0.5641895835
  END FUNCTION psi_func
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION voigt_function2(a,vs,psi) result(h)
    !---------------------------------------------------------------------
    !  computes a voigt function  h = h(a,v)
    !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
    !  traving (landolt-b\rnstein, p. 449).
    !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
    !---------------------------------------------------------------------
    real(4), intent(in) :: a,vs,psi
    real(4) :: h
    real(4) :: v,u,v2,v4,ex,h1,h1p,h2p,h3p,h4p
    real(4) :: a2,u2,v2ex

    real(4), parameter :: pqs  = 1.128379167  ! 2./pi**0.5
    real(4), parameter :: isqp = 0.5641895835 ! 1/pi**0.5
    real(4), parameter :: ipi  = 0.3183098862 ! 1/pi
    !----------------------------------------------------------------------------------
    v = abs(vs)
    u = a + v
    v2 = v*v

    if (a == 0.0) then
      if (v2 < 100.) then
        h = exp(-v2) * isqp
      else
        h = 0.
      endif

    elseif (a > 0.2) then
      if (a > 1.4 .or. u > 3.2) then
        a2 = a*a
        u = 1 / (a2 + v2)
        u2 = 0.5 * u*u
        ! a gt 1.4  or  a + v gt 3.2
        h = ipi*a*u * (1. + u2*((3.*v2 - a2) + u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2)))

      else
        if (v >= 2.4) then
          h1 = 1./(v2 - 1.5) * (0.554153432 + v*(0.278711796 + v*(-0.188325687 + v*(0.042991293 + v*-0.003278278))))
        else
          if (v >= 1.3) then
            h1 = -4.48480194 + v*(9.39456063 + v*(-6.61487486 + v*(1.98919585 + v*-0.22041650)))
          else
            h1 = -1.12470432 + v*(-0.15516677 + v*(3.28867591 + v*(-2.34357915 + v*0.42139162)))
          endif
        endif

        ex = exp(-v2)
        v2ex = 2.*v2*ex
        h1p = h1 + pqs*ex
        h2p = pqs*h1p - v2ex
        h3p = (pqs*(1. - ex + v2ex) - 2.*v2*h1p)/3. + pqs*h2p
        h4p = (v2ex*v2- pqs*h1p)/3. + pqs*h3p
        !psi = 0.979895023 + a*(-0.962846325 + a*(0.532770573 + a*-0.122727278)) / sqp
        ! 0.2 < a <= 1.4  and  a + v <= 3.2
        h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
      endif

    elseif (v >= 5.0) then
      v4 = 4*v2*v2
      h = ipi * a*(15. + 6.*v2 + v4) / (v4*v2)

    else
      if (v >= 2.4) then
        h1 = 1./(v2 - 1.5) * (0.554153432 + v*(0.278711796 + v*(-0.188325687 + v*(0.042991293 + v*-0.003278278))))
      else
        if (v >= 1.3) then
          h1 = -4.48480194 + v*(9.39456063 + v*(-6.61487486 + v*(1.98919585 + v*-0.22041650)))
        else
          h1 = -1.12470432 + v*(-0.15516677 + v*(3.28867591 + v*(-2.34357915 + v*0.42139162)))
        endif
      endif

      ! a le 0.2  and v lt 5.
      h = (h1*a + exp(-v2)*(1. + a*a*(1. - 2.*v2))) * isqp
    endif
  END FUNCTION voigt_function2
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION voigt_function(a,vs) result(h)
    !---------------------------------------------------------------------
    !  computes a voigt function  h = h(a,v)
    !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
    !  traving (landolt-b\rnstein, p. 449).
    !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
    !---------------------------------------------------------------------
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
    !----------------------------------------------------------------------------------
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
  END FUNCTION voigt_function
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE line_damping (self, adamp, temp, ne, tothi, unsold_temp, stark_pe)
    class(bbtrans_type), intent(in)  :: self
    real(4), intent(in)  :: temp, ne, tothi, unsold_temp, stark_pe
    real(4), intent(inout) :: adamp
    !----------------------------------------------------------------------------------
    ! Van der Waals Broadening
    if(self%unsold) then
      adamp = adamp + self%gw1 * unsold_temp**0.3 * tothi
    elseif (self%gw2 == 0) then
      ! recipe 'W' Turbospectrum
      adamp = adamp + self%gw1 * temp**0.3 * tothi
    else
      adamp = adamp + self%gw1 * unsold_temp**self%gw2 * tothi
    endif

    ! Stark Broadening
    if (self%gq >= 0.0) then
      adamp = adamp + self%gq * ne
    else! if line%gq < 0
      adamp = adamp - self%gq * stark_pe
    end if
    
    adamp = adamp + self%ga
  END SUBROUTINE line_damping
  !------------------------------------------------------------------------------------ 
  
  !====================================================================================
  PURE SUBROUTINE line_dampingH (self, adamp, temp, ne, nh6)
    class(bbtrans_type), intent(in)  :: self
    real(4), intent(in), dimension(:,:,:)  :: temp, ne, nh6
    real(4), intent(inout), dimension(:,:,:) :: adamp
    integer :: i, j
    !----------------------------------------------------------------------------------
    j = nint(sqrt(self%upper_g / 2))
    i = nint(sqrt(1 / (1./i**2 - self%nu0 * H_E / RYD)))
  
    call self%LinearStarkH(adamp,ne,i,j)
    call self%IonicH(adamp,temp,nh6,ne,i,j)
  END SUBROUTINE line_dampingH
  !------------------------------------------------------------------------------------

  !====================================================================================
  ! Linear Stark broadening for hydrogen, see
  ! Sutton K. (1978), JQSRT, 20, 333
  !====================================================================================
  PURE ELEMENTAL SUBROUTINE LinearStarkH(self, gamma, ne, i, j)
    class(bbtrans_type), intent(in)  :: self
    real(4), intent(inout) :: gamma
    real(4), intent(in)    :: ne
    integer, intent(in) :: i, j
    real(4) :: a1
    !----------------------------------------------------------------------------------
    if (j-i == 1) then
      a1=0.642
    else
      a1=1.000
    endif

    gamma = gamma + a1 * 0.6 * (j**2-i**2) * ne**two_third
  END SUBROUTINE LinearStarkH
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE prepare_Stark(self, gq)
    class(bbtrans_type), intent(inout)  :: self
    real(dp), intent(in) :: gq
    !----------------------------------------------------------------------------------
    if (gq < 0) then
      ! make negative to mark this recipe
      self%gq = - 10**(19.4 + 2./3.*self%gq)
    else
      self%gq = gq
    endif
  END SUBROUTINE prepare_Stark
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE IonicH (self, gamma, temp, nh, ne, il, jl) 
    class(bbtrans_type), intent(in)  :: self
    real(4),intent(in)         :: temp,nh,ne
    real(4),intent(inout)      :: gamma
    integer, intent(in) :: il, jl
    integer :: i,j,k
    real(8) :: fn,fn1,rd,rmax,rmin,dohm,gprot,xmr,xmp,xm,rratio,xn
    real(8) :: gion, lam0
    !----------------------------------------------------------------------------------
    lam0 = CC / self%nu0 * 1e8

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
      rmax=min( rd , lam0 / 2. / pi )

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
        rmax=min(rd,lam0/2./pi)
        rmin=sqrt(2./3.)*sqrt(fn)/2./pi/em*hh*sqrt(2*xm*uu/pi/bk/temp) 
        rratio=max(rmax/rmin,one)
        dohm=2.7e8*sqrt(uu*2*xm/temp)*xn*fn*(0.5+log(rratio))
        gion=dohm
      else
        gion=0.0
      endif

      gamma = gamma + gprot + gion
    endif
  END SUBROUTINE IonicH
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE prepare_VanderWaals(self, gv, lower_ion, cnt_ev)
    class(bbtrans_type), intent(inout)  :: self
    integer, intent(in) :: lower_ion
    real(dp), intent(in) :: gv, cnt_ev
    real(dp) :: upper_ev, sigma, alpha, gx, gammaf, c1
    !----------------------------------------------------------------------------------
    ! compute van der waals broadening due to hydrogen
    ! following Barklem
    self%unsold = .false.
    if(gv.ge.20.0) then
      sigma=int(gv)*2.80028e-21
      alpha=gv-int(gv)
      gx=2.0-alpha*0.5
      gx=gx-1.0
      gammaf=1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*gx)*gx)*gx)*gx)*gx

      self%gw2 = 0.5*(1 - alpha)
      self%gw1 = (4.0/pi)**(alpha*0.5) * gammaf * 2.e10 * sigma * 10**(-12*self%gw2)
    
    elseif (gv < 0.0) then
      ! recipe 'W' Turbospectrum
      self%gw1 = 10**gv / 10000.**0.3
      self%gw2 = 0.0

    elseif(gv.ne.0.0) then
      self%unsold = .true.
      upper_ev = self%lower_ev + H_E * self%nu0
      ! Unsold approximation
      ! Jul  5 2010 jorritl: I'm sceptical about the numerical 
      ! constant  1.283984e-12. According to Gray it should be 
      ! 0.3d-30^0.4 = 6.1780059e-13 and  8.411 should be 17
      ! I haven't found the original reference to Unsold 1955 yet...
      c1 = (1./(cnt_ev - upper_ev)**2 - 1./(cnt_ev - self%lower_ev)**2)**.4
      self%gw1 = gv * 8.411 * 1.283984e-12 * lower_ion**0.8 * c1
    endif
  END SUBROUTINE prepare_VanderWaals
  !------------------------------------------------------------------------------------

  !====================================================================================
  FUNCTION estimate_qmax(self, teff) result(qmax)
    class(bbtrans_type), intent(in)  :: self
    real(8), intent(in) :: teff
    real(8) :: blz, fac, qmax
    real(8), parameter :: a=3.0, b=2.8
    !----------------------------------------------------------------------------------
    blz = self%lower_ev * EK / teff
    fac = log(self%Aji) - blz
    !qmax = 0.625 + fac*( 0.73 + fac*(-0.24 + fac*0.03))
    qmax = exp(fac/b)/a + (1 - 1/a)
    if (qmax < 1) then
      qmax = 1
    endif
  END FUNCTION estimate_qmax
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE get_range(self, nu, dlam, ired, iblue)
    class(bbtrans_type), intent(in)  :: self
    real(dp), dimension(:), intent(in) :: nu
    real(dp), intent(in) :: dlam
    integer, intent(out) :: ired, iblue
    real(dp) :: lam0, lam_red, lam_blue, nu_red, nu_blue
    !----------------------------------------------------------------------------------
    lam0 = CC / self%nu0 * 1e8
    lam_red  = lam0 + dlam
    lam_blue = lam0 - dlam
    nu_red  = CC / lam_red  * 1e8
    nu_blue = CC / lam_blue * 1e8

    ired  = searchsorted(nu, nu_red,  outbounds_ok=.true.)
    iblue = searchsorted(nu, nu_blue, outbounds_ok=.true.)
    if (iblue==ired) then
     ired  = 0
     iblue = 0
    elseif (ired == 0) then
      ired = 1
    endif
  END SUBROUTINE get_range
  !------------------------------------------------------------------------------------
END MODULE bbtrans_mod