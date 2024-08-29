MODULE bbtrans_mod
  USE prec_mod
  USE const_mod
  USE op_math_mod, only: searchsorted
  USE trace_mod
  implicit none

  type :: bbtrans_type
    real(dp) :: nu0
    integer :: nnu ! number of frequencies in the line
    real(dp) :: ga=0, gw1=0, gw2=0, gw3=0, gq=0, gv=0 ! line broadening (damping) constants
    real(dp) :: lower_ev, upper_g, loggf
    real(dp) :: Aji, Aij, Bji, Bij ! Einstein coefficients, Bij is the
                              ! transition from the lower to the upper
                              ! level
    logical :: active=.true., unsold
    logical, pointer :: used => null()
    integer :: iso_idx=0 ! isotope index for lines in line list
    ! Added to deal with multi component lines
    ! IO > 1 leads to the following being read from the atom file
    integer :: io ! number of components in term (lifted from Multi 2.3)
    real(8), dimension(:), allocatable :: dq
    real(8), dimension(:), allocatable :: wt 
    real(8), dimension(:), allocatable :: da
    integer :: j=0 ! upper level  
    integer :: i=0 ! lower level
  contains
    procedure :: init_einstein_coefficients
    procedure :: damping => line_damping
    procedure :: dampingH => line_dampingH
    procedure :: add_lte_opacities
    procedure :: LinearStarkH
    procedure :: IonicH
    procedure :: prepare_Stark
    procedure :: prepare_VanderWaals
    procedure :: estimate_qmax
    procedure :: get_range
    procedure :: deactivate
  end type

  integer, save :: count1, count2, count3, count4, count5, count6
  logical, public :: disable_stark=.false., disable_vdw=.false., disable_ga=.false.
  real(4), public, save :: tolerance=1e-5
    
contains

  !====================================================================================
  SUBROUTINE init_einstein_coefficients(self, f, gij)
    class(bbtrans_type), intent(inout) :: self
    real(8), intent(in) :: f, gij
    real(8), parameter :: f_to_Aji = 8.0d0 * pi**2 * ec**2 / em / cc**3
    !----------------------------------------------------------------------------------
    ! Einstein coefficients include an additional factor of h_4pi for convenience
    self%Aij = f_to_Aji * f * self%nu0**2
    self%Aij = self%Aij * h_4pi
    self%Bij = c2_2h / self%nu0**3 * self%Aij
    self%Aji = self%Aij * gij
    self%Bji = self%Bij * gij
  END SUBROUTINE init_einstein_coefficients

  !====================================================================================
  SUBROUTINE add_lte_opacities(self, ired, iblue, nu, vproj, inv_df, EK_over_temp, dxij, inv_part, &
    ne, tothi, vbar2, stark_pe, unsold, chi, cont)
    class(bbtrans_type), intent(inout) :: self
    integer, intent(in) :: ired, iblue
    real(8), dimension(:), intent(in), contiguous :: nu
    real(4), dimension(:,:,:), intent(in), contiguous :: vproj, inv_df, inv_part
    real(4), dimension(:,:,:), intent(in), contiguous :: EK_over_temp, dxij, ne, tothi, vbar2, stark_pe, unsold
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    real(4), dimension(:,:,:,:), intent(in) :: cont
    real(8), allocatable :: qq(:)
    real(4) :: ex, ni, nj, nij, upper_ev, nu0_ev, scalar, vs, lchi
    real(4) :: a, idf, cf, v, ref, ekt, Bji
    integer :: i,j,k,f, f0, m(3), itimer=0
    logical :: used
    !----------------------------------------------------------------------------------
    !call trace%begin('add_lte_opacities_1', itimer=itimer)
    m = shape(vproj)
    scalar = over_4pi / self%nu0
    nu0_ev = H_E * self%nu0
    upper_ev = self%lower_ev + nu0_ev
    Bji = self%Bji

    allocate(qq(ired:iblue))
    qq(ired:iblue) = nu(ired:iblue) / self%nu0 - 1
    f0 = minloc(abs(qq),1) + ired - 1
    used = .false.
    
    do k=1, m(3)
    do j=1, m(2)
    do i=1, m(1)
      ! Fixed Broadening
      a = self%ga
      ! Van der Waals Broadening
      if (self%unsold) then
        a = a + self%gw1 * unsold(i,j,k)
      else
        a = a + self%gw1 * vbar2(i,j,k)**self%gw2 * tothi(i,j,k)
      endif
      ! Stark Broadening
      if (self%gq >= 0.0) then
        a = a + self%gq * ne(i,j,k)
      else! if line%gq < 0
        a = a - self%gq * stark_pe(i,j,k)
      end if

      ! fetch variables
      idf = inv_df(i,j,k)
      v   = vproj(i,j,k)
      ekt = EK_over_temp(i,j,k)
      ex  = (upper_ev - dxij(i,j,k)) * ekt

      if (ex < 85.) then
        ! compute level populations
        a = a * scalar * idf
        nij = exp(nu0_ev * ekt)             ! this is ni/nj gj/gi
        nj = inv_part(i,j,k) * exp(-ex)     ! this is nj/gj
        cf = nj * Bji * (nij - 1.)          ! this is (ni*Bij - nj*Bji)*h/4pi
        cf = cf * idf

        ! check if line peak is larger than continuum tolerance
        ref = tolerance * cont(f0,i,j,k)
        !call scater_voigt(f0,ired,iblue,a,qq,v,cf,idf,ref,chi(:,i,j,k))
        
        vs = (qq(f0) - v) * idf
        lchi = voigt_function2(a, vs) * cf
        if (lchi > ref) then
          used = .true.
          chi(f0,i,j,k) = chi(f0,i,j,k) + lchi
          do f=f0+1, iblue
            vs = (qq(f) - v) * idf
            lchi = voigt_function2(a, vs) * cf
            chi(f,i,j,k) = chi(f,i,j,k) + lchi
            if (lchi < ref) exit ! if line opacity below tolerance exit
          end do
          do f=f0-1, ired, -1
            vs = (qq(f) - v) * idf
            lchi = voigt_function2(a, vs) * cf
            chi(f,i,j,k) = chi(f,i,j,k) + lchi
            if (lchi < ref) exit ! if line opacity below tolerance exit
          end do
        endif
      endif
    end do
    end do
    end do
    if (used) self%used = used
    !call trace%end(itimer)
  END SUBROUTINE add_lte_opacities
  !------------------------------------------------------------------------------------

  !====================================================================================
  !DIR$ ATTRIBUTES FORCEINLINE :: voigt_function2
  ELEMENTAL FUNCTION voigt_function2(a,vs) result(h)
    !---------------------------------------------------------------------
    !  computes a voigt function  h = h(a,v)
    !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
    !  traving (landolt-b\rnstein, p. 449).
    !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
    !---------------------------------------------------------------------
    real(4), intent(in) :: a,vs
    real(4) :: h
    real(4) :: v,u,v2,v4,ex
    real(4) :: a2,u2
    complex(8) :: tav,uav,w4,u4 
    real(4), parameter :: pqs  = 1.128379167  ! 2./pi**0.5
    real(4), parameter :: isqp = 0.5641895835 ! 1/pi**0.5
    real(4), parameter :: ipi  = 0.3183098862 ! 1/pi
    !----------------------------------------------------------------------------------
    v = abs(vs)
    u = a + v
    v2 = v*v

    if (a > 0.2) then
      if (a > 1.4 .or. u > 3.2) then
        a2 = a*a
        u = 1 / (a2 + v2)
        u2 = 0.5 * u*u
        ! a gt 1.4  or  a + v gt 3.2
        h = ipi*a*u * (1. + u2*((3.*v2 - a2) + u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2)))

      else
        tav = cmplx(a,-v)
        if (a >= 0.195e0*v - 0.176e0) then
          w4 = (16.4955e0 + tav*(20.20933e0 + tav*(11.96482e0      &
            +  tav*(3.778987e0 + tav*0.5642236e0)))) / (16.4955e0  &
            +  tav*(38.82363e0 + tav*(39.27121e0 + tav*(21.69274e0 &
            +  tav*(6.699398e0 + tav)))))
        else
          uav = tav*tav 
          w4 = tav*(36183.31e0 - uav*(3321.9905e0 - uav*(1540.787e0 &
             - uav*(219.0313e0 - uav*( 35.76683e0 - uav*(1.320522e0 &
             - uav*0.56419e0))))))
          u4 = (32066.6e0 - uav*(24322.84e0 - uav*(9022.228e0       &
             - uav*(2186.181e0 - uav*(364.2191e0                    &
             - uav*(61.57037e0 - uav*(1.841439e0 - uav)))))))
          w4 = exp(uav) - w4/u4
        endif
        h = real(w4) * isqp
      endif

    elseif (v >= 5.0) then
      v4 = 4*v2*v2
      h = ipi * a*(15. + 6.*v2 + v4) / (v4*v2)

    else
      if (v >= 2.4) then
        h = 1./(v2 - 1.5) * (0.554153432 + v*(0.278711796 + v*(-0.188325687 + v*(0.042991293 + v*-0.003278278))))
      elseif (v >= 1.3) then
        h = -4.48480194 + v*(9.39456063 + v*(-6.61487486 + v*(1.98919585 + v*-0.22041650)))
      else
        h = -1.12470432 + v*(-0.15516677 + v*(3.28867591 + v*(-2.34357915 + v*0.42139162)))
      endif

      ! a le 0.2  and v lt 5.
      h = (h*a + exp(-v2)*(1. + a*a*(1. - 2.*v2))) * isqp
    endif
  END FUNCTION voigt_function2
  !------------------------------------------------------------------------------------

  !====================================================================================
  !DIR$ ATTRIBUTES FORCEINLINE :: voigt_function2
  ELEMENTAL FUNCTION voigt_profile(a,q,vpr,idf) result(h)
    !---------------------------------------------------------------------
    !  computes a voigt function  h = h(a,v)
    !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
    !  traving (landolt-b\rnstein, p. 449).
    !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
    !---------------------------------------------------------------------
    real(4), intent(in) :: a,q,vpr,idf
    real(4) :: h
    real(4) :: v,u,v2,v4,ex
    real(4) :: a2,u2
    complex(8) :: tav,uav,w4,u4,num,den
    real(4), parameter :: pqs  = 1.128379167  ! 2./pi**0.5
    real(4), parameter :: isqp = 0.5641895835 ! 1/pi**0.5
    real(4), parameter :: ipi  = 0.3183098862 ! 1/pi
    !----------------------------------------------------------------------------------
    v = abs(q - vpr) * idf
    u = a + v
    v2 = v*v

    if (a > 0.2) then
      if (a > 1.4 .or. u > 3.2) then
        a2 = a*a
        u = 1 / (a2 + v2)
        u2 = 0.5 * u*u
        ! a gt 1.4  or  a + v gt 3.2
        h = ipi*a*u * (1. + u2*((3.*v2 - a2) + u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2)))

      else
        tav = cmplx(a,-v)
        if (a >= 0.195e0*v - 0.176e0) then
          num = 16.4955e0 + tav*(20.20933e0 + tav*(11.96482e0      &
                          + tav*(3.778987e0 + tav*0.5642236e0)))
          den = 16.4955e0 + tav*(38.82363e0 + tav*(39.27121e0 &
                          + tav*(21.69274e0 + tav*(6.699398e0 + tav))))
          w4 = num/den
        else
          uav = tav*tav 
          num = tav*(36183.31e0 - uav*(3321.9905e0 - uav*(1540.787e0 &
             - uav*(219.0313e0 - uav*( 35.76683e0 - uav*(1.320522e0 &
             - uav*0.56419e0))))))
          den = (32066.6e0 - uav*(24322.84e0 - uav*(9022.228e0       &
             - uav*(2186.181e0 - uav*(364.2191e0                    &
             - uav*(61.57037e0 - uav*(1.841439e0 - uav)))))))
          w4 = num/den
          w4 = exp(uav) - w4
        endif
        h = real(w4) * isqp
      endif

    elseif (v >= 5.0) then
      v4 = 4*v2*v2
      h = ipi * a*(15. + 6.*v2 + v4) / (v4*v2)

    else
      if (v >= 2.4) then
        h = 1./(v2 - 1.5) * (0.554153432 + v*(0.278711796 + v*(-0.188325687 + v*(0.042991293 + v*-0.003278278))))
      elseif (v >= 1.3) then
        h = -4.48480194 + v*(9.39456063 + v*(-6.61487486 + v*(1.98919585 + v*-0.22041650)))
      else
        h = -1.12470432 + v*(-0.15516677 + v*(3.28867591 + v*(-2.34357915 + v*0.42139162)))
      endif

      ! a le 0.2  and v lt 5.
      ex = exp(-v2)
      h = (h*a + ex*(1. + a*a*(1. - 2.*v2))) * isqp
    endif
    h = h * idf * (q + 1)
  END FUNCTION voigt_profile
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
  ELEMENTAL FUNCTION new_voigt(a,v,f) result(h)
    ! Champion Voigt routine in speed and accuracy test by BE+NP 97-11-18
    !  VOIGT FUNCTION CALCULATION: HUMLICEK'S APPROXIMATION 
    ! The results are normalized to an area of 1/sqrt(pi).
    real(4), intent(in) :: f
    real(4), intent(in) :: a, v
    complex(8) :: tav,uav,w4,v4 
    real(4)    :: sav,aa,vv,h
    real(4), parameter :: isqp = 0.5641895835 ! 1/pi**0.5
    integer :: i,j,k
    !----------------------------------------------------------------------------------
    aa = a
    vv = abs(v)
    tav = cmplx(aa,-vv) 
    sav = vv + aa 
    uav = tav*tav 
    if (sav >= 15.e0) then
      w4 = tav*0.5641896e0 / (0.5e0 + uav)
    elseif (sav >= 5.5e0) then
      w4 = tav*(1.410474e0 + uav*0.5641896e0) / (0.75e0 + uav*(3.e0 + uav))
    elseif (aa >= 0.195e0*vv - 0.176e0) then
      w4 = (16.4955e0 + tav*(20.20933e0 + tav*(11.96482e0       &
        +  tav*(3.778987e0 + tav*0.5642236e0)))) / (16.4955e0  &
        +  tav*(38.82363e0 + tav*(39.27121e0 + tav*(21.69274e0 &
        +  tav*(6.699398e0 + tav)))))
    else
      w4 = tav*(36183.31e0 - uav*(3321.9905e0 - uav*(1540.787e0 &
        - uav*(219.0313e0 - uav*( 35.76683e0 - uav*(1.320522e0 &
        - uav*0.56419e0))))))
      v4 = (32066.6e0 - uav*(24322.84e0 - uav*(9022.228e0       &
        - uav*(2186.181e0 - uav*(364.2191e0                    &
        - uav*(61.57037e0 - uav*(1.841439e0 - uav)))))))
      w4 = exp(uav) - w4/v4
    endif
    h = real(w4) * f ! * isqp
  END FUNCTION new_voigt
  !------------------------------------------------------------------------------------

  !====================================================================================
  ELEMENTAL FUNCTION dirty_voigt(a, x) result(h)
    real(4), intent(in) :: a, x
    real(4) :: h, z
    ! Coefficients for the approximation function (Table A1, Smith+15)
    real(4), parameter :: A0 = 15.75328153963877
    real(4), parameter :: A1 = 286.9341762324778
    real(4), parameter :: A2 = 19.05706700907019
    real(4), parameter :: A3 = 28.22644017233441
    real(4), parameter :: A4 = 9.526399802414186
    real(4), parameter :: A5 = 35.29217026286130
    real(4), parameter :: A6 = 0.8681020834678775
    real(4), parameter :: B0 = 0.0003300469163682737
    real(4), parameter :: B1 = 0.5403095364583999
    real(4), parameter :: B2 = 2.676724102580895
    real(4), parameter :: B3 = 12.82026082606220
    real(4), parameter :: B4 = 3.21166435627278
    real(4), parameter :: B5 = 32.032981933420
    real(4), parameter :: B6 = 9.0328158696
    real(4), parameter :: B7 = 23.7489999060
    real(4), parameter :: B8 = 1.82106170570
    real(4), parameter :: isqp = 0.5641895835 ! 1/pi**0.5
    !----------------------------------------------------------------------------------
    z = x * x
    if (z <= 3.0) then
        h = exp(-z) * (1.0 - a * (A0 + A1 / (z - A2 + A3 / (z - A4 + A5 / (z - A6)))))
    else if (z < 25.0) then
        h = exp(-z) + a * (B0 + B1 / (z - B2 + B3 / (z + B4 + B5 / (z - B6 + B7 / (z - B8)))))
    else
        h = 0.5 * sqrt(2.0 / pi) * a / (z - 1.5 - 1.5 / (z - 3.5 - 5.0 / (z - 5.5)))
    end if
    h = h * isqp
  END FUNCTION dirty_voigt

  !====================================================================================
  SUBROUTINE scater_voigt(f0,ired,iblue,a,qq,vel,cf,idf,ref,chi)
    integer, intent(in) :: f0, ired, iblue
    real(4), intent(in) :: a, vel, cf, idf, ref
    real(8), intent(in), dimension(:), contiguous :: qq
    real(4), intent(inout), dimension(:), contiguous :: chi
    ! this vectorizable voigt function is based on the paper by
    ! hui, armstrong and wray, jqsrt 19, 509 (1977). it has been
    ! checked agains the old landolt & boernstein standard voigt. errors
    ! become significant (at the < 1 % level) around the "knee" between
    ! the dopler core and the damping wings for a smaller than 1.e-3.
    ! note that, as written here, the routine has no provision for
    ! the trivial case a=0. the normalization is such that the integral
    ! is sqrt(pi); i.e., a plain exp(-x**2) should be used at a=0.
    ! the new landolt & boernstein gives a useful expression for the
    ! small but finite a case.

    ! coded by: a. nordlund/16-dec-82.
    
    ! WH: reprogrammed the algorithm:
    ! - expanding the complex expression allows the loop to vectorize
    ! - interleaving the evaluation of real and imaginary parts, as well
    !   as the numerator and denominator reduces dependencies and thus
    !   improves pipelining
    !
    ! The original expression was:
    !
    !    z=cmplx(a(i),-abs(v(i)))
    !    h(i)=real(((((((   a6 *z+a5)*z+a4)*z+a3)*z+a2)*z+a1)*z+a0)/
    !            & (((((((z+b6)*z+b5)*z+b4)*z+b3)*z+b2)*z+b1)*z+b0))
    !

    ! Single precision coefficients for the Hui et al. recipe
    real(4), parameter :: a0 = 122.6079
    real(4), parameter :: a1 = 214.3824
    real(4), parameter :: a2 = 181.9285
    real(4), parameter :: a3 = 93.15558
    real(4), parameter :: a4 = 30.18014
    real(4), parameter :: a5 = 5.912626
    real(4), parameter :: a6 = 0.5641896
    real(4), parameter :: b0 = 122.6079
    real(4), parameter :: b1 = 352.7306
    real(4), parameter :: b2 = 457.3345
    real(4), parameter :: b3 = 348.7039
    real(4), parameter :: b4 = 170.3540
    real(4), parameter :: b5 = 53.99291
    real(4), parameter :: b6 = 10.47986
    real(4), parameter :: isqp = 0.5641895835 !1/sqrt(pi)
    ! Compute denominator and numerator of the complex
    ! expression in an alternating manner to improve the
    ! pipelining:
    ! c-prefixes: numerator
    ! d-prefixes: denominator
    real(4) :: vabs,zre,zim,cre1,cim1,cre2,cim2,dre1,dim1,dre2,dim2
    real(4) :: fac, lchi
    integer :: f, f1
    !----------------------------------------------------------------------------------
    fac = cf * isqp
    lchi = voigt(a, qq(f0))
    if (lchi > ref) then
      ! if so, add and continue with blue wing
      chi(f0) = chi(f0) + lchi
      f1 = ired
      do f=f0+1, iblue
        lchi = voigt(a, qq(f))
        ! if line opacity below tolerance exit and save f
        if (lchi < ref) then
          f1 = max(ired, 2*f0-f)
          exit
        endif
        chi(f) = chi(f) + lchi
      end do
      ! now do red wing, starting where blue wing left off
      do f=f1, f0-1
        lchi = voigt(a, qq(f))
        chi(f) = chi(f) + lchi
      end do
    endif

  contains
    function voigt(a,q) result(h)
      real(4), intent(in) :: a
      real(8), intent(in) :: q
      real(4) :: v, h
      !----------------------------------------------------------------------------------
      v = (q - vel) * idf
      vabs=abs(v)
      zre=a
      zim=-vabs

      dre2=zre+b6

      cim1=a6  *zim
      dim1=(zre+dre2)*zim

      cre1=a6  *zre        +a5
      dre1=dre2*zre-zim*zim+b5

      cim2=cim1*zre+cre1*zim
      dim2=dim1*zre+dre1*zim

      cre2=cre1*zre-cim1*zim+a4
      dre2=dre1*zre-dim1*zim+b4

      cim1=cim2*zre+cre2*zim
      dim1=dim2*zre+dre2*zim

      cre1=cre2*zre-cim2*zim+a3
      dre1=dre2*zre-dim2*zim+b3

      cim2=cim1*zre+cre1*zim
      dim2=dim1*zre+dre1*zim

      cre2=cre1*zre-cim1*zim+a2
      dre2=dre1*zre-dim1*zim+b2

      cim1=cim2*zre+cre2*zim
      dim1=dim2*zre+dre2*zim

      cre1=cre2*zre-cim2*zim+a1
      dre1=dre2*zre-dim2*zim+b1

      cim2=cim1*zre+cre1*zim
      dim2=dim1*zre+dre1*zim

      cre2=cre1*zre-cim1*zim+a0
      dre2=dre1*zre-dim1*zim+b0

      ! Gather numbers
      cim2=cim2*dim2
      dim2=dim2*dim2

      cre2=cre2*dre2
      dre2=dre2*dre2

      ! Build the denominator and numerator
      dre2=dre2+dim2
      cre2=cre2+cim2

      ! Divide
      h=cre2/dre2
      h = h * fac
    end function voigt
  END SUBROUTINE scater_voigt

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE line_damping (self, adamp, ne, tothi, vbar2, stark_pe, unsold)
    class(bbtrans_type), intent(in)  :: self
    real(4), intent(in)  :: ne, tothi, vbar2, stark_pe, unsold
    real(4), intent(inout) :: adamp
    !----------------------------------------------------------------------------------
    ! Van der Waals Broadening
    if (self%unsold) then
      adamp = adamp + self%gw1 * unsold
    else
      adamp = adamp + self%gw1 * vbar2**self%gw2 * tothi
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
    i = nint(sqrt(1 / (1./j**2 + self%nu0 * H_E / RYD)))
  
    if (.not. disable_stark) call self%LinearStarkH(adamp,ne,i,j)
    if (.not. disable_vdw) call self%IonicH(adamp,temp,nh6,ne,i,j)
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
      ! formula from Gray
      !self%gq = - BK * 10**(19.4 + 2./3.*gq)
      self%gq = - 10**gq * 1e-4**one_sixth
    else
      self%gq = gq
    endif
    if (disable_ga) self%ga = 0
    if (disable_stark) self%gq = 0
  END SUBROUTINE prepare_Stark
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE ELEMENTAL SUBROUTINE IonicH (self, gamma, temp, nh, ne, i, j) 
    class(bbtrans_type), intent(in)  :: self
    real(4),intent(in)         :: temp,nh,ne
    real(4),intent(inout)      :: gamma
    integer, intent(in) :: i, j
    real(8) :: fn,fn1,rd,rmax,rmin,dohm,gprot,xmr,xmp,xm,rratio,xn
    real(8) :: gion, lam0
    !----------------------------------------------------------------------------------
    lam0 = CC / self%nu0 * 1e8

    ! ionic broadening
    if( (j-i >= 1) .and. (j-i <= 3) .and. (i >= 4) .and. (i <= 20) ) then

      fn1 = 2.25 * i**2 * ( 1.5 + 0.271*log(2.0/3.0*i)+  1.92/i**0.422 )
      if(j-i.eq.1) then
        fn=fn1
      else if(j-i.eq.2) then
        fn=fn1*(9.1333/i+5.57333)
      else if(j-i.eq.3) then
        fn=fn1*(38.6667/i+12.9667)
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
  PURE SUBROUTINE prepare_VanderWaals(self, gv, lower_ion, cnt_ev, vbar2e4)
    class(bbtrans_type), intent(inout)  :: self
    integer, intent(in) :: lower_ion
    real(dp), intent(in) :: gv, cnt_ev, vbar2e4
    real(dp) :: upper_ev, sigma, alpha, gx, gammaf, balfac, c1, c625, vbar2
    !----------------------------------------------------------------------------------
    ! compute van der waals broadening due to hydrogen
    ! following Barklem
    self%unsold = .false.
    if(gv >= 20.0) then
      sigma = 2.80028e-21*int(gv)
      alpha = gv - int(gv)
      gx = 1.0 - 0.5*alpha
      gammaf = 1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*gx)*gx)*gx)*gx)*gx
      balfac = (4.0/pi)**(0.5*alpha) * gammaf * 2.e10

      self%gw2 = 0.5*(1 - alpha)
      self%gw1 = balfac * sigma * 10**(-12*self%gw2)
    
    elseif (gv < 0.0) then
      ! recipe 'W' Turbospectrum
      self%unsold = .true.
      self%gw1 = 10**gv / (vbar2e4)**0.3
      self%gw2 = 0.0

    elseif(gv.ne.0.0) then
      self%unsold = .true.
      upper_ev = self%lower_ev + H_E * self%nu0
      ! Unsold approximation
      ! Jul  5 2010 jorritl: I'm sceptical about the numerical 
      ! constant  1.283984e-12. According to Gray it should be 
      ! 0.3d-30^0.4 = 6.1780059e-13 and  8.411 should be 17
      ! I haven't found the original reference to Unsold 1955 yet...
      if (self%lower_ev < cnt_ev) then
        c1 = 1./(cnt_ev - upper_ev)**2 - 1./(cnt_ev - self%lower_ev)**2
        if (c1 > 0) then
          c625 = (lower_ion**2 * c1)**0.4 ! c6^(2/5) impact approximation
          self%gw1 = gv * 8.411 * 1.283984e-12 * c625
        else
          self%gw1 = 0
        endif
      endif
    endif
    if (disable_vdw) then
      self%unsold = .true.
      self%gw1 = 0
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
    if (dlam == 0) then
      ired = findloc(nu, self%nu0, 1)
      iblue = ired
    else
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
    endif
  END SUBROUTINE get_range
  !------------------------------------------------------------------------------------

  !====================================================================================
  PURE SUBROUTINE deactivate(self)
    class(bbtrans_type), intent(inout) :: self
    !----------------------------------------------------------------------------------
    self%active = .false.
  END SUBROUTINE deactivate
END MODULE bbtrans_mod