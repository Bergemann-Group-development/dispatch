MODULE level_mod
    USE prec_mod
    implicit none

    type :: state_type
        real(dp) :: xion
        real(dp), dimension(:,:,:), allocatable :: n
    contains
        procedure :: boltzmann
    end type

    type, extends(state_type) :: level_type
        integer :: i
        real(dp) :: xl, g=0
        real(dp), dimension(:,:,:), allocatable :: w
    contains
        procedure :: fischel_sparks
    end type

    public :: wcalc
    
contains

    !===============================================================================
    elemental real(dp) FUNCTION boltzmann(self, EK_over_temp, dxi, natural_log)
        class(state_type), intent(in) :: self
        real(dp), intent(in) :: EK_over_temp
        real(dp), intent(in), optional :: dxi
        logical, intent(in), optional :: natural_log
        real(dp) :: dum
        !---------------------------------------------------------------------------
        if (present(dxi)) then
            dum = -EK_over_temp*(self%xion - dxi)
        else
            dum = -EK_over_temp*self%xion
        endif 

        if (dum < -100.d0) then
            boltzmann = 0
        else
            if (.not. present(natural_log)) then
                boltzmann = exp(dum * prec_err) ! keep error for consistency here
            elseif(natural_log) then
                boltzmann = exp(dum)
            else
                boltzmann = exp(dum * prec_err)
            endif
        endif
    end FUNCTION boltzmann

    !===============================================================================
    elemental real(dp) FUNCTION fischel_sparks(self, ion, EK_over_temp, ne)
        class(level_type), intent(in) :: self
        real(dp), intent(in) :: EK_over_temp, ne
        integer, intent(in) :: ion
        real(dp) :: a, a2, r, r2, r3, axl2, p, p2, p3, p4
        real(dp) :: bk_err
        !---------------------------------------------------------------------------
        bk_err = BK / uop_BK ! keep around for consistency

        a = 31.321*ion**2 * log10e * EK_over_temp * prec_err
        p = ion * 4.2E3 / (ne*bk_err)**0.166666667
    
        if (p <= self%xl) then
            axl2 = a / self%xl**2
            ! FISCHEL AND SPARKS, EQ. (27)
            fischel_sparks = p**4 / self%xl * (1 + axl2*(0.33333333 + 0.1*axl2))
        else
            r = self%xl - 1
            r2 = r*r
            r3 = r2*r
            p2 = p*p
            p3 = p2*p
            p4 = p3*p
            a2 = a*a
        
            ! FISCHEL AND SPARKS, EQ. (26)
            fischel_sparks = 1.33333333*p3 + 0.5*p2 + 0.16666667*p + 1.33333333*a*p - 0.4*a2/p &
                            - 0.33333333*r3 - 0.5*r2 - 0.16666667*r - a*r + 0.5*a2/self%xl
        endif
    end FUNCTION fischel_sparks

    !===============================================================================
    ! Computes the occupation probability for a H atom in 
    ! state with effective principle quantum number NS in a plasma
    ! enviroment with NH, NE, NHE number densities of H, ions,
    ! and He atoms respectively.  This code assumes the perturbing
    ! neutral H and He are in the ground state, (noting the hard
    ! sphere approximation used is quite crude anyway) and that ions 
    ! are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
    ! See eqn 4.71 Hummer & Milhalas (1988) 
    ! Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
    ! which is a fit to detailed calculation including correlations,
    ! thus the temperature dependence
    ! Sizes of neutral atoms adopted are sqrt(<r^2>)
    ! Coded by Paul Barklem and Kjell Eriksson, Aug 2003
    ! Optimized by Richard Hoppe, Jul 2023
    !===============================================================================
    PURE SUBROUTINE wcalc(nh, nhe, ne, shield, ns, w)
        real(dp), intent(in) , dimension(:,:,:) :: nh, nhe, ne, shield
        real(dp), intent(in) :: ns
        real(dp), intent(out), dimension(:,:,:) :: w
        real(dp) :: ns2, rih, x1, x2, neutr, k1, k2, x, f, charg, res
        real(dp), parameter :: a0=5.29177d-9
        real(dp), parameter :: pi43=4.18879d0  ! 4.18879 is 4*pi/3
        integer :: i,j,k, m(3)
        !---------------------------------------------------------------------------
        ns2 = ns*ns
        rih = ns * sqrt(2.5*ns2 + 0.5)
        x1 = (a0*(rih + 1.73))**3
        x2 = (a0*(rih + 1.02))**3

        k1 = 8.3e14 / (ns2*ns2)
        if(ns > 3.) then
          k1 = k1 * 5.33333333*ns / (ns + 1.d0)**2
        endif
        k1 = k1*sqrt(k1)
        k2 = k1*k1

        m = shape(ne)
        do k=1, m(3)
        do j=1, m(2)
        do i=1, m(1)
          ! neutral perturbers
          neutr = pi43 * (nh(i,j,k)*x1 + nhe(i,j,k)*x2)
          if (neutr > 1e-5) then
            res = exp(-neutr)
          else
            res = 1.0
          endif
          
          ! charged perturbers
          x = shield(i,j,k)
          if (k1*x < 1e4) then
            f = 0.1402d0*x*k2 / (ne(i,j,k) * (1. + 0.1285*x*k1))
            charg = f / (1.d0 + f)
            res = res * charg
          endif

          w(i,j,k) = res
        end do
        end do
        end do
    END SUBROUTINE wcalc
    
END MODULE level_mod