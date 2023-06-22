MODULE level_mod
    USE prec_mod
    implicit none

    type :: level_type
        real(dp) :: xion, xl, g2=0
    contains
        procedure :: boltzmann
        procedure :: fischel_sparks
    end type
    
contains

    elemental real(dp) FUNCTION boltzmann(self, EK_over_temp, dxi, natural_log)
        class(level_type), intent(in) :: self
        real(dp), intent(in) :: EK_over_temp
        real(dp), intent(in), optional :: dxi
        logical, intent(in), optional :: natural_log
        real(dp) :: dum

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

    elemental real(dp) FUNCTION fischel_sparks(self, ion, EK_over_temp, ne)
        class(level_type), intent(in) :: self
        real(dp), intent(in) :: EK_over_temp, ne
        integer, intent(in) :: ion
        real(dp) :: a, a2, r, r2, r3, axl2, p, p2, p3, p4
        real(dp) :: bk_err

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
    
END MODULE level_mod