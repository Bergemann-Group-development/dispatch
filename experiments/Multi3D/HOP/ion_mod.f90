module ion_mod
    use prec_mod
    use level_mod
    implicit none

    type particle_type
        real(dp), dimension(:,:,:), allocatable :: f, rho_per_h
    end type

    type ion_pointer_t
        type(ion_type), pointer :: p
    end type ion_pointer_t
    
    type, extends(particle_type) :: ion_type
        character(len=6) :: id
        integer :: g0, nk, j, iqfix
        type(level_type), allocatable :: levels(:)
        type(ion_type), pointer :: above
        type(ion_type), pointer :: below
        integer :: active
        real(dp) :: parco, awgt
        real(dp), pointer :: abund
        real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
        real(dp), dimension(:,:,:), allocatable :: anjon, saha
        real(dp), dimension(:,:,:), allocatable :: pk, part
    contains
        procedure :: partition
        procedure :: saha_eq
        !procedure :: absorption_factor => ion_absorption_factor
    end type

contains


elemental real(dp) function partition(self, EK_over_temp, ne)
    class(ion_type), intent(in) :: self
    real(dp), intent(in) :: EK_over_temp, ne
    real(dp), dimension(self%nk) :: blz, fis, g

    fis = 1

    if (self%iqfix == 0) then
        partition = self%parco
    else
        g   = self%levels%g2
        blz = self%levels%boltzmann(EK_over_temp)
        where(self%levels%xl > 0) fis = self%levels%fischel_sparks(self%j, EK_over_temp, ne)
        partition = self%g0 + sum(g*fis*blz)
    endif
end function partition

elemental real(dp) function saha_eq(self, EK_over_temp, ne, saha, dxi, partition_ratio)
    class(ion_type), intent(in) :: self
    real(dp), intent(in) :: EK_over_temp, ne, saha, dxi
    real(dp), intent(in), optional :: partition_ratio
    real(dp) :: dxi_j

    if (self%j == -2) then
        dxi_j = 2*dxi
    else
        dxi_j = dxi*(self%j*(self%j+1)*0.5)
    endif

    if (.not. present(partition_ratio)) then
        saha_eq = saha * self%levels(1)%boltzmann(EK_over_temp, dxi=dxi_j) * &
            self%above%partition(EK_over_temp, ne) / self%partition(EK_over_temp, ne)
    else
        saha_eq = saha * self%levels(1)%boltzmann(EK_over_temp, dxi=dxi_j) * partition_ratio
    endif
end function saha_eq
    
end module ion_mod