module molecule_mod
    use prec_mod
    use op_math_mod
    use ion_mod
    implicit none

    type, extends(particle_type) :: molec_type
        integer :: nions
        class(ion_pointer_t), allocatable :: ions(:)
        real(dp) :: dc(5) ! dissociation_constants
        real(dp), dimension(:), allocatable :: de ! dissociation_energies
        real(dp), dimension(:,:,:), allocatable :: dis
    contains
        procedure :: pk
        procedure :: calc_f
        procedure :: calc_e
        procedure :: allocate => allocate_molec
    end type
    
contains

! ----------------------------------------------------------------------------
! Computes dissociation constants N(AB) / (N(A)*N(B)) for molecules
! ----------------------------------------------------------------------------
elemental real(dp) FUNCTION pk(self, teta)
    class(molec_type), intent(in) :: self
    real(dp), intent(in) :: teta

    pk = 10**(-poly1d(teta, self%dc(1), self%dc(2), self%dc(3), self%dc(4), self%dc(5)))
end FUNCTION

!===============================================================================
FUNCTION calc_f(self, teta, pe_over_fe) result(f)
    class(molec_type), intent(in) :: self
    real(dp), dimension(:,:,:) :: teta, pe_over_fe
    real(dp), dimension(self%nx,self%ny,self%nz) :: f
    type(ion_type), pointer :: ion
    integer :: j
    !---------------------------------------------------------------------------
    f = self%dis

    do j=1, self%nions
        ion => self%ions(j)%p
        f = f * ion%f
        if (j > 1) f = f * pe_over_fe
    end do
end FUNCTION calc_f

!===============================================================================
FUNCTION calc_e(self, teta, temp, dxi) result(e)
    class(molec_type), intent(in) :: self
    class(ion_type), pointer :: ion
    real(dp), dimension(:,:,:) :: teta, temp, dxi
    real(dp), dimension(self%nx,self%ny,self%nz) :: dis, e
    integer :: i
    !---------------------------------------------------------------------------
    if (size(self%de) > 1) then
        dis = poly1d(teta, self%de(1), self%de(2), self%de(3), self%de(4), self%de(5))*temp/EK
        
        if (trim(self%id) == 'H2' ) dis = - (dis - 4.476)
        if (trim(self%id) == 'H2p') dis = - (dis - 2.648)
    else
        dis = self%de(1)
    endif

    do i=1, self%nions
        ion => self%ions(i)%p

        ! why only consider ionisation energy from H?
        if (trim(ion%id) == 'H  I') then
            dis = dis + ion%levels(1)%xion - dxi
        endif
    end do

    e = -dis * self%f
end FUNCTION

!===============================================================================
SUBROUTINE allocate_molec (self, nx, ny, nz)
    class(molec_type) :: self
    integer :: nx, ny, nz
    !---------------------------------------------------------------------------
    if (.not. self%allocated) then
        self%nx = nx
        self%ny = ny
        self%nz = nz
        allocate(self%f(nx,ny,nz))
        allocate(self%e(nx,ny,nz))
        allocate(self%dis(nx,ny,nz))
        allocate(self%rho_per_h(nx,ny,nz))
        self%allocated = .true.
    endif
END SUBROUTINE allocate_molec

end module molecule_mod