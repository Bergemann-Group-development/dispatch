module element_mod
    use prec_mod
    use op_math_mod
    use ion_mod
    implicit none

    type :: element_type
        integer :: np, nion, nmolec, mcount, nx, ny, nz
        real(dp) :: awgt, abund
        type(ion_type), allocatable :: ions(:)
        !type(mol_pointer_t), allocatable :: molecules(:)
        integer :: active
        real(dp), dimension(:,:,:), allocatable :: fil, n
        character(len=3) :: id
        contains
        procedure :: compute_ion_equilibrium
    end type element_type
    
contains

!===============================================================================
SUBROUTINE compute_ion_equilibrium(self, EK_over_temp, ne, saha, dxi)
    class(element_type), target :: self
    type(ion_type), pointer :: ions(:), ion
    real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, ne, dxi, saha
    real(dp), dimension(self%nx, self%ny, self%nz) :: fil, partition_ratio
    integer :: j
    !---------------------------------------------------------------------------
    ions => self%ions

    do j=1, self%nion
        ion => self%ions(j)
        call ion%allocate(self%nx, self%ny, self%nz)
        ion%part = 1 / ion%partition(EK_over_temp, ne)
    end do

    fil = 1
    do j=self%nion-1, 1, -1
        partition_ratio = ions(j)%part / ions(j+1)%part
        ions(j)%saha = ions(j)%saha_eq(EK_over_temp, ne, saha, dxi, partition_ratio)
        fil = 1 + fil * ions(j)%saha
    end do

    fil = 1 / fil
    ions(1)%anjon = fil

    do j=2, self%nion
        fil = fil*ions(j-1)%saha
        ions(j)%anjon = fil
    end do
END SUBROUTINE compute_ion_equilibrium

end module element_mod