MODULE element_mod
    USE prec_mod
    USE op_math_mod
    USE ion_mod
    implicit none
    type isotope_type
        integer  :: iso
        real(dp) :: mass, f
    end type isotope_type

    type :: element_type
        integer :: np, nion, nmolec, mcount, nx, ny, nz
        real(dp) :: awgt, abund
        type(ion_type), allocatable :: ions(:)
        class(isotope_type), allocatable :: isotope(:)
        !type(mol_pointer_t), allocatable :: molecules(:)
        integer :: active
        real(dp), dimension(:,:,:), allocatable :: fil, n
        character(len=3) :: id
        contains
        procedure :: compute_ion_equilibrium
        procedure :: print_info
    end type element_type
    
contains

!===============================================================================
SUBROUTINE compute_ion_equilibrium(self, EK_over_temp, ne, saha, dxi)
    class(element_type), target :: self
    type(ion_type), pointer :: ions(:), ion
    real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, ne, dxi, saha
    real(dp), dimension(self%nx, self%ny, self%nz) :: fil, partition_ratio, partp
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
        partp = ions(j+1)%part
        ! --- this only shows up in the Oslo legacy code ---!
        !partp = max(1d-3, partp)
        ! --------------------------------------------------!
        partition_ratio = ions(j)%part / partp
        ions(j)%saha = ions(j)%saha_eq(EK_over_temp, ne, saha, dxi, partition_ratio)
        fil = 1 + fil * ions(j)%saha
    end do

    ions(1)%f = self%abund / fil

    do j=2, self%nion
        ions(j)%f = ions(j-1)%f * ions(j-1)%saha
    end do
END SUBROUTINE compute_ion_equilibrium

!===============================================================================
SUBROUTINE print_info(self, lunit, p)
    class(element_type), target :: self
    integer, intent(in) :: lunit, p(:)
    type(ion_type), pointer :: ion
    character(32) :: info
    character(32) :: f_string
    character(48) :: saha_string
    character(64) :: part_string
    integer :: j
    !---------------------------------------------------------------------------
    write(info,'(A4, e12.4)') self%id, self%abund
    f_string=''
    part_string=''
    saha_string=''
    do j=1, self%nion
        ion => self%ions(j)
        write(f_string, '(A, f8.4)') trim(f_string), ion%f(p(1),p(2),p(3))
    end do

    do j=1, self%nion
        ion => self%ions(j)
        write(saha_string, '(A, es10.2)') trim(saha_string), ion%saha(p(1),p(2),p(3))
    end do

    do j=1, self%nion
        ion => self%ions(j)
        write(part_string, '(A, es10.2)') trim(part_string), 1/ion%part(p(1),p(2),p(3))
    end do
    write(lunit, '(A,A,A,A)'), trim(info), f_string, saha_string, part_string
END SUBROUTINE print_info

END MODULE element_mod