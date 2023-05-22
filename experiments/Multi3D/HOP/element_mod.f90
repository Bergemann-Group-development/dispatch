MODULE element_mod
    USE prec_mod
    USE op_math_mod
    USE ion_mod
    implicit none

    type isotope_type
        integer  :: iso
        real(dp) :: awgt, f
    contains
        procedure :: compute_doppler_factor
        procedure :: compute_vbar2
    end type isotope_type

    type, extends(isotope_type) :: element_type
        integer :: np, nion, nmolec, mcount, nx, ny, nz
        real(dp) :: abund
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

!====================================================================================
pure elemental function compute_doppler_factor(self, temp, vmic) result(doppler_factor)
  class(isotope_type), intent(in) :: self
  real(4),  intent(in)  :: temp, vmic
  real(4) :: doppler_factor
  !---------------------------------------------------------------------------------
  doppler_factor = sqrt(  2.0d0*bk*temp/(self%awgt*UU) + (vmic*1e5)**2 ) / cc
end function compute_doppler_factor

!====================================================================================
pure elemental function compute_vbar2(self, temp) result(vbar2)
  class(isotope_type), intent(in) :: self
  real(4),  intent(in)  :: temp
  real(4) :: vbar2
  !---------------------------------------------------------------------------------
  vbar2 = 8.0*BK*temp/PI * (1.0/MH + 1.0/(self%awgt*UU))
end function compute_vbar2

END MODULE element_mod