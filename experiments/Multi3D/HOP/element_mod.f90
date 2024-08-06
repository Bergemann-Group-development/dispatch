MODULE element_mod
    USE prec_mod
    USE op_math_mod
    USE ion_mod
    implicit none

    type, extends(massive_type) :: isotope_type
        integer  :: iso
        real(dp) :: frac
    end type isotope_type

    type, extends(isotope_type) :: element_type
        integer :: np, nion, nmolec, mcount
        real(dp) :: abund=-100
        type(ion_type), allocatable :: ions(:)
        class(isotope_type), allocatable :: isotope(:)
        !type(mol_pointer_t), allocatable :: molecules(:)
        integer :: active
        character(len=3) :: id
    contains
        procedure :: compute_partition_functions
        procedure :: compute_ion_equilibrium
        procedure :: print_info
    end type element_type
    
contains

!===============================================================================
SUBROUTINE compute_partition_functions(self, temp, EK_over_temp, ne)
    class(element_type), target :: self
    type(ion_type), pointer :: ion
    type(level_type), pointer :: lvl
    real(dp), dimension(:,:,:), intent(in) :: temp, EK_over_temp, ne
    real(dp), dimension(:,:,:), allocatable :: blz
    integer :: i, j, l, u, x, y, z
    !---------------------------------------------------------------------------
    l = lbound(self%ions,1)
    u = ubound(self%ions,1)

    allocate(blz, mold=temp)

    do j=l, u
        ion => self%ions(j)
        call ion%allocate(shape(ne))

        if (ion%fixed_part) then
            ion%part = ion%parco
        elseif (allocated(ion%irwin)) then ! IRWIN FORMAT
            ion%part = ion%irwin_partf(temp)
        else
            ion%part = 0
            do i=1, ubound(ion%levels, 1)
                lvl => ion%levels(i)
                blz = lvl%g * lvl%boltzmann(EK_over_temp)
                if (lvl%xl > 0) then
                    blz = blz * lvl%fischel_sparks(ion%j, EK_over_temp, ne)
                endif
                ion%part = ion%part + blz
            end do
        endif

        ion%part = 1 / ion%part
    end do
    deallocate(blz)
END SUBROUTINE compute_partition_functions

!===============================================================================
SUBROUTINE compute_ion_equilibrium(self, EK_over_temp, saha, dxi)
    class(element_type), target :: self
    type(ion_type), pointer :: ions(:), ion
    real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, dxi, saha
    real(dp), dimension(:,:,:), allocatable :: fil, partition_ratio, partp
    integer :: j, l, u, x, y, z
    !---------------------------------------------------------------------------
    allocate(fil, partition_ratio, partp, mold=dxi)
    ions => self%ions
    l = lbound(ions,1)
    u = ubound(ions,1)

    fil = 1
    do j=u-1, l, -1
        partp = ions(j+1)%part
        ! --- this only shows up in the Oslo legacy code ---!
        !partp = max(1d-3, partp)
        ! --------------------------------------------------!
        partition_ratio = ions(j)%part / partp
        ions(j)%saha = ions(j)%saha_eq(EK_over_temp, saha, dxi, partition_ratio)
        fil = 1 + fil * ions(j)%saha
    end do

    ions(l)%f = self%abund / fil

    do j=l+1, u
        ions(j)%f = ions(j-1)%f * ions(j-1)%saha
    end do
    deallocate(fil, partition_ratio, partp)
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