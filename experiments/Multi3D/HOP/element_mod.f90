module element_mod
    use prec_mod
    use op_math_mod
    use ion_mod
    implicit none

    type :: element_pointer_t
        type(element_type), pointer :: p
    end type element_pointer_t

    type :: element_type
        character(len=3) :: id
        integer :: np, nion, nmolec, mcount
        real(dp) :: awgt, abund
        type(ion_type), allocatable :: ions(:)
        type(mol_pointer_t), allocatable :: molecules(:)
        integer :: active
        real(dp), dimension(:,:,:), allocatable :: f, fil, n
        contains
        procedure :: compute_ion_equilibrium
        procedure :: f_of_own_molecules
        procedure :: f_of_molecules_with
        procedure :: e_of_molecules_with
        procedure :: molecules_with
    end type element_type


    type mol_pointer_t
        type(molec_type), pointer :: p
    end type mol_pointer_t

    type, extends(particle_type) :: molec_type
        character(len=6) :: id
        integer :: nel
        type(element_pointer_t), allocatable :: elements(:)
        type(ion_type), pointer :: ion => null()
        real(dp) :: dc(5) ! dissociation_constants
        real(dp), dimension(:), allocatable :: de ! dissociation_energies
        real(dp), dimension(:,:,:), allocatable :: e, dis
    contains
        procedure :: pk
        procedure :: calc_f
        procedure :: calc_e
    end type
    
contains

subroutine compute_ion_equilibrium(self, EK_over_temp, ne, saha, dxi)
    class(element_type), target :: self
    type(ion_type), pointer :: ions(:)
    real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, ne, dxi, saha
    real(dp), dimension(:,:,:), allocatable :: fil, partition_ratio
    integer :: j

    allocate(fil, mold=ne)
    
    ions => self%ions

    do j=1, self%nion
        ions(j)%part = ions(j)%partition(EK_over_temp, ne)
    end do

    fil = 1
    do j=self%nion-1, 1, -1
        partition_ratio = ions(j+1)%part / ions(j)%part
        ions(j)%saha = ions(j)%saha_eq(EK_over_temp, ne, saha, dxi, partition_ratio)
        fil = 1 + fil * ions(j)%saha
    end do

    fil = 1 / fil
    ions(1)%anjon = fil

    do j=2, self%nion
        fil = fil*ions(j-1)%saha
        ions(j)%anjon = fil
    end do

end subroutine compute_ion_equilibrium


! ----------------------------------------------------------------------------
! Computes dissociation constants N(AB) / (N(A)*N(B)) for molecules
! ----------------------------------------------------------------------------
elemental real(dp) function pk(self, teta)
    class(molec_type), intent(in) :: self
    real(dp), intent(in) :: teta

    pk = 10**(-poly1d(teta, self%dc(1), self%dc(2), self%dc(3), self%dc(4), self%dc(5)))
end function


function calc_f(self, teta, pe_over_fe) result(f)
    class(molec_type), intent(in) :: self
    real(dp), dimension(:,:,:) :: teta, pe_over_fe
    real(dp), dimension(:,:,:), allocatable :: f
    type(element_type), pointer :: element
    integer :: j

    f = self%dis

    if (associated(self%ion)) then
        f = f * self%ion%saha
    endif

    do j=1, self%nel
        element => self%elements(j)%p

        if (.not. allocated(element%f)) then
            allocate(element%f, mold=f)
            element%f = 1
        endif

        f = f * element%f
        if (j > 1) f = f * pe_over_fe
    end do
end function calc_f


function calc_e(self, teta, temp, dxi) result(e)
    class(molec_type), intent(in) :: self
    type(element_type), pointer :: element
    real(dp), dimension(:,:,:) :: teta, temp, dxi
    real(dp), dimension(:,:,:), allocatable :: dis, e
    integer :: i

    allocate(dis, mold=temp)

    if (size(self%de) > 1) then
        dis = poly1d(teta, self%de(1), self%de(2), self%de(3), self%de(4), self%de(5))*temp/EK
        
        if (trim(self%id) == 'H2' ) dis = - (dis - 4.476)
        if (trim(self%id) == 'H2p') dis = - (dis - 2.648)
    else
        dis = self%de(1)
    endif

    do i=1, self%nel
        element => self%elements(i)%p

        ! why only consider ionisation energy from H?
        if (trim(element%id) == 'H') then
            dis = dis + element%ions(1)%levels(1)%xion - dxi
            if (trim(self%id) == 'H2p') exit ! only one H contribution from H2p
        endif
    end do

    e = -dis * self%f
end function


function f_of_own_molecules(self) result(fsum)
    class(element_type) :: self
    type(molec_type), pointer :: molecule
    integer :: i
    real(dp), dimension(:,:,:), allocatable :: fsum

    allocate(fsum, mold=self%f)

    fsum = 0
    do i=1, self%nmolec
        molecule => self%molecules(i)%p
        fsum = fsum + molecule%f
    end do
end function


function f_of_molecules_with(self, other, count_double) result(fsum)
    class(element_type) :: self
    type(element_type) :: other
    type(molec_type), pointer :: molecule
    type(mol_pointer_t), allocatable :: mol_list(:)
    integer :: i
    real(dp), dimension(:,:,:), allocatable :: fsum
    logical, optional :: count_double

    allocate(fsum, mold=self%f)

    call self%molecules_with(other, mol_list, count_double)

    fsum = 0
    do i=1, size(mol_list)
        molecule => mol_list(i)%p
        fsum = fsum + molecule%f
    end do
end function

function e_of_molecules_with(self, other, count_double) result(esum)
    class(element_type) :: self
    type(element_type) :: other
    type(molec_type), pointer :: molecule
    type(mol_pointer_t), allocatable :: mol_list(:)
    integer :: i
    real(dp), dimension(:,:,:), allocatable :: esum
    logical, optional :: count_double

    allocate(esum, mold=self%f)

    call self%molecules_with(other, mol_list, count_double)

    esum = 0
    do i=1, size(mol_list)
        molecule => mol_list(i)%p
        esum = esum + molecule%e
    end do
end function

subroutine molecules_with (self, other, mol_list, count_double)
    class(element_type), target, intent(in) :: self
    type(element_type), intent(in) :: other
    logical, intent(in), optional :: count_double
    type(mol_pointer_t), intent(out), allocatable :: mol_list(:)
    type(molec_type), pointer :: molecule
    logical :: cd
    integer :: mi, mj, count

    if (present(count_double)) then
        cd = count_double
    else
        cd = .false.
    endif

    count = 0

    do mi=1, self%nmolec
        molecule => self%molecules(mi)%p

        if (mi>1 .and. .not. cd) then
            if(molecule%id == self%molecules(mi-1)%p%id) cycle
        endif

        do mj=1, other%nmolec
            if (molecule%id == other%molecules(mj)%p%id) count = count + 1
        end do
    end do

    allocate(mol_list(count))

    count = 0
    do mi=1, self%nmolec
        molecule => self%molecules(mi)%p

        if (mi>1 .and. .not. cd) then
            if(molecule%id == self%molecules(mi-1)%p%id) cycle
        endif

        do mj=1, other%nmolec
            if (molecule%id == other%molecules(mj)%p%id) then
                count = count + 1
                mol_list(count)%p => molecule
            endif
        end do
    end do

end subroutine molecules_with

end module element_mod