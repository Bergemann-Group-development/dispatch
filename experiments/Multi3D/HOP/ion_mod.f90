module ion_mod
    use prec_mod
    use level_mod
    implicit none

    type particle_type
        real(dp), dimension(:,:,:), allocatable :: f, rho_per_h, e
        character(len=6) :: id
    end type

    type particle_pointer_t
        class(particle_type), pointer :: p
    end type

    type ion_pointer_t
        type(ion_type), pointer :: p
    end type ion_pointer_t
    
    type, extends(particle_type) :: ion_type
        integer :: g0, nk, j, iqfix, nmolec=0, mcount
        type(level_type), allocatable :: levels(:)
        type(ion_type), pointer :: above
        type(ion_type), pointer :: below
        class(particle_pointer_t), allocatable :: molecules(:)
        integer :: active
        real(dp) :: parco, awgt
        real(dp), pointer :: abund
        real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
        real(dp), dimension(:,:,:), allocatable :: anjon, saha
        real(dp), dimension(:,:,:), allocatable :: pk, part, fil
    contains
        procedure :: partition
        procedure :: saha_eq
        procedure :: molecules_with
        procedure :: f_of_own_molecules
        procedure :: f_of_molecules_with
        procedure :: e_of_molecules_with
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

function f_of_own_molecules(self) result(fsum)
    class(ion_type) :: self
    type(particle_type), pointer :: molecule
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
    class(ion_type) :: self
    type(ion_type) :: other
    type(particle_type), pointer :: molecule
    type(particle_pointer_t), allocatable :: mol_list(:)
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
    class(ion_type) :: self
    type(ion_type) :: other
    class(particle_type), pointer :: molecule
    type(particle_pointer_t), allocatable :: mol_list(:)
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
    class(ion_type), target, intent(in) :: self
    type(ion_type), intent(in) :: other
    logical, intent(in), optional :: count_double
    type(particle_pointer_t), intent(out), allocatable :: mol_list(:)
    class(particle_type), pointer :: molecule
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
    
end module ion_mod