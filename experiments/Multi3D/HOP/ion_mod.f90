module ion_mod
    use prec_mod
    use level_mod
    use bbtrans_mod, only: bbtrans_type
    implicit none

    type particle_type
        real(dp), dimension(:,:,:), allocatable :: f, rho_per_h, e, part
        integer :: nx, ny, nz, nline=0
        character(len=6) :: id
        logical :: allocated=.false.
        class(bbtrans_type), dimension(:), allocatable :: line
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
        real(dp) :: parco, awgt, iso_awgt
        real(dp), pointer :: abund
        real(dp) :: acolc, tcolc, aradc, xradc, adic, bdic, t0c, t1c
        real(dp), dimension(:,:,:), allocatable :: anjon, saha
        real(dp), dimension(:,:,:), allocatable :: pk, fil
    contains
        procedure :: partition
        procedure :: saha_eq
        procedure :: molecules_with
        procedure :: f_of_own_molecules
        procedure :: f_of_molecules_with
        procedure :: e_of_molecules_with
        procedure :: allocate => allocate_ion
        !procedure :: absorption_factor => ion_absorption_factor
    end type

contains

!===============================================================================
elemental real(dp) FUNCTION partition(self, EK_over_temp, ne)
    class(ion_type), intent(in) :: self
    real(dp), intent(in) :: EK_over_temp, ne
    real(dp), dimension(self%nk) :: blz, fis, g
    !---------------------------------------------------------------------------
    fis = 1

    if (self%iqfix == 0) then
        partition = self%parco
    else
        g   = self%levels%g2
        blz = self%levels%boltzmann(EK_over_temp)
        where(self%levels%xl > 0) fis = self%levels%fischel_sparks(self%j, EK_over_temp, ne)
        partition = self%g0 + sum(g*fis*blz)
    endif
end FUNCTION partition

!===============================================================================
elemental real(dp) FUNCTION saha_eq(self, EK_over_temp, ne, saha, dxi, partition_ratio)
    class(ion_type), intent(in) :: self
    real(dp), intent(in) :: EK_over_temp, ne, saha, dxi
    real(dp), intent(in), optional :: partition_ratio
    real(dp) :: dxi_j, partp
    !---------------------------------------------------------------------------
    if (self%j == -2) then
        dxi_j = 2*dxi
    else
        dxi_j = dxi*self%j
    endif

    saha_eq = saha * self%levels(1)%boltzmann(EK_over_temp, dxi=dxi_j) * partition_ratio
end FUNCTION saha_eq

!===============================================================================
FUNCTION f_of_own_molecules(self) result(fsum)
    class(ion_type) :: self
    type(particle_type), pointer :: molecule
    integer :: i
    real(dp), dimension(self%nx, self%ny, self%nz) :: fsum
    !---------------------------------------------------------------------------
    fsum = 0
    do i=1, self%nmolec
        molecule => self%molecules(i)%p
        fsum = fsum + molecule%f
    end do
end FUNCTION

!===============================================================================
FUNCTION f_of_molecules_with(self, other, count_double) result(fsum)
    class(ion_type) :: self
    type(ion_type) :: other
    type(particle_type), pointer :: molecule
    type(particle_pointer_t), allocatable :: mol_list(:)
    integer :: i
    real(dp), dimension(self%nx, self%ny, self%nz) :: fsum
    logical, optional :: count_double
    !---------------------------------------------------------------------------
    call self%molecules_with(other, mol_list, count_double)

    fsum = 0
    do i=1, size(mol_list)
        molecule => mol_list(i)%p
        fsum = fsum + molecule%f
    end do
end FUNCTION

!===============================================================================
FUNCTION e_of_molecules_with(self, other, count_double) result(esum)
    class(ion_type) :: self
    type(ion_type) :: other
    class(particle_type), pointer :: molecule
    type(particle_pointer_t), allocatable :: mol_list(:)
    integer :: i
    real(dp), dimension(self%nx, self%ny, self%nz) :: esum
    logical, optional :: count_double
    !---------------------------------------------------------------------------
    call self%molecules_with(other, mol_list, count_double)

    esum = 0
    do i=1, size(mol_list)
        molecule => mol_list(i)%p
        esum = esum + molecule%e
    end do
end FUNCTION

!===============================================================================
SUBROUTINE molecules_with (self, other, mol_list, count_double)
    class(ion_type), target, intent(in) :: self
    type(ion_type), intent(in) :: other
    logical, intent(in), optional :: count_double
    type(particle_pointer_t), intent(out), allocatable :: mol_list(:)
    class(particle_type), pointer :: molecule
    logical :: cd
    integer :: mi, mj, count
    !---------------------------------------------------------------------------
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
END SUBROUTINE molecules_with

SUBROUTINE allocate_ion (self, nx, ny, nz)
    class(ion_type) :: self
    integer :: nx, ny, nz
    !---------------------------------------------------------------------------
    if (self%allocated) then
        deallocate(self%f, self%saha, self%part, self%anjon, self%rho_per_h)
    endif
    
    self%nx = nx
    self%ny = ny
    self%nz = nz
    allocate(self%f(nx,ny,nz))
    allocate(self%saha(nx,ny,nz))
    allocate(self%part(nx,ny,nz))
    allocate(self%anjon(nx,ny,nz))
    allocate(self%rho_per_h(nx,ny,nz))
    self%allocated = .true.
END SUBROUTINE allocate_ion
    
end module ion_mod