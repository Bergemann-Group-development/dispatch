MODULE molecule_mod
    USE prec_mod
    USE op_math_mod
    USE ion_mod
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
        procedure :: init_mol_dissociation
    end type
    
contains

!===============================================================================
! Computes dissociation constants N(AB) / (N(A)*N(B)) for molecules
!===============================================================================
elemental real(dp) FUNCTION pk(self, teta)
    class(molec_type), intent(in) :: self
    real(dp), intent(in) :: teta
    !---------------------------------------------------------------------------
    pk = 10**(-poly1d(teta, self%dc(1), self%dc(2), self%dc(3), self%dc(4), self%dc(5)))
END FUNCTION

!===============================================================================
FUNCTION calc_f(self, pe_over_fe) result(f)
    class(molec_type), intent(in) :: self
    real(dp), dimension(:,:,:) :: pe_over_fe
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
END FUNCTION calc_f

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
        if (trim(self%id) == 'H2+') dis = - (dis - 2.648)
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
END FUNCTION calc_e

!===============================================================================
SUBROUTINE allocate_molec (self, nx, ny, nz)
    class(molec_type) :: self
    integer :: nx, ny, nz
    !---------------------------------------------------------------------------
    if (self%allocated) then
        deallocate(self%f, self%e, self%dis, self%n)
    endif

    self%nx = nx
    self%ny = ny
    self%nz = nz
    allocate(self%f(nx,ny,nz))
    allocate(self%e(nx,ny,nz))
    allocate(self%n(nx,ny,nz))
    allocate(self%dis(nx,ny,nz))
    self%f = 0
    self%e = 0
    self%dis = 0
    self%n = 0
    self%allocated = .true.
END SUBROUTINE allocate_molec

!===============================================================================
SUBROUTINE init_mol_dissociation (self)
    class(molec_type), target :: self
    integer :: i
    !---------------------------------------------------------------------------
    select case (trim(self%id))
    case ('H2', 'H2+')
        allocate(self%de(5))
    case default
        allocate(self%de(1))
    end select

    select case (trim(self%id))
    case ('H2')
        self%dc = [ 6.3021E-4, -1.4149E-2,  1.2572E-1, -5.1172E0,  12.739E0 ]
        self%de = [  0.00975,  -0.12427,   0.60602,   -1.4772,    2.6757 ]
    case ('H2+')
        self%dc = [ 0.0000e+00,  0.024790, -7.9196803E-2, -2.7942767E0,  11.206998E0 ]
        self%de = [  0.15253,  -0.82685,    1.7231,   -2.0036,    2.9216 ]
    case ('H2O')
        self%dc = [ 8.1730E-04, -1.8368E-02,  1.6939E-1, -10.522E0,  25.420E0 ]
        self%de = [       9.50 ]
    case ('OH')
        self%dc = [ 7.7224E-04, -1.6547E-02,  1.3822E-01, -5.0578E0,  12.371E0 ]
        self%de = [       4.38 ]
    case ('CH')
        self%dc = [ 7.2661E-04, -1.5473E-02,  1.2768E-01, -4.0760E0,  12.135E0 ]
        self%de = [       3.47 ]
    case ('CO')
        self%dc = [ 1.1349E-03, -2.2888E-02,  1.7217E-01, -11.795E0,  13.820E0 ]
        self%de = [      11.11 ]
    case ('CN')
        self%dc = [ 3.4666E-04, -7.3627E-03,  6.4162E-02, -8.2793E0,  12.805E0 ]
        self%de = [       7.90 ]
    case ('C2')
        self%dc = [ 6.2603E-04, -1.2739E-02,  9.7719E-02, -6.5178E0,  12.804E0 ]
        self%de = [       6.12 ]
    case ('N2')
        self%dc = [ 1.4993E-03, -2.9997E-02,  2.2067E-01, -10.585E0,  13.590E0 ]
        self%de = [       9.76 ]
    case ('O2')
        self%dc = [ 3.7970E-04, -8.1511E-03,  6.9935E-02, -5.5181E0,  13.228E0 ]
        self%de = [       5.12 ]
    case ('NO')
        self%dc = [ 1.1380E-03, -2.3065E-02,  1.7349E-01, -7.1964E0,  12.831E0 ]
        self%de = [       6.51 ]
    case ('NH')
        self%dc = [ 7.8691E-04, -1.6643E-02,  1.3629E-01, -3.8435E0,  12.033E0 ]
        self%de = [       3.21 ]
    end select
END SUBROUTINE init_mol_dissociation
END MODULE molecule_mod