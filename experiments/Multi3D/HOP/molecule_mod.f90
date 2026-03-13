MODULE molecule_mod
  USE os_mod
  USE prec_mod
  USE op_math_mod
  USE ion_mod
  USE handyfuncs_mod
  implicit none

  type, extends(particle_type) :: molec_type
    integer :: nions
    class(ion_pointer_t), allocatable :: ions(:)
    real(dp) :: dc(5) ! dissociation_constants
    real(dp), dimension(:), allocatable :: de ! dissociation_energies
    real(dp), dimension(:,:,:), allocatable :: dis
  contains
    procedure :: read_molec_data
    procedure :: get_dis
    procedure :: get_dis2
    procedure :: calc_f
    procedure :: calc_e
    procedure :: allocate => allocate_molec
    procedure :: deallocate => deallocate_molec
    procedure :: generate_isotopic_combinations
  end type

contains

  !===============================================================================
  SUBROUTINE read_molec_data(self, lunit, iostat)
    class(molec_type), intent(inout) :: self
    integer, intent(in) :: lunit
    integer, intent(out) :: iostat
    integer :: defix, n_irwin, n_de
    character(64) :: text
    !---------------------------------------------------------------------------
    read(lunit, *, iostat=iostat) text
    if (iostat < 0) return
    do while (text(1:1) == '!')
      read(lunit, *, iostat=iostat) text
    end do
    backspace lunit
    if (allocated(self%de)) deallocate(self%de)
    if (allocated(self%irwin)) deallocate(self%irwin)
    read(lunit, *, iostat=iostat) self%id, defix, n_irwin
    n_de = merge(6, 1, defix==0)
    allocate(self%de(n_de))
    allocate(self%irwin(n_irwin))
    read(lunit, *) self%de
    read(lunit, *) self%dc ! always 5 numbers
    read(lunit, *) self%irwin
    self%irwin(:) = self%irwin(n_irwin:1:-1)
  END SUBROUTINE read_molec_data

  !===============================================================================
  ! Computes the inverse dissociation constant P(AB) / (P(A)*P(B)) for molecules
  !===============================================================================
  elemental real(dp) FUNCTION get_dis(self, teta)
    class(molec_type), intent(in) :: self
    real(dp), intent(in) :: teta
    !---------------------------------------------------------------------------
    get_dis = 10**(-poly1d(teta, self%dc(1), self%dc(2), self%dc(3), self%dc(4), self%dc(5)))
  END FUNCTION

  !===============================================================================
  ! Computes the inverse dissociation constant P(AB) / (P(A)*P(B)) for molecules
  !===============================================================================
  SUBROUTINE get_dis2(self, temp, dis)
    class(molec_type), intent(in) :: self
    real(dp), intent(in), dimension(:,:,:) :: temp
    real(dp), intent(out), dimension(:,:,:) :: dis
    type(ion_type), pointer :: ion
    real(dp) :: awgt15, scalar, N_tot, const_term
    real(dp), parameter :: const_base = (2*PI*UU/HH**2)**(-1.5) * BK**(-2.5)
    integer :: j
    !---------------------------------------------------------------------------
    const_term = const_base**(self%nions - 1.0)

    dis = 1.d0
    do j=1, self%nions
      ion => self%ions(j)%p
      awgt15 = ion%awgt**(-1.5)
      dis = dis * awgt15 * ion%part
    end do

    scalar = const_term * self%awgt**(1.5)
    dis = dis * scalar
    dis = dis * exp(self%de(1)*EK/temp)
    dis = dis / (self%part * temp**(2.5 * (self%nions - 1.0)))
  END SUBROUTINE get_dis2

  !===============================================================================
  FUNCTION calc_f(self, pefe) result(f)
    class(molec_type), intent(in) :: self
    real(dp), dimension(:,:,:) :: pefe
    real(dp), dimension(:,:,:), allocatable :: f
    type(ion_type), pointer :: ion
    integer :: j
    !---------------------------------------------------------------------------
    allocate(f, source=self%dis)
    do j=1, self%nions
      ion => self%ions(j)%p
      f = f * ion%f
      ! pefe needed to convert from pressure to f, first pefe cancels out
      if (j > 1) f = f * pefe
    end do
  END FUNCTION calc_f

  !===============================================================================
  FUNCTION calc_e(self, teta, temp, dxi) result(e)
    class(molec_type), intent(in) :: self
    class(ion_type), pointer :: ion
    real(dp), dimension(:,:,:) :: teta, temp, dxi
    real(dp), dimension(:,:,:), allocatable :: de, e
    integer :: i
    !---------------------------------------------------------------------------
    allocate(de, e, mold=self%n)

    if (size(self%de) > 1) then
      ! rotational and vibrational energy has to be subtracted
      de = self%de(1) - poly1d(teta, self%de(2), self%de(3), self%de(4), self%de(5), self%de(6))*temp/EK
    else
      de = self%de(1)
    endif

    do i=1, self%nions
      ion => self%ions(i)%p

      ! why only consider ionisation energy from H?
      if (trim(ion%id) == 'H I') then
        de = de + ion%xion - dxi
      endif
    end do

    e = - de * self%n
  END FUNCTION calc_e

  !===============================================================================
  SUBROUTINE allocate_molec (self, m)
    class(molec_type) :: self
    integer :: m(3), mold(3), s, p
    !---------------------------------------------------------------------------
    if (all(m == self%m)) return
    self%m = m
    s = storage_size(self%f)
    if (self%allocated) then
      mold = shape(self%f)
      p = product(mold)
      deallocate(self%f, self%dis, self%part, self%n)
      call os%mem_register(s, -p, label='mol%f')
      call os%mem_register(s, -p, label='mol%n')
      call os%mem_register(s, -p, label='mol%dis')
      call os%mem_register(s, -p, label='mol%part')
    endif
    p = product(m)
    allocate(self%f(m(1),m(2),m(3)), source=0.d0)
    allocate(self%n, self%dis, self%part, source=self%f)
    call os%mem_register(s, p, label='mol%f')
    call os%mem_register(s, p, label='mol%n')
    call os%mem_register(s, p, label='mol%dis')
    call os%mem_register(s, p, label='mol%part')
    self%allocated = .true.
  END SUBROUTINE allocate_molec

  !===============================================================================
  SUBROUTINE deallocate_molec (self)
    class(molec_type) :: self
    integer :: m(3), s, p
    !---------------------------------------------------------------------------
    m = self%m
    p = product(m)
    s = storage_size(self%f)
    if (self%allocated) then
      deallocate(self%f, self%dis)
      call os%mem_register(s, -p, label='mol%f')
      call os%mem_register(s, -p, label='mol%dis')
      if (.not. (self%in_opac .or. self%in_linelist)) then
        deallocate(self%n, self%part)
        call os%mem_register(s, -p, label='mol%n')
        call os%mem_register(s, -p, label='mol%part')
      endif
    endif
    self%allocated = .false.
  END SUBROUTINE deallocate_molec

  !===============================================================================
  SUBROUTINE generate_isotopic_combinations(self)
    class(molec_type), intent(inout) :: self
    class(ion_type), pointer :: ion
    class(isotope_type), pointer :: iso, iso2
    integer :: i, j, k
    integer, allocatable :: current_iso(:), order(:)
    real(8) :: total_mass, total_frac
    character(16) :: mass_id
    logical :: done
    !---------------------------------------------------------------------------
    ! Estimate number of isotopic combinations (upper bound)
    self%niso = 1
    do i=1, self%nions
      ion => self%ions(i)%p
      if (associated(ion%isotope)) then
        self%niso = self%niso * (ion%niso+1)
      endif
    end do

    allocate(self%isotope(self%niso))

    ! Initialize index tracking for each ion
    allocate(current_iso(self%nions))
    allocate(order(self%nions))
    call argsort_int([(nint(self%ions(i)%p%awgt), i=1, self%nions)], order)
    current_iso = [(merge(0,1,associated(self%ions(order(i))%p%isotope)), i=1, self%nions)]
    k = 1
    done = .false.

    ! Generate all isotopic combinations
    do while (.not. done)
      ! Compute total mass and fraction
      total_mass = 0.0d0
      total_frac = 1.0d0
      mass_id = ''

      do i=1, self%nions
        ion => self%ions(order(i))%p
        if (.not. associated(ion%isotope)) then
          total_mass = total_mass + ion%awgt
          write(mass_id, '(A, I0)') trim(mass_id), nint(ion%awgt)
        else
          if (current_iso(i)==0) then
            total_mass = total_mass + ion%awgt
            write(mass_id, '(A, I0)') trim(mass_id), 0
          else
            iso => ion%isotope(current_iso(i))
            total_mass = total_mass + iso%awgt
            total_frac = total_frac * iso%frac
            write(mass_id, '(A, I0)') trim(mass_id), nint(iso%awgt)
          endif
        endif
      end do

      ! Store computed isotope in molecule
      iso => self%isotope(k)
      iso%awgt = total_mass
      iso%frac = total_frac
      read(mass_id, *) iso%iso
      k = k + 1

      ! Generate next combination
      do i=self%nions, 1, -1
        ion => self%ions(order(i))%p
        current_iso(i) = current_iso(i) + 1
        if (current_iso(i) <= ion%niso) exit
        current_iso(i) = merge(0,1,associated(ion%isotope))
        if (i == 1) done = .true.
      end do
    end do
    ! sum up duplicates
    ! this isn't quite perfect yet as it generates extra isotopes
    ! for H2, C2, etc. with frac > 1, as e.g. (0)12 and 120 are summed.
    ! It's fine as long as we don't use them though.
    ! Also haven't tested this for triatomic molecules yet which do not contain H.
    do i=1, self%niso
      iso => self%isotope(i)
      do k=i+1, self%niso
        iso2 => self%isotope(k)
        if (iso%awgt == iso2%awgt .and. iso%frac == iso2%frac) then
          iso%frac = iso%frac + iso2%frac
          iso2%frac = iso%frac
        endif
      end do
    end do
    deallocate(current_iso)
  END SUBROUTINE generate_isotopic_combinations

END MODULE molecule_mod
