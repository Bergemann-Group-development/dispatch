module transition_mod
    use prec_mod
    use op_math_mod
    use level_mod
    use ion_mod
    use element_mod
    use molecule_mod
    implicit none
    
    type :: transition_type
        integer :: j, k, ntt, nnu, fs, fe, di
        type(level_type) :: level
        class(particle_type), pointer :: ion
        real(dp) :: mint, maxt, minl, maxl
        logical :: ilogl, kvadl, minex, maxex
        logical :: ilogt, kvadt, minet, maxet, iteta
        integer, dimension(:), allocatable :: iw
        real(dp), dimension(:), allocatable :: lam, temp, teta, tt
        real(dp), dimension(:,:), allocatable :: xkap, yy2d
        real(dp), dimension(:,:), allocatable :: log_xkap
        real(dp), dimension(:,:), allocatable :: lagrange
        logical, dimension(:), allocatable :: quad, ldisk!, tdisk
        logical :: bound_free, regular, active
        character(len=5) :: key
        character(len=20) :: label

    contains
        procedure :: get_temp_index
        procedure :: interp_xkap
        procedure :: check_temp_quadratic
        procedure :: check_lam_quadratic
        procedure :: boltzmann_weight => trans_boltzmann_weight
        procedure :: absorption_factor
        procedure :: add_absorption
        procedure :: set_xweight
        procedure :: load_lam_interp_weights
        procedure :: interp_temp_quadratic
        procedure :: interp_temp_linear
    end type

contains

!===============================================================================
subroutine load_lam_interp_weights(self,lam)
    class(transition_type) :: self
    real(dp), dimension(:), intent(in) :: lam
    real(dp), dimension(:), allocatable :: test
    integer :: nf, f, fs, fe, di
    integer :: nx, j, i
    logical :: quadratic
    !---------------------------------------------------------------------------
    self%active = .true.

    nx = size(self%lam)
    nf = size(lam)

    if (lam(1) > lam(nf)) then
        di = -1
        fs = nf
        fe = 1
    else
        di = 1
        fs = 1
        fe = nf
    endif

    if (.not. self%minex .and. self%minl > lam(fs)) then
        if (lam(fe) < self%minl) then
            self%active = .false.
            return
        endif
        do f=fs,fe,di
            if (lam(f) > self%minl) then
                fs = f
                exit
            endif
        end do
    endif

    if (.not. self%maxex .and. self%maxl < lam(fe)) then
        if (lam(fs) > self%maxl) then
            self%active = .false.
            return
        endif
        do f=fe,fs,-di
            if (lam(f) < self%maxl) then
                fe = f
                exit
            endif
        end do
    endif

    j = 1
    f = fs

    quadratic = self%check_lam_quadratic(j)
    self%quad(f) = quadratic

    do i=1, nf+nx ! i is not actually needed
        if (j >= nx) then
            self%iw(f) = nx

        elseif (self%lam(j+1) < lam(f)) then
            j = j+1
            quadratic = self%check_lam_quadratic(j)
            cycle

        elseif (self%lam(j) < lam(f)) then
            self%iw(f) = j

        else
            self%iw(f) = 1
        endif

        call self%set_xweight(lam, f)

        if (f == fe) exit
        f = f+di
        self%quad(f) = quadratic
    end do

    self%fs = fs
    self%fe = fe
    self%di = di
end subroutine

!===============================================================================
subroutine set_xweight(self, lam, f)
    class(transition_type) :: self
    real(dp), dimension(:), intent(in) :: lam
    integer, intent(in) :: f
    integer :: j
    !---------------------------------------------------------------------------
    j = self%iw(f)

    if (self%quad(f)) then
        if (size(self%lam) - j < 2) then
            ! if there are not enough points, step back one index
            self%iw(f) = self%iw(f) - 1
            j = self%iw(f)
        endif

        self%lagrange(1:3,f) = quadratic_weight(self%lam(j:j+2), lam(f))
    else
        self%lagrange(1:2,f) = linear_weight(self%lam(j:j+1), lam(f))
    endif
end subroutine set_xweight

!===============================================================================
function quadratic_weight(xxx,newx) result(w)
    real(dp), dimension(3), intent(in) :: xxx
    real(dp), intent(in) :: newx
    real(dp) :: w(3)
    !---------------------------------------------------------------------------
    w(1) = (newx-xxx(2)) * (newx-xxx(3)) / ((xxx(1)-xxx(2)) * (xxx(1)-xxx(3)))
    w(2) = (newx-xxx(1)) * (newx-xxx(3)) / ((xxx(2)-xxx(1)) * (xxx(2)-xxx(3)))
    w(3) = (newx-xxx(1)) * (newx-xxx(2)) / ((xxx(3)-xxx(1)) * (xxx(3)-xxx(2)))
end function quadratic_weight

!===============================================================================
function linear_weight(xx,newx) result(w)
    real(dp), dimension(2), intent(in) :: xx
    real(dp), intent(in) :: newx
    real(dp) :: w(2)
    !---------------------------------------------------------------------------
    w(1) = (xx(2) - newx) / (xx(2) - xx(1))
    w(2) = (newx - xx(1)) / (xx(2) - xx(1))
end function linear_weight

!===============================================================================
subroutine interp_xkap(self, abslte, lam, temp, teta, fakt)
    class(transition_type) :: self
    real(dp), dimension(:), intent(in) :: lam
    real(dp), dimension(:,:,:,:), intent(inout) :: abslte
    real(dp), intent(in), dimension(:,:,:), target :: temp, teta
    real(dp), intent(in), dimension(:,:,:) :: fakt
    real(dp), dimension(:,:,:), pointer :: newt
    real(dp), dimension(:,:,:,:), allocatable :: yyt
    real(dp), dimension(:,:,:), allocatable :: ab, yt
    real(dp) :: ab1
    real(dp), dimension(self%nnu) :: yy
    real(dp), dimension(self%ntt) :: tt
    real(dp), dimension(:,:), allocatable :: yy2d, dum
    integer :: i, fs, fe, f, di, nx, ny, nz, nf, iw
    integer, dimension(:,:,:), allocatable :: it
    logical, dimension(:,:,:), allocatable :: quadratic, in_bounds
    !---------------------------------------------------------------------------
    fs = self%fs
    fe = self%fe
    di = self%di

    nx = size(temp, 1)
    ny = size(temp, 2)
    nz = size(temp, 3)
    nf = size(lam)
    
    if (self%ntt > 1) then
        allocate(yyt(nx,ny,nz,self%nnu))
        allocate(ab, yt, mold=temp)
        allocate(quadratic(nx,ny,nz))
        allocate(in_bounds(nx,ny,nz))
        allocate(it(nx,ny,nz))

        ab = 0

        in_bounds = .true.

        if (.not. self%minet) then
            where (temp < self%mint)
                in_bounds = .false.
            end where
        endif

        if (.not. self%maxet) then
            where (temp > self%maxt)
                in_bounds = .false.
            end where
        endif

        if (.not. all(in_bounds)) then
            print *, 'temperature out of bounds in continuum ', self%key
            print *, 'this could be a problem, please check'
        endif

        if (self%ilogt) then
            self%yy2d = self%log_xkap
        else
            self%yy2d = self%xkap
        endif

        if (self%iteta) then
            newt => teta
        else
            newt => temp
        endif

        call self%get_temp_index(newt, it, quadratic)

        where (quadratic) 
            quadratic = self%check_temp_quadratic(it)
        end where

        do f=1, self%nnu
            where (quadratic .and. in_bounds)
                yt = self%interp_temp_quadratic(f, it, newt)
            elsewhere (in_bounds)
                yt = self%interp_temp_linear(f, it, newt)
            elsewhere
                yt = 0
            end where

            if (self%ilogt .and. .not. self%ilogl) then
                where (in_bounds) yt = exp(yt)
            elseif(self%ilogl .and. .not. self%ilogt) then
                where (in_bounds) yt = log(yt)
            endif

            yyt(:,:,:,f) = yt
        end do

        do f=fs,fe,di
            iw = self%iw(f)
            if (self%quad(f)) then
                where (in_bounds)
                    ab = yyt(:,:,:,iw)   * self%lagrange(1,f) + &
                            yyt(:,:,:,iw+1) * self%lagrange(2,f) + &
                            yyt(:,:,:,iw+2) * self%lagrange(3,f)
                end where
            else
                where (in_bounds)
                    ab = yyt(:,:,:,iw)   * self%lagrange(1,f) + & 
                            yyt(:,:,:,iw+1) * self%lagrange(2,f)
                end where
            endif

            if (self%ilogl) where (in_bounds) ab = exp(ab)

            where (in_bounds)
                abslte(:,:,:,f) = abslte(:,:,:,f) + fakt*ab
            end where

            !if (f == 1) print '(A20,4(E18.10))', self%label, fakt(1,1,230)*ab(1,1,230), fakt(1,1,230), ab(1,1,230), abslte(1,1,230,f)
        end do
    else
        if(self%ilogl) then
            yy = self%log_xkap(:,1)
        else
            yy = self%xkap(:,1)
        endif

        do f=fs,fe,di
            iw = self%iw(f)
            if (self%quad(f)) then
                ab1 = sum(yy(iw:iw+2) * self%lagrange(1:3,f))
            else
                ab1 = sum(yy(iw:iw+1) * self%lagrange(1:2,f))
            endif

            if (self%ilogl) ab1 = exp(ab1)
            abslte(:,:,:,f) = abslte(:,:,:,f) + fakt*ab1

            !if (f == 1) then
            !    print '(A20,4(E18.10))', self%label, fakt(1,1,230)*ab1, fakt(1,1,230), ab1, abslte(1,1,230,f)
            !endif
        end do
    endif
end subroutine interp_xkap

!===============================================================================
elemental function interp_temp_quadratic(self, f, it, newt) result(newy)
    class(transition_type), intent(in) :: self
    real(dp), intent(in) :: newt
    integer, intent(in) :: f, it
    real(dp) :: newy, w(3)
    integer :: j1, j2, j3
    !---------------------------------------------------------------------------
    if (self%ntt - it > 1) then
        j1 = it
        j2 = it+1
        j3 = it+2
    else
        j1 = it-1
        j2 = it
        j3 = it+1
    endif

    w(1) = (newt-self%tt(j2)) * (newt-self%tt(j3)) / ((self%tt(j1)-self%tt(j2)) * (self%tt(j1)-self%tt(j3)))
    w(2) = (newt-self%tt(j1)) * (newt-self%tt(j3)) / ((self%tt(j2)-self%tt(j1)) * (self%tt(j2)-self%tt(j3)))
    w(3) = (newt-self%tt(j1)) * (newt-self%tt(j2)) / ((self%tt(j3)-self%tt(j1)) * (self%tt(j3)-self%tt(j2)))

    newy = sum(self%yy2d(f,j1:j3) * w(1:3))
end function interp_temp_quadratic 

!===============================================================================
elemental function interp_temp_linear(self, f, it, newt) result(newy)
    class(transition_type), intent(in) :: self
    real(dp), intent(in) :: newt
    integer, intent(in) :: f, it
    real(dp) :: m, newy
    !---------------------------------------------------------------------------
    m = (self%yy2d(f,it+1) - self%yy2d(f,it)) / (self%tt(it+1) - self%tt(it))
    newy = m * (newt - self%tt(it)) + self%yy2d(f,it)
end function interp_temp_linear

!===============================================================================
elemental function check_lam_quadratic(self, iw) result(quadratic)
    class(transition_type), intent(in) :: self
    integer, intent(in) :: iw
    logical :: quadratic
    !---------------------------------------------------------------------------
    quadratic = .false.

    if (self%kvadl) then
        if (self%nnu - iw > 1) then
            if (.not. self%ldisk(iw)) then
                quadratic = .true.
            endif
        endif
    endif
end function check_lam_quadratic

!===============================================================================
elemental function check_temp_quadratic(self, it) result(quadratic)
    class(transition_type), intent(in) :: self
    integer, intent(in) :: it
    logical :: quadratic
    !---------------------------------------------------------------------------
    quadratic = .false.

    if (self%kvadt) then
        if (self%ntt - it > 1) then
            !if (.not. self%tdisk(it)) then
            quadratic = .true.
            !endif
        endif
    endif
end function check_temp_quadratic

!===============================================================================
elemental subroutine get_temp_index(self, newt, idx, no_extrap)
    class(transition_type), intent(in) :: self
    real(dp), intent(in) :: newt
    integer,  intent(out) :: idx
    logical,  intent(out) :: no_extrap
    !---------------------------------------------------------------------------
    no_extrap = .false.

    if (newt < self%tt(1)) then
        idx = 1
    elseif (newt > self%tt(self%ntt)) then
        idx = self%ntt - 1
    else
        idx = BinarySearch(self%tt, 1, self%ntt, newt)
        no_extrap = .true.
    endif
end subroutine get_temp_index

!===============================================================================
pure function trans_boltzmann_weight(self, EK_over_temp) result(frac)
    class(transition_type), intent(in) :: self
    real(dp), dimension(:,:,:), intent(in) :: EK_over_temp
    real(dp), dimension(:,:,:), allocatable :: frac
    !---------------------------------------------------------------------------
    frac = self%level%g2 * self%level%boltzmann(EK_over_temp, natural_log=.true.)
end function trans_boltzmann_weight

!===============================================================================
function absorption_factor(self, EK_over_temp, teta, pe, rho, dxi, xfakh) result(fakt)
    class(transition_type), intent(in) :: self
    class(particle_type), pointer :: ion
    real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, teta, pe, rho, dxi, xfakh
    real(dp), dimension(:,:,:), allocatable :: teta31, fakt, h1
    integer, dimension(:,:,:), allocatable :: nniv
    !---------------------------------------------------------------------------
    ion => self%ion

    if (self%regular) then
        fakt = self%boltzmann_weight(EK_over_temp) / ion%part * ion%rho_per_h
    else
        select case(trim(self%key))
        case('H- BF')
            fakt = 1.e-17 * ion%rho_per_h 
        case('H- FF')
            fakt = 2.E-26 / ion%part * pe * ion%rho_per_h
        case('H BF')
            if (self%k <= 5) then
                fakt = self%boltzmann_weight(EK_over_temp)
                fakt = 1.045E-26 / ion%part * fakt * ion%rho_per_h / self%k**5
            else
                teta31 = 31.30364*teta
                fakt = 2.0898E-26 / ion%part * xfakh
                fakt = fakt * exp(-teta31*(1. - 1./self%k**2)) / self%k**3
            endif
            h1 = sqrt(13.595/dxi)
            nniv = int(h1)
            where (nniv < self%k)
                fakt = 0
            elsewhere (nniv == self%k)
                fakt = fakt * min(h1 - float(nniv), 1.0)
            end where

        case('H FF')
            teta31 = 31.30364*teta
            fakt = 2.0898E-26 / ion%part * xfakh
            fakt = fakt * exp(-teta31 + 2.3026 * dxi * teta) / (2*teta31)

        case('(H+H)')
            fakt = (1e-25 * ion%rho_per_h)**2 * rho
        case('H2+')
            fakt = (1e-20 * ion%rho_per_h)**2 * rho
        case('H2-', 'HE-')
            fakt = pe * ion%rho_per_h
        end select
    endif
    !print *, self%key, fakt
end function absorption_factor

!===============================================================================
subroutine add_absorption(self, abslte, lam, temp, EK_over_temp, teta, pe, rho, dxi, xfakh, fakt)
    use timer_mod
    class(transition_type) :: self
    real(dp), dimension(:,:,:,:), intent(inout) :: abslte
    real(dp), dimension(:,:,:), intent(in) :: temp, EK_over_temp, teta, pe, rho, dxi, xfakh
    real(dp), dimension(:), intent(in) :: lam
    real(dp), dimension(:), allocatable :: ab, abkof
    real(dp), dimension(:,:,:), allocatable, intent(out) :: fakt
    integer :: i,j,k,nf
    !---------------------------------------------------------------------------

    nf = size(lam)
    allocate(self%quad(nf))
    allocate(self%iw(nf))
    allocate(self%lagrange(3,nf))

    self%lagrange = 0

    call self%load_lam_interp_weights(lam)

    if (self%active) then
        fakt = self%absorption_factor(EK_over_temp, teta, pe, rho, dxi, xfakh)
        call self%interp_xkap(abslte, lam, temp, teta, fakt)
    endif

    !print '(A20, E18.10)', self%label, abslte(1,1,1,1)

    deallocate(self%quad)
    deallocate(self%iw)
    deallocate(self%lagrange)
end subroutine add_absorption
end module