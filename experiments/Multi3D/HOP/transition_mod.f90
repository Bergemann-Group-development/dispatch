module transition_mod
    use prec_mod
    use op_math_mod
    use level_mod
    use ion_mod
    use element_mod
    implicit none
    
    type :: transition_type
        integer :: j, k, ntt, nnu, fs, fe, di
        type(level_type) :: level
        class(particle_type), pointer :: ion
        real(dp) :: mint, maxt, minl, maxl
        logical :: ilogl, kvadl, minex, maxex
        logical :: ilogt, kvadt, minet, maxet, iteta
        integer, dimension(:), allocatable :: iw
        real(dp), dimension(:), allocatable :: lam, temp, teta
        real(dp), dimension(:,:), allocatable :: xkap
        real(dp), dimension(:,:), allocatable :: log_xkap
        real(dp), dimension(:,:), allocatable :: lagrange
        logical, dimension(:), allocatable :: quad, ldisk!, tdisk
        logical :: bound_free, regular, active
        character(len=5) :: key

    contains
        procedure :: interp_xkap
        procedure :: check_temp_quadratic
        procedure :: check_lam_quadratic
        procedure :: boltzmann_weight => trans_boltzmann_weight
        procedure :: absorption_factor
        procedure :: add_absorption
        procedure :: set_xweight
        procedure :: load_lam_interp_weights
    end type

contains

    subroutine load_lam_interp_weights(self,lam)
        use io_mod
        class(transition_type) :: self
        real(dp), dimension(:), intent(in) :: lam
        real(dp), dimension(:), allocatable :: test
        integer :: nf, f, fs, fe, di
        integer :: nx, j, i
        logical :: quadratic
        
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

    subroutine set_xweight(self, lam, f)
        class(transition_type) :: self
        real(dp), dimension(:), intent(in) :: lam
        integer, intent(in) :: f
        integer :: j

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

    function quadratic_weight(xxx,newx) result(w)
        real(dp), dimension(3), intent(in) :: xxx
        real(dp), intent(in) :: newx
        real(dp) :: w(3)

        w(1) = (newx-xxx(2)) * (newx-xxx(3)) / ((xxx(1)-xxx(2)) * (xxx(1)-xxx(3)))
        w(2) = (newx-xxx(1)) * (newx-xxx(3)) / ((xxx(2)-xxx(1)) * (xxx(2)-xxx(3)))
        w(3) = (newx-xxx(1)) * (newx-xxx(2)) / ((xxx(3)-xxx(1)) * (xxx(3)-xxx(2)))
    end function quadratic_weight

    function linear_weight(xx,newx) result(w)
        real(dp), dimension(2), intent(in) :: xx
        real(dp), intent(in) :: newx
        real(dp) :: w(2)

        w(1) = (xx(2) - newx) / (xx(2) - xx(1))
        w(2) = (newx - xx(1)) / (xx(2) - xx(1))
    end function linear_weight

    pure function interp_xkap(self, lam, temp) result(xkap)
        class(transition_type), intent(in) :: self
        real(dp), dimension(:), intent(in) :: lam
        real(dp), intent(in) :: temp
        real(dp) :: newt
        real(dp), dimension(size(lam)) :: xkap
        real(dp), dimension(self%nnu) :: yy
        real(dp), dimension(self%ntt) :: tt
        real(dp), dimension(:,:), allocatable :: yy2d
        integer :: i, it, iw, fs, fe, f, di
        logical :: no_extrap, quadratic

        xkap = 0

        if (self%ntt > 1) then
            if (     (.not. self%minet .and. temp < self%mint)      &
                .or. (.not. self%maxet .and. temp > self%maxt)) then
                return
            endif
        endif
        
        if (self%ntt > 1) then
            if (self%ilogt) then
                yy2d = self%log_xkap
            else
                yy2d = self%xkap
            endif

            if (self%iteta) then
                tt = self%teta
                newt = temp_to_theta/temp
            else
                tt = self%temp
                newt = temp
            endif

            call get_interp_index(tt, newt, it, no_extrap)

            quadratic = .false.

            if (no_extrap) then
                quadratic = self%check_temp_quadratic(it)
                if (quadratic) then
                    do i=1, self%nnu
                        yy(i) = interp1d_quadratic(it, tt, yy2d(i,:), newt)
                    end do
                endif
            endif

            if (.not. quadratic) then
                do i=1, self%nnu
                    yy(i) = interp1d_linear(it, tt, yy2d(i,:), newt)
                end do
            endif

            if (self%ilogt .and. .not. self%ilogl) then
                yy = exp(yy)
            elseif(self%ilogl .and. .not. self%ilogt) then
                yy = log(yy)
            endif
        else
            if(self%ilogl) then
                yy = self%log_xkap(:,1)
            else
                yy = self%xkap(:,1)
            endif
        endif

        fs = self%fs
        fe = self%fe
        di = self%di

        do f=fs,fe,di
            iw = self%iw(f)
            if (self%quad(f)) then
                xkap(f) = sum(yy(iw:iw+2) * self%lagrange(1:3,f))
            else
                xkap(f) = sum(yy(iw:iw+1) * self%lagrange(1:2,f))
            endif
        end do

        if (self%ilogl) xkap(fs:fe:di) = exp(xkap(fs:fe:di))
    end function interp_xkap

    pure function check_lam_quadratic(self, iw) result(quadratic)
        class(transition_type), intent(in) :: self
        integer, intent(in) :: iw
        logical :: quadratic

        quadratic = .false.

        if (self%kvadl) then
            if (self%nnu - iw > 1) then
                if (.not. self%ldisk(iw)) then
                    quadratic = .true.
                endif
            endif
        endif
    end function check_lam_quadratic

    pure function check_temp_quadratic(self, it) result(quadratic)
        class(transition_type), intent(in) :: self
        integer, intent(in) :: it
        logical :: quadratic

        quadratic = .false.

        if (self%kvadt) then
            if (self%ntt - it > 1) then
                !if (.not. self%tdisk(it)) then
                quadratic = .true.
                !endif
            endif
        endif
    end function check_temp_quadratic

    pure subroutine get_interp_index(xx, newx, idx, no_extrap)
        real(dp), dimension(:), intent(in) :: xx
        real(dp), intent(in) :: newx
        integer,  intent(out) :: idx
        logical,  intent(out) :: no_extrap
        integer :: nx

        nx = size(xx)
        no_extrap = .false.

        if (newx < xx(1)) then
            idx = 1
        elseif (newx > xx(nx)) then
            idx = nx - 1
        else
            idx = BinarySearch(xx, 1, nx, newx)
            no_extrap = .true.
        endif
    end subroutine get_interp_index

    pure function trans_boltzmann_weight(self, EK_over_temp) result(frac)
        class(transition_type), intent(in) :: self
        real(dp), dimension(:,:,:), intent(in) :: EK_over_temp
        real(dp), dimension(:,:,:), allocatable :: frac

        frac = self%level%g2 * self%level%boltzmann(EK_over_temp, natural_log=.true.)
    end function trans_boltzmann_weight

    function absorption_factor(self, EK_over_temp, teta, pe, rho, dxi, xfakh) result(fakt)
        class(transition_type), intent(in) :: self
        type(ion_type), pointer :: ion
        type(molec_type), pointer :: H2
        real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, teta, pe, rho, dxi, xfakh
        real(dp), dimension(:,:,:), allocatable :: teta31, fakt

        select type(ptr => self%ion)
        class is(ion_type)
            ion => ptr
        class is(molec_type)
            H2 => ptr
        end select

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
                    !fakt = 2.0898E-26 / ion%part * ion%rho_per_h
                    fakt = 2.0898E-26 / ion%part * xfakh
                    fakt = fakt * exp(-teta31*(1. - 1./self%k**2)) / self%k**3
                endif
            case('H FF')
                teta31 = 31.30364*teta
                !fakt = 2.0898E-26 / ion%part * ion%rho_per_h
                fakt = 2.0898E-26 / ion%part * xfakh
                fakt = fakt * exp(-teta31 + 2.3026 * dxi * teta) / (2*teta31)

            case('(H+H)')
                fakt = (1e-25 * ion%rho_per_h)**2 * rho
            case('H2+')
                fakt = (1e-20 * ion%rho_per_h)**2 * ion%saha * rho
            case('H2-')
                fakt = pe * H2%rho_per_h
            case('HE-')
                fakt = pe * ion%rho_per_h
            end select
        endif
        !print *, self%key, fakt
    end function absorption_factor

    subroutine add_absorption(self, abslte, lam, temp, EK_over_temp, teta, pe, rho, dxi, xfakh)
        use timer_mod
        class(transition_type) :: self
        real(dp), dimension(:,:,:,:), intent(inout) :: abslte
        real(dp), dimension(:,:,:), intent(in) :: temp, EK_over_temp, teta, pe, rho, dxi, xfakh
        real(dp), dimension(:), intent(in) :: lam
        real(dp), dimension(:), allocatable :: ab, abkof
        real(dp), dimension(:,:,:), allocatable :: fakt
        integer :: i,j,k,nf

        fakt = self%absorption_factor(EK_over_temp, teta, pe, rho, dxi, xfakh)

        nf = size(lam)
        allocate(self%quad(nf))
        allocate(self%iw(nf))
        allocate(self%lagrange(3,nf))

        call self%load_lam_interp_weights(lam)

        if (self%active) then
            do k=1, size(temp,3)
            do j=1, size(temp,2)
            do i=1, size(temp,1)
                abkof = self%interp_xkap(lam, temp(i,j,k))
                ab = fakt(i,j,k) * abkof
                abslte(i,j,k,:) = abslte(i,j,k,:) + ab
                !if (k==1) then
                !    print *, self%key, ab(1), fakt(i,j,k), abkof(1)
                !endif
            end do
            end do
            end do
        endif

        deallocate(self%quad)
        deallocate(self%iw)
        deallocate(self%lagrange)
    end subroutine add_absorption
end module