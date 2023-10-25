MODULE bftrans_mod
  USE prec_mod
  USE op_math_mod
  USE level_mod, only: level_type, wcalc, wcalc2
  USE ion_mod
  USE element_mod
  USE molecule_mod
  USE os_mod
  USE handyfuncs_mod
  USE trace_mod
  implicit none
  private
  
  type, public :: bftrans_type
      integer :: k, ntt, nnu, verbose
      type(level_type) :: level
      class(particle_type), pointer :: ion
      real(dp) :: mint, maxt, minl, maxl
      logical :: ilogl, kvadl, minex, maxex, loglam=.false.
      logical :: ilogt, kvadt, minet, maxet, iteta, pedep, istim
      real(dp), dimension(:), allocatable :: lam, temp, teta, tt
      real(dp), dimension(:,:), allocatable :: xkap
      real(dp), dimension(:,:), allocatable :: log_xkap
      logical, dimension(:), allocatable :: ldisk!, tdisk
      logical :: free_free=.false., regular=.false., molecular=.false.
      logical :: active=.true.
      character(len=5) :: key
      character(len=20) :: label

  contains
      procedure :: get_temp_index
      procedure :: interp_xkap
      procedure :: check_temp_quadratic
      procedure :: check_lam_quadratic
      procedure :: absorption_factor
      procedure :: add_absorption
      procedure :: set_xweight
      procedure :: load_lam_interp_weights
      procedure :: load_unsorted_lam_interp_weights
      procedure :: interp_temp_quadratic
      procedure :: interp_temp_linear
      procedure :: H_occupation
  end type

  logical :: active
  integer :: fs, fe
  logical, public :: H_occupation_probability=.true.
  integer, dimension(:), pointer, contiguous :: iw => null()
  logical, dimension(:), pointer, contiguous :: quad => null()
  real(dp), dimension(:,:), pointer, contiguous :: lagrange => null()
  real(dp), dimension(:,:), pointer, contiguous :: pvt_xkap => null()
  !$omp threadprivate(active, fs, fe, iw, quad, lagrange, pvt_xkap)

contains

  !===============================================================================
  SUBROUTINE load_lam_interp_weights(self, lam)
    class(bftrans_type) :: self
    real(dp), dimension(:), intent(in) :: lam
    integer :: nf, f, di
    integer :: nx, j, i
    logical :: quadratic
    !---------------------------------------------------------------------------
    active = .true.

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
        active = .false.
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
        active = .false.
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
    quad(f) = quadratic

    do i=1, nf+nx ! i is not actually needed
      if (j >= nx) then
          iw(f) = nx

      elseif (self%lam(j+1) < lam(f)) then
          j = j+1
          quadratic = self%check_lam_quadratic(j)
          cycle

      elseif (self%lam(j) < lam(f)) then
          iw(f) = j

      else
          iw(f) = 1
      endif

      call self%set_xweight(lam, f)

      if (f == fe) exit
      f = f+di
      quad(f) = quadratic
    end do

    i = fs
    fs = min(fs, fe)
    fe = max(fe, i)
  END SUBROUTINE load_lam_interp_weights

    !===============================================================================
  SUBROUTINE load_unsorted_lam_interp_weights(self, lam)
    class(bftrans_type) :: self
    real(dp), dimension(:), intent(in) :: lam
    real(dp) :: minl, maxl
    integer :: nf, f
    integer :: nx, j, i
    logical :: quadratic
    !---------------------------------------------------------------------------
    active = .true.

    nx = size(self%lam)
    nf = size(lam)

    minl = minval(lam, 1)
    maxl = maxval(lam, 1)

    if (.not. self%minex .and. self%minl > maxl) then
      active = .false.
      return
    endif

    if (.not. self%maxex .and. self%maxl < minl) then
      active = .false.
      return
    endif

    fs = nf
    fe = 1
    iw(:) = 0

    do i=1, nf
      if (.not. self%maxex .and. lam(i) > self%maxl) cycle
      if (.not. self%minex .and. lam(i) < self%minl) cycle
      fs = min(i, fs)
      fe = max(i, fe)
      iw(i) = searchsorted(self%lam, lam(i), outbounds_ok=.true.)
      quad(i) = self%check_lam_quadratic(iw(i))
      call self%set_xweight(lam, i)
    end do
  END SUBROUTINE load_unsorted_lam_interp_weights

  !===============================================================================
  SUBROUTINE set_xweight(self, lam, f)
      class(bftrans_type) :: self
      real(dp), dimension(:), intent(in) :: lam
      real(dp) :: xx(3), x
      integer, intent(in) :: f
      integer :: j, k
      !---------------------------------------------------------------------------
      ! if there are not enough points, step back
      k = merge(2, 1, quad(f))
      iw(f) = min(iw(f), self%nnu-k)
      iw(f) = max(iw(f), 1)
      j = iw(f)
      k = k+1

      associate(xx => xx(1:k))
        x = lam(f)
        xx = self%lam(j:j+k-1)

        if (self%loglam) then
          x = log(x)
          xx = log(xx)
        endif

        if (quad(f)) then
            lagrange(1:k,f) = quadratic_weight(xx, x)
        else
            lagrange(1:k,f) = linear_weight(xx, x)
        endif
      end associate
  END SUBROUTINE set_xweight

  !===============================================================================
  FUNCTION quadratic_weight(xxx,newx) result(w)
      real(dp), dimension(3), intent(in) :: xxx
      real(dp), intent(in) :: newx
      real(dp) :: diff(3), w(3)
      !---------------------------------------------------------------------------
      diff = newx-xxx
      w(1) = diff(2) * diff(3) / ((xxx(1)-xxx(2)) * (xxx(1)-xxx(3)))
      w(2) = diff(1) * diff(3) / ((xxx(2)-xxx(1)) * (xxx(2)-xxx(3)))
      w(3) = diff(1) * diff(2) / ((xxx(3)-xxx(1)) * (xxx(3)-xxx(2)))
  end FUNCTION quadratic_weight

  !===============================================================================
  FUNCTION linear_weight(xx,newx) result(w)
      real(dp), dimension(2), intent(in) :: xx
      real(dp), intent(in) :: newx
      real(dp) :: w(2)
      !---------------------------------------------------------------------------
      w(1) = (xx(2) - newx) / (xx(2) - xx(1))
      w(2) = (newx - xx(1)) / (xx(2) - xx(1))
  end FUNCTION linear_weight

  !===============================================================================
  SUBROUTINE interp_xkap(self, chi, temp, teta, fakt, lam, ne, nhe, shield)
    class(bftrans_type) :: self
    class(level_type), pointer :: lvl
    real(4), dimension(:,:,:,:), intent(inout) :: chi
    real(dp), intent(in), dimension(:,:,:), target :: temp, teta
    real(dp), intent(inout), dimension(:,:,:) :: fakt
    real(dp), intent(in), dimension(:) :: lam
    real(dp), intent(in), dimension(:,:,:), optional :: ne, nhe, shield
    real(dp), dimension(:,:,:), pointer :: newt
    real(dp), dimension(:), allocatable :: ab
    real(dp), dimension(:,:), pointer :: xkap
    real(dp), dimension(3,self%ntt-1) :: dtdt
    real(dp), dimension(self%ntt-1) :: dtdt1
    real(dp) :: ab1, fakt1, t, diff(3), w(3), t3(3), dt, a, b, y1,y2,y3
    real(dp), dimension(self%nnu) :: yy
    real(dp), dimension(self%ntt) :: tt
    integer :: i, f, nx, ny, nz, x, y, z, it, it2, j1, j2, j3
    logical :: skip, do_log, do_exp, in_bounds, quadratic
    integer, save :: itimer1=0, itimer2=0, itimer3=0, itimer4=0, itimer5=0
    !---------------------------------------------------------------------------
    nx = size(temp, 1)
    ny = size(temp, 2)
    nz = size(temp, 3)
    allocate(ab(fs:fe), mold=lam(fs:fe))
    
    if (self%ntt > 1) then
      call trace%begin('interp_xkap_1', itimer=itimer1)
      do_exp = .false.
      do_log = .false.
      if (self%ilogt) then
          pvt_xkap = self%log_xkap
          do_exp = .not. self%ilogl
      else
          pvt_xkap = self%xkap
          do_log = self%ilogl
      endif

      if (self%iteta) then
          newt => teta
      else
          newt => temp
      endif

      if (self%kvadt) then
        do it=1, self%ntt-1
          if (it+2 <= self%ntt) then
            t3 = self%tt(it:it+2)
          else
            t3 = self%tt(it-1:it+1)
          endif
          dtdt(1,it) = 1. / ((t3(1)-t3(2)) * (t3(1)-t3(3)))
          dtdt(2,it) = 1. / ((t3(2)-t3(1)) * (t3(2)-t3(3)))
          dtdt(3,it) = 1. / ((t3(3)-t3(1)) * (t3(3)-t3(2)))
        end do
      endif
      
      do it=1, self%ntt-1
        dtdt1(it) = 1. / (self%tt(it+1) - self%tt(it))
      end do

      do z=1, nz
      do y=1, ny
      do x=1, nx
        t = temp(x,y,z)
        in_bounds = .true.
        quadratic = .false.
        if (.not. self%minet) in_bounds = t > self%mint
        if (.not. self%maxet) in_bounds = t < self%maxt

        if (in_bounds) then
          t = newt(x,y,z)
          if (t < self%tt(1)) then
            it = 1
          elseif (t > self%tt(self%ntt)) then
            it = self%ntt - 1
          else
            do it=1, self%ntt-1
              if (self%tt(it+1) > t) exit
            end do
            quadratic = self%kvadt
          endif

          if (quadratic) then
            ! quadratic interpolation
            if (it+2 <= self%ntt) then
              i = it
              xkap => pvt_xkap(:,it:it+2)
            else
              i = it-1
              xkap => pvt_xkap(:,it-1:it+1)
            endif
            diff = t - self%tt(i:i+2)
            w(1) = diff(2) * diff(3) * dtdt(1,it)
            w(2) = diff(1) * diff(3) * dtdt(2,it)
            w(3) = diff(1) * diff(2) * dtdt(3,it)
            yy(:) = 0
            yy(:) = yy(:) + xkap(:,1) * w(1)
            yy(:) = yy(:) + xkap(:,2) * w(2)
            yy(:) = yy(:) + xkap(:,3) * w(3)
          else
            ! linear interpolation
            dt = (t - self%tt(it)) * dtdt1(it)
            do f=1, self%nnu
              a = pvt_xkap(f,it)
              b = pvt_xkap(f,it+1)
              yy(f) = (b - a)*dt + a
            end do
          endif

          if (do_exp) then
            yy = exp(yy)
          elseif(do_log) then
            yy = log(yy)
          endif

          fakt1 = fakt(x,y,z)
          i = 0
          do f=fs,fe
            if (i /= iw(f)) then
              i = iw(f)
              y1 = yy(i)
              y2 = yy(i+1)
              if (quad(f)) y3 = yy(i+2)
            endif
            if (i /= 0) then
              ab1 =       y1*lagrange(1,f)
              ab1 = ab1 + y2*lagrange(2,f)
              if (quad(f)) then
                ab1 = ab1 + y3*lagrange(3,f)
              endif
              if (self%ilogl) ab1 = exp(ab1)
              chi(f,x,y,z) = chi(f,x,y,z) + fakt1*ab1
            endif
          end do
        endif
      end do
      end do
      end do
      call trace%end(itimer1)
    else
      call trace%begin('interp_xkap_2', itimer=itimer2)

      if(self%ilogl) then
        yy = self%log_xkap(:,1)
      else
        yy = self%xkap(:,1)
      endif

      do f=fs,fe
        i = iw(f)
        if (i==0) cycle
        if (quad(f)) then
            ab(f) = sum(yy(i:i+2) * lagrange(1:3,f))
        else
            ab(f) = sum(yy(i:i+1) * lagrange(1:2,f))
        endif
      end do
      if (self%ilogl) ab = exp(ab)

      if (H_occupation_probability .and. self%key == 'H BF') then
        lvl => self%ion%levels(self%k)
        lvl%i = self%k
        call trace%begin('wcalc2', itimer=itimer3)
        call wcalc2(lvl, chi, lam, ab, fakt, self%ion%n, nhe, ne, shield)
        call trace%end(itimer3)
      else
        do z=1, nz
        do y=1, ny
        do x=1, nx
          do f=fs,fe
            chi(f,x,y,z) = chi(f,x,y,z) + ab(f)*fakt(x,y,z)
          end do
        end do
        end do
        end do
      endif
      call trace%end(itimer2)
    endif
  END SUBROUTINE interp_xkap

  !===============================================================================
  PURE SUBROUTINE H_occupation(self, lam, nhe, ne, shield, temp, fakt, skip)
      class(bftrans_type), intent(in) :: self
      real(dp), intent(in) :: lam
      real(dp), intent(in), dimension(:,:,:) :: nhe, ne, shield, temp
      real(dp), intent(out), dimension(:,:,:) :: fakt
      logical, intent(out) :: skip
      real(dp) :: redcut, nu, x, scalar
      !---------------------------------------------------------------------------
      associate(lvl => self%ion%levels(self%k))
        redcut = hc_e / ((self%k**(-2.d0) - (self%k+1)**(-2.d0))*RYD) * 1d8
        if (lam > redcut) then
          redcut = exp((redcut - lam)*1.d-2)
        else
          redcut = 1.0
        endif
        skip = redcut < 1d-1
        if (.not. skip) then
          nu = cc/lam*1e8
          scalar = 2.0898E-26 / (2*self%k**5) * redcut

          x = 1. - (lvl%xion + h_e*nu)/RYD
          if(x > 0.d0) then
            x = 1.d0/sqrt(x)
            call wcalc(self%ion%n, nhe, ne, shield, x, fakt)
            fakt = lvl%w - fakt
          else
            fakt = lvl%w
          endif

          fakt = fakt * scalar * lvl%n
        endif
      end associate
  END SUBROUTINE H_occupation

  !===============================================================================
  elemental FUNCTION interp_temp_quadratic(self, f, it, newt) result(newy)
      class(bftrans_type), intent(in) :: self
      real(dp), intent(in) :: newt
      integer, intent(in) :: f, it
      real(dp) :: newy, w(3), diff(3)
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

      diff = newt - self%tt(j1:j3)
      w(1) = diff(2) * diff(3) / ((self%tt(j1)-self%tt(j2)) * (self%tt(j1)-self%tt(j3)))
      w(2) = diff(1) * diff(3) / ((self%tt(j2)-self%tt(j1)) * (self%tt(j2)-self%tt(j3)))
      w(3) = diff(1) * diff(2) / ((self%tt(j3)-self%tt(j1)) * (self%tt(j3)-self%tt(j2)))

      newy = sum(pvt_xkap(f,j1:j3) * w(1:3))
  end FUNCTION interp_temp_quadratic 

  !===============================================================================
  elemental FUNCTION interp_temp_linear(self, f, it, newt) result(newy)
      class(bftrans_type), intent(in) :: self
      real(dp), intent(in) :: newt
      integer, intent(in) :: f, it
      real(dp) :: m, newy
      !---------------------------------------------------------------------------
      m = (pvt_xkap(f,it+1) - pvt_xkap(f,it)) / (self%tt(it+1) - self%tt(it))
      newy = m * (newt - self%tt(it)) + pvt_xkap(f,it)
  end FUNCTION interp_temp_linear

  !===============================================================================
  elemental FUNCTION check_lam_quadratic(self, iw) result(quadratic)
      class(bftrans_type), intent(in) :: self
      integer, intent(in) :: iw
      logical :: quadratic
      !---------------------------------------------------------------------------
      quadratic = .false.

      if (self%kvadl) then
        if (self%nnu - iw > 1) then
          if (iw <= 0) then
            quadratic = .false.
          elseif (.not. self%ldisk(iw)) then
            quadratic = .true.
          endif
        endif
      endif
  end FUNCTION check_lam_quadratic

  !===============================================================================
  elemental FUNCTION check_temp_quadratic(self, it) result(quadratic)
      class(bftrans_type), intent(in) :: self
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
  end FUNCTION check_temp_quadratic

  !===============================================================================
  elemental SUBROUTINE get_temp_index(self, newt, idx, no_extrap)
      class(bftrans_type), intent(in) :: self
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
  END SUBROUTINE get_temp_index

  !===============================================================================
  SUBROUTINE absorption_factor(self, EK_over_temp, teta, pe, dxi, fakt)
      class(bftrans_type), intent(in) :: self
      class(particle_type), pointer :: ion
      real(dp), dimension(:,:,:), intent(in) :: EK_over_temp, teta, pe, dxi
      real(dp), dimension(:,:,:), intent(out) :: fakt
      real(dp), dimension(:,:,:), allocatable :: h1
      integer, dimension(:,:,:), allocatable :: nniv
      real(dp) :: scalar
      !---------------------------------------------------------------------------
      ion => self%ion
      !if (self%key /= 'H- BF') then
      !  fakt = 0
      !  return
      !endif
      select case(trim(self%key))
      case('H- BF')
          fakt = ion%n
      case('H- FF')
          fakt = ion%n * pe
      case('H BF')
          if (.not. H_occupation_probability) then
            scalar = 2.0898E-26 / (2*self%k**5)
            fakt = scalar * ion%levels(self%k)%n
          endif
      case('H FF')
          !fakt = 2.0898E-26 / 2 * ion%part * ion%n * exp(-31.30364*teta + 2.3026*dxi*teta) / (31.30364*teta)
          fakt = ion%boltzmann(EK_over_temp, dxi, natural_log=.true.)
          fakt = 2.0898E-26 / (2*RYD*EK_over_temp) * fakt * ion%part * ion%n
      case('(H+H)')
          ! see Doyle 1968 "The astrophysical significance of the continuous spectrum of the hydrogen quasi-molecule"
          fakt = (1d-25 * ion%n)**2
      case('H2+')
          select type(ion)
          class is (ion_type)
              fakt = 1d-40 * ion%n * ion%above%n
          end select
      case default
        if (self%pedep) then
          fakt = pe * ion%n
          if (self%free_free) then
            fakt = (1d-40 / EE) * EK_over_temp * fakt
          endif
        elseif (self%molecular) then
            fakt = ion%n * ion%awgt / 6.023d23
        elseif (self%regular) then
            fakt = self%level%g * self%level%boltzmann(EK_over_temp, natural_log=.true.)
            fakt = fakt * ion%part * ion%n
        endif
      end select
  END SUBROUTINE absorption_factor

  !===============================================================================
  SUBROUTINE add_absorption(self, chi, lam, temp, EK_over_temp, teta, pe, dxi, ne, nhe, shield, sorted)
      class(bftrans_type), intent(in) :: self
      real(4), dimension(:,:,:,:), intent(inout) :: chi
      real(dp), dimension(:,:,:), intent(in) :: temp, EK_over_temp, teta, pe, dxi
      real(dp), dimension(:,:,:), intent(in), optional :: ne, nhe, shield
      real(dp), dimension(:), intent(in) :: lam
      real(dp), dimension(:,:,:), allocatable :: fakt
      logical, optional :: sorted
      real(dp) :: prev
      integer :: i,j,k,nf, m(2), p(4)=pp
      !---------------------------------------------------------------------------
      if (.not. self%active) return
      nf = size(lam)
      call os%allocate(quad, nf)
      call os%allocate(iw, nf)
      call os%allocate(lagrange, [3,nf], zero=.true.)
      call os%allocate(pvt_xkap, shape(self%xkap))
      allocate(fakt, mold=temp)

      if (option(sorted, .true.)) then
        call self%load_lam_interp_weights(lam)
      else
        call self%load_unsorted_lam_interp_weights(lam)
      endif

      if (active) then
          !if (nf > 1) prev = chi(p(1), p(2), p(3), p(4))
          call self%absorption_factor(EK_over_temp, teta, pe, dxi, fakt)
          call self%interp_xkap(chi, temp, teta, fakt, lam, ne, nhe, shield)
          !if (nf > 1) print '(A8, 2x, A12, 2x, A, 8(ES14.6))', self%key, self%label, self%ion%id, &
          ! fakt(p(1),p(2),p(3)), chi(p(1),p(2),p(3),p(4)) - prev, chi(p(1),p(2),p(3),p(4))
      endif

      call os%deallocate(quad)
      call os%deallocate(iw)
      call os%deallocate(lagrange)
      call os%deallocate(pvt_xkap)
  END SUBROUTINE add_absorption


END MODULE
