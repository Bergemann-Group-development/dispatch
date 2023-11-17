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
      procedure :: interp_xkap
      procedure :: check_lam_quadratic
      procedure :: absorption_factor
      procedure :: add_absorption
      procedure :: set_xweight
      procedure :: load_lam_interp_weights
  end type

  logical :: active
  integer :: fs, fe
  logical, public :: H_occupation_probability=.true.
  integer, dimension(:), allocatable :: iw
  logical, dimension(:), allocatable :: quad
  real(dp), dimension(:,:), allocatable :: lagrange
  real(dp), dimension(:,:), pointer, contiguous :: pvt_xkap => null()
  !$omp threadprivate(active, fs, fe, iw, quad, lagrange, pvt_xkap)

contains

    !===============================================================================
  SUBROUTINE load_lam_interp_weights(self, lam)
    class(bftrans_type) :: self
    real(dp), dimension(:), intent(in) :: lam
    real(dp) :: minl, maxl
    integer :: nf, f
    integer :: nx, j, i
    logical :: quadratic, above, below, skip
    integer, save :: itimer=0
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
    call trace%begin('load_lam_interp_weights', itimer=itimer)

    fs = nf
    fe = 1
    iw(:) = 0

    do i=1, nf
      above = lam(i) > self%maxl
      below = lam(i) < self%minl
      skip = (.not. self%maxex .and. above) .or. (.not. self%minex .and. below)
      if (.not. skip) then
        fs = min(i, fs)
        fe = max(i, fe)
        j = -1
        if (below) then
          j = 0
        elseif (above) then
          j = self%nnu
        elseif (i > 1) then
          if (fits_in(iw(i-1))) then
            j = iw(i-1)
          elseif (fits_in(iw(i-1)-1)) then
            j = iw(i-1)-1
          endif
        endif
        if (j<0) then
          do j=1, self%nnu-1
            if (fits_in(j)) exit
          end do
        endif
        iw(i) = j
        quad(i) = self%check_lam_quadratic(j)
        call self%set_xweight(lam, i)
      endif
    end do
    call trace%end(itimer)
  contains
    pure logical function fits_in(idx)
      integer, intent(in) :: idx
      fits_in = idx > 0 .and. idx < self%nnu
      if (fits_in) then
        fits_in = lam(i) > self%lam(idx) .and. lam(i) < self%lam(idx+1)
      endif
    end function fits_in
  END SUBROUTINE load_lam_interp_weights

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
    real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
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

          i = 0
          ab = 0
          do f=fs,fe
            if (i /= iw(f)) then
              i = iw(f)
              y1 = yy(i)
              y2 = yy(i+1)
              if (quad(f)) then
                y3 = yy(i+2)
              else
                y3 = 0
              endif
            endif
            if (i /= 0) then
              ab(f) = sum([y1, y2, y3]*lagrange(:,f))
            endif
          end do
          fakt1 = fakt(x,y,z)
          if (self%ilogl) ab = exp(ab)
          chi(fs:fe,x,y,z) = chi(fs:fe,x,y,z) + fakt1*ab(fs:fe)
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
        call wcalc2(lvl, chi, lam(fs:fe), ab, fakt, self%ion%n, nhe, ne, shield)
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
      real(4), dimension(:,:,:,:), contiguous, intent(inout) :: chi
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
      allocate(quad(nf))
      allocate(iw(nf))
      allocate(lagrange(3,nf), source=0.d0)
      call os%allocate(pvt_xkap, shape(self%xkap))
      allocate(fakt, mold=temp)

      call self%load_lam_interp_weights(lam)

      if (active) then
          !if (nf > 1) prev = chi(p(1), p(2), p(3), p(4))
          call self%absorption_factor(EK_over_temp, teta, pe, dxi, fakt)
          call self%interp_xkap(chi, temp, teta, fakt, lam, ne, nhe, shield)
          !if (nf > 1) print '(A8, 2x, A12, 2x, A, 8(ES14.6))', self%key, self%label, self%ion%id, &
          ! fakt(p(1),p(2),p(3)), chi(p(1),p(2),p(3),p(4)) - prev, chi(p(1),p(2),p(3),p(4))
      endif

      deallocate(quad)
      deallocate(iw)
      deallocate(lagrange)
      call os%deallocate(pvt_xkap)
  END SUBROUTINE add_absorption


END MODULE
