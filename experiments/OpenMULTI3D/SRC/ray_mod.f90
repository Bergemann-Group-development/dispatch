MODULE ray_mod
  USE const_mod
  USE radau_mod
  USE console_mod
  implicit none

  ! *******************************************************************

  type ray_type
    ! cubic-convolution weights in x and y for the upwind intersection with the
    ! first XY plane, and the offset of the first of the four points relative to
    ! the middle point of the short characteristic
    real(8), dimension(4) :: xkern, ykern
    integer :: upo1, upo2
    ! ray length to the upwind XY-plane intersection
    real(8) :: dsup
  contains
    procedure :: set_stencil
  end type ray_type

  logical, public :: one_dimensional=.false.

contains

  !===============================================================================
  SUBROUTINE set_carlson_quadratures (quad_scheme, vec, wts, Jnu_weights)
    !---------------------------------------------------------------------
    !  Routine is modified from set_ang by AJG for output ray quadrature
    !  schemes
    !---------------------------------------------------------------------
    character(len=*), intent(in) :: quad_scheme
    integer :: m, i, nro, nmu
    real(8), dimension(:), allocatable :: mux,  muy,  wmu
    real(8), dimension(:), allocatable, intent(out) :: wts
    real(8), dimension(:,:), allocatable, intent(out) :: vec
    logical, optional :: Jnu_weights
    logical :: mirror=.true.
    ! ----------------------------------------------------------------------------
    if (present(Jnu_weights)) then
      mirror = Jnu_weights
    endif

    ! define the different angle sets
    select case (trim(quad_scheme))

     case('disk_center')
      if (mirror) then
        allocate(vec(3, 2), wts(2))
        vec(:,1) = [0.0, 0.0, 1.0]
        vec(:,2) = [0.0, 0.0,-1.0]
        wts = 0.5
      else
        allocate(vec(3, 1), wts(1))
        vec(:,1) = [0.0, 0.0, 1.0]
        wts = 1.0
      endif
      return

     case('set_a2')
      nro    = 1
      allocate(mux(nro), muy(nro), wmu(nro))
      mux    = (/0.57735026/)
      muy    = (/0.57735026/)
      wmu    = (/1.00000000/)

     case('set_a4')
      nro = 3
      allocate(mux(nro), muy(nro), wmu(nro))
      mux = (/0.88191710, 0.33333333, 0.33333333/)
      muy = (/0.33333333, 0.88191710, 0.33333333/)
      wmu = (/0.33333333, 0.33333333, 0.33333333/)

     case('set_a6')
      nro = 6
      allocate(mux(nro), muy(nro), wmu(nro))
      mux = (/0.93094934, 0.68313005, 0.25819889,&
        0.68313005, 0.25819889, 0.25819889/)
      muy = (/0.25819889, 0.68313005, 0.93094934,&
        0.25819889, 0.68313005, 0.25819889/)
      wmu = (/0.18333333, 0.15000000, 0.18333333,&
        0.15000000, 0.15000000, 0.18333333/)

     case('set_a8')
      nro = 10
      allocate(mux(nro), muy(nro), wmu(nro))
      mux = (/0.95118973, 0.78679579, 0.57735027,&
        0.21821789, 0.78679579, 0.57735027,&
        0.21821789, 0.57735027, 0.21821789, 0.21821789/)
      muy = (/0.21821789, 0.57735027, 0.78679579,&
        0.95118973, 0.21821789, 0.57735027,&
        0.78679579, 0.21821789, 0.57735027, 0.21821789/)
      wmu = (/0.12698138, 0.09138353, 0.09138353,&
        0.12698138, 0.09138353, 0.07075469,&
        0.09138353, 0.09138353, 0.09138353, 0.12698138/)

     case('set_b4')
      nro = 3
      allocate(mux(nro), muy(nro), wmu(nro))
      mux = (/0.70412415, 0.09175171, 0.09175171/)
      muy = (/0.09175171, 0.70412415, 0.09175171/)
      wmu = (/0.33333333, 0.33333333, 0.33333333/)

     case('set_b6')
      nro = 6
      allocate(mux(nro), muy(nro), wmu(nro))
      mux = (/0.80847426, 0.57735027, 0.11417547,&
        0.57735027, 0.11417547, 0.11417547/)
      muy = (/0.11417547, 0.57735027, 0.80847426,&
        0.11417547, 0.57735027, 0.11417547/)
      wmu = (/0.22222222, 0.11111111, 0.22222222,&
        0.11111111, 0.11111111, 0.22222222/)

     case('set_b8')
      nro = 10
      allocate(mux(nro), muy(nro), wmu(nro))
      mux = (/0.85708018, 0.70273364, 0.50307327,&
        0.11104440, 0.70273364, 0.50307327,&
        0.11104440, 0.50307327, 0.11104440, 0.11104440/)
      muy = (/0.11104440, 0.50307327, 0.70273364,&
        0.85708018, 0.11104440, 0.50307327,&
        0.70273364, 0.11104440, 0.50307327, 0.11104440/)
      wmu = (/0.17152232, 0.06900520, 0.06900520,&
        0.17152232, 0.06900520, 0.07140185,&
        0.06900520, 0.06900520, 0.06900520, 0.17152232/)

     case default
      call console_error ('quadrature keyword "'// trim(quad_scheme) // '" not recognized')
      error stop
    end select

    if (mirror) then
      nmu = 8*nro
    else
      nmu = 4*nro
    endif

    allocate(vec(3, nmu), wts(nmu))
    vec = 0
    wts = 0

    vec(1, :nro) = -mux
    vec(2, :nro) = -muy
    vec(3, :nro) = sqrt(1.0d0 - mux**2 - muy**2)
    wts(:nro) = wmu

    ! Fill the second octant with rays reflected in the y-z plane
    wts(  nro+1:2*nro) =  wts(  1:nro)
    vec(1,nro+1:2*nro) = -vec(1,1:nro)
    vec(2,nro+1:2*nro) =  vec(2,1:nro)
    vec(3,nro+1:2*nro) =  vec(3,1:nro)

    ! Fill the third and fourth octant with rays reflected in the x-z plane
    wts(  2*nro+1:4*nro) =  wts(  1:2*nro)
    vec(1,2*nro+1:4*nro) = -vec(1,1:2*nro)
    vec(2,2*nro+1:4*nro) = -vec(2,1:2*nro)
    vec(3,2*nro+1:4*nro) =  vec(3,1:2*nro)

    ! Normalize, so that integration over 8 octants yields 1
    ! the above sets are normalized so that integration over one octant yields 1.

    if (mirror) then
      ! Fill the other octants with inverted rays
      wts(  4*nro+1:8*nro) =  wts(  1:4*nro)
      vec(:,4*nro+1:8*nro) = -vec(:,1:4*nro)
      wts = wts / 8
    else
      wts = wts / 4
      wts = wts * vec(:,3)
    endif
  END SUBROUTINE set_carlson_quadratures
!---------------------------------------------------------------------

  !===============================================================================
  SUBROUTINE set_grid_quadratures (quad_scheme, ds, vec, wts)
    character(len=*), intent(in) :: quad_scheme
    real(8), intent(in) :: ds(3)
    real(8), dimension(:), allocatable :: mu, zrot
    real(8), dimension(:), allocatable, intent(out) :: wts
    real(8), dimension(:,:), allocatable, intent(out) :: vec
    real(8) :: plus, cross
    ! ----------------------------------------------------------------------------
    plus = ds(1)/ds(3)
    cross = sqrt(ds(1)**2 + ds(2)**2)/ds(3)

    ! define the different angle sets
    select case (trim(quad_scheme))
     case('grid+')
      allocate(mu, source=[0.d0, atan(plus)])
      allocate(zrot, source=[0.d0, 0.d0])
     case('gridx')
      allocate(mu, source=[0.d0, atan(cross)])
      allocate(zrot, source=[0.d0, PI/4.])
     case('grid*')
      allocate(mu, source=[0.d0, atan(plus), atan(cross)])
      allocate(zrot, source=[0.d0, 0.d0, PI/4.])
    end select

    mu = cos(mu)
    call set_custom_quadratures(mu, nphi=4, vec=vec, wts=wts, zrot=zrot)
  END SUBROUTINE set_grid_quadratures

  !===============================================================================
  SUBROUTINE set_lobatto_quadratures(ntheta, nphi, vec, wts, Jnu_weights)
    ! -------------------------------------------------------------------
    ! Computes Lobatto quadratures for a given number of theta and phi
    ! angles. It then converts the spherical geometry into cartesian geometry
    ! -------------------------------------------------------------------
    integer, intent(in) :: ntheta, nphi
    real(8), intent(out), dimension(:), allocatable :: wts
    real(8), intent(out), dimension(:,:), allocatable :: vec
    real(8), dimension(:), allocatable :: mu, wmu
    logical, optional :: Jnu_weights
    ! -------------------------------------------------------------------
    allocate(mu(ntheta))
    allocate(wmu(ntheta))

    ! setting indexes and hard-coded Lobatto coefficients
    select case (ntheta)
     case (1)
      mu  = 1.0
      wts = 1.0

     case (2)
      mu  = (/1.0,        0.44721360/)
      wmu = (/0.16666667, 0.83333333/)

     case (3)
      mu  = (/1.0,        0.76505532, 0.28523152/)
      wmu = (/0.06666667, 0.37847496, 0.55485838/)

     case (4)
      mu  = (/1.0,        0.87174015, 0.59170018, 0.20929922/)
      wmu = (/0.03571428, 0.21070422, 0.34112270, 0.41245880/)

     case (5)
      mu  = (/1.0,          0.9195339082, 0.7387738651, &
        0.4779249498, 0.1652789577/)
      wmu = (/0.0222222222, 0.1333059908, 0.2248893420, &
        0.2920426836, 0.3275397612/)

     case (6)
      mu  = (/1.0           , 0.944899272223, 0.819279321644, &
        0.632876153032, 0.399530940965, 0.136552932855/)
      wmu = (/0.015151515153, 0.091684517440, 0.157974705655, &
        0.212508417834, 0.251275603210, 0.271405240910/)

     case default
      error stop 'Error: ntheta is set too high. *ntheta<=5*'
    end select

    call nphi_copies(nphi, mu, wmu, vec, wts, Jnu_weights=Jnu_weights)
  END SUBROUTINE set_lobatto_quadratures

  !===============================================================================
  SUBROUTINE set_gauss_quadratures(ntheta, nphi, vec, wts, Jnu_weights)
    ! -------------------------------------------------------------------
    ! Computes Gauss-Legendre quadratures for a given number of theta and phi
    ! angles. It then converts the spherical geometry into cartesian geometry
    ! -------------------------------------------------------------------
    integer, intent(in) :: ntheta, nphi
    real(8), intent(out), dimension(:), allocatable :: wts
    real(8), intent(out), dimension(:,:), allocatable :: vec
    real(8), dimension(:), allocatable :: mu, wmu
    logical, optional :: Jnu_weights
    ! -------------------------------------------------------------------
    allocate(mu(ntheta))
    allocate(wmu(ntheta))

    ! setting indexes and hard-coded Lobatto coefficients
    select case (ntheta)
     case (1)
      mu  = 1.0
      wts = 1.0

     case (2)
      mu  = (/0.78867513, 0.21132487/)
      wmu = (/0.50000000, 0.50000000/)

     case (3)
      mu  = (/0.88729833, 0.50000000, 0.11270167/)
      wmu = (/0.27777778, 0.44444444, 0.27777778/)

     case (4)
      mu  = (/0.93056816, 0.66999052, 0.33000948, 0.06943184/)
      wmu = (/0.17392742, 0.32607258, 0.32607258, 0.17392742/)

     case (5)
      mu  = (/0.95308992, 0.76923466, 0.50000000, 0.23076534, 0.04691008/)
      wmu = (/0.11846344, 0.23931434, 0.28444444, 0.23931434, 0.11846344/)

     case (6)
      mu  = (/0.96623476, 0.83060469, 0.61930959, 0.38069041, 0.16939531, 0.03376524/)
      wmu = (/0.08566225, 0.18038079, 0.23395697, 0.23395697, 0.18038079, 0.08566225/)

     case default
      error stop 'Error: ntheta is set too high. *ntheta<=6*'
    end select

    call nphi_copies(nphi, mu, wmu, vec, wts, Jnu_weights=Jnu_weights)
  END SUBROUTINE set_gauss_quadratures

  !===============================================================================
  SUBROUTINE set_custom_quadratures(mu, nphi, vec, wts, zrot)
    real(8), intent(inout), dimension(:) :: mu
    integer, intent(in) :: nphi
    real(8), intent(out), dimension(:), allocatable :: wts
    real(8), intent(out), dimension(:,:), allocatable :: vec
    real(8), dimension(:), allocatable :: wmu, smus, swts, swts2
    integer, dimension(:), allocatable :: order
    real(8), dimension(:), optional :: zrot
    integer :: nrays
    ! -------------------------------------------------------------------
    nrays = size(mu)
    allocate(wmu(nrays))
    allocate(order(nrays))
    allocate(smus(nrays+1))
    allocate(swts(nrays+1))

    if (nrays==1) then
      wmu = 1
    else
      call sort_array(mu, mu, order)
      smus(1) = 0.d0
      smus(2:nrays+1) = mu

      if (mod(size(smus), 2) == 1) then
        call simps_weights(smus, swts)
      else
        allocate(swts2(nrays))
        call simps_weights(smus(1:nrays), swts(1:nrays))
        call simps_weights(smus(2:nrays+1), swts2)
        swts(1:2) = swts(1:2) + 0.5*(smus(2)-smus(1))
        swts2(nrays-1:nrays) = swts2(nrays-1:nrays) + 0.5*(smus(nrays+1)-smus(nrays))
        swts(2:nrays+1) = swts(2:nrays+1) + swts2
        swts = swts / 2
      endif

      wmu = swts(nrays+1:2:-1)
      mu = mu(nrays:1:-1)
    endif

    if (present(zrot)) zrot = zrot(order(nrays:1:-1))
    call nphi_copies(nphi, mu, wmu, vec, wts, zrot, Jnu_weights=.false.)
  END SUBROUTINE set_custom_quadratures

  !===============================================================================
  SUBROUTINE set_radau_quadratures(ntheta, nphi, vec, wts, Jnu_weights)
    integer, intent(in) :: ntheta, nphi
    real(8), intent(out), dimension(:), allocatable :: wts
    real(8), intent(out), dimension(:,:), allocatable :: vec
    logical, intent(in), optional :: Jnu_weights
    real(4), dimension(:), allocatable :: wmu2, mu2
    real(8), dimension(:), allocatable :: wmu, mu
    real(8), dimension(:), allocatable :: zrot
    ! -------------------------------------------------------------------
    allocate(mu(ntheta), wmu(ntheta))
    allocate(mu2(ntheta), wmu2(ntheta))
    allocate(zrot(ntheta))
    call radau%init(ntheta, mu2, wmu2)
    mu = mu2(ntheta:1:-1)
    wmu = wmu2(ntheta:1:-1)
    zrot = PI/4.
    call nphi_copies(nphi, mu, wmu, vec, wts, zrot, Jnu_weights=Jnu_weights)
  END SUBROUTINE set_radau_quadratures

  !===============================================================================
  SUBROUTINE nphi_copies (nphi_in, mu, wmu, vec, wts, zrot, Jnu_weights)
    integer, intent(in) :: nphi_in
    real(8), intent(inout), dimension(:) :: mu, wmu
    real(8), intent(out), dimension(:), allocatable :: wts
    real(8), intent(out), dimension(:,:), allocatable :: vec
    real(8), intent(in), optional, dimension(:) :: zrot
    real(8) :: xphi, theta, xoff
    integer :: nrays, ntheta, i, nt, np, nmuphi, n1, nphi
    logical :: includes_mu1, mirror=.false.
    logical, optional :: Jnu_weights
    ! -------------------------------------------------------------------
    if (present(Jnu_weights)) then
      mirror = Jnu_weights
    endif

    nrays = size(mu)
    includes_mu1 = findloc(mu, 1.0, 1) > 0

    ! no need for nphi copies if atmosphere is plane parallel
    nphi = merge(1, nphi_in, one_dimensional)

    if (includes_mu1) then
      nmuphi = (nrays-1)*nphi + 1
    else
      nmuphi = nrays*nphi
    endif

    if (mirror) nmuphi = 2*nmuphi
    allocate(wts(nmuphi))
    allocate(vec(3, nmuphi))

    ! storing first index of vec and wts
    if (includes_mu1) then
      wts(1) = wmu(1)
      vec(:,1) = [0.,0.,1.]
      n1 = 2
    else
      n1 = 1
    endif

    ! storing the rest of vec and wts indexes
    do nt = n1, nrays
      do np = 1, nphi
        i = (nt-n1) * nphi + np
        if (includes_mu1) i = i+1

        if (present(zrot)) then
          xoff = zrot(nt) ! zrot assumed in radians
        else
          xoff = 0.d0
        endif

        ! Computing cartesian geometries
        xphi = 2.d0*PI / dble(nphi) * (np-1) + xoff
        theta = acos(mu(nt))
        vec(1,i) = sin(theta) * cos(xphi)
        vec(2,i) = sin(theta) * sin(xphi)
        vec(3,i) = mu(nt)
        wts(i)   = wmu(nt) / dble(nphi)
      end do
    end do

    ! Remove computational anomalies
    where (abs(vec) < 1d-3) vec = 0d0

    if (mirror) then
      wts(  nmuphi/2+1:nmuphi) =  wts(1:nmuphi/2)
      vec(:,nmuphi/2+1:nmuphi) = -vec(:,1:nmuphi/2)
      wts = wts / 2
    else
      ! Flux weights
      wts = wts * vec(3,:)
    endif
  END SUBROUTINE nphi_copies

!---------------------------------------------------------------------
  SUBROUTINE set_stencil(ray, ds, vec)
    !---------------------------------------------------------------------
    !     Calculates the interpolation weights, the physical length of the
    !     ray, and the side the ray hits.
    !
    !     Positive angles point in the direction of increasing index
    !
    !     The upwind point is the point where the ray intersects
    !     the first grid plane, diagram below is for positive mux and muz
    !
    !                        downwind
    !                        xdn
    !z9   -------------------0----------
    !     |              |  /          |
    !     |              | /           |
    !     |              |/            |
    !z8   ---------------0--------------
    !     |            / |             |
    !     |           /  |             |
    !     |          /   |             |
    !z7   ----------0-------------------
    !     x3        xup  x4            x5
    !               upwind
    !
    !---------------------------------------------------------------------
    class(ray_type) :: ray
    real(8), intent(in) :: ds(3), vec(3)
    real(8) :: wt1, wt2, h1, h2
    !---------------------------------------------------------------------
    ray%dsup = ds(3) / abs(vec(3))
    wt1 = -vec(1) / ds(1) * ray%dsup
    wt2 = -vec(2) / ds(2) * ray%dsup
    h1 = merge(wt1, wt1+1, wt1>0)
    h2 = merge(wt2, wt2+1, wt2>0)
    ray%upo1 = merge(-1, -2, wt1>0)
    ray%upo2 = merge(-1, -2, wt2>0)
    call cc_kernel(h1,ray%xkern)
    call cc_kernel(h2,ray%ykern)
  END SUBROUTINE set_stencil
!---------------------------------------------------------------------


!===============================================================================
! Returns simpson integration weights for odd numbered xx-grid.
  !===============================================================================
  SUBROUTINE simps_weights(xx, wts)
    real(8), intent(in), dimension(:) :: xx
    real(8), intent(out), dimension(:) :: wts
    real(8), allocatable, target, dimension(:) :: h
    real(8), pointer, dimension(:) :: h0, h1
    integer :: n
    !-----------------------------------------------------------------------------
    n = size(xx)

    if (n == 1) then
      wts = 1.
      return
    elseif (n < 3) then
      error stop 'less than 3 weights in simps_weights'
    endif

    allocate(h(n-1))

    wts = 0.d0
    h = xx(2:n) - xx(1:n-1)
    h0 => h(1:n-2:2)
    h1 => h(2:n-1:2)

    wts(1:n-2:2) = wts(1:n-2:2) +  (h0 + h1)/6 * (2 - h1 / h0)
    wts(2:n-1:2) = wts(2:n-1:2) +  (h0 + h1)**3 / (6*h0*h1)
    wts(3:n-0:2) = wts(3:n-0:2) +  (h0 + h1)/6 * (2 - h0 / h1)
  END SUBROUTINE simps_weights

  !===============================================================================
  SUBROUTINE sort_array(array, sorted_array, permutation_order)
    real(kind=8), intent(in) :: array(:)
    real(kind=8), intent(out) :: sorted_array(:)
    integer, intent(out) :: permutation_order(:)
    real(kind=8), allocatable :: temp_array(:)
    integer :: n, i, j
    !-----------------------------------------------------------------------------
    n = size(array)

    ! Copy the original array to a temporary array
    allocate(temp_array, source=array)

    ! Initialize the permutation order
    permutation_order = [(i, i = 1, n)]

    ! Sort the temporary array and permutation order using bubble sort
    do i = 1, n - 1
      do j = 1, n - i
        if (temp_array(j) > temp_array(j+1)) then
          ! Swap elements in the temporary array and permutation order
          temp_array(j:j+1) = [temp_array(j+1), temp_array(j)]
          permutation_order(j:j+1) = [permutation_order(j+1), permutation_order(j)]
        end if
      end do
    end do

    ! Return the sorted array and permutation order
    sorted_array = temp_array
  END SUBROUTINE sort_array

  !---------------------------------------------------------------------
  SUBROUTINE cc_kernel(s,k)
    !---------------------------------------------------------------------
    !
    !  Computes the 1D cubic convolution kernel. For an equidistant
    !  grid x_i (i={1,2,3,4}) with spacing h, a coordinate x_2<x<x_3
    !  where the function value is desired, and s= (x-x_2)/h, the
    !  kernel k is computed, so that y(x) approx sigma_1^4 k[i]*y_i
    !  for tabulated values y_i at position x_i
    !
    !---------------------------------------------------------------------
    implicit none

    real(8),intent(in)               :: s
    real(8),intent(out),dimension(4) :: k
    real(8)                          :: s2,s3

    s2=s*s
    s3=s2*s

    k(1) = 0.5*(-s3+2*s2-s)
    k(2) = 0.5*( 3*s3-5*s2+2)
    k(3) = 0.5*(-3*s3+4*s2+s)
    k(4) = 0.5*(s3-s2)

  END SUBROUTINE cc_kernel
  !---------------------------------------------------------------------

END MODULE ray_mod
