MODULE ray_mod
  USE const_mod
  USE radau_mod
  implicit none

  ! *******************************************************************
 
  type intersect_type
    ! intersection plane 1:XZ;2:YZ;
    integer :: plane
    ! linear interpolation weights for four points surrounding the intersection
    real(8),dimension(2,2) :: wght
    ! offsets of the four points relative to where the ray originates from
    integer,dimension(2,2) :: iof,jof,kof
    ! distance to previous plane intersection (the last ds is stored in the ray struct)
    real(8) :: ds 
  end type intersect_type
 
! *******************************************************************

  type ray_type

    ! ange cosines
    real(8) :: mux,muy,muz
    ! quadrature weights
    real(8) :: wmu

    ! iu is the ofset of the upwind point, id offset of the downwind point
    integer :: iup,idn
    integer :: jup,jdn
    integer :: kup,kdn

    ! interpolation weights for the cubic convolution in the XY plane
    ! for the upwind point upwind and the first intersection for the
    ! downwind point
    real(8),dimension(4,4) :: upkern,dnkern
    ! offset of the point closest to the origin of the 16 points
    ! surrounding the first plane intersection, relative to the middle
    ! point of the short characteristic
    integer :: upo1,upo2,dno1,dno2

    ! distance to upwind XY-plane intersection in case ray is long
    real(8) :: longds
    ! interpolation weights for the cubic convolution in the upwind
    ! XY plane in case the ray is long
    real(8),dimension(4,4) :: longkern
    ! offset of the point closest to the origin of the 16 points
    ! surrounding this plane intersection, relative to the middle
    ! point of the short characteristic
    integer :: longo1,longo2

    ! which plane is hit first by downwind characteristic
    ! 1: YZ; 2:XZ; 3:XY
    integer :: dnplane
    integer :: upplane

    ! ray length to next XY plane intersection for the upwind point
    ! and next intersection for the downwind point
    real(8) :: dsup,dsdn

    ! The number of intersections of this ray with cell walls
    integer :: nintersect
    ! is this a long ray
    logical :: long
    ! Indicates that a ray is cut off at the boundary
    logical :: cut

    type(intersect_type),dimension(8) :: int

  contains
    procedure :: init => ray_init
    procedure :: set_stencil

  end type ray_type
  
  logical, public :: one_dimensional=.false.
  real(4), save :: rstar=0, rorigin=0
  
contains

!===============================================================================
SUBROUTINE set_carlson_quadratures (quad_scheme, vec, wts, up_and_down)
  !---------------------------------------------------------------------
  !  Routine is modified from set_ang by AJG for output ray quadrature
  !  schemes
  !---------------------------------------------------------------------
  character(len=*), intent(in) :: quad_scheme
  integer :: m, i, nro, nmu
  real(8), dimension(:), allocatable :: mux,  muy,  wmu
  real(8), dimension(:), allocatable, intent(out) :: wts
  real(8), dimension(:,:), allocatable, intent(out) :: vec
  logical, optional :: up_and_down
  logical :: mirror=.true.
  ! ----------------------------------------------------------------------------
  if (present(up_and_down)) then
    mirror = up_and_down
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
    print *, 'Quadrature keyword "'// trim(quad_scheme) // '" not recognized. Stopping...'
    stop
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
SUBROUTINE set_lobatto_quadratures(ntheta, nphi, vec, wts, up_and_down)
  ! -------------------------------------------------------------------
  ! Computes Lobatto quadratures for a given number of theta and phi
  ! angles. It then converts the spherical geometry into cartesian geometry
  ! -------------------------------------------------------------------
  integer, intent(in) :: ntheta, nphi
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  real(8), dimension(:), allocatable :: mu, wmu
  logical, optional :: up_and_down
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
     stop 'Error: ntheta is set too high. *ntheta<=5*'
  end select

  call nphi_copies(nphi, mu, wmu, vec, wts, up_and_down=up_and_down)
END SUBROUTINE set_lobatto_quadratures

!===============================================================================
SUBROUTINE set_gauss_quadratures(ntheta, nphi, vec, wts, up_and_down)
  ! -------------------------------------------------------------------
  ! Computes Gauss-Legendre quadratures for a given number of theta and phi
  ! angles. It then converts the spherical geometry into cartesian geometry
  ! -------------------------------------------------------------------
  integer, intent(in) :: ntheta, nphi
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  real(8), dimension(:), allocatable :: mu, wmu
  logical, optional :: up_and_down
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
     stop 'Error: ntheta is set too high. *ntheta<=6*'
  end select

  call nphi_copies(nphi, mu, wmu, vec, wts, up_and_down=up_and_down)
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
  call nphi_copies(nphi, mu, wmu, vec, wts, zrot, up_and_down=.false.)
END SUBROUTINE set_custom_quadratures

!===============================================================================
SUBROUTINE set_radau_quadratures(ntheta, nphi, vec, wts, up_and_down)
  integer, intent(in) :: ntheta, nphi
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  logical, intent(in), optional :: up_and_down
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
  call nphi_copies(nphi, mu, wmu, vec, wts, zrot, up_and_down=up_and_down)
END SUBROUTINE set_radau_quadratures

!===============================================================================
SUBROUTINE nphi_copies (nphi_in, mu, wmu, vec, wts, zrot, up_and_down)
  integer, intent(in) :: nphi_in
  real(8), intent(inout), dimension(:) :: mu, wmu
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  real(8), intent(in), optional, dimension(:) :: zrot
  real(8) :: xphi, theta, xoff
  integer :: nrays, ntheta, i, nt, np, nmuphi, n1, nphi
  logical :: includes_mu1, mirror=.false.
  logical, optional :: up_and_down
  ! -------------------------------------------------------------------
  if (present(up_and_down)) then
    mirror = up_and_down
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
    !wts = wts * 2*PI
    wts = wts / 2
  else
    wts = wts * vec(3,:)
  endif
END SUBROUTINE nphi_copies

!===============================================================================
SUBROUTINE ray_init(ray, vec)
  class(ray_type), intent(inout) :: ray
  real(8), intent(in) :: vec(3)
  ! ----------------------------------------------------------------------------
  ray%mux = vec(1)
  ray%muy = vec(2)
  ray%muz = vec(3)

  if (ray%mux.gt.0) then
    ! x sweep in direction of increasing x index
    ray%iup = -1
    ray%idn =  1
  else
    ! x sweep in direction of decreasing x index
    ray%iup =  1
    ray%idn = -1
  end if

  if (ray%muy.gt.0) then
    ! y sweep in direction of increasing y index
    ray%jup = -1
    ray%jdn =  1
  else
    ! y sweep in direction of decreasing y index
    ray%jup = 1
    ray%jdn = -1
  end if

  if (ray%muz.gt.0) then
    ! z sweep in direction of increasing z index
    ray%kup = -1
    ray%kdn =  1
  else
    ! z sweep in direction of decreasing z index
    ray%kup =  1
    ray%kdn = -1
  end if
END SUBROUTINE ray_init

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
  call SetCcKern(ray%upkern,h1,h2)
END SUBROUTINE set_stencil
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE GetIntersect(rx,ry,rz,plane)
  !---------------------------------------------------------------------
  !
  ! Determines which plane is hit first by a ray
  !
  ! ds[x,y,z] : distance to first intersection with [yz,xz,xy]= plane
  !
  ! thus, if dsx < dsy and dsx < dsz then first intersection is yz-plane
  ! 
  ! ds[x,y,z]= 1/abs(r[x,y,z]), so this inequality becomes 
  !
  ! rx > ry and rx > rz
  !
  !---------------------------------------------------------------------
  real(8),intent(in) :: rx,ry,rz
  integer,intent(out) :: plane
  real(8) :: rxt,ryt,rzt

  rxt=abs(rx)
  ryt=abs(ry)
  rzt=abs(rz)

  if (rzt.ge.rxt.and. rzt.ge.ryt) then
     plane = 3 !xy plane
  else if (ryt.ge.rxt.and.ryt.ge.rzt) then
     plane = 2 !xz plane
  else
     plane = 1 !yz plane
  endif

END SUBROUTINE getintersect
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE CalcFrac(wt1,wt2,rx,ry,rz,plane)
  !---------------------------------------------------------------------
  ! Calculate coordinate of point where the characteristic hits cell 
  ! boundary in units of cell size, relative to the central point of the 
  ! chracteristic.
  !
  ! Ds is length of ray to first plane intersection. 
  !
  ! Example if first plane is yz-plane:
  ! dy is y-coordinate of intersection point if the origin is placed
  ! in the center point of the short characteristic.
  ! Dy is the y-size of the cell
  !
  !  Ds = dx/abs(mux) = 1 / abs(rx)
  ! 
  !  dy = muy * Ds = muy / abs(rx)
  ! 
  ! in units of cell size (wt1) this is:
  !
  ! wt1 = dy/Dy = muy/Dy * 1/abs(rx) = ry/abs(rx)
  !
  ! Note: -1 < wt[12] < 1
  !---------------------------------------------------------------------
  real(8), intent(out) :: wt1,wt2
  real(8), intent(in)  :: rx,ry,rz
  integer, intent(in)  :: plane

  select case(plane) 
  case(1) ! yz plane
     wt1 = ry/abs(rx)
     wt2 = rz/abs(rx)
  case(2) ! xz plane
     wt1 = rx/abs(ry)
     wt2 = rz/abs(ry)
  case(3) ! xy plane
     wt1 = rx/abs(rz)
     wt2 = ry/abs(rz)
  end select

END SUBROUTINE CalcFrac
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE CalcDistance(ds,plane,rx,ry,rz)
  !---------------------------------------------------------------------
  !
  ! Computes distance ds between center point of the short charactersitic
  ! and the first plane intersection
  !
  ! ds = abs(mu[x,y,z]) / d[x,y,z] = 1/abs(r[x,y,z]) 
  !
  !---------------------------------------------------------------------
  real(8), intent(in)  :: rx,ry,rz
  integer, intent(in)  :: plane
  real(8), intent(out) :: ds

  select case(plane)
  case(1) 
     !YZ-plane
     ds = 1.0d0/abs(rx)
  case(2) 
     ! XZ-plane
     ds= 1.0d0/abs(ry)
  case(3) 
     ! XY-plane
     ds = 1.0d0/abs(rz)
  end select

END SUBROUTINE CalcDistance
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE SetShortKern(ray,wt1,wt2,isup)
  !---------------------------------------------------------------------
  !
  ! Set up the interpolation kernel for short characteristics.
  ! We assume a 4x4 kernel with the following geometry 
  ! (example uses xy-plane intersection of the upwind point with 
  ! mux<0 and muy<0)
  !
  ! 0 is the (x,y) coordinate of the center point of the short
  !   charactersitic
  ! + is the point where the ray crosses the xy-plane
  ! (upo1,upo2) are the indices of point of the kernel closest
  !   to the origin
  !
  !
  !           y
  !           ^
  !           | 
  !
  !       1   *-----*-----*-----*
  !           |     |     |     |
  !           |     |     |     |
  !       0   *-----*-----0-----*
  !           |     |     |     |
  !           |     |  +  |     |
  !      -1   *-----*-----*-----*
  !           |     |     |     |
  !           |     |     |     |
  ! upo2 -2   *-----*-----*-----*   --->x
  ! 
  !         -2    -1    0     1
  !        upo1
  !   
  !---------------------------------------------------------------------
  real(8),intent(in) :: wt1,wt2
  logical,intent(in) :: isup
  type(ray_type),intent(inout) :: ray
  real(8) :: h1,h2

  ! set weights such it is relative to the surrounding point closest
  ! to the origin:  0  <  h[12] < 1
  if (wt1.lt.0) then
     h1=1.0d0+wt1
     if (isup) then
        ray%upo1=-2
     else
        ray%dno1=-2
     endif
  else
     h1=wt1
     if (isup) then
        ray%upo1=-1
     else
        ray%dno1=-1
     endif
  end if

  if (wt2.lt.0) then 
     h2=1.0d0+wt2
     if (isup) then
        ray%upo2=-2
     else
        ray%dno2=-2
     endif
  else
     h2=wt2
     if (isup) then
        ray%upo2=-1
     else
        ray%dno2=-1
     endif
  end if

  if (isup) then 
     select case(ray%upplane) 
     case(1:2) !xz and yz
        ! z grid is nonequidistant, so we set a linear interpolation
        ! kernel on a 2x2 grid
        call SetLinKern(ray%upkern,h1,h2)  
     case(3) !xy   
        ! xy grids are equidistant, so we can do cubic convolution on  a
        ! 4x4 grid 
        call SetCcKern(ray%upkern,h1,h2)
     end select
  else
     select case(ray%dnplane) 
     case(1:2) !xz and yz
        ! z grid is nonequidistant, so we set a linear interpolation
        ! kernel on a 2x2 grid
        call SetLinKern(ray%dnkern,h1,h2)  
     case(3) !xy   
        ! xy grids are equidistant, so we can do cubic convolution on  a
        ! 4x4 grid 
        call SetCcKern(ray%dnkern,h1,h2)
     end select
  endif

END SUBROUTINE SetShortKern
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE SetLongKern(ray,wt1,wt2,longfracX,longfracY)
  type(ray_type),intent(inout) :: ray
  real(8)       ,intent(in)    :: wt1,wt2,longfracX,longfracY
  real(8) :: h1,h2

  ! set weights such that it is relative to the surrounding point
  ! closest to the origin: 0 < h[12] < 1
  if (wt1.lt.0) then
     h1=1.0d0+wt1
     ray%longo1=-2+int(longfracX)
  else
     h1=wt1
     ray%longo1=-1+int(longfracX)
  end if
  if (wt2.lt.0) then 
     h2=1.0d0+wt2
     ray%longo2=-2+int(longfracY)
  else
     h2=wt2
     ray%longo2=-1+int(longfracY)
  end if
  ! xy grid is equidistant, so we can do cubic convolution    
  call SetCcKern(ray%longkern,h1,h2)

END SUBROUTINE SetLongKern
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE SetCcKern(kern,h1,h2)
  USE mmath_mod,only:cc_kernel
  !---------------------------------------------------------------------
  !
  ! set a 2D cubic convulution kernel as the multiplication
  ! of two 1D cubic convolutions
  !
  !---------------------------------------------------------------------
  real(8),intent(in) :: h1,h2
  real(8),dimension(4,4),intent(out) :: kern
  real(8),dimension(4) :: xkern,ykern

  integer :: i,j

  call cc_kernel(h1,xkern)
  call cc_kernel(h2,ykern)
  forall (i=1:4,j=1:4)
     kern(i,j) = xkern(i)*ykern(j)
  endforall

END SUBROUTINE SetCcKern
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE SetLinKern(kern,h1,h2)
  !---------------------------------------------------------------------
  ! set the inner 2x2 grid of a 4x4 kernel to a linear interpolation
  ! kernel
  !
  !          dir2
  !           ^
  !           | 
  !
  !       4   *-----*-----*-----*
  !           |     |     |     |
  !           |     |     |     |
  !       3   *-----*-----0-----*
  !           |     |     |     |
  !           |     |  +  |     |
  !       2   *-----*-----*-----*
  !           |     |     |     |
  !           |     |     |     |
  !       1   *-----*-----*-----*   ---> dir1
  ! 
  !           1     2     3     4
  !        
  !---------------------------------------------------------------------
  real(8),intent(in) :: h1,h2
  real(8),dimension(4,4),intent(out) :: kern

  kern=0.0d0
  kern(2,2)=(1.0-h1)*(1.0-h2)
  kern(3,2)=     h1 *(1.0-h2)
  kern(2,3)=(1.0-h1)*     h2
  kern(3,3)=     h1 *     h2  

END SUBROUTINE SetLinKern
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
    stop 'less than 3 weights in simps_weights'
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

END MODULE ray_mod