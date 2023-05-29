MODULE ray_mod
  USE mconst
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
      wts = 1.0
    else
      allocate(vec(3, 1), wts(1))
      vec(:,1) = [0.0, 0.0, 1.0]
      wts = 1.0
    endif
    wts = wts * 2 * pi
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

  ! Normalize, so that integration over 8 octants yields 4 pi
  ! the above sets are normalized so that integration over one octant yields 1.
  wmu = wmu * 0.5 * pi

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

  if (mirror) then
    ! Fill the other octants with inverted rays 
    wts(  4*nro+1:8*nro) =  wts(  1:4*nro)
    vec(:,4*nro+1:8*nro) = -vec(:,1:4*nro)
  endif
END SUBROUTINE set_carlson_quadratures
!---------------------------------------------------------------------

!===============================================================================
SUBROUTINE set_grid_quadratures (quad_scheme, ds, vec, wts, up_and_down)
  character(len=*), intent(in) :: quad_scheme
  real(8), intent(in) :: ds(3)
  integer :: m, i, nro, nmu
  real(8), dimension(:), allocatable :: mu, zrot
  real(8), dimension(:), allocatable, intent(out) :: wts
  real(8), dimension(:,:), allocatable, intent(out) :: vec
  real(8) :: plus, cross
  logical, optional :: up_and_down
  logical :: mirror=.true.
  ! ----------------------------------------------------------------------------
  if (present(up_and_down)) then
    mirror = up_and_down
  endif

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
SUBROUTINE set_lobatto_quadratures(ntheta, nphi, vec, wts)
  ! -------------------------------------------------------------------
  ! Computes Lobatto quadratures for a given number of theta and phi
  ! angles. It then converts the spherical geometry into cartesian geometry
  ! -------------------------------------------------------------------
  integer, intent(in) :: ntheta, nphi
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  real(8), dimension(:), allocatable :: mu, wmu
  integer :: nrays
  ! -------------------------------------------------------------------
  nrays  = ntheta + 1

  allocate(mu(nrays))
  allocate(wmu(nrays))

  ! setting indexes and hard-coded Lobatto coefficients
  select case (nrays)
  case (1)
     vec(:,1) = (/0.0, 0.0, 1.0/)
     wts = 1.0
     return
     
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
     stop 'Error: ntheta is set too high. *ntheta<5*'

  end select

  wmu = wmu * mu * 2 * PI
  call nphi_copies(nphi, mu, wmu, vec, wts)
END SUBROUTINE set_lobatto_quadratures

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
  wmu = wmu * mu * 2 * PI
  if (present(zrot)) zrot = zrot(order(nrays:1:-1))

  call nphi_copies(nphi, mu, wmu, vec, wts, zrot)
END SUBROUTINE set_custom_quadratures

!===============================================================================
SUBROUTINE set_radau_quadratures(quad_scheme, vec, wts, up_and_down)
  character(len=*), intent(in) :: quad_scheme
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  real(4), dimension(:), allocatable :: wmu2, mu2
  real(8), dimension(:), allocatable :: wmu, mu
  real(8), dimension(:), allocatable :: wts2, zrot
  real(8), dimension(:,:), allocatable :: vec2
  integer :: nrays, ntheta, nphi
  logical, optional :: up_and_down
  logical :: mirror=.true.
  ! -------------------------------------------------------------------
  select case (quad_scheme)
  case ('radau1')
    ntheta = 1
    nphi = 1
  case ('radau2')
    ntheta = 2
    nphi = 4
  case ('radau4')
    ntheta = 3
    nphi = 4
  end select

  allocate(mu(ntheta), wmu(ntheta))
  allocate(mu2(ntheta), wmu2(ntheta))
  allocate(zrot(ntheta))
  call radau%init(ntheta, mu2, wmu2)
  mu = mu2(ntheta:1:-1)
  wmu = wmu2(ntheta:1:-1)
  zrot = PI/4.
  call nphi_copies(nphi, mu, wmu, vec2, wts2, zrot=zrot)
  nrays = size(wts2)

  if (present(up_and_down)) then
    mirror = up_and_down
  endif

  if (mirror) then
    allocate(wts(2*nrays))
    allocate(vec(3, 2*nrays))
    wts(  1:nrays) = wts2
    vec(:,1:nrays) = vec2
    wts(  nrays+1:2*nrays) =  wts2
    vec(:,nrays+1:2*nrays) = -vec2
    wts = wts * 2*PI
  else
    allocate(wts(nrays))
    allocate(vec(3, nrays))
    wts = wts2
    vec = vec2
  endif
END SUBROUTINE set_radau_quadratures

!===============================================================================
SUBROUTINE nphi_copies (nphi, mu, wmu, vec, wts, zrot)
  integer, intent(in) :: nphi
  real(8), intent(inout), dimension(:) :: mu, wmu
  real(8), intent(out), dimension(:), allocatable :: wts
  real(8), intent(out), dimension(:,:), allocatable :: vec
  real(8), intent(in), optional, dimension(:) :: zrot
  real(8) :: xphi, theta, xoff
  integer :: nrays, ntheta, i, nt, np, nmuphi, n1
  logical :: includes_mu1
  ! -------------------------------------------------------------------
  nrays = size(mu)
  includes_mu1 = findloc(mu, 1.0, 1) > 0
  if (includes_mu1) then
    nmuphi = (nrays-1)*nphi + 1
  else
    nmuphi = nrays*nphi
  endif

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
END SUBROUTINE nphi_copies

!===============================================================================
SUBROUTINE ray_init(ray, vec, wmu)
  class(ray_type), intent(inout) :: ray
  real(8), intent(in) :: vec(3), wmu
  ! ----------------------------------------------------------------------------
  ray%mux = vec(1)
  ray%muy = vec(2)
  ray%muz = vec(3)
  ray%wmu = wmu

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
SUBROUTINE set_stencil(ray, ds)
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
  USE mmath, only:qsortd

  class(ray_type) :: ray
  real(8), intent(in) :: ds(3)
  integer :: i,j,ncrossXZ,ncrossYZ,slct
  integer :: count,Nvalid,nmc
  integer, parameter :: mb=3
  real(8) :: dx,dy,dz, mux,muy,muz
  real(8) :: wt1, wt2, rx, ry, rz
  real(8) :: longfracx,longfracy 
  real(8),dimension(:),allocatable :: ts,tx,ty,tz
  integer,dimension(:),allocatable :: tp,order
  logical,dimension(:),allocatable :: valid
  character(len=130) :: txt

  ! maximum number of cell boundaries a long ray may cross
  nmc=mb-2

  ! temporary storage for 
  allocate(tx(2*(mb-2))) ! x distance to intersection
  allocate(ty(2*(mb-2))) ! y distance to intersection
  allocate(ts(2*(mb-2))) ! abs distance to plane intersection
  allocate(tz(2*(mb-2))) ! z distance to intersection
  allocate(tp(2*(mb-2))) ! which plane is intersected
  allocate(order(2*(mb-2))) ! order in which planes are intersected
  allocate(valid(2*(mb-2))) ! .true. if intersection is used: in case
  ! two intersections occur very close together
  ! one of the two is discarded)

  ! first we determine the downwind quantities

  ! cell sizes in the downwind direction
  dx=ds(1)
  dy=ds(2)
  dz=ds(3)

  ! directional cosines divided by cell size
  mux = ray%mux
  muy = ray%muy
  muz = ray%muz
  rx = mux / dx
  ry = muy / dy
  rz = muz / dz          
  
  call GetIntersect(rx,ry,rz,ray%dnplane)
  call CalcFrac(wt1,wt2,rx,ry,rz,ray%dnplane)
  call CalcDistance(ray%dsdn,ray%dnplane,rx,ry,rz)       
  call SetShortKern(ray,wt1,wt2,.false.)

  ! now we determine the upwind quantities
  ! the upwind point has the angles reversed
  mux=-mux
  muy=-muy
  muz=-muz
  rx=mux / dx
  ry=muy / dy
  rz=muz / dz          

  call GetIntersect(rx,ry,rz,ray%upplane)
  call CalcFrac(wt1,wt2,rx,ry,rz,ray%upplane)
  call CalcDistance(ray%dsup,ray%upplane,rx,ry,rz) 
  call SetShortKern(ray,wt1,wt2,.true.)
  ray%long=.false.

  if (ray%upplane.eq.1.or.ray%upplane.eq.2) then
    ! The upwind point of the short characteristic is lies on the           
    ! xz or the yz-plane. Because of the periodic boundary conditions we do
    ! not know the incoming intensity on the kernel points.
    ! We do know however the intensity on all the points in the xy-plane
    ! in the upwind direction.
    ! In order to get the intensity on the upwind point of the short char,
    ! we have to trace the shortchar upwind to the first intersection with 
    ! the xy-plane

    ! found one more long ray, and flag this ray is long
    ray%long=.true.
    if(ray%long) write(txt,*) 'long rays'

    ! distance to next XY-plane intersection      
    ray%longds=abs(1./rz)

    ! offset of ray to next XY plane crossing in units of cell size, relative
    ! to center point of the short char, i.e., if mux<0 then longfracX < 0
    longfracX = rx / abs(rz)
    longfracY = ry / abs(rz)

    ! number of XZ or YZ plane crossings 
    ! example: if longfracX = -3.2 then ncrossYZ = 3
    ncrossYZ = abs(int(longfracX))
    ncrossXZ = abs(int(longfracY))

    ! set the cubic convolution kernel for the XY crossing, but
    ! make sure the ray does not cross more than mb-2
    ! vertical planes in each direction otherwise the
    ! cubic-convolution will not work

    ! fractional part of longfrac[X,Y]
    wt1=longfracX-int(longfracX)
    wt2=longfracY-int(longfracY)

    ! if long ray extends to far in the x or y direction we cheat and truncate the ray.
    if (ncrossYZ.gt.nmc.or.ncrossXZ.gt.nmc) then
        ! one process per plane of processes issues a warning
        print *, 'warning:long ray extends upwind beyond boundary zones, truncating at'
        write(txt,fmt='(A,3(F8.4),A,2(F8.4))') 'mux,muy,muz=',mux,muy,muz,',  lf1,lf2=',longfracX,longfracY
        print *, txt
    endif

    ! check whether long ray extends more than mb-1 cell crossing in the x or y-direction.
    slct=0
    ! ray extends too far in the x-direction
    if (ncrossYZ.gt.nmc.and.ncrossXZ.le.nmc) slct=1
    ! ray extends too far in the y-direction                 
    if (ncrossYZ.le.nmc.and.ncrossXZ.gt.nmc) slct=2
    ! ray extends too far in both the x and y-direction 
    if (ncrossYZ.gt.nmc.and.ncrossXZ.gt.nmc) then
        if ( abs(rx).gt.abs(ry) ) then
          ! ray hits max number of YZ crossings first
          slct=1
        else
          ! ray hits max number of XZ crossings first
          slct=2
        endif
    end if

    ! if ray extends to far, then truncate
    select case(slct)
    case(1)
        ncrossYZ = nmc
        if (mux.lt.0) then
          longfracX =-nmc-.99
          wt1=-0.99
        else
          longfracX = nmc+.99
          wt1= 0.99                     
        endif
        ! 'lenght' of truncated ray
        ray%longds=abs(longfracX*dx/mux)
        ! new number of YZ-plane crossings
        longfracY= ray%longds*muy/dy                   
        ncrossXZ = abs(int(longfracY))
        wt2=longfracY-int(longfracY)

    case(2)
        ncrossXZ = nmc
        if (muy.lt.0) then
          longfracY =-nmc-.99
          wt1=-0.99
        else
          longfracY = nmc+.99
          wt1= 0.99                     
        endif
        ! 'lenght' of truncated ray
        ray%longds=abs(longfracY*dy/muy)
        ! new number of XZ-plane crossings
        longfracX= ray%longds*mux/dx                   
        ncrossYZ = abs(int(longfracX))
        wt2=longfracX-int(longfracX)
    case default
        ! don't do anything
    end select


    ! set up the interpolation kernel for the intersection with the XY plane
    call SetLongKern(ray,wt1,wt2,longfracX,longfracY)

    ! now put some information about the different
    ! vertical plane intersections in temporary arrays

    ! first zero arrays from previous iteration
    tx=0.0
    ty=0.0
    tz=0.0
    ts=0.0
    tp=0
    order=0
    valid=.false.

    count=0
    do i=1,ncrossYZ
        count = count+1
        tp(count)=1             ! YZ plane
        ts(count)=i*dx/abs(mux) ! abs distance to plane intersection from center point of short char
        tx(count)=mux*ts(count)  ! x distance to intersection
        ty(count)=muy*ts(count)  ! y distance to intersection
        tz(count)=muz*ts(count)  ! z distance to intersection
    end do
    do i=1,ncrossXZ
        count = count+1 
        tp(count)=2             ! XZ plane
        ts(count)=i*dy /abs(muy) ! abs distance to plane intersection
        tx(count)=mux*ts(count)  ! x distance to intersection
        ty(count)=muy*ts(count)  ! y distance to intersection
        tz(count)=muz*ts(count)  ! z distance to intersection
    end do

    ! now we sort the intersections according to the
    ! distance to the point of origin.
    call qsortd(ts,order,count)
    tp(1:count)=tp(order(1:count))
    ts(1:count)=ts(order(1:count))
    tx(1:count)=tx(order(1:count))
    ty(1:count)=ty(order(1:count))
    tz(1:count)=tz(order(1:count))

    ! Determine if any of the intersections appear too close to one
    ! another, for instance when they occur on an intersection between YZ
    ! and XZ planes 
    Nvalid=count
    valid(1)=.true. ! the first intersection (as measured from center point of short char)
    ! is always valid
    do i=2,count
        if ( (ts(i)-ts(i-1)) .lt. 0.05*min(dx,dy) ) then
          valid(i)=.false.
          Nvalid=Nvalid-1
        else
          valid(i) = .true.
        end if
    end do

    ! store number of valid intersections
    ray%Nintersect=Nvalid

    ! fill the vertical wall intersections in the ray struct; storage in
    ! order of decreasing distance ds to point of origin
    j=0

    do i=count,1,-1

        ! don't use this intersection
        if (.not.valid(i)) cycle

        ! valid intersection found, add 1 to counter 
        j=j+1

        ! distance to central point
        ray%int(j)%ds=ts(i)

        if (tp(i).eq.1) then ! YZ plane intersection

          ray%int(j)%plane=1

          ! set linear interpolation weights
          ! wt1 and wt2 are x and y distances of intersection
          ! to the surrounding point that lies closest to the origin
          ! in units of cell size
          wt1 = (ty(i)/dy)-aint(ty(i)/dy)
          if (muy.lt.0) wt1=1.0+wt1
          wt2 = (tz(i)/dz)-aint(tz(i)/dz)
          if (muz.lt.0) wt2=1.0+wt2

          ray%int(j)%wght(1,1)=(1.0-wt1)*(1.0-wt2)  !lower left
          ray%int(j)%wght(2,1)=     wt1 *(1.0-wt2)  !lower right
          ray%int(j)%wght(1,2)=(1.0-wt1)*     wt2   !upper left
          ray%int(j)%wght(2,2)=     wt1 *     wt2   !upper right

          ! set the offset of the x coordinate of the 4
          ! points surrounding the intersection relative to
          ! the central shortchar point
          ray%int(j)%iof(:,:)=nint(tx(i)/dx) 

          ! the same for the y coordinates
          ray%int(j)%jof(1,1)=  floor(ty(i)/dy) !lower left
          ray%int(j)%jof(2,1)=ceiling(ty(i)/dy) !lower right
          ray%int(j)%jof(1,2)=  floor(ty(i)/dy) !upper left
          ray%int(j)%jof(2,2)=ceiling(ty(i)/dy) !upper right

        else !XZ plane intersection

          wt1=tx(i)/dx-aint(tx(i)/dx)
          if (mux.lt.0) wt1=1.0+wt1
          wt2=tz(i)/dz-aint(tz(i)/dz)
          if (muz.lt.0) wt2=1.0+wt2

          ray%int(j)%plane=2

          ray%int(j)%wght(1,1)=(1.0-wt1)*(1.0-wt2)  !lower left
          ray%int(j)%wght(2,1)=     wt1 *(1.0-wt2)  !lower right
          ray%int(j)%wght(1,2)=(1.0-wt1)*     wt2   !upper left
          ray%int(j)%wght(2,2)=     wt1 *     wt2   !upper right

          ray%int(j)%jof(:,:)=nint(ty(i)/dy) 

          ray%int(j)%iof(1,1)=  floor(tx(i)/dx) !lower left
          ray%int(j)%iof(2,1)=ceiling(tx(i)/dx) !lower right
          ray%int(j)%iof(1,2)=  floor(tx(i)/dx) !upper left
          ray%int(j)%iof(2,2)=ceiling(tx(i)/dx) !upper right

        end if

        ! set the offset of the z coordinates of the 4
        ! points surrouding the intersection relative to
        ! central ray point
        if (ray%muz>0) then
          ! ray points in direction of increasing z index
          ray%int(j)%kof(1,1)=-1  !lower left
          ray%int(j)%kof(2,1)=-1  !lower right
          ray%int(j)%kof(1,2)= 0  !upper left
          ray%int(j)%kof(2,2)= 0  !upper right
        else
          ray%int(j)%kof(1,1)= 0  !lower left
          ray%int(j)%kof(2,1)= 0  !lower right
          ray%int(j)%kof(1,2)= 1  !upper left
          ray%int(j)%kof(2,2)= 1  !upper right
        end if

    end do

    ! now set distances between intersections
    do i=ray%Nintersect,2,-1
        ray%int(i)%ds= ray%int(i-1)%ds - ray%int(i)%ds
    end do
    ray%int(1)%ds=ray%longds-ray%int(1)%ds

  end if 
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
  USE mmath,only:cc_kernel
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

  if (n < 3) then
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