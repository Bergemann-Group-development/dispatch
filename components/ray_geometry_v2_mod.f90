!> =====================================================================
!> A module to set up ray geometry. Several different geometries can be
!> created if needed for different experiment geometries (e.g. cartesian,
!> spherical, etc.). As an input the module needs the geometry index,
!> lower and upper indices of the mesh (li, ui) and optionally cell
!> sizes (otherwse they are assumed to be cubic and with unit size walls).
!> The resulting heararchy is such:
!> ray ⊂ bundle ⊂ direction ⊂ geometry (⊂ indicates that A is a
!> subset of B).
!> to access levels: ray_geometry%ray_list%bundle%rays

!> ray_geometry_t:  has pointers to optional different geometries, and contains
!> all the information about itself - li, ui, cell dimensions, n_ray_dirs
!> and a pointer to a list of ray directions

!> ray_dirs_t: is a linked list of directions. Every direction has an unique ID
!> for particular geometry, number of bundles and a directionality parameter,
!> which indicates, whether rays should go both ways or to only some particular
!> direction and a pointer to bundles.

!> ray_bndl_t: is a linked list of ray bundles. Every bundle has an unique ID
!> for particular direction, direction ID for safekeeping, number of ray packets
!> it contains, which face the bundle originates/terminates at (lwall and uwall
!> accordingly), axis, along which the bundle is oriented (e.g. x|y|z in cartesian),
!> theta and phi values of slanted rays.

!> ray_t: is a ray packet, containing n_rays of rays of same length (n_points).
!> it has all the necessary info about the ray - where it originates, where it
!> terminates (fwd_coords/rev_coords, depending on direction), which cells
!> it crosses (cell_coords), distance it traves through particular cell(ds),
!> and which cell face it crosses when it enters and leaves the cell (pt_wall)
!> =====================================================================
MODULE ray_geometry_v2_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  implicit none
  private
! ----------------------------------------------------------------------
  type, public :: ray_t
    class(ray_t), pointer :: head => null()
    class(ray_t), pointer :: next => null()
    class(ray_t), pointer :: prev => null()
    class(ray_t), pointer :: tail => null()
    integer :: ID
    integer :: n_rays
    integer :: n_points
    integer, dimension(:,:,:), pointer   :: cell_coords => null()       ! x|y|z,        ray point, ray #
    integer, dimension(:,:  ), pointer   :: fwd_coords  => null()       ! x|y|z,                   ray #
    integer, dimension(:,:  ), pointer   :: rev_coords  => null()       ! x|y|z,                   ray #
    real,    dimension(:,:  ), pointer   :: ds          => null()       ! ds   ,        ray point, ray #
    integer, dimension(:,:,:,:), pointer :: pt_wall     => null()       ! x|y|x, lw|uw, ray point, ray #
  contains
    procedure :: fill_ray_list_data
    procedure :: print => print_ray
  end type
! ----------------------------------------------------------------------
  type, public :: ray_bndl_t
    class(ray_bndl_t),  pointer :: head => null()
    class(ray_bndl_t),  pointer :: next => null()
    class(ray_bndl_t),  pointer :: prev => null()
    class(ray_bndl_t),  pointer :: tail => null()
    class(ray_t), pointer :: rays => null()
    integer :: ID
    integer :: dir_ID
    integer :: n_packets                                                ! number of packets in this bundle
    integer :: lwall(3)
    integer :: uwall(3)
    integer :: axis                                                     ! identificator along which axis ray goes 1-3
    real    :: theta                                                    ! angle between ray and axis [0, 45], stored in radians
    real    :: phi                                                      ! angle of the ray in the axis plane [0, 360), stored in radians
  contains
    procedure :: print => print_bundle
  end type
! ----------------------------------------------------------------------
  type, public :: ray_dirs_t
    class(ray_dirs_t), pointer :: head => null()
    class(ray_dirs_t), pointer :: next => null()
    class(ray_dirs_t), pointer :: prev => null()
    class(ray_dirs_t), pointer :: tail => null()
    class(ray_bndl_t), pointer :: bundle => null()
    integer :: ID
    integer :: n_bundles                                                ! number of bundles in this direction
    integer :: directionality(3)                                        ! 1 - two directions, 2 - only forward, 3 - only reverse (x|y|z)
  contains
    procedure :: print => print_dirs
  end type
! ----------------------------------------------------------------------
  type, public :: ray_geometry_t
    class(ray_geometry_t), pointer :: head => null()
    class(ray_geometry_t), pointer :: next => null()
    class(ray_geometry_t), pointer :: prev => null()
    class(ray_geometry_t), pointer :: tail => null()
    class(ray_dirs_t), pointer :: ray_list => null()
    integer :: n_ray_dirs
    integer :: li(3), ui(3)
    integer :: cell_dims(3)
    integer :: geometry = 1 !cartesian = 1, spherical = 2, cylindrical = 3
    logical :: init_done = .false.
    integer :: axis, dr(3)                                              ! convenience
    real    :: theta, phi                                               ! convenience
  contains
    procedure :: add_axis
    procedure :: add_diag
    procedure :: add_ray
    procedure :: finish_single_ray_list
    procedure :: get_my_bundle
    procedure :: get_my_geometry
    procedure :: init
    procedure :: init_cartesian
    procedure :: init_spherical
    procedure :: init_cylindrical
    procedure :: test_ray_geometry
    procedure :: traverse_2d
    procedure :: traverse_3d
    procedure :: print => print_geometry
  end type
! ----------------------------------------------------------------------
  type :: single_ray_t
    class(single_ray_t), pointer :: head => null()
    class(single_ray_t), pointer :: next => null()
    class(single_ray_t), pointer :: prev => null()
    class(single_ray_t), pointer :: tail => null()
    integer                          :: lwall(3)
    integer                          :: uwall(3)
    integer                          :: n_points
    integer, dimension(:,:), pointer :: cell_coords => null()  ! x|y|z, ray point
    integer, dimension(3)            :: fwd_coords             ! x|y|z,
    integer, dimension(3)            :: rev_coords             ! x|y|z,
    real,    dimension(:  ), pointer :: ds          => null()  ! ds   , ray point
    integer, dimension(:,:,:), pointer :: cell_walls => null()  ! lw|uw(x,y,z), ray point
  contains
    procedure :: fill_single_ray_data
  end type
! ----------------------------------------------------------------------
  type :: single_ray_packet_t
    class(single_ray_packet_t), pointer :: head => null()
    class(single_ray_packet_t), pointer :: next => null()
    class(single_ray_packet_t), pointer :: prev => null()
    class(single_ray_packet_t), pointer :: tail => null()
    class(single_ray_t), pointer:: rays => null()
    integer :: n_rays, n_points
  end type
! ----------------------------------------------------------------------
  type :: single_ray_bundle_t
    class(single_ray_bundle_t), pointer :: head => null()
    class(single_ray_bundle_t), pointer :: next => null()
    class(single_ray_bundle_t), pointer :: prev => null()
    class(single_ray_bundle_t), pointer :: tail => null()
    class(single_ray_packet_t), pointer :: packets  => null()
    integer :: lwall(3)
    integer :: uwall(3)
  contains
    procedure :: finish_ray
  end type
! ----------------------------------------------------------------------
  type :: ray_point_t
    class(ray_point_t), pointer :: head => null()
    class(ray_point_t), pointer :: next => null()
    class(ray_point_t), pointer :: prev => null()
    class(ray_point_t), pointer :: tail => null()
    integer                     :: n_points
    integer                     :: coords(3)
    integer                     :: fwd(3)
    integer                     :: rev(3)
    integer                     :: uwall(3)
    integer                     :: lwall(3)
    real                        :: ds
  contains
    procedure :: add_ray_point
  end type
! ----------------------------------------------------------------------
  type :: traverse_2d_t
    integer, dimension(3) :: cell_max, cell_min, cell_step, cell_pos
    real, dimension(2) :: dstep, ray_step, ray_dist, ray_p_init, ray_p_last
    real               :: r_tot, r_trace, ds
    integer            :: a, b, fc, uw(3), lw(3)
    logical            :: oob
  contains
    procedure :: init           => init_traverse_2D
    procedure :: posit          => posit_2D
    procedure :: hop            => hop_2D
    procedure :: restart        => restart_2D
    procedure :: get_downstream => get_downstream_2D
  end type
! ----------------------------------------------------------------------
  type :: traverse_3d_t
    integer, dimension(3) :: cell_max, cell_min, cell_step, cell_pos
    real, dimension(3) :: dstep, ray_step, ray_dist, ray_p_init, ray_p_last
    real               :: r_tot, r_trace, ds
    integer            :: a, b, c, uw(3), lw(3)
    logical            :: oob, reset_a, reset_b
  contains
    procedure :: init           => init_traverse_3D
    procedure :: posit          => posit_3D
    procedure :: hop            => hop_3D
    procedure :: restart        => restart_3D
    procedure :: get_downstream => get_downstream_3D
    procedure :: wall           => cw_3d
  end type
! ----------------------------------------------------------------------
  type(ray_geometry_t), public :: ray_geometry
  real(8), public, parameter:: pi = asin(1.0d0) * 2.0d0, radian = 180./pi
CONTAINS

!> =====================================================================
!> Initialisation subroutine for ray geometry
!> =====================================================================
SUBROUTINE init (self, geometry, li, ui, cell_dims)
  class(ray_geometry_t) :: self
  class(ray_geometry_t), pointer :: new
  integer, optional :: cell_dims
  integer :: li(3), ui(3)
  integer :: geometry
! ......................................................................
  call trace_begin('ray_geometry_t%init')

  allocate (new)
  new%geometry = geometry
  new%li = li
  new%ui = ui

  if(present(cell_dims)) then
      new%cell_dims = cell_dims
    else
      new%cell_dims(:) = 1                                              ! default - cubic cells
  end if

  if(.not.associated(self%head)) then
    self%head => new
    self%tail => new
  else
    new%prev      => self%tail
    new%head      => self%head
    new%prev%next => new
    self%tail     => new
  end if

  select case(geometry)
    case(1)
      call new%init_cartesian
    case(2)
      call new%init_spherical
    case(3)
      call new%init_cylindrical
  end select

  call new%print
  call trace_end
END SUBROUTINE init

SUBROUTINE get_my_geometry(self, geometry, li, ui, cell_dims, my_geom)
  class(ray_geometry_t) :: self
  class(ray_geometry_t), pointer :: this, my_geom
  integer, optional :: cell_dims
  integer :: li(3), ui(3)
  integer :: geometry
  logical :: found
! ......................................................................
  call trace_begin('ray_geometry_t%get_my_geometry')

  if(.not.associated(self%head)) then
    call self%init(geometry, li, ui, cell_dims)
  else
    this => self%head
    found = .false.
    do while (associated(this).and..not.found)
      if ((this%geometry == geometry).and.all(this%li==li).and.&
        all(this%li==li)) then
        my_geom => this
        found = .true.
      else
        this => this%next
      end if
    end do
  end if

  nullify(this)
  call trace_end
END SUBROUTINE get_my_geometry

!> =====================================================================
!> Initialisation subroutine for cartesian geometry.
!> It is important to note, that although rays are created independently
!> for arbitrary directions (or can be not created, if user doesn't need
!> them), 45 degree rays inevitably cover all tree directions, so they
!> should be treated carefully if one doesn't need some directions.
!> Namelist parameters:
!>   nx, ny, nz     : number of angles relative to x-, y-, z-axis
!>   xti, yti, zti  : degrees relative to x-, y-, and z-axis
!>   dx, dy, dz     : currently not used
!> =====================================================================
SUBROUTINE init_cartesian(self)
  class(ray_geometry_t) :: self
  integer                         :: i, d45
  real, dimension(:), allocatable :: xt, yt, zt                         ! FIXME: These seem to not be used
  real, dimension(5)              :: xti=0.0, yti=0.0, zti=0.0
  integer, save                   :: nx=1, ny=0, nz=0
  integer, save                   :: dx=0, dy=0, dz=0                   ! FIXME: These seem to not be used
! ......................................................................
  namelist /cartesian_ray_geometry_params/ nx, ny, nz, xti, yti, zti, dx, dy, dz
! ......................................................................
  call trace_begin('ray_geometry_t%init_cartesian',2)
  !> TO DO: add an option for 2D simulation setup (i.e. any(cell_dims) == 0)

  self%init_done = .true.
  allocate(self%ray_list)
  self%n_ray_dirs = 0

  rewind (io%input)
  read (io%input, cartesian_ray_geometry_params)
  self%dr(1) = dx
  self%dr(2) = dy
  self%dr(3) = dz
  d45 = 0                                                               ! flag for 45 degree angles
  if (nx > 0) then
    allocate (xt(nx)); xt = xti(1:nx)
    do i = 1, nx
      if (xti(i) == 0) then
        call self%add_axis(1)
      else if (xti(i) == 45.) then
        if (d45 < 1) then
          call self%add_diag(1,1)
          d45 = d45 + 1
        else if (d45 < 2) then
          call self%add_diag(2,1)
          d45 = d45 + 1
        end if
      else
        call self%add_ray(xti(i),1)
      end if
    end do
    deallocate (xt)
  end if

  if (ny > 0) then
    allocate (yt(ny)); yt = yti(1:ny)
    do i = 1, ny
      if (yti(i) == 0) then
        call self%add_axis(2)
      else if (yti(i) == 45.) then
       if (d45 < 1) then
          call self%add_diag(1,2)
          d45 = d45 + 1
        else if (d45 < 2) then
          call self%add_diag(2,2)
          d45 = d45 + 1
        end if
      else
        call self%add_ray(yti(i),2)
      end if
    end do
    deallocate (yt)
  end if

  if (nz > 0) then
    allocate (zt(nz)); zt = zti(1:nz)
    do i = 1, nz
      if (zti(i) == 0) then
        call self%add_axis(3)
      else if (zti(i) == 45.) then
        if (d45 < 1) then
          call self%add_diag(1,3)
          d45 = d45 + 1
        else if (d45 < 2) then
          call self%add_diag(2,3)
          d45 = d45 + 1
        end if
      else
        call self%add_ray(zti(i),3)
      end if
    end do
    deallocate (zt)
  end if

  call trace_end
END SUBROUTINE init_cartesian

!> =====================================================================
!> Initialisation subroutine for cylindrical geometry
!> =====================================================================
SUBROUTINE init_cylindrical(self)
  class(ray_geometry_t)      :: self
  integer                         :: i, nz, nR, nPhi
  real, dimension(:), allocatable :: zt, Rt, Phit
  real, dimension(5)              :: zti, Rti, Phiti
! ......................................................................
  namelist /cylindrical_ray_geometry_params/ nz, nR, nPhi, zti, Rti, Phiti
! ......................................................................
  error stop "Cylindrical ray geometry not implemented yet"

END SUBROUTINE init_cylindrical

!> =====================================================================
!> Initialisation subroutine for spherical geometry
!> =====================================================================
SUBROUTINE init_spherical(self)
  class(ray_geometry_t)      :: self
  integer                         :: i, nr, ntheta, nphi
  real, dimension(:), allocatable :: rt, thetat, phit
  real, dimension(5)              :: rti, thetati, phiti
! ......................................................................
  namelist /spherical_ray_geometry_params/ nr, ntheta, nphi, rti, thetati, phiti
! ......................................................................
  error stop "Spherical ray geometry not implemented yet"

END SUBROUTINE init_spherical

!> =====================================================================
!> Subroutine to add axisymmetric rays. It is done separately, because
!> it is much simpler and therefore faster than with rays with other
!> angles. Adds stiff ray arrays instantly.
!> =====================================================================
SUBROUTINE add_axis(self, dir)
  class(ray_geometry_t) :: self
  class(ray_dirs_t), pointer :: ndir
  class(ray_bndl_t), pointer :: bndl
  class(ray_t),      pointer :: ray
  integer :: dir, ia, ib, ic, ray_n, a, b, c, rp, p(3), pmax(3)
  real    :: ds
! ......................................................................
  call trace_begin('ray_geometry_t%add_axis', 4)
  pmax = (self%ui - self%li) + 1                                        ! Total number of points
  self%n_ray_dirs = self%n_ray_dirs + 1                                ! Adding new dir
  SELECT case(dir)
  case(1)
    a = 2; b = 3; c = 1                                                 ! along x axis
  case(2)
    a = 3; b = 1; c = 2                                                 ! along y axis
  case(3)
    a = 1; b = 2; c = 3                                                 ! along z axis
  end SELECT
  ! Create new ray bundle and add it to ray_geometry
  allocate (ndir, bndl, ray)                                            ! create new ray
  if (associated(self%ray_list%head)) then                              ! if ray list is not empty
    ndir%prev => self%ray_list%tail
    ndir%prev%next => ndir
    self%ray_list%tail => ndir
    ndir%head => self%ray_list%head
  else
    self%ray_list => ndir
    self%ray_list%head => ndir                                           ! if empty, point to itself
    self%ray_list%tail => ndir
  end if

  ! --- Fill in direction info ---
  ndir%bundle => bndl
  ndir%n_bundles = 1
  ndir%directionality = self%dr
  ndir%ID = self%n_ray_dirs
  ! --- Fill in bundle info ---
  bndl%head   => bndl
  bndl%tail   => bndl
  bndl%rays   => ray
  bndl%axis  = c
  bndl%theta = 0
  bndl%phi   = 0
  bndl%ID = 1
  bndl%dir_ID = ndir%ID
  bndl%n_packets = 1
  bndl%lwall(a) =  0;  bndl%uwall(a) = 0
  bndl%lwall(b) =  0;  bndl%uwall(b) = 0
  bndl%lwall(c) = -1;  bndl%uwall(c) = 1                                ! lower/upper wall for this direction

  ! --- Fill in rays info ---
  ray%head  => ray
  ray%tail  => ray
  ray%ID = 1                                                            ! it's packet ID
  ray%n_rays    = pmax(a)*pmax(b)                                       ! number of rays in this direction
  ray%n_points  = pmax(c)                                               ! number of points per ray corresponds to n_cells in that direction

  allocate (ray%cell_coords(3,pmax(c),ray%n_rays))                      ! create ray table, containing x|y|z loc, for all n_cells and n_rays
  allocate (ray%ds(ray%n_rays,pmax(c)))                                 ! create ds table, containing ds for all points in all rays
  allocate (ray%fwd_coords(3,ray%n_rays))                               ! create table, containing x|y|z coords for upstream position for particular ray, for all rays
  allocate (ray%rev_coords(3,ray%n_rays))                               ! create table, containing x|y|z coords for downstream position for particular ray, for all rays
  allocate (ray%pt_wall(3,2,pmax(c), ray%n_rays))                       ! lower and upper walls the ray crosses in the cell

  ray_n = 1                                                             ! ray number
  do ia = self%li(a), self%ui(a)                                        ! iterate over cells in first axis
    do ib = self%li(b), self%ui(b)                                      ! iterate over cells in second axis
      p(a) = ia; p(b) = ib                                              ! point position in first two axes
      ray%fwd_coords(a,ray_n) = p(a)                                    ! upstream position in terms of first axis
      ray%fwd_coords(b,ray_n) = p(b)                                    ! upstream position in terms of second axis
      ray%fwd_coords(c,ray_n) = self%li(c) - 1                          ! upstream position in terms of ray direcion
      rp = 1
      do ic = self%li(c), self%ui(c)                                    ! move along ray direction
        p(c) = ic                                                       ! current cell position in terms of third axis
        ray%cell_coords(:,rp,ray_n) = p                                 ! cell position in particular ray
        ray%ds(ray_n,rp)            = self%cell_dims(c)                 ! ray patch length in particular cell
        ray%pt_wall(:,1,rp,ray_n) = bndl%lwall                          ! ray enters the cell through this wall
        ray%pt_wall(:,2,rp,ray_n) = bndl%uwall                          ! ray leaves the cell through this wall
        rp = rp+1                                                       ! increment ray point
      end do
      ray%rev_coords(a,ray_n) = p(a)                                    ! downstream position in terms of first axis
      ray%rev_coords(b,ray_n) = p(b)                                    ! downstream position in terms of second axis
      ray%rev_coords(c,ray_n) = self%ui(c) + 1                          ! downstream position in terms of ray direction
      ray_n = ray_n + 1                                                 ! move to next ray
    end do
  end do
  nullify(ndir, bndl, ray)
  call trace_end
END SUBROUTINE add_axis

!> =====================================================================
!> Subroutine to add rays (except 0 and 45 degree angles).
!> Done using modified "fast voxel traversal" algorithm
!> (original, not modified, one is from Amanatides & Woo, 1987)
!> =====================================================================
SUBROUTINE add_ray(self, theta, dir)
  class(ray_geometry_t) :: self
  real    :: theta, t, phi
  integer :: dir, nphi, a, b, c, pdir, round
  real, dimension(3) :: dp, init_step
  real, dimension(2) :: dp2d
  real, dimension(:), allocatable :: phis
  logical :: oob_dir
  ! ....................................................................
  call trace_begin('ray_geometry_t%add_ray', 4)

  nphi=floor(360./theta)                                                ! equispaced solid angles for all rays
  allocate(phis(nphi))
  phis(1) = 0; t = 360./nphi                                            ! phi = 0 direction

  do a = 2, nphi-1                                                      ! create phis array, avoiding phi = 360
    phis(a) = phis(a-1) + t
  end do

  phis(nphi) = phis(nphi-1) + t
  phis = phis/radian
  t = theta/radian

  SELECT case(dir)
  case(1)
    a = 2; b = 3; c = 1                                                 ! along x axis
  case(2)
    a = 3; b = 1; c = 2                                                 ! along y axis
  case(3)
    a = 1; b = 2; c = 3                                                 ! along z axis
  end SELECT

  do pdir = 1, nphi                                                     ! loop over phis
    phi = phis(pdir)
    self%n_ray_dirs = self%n_ray_dirs + 1                               ! increment ray dirs

    self%axis = dir
    self%theta = t
    self%phi = phis(pdir)

    dp(a) = sin(t)*cos(phi)*self%cell_dims(a)                           ! step size along a direction
    dp(b) = sin(t)*sin(phi)*self%cell_dims(b)                           ! step size along b direction
    dp(c) = cos(t)         *self%cell_dims(c)                           ! step size along c direction

    round = anint(dp(a)*100.0)                                          ! removing computer uncertainties
    dp(a) = round/100.0                                                 ! and giving angle error of around 1 degree
                                                                          ! this is necessary in order to avoid infinite loops
    round = anint(dp(b)*100.0)                                          ! and overabundance of ray points with machine epsilon
    dp(b) = round/100.0                                                 ! ds

    round = anint(dp(c)*100.0)
    dp(c) = round/100.0

    if (dp(a)==0) then                                                  ! if a direction == 0
      dp2d(1) = dp(b)                                                   ! prepare b direction for 2D traversal
      dp2d(2) = dp(c)                                                   ! prepare c direction for 2D traversal
      init_step(1) = 0.5                                                ! cell center for b direction
      init_step(2) = 1                                                  ! cell wall for c direction (ray propagation direction)
      init_step(3) = 0.5                                                ! passive cell step
      call self%traverse_2d(dp2d, b, c, a, init_step)                   ! proceed with 2D traversal

    else if (dp(b) == 0) then                                           ! if b direction == 0
      dp2d(1) = dp(a)                                                   ! prepare a direction for 2D traversal
      dp2d(2) = dp(c)                                                   ! prepare c direction for 2D traversal
      init_step(1) = 0.5                                                ! cell center for a direction
      init_step(2) = 1                                                  ! cell wall for c direction (ray propagation direction)
      init_step(3) = 0.5                                                ! passive cell step
      call self%traverse_2d(dp2d, a, c, b, init_step)                   ! proceed with 2D traversal

    else                                                                ! movement to all directions
      init_step(1) = 0.5                                                ! cell center for a direction
      init_step(2) = 0.5                                                ! cell center for b direction
      init_step(3) = 1                                                  ! cell wall for c direction (ray propagation direction)
      call self%traverse_3d (dp, a, b, c, init_step)                    ! proceed with 3D traversal
    end if
  end do

  deallocate(phis)
  call trace_end
END SUBROUTINE add_ray

!> =====================================================================
!> Subroutine to add diagonal (45 degree) rays. Done separately, because
!> it is a special case - rays overlap when done from 3 directions.
!> Otherwise it is still simpler and therefore faster than other angles -
!> ds are known apriori.
!> Same as add_ray for one direction and only 2 rays for second direction.
!> due to ray alignment, third direction is omitted.
!> In this setup rays DO NOT go through cell centres, therefore they DO NOT
!> cross cell corners - avoiding diagonal nbor problem. To set rays to go
!> through cell centers, move ray from cell center to wall and vice-versa
!> for 3D movement.
!> =====================================================================
SUBROUTINE add_diag(self, rep, dir)
  class(ray_geometry_t) :: self
  real, dimension(3) :: init_step, dp
  real, dimension(8) :: phis
  real, dimension(2) :: dp2d
  real    :: theta, t, phi
  integer :: rep, dir, nphi, a, b, c, pdir, round
! ......................................................................
  call trace_begin('ray_geometry_t%add_diag',4)

  t = 45./radian

  SELECT case(dir)
  case(1)
    a = 2; b = 3; c = 1                                                 ! along x axis
  case(2)
    a = 3; b = 1; c = 2                                                 ! along y axis
  case(3)
    a = 1; b = 2; c = 3                                                 ! along z axis
  end SELECT

  SELECT case(rep)                                                      ! 45 degree rays can be added only twice
  case(1)
    nphi = 8
    phis(1:8) = (/0., 45., 90., 135., 180., 225., 270., 315./)
  case(2)
    nphi = 2
    phis(1:2) = (/45., 315./)
  end SELECT

  phis = phis/radian
  do pdir = 1, nphi                                                     ! loop over phis
    self%n_ray_dirs = self%n_ray_dirs + 1                               ! increment ray dirs
    self%axis  = dir
    self%theta = t
    self%phi   = phis(pdir)
    phi=phis(dir)
    dp(a) = sin(t)*cos(phi)*self%cell_dims(a)                           ! step size along a direction
    dp(b) = sin(t)*sin(phi)*self%cell_dims(b)                           ! step size along b direction
    dp(c) = cos(t)         *self%cell_dims(c)                           ! step size along c direction

    round = anint(dp(a)*100.0)                                          ! removing computer uncertainties
    dp(a) = round/100.0                                                 ! and giving angle error of around 1 degree
    round = anint(dp(b)*100.0)                                          ! this is necessary in order to avoid infinite loops
    dp(b) = round/100.0                                                 ! and overabundance of ray points with machine epsilon
    round = anint(dp(c)*100.0)                                          ! ds
    dp(c) = round/100.0

    if (dp(a)==0) then                                                  ! if a direction == 0
      dp2d(1) = dp(b)                                                   ! prepare b direction for 2D traversal
      dp2d(2) = dp(c)                                                   ! prepare c direction for 2D traversal
      init_step(1) = 0.5                                                ! cell center for a direction
      init_step(2) = 1                                                  ! cell wall for c direction (ray propagation direction)
      init_step(3) = 0.5                                                ! passive cell step
      call self%traverse_2d(dp2d, b, c, a, init_step)                   ! proceed with 2D traversal

    else if (dp(b) == 0) then                                           ! if b direction == 0
      dp2d(1) = dp(a)                                                   ! prepare a direction for 2D traversal
      dp2d(2) = dp(c)                                                   ! prepare c direction for 2D traversal
      init_step(1) = 0.5                                                ! cell center for a direction
      init_step(2) = 1                                                  ! cell wall for c direction (ray propagation direction)
      init_step(3) = 0.5                                                ! passive cell step
      call self%traverse_2d(dp2d, a, c, b, init_step)                   ! proceed with 2D traversal

    else                                                                ! movement to all directions
      init_step(1) = 0.5                                                ! cell center for a direction
      init_step(2) = 0.5                                                ! cell center for b direction
      init_step(3) = 1                                                  ! cell wall for c direction (ray propagation direction)
      call self%traverse_3d (dp, a, b, c, init_step)                    ! proceed with 3D traversal
    end if

  end do
  call trace_end
END SUBROUTINE add_diag

!> =====================================================================
!> Subroutine to traverse rays in 2D. Similar to "fast voxel traversal
!> algorithm" (Amanatides & Woo, 1987) and appendix A in Rijkhorst et al. (2006),
!> but without their serious drawbacks and is (probably) faster.
!> input: self, step size inside the cells (dp), nominal direction (a),
!> ray direction (b), fixed coordinate (fc), starting position (init_step)
!> =====================================================================
SUBROUTINE traverse_2d(self, dp, a, b, fc, init_step)
  class(ray_geometry_t) :: self
  class(ray_point_t), pointer :: ray
  class(single_ray_bundle_t), pointer :: ray_bundle
  class(traverse_2d_t), pointer :: t
  integer :: a, b, fc, ia, ib
  real    :: dp(2), init_step(3), I_up(3), lw(3)
  logical :: first_point
! ......................................................................
call trace_begin('traverse_2d_t%traverse_2d',4)
  allocate(t)                                                           ! allocate traversal object
  allocate(ray_bundle)                                                  ! create ray tree to store rays

  call t%init(dp, self%ui, self%li, a, b, fc)                           ! initialize traversal object

  do ia = t%cell_min(1), t%cell_max(1), t%cell_step(1)                  ! loop over nominal direction
    do ib = t%cell_min(3), t%cell_max(3), t%cell_step(3)                ! loop over passive direction
      first_point = .true.
      allocate(ray)                                                     ! create new ray (point)
      call t%posit(ia, ib, init_step, I_up)                             ! position traversal object
      do while (.not.t%oob)                                             ! loop until ray reaches wall at ray direction
        call ray%add_ray_point(t%cell_pos)                              ! add a ray point to the ray
        if (first_point) then                                           ! if it is the starting point, link stuff
          ray%head%lwall = t%lw                                         ! ray enters the cell through this wall
          first_point = .false.
        else                                                            ! add it to the end
          ray%tail%lwall = t%lw                                         ! ray enters the cell through this wall
          ray%tail%prev%uwall = -ray%tail%lwall                         ! ray leaves the cell through this wall
        end if

        call t%hop                                                      ! hop to the next cell
        ray%tail%ds = t%ds                                              ! save distance traveled
        if (t%r_trace >= t%r_tot) then                                  ! if the ray is traced
          ray%tail%uwall = -t%lw
          call t%get_downstream(ray%head%rev)                           ! get ray downstream cell
          ray%head%fwd = I_up                                           ! save ray upstream cell
          call ray_bundle%finish_ray(ray)                               ! and save ray to the ray tree
          if (t%cell_pos(b) > t%cell_max(2)) then                       ! if cell wall along the ray direction is reached
            t%oob = .true.                                              ! this ray is done
          else                                                          ! otherwise
            allocate(ray)                                               ! get new ray
            first_point = .true.
            call t%restart(I_up)                                        ! and restart the nominal direction
          end if
        end if

      end do
    end do
  end do
  call self%finish_single_ray_list(ray_bundle)                          ! save the ray tree to ray_geometry
  deallocate(t)
  call trace_end
END SUBROUTINE traverse_2D

!> =====================================================================
!> Subroutine to traverse rays in 3D. Similar to "fast voxel traversal
!> algorithm" (Amanatides & Woo, 1987) and appendix A in Rijkhorst et al. (2006),
!> but without their serious drawbacks and is faster.
!> =====================================================================
SUBROUTINE traverse_3d(self, dp, a, b, c, init_step)
  class(ray_geometry_t)             :: self
  class(ray_point_t),       pointer :: ray
  class(single_ray_bundle_t), pointer :: ray_bundle
  class(traverse_3d_t),     pointer :: t
  integer            :: a, b, c, ia, ib
  real, dimension(3) :: init_step, dp, I_up
  logical :: first_point
! ......................................................................
  call trace_begin('traverse_3d_t%traverse_3d')
  allocate(t)                                                           ! allocate traversal object
  allocate(ray_bundle)                                                  ! create ray tree to store rays

  call t%init(dp, self%ui, self%li, a, b, c)                            ! initialize traversal object

  do ia = t%cell_min(1), t%cell_max(1), t%cell_step(1)                  ! loop over 1st direction
    do ib = t%cell_min(2), t%cell_max(2), t%cell_step(2)                ! loop over 2nd direction
      first_point = .true.
      allocate(ray)                                                     ! create new ray
      call t%posit(ia, ib, init_step, I_up)                             ! position traversal object
      do while(.not.t%oob)                                              ! loop until ray reaches wall at ray direction
        call ray%add_ray_point(t%cell_pos)                              ! add a ray point to the ray
        if (first_point) then
          ray%head%lwall = t%lw                                         ! ray enters the cell through this wall
          first_point = .false.
        else
          ray%tail%lwall = t%lw                                         ! ray enters the cell through this wall
          ray%tail%prev%uwall = -ray%tail%lwall                         ! ray leaves the cell through this wall
        end if
        call t%hop                                                      ! hop to the next cell
        ray%tail%ds = t%ds                                              ! save distance traveled

        if (t%r_trace >= t%r_tot) then                                  ! if the ray is fully traced
          ray%tail%uwall = -t%lw
          call t%get_downstream(ray%head%rev)                           ! get ray downstream cell
          ray%head%fwd = I_up                                           ! save ray upstream cell
          call ray_bundle%finish_ray(ray)                               ! and save ray to the ray tree
          if (t%cell_pos(c) >= t%cell_max(3)) then                      ! if cell wall along the ray direction is reached
              t%oob = .true.                                            ! this ray is done
          else                                                          ! otherwise
            allocate(ray)                                               ! get new ray
            call t%restart(I_up)                                        ! and restart the 1st or 2nd direction
            first_point = .true.
          end if
        end if

      end do
    end do
  end do
  call self%finish_single_ray_list(ray_bundle)                          ! save the ray tree to ray_geometry
  deallocate(t)
END SUBROUTINE traverse_3D

!> =====================================================================
!> INPUT: self, step size along ray direction (dp), smallest cell coord (li)
!> highest cell coord (ui), nominal direction (a), ray direction (b),
!> fixed coord (fc)
!> =====================================================================
SUBROUTINE init_traverse_2D(self, dp, ui, li, a, b, fc)
  class(traverse_2d_t)  :: self
  integer, dimension(3) :: li, ui
  real, dimension(2)    :: dp
  integer               :: a, b, fc
! ......................................................................
  call trace_begin('traverse_2d_t%init')
  self%a = a; self%b = b; self%fc = fc
  self%cell_max(1) = ui(a)                                              ! last cell point along a axis
  self%cell_max(2) = ui(b)                                              ! last cell point along b axis
  self%cell_max(3) = ui(fc)                                             ! last cell point along fixed coord axis
  if (dp(1) < 0) then                                                   ! if ray goes from "top" to "bottom"
    self%cell_min (1) =  self%cell_max(1)                               ! swap cell_max and cell_min
    self%cell_max (1) =  li(a)                                          ! cell_max now is first cell point
    self%cell_step(1) = -1                                              ! stepper now is negative
    dp(1) = -dp(1)                                                      ! step size always has to be positive
  else
    self%cell_step(1) = 1                                               ! step through one cell to another
    self%cell_min (1) = li(a)
  end if
  self%dstep = 1/dp                                                     ! distance the ray travels to go one step
  self%cell_min(2)  = li(b);  self%cell_step(2) = 1                     ! ray direction is always positive
  self%cell_min(3)  = li(fc); self%cell_step(3) = 1                     ! fixed coordinate
  call trace_end
END SUBROUTINE init_traverse_2D



!> =====================================================================
!> INPUT: self, position in nominal direction (ia) and passive direction (ib),
!> distance to cell walls (init_step), upstream cell position for particular
!> ray (upstream)
!> =====================================================================
SUBROUTINE posit_2D(self, ia, ib, init_step, upstream)
  class (traverse_2d_t) :: self
  real, dimension(3) :: init_step
  real, dimension(3) :: upstream
  integer               :: ia, ib
! ......................................................................
  self%cell_pos(self%a ) = ia                                           ! position in a axis
  self%cell_pos(self%b ) = self%cell_min(2)                             ! position in ray direction axis
  self%cell_pos(self%fc) = ib                                           ! position in passive axis
  self%lw = 0; self%uw = 0

  if (init_step(1) == 1) then                                           ! initial ray position in a axis
    self%ray_p_init(1) = ia - self%cell_step(1)/2                       ! it starts from the wall
  else
    self%ray_p_init(1) = ia                                             ! or at predetermined position
  end if

  if (init_step(2) == 1) then                                           ! initial ray position in ray direction
    self%ray_p_init(2) = self%cell_min(2) - self%cell_step(2)/2         ! it starts from the wall
  else
    self%ray_p_init(2) = self%cell_min(2)                               ! or at predetermined position
  end if

  self%ray_p_last(1) = self%cell_max(1) + self%cell_step(1)/2           ! ray ends at cell wall
  self%ray_p_last(2) = self%cell_max(2) + self%cell_step(2)/2           ! ray ends at cell wall

  if (self%cell_step(1) > 0) then                                       ! distance ray travels in a direction
    self%ray_dist(1) = self%ray_p_last(1) - self%ray_p_init(1)
  else
    self%ray_dist(1) = self%ray_p_init(1) - self%ray_p_last(1)
  end if

  self%ray_dist(2) = self%ray_p_last(2) - self%ray_p_init(2)            ! distance ray travels in ray direction axis
  self%ray_dist    = self%ray_dist * self%dstep                         ! total distance ray travels in each direction to reach the end of patch
  self%r_tot       = minval(self%ray_dist)                              ! determine the smallest distance
  self%ray_step(1) = init_step(1)
  self%ray_step(2) = init_step(2)
  self%ray_step    = self%ray_step * self%dstep                         ! distance along the ray to cross to another cell

  if (init_step(1) == init_step(2)) then                                ! find upstream cell position
    upstream(self%a)  = self%cell_pos(self%a ) - self%cell_step(1)
    upstream(self%b)  = self%cell_pos(self%b ) - self%cell_step(2)
    upstream(self%fc) = self%cell_pos(self%fc)
    if (self%cell_step(1) > 0) then
      self%lw(self%a) = -1
    else
      self%lw(self%a) = 1
    end if
    if (self%cell_step(2) > 0) then
      self%lw(self%b) = -1
    else
      self%lw(self%b) = 1
    end if
  else if (init_step(1) > init_step(2)) then
    upstream(self%a)  = self%cell_pos(self%a ) - self%cell_step(1)
    upstream(self%b)  = self%cell_pos(self%b )
    upstream(self%fc) = self%cell_pos(self%fc)
    if (self%cell_step(1) > 0) then
      self%lw(self%a) = -1
    else
      self%lw(self%a) = 1
    end if
  else
    upstream(self%a)  = self%cell_pos(self%a )
    upstream(self%b)  = self%cell_pos(self%b ) - self%cell_step(2)
    upstream(self%fc) = self%cell_pos(self%fc)
    if (self%cell_step(2) > 0) then
      self%lw(self%b) = -1
    else
      self%lw(self%b) = 1
    end if
  end if
  self%lw(self%fc) = 0
  self%oob = .false.
  self%r_trace = 0.                                                     ! distance, the ray has traveled
END SUBROUTINE posit_2D

!> =====================================================================
!> Subroutine, which handles moving through cells in 2D
!> =====================================================================
SUBROUTINE hop_2D (self)
  class (traverse_2d_t) :: self
! ......................................................................
  associate (dstep => self%dstep, ray_step => self%ray_step, &
  cell_step => self%cell_step, cell_pos => self%cell_pos, &
  r_trace => self%r_trace, ds => self%ds, a => self%a, b => self%b)
  self%lw = 0
  if (ray_step(1) < ray_step(2)) then                                   ! if cell wall in direction 1 is closer than in
    cell_pos(a) = cell_pos(a) + cell_step(1)                            ! direction 2, hop to cell towards direction 1
    ray_step(2) = ray_step(2) - ray_step(1)                             ! reduce the distance to wall in 2 with distance traveled
             ds = ray_step(1)                                           ! save distance traveled
    r_trace     = r_trace + ray_step(1)                                 ! save total distance traveled (for sanity check)
    ray_step(1) = dstep(1)                                              ! reset distance to wall for 1
    self%lw(a)  = -1
  else if (ray_step(1) == ray_step(2)) then                             ! if both equal, move diagonally
    cell_pos(a) = cell_pos(a) + cell_step(1)                            ! cell position is diagonal to the old one
    cell_pos(b) = cell_pos(b) + cell_step(2)
             ds = ray_step(1)                                           ! distance of both 1 and 2 is the same
    r_trace     = r_trace + ray_step(1)                                 ! save total distance traveled
    ray_step    = dstep                                                 ! reset distance to wall for 1 and 2
    self%lw(a)  = -1
    self%lw(b)  = -1
  else
    cell_pos(b) = cell_pos(b) + cell_step(2)                            ! if wall is closer in 2, move towards there
    ray_step(1) = ray_step(1) - ray_step(2)                             ! reduce the distance to 1
             ds = ray_step(2)                                           ! save distance traveled
    r_trace     = r_trace + ray_step(2)                                 ! save total distance traveled
    ray_step(2) = dstep(2)                                              ! reset distance to wall for 2
    self%lw(b)  = -1
  end if
  end associate
END SUBROUTINE hop_2D

!> =====================================================================
!> Get downstream coordinates for 2D case
!> =====================================================================
SUBROUTINE get_downstream_2D (self, downstream)
  class (traverse_2d_t) :: self
  real                  :: r_trace
  real, dimension(2)    :: ray_step
  integer, dimension(3) :: cell_pos, downstream
! ......................................................................
  ray_step = self%ray_step                                              ! save current ray steps
  cell_pos = self%cell_pos                                              ! save current cell position
  r_trace  = self%r_trace                                               ! save current distance
  self%uw = 0
  call hop_2D(self)                                                     ! make a hop
  downstream = self%cell_pos                                            ! that is your downstream cell
  ! determine which wall is crossed
  if (cell_pos(self%a) /= downstream(self%a)) then
    if (self%cell_step(1) > 0) then
      self%uw(self%a) = 1
    else
      self%uw(self%a) = -1
    end if
  else
    self%uw(self%a) = 0
  end if

  if (cell_pos(self%b) /= downstream(self%b)) then
    if (self%cell_step(2) > 0) then
      self%uw(self%b) = 1
    else
      self%uw(self%b) = -1
    end if
  else
    self%uw(self%b) = 0
  end if
  self%uw(self%fc) = 0

  self%cell_pos = cell_pos                                              ! reset cell position
  self%ray_step = ray_step                                              ! reset ray steps
  self%r_trace  = r_trace                                               ! reset total distance
END SUBROUTINE get_downstream_2D

!> =====================================================================
!> INPUT: self, step size along ray direction (dp), smallest cell coord (li)
!> highest cell coord (ui), nominal direction (a, b), ray direction (c)
!> =====================================================================
SUBROUTINE init_traverse_3D(self, dp, ui, li, a, b, c)
  class(traverse_3d_t) :: self
  integer, dimension(3)   :: li, ui
  real, dimension(3)   :: dp
  integer              :: a, b, c
! ......................................................................
  call trace_begin('traverse_3d_t%init')
  self%a = a; self%b = b; self%c = c                                    ! convenience
  self%cell_max(1) = ui(a)                                              ! last cell point along a axis
  self%cell_max(2) = ui(b)                                              ! last cell point along b axis
  self%cell_max(3) = ui(c)                                              ! last cell point along c axis
  if (dp(a) < 0) then                                                   ! if ray goes from "top" to "bottom"
    self%cell_min (1) =  self%cell_max(1)                               ! swap cell_max and cell_min
    self%cell_max (1) =  li(a)                                          ! cell_max now is first cell point
    self%cell_step(1) = -1                                              ! stepper now is negative
    dp(a) = -dp(a)                                                      ! step size always has to be positive
  else
    self%cell_step(1) = 1                                               ! step through one cell to another
    self%cell_min (1) = li(a)
  end if
  if (dp(b) < 0) then                                                   ! if ray goes from "top" to "bottom"
    self%cell_min (2) = self%cell_max(2)                                ! swap cell_max and cell_min
    self%cell_max (2) = li(b)                                           ! cell_max now is first cell point
    self%cell_step(2) = -1                                              ! stepper now is negative
    dp(b)        = -dp(b)                                               ! step size always has to be positive
  else
    self%cell_step(2) = 1
    self%cell_min (2) = li(b)
  end if
  self%dstep(1) = 1/dp(a)                                               ! distance the ray travels to go one step
  self%dstep(2) = 1/dp(b)                                               ! distance the ray travels to go one step
  self%dstep(3) = 1/dp(c)                                               ! distance the ray travels to go one step
  self%cell_min(3)  = li(c); self%cell_step(3) = 1                      ! ray direction is always positive
  call trace_end
END SUBROUTINE init_traverse_3D

!> =====================================================================
!> INPUT: self, position in nominal directions (ia, ib),
!> distance to cell walls (init_step), upstream cell position for particular
!> ray (upstream)
!> =====================================================================
SUBROUTINE posit_3D(self, ia, ib, init_step, upstream)
  class (traverse_3d_t) :: self
  real, dimension(3) :: init_step
  real, dimension(3) :: upstream
  integer :: ia, ib
! ......................................................................
  self%cell_pos(self%a) = ia                                            ! position in a axis
  self%cell_pos(self%b) = ib                                            ! position in b axis
  self%cell_pos(self%c) = self%cell_min(3)                              ! starting position for c axis
  self%lw = 0; self%uw = 0
  if (init_step(1) == 1) then                                           ! determine whether ray starts
    self%ray_p_init(1) = ia - self%cell_step(1)/2                       ! at the wall or at the centre
  else                                                                  ! in terms of thjs direction
    self%ray_p_init(1) = ia
  end if

  if (init_step(2) == 1) then
    self%ray_p_init(2) = ib - self%cell_step(2)/2
  else
    self%ray_p_init(2) = ib
  end if

  if (init_step(3) == 1) then
    self%ray_p_init(3) = self%cell_min(3) - self%cell_step(3)/2
  else
    self%ray_p_init(3) = self%cell_min(3)
  end if


  self%ray_p_last = self%cell_max + self%cell_step/2                    ! ray finishes at the wall

  if (self%cell_step(1) > 0) then                                       ! estimate ray distance for a dirrection
    self%ray_dist(1) = self%ray_p_last(1) - self%ray_p_init(1)
  else
    self%ray_dist(1) = self%ray_p_init(1) - self%ray_p_last(1)
  end if

  if (self%cell_step(2) > 0) then                                       ! estimate ray distance for b dirrection
    self%ray_dist(2) = self%ray_p_last(2) - self%ray_p_init(2)
  else
    self%ray_dist(2) = self%ray_p_init(2) - self%ray_p_last(2)
  end if

  self%ray_dist(3) = self%ray_p_last(3) - self%ray_p_init(3)            ! estimate ray distance for c dirrection
  self%ray_dist    = self%ray_dist * self%dstep
  self%r_tot       = minval(self%ray_dist)                              ! determine the shortest path

  if (self%r_tot< 0) then                                               ! never happens, but just in case
    error stop "ray_geometry:posit_3d: crap"
  end if

  self%ray_step    = init_step                                          ! set initial ray steps
  self%ray_step    = self%ray_step * self%dstep                         ! normalise in axis point of view

  if (init_step(1) == init_step(2)) then                                ! determine upstream cell position
    if  (init_step(2) == init_step(3)) then
      upstream(self%a) = self%cell_pos(self%a) - self%cell_step(1)
      upstream(self%b) = self%cell_pos(self%b) - self%cell_step(2)
      upstream(self%c) = self%cell_pos(self%c) - self%cell_step(3)
      call self%wall(1); call self%wall(2); call self%wall(3)           ! and the lower face

    else if (init_step(2) < init_step(3)) then
      upstream(self%a) = self%cell_pos(self%a)
      upstream(self%b) = self%cell_pos(self%b)
      upstream(self%c) = self%cell_pos(self%c) - self%cell_step(3)
      call self%wall(1,0); call self%wall(2,0); call self%wall(3)
    else
      upstream(self%a) = self%cell_pos(self%a) - self%cell_step(1)
      upstream(self%b) = self%cell_pos(self%b) - self%cell_step(2)
      upstream(self%c) = self%cell_pos(self%c)
      call self%wall(1); call self%wall(2); call self%wall(3,0)
    end if
  else if (init_step(2) == init_step(3)) then
    if (init_step(1) < init_step(3)) then
      upstream(self%a) = self%cell_pos(self%a)
      upstream(self%b) = self%cell_pos(self%b) - self%cell_step(2)
      upstream(self%c) = self%cell_pos(self%c) - self%cell_step(3)
      call self%wall(1,0); call self%wall(2); call self%wall(3)
    else
      upstream(self%a) = self%cell_pos(self%a) - self%cell_step(1)
      upstream(self%b) = self%cell_pos(self%b)
      upstream(self%c) = self%cell_pos(self%c)
      call self%wall(1); call self%wall(2,0); call self%wall(3,0)
    end if
  else if (init_step(1) == init_step(3)) then
    if (init_step(1) < init_step(2)) then
      upstream(self%a) = self%cell_pos(self%a)
      upstream(self%b) = self%cell_pos(self%b) - self%cell_step(2)
      upstream(self%c) = self%cell_pos(self%c)
      call self%wall(1,0); call self%wall(2); call self%wall(3,0)
    else
      upstream(self%a) = self%cell_pos(self%a) - self%cell_step(1)
      upstream(self%b) = self%cell_pos(self%b)
      upstream(self%c) = self%cell_pos(self%c) - self%cell_step(3)
      call self%wall(1); call self%wall(2,0); call self%wall(3)
    end if
  else if (init_step(1) < init_step(3)) then
    if (init_step(2) < init_step(3)) then
      upstream(self%a) = self%cell_pos(self%a)
      upstream(self%b) = self%cell_pos(self%b)
      upstream(self%c) = self%cell_pos(self%c) - self%cell_step(3)
      call self%wall(1,0); call self%wall(2,0); call self%wall(3)
    else
      upstream(self%a) = self%cell_pos(self%a)
      upstream(self%b) = self%cell_pos(self%b) -  self%cell_step(2)
      upstream(self%c) = self%cell_pos(self%c)
      call self%wall(1,0); call self%wall(2); call self%wall(3,0)
    end if
  else if (init_step(1) > init_step(2)) then
    upstream(self%a) = self%cell_pos(self%a) - self%cell_step(1)
    upstream(self%b) = self%cell_pos(self%b)
    upstream(self%c) = self%cell_pos(self%c)
    call self%wall(1); call self%wall(2,0); call self%wall(3,0)
  else
    upstream(self%a) = self%cell_pos(self%a)
    upstream(self%b) = self%cell_pos(self%b) -  self%cell_step(2)
    upstream(self%c) = self%cell_pos(self%c)
    call self%wall(1,0); call self%wall(2); call self%wall(3,0)
  end if
  self%oob     = .false.
  self%reset_a = .false.
  self%reset_b = .false.
  self%r_trace = 0.
END SUBROUTINE posit_3D

!> =====================================================================
!> Convenience. Sets the lower wall index. it is self%wall
!> =====================================================================
SUBROUTINE cw_3d(self, axis, zero)
  class (traverse_3d_t) :: self
  integer               :: axis, c
  integer, optional     :: zero
  !---------------------------------------------------------------------
  SELECT case(axis)
  case(1)
    c = self%a
  case(2)
    c = self%b
  case(3)
    c = self%c
  end SELECT
  if(present(zero)) then
    self%lw(c) = 0
  else
    if (self%cell_step(axis) > 0) then
      self%lw(c) = -1
    else
      self%lw(c) = 1
    end if
  end if
END SUBROUTINE cw_3d

!> =====================================================================
!> Subroutine, which handles moving through cells in 3D
!> =====================================================================
SUBROUTINE hop_3D (self)
  class (traverse_3d_t) :: self
! ......................................................................
  associate (dstep => self%dstep, ray_step => self%ray_step, &
  cell_step => self%cell_step, cell_pos => self%cell_pos, &
  r_trace => self%r_trace, ds => self%ds, a => self%a, b => self%b, c => self%c)
  if (ray_step(1) == ray_step(2)) then
    if (ray_step(2) == ray_step(3)) then
      cell_pos(a) = cell_pos(a) + cell_step(1)
      cell_pos(b) = cell_pos(b) + cell_step(2)
      cell_pos(c) = cell_pos(c) + cell_step(3)
               ds = ray_step(1)
      r_trace     = r_trace + ray_step(1)
      ray_step    = dstep
    else if (ray_step(2) < ray_step(3)) then
      cell_pos(a) = cell_pos(a) + cell_step(1)
      cell_pos(b) = cell_pos(b) + cell_step(2)
               ds = ray_step(1)
      ray_step(3) = ray_step(3) - ray_step(1)
      r_trace     = r_trace + ray_step(1)
      ray_step(1) = dstep(1)
      ray_step(2) = dstep(2)
    else
      cell_pos(c) = cell_pos(c) + cell_step(3)
               ds = ray_step(3)
      ray_step(1) = ray_step(1) - ray_step(3)
      ray_step(2) = ray_step(2) - ray_step(3)
      r_trace     = r_trace + ray_step(3)
      ray_step(3) = dstep(3)
    end if
  else if (ray_step(2) == ray_step(3)) then
    if (ray_step(1) < ray_step(3)) then
      cell_pos(a) = cell_pos(a) + cell_step(1)
               ds = ray_step(1)
      ray_step(2) = ray_step(2) - ray_step(1)
      ray_step(3) = ray_step(3) - ray_step(1)
      r_trace     = r_trace + ray_step(1)
      ray_step(1) = dstep(1)
    else
      cell_pos(b) = cell_pos(b) + cell_step(2)
      cell_pos(c) = cell_pos(c) + cell_step(3)
               ds = ray_step(2)
      ray_step(1) = ray_step(1) - ray_step(2)
      r_trace     = r_trace + ray_step(2)
      ray_step(2) = dstep(2)
      ray_step(3) = dstep(3)
    end if
  else if (ray_step(1) == ray_step(3)) then
    if (ray_step(2) < ray_step(3)) then
      cell_pos(b) = cell_pos(b) + cell_step(2)
               ds = ray_step(2)
      ray_step(1) = ray_step(1) - ray_step(2)
      ray_step(3) = ray_step(3) - ray_step(2)
      r_trace     = r_trace + ray_step(2)
      ray_step(2) = dstep(2)
    else
      cell_pos(a) = cell_pos(a) + cell_step(1)
      cell_pos(c) = cell_pos(c) + cell_step(3)
               ds = ray_step(1)
      ray_step(2) = ray_step(2) - ray_step(1)
      r_trace     = r_trace + ray_step(1)
      ray_step(1) = dstep(1)
      ray_step(3) = dstep(3)
    end if
  else if (ray_step(1) < ray_step(2)) then
    if (ray_step(1) < ray_step(3)) then
      cell_pos(a) = cell_pos(a) + cell_step(1)
               ds = ray_step(1)
      ray_step(2) = ray_step(2) - ray_step(1)
      ray_step(3) = ray_step(3) - ray_step(1)
      r_trace     = r_trace + ray_step(1)
      ray_step(1) = dstep(1)
    else
      cell_pos(c) = cell_pos(c) + cell_step(3)
               ds = ray_step(3)
      ray_step(1) = ray_step(1) - ray_step(3)
      ray_step(2) = ray_step(2) - ray_step(3)
      r_trace     = r_trace + ray_step(3)
      ray_step(3) = dstep(3)
    end if
  else if (ray_step(2) < ray_step(3)) then
    cell_pos(b) = cell_pos(b) + cell_step(2)
             ds = ray_step(2)
    ray_step(1) = ray_step(1) - ray_step(2)
    ray_step(3) = ray_step(3) - ray_step(2)
    r_trace     = r_trace + ray_step(2)
    ray_step(2) = dstep(2)
  else
    cell_pos(c) = cell_pos(c) + cell_step(3)
             ds = ray_step(3)
    ray_step(1) = ray_step(1) - ray_step(3)
    ray_step(2) = ray_step(2) - ray_step(3)
    r_trace     = r_trace + ray_step(3)
    ray_step(3) = dstep(3)
  end if
  end associate
END SUBROUTINE hop_3D

!> =====================================================================
!> Obtain downstream cell position
!> =====================================================================
SUBROUTINE get_downstream_3D (self, downstream)
  class (traverse_3d_t) :: self
  real                  :: r_trace
  integer, dimension(3) :: cell_pos, downstream, ray_step
! ......................................................................
  ray_step = self%ray_step
  cell_pos = self%cell_pos
  r_trace  = self%r_trace

  call hop_3D(self)
  downstream = self%cell_pos
  self%uw = 0
  if (cell_pos(self%a) /= downstream(self%a)) then
    if (self%cell_step(1) > 0) then
      self%uw(self%a) =  1
    else
      self%uw(self%a) = -1
    end if
  else
    self%uw(self%a) = 0
  end if

  if (cell_pos(self%b) /= downstream(self%b)) then
    if (self%cell_step(2) > 0) then
      self%uw(self%b) =  1
    else
      self%uw(self%b) = -1
    end if
  else
    self%uw(self%b) = 0
  end if

  if (cell_pos(self%c) /= downstream(self%c)) then
    if (self%cell_step(3) > 0) then
      self%uw(self%c) =  1
    else
      self%uw(self%c) = -1
    end if
  else
    self%uw(self%c) = 0
  end if

  self%cell_pos = cell_pos
  self%ray_step = ray_step
  self%r_trace  = r_trace
END SUBROUTINE get_downstream_3D

!> =====================================================================
!> restart the ray when it hits a side-wall (invoke pseudo-periodic
!> boundaries). This is for 2D case (less steps needed than in 3D ->
!> is faster
!> =====================================================================
SUBROUTINE restart_2D (self, upstream)
  class(traverse_2d_t) :: self
  real, dimension(3)   :: upstream
! ......................................................................
  self%ray_p_init(1)      = self%cell_min(1)   - self%cell_step(1)/2    ! ray starts at the wall
  self%ray_p_last(1)      = self%cell_max(1)   + self%cell_step(1)/2    ! ray ends at the wall
  if (self%cell_step(1) > 0) then                                       ! scheck ray direction
    self%ray_dist(1)      = self%ray_p_last(1) - self%ray_p_init(1)     ! and calculate distance to the wall
  else
    self%ray_dist(1)      = self%ray_p_init(1) - self%ray_p_last(1)
  end if
  self%ray_dist(1)        = self%ray_dist(1)   * self%dstep(1)          ! determine new ray distance
  self%ray_dist(2)        = self%ray_dist(2)   - self%r_trace           ! but save the distance, traveled in ray direction
  self%r_tot              = minval(self%ray_dist)                       ! find the minimal value
  self%r_trace            = 0                                           ! restart the ray total distance
  self%cell_pos(self%a)   = self%cell_min(1)                            ! restart the cell position in nominal axis
  self%ray_step(1)        = self%dstep(1)                               ! restart the distance, the ray has to travel in nominal axis

  upstream(self%a)  = self%cell_pos(self%a) - self%cell_step(1)         ! get upstream cell position
  upstream(self%b)  = self%cell_pos(self%b)
  upstream(self%fc) = self%cell_pos(self%fc)
  self%uw = 0
  self%lw(self%a)  = -1; self%lw(self%b)  =  0; self%lw(self%fc) =  0   ! save wall data
END SUBROUTINE restart_2D

!> =====================================================================
!> restart the ray when it hits a side-wall (invoke pseudo-periodic
!> boundaries). This is for the general, 3D case.
!> =====================================================================
SUBROUTINE restart_3D (self, upstream)
  class(traverse_3d_t) :: self
  real, dimension(3), intent(inout) :: upstream
! ......................................................................
  if (((self%cell_pos(self%a) >= self%cell_max(1)).and.(self%cell_step(1)>0)).or.(&
  (self%cell_pos(self%a) <= self%cell_max(1)).and.(self%cell_step(1)<0))) then
    self%reset_a = .true.
  end if
  if (((self%cell_pos(self%b) >= self%cell_max(2)).and.(self%cell_step(2)>0)).or.(&
     (self%cell_pos(self%b) <= self%cell_max(2)).and.(self%cell_step(2)<0))) then
    self%reset_b = .true.
  end if
  if((self%reset_a).and.(.not.self%reset_b)) then
    self%ray_p_init(1)    = self%cell_min(1)   - self%cell_step(1)/2
    self%ray_p_last(1)    = self%cell_max(1)   + self%cell_step(1)/2
    if (self%cell_step(1) > 0) then
      self%ray_dist(1)    = self%ray_p_last(1) - self%ray_p_init(1)
    else
      self%ray_dist(1)    = self%ray_p_init(1) - self%ray_p_last(1)
    end if
    self%ray_dist(1)      = self%ray_dist(1)   * self%dstep(1)
    self%ray_dist(2)      = self%ray_dist(2)   - self%r_trace
    self%ray_dist(3)      = self%ray_dist(3)   - self%r_trace
    self%r_tot            = minval(self%ray_dist)
    self%r_trace          = 0
    self%cell_pos(self%a) = self%cell_min(1)
    self%ray_step(1)      = self%dstep(1)
    upstream(self%a)      = self%cell_pos(self%a) - self%cell_step(1)
    upstream(self%b)      = self%cell_pos(self%b)
    upstream(self%c)      = self%cell_pos(self%c)
    self%lw(self%a)  = -1; self%lw(self%b)  =  0; self%lw(self%c) =  0

  else if ((self%reset_b).and.(.not.self%reset_a)) then
    self%ray_p_init(2)    = self%cell_min(2)      - self%cell_step(2)/2
    self%ray_p_last(2)    = self%cell_max(2)      + self%cell_step(2)/2
    if (self%cell_step(2) > 0) then
      self%ray_dist(2)    = self%ray_p_last(2)    - self%ray_p_init(2)
    else
      self%ray_dist(2)    = self%ray_p_init(2)    - self%ray_p_last(2)
    end if
    self%ray_dist(2)      = self%ray_dist(2)      * self%dstep(2)
    self%ray_dist(1)      = self%ray_dist(1)      - self%r_trace
    self%ray_dist(3)      = self%ray_dist(3)      - self%r_trace
    self%r_tot            = minval(self%ray_dist)
    self%r_trace          = 0
    self%cell_pos(self%b) = self%cell_min(2)
    self%ray_step(2)      = self%dstep(2)
    upstream(self%a)      = self%cell_pos(self%a)
    upstream(self%b)      = self%cell_pos(self%b) - self%cell_step(2)
    upstream(self%c)      = self%cell_pos(self%c)
    self%lw(self%a)  = 0; self%lw(self%b)  =  -1; self%lw(self%c) =  0
  else if ((self%reset_a).and.(self%reset_b)) then
    self%ray_p_init(1)    = self%cell_min(1)      - self%cell_step(1)/2
    self%ray_p_last(1)    = self%cell_max(1)      + self%cell_step(1)/2
    self%ray_p_init(2)    = self%cell_min(2)      - self%cell_step(2)/2
    self%ray_p_last(2)    = self%cell_max(2)      + self%cell_step(2)/2
    if (self%cell_step(1) > 0) then
      self%ray_dist(1)    = self%ray_p_last(1)    - self%ray_p_init(1)
    else
      self%ray_dist(1)    = self%ray_p_init(1)    - self%ray_p_last(1)
    end if
    if (self%cell_step(2) > 0) then
      self%ray_dist(2)    = self%ray_p_last(2)    - self%ray_p_init(2)
    else
      self%ray_dist(2)    = self%ray_p_init(2)    - self%ray_p_last(2)
    end if
    self%ray_dist(1)      = self%ray_dist(1)      * self%dstep(1)
    self%ray_dist(2)      = self%ray_dist(2)      * self%dstep(2)
    self%ray_dist(3)      = self%ray_dist(3)      - self%r_trace
    self%r_tot            = minval(self%ray_dist)
    self%r_trace          = 0
    self%cell_pos(self%a) = self%cell_min(1)
    self%cell_pos(self%b) = self%cell_min(2)
    self%ray_step(1)      = self%dstep(1)
    self%ray_step(2)      = self%dstep(2)
    upstream(self%a)      = self%cell_pos(self%a) - self%cell_step(1)
    upstream(self%b)      = self%cell_pos(self%b) - self%cell_step(2)
    upstream(self%c)      = self%cell_pos(self%c)
    self%lw(self%a)  = -1; self%lw(self%b)  =  -1; self%lw(self%c) =  0
  else
    error stop "ray_geometry: in reset 3D resets went haywire"
  end if
  self%uw = 0
END SUBROUTINE restart_3D

!> =====================================================================
!> Subroutine to create linked-list of ray points. Later this list is
!> redone into a stiff (array) ray.
!> INPUT: self, cell coordinates where the ray point have moved to (p)
!> =====================================================================
SUBROUTINE add_ray_point(self, p)
  class(ray_point_t)          :: self
  class(ray_point_t), pointer :: point
  integer :: p(3)
! ......................................................................
  call trace_begin('ray_geometry_t%add_ray_point')
  allocate(point)
  if (associated(self%tail)) then
    point%prev => self%tail
    point%prev%next => point
    point%head => self%head
    self%tail => point
    self%head%n_points = self%head%n_points + 1                         ! increment number of points in this ray
  else
    self%head => point
    self%tail => point
    self%head%n_points = 1
  end if
  point%coords   = p                                                    ! set coordinates
  call trace_end
END SUBROUTINE add_ray_point

!> =====================================================================
!> Subroutine to transform linked-list of ray points into a stiff (array)
!> ray. Arranging is done based on the walls the ray hits AND the number
!> of ray points.
!> INPUT: self, linked list of ray points (ray_point)
!> =====================================================================

SUBROUTINE fill_single_ray_data(self, ray)
  class(single_ray_t) :: self
  class(ray_point_t), pointer :: ray, prev
  integer :: pt
! ......................................................................
  prev => ray%tail
  ray  => ray%head

  self%lwall    = ray%lwall                                          ! store lower wall
  self%uwall    = prev%uwall                                            ! store upper wall
  self%n_points = ray%n_points                                          ! store number of points in this ray

  allocate(self%cell_coords(3, ray%n_points))                           ! create a stiff array for this ray with cell cords
  allocate(self%ds (ray%n_points))                                      ! create a stiff array for this ray with ds
  allocate(self%cell_walls (3,2,ray%n_points))                          ! create a stiff array for this ray with cell uw and lw
  self%fwd_coords = ray%fwd                                             ! save upstream cell coords
  self%rev_coords = ray%rev                                             ! save downstream cell coords

  pt=1
  do while (associated(ray))                                            ! loop over ray points
    self%cell_coords(:,pt)  = ray%coords                                ! save ray point coords
    self%ds(pt)             = ray%ds                                    ! save ray point ds
    self%cell_walls(:,1,pt) = ray%lwall                                 ! save ray entry wall
    self%cell_walls(:,2,pt) = ray%uwall                                 ! save ray exit wall
    prev => ray                                                         ! save previous point
    ray => ray%next                                                     ! go to next point
    deallocate(prev)                                                    ! deallocate previous point
    pt = pt+1
  end do
END SUBROUTINE fill_single_ray_data

SUBROUTINE finish_ray(self, ray_point)
  class(single_ray_bundle_t) :: self
  class(single_ray_bundle_t), pointer :: bundle
  class(single_ray_packet_t), pointer :: packet
  class(single_ray_t), pointer  :: ray
  class(ray_point_t),  pointer  :: ray_point
  logical :: fnd_p, fnd_b
! ......................................................................
  call trace_begin('ray_geometry_t%finish_ray')

  allocate (ray)                                                        ! allocate ray list
  call ray%fill_single_ray_data(ray_point)

  if (.not.associated(self%head)) then                                  ! if the ray tree is empty
    allocate(bundle)
    self%head => bundle
    self%tail => bundle
    bundle%lwall = ray%lwall
    bundle%uwall = ray%uwall
  else
    bundle => self%head
    fnd_b = .false.
    do while (associated(bundle).and..not.fnd_b)
      if (all(bundle%lwall == ray%lwall).and. &                         ! bundles are defined by which walls
          all(bundle%uwall == ray%uwall)) then                          ! the ray enters and leaves the patch
        fnd_b = .true.
      else
        bundle => bundle%next
      end if
    end do
    if(.not.fnd_b) then
      allocate(bundle)
      bundle%head => self%head
      bundle%prev => self%tail
      bundle%prev%next => bundle
      self%tail => bundle
      bundle%lwall = ray%lwall
      bundle%uwall = ray%uwall
    end if
  end if

  if (.not.associated(bundle%packets)) then
    allocate(packet)
    bundle%packets => packet
    packet%head => packet
    packet%tail => packet
    packet%n_rays = 1
    packet%n_points = ray%n_points
  else
    packet => bundle%packets
    fnd_p = .false.
    do while (associated(packet).and..not.fnd_p)
      if (ray%n_points == packet%n_points) then                         ! rays with same number of points are packed together
        fnd_p = .true.
        packet%n_rays = packet%n_rays + 1
      else
        packet => packet%next
      end if
    end do
    if(.not.fnd_p) then
      allocate(packet)
      packet%head => bundle%packets
      packet%prev => bundle%packets%tail
      packet%prev%next => packet
      bundle%packets%tail => packet
      packet%n_rays = 1
      packet%n_points = ray%n_points
    end if
  end if

  if (.not.associated(packet%rays)) then
    packet%rays => ray
    ray%head => ray
    ray%tail => ray
  else
    ray%head => packet%rays
    ray%prev => packet%rays%tail
    ray%prev%next => ray
    packet%rays%tail => ray
  end if

  nullify (ray, bundle, packet, ray_point)
  call trace_end
END SUBROUTINE finish_ray

SUBROUTINE fill_ray_list_data(self, sray)
  class(ray_t) :: self
  class(single_ray_t), pointer :: sray, prev_sray
  integer :: iray
! ......................................................................
  allocate (self%cell_coords(3,self%n_points,self%n_rays))              ! cell coords per point per ray in the bundle
  allocate (self%ds(self%n_rays,self%n_points))                         ! ds per point per ray in the bundle
  allocate (self%fwd_coords(3,self%n_rays))                             ! upstream cell coords per ray
  allocate (self%rev_coords(3,self%n_rays))                             ! downstream cell coords per ray
  allocate (self%pt_wall(3,2,self%n_points,self%n_rays))                ! wall info per point per ray

  sray => sray%head
  iray = 0
  do while (associated(sray))                                           ! loop over all rays in the tree
    iray = iray + 1
    self%cell_coords(:,:,iray) = sray%cell_coords                       ! save cell cords
    self%ds(iray,:)            = sray%ds(:)                             ! save distances, traveled in those cells
    self%pt_wall(:,:,:, iray)  = sray%cell_walls                        ! save wall crossing info per point per ray
    self%fwd_coords(:,iray)    = sray%fwd_coords                        ! save upstream coords
    self%rev_coords(:,iray)    = sray%rev_coords                        ! save downstream coords
    prev_sray   => sray                                                 ! hold pointer to this ray
    deallocate(sray%cell_coords, sray%ds, sray%cell_walls)
    sray => sray%next                                                   ! go to next ray
    deallocate (prev_sray)                                              ! deallocate previous ray
  end do
  nullify (sray, prev_sray)
END SUBROUTINE fill_ray_list_data

!> =====================================================================
!> Subroutine to restructure ray tree into 2D sorted ray tables.
!> =====================================================================
SUBROUTINE finish_single_ray_list(self, sbundle)
  class(ray_geometry_t) :: self
  class(ray_dirs_t), pointer :: dir
  class(ray_bndl_t), pointer :: bundle
  class(ray_t)     , pointer :: ray
  ! ---
  class(single_ray_bundle_t), pointer :: sbundle, sb_prev
  class(single_ray_packet_t), pointer :: spacket, sp_prev
  class(single_ray_t),   pointer :: sray
! ......................................................................

  ! --- Sort out direction issues ---
  allocate (dir)                                                        ! create new ray direction
  if (associated(self%ray_list%head)) then                              ! if ray list is not empty
    dir%prev => self%ray_list%tail
    dir%prev%next => dir
    self%ray_list%tail => dir
    dir%head => self%ray_list%head
  else
    self%ray_list => dir
    self%ray_list%head => dir                                           ! if empty, point to itself
    self%ray_list%tail => dir
  end if

  dir%ID = self%n_ray_dirs
  dir%n_bundles = 0
  dir%directionality = self%dr

  ! --- Go through single bundles ---
  sbundle => sbundle%head
  do while (associated(sbundle))                                        ! loop over ray bundles to be saved
    dir%n_bundles = dir%n_bundles + 1
    allocate(bundle)
    bundle%dir_ID = dir%ID
    bundle%axis  = self%axis
    bundle%theta = self%theta
    bundle%phi   = self%phi
    if (.not.associated(dir%bundle)) then
      dir%bundle => bundle
      bundle%head => bundle
      bundle%tail => bundle
      bundle%ID = 1
    else
      bundle%head => dir%bundle
      bundle%prev => dir%bundle%tail
      bundle%prev%next => bundle
      dir%bundle%tail => bundle
      bundle%ID = bundle%prev%ID
    end if

    bundle%lwall = sbundle%lwall
    bundle%uwall = sbundle%uwall
    bundle%n_packets = 0

    ! --- Go through single packets ---
    spacket => sbundle%packets
    do while (associated(spacket))
      bundle%n_packets = bundle%n_packets + 1
      allocate(ray)
      if (.not.associated(bundle%rays)) then
        bundle%rays => ray
        ray%head => ray
        ray%tail => ray
        ray%ID = 1
      else
        ray%head => bundle%rays
        ray%prev => bundle%rays%tail
        ray%prev%next => ray
        bundle%rays%tail => ray
        ray%ID = ray%prev%ID
      end if

      ray%n_rays = spacket%n_rays
      ray%n_points = spacket%n_points
      call ray%fill_ray_list_data(spacket%rays)
      sp_prev => spacket
      spacket => spacket%next
      deallocate(sp_prev)
    end do
    sb_prev => sbundle
    sbundle => sbundle%next
    deallocate(sb_prev)
  end do

END SUBROUTINE finish_single_ray_list

!> =====================================================================
!> Sort through the ray geometry list and select appropriate ray bundles
!> for particular direction
!> =====================================================================
SUBROUTINE get_my_bundle (self, bundle, d_ID, b_ID)
  class(ray_geometry_t) :: self
  class(ray_bndl_t), pointer :: bundle
  class(ray_dirs_t), pointer :: dir
  integer :: b_ID, d_ID
  logical :: df, bf
!-----------------------------------------------------------------------
  df = .false.
  dir => self%ray_list%head
  do while (associated(dir).and.(.not.df))
    if (dir%id == d_ID) then
      df = .true.
    else
      dir => dir%next
    end if
  end do

  if (df) then
    bf = .false.
    bundle => dir%bundle
    do while(associated(bundle).and.(.not.bf))
      if (bundle%id == b_ID) then
        bf = .true.
      else
        bundle => bundle%next
      end if
    end do
    if (.not.bf) error stop "Appropriate bundle not found"
  else
    error stop "Appropriate ray direction not found"
  end if

END SUBROUTINE get_my_bundle


SUBROUTINE test_ray_geometry(self, dir, bundle_ID, ray)
  class(ray_geometry_t) :: self
  class(ray_bndl_t), pointer :: bundle
  integer :: dir, bundle_ID, ray, pt, p(3)

  call self%get_my_bundle(bundle, dir, bundle_ID)

  print *, 'Testing ray_geom for dir_ID, bundle_ID, ray:',bundle%dir_ID, bundle_ID, ray
  print *, 'packets in bundle:',bundle%n_packets
  print *, 'lwall, uwall', bundle%lwall, bundle%uwall
  print *, 'n_points in test ray:',bundle%rays%n_points

END SUBROUTINE test_ray_geometry

!> =============================================================================
!> Print summary info for ray_t
!> =============================================================================
SUBROUTINE print_ray (self, i)
  class(ray_t) :: self
  integer :: i
  class(ray_t), pointer :: ray
  !.............................................................................
  if (io%verbose < 3) return
  write (*,*) '      RAY_PACKET:', i
  write (*,*) '      number of rays:', self%n_rays
  write (*,*) '      number of points:', self%n_points
END SUBROUTINE print_ray

!> =============================================================================
!> Print summary info for ray_bndl_t
!> =============================================================================
SUBROUTINE print_bundle (self, i)
  class(ray_bndl_t) :: self
  integer :: i, j
  class(ray_t), pointer :: ray
  !.............................................................................
  if (io%verbose < 3) return
  write (*,*) '    RAY_BUNDLE:', i
  write (*,*) '    number of packets:', self%n_packets
  ray => self%rays
  j = 1
  do while (associated(ray))
    call ray%print (j)
    j = j+1
    ray => ray%next
  end do
END SUBROUTINE print_bundle

!> =============================================================================
!> Print summary info for ray_dirs_t
!> =============================================================================
SUBROUTINE print_dirs (self, i)
  class(ray_dirs_t) :: self
  integer :: i, j
  class(ray_bndl_t), pointer :: ray_bndl
  !.............................................................................
  if (io%verbose < 2) return
  write (*,*) '  RAY_DIR:', i
  write (*,*) '  number of bundles:', self%n_bundles
  ray_bndl => self%bundle
  j = 1
  do while (associated(ray_bndl))
    call ray_bndl%print (j)
    j = j+1
    ray_bndl => ray_bndl%next
  end do
END SUBROUTINE print_dirs

!> =============================================================================
!> Print summary info for ray_geometry_t
!> =============================================================================
SUBROUTINE print_geometry (self)
  class(ray_geometry_t) :: self
  integer :: j
  class(ray_dirs_t), pointer :: ray_dir
  !.............................................................................
  if (io%verbose < 1) return
  write (*,*) 'RAY_GEOMETRY:'
  write (*,*) 'number of directions:', self%n_ray_dirs
  ray_dir => self%ray_list
  j = 1
  do while (associated(ray_dir))
    call ray_dir%print (j)
    j = j+1
    ray_dir => ray_dir%next
  end do
END SUBROUTINE print_geometry

END MODULE ray_geometry_v2_mod
