!=======================================================================
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

!> packet_t: is a ray packet, containing n_rays of rays of same length (n_points).
!> it has all the necessary info about the ray - where it originates, where it
!> terminates (fwd_coords/rev_coords, depending on direction), which cells
!> it crosses (cell_coords), distance it traves through particular cell(ds),
!> and which cell face it crosses when it enters and leaves the cell (pt_wall)
!=======================================================================
MODULE ray_geometry_v3_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE math_mod
  implicit none
  private
! ----------------------------------------------------------------------
  type, public :: packet_t
    integer :: ID                                                       ! packet ID
    integer :: n_rays                                                   ! number of rays in this packet
    integer :: n_points                                                 ! ray lengths in this packet
    real    :: dstep(3)                                                 ! traverse index increment
    integer, dimension(:)    , pointer   :: p_id_fw => null()           ! packet ID for ray connection in forward direction
    integer, dimension(:)    , pointer   :: p_id_rev => null()          ! packet ID for ray connection in reverse direction
    integer, dimension(:)    , pointer   :: r_id_fw => null()           ! ray    ID for ray connection in forward direction
    integer, dimension(:)    , pointer   :: r_id_rev => null()          ! ray    ID for ray connection in reverse direction
    integer, dimension(:)    , pointer   :: b_id => null()              ! ray    ID  in the bundle
    integer, dimension(:,:,:), pointer   :: coords => null()            ! x|y|z,        ray point, ray #
    integer, dimension(:,:  ), pointer   :: fwd_coords  => null()       ! x|y|z,                   ray #, next point after the ray left the patch
    integer, dimension(:,:  ), pointer   :: rev_coords  => null()       ! x|y|z,                   ray #, point before the ray enters the patch
  contains
    procedure :: fill_ray_list_data
    procedure :: print => print_ray
  end type
! ----------------------------------------------------------------------
  type, public :: ray_bndl_t
    class(packet_t), dimension(:), pointer :: packet => null()             ! packet list in array mode
    integer, dimension(:), pointer   :: rays_fw => null()               ! ray locations in the next bundle
    integer, dimension(:), pointer   :: rays_rev => null()              ! ray locations in the previous bundle
    integer :: ID                                                       ! bundle ID
    integer :: dir_ID                                                   ! direction ID
    integer :: b_id_fw                                                  ! bundle ID for ray connection in forward direction
    integer :: b_id_rev                                                 ! bundle ID for ray connection in forward direction
    integer :: n_packets                                                ! number of packets in this bundle
    integer :: n_rays                                                   ! total number of rays in this bundle
    integer :: lwall(3)
    integer :: uwall(3)
    real    :: dstep(3)                                                 ! traverse index increment
    integer :: axis                                                     ! identificator along which axis ray goes 1-3
    real    :: theta                                                    ! angle between ray and axis [0, 45], stored in radians
    real    :: phi                                                      ! angle of the ray in the axis plane [0, 360), stored in radians
    real    :: mu
    real    :: weight
  contains
    procedure :: print => print_bundle
  end type
! ----------------------------------------------------------------------
  type, public :: ray_dirs_t
    class(ray_dirs_t), pointer :: head => null()
    class(ray_dirs_t), pointer :: next => null()
    class(ray_dirs_t), pointer :: prev => null()
    class(ray_dirs_t), pointer :: tail => null()
    class(ray_bndl_t), dimension(:), pointer :: bundle => null()
    integer :: ID
    integer :: n_bundles = 0                                            ! number of bundles in this direction
    integer :: directionality(3)                                        ! 0 - two directions, 1 - only forward, -1 - only reverse (x|y|z)
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
    integer :: li(3), ui(3), n(3), ng(3)
    real(kind=8) :: cell_dims(3)
    integer    :: geometry = 1 !cartesian = 1, spherical = 2, cylindrical = 3
    integer    :: wall_n_dirs(3) = 0                                    ! number of directions entering/leaving a particular wall
    logical    :: init_done = .false.
    logical    :: space_diagonals = .true.
    logical    :: do_chains = .true.
    integer    :: n_warm_up = 10
    integer    :: axis, dr(3)                                           ! convenience
    integer    :: total_n_bundles = 0
    integer(8) :: total_n_packets = 0
    real       :: theta, phi                                            ! convenience
    logical    :: no_mans_land
  contains
    procedure :: add_axis
    procedure :: add_diag
    procedure :: add_ray
    procedure :: finish_single_ray_list
    procedure :: finish_direction
    procedure :: get_my_geometry
    procedure :: init
    procedure :: init_cartesian
    procedure :: init_spherical
    procedure :: init_cylindrical
    procedure :: traverse
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
    integer                          :: ID
    integer, dimension(:,:), pointer :: coords => null()  ! x|y|z, ray point
    integer, dimension(3)            :: fwd_coords             ! x|y|z,
    integer, dimension(3)            :: rev_coords             ! x|y|z,
    real :: ds
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
    integer :: n_rays, n_points, packet_ID
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
    integer :: dstep(3)
    integer :: l_chain = 0
    integer :: bundle_ID
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
    integer                     :: uwall(3), uw(3)
    integer                     :: lwall(3), lw(3)
    real                        :: ds
  contains
    procedure :: add_ray_point
  end type
! ----------------------------------------------------------------------
  type :: traverse_t
    real :: ds                                                          ! distance a ray travels between points
    integer :: dstep(3)                                                 ! distance the ray covers by moving one cell size
    integer :: cell_max(3), cell_min(3)
    integer :: p(3), upstream(3), downstream(3)
    integer :: lwall(3)
    integer :: pos(3)
    logical :: oob, reset(2)
  contains
    procedure :: init           => init_traverse
    procedure :: posit
    procedure :: restart
  end type
! ----------------------------------------------------------------------
  type(ray_geometry_t), public :: ray_geometry
  logical :: first_time = .true.
  real(8), public, parameter:: pi = asin(1.0d0) * 2.0d0, radian = 180./pi
CONTAINS

!=======================================================================
!> Initialisation subroutine for ray geometry
!=======================================================================
SUBROUTINE init (self, geometry, n, ng, cell_dims, no_mans_land)
  class(ray_geometry_t) :: self
  class(ray_geometry_t), pointer :: new
  real(kind=8), optional :: cell_dims(3)
  integer :: n(3), ng(3)
  integer :: geometry
  logical :: no_mans_land
  real(kind=8) :: maxdim
! ......................................................................
  call trace_begin('ray_geometry_t%init')
  if (first_time) then
    first_time = .false.
    allocate (new)                                                        ! allocate new geometry object
    new%geometry = geometry
    new%n  = n
    new%li = ng+1                                                         ! cell domain lower point
    new%ui = new%li+n - 1                                                 ! cell domain upper point
    new%no_mans_land = no_mans_land
    if(present(cell_dims)) then
       if ((cell_dims(1) == cell_dims(2)).and.(cell_dims(2) == cell_dims(3))) then
         new%cell_dims = 1.
       else
         maxdim = maxval(cell_dims)
         new%cell_dims = cell_dims/maxdim
       end if
    else
        new%cell_dims(:) = 1.                                              ! default - cubic cells
    end if

    if(.not.associated(self%head)) then                                   ! connect to geometry list
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
  end if
  call trace_end
END SUBROUTINE init

SUBROUTINE get_my_geometry(self, geometry, n, ng, cell_dims, my_geom)
  class(ray_geometry_t) :: self
  class(ray_geometry_t), pointer :: this, my_geom
  real(kind=8), optional :: cell_dims(3)
  integer :: n(3), ng(3)
  integer :: geometry
  logical :: found
! ......................................................................
  !call trace_begin('ray_geometry_t%get_my_geometry')

  if(.not.associated(self%head)) then
    call self%init(geometry, n, ng,  cell_dims, self%no_mans_land)
    this => self%head
    found = .false.
    do while (associated(this).and..not.found)
      if ((this%geometry == geometry).and.all(this%n == n).and.&
           all(this%ng==ng)) then
        my_geom => this
        found = .true.
      else
        this => this%next
      end if
    end do
  end if

  nullify(this)
  !call trace_end
END SUBROUTINE get_my_geometry

!=======================================================================
!> Initialisation subroutine for cartesian geometry.
!> It is important to note, that although rays are created independently
!> for arbitrary directions (or can be not created, if user doesn't need
!> them), 45 degree rays inevitably cover all tree directions, so they
!> should be treated carefully if one doesn't need some directions.
!> Namelist parameters:
!>   nx, ny, nz     : number of angles relative to x-, y-, z-axis
!>   xti, yti, zti  : degrees relative to x-, y-, and z-axis
!>   dx, dy, dz     : ray directionality
!=======================================================================
SUBROUTINE init_cartesian(self)
  class(ray_geometry_t) :: self
  integer                         :: i, d45
  real, dimension(5)              :: xti=0.0, yti=0.0, zti=0.0
  integer, save                   :: nx=1, ny=0, nz=0
  integer, save                   :: dx=0, dy=0, dz=0
  logical, save                   :: space_diagonals=.true.
  logical, save                   :: do_chains = .true.
  integer, save                   :: n_warm_up=10
  logical, save                   :: first_time=.true.
  ! ----------------------------------------------------------------------------
  class(ray_bndl_t), pointer :: bnd
  class(ray_dirs_t), pointer :: dir
  integer:: ibn
  real:: wgt
! ......................................................................
  namelist /cartesian_ray_geometry_params/ do_chains, &
    nx, ny, nz, xti, yti, zti, dx, dy, dz, space_diagonals, n_warm_up
! ......................................................................
  call trace_begin('ray_geometry_t%init_cartesian',2)
  !> TO DO: add an option for 2D simulation setup (i.e. any(cell_dims) == 0)

  self%init_done = .true.
  allocate(self%ray_list)
  self%n_ray_dirs = 0

  call io%lock%set ('ray_geometry_t%init_cartesian')
    first_time = .false.
    rewind (io%input)
    read (io%input, cartesian_ray_geometry_params)
    if (io%master) write (io%output, cartesian_ray_geometry_params)
  call io%lock%unset ('ray_geometry_t%init_cartesian')
  self%dr(1) = dx
  self%dr(2) = dy
  self%dr(3) = dz
  self%space_diagonals = space_diagonals
  self%do_chains = do_chains                                            ! do connection between packets and separate rays, leaving to different walls
  self%n_warm_up = n_warm_up                                            ! number of run-ins to do rt in no-chains mode, before hydro updates are started
  d45 = 0                                                               ! flag for 45 degree angles
  if (nx > 0) then
    do i = 1, nx
      if (xti(i) == 0) then                                             ! 0 degree angle == axisymmetric
        call self%add_axis(1)
      else if (xti(i) == 45.) then                                      ! diagonal
        if (d45 < 1) then
          call self%add_diag(1,1)                                       ! first issue of diagonals
          d45 = d45 + 1
        else if (d45 < 2) then
          call self%add_diag(2,1)                                       ! second issue of diagonals
          d45 = d45 + 1
        end if
      else
        call self%add_ray(xti(i),1)                                     ! any other angle
      end if
    end do
  end if

  if (ny > 0) then
    do i = 1, ny
      if (yti(i) == 0) then
        call self%add_axis(2)                                           ! axisymmetric
      else if (yti(i) == 45.) then
       if (d45 < 1) then
          call self%add_diag(1,2)                                       ! first issue of diagonals
          d45 = d45 + 1
        else if (d45 < 2) then
          call self%add_diag(2,2)                                       ! second issue of diagonals
          d45 = d45 + 1
        end if
      else
        call self%add_ray(yti(i),2)
      end if
    end do
  end if

  if (nz > 0) then
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
  end if
  
  if (self%ray_list%directionality(1) == 0) then
    wgt = math%pi4/(self%n_ray_dirs * 2)
  else
    wgt = math%pi4/(self%n_ray_dirs * 2)
  end if
  dir => self%ray_list%head
  do while (associated(dir))
    do ibn = 1, dir%n_bundles
      dir%bundle%weight = wgt
    end do
    dir => dir%next
  end do
  call trace_end
END SUBROUTINE init_cartesian

!=======================================================================
!> Initialisation subroutine for cylindrical geometry
!=======================================================================
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

!=======================================================================
!> Initialisation subroutine for spherical geometry
!=======================================================================
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

!=======================================================================
!> Subroutine to add axisymmetric rays. It is done separately, because
!> it is much simpler and therefore faster than with rays with other
!> angles. Adds stiff ray arrays instantly.
!=======================================================================
SUBROUTINE add_axis(self, dir)
  class(ray_geometry_t) :: self
  class(ray_dirs_t), pointer :: ndir
  integer :: dir, ia, ib, ic, ray_n, a, b, c, rp, p(3), pmax(3), dstep(3), nr
  real    :: ds
  class(ray_bndl_t), pointer :: bundle
  class(packet_t)  , pointer :: packet
  integer :: ibn, jbn, jpck, ray1, ray2
! ......................................................................
 ! call trace_begin('ray_geometry_t%add_axis', 4)
  pmax = (self%ui - self%li) + 1                                        ! Total number of points
  self%n_ray_dirs = self%n_ray_dirs + 1                                 ! Adding new direction to the ray geometry
  dstep = 0
  SELECT case(dir)                                                      ! select the direction
  case(1)
    a = 2; b = 3; c = 1                                                 ! along x axis
  case(2)
    a = 3; b = 1; c = 2                                                 ! along y axis
  case(3)
    a = 1; b = 2; c = 3                                                 ! along z axis
  end SELECT
  dstep(c) = 1
  ! Create new ray bundle and add it to ray_geometry
  allocate (ndir)                                                       ! create new direction
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
  nr = pmax(a)*pmax(b)
  ! --- Fill in direction info ---
  allocate(ndir%bundle(1))
  bundle => ndir%bundle(1)
  ndir%n_bundles = 1                                                    ! there can only be one bundle for this direction
  self%total_n_bundles = self%total_n_bundles + 1                       ! increment total number of bundles by 1
  self%total_n_packets = self%total_n_packets + 1
  ndir%directionality = self%dr                                         ! one/two directional
  ndir%ID = self%n_ray_dirs                                             ! get the dir_id
  self%wall_n_dirs(c) = self%wall_n_dirs(c)+1
  ! --- Fill in bundle info ---
  allocate(bundle%packet(1))                                            ! packet for this bundle
  bundle%axis      = c                                                  ! direction axis of this bundle
  bundle%theta     = 0
  bundle%mu        = cos(bundle%theta)
  bundle%phi       = 0
  bundle%id        = self%total_n_bundles                               ! bundle ID in this direction
  bundle%dir_id    = ndir%id                                            ! dir_ID of this bundle
  bundle%n_packets = 1                                                  ! there can only be one packet for this bundle
  bundle%n_rays    = nr                                                 ! number of rays in this direction
  bundle%dstep     = dstep                                              ! step along direction axis
  bundle%lwall(a)  = 0                                                  ! lower wall along A axis
  bundle%uwall(a)  = 0                                                  ! upper wall along A axis
  bundle%lwall(b)  = 0                                                  ! lower wall along B axis
  bundle%uwall(b)  = 0                                                  ! upper wall along B axis
  bundle%lwall(c)  =-1                                                  ! lower wall for this direction
  bundle%uwall(c)  = 1                                                  ! upper wall for this direction
  ! --- Fill in rays info ---
  packet => bundle%packet(1)
  packet%id = self%total_n_packets                                      ! it's packet ID
  packet%n_rays   = nr                                                  ! number of rays in this direction
  packet%dstep    = dstep                                               ! step along direction axis
  packet%n_points = pmax(c)                                             ! number of points per ray corresponds to n_cells in that direction


  allocate (packet%coords(3,pmax(c),nr))                                ! create ray table, containing x|y|z loc, for all n_cells and n_rays
  allocate (packet%fwd_coords(3,nr))                                    ! create table, containing x|y|z coords for upstream position for particular ray, for all rays
  allocate (packet%rev_coords(3,nr))                                    ! create table, containing x|y|z coords for downstream position for particular ray, for all rays

  bundle%b_id_fw  = self%total_n_bundles                                ! bundle ID for ray connection in forward direction
  bundle%b_id_rev = self%total_n_bundles                                ! same for bundle ID for ray connection in forward direction
  allocate(packet%p_id_fw(nr))                                          ! packet  ID for ray connection in forward direction
  packet%p_id_rev =>   packet%p_id_fw(:)                                ! same for packet ID for ray connection in reverse direction
  allocate(packet%r_id_fw(nr))                                          ! ray ID for ray connection in forward direction
  packet%r_id_fw(:) = 1
  packet%p_id_fw(:) = self%total_n_packets
  packet%r_id_rev => packet%r_id_fw(:)                                  ! same for ray ID for ray connection in reverse direction
  bundle%rays_fw  => packet%r_id_fw(:)
  bundle%rays_rev => packet%r_id_fw(:)
  packet%b_id     => bundle%rays_rev(:)
  ray_n = 1                                                             ! ray number
  do ia = self%li(a), self%ui(a)                                        ! iterate over cells in first axis
    do ib = self%li(b), self%ui(b)                                      ! iterate over cells in second axis
      p(a) = ia; p(b) = ib                                              ! point position in first two axes
      packet%rev_coords(a,ray_n) = p(a)                                 ! upstream position in terms of first axis
      packet%rev_coords(b,ray_n) = p(b)                                 ! upstream position in terms of second axis
      packet%rev_coords(c,ray_n) = self%li(c) - 1                       ! upstream position in terms of ray direcion
      rp = 1
      do ic = self%li(c), self%ui(c)                                    ! move along ray direction
        p(c) = ic                                                       ! current cell position in terms of third axis
        packet%coords(:,rp,ray_n) = p                                   ! cell position in particular ray
        rp = rp+1                                                       ! increment ray point
      end do
      packet%fwd_coords(a,ray_n) = p(a)                                 ! downstream position in terms of first axis
      packet%fwd_coords(b,ray_n) = p(b)                                 ! downstream position in terms of second axis
      packet%fwd_coords(c,ray_n) = self%ui(c) + 1                       ! downstream position in terms of ray direction
      packet%r_id_fw(ray_n) = ray_n
      ray_n = ray_n + 1                                                 ! move to next ray
    end do
  end do
!  call trace_end
END SUBROUTINE add_axis

!=======================================================================
!> Subroutine to add rays (except 0 and 45 degree angles).
!> Done using modified "fast voxel traversal" algorithm
!> (original, not modified, one is from Amanatides & Woo, 1987)
!=======================================================================
SUBROUTINE add_ray(self, theta, dir)
  class(ray_geometry_t) :: self
  integer :: dir, nphi, a, b, c, pdir
  real :: dp(3), init_step(3), theta, t, phi
  real, dimension(:), allocatable :: phis
  logical :: oob_dir
  ! ....................................................................
 ! call trace_begin('ray_geometry_t%add_ray', 4)

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

    dp(a) = sin(t)*cos(phi)                                             ! step size multiplier along a direction
    dp(b) = sin(t)*sin(phi)                                             ! step size ,ultiplier along b direction
    dp(c) = self%cell_dims(c) / cos(t)                                  ! step size along c direction


                                                                        ! removing computer uncertainties
    dp(a) = round(dp(a))                                                ! and giving angle error of around 1 degree
    dp(b) = round(dp(b))                                                ! this is necessary in order to avoid infinite loops
    dp(c) = round(dp(c))                                                ! and overabundance of ray points with machine epsilon ds

    init_step(1) = 0                                                    ! cell center for a direction
    init_step(2) = 0                                                    ! cell center for b direction
    init_step(3) = 0.5                                                  ! cell wall for c direction (ray propagation direction)
    call self%traverse (dp, a, b, c, init_step)                         ! proceed with 3D traversal
  end do

  deallocate(phis)
  !call trace_end
END SUBROUTINE add_ray

!=======================================================================
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
!=======================================================================

SUBROUTINE add_diag(self, rep, dir)
  class(ray_geometry_t) :: self
  real :: init_step(3), dp(3), theta, t, phi
  real, dimension(8) :: phis
  integer :: rep, dir, nphi, a, b, c, pdir
! ......................................................................
  !call trace_begin('ray_geometry_t%add_diag',4)

  t = 45./radian

  SELECT case(dir)
  case(1)
    a = 1; b = 2; c = 3                                                 ! along z axis
  case(3)
    a = 2; b = 3; c = 1                                                 ! along x axis
  case(2)
    a = 3; b = 1; c = 2                                                 ! along y axis
  end SELECT

  SELECT case(rep)                                                      ! 45 degree rays can be added only twice
  case(1)
    if (self%space_diagonals) then
      nphi = 8
      phis(1:8) = (/0., 45., 90., 135., 180., 225., 270., 315./)/radian
    else
      nphi = 4
      phis(1:4) = (/0., 90., 180., 270./)/radian
    end if
  case(2)
    nphi = 2
    phis(1:2) = (/90., 270./)/radian
  end SELECT

  if (rep < 3) then
    do pdir = 1, nphi                                                     ! loop over phis
      self%n_ray_dirs = self%n_ray_dirs + 1                               ! increment ray dirs
      self%axis  = dir
      self%theta = t
      self%phi   = phis(pdir)
      phi=phis(pdir)
      dp(a) = sin(t)*cos(phi)                                           ! step size multiplier along a direction
      dp(b) = sin(t)*sin(phi)                                           ! step size ,ultiplier along b direction
      dp(c) = self%cell_dims(c) / cos(t)                                ! step size along c direction

                                                                        ! removing computer uncertainties
      dp(a) = round(dp(a))                                                ! and giving angle error of around 1 degree
      dp(b) = round(dp(b))                                                ! this is necessary in order to avoid infinite loops
      dp(c) = round(dp(c))                                                ! and overabundance of ray points with machine epsilon ds

      if (dp(a)==0) then                                                  ! if a direction == 0
        init_step(1) = 0                                                  ! cell center for a direction
        init_step(2) = 0.5                                                ! cell wall for b direction
        init_step(3) = 0.5                                                ! cell wall for c direction (ray propagation direction)

      else if (dp(b) == 0) then                                           ! if b direction == 0
        init_step(1) = 0.5                                                ! cell wall for a direction
        init_step(2) = 0                                                  ! cell center for b direction
        init_step(3) = 0.5                                                ! cell wall for c direction (ray propagation direction)

      else                                                                ! movement to all directions
        init_step(1) = 0.5                                                ! cell wall for a direction
        init_step(2) = 0.5                                                ! cell wall for b direction
        init_step(3) = 0.5                                                ! cell wall for c direction (ray propagation direction)
      end if

      call self%traverse (dp, a, b, c, init_step)                         ! proceed with 3D traversal

    end do
  end if
  !call trace_end
END SUBROUTINE add_diag

!=======================================================================
!> Subroutine to traverse rays in 3D. Similar to "fast voxel traversal
!> algorithm" (Amanatides & Woo, 1987) and appendix A in Rijkhorst et al. (2006),
!> but without their serious drawbacks and is faster.
!=======================================================================
SUBROUTINE traverse(self, dp, a, b, c, init_step)
  class(ray_geometry_t) :: self
  class(ray_point_t), pointer :: ray
  class(single_ray_bundle_t), pointer :: ray_bundle
  class(traverse_t), pointer :: t
  integer :: a, b, c, ia, ib, x, r_id
  real :: init_step(3), dp(3)
  logical :: first_point
! ......................................................................
  !call trace_begin('traverse_t%traverse')
  r_id = 0
  allocate(t)                                                           ! allocate traversal object
  allocate(ray_bundle)                                                  ! create ray tree to store rays
  call t%init(dp, self%ui, self%li, a, b, c)                            ! initialize traversal object
  t%dstep(c) = self%cell_dims(c)
  ray_bundle%dstep = t%dstep                                            ! save for fast access later
  do ia = t%cell_min(1), t%cell_max(1)                                  ! loop over 1st direction
    do ib = t%cell_min(2), t%cell_max(2)                                ! loop over 2nd direction
      first_point = .true.
      allocate(ray)                                                     ! create new ray
      r_id = r_id + 1                                                   ! increment ray ID
      call t%posit(ia, ib, init_step)                                   ! position traversal object
      do while(.not.t%oob)                                              ! loop until ray reaches wall at ray direction
        call ray%add_ray_point(t%pos)                                   ! add a ray point to the ray
        if (first_point) then
          ray%head%lwall = t%lwall                                      ! ray enters the patch through this wall
          ray%head%ds = t%ds
          first_point = .false.
        end if

        t%pos = t%pos + t%dstep                                         ! move to next point

        if (any((t%pos > t%cell_max).or. &
                (t%pos < t%cell_min))) then                             ! if the ray is fully traced
          ray%head%fwd = t%pos
          ray%head%rev = t%upstream                                     ! save ray upstream cell

          do x = 1,3
            if (t%dstep(t%p(x)) < 0 ) then
              if (t%pos(t%p(x)) < t%cell_min(x)) ray%head%uwall(t%p(x)) = -1
            else
              if (t%pos(t%p(x)) > t%cell_max(x)) ray%head%uwall(t%p(x)) = 1
            end if
          end do
          call ray_bundle%finish_ray(ray, r_id)                         ! and save ray to the ray tree
          if (t%pos(t%p(3)) > t%cell_max(3)+0.5) then                   ! if cell wall along the ray direction is reached
              t%oob = .true.                                            ! this ray is done

          else                                                          ! otherwise
            allocate(ray)                                               ! get new ray
            r_id = r_id + 1                                             ! increment ray ID
            call t%restart                                              ! and restart the 1st or 2nd direction
            first_point = .true.
          end if
        end if

      end do
    end do
  end do
  call self%finish_single_ray_list(ray_bundle)                          ! save the ray tree to ray_geometry
  deallocate(t)
  !call trace_end
END SUBROUTINE traverse

!=======================================================================
!> INPUT: self, step size along ray direction (dp), smallest cell coord (li)
!> highest cell coord (ui), nominal direction (a, b), ray direction (c)
!=======================================================================
SUBROUTINE init_traverse(self, dp, ui, li, a, b, c)
  class(traverse_t) :: self
  integer :: a, b, c, li(3), ui(3)
  real :: dp(3)
! ......................................................................
  !call trace_begin('traverse_t%init')
  self%p(1) = a; self%p(2) = b; self%p(3) = c                           ! convenience
  self%cell_max(1) = ui(a)                                              ! last cell point along a axis
  self%cell_max(2) = ui(b)                                              ! last cell point along b axis
  self%cell_max(3) = ui(c)                                              ! last cell point along c axis
  self%cell_min(1) = li(a)                                              ! first cell point along a axis
  self%cell_min(2) = li(b)                                              ! first cell point along b axis
  self%cell_min(3) = li(c)                                              ! first cell point along c axis
  self%ds = dp(c)
  self%dstep = 0
  if (dp(a) /= 0) then
    self%dstep(a) = nint(dp(a)*self%ds)                                 ! distance between points along a axis
  end if
  if (dp(b) /= 0) then
    self%dstep(b) = nint(dp(b)*self%ds)                                 ! distance between points along b axis
  end if
  self%dstep(c) = 1                                                     ! distance between points along c axis
  !call trace_end
END SUBROUTINE init_traverse

!=======================================================================
!> INPUT: self, position in nominal directions (ia, ib),
!> distance to cell walls (init_step), upstream cell position for particular
!> ray (upstream)
!=======================================================================
SUBROUTINE posit(self, ia, ib, init_step)
  class (traverse_t) :: self
  real, dimension(3) :: init_step
  integer :: ia, ib, x
! ......................................................................
  self%pos(self%p(1)) = ia                                              ! position in a axis
  self%pos(self%p(2)) = ib                                              ! position in b axis
  self%pos(self%p(3)) = self%cell_min(3)                                ! starting position for c axis
  self%lwall = 0


  do x = 1,3
    self%upstream(self%p(x)) = self%pos(self%p(x)) - self%dstep(self%p(x))

    if (self%dstep(self%p(x)) > 0 ) then
      if (self%upstream(self%p(x)) < self%cell_min(x)) self%lwall(self%p(x)) = -1
    else
      if (self%upstream(self%p(x)) > self%cell_max(x)) self%lwall(self%p(x)) = 1
    end if
  end do

  self%oob     = .false.
  self%reset(1) = .false.
  self%reset(2) = .false.
END SUBROUTINE posit

!=======================================================================
!> restart the ray when it hits a side-wall (invoke pseudo-periodic
!> boundaries). This is for the general, 3D case.
!=======================================================================
SUBROUTINE restart (self)
  class(traverse_t) :: self
  integer :: x
! ......................................................................
  self%lwall = 0
  do x = 1,2
    if (self%dstep(self%p(x)) /= 0) then
      if (self%dstep(self%p(x)) > 0) then
        if (self%pos(self%p(x)) > self%cell_max(x)) self%reset(x) = .true.
      else
        if (self%pos(self%p(x)) < self%cell_min(x)) self%reset(x) = .true.
      end if
    end if
  end do

  do x = 1,2
    if(self%reset(x)) then
      if (self%dstep(self%p(x)) > 0) then
        self%pos(self%p(x)) = self%cell_min(x)
      else
        self%pos(self%p(x)) = self%cell_max(x)
      end if
    end if
  end do

  do x = 1,3
    self%upstream(self%p(x)) = self%pos(self%p(x)) - self%dstep(self%p(x))

    if (self%dstep(self%p(x)) > 0 ) then
      if (self%upstream(self%p(x)) < self%cell_min(x)) self%lwall(self%p(x)) = -1
    else
      if (self%upstream(self%p(x)) > self%cell_max(x)) self%lwall(self%p(x)) = 1
    end if
  end do

  self%reset(1) = .false.
  self%reset(2) = .false.
END SUBROUTINE restart

!=======================================================================
!> Subroutine to create linked-list of ray points. Later this list is
!> redone into a stiff (array) ray.
!> INPUT: self, cell coordinates where the ray point have moved to (p)
!=======================================================================
SUBROUTINE add_ray_point(self, p)
  class(ray_point_t)          :: self
  class(ray_point_t), pointer :: point
  integer :: p(3)
! ......................................................................

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

END SUBROUTINE add_ray_point

!=======================================================================
!> Subroutine to transform linked-list of ray points into a stiff (array)
!> ray. Arranging is done based on the walls the ray hits AND the number
!> of ray points.
!> INPUT: self, linked list of ray points (ray_point)
!=======================================================================

SUBROUTINE fill_single_ray_data(self, ray)
  class(single_ray_t) :: self
  class(ray_point_t), pointer :: ray, prev
  integer :: pt, li(3), ui(3), w
! ......................................................................
  prev => ray%tail
  ray  => ray%head

  self%lwall = ray%lwall
  self%uwall = ray%uwall

  self%n_points = ray%n_points                                          ! store number of points in this ray

  allocate(self%coords(3, ray%n_points))                                ! create a stiff array for this ray with cell cords
  self%fwd_coords = ray%fwd                                             ! save upstream cell coords
  self%rev_coords = ray%rev                                             ! save downstream cell coords
  pt=1
  do while (associated(ray))                                            ! loop over ray points
    self%coords(:,pt)  = ray%coords                                     ! save ray point coords
    prev => ray                                                         ! save previous point
    ray => ray%next                                                     ! go to next point
    deallocate(prev)                                                    ! deallocate previous point
    pt = pt+1
  end do
END SUBROUTINE fill_single_ray_data

!=======================================================================
!> finish a single ray - put it into a packet and a bundle. For the time
!> being, rays are not merged into arrays.
!=======================================================================
SUBROUTINE finish_ray(self, ray_point, id)
  class(single_ray_bundle_t) :: self
  class(single_ray_bundle_t), pointer :: bundle
  class(single_ray_packet_t), pointer :: packet
  class(single_ray_t), pointer  :: ray
  class(ray_point_t),  pointer  :: ray_point
  integer :: id
  logical :: fnd_p, fnd_b
! ......................................................................
  !call trace_begin('ray_geometry_t%finish_ray')

  allocate (ray)                                                        ! allocate ray list
  ray%id = id                                                           ! unique ray ID for this direction
  call ray%fill_single_ray_data(ray_point)                              ! store ray information of this ray

  ! ----- choose an appropriate bundle -----
  if (.not.associated(self%head)) then                                  ! if the ray tree is empty
    allocate(bundle)                                                    ! allocate a new bundle
    bundle%dstep = self%dstep                                           ! store basic data - dstep
    self%head => bundle                                                 ! point to self
    self%tail => bundle                                                 ! be such self-centered...
    bundle%lwall = ray%lwall                                            ! lwall of this bundle
    bundle%uwall = ray%uwall                                            ! uwall of this bundle
    bundle%bundle_ID = 1                                                ! this is the first bundle in the bundle list
  else
    bundle => self%head
    fnd_b = .false.
    do while (associated(bundle).and..not.fnd_b)                        ! look for a bundle, which match the criterion
      if (all(bundle%lwall == ray%lwall).and. &                         ! bundles are defined by which walls
          all(bundle%uwall == ray%uwall)) then                          ! the ray enters and leaves the patch
        fnd_b = .true.
      else
        bundle => bundle%next
      end if
    end do
    if(.not.fnd_b) then
      allocate(bundle)                                                  ! if appropriate bundle is not found, create a new one
      bundle%dstep = self%dstep                                         ! fill in basic data - dstep
      bundle%head => self%head                                          ! point to first bundle
      bundle%prev => self%tail                                          ! add itself to bundle list
      bundle%prev%next => bundle
      self%tail => bundle
      bundle%lwall = ray%lwall                                          ! lwall of this bundle
      bundle%uwall = ray%uwall                                          ! uwall of this bundle
      bundle%bundle_ID = bundle%prev%bundle_ID + 1                      ! save it's bundle ID
    end if
  end if
  ! ----- choose an appropriate packet -----
  if (.not.associated(bundle%packets)) then                             ! check if there are packets for selected bundle
    allocate(packet)                                                    ! allocate a new packet
    bundle%packets => packet
    packet%head => packet
    packet%tail => packet
    packet%n_rays = 1                                                   ! this will be the first ray in the packet
    packet%n_points = ray%n_points                                      ! packet will contain rays of such length
    packet%packet_ID = 1                                                ! new packet ID. Since it is the first packet, it is 1
  else                                                                  ! if packet list is not empty
    packet => bundle%packets                                            ! browse through packets
    fnd_p = .false.
    do while (associated(packet).and..not.fnd_p)
      if (ray%n_points == packet%n_points) then                         ! rays with same number of points are packed together
        fnd_p = .true.
        packet%n_rays = packet%n_rays + 1                               ! number of rays in this packet will increase
      else
        packet => packet%next
      end if
    end do
    if(.not.fnd_p) then                                                 ! if a packet with needed number of ray points is not found
      allocate(packet)                                                  ! create a new one
      packet%head => bundle%packets                                     ! add it to the list
      packet%prev => bundle%packets%tail
      packet%prev%next => packet
      bundle%packets%tail => packet
      packet%n_rays = 1                                                 ! this will be the first ray
      packet%n_points = ray%n_points                                    ! this packet will contain rays with such number of points
      packet%packet_ID = packet%prev%packet_ID + 1                      ! get an ID for this packet
    end if
  end if
  ! ----- add the ray to the packet -----
  if (.not.associated(packet%rays)) then                                ! if there are no rays in the packet
    packet%rays => ray
    ray%head => ray
    ray%tail => ray
  else                                                                  ! otherwise add it to the end
    ray%head => packet%rays
    ray%prev => packet%rays%tail
    ray%prev%next => ray
    packet%rays%tail => ray
  end if

  nullify (ray, bundle, packet, ray_point)
  !call trace_end
END SUBROUTINE finish_ray

SUBROUTINE fill_ray_list_data(self, sray, rid)
  class(packet_t) :: self
  class(single_ray_t), pointer :: sray, prev_sray
  integer :: iray,rid
! ......................................................................
  allocate (self%coords(3,self%n_points,self%n_rays))                   ! cell coords per point per ray in the bundle
  allocate (self%fwd_coords(3,self%n_rays))                             ! upstream cell coords per ray
  allocate (self%rev_coords(3,self%n_rays))                             ! downstream cell coords per ray
  allocate (self%b_id(self%n_rays))
  sray => sray%head
  iray = 0
  do while (associated(sray))                                           ! loop over all rays in the tree
    iray = iray + 1
    self%coords    (:,:,iray)    = sray%coords                          ! save cell cords
    self%fwd_coords(:,  iray)    = sray%fwd_coords                      ! save upstream coords
    self%rev_coords(:,  iray)    = sray%rev_coords                      ! save downstream coords
    self%b_id(iray) = rid
    rid = rid+1
    prev_sray   => sray                                                 ! hold pointer to this ray
    deallocate(sray%coords)
    sray => sray%next                                                   ! go to next ray
    deallocate (prev_sray)                                              ! deallocate previous ray
  end do
  nullify (sray, prev_sray)
END SUBROUTINE fill_ray_list_data

!=======================================================================
!> Subroutine to restructure ray tree into 2D sorted ray tables.
!=======================================================================
SUBROUTINE finish_single_ray_list(self, sbundle)
  class(ray_geometry_t) :: self
  class(ray_dirs_t), pointer :: dir
  class(ray_bndl_t), pointer :: bundle
  class(packet_t)  , pointer :: packet
  ! ---
  class(single_ray_bundle_t), pointer :: sbundle, sb_prev, tbundle
  class(single_ray_packet_t), pointer :: spacket, sp_prev, tpacket, sp_next
  class(single_ray_t),   pointer :: sray
  logical :: found
  integer :: r_id, b_ID, p_id
! ......................................................................
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

  dir%ID = self%n_ray_dirs                                              ! dir_ID for this direction
  dir%n_bundles = 0                                                     ! initially, there are no bundles
  dir%directionality = self%dr                                          ! set directionality

  ! =============== IF CHAINS ARE DONE ORGANIZE THIS WAY ===============
  if (self%do_chains) then

  ! --- Go through single ray bundles, count them ---
  sbundle => sbundle%head                                               ! take the first bundle
  tbundle => sbundle                                                    ! save the head for later
  b_ID = 0                                                              ! set the counter to 0
  do while (associated(sbundle))                                        ! go through bundles
    b_id = b_id+1
    sbundle => sbundle%next
  end do
  ! --- allocate the bundle array, corresponding to the number of bundles
  dir%n_bundles = b_id                                                  ! now we know how many bundles there are in this direction
  allocate(dir%bundle(b_id))                                            ! allocate the bundles

  ! --- Now go through single ray bundles and organize data
  sbundle => tbundle
  b_ID = 1
  do while (associated(sbundle))                                        ! loop over the bundles again
    self%total_n_bundles = self%total_n_bundles + 1                     ! increment the total number of bundles in ray geometry
    bundle => dir%bundle(b_id)
    bundle%id     = self%total_n_bundles                                ! bundle ID
    bundle%dir_ID = dir%ID                                              ! dir ID
    bundle%axis   = self%axis                                           ! axis of ray flow
    bundle%theta  = self%theta
    bundle%mu     = cos(bundle%theta)
    bundle%phi    = self%phi
    bundle%lwall  = sbundle%lwall                                       ! wall where bundle enters the box
    bundle%uwall  = sbundle%uwall                                       ! wall where bundle leaves the box
    bundle%dstep  = sbundle%dstep                                       ! step through the cells

    ! --- go through packets in this single ray bundle and count them
    p_id = 0
    spacket => sbundle%packets                                          ! take the first packet
    do while(associated(spacket))
      p_id=p_id+1
      spacket => spacket%next
    end do

    bundle%n_packets = p_id                                             ! set the number of packets in this bundle
    allocate(bundle%packet(p_id))                                       ! allocate the packets

    spacket => sbundle%packets

    p_id=1
    r_id = 1
    do while(associated(spacket))                                       ! go through packets again
      self%total_n_packets = self%total_n_packets + 1
      packet => bundle%packet(p_id)
      packet%id       = self%total_n_packets                            ! save the packet ID
      packet%n_rays   = spacket%n_rays                                  ! pick up the number of rays in this packet
      packet%n_points = spacket%n_points                                ! how long the rays are in this packet
      packet%dstep    = bundle%dstep
      bundle%n_rays   = bundle%n_rays + spacket%n_rays                  ! increment the number of rays in this bundle
      call packet%fill_ray_list_data(spacket%rays,r_id)                 !fill in the ray data


      p_id = p_id+1                                                     ! increment the packet ID
      sp_prev => spacket
      spacket => spacket%next
      deallocate(sp_prev)                                               ! deallocate the single ray packet
    end do
    b_ID = b_ID + 1
    sb_prev => sbundle
    sbundle => sbundle%next
    deallocate(sb_prev)                                                 ! deallocate the single ray bundle
  end do
  nullify(tbundle)
  call self%finish_direction(dir)                                       ! finalize the direction by connecting the ray chains

  ! ============= IF CHAINS ARE NOT DONE ORGANIZE THIS WAY =============
  else
    sbundle => sbundle%head
    dir%n_bundles=1                                                     ! there can only be one
    self%total_n_bundles = self%total_n_bundles + 1                     ! increment total number of bundles for this geometry
    allocate(dir%bundle(1))                                             ! allocate the bundle
    bundle => dir%bundle(1)
    bundle%dir_id = dir%id                                              ! give its dir ID
    bundle%axis   = self%axis                                           ! save axis
    bundle%theta  = self%theta                                          ! save theta
    bundle%mu     = cos(bundle%theta)
    bundle%phi    = self%phi                                            ! save phi
    bundle%id     = self%total_n_bundles                                ! it is the only bundle
    bundle%dstep  = sbundle%dstep

   ! -- loop over the bundles and merge all packets of same ray lenght
   ! by placing all of them into a temporary sbundle
    allocate (tbundle)
    do while (associated(sbundle))                                      ! loop over ray bundles
    ! --- Go through single packets ---
      spacket => sbundle%packets
      do while (associated(spacket))
        sp_next => spacket%next
        if (.not.associated(tbundle%packets)) then
          tbundle%packets => spacket                                    ! pop the first packet
          nullify(spacket%next, spacket%prev)
          spacket%tail => spacket
          spacket%head => spacket
        else
          found = .false.
          tpacket => tbundle%packets

          find_loop: do while (associated(tpacket))

            if (tpacket%n_points == spacket%n_points) then              ! merge
               tpacket%n_rays = tpacket%n_rays + spacket%n_rays
               spacket%rays%prev => tpacket%rays%tail
               tpacket%rays%tail%next => spacket%rays
               spacket%rays%head => tpacket%rays
               tpacket%rays%tail => spacket%rays%tail
               deallocate(spacket)
               found = .true.
               exit find_loop
            end if
            if (tpacket%n_points < spacket%n_points) then               ! prepend
              if (.not.associated(tpacket%prev)) then
                tbundle%packets => spacket
                nullify(spacket%prev)
              else
                spacket%prev => tpacket%prev
              end if
              tpacket%prev => spacket
              spacket%next => tpacket
              found = .true.
              exit find_loop
            end if
            tpacket => tpacket%next
          end do find_loop
          if (.not.found) then                                          ! append
            spacket%prev => tbundle%packets%tail
            tbundle%packets%tail%next=> spacket
            tbundle%packets%tail => spacket
            spacket%head => tbundle%packets
            nullify(spacket%next)
          end if
        end if
        spacket => sp_next
      end do
    sb_prev => sbundle
    sbundle => sbundle%next
    deallocate(sb_prev)
    end do

   spacket => tbundle%packets
   p_id = 0
   do while (associated(spacket))
     p_id = p_id + 1
     spacket => spacket%next
   end do

   bundle%n_packets = p_id
   allocate(bundle%packet(p_id))
    r_id = 1
    p_id = 1
    spacket => tbundle%packets
    do while(associated(spacket))
      self%total_n_packets = self%total_n_packets + 1
      sp_next => spacket%next
      packet => bundle%packet(p_id)
      packet%id        = self%total_n_packets                           ! unique ID for this packet
      packet%n_rays    = spacket%n_rays                                 ! number of rays in this packet
      packet%n_points  = spacket%n_points                               ! how long are the rays in this packet
      packet%dstep     = bundle%dstep
      bundle%n_rays = bundle%n_rays + spacket%n_rays                    ! number of rays in the bundle
      call packet%fill_ray_list_data(spacket%rays,r_id)                 !fill in the ray data
      p_id = p_id+1
      deallocate(spacket)
      spacket => sp_next
    end do
  end if
END SUBROUTINE finish_single_ray_list

!=======================================================================
!> Finish ray direction by sorting ray chains, or, if chains are not done,
!> by merging all bundles into one
!=======================================================================
SUBROUTINE finish_direction(self, dir)
  class(ray_geometry_t) :: self
  class(ray_dirs_t), pointer :: dir
  !.....................................................................
  class(ray_bndl_t), pointer :: bnd1, bnd2
  class(packet_t)  , pointer :: packet1, packet2
  integer :: ibn, jbn, ipck, jpck, ray1, ray2, coords(3), i
  integer :: neq_fw, neq_rv, diff_fw,diff_rv, d_bn
  logical :: wallupd(3)
  !---------------------------------------------------------------------
  wallupd(:) = .false.


  do ibn = 1, dir%n_bundles                                             ! loop over bundles
    bnd1 => dir%bundle(ibn)
    if (.not. associated(bnd1%rays_fw))  allocate(bnd1%rays_fw (bnd1%n_rays))
    if (.not. associated(bnd1%rays_rev)) allocate(bnd1%rays_rev(bnd1%n_rays))

    do ipck = 1, bnd1%n_packets                                         ! loop over packets
      packet1 => bnd1%packet(ipck)
      if (.not.associated(packet1%r_id_fw)) then                        ! if chain connection data is not allocated
        allocate(packet1%p_id_fw(packet1%n_rays))                       ! packet ids forward
        allocate(packet1%p_id_rev(packet1%n_rays))                      ! packet ids reverse
        allocate(packet1%r_id_fw(packet1%n_rays))                       ! ray ids forward
        allocate(packet1%r_id_rev(packet1%n_rays))                      ! ray ids reverse
        bnd1%b_id_fw    = 0; bnd1%b_id_rev    = 0
        packet1%p_id_fw = 0; packet1%p_id_rev = 0
        packet1%r_id_fw = 0; packet1%r_id_rev = 0
      end if

      do ray1=1, packet1%n_rays
        if (self%no_mans_land) then
          coords = packet1%fwd_coords(:,ray1)
        else
          coords = packet1%coords(:,packet1%n_points,ray1)
        end if
        do i=1,3
          if (bnd1%uwall(i) > 0) then
            coords(i) = self%li(i)
            if (.not.wallupd(i)) then                                   ! update the leaving wall, if needed
              self%wall_n_dirs(i) = self%wall_n_dirs(i)+1
              wallupd(i) = .true.
            end if
          else if (bnd1%uwall(i)<0) then
            coords(i) = self%ui(i)
            if (.not.wallupd(i)) then                                   ! update the leaving wall, if needed
              self%wall_n_dirs(i) = self%wall_n_dirs(i)+1
              wallupd(i) = .true.
            end if
          end if
        end do
        do jbn = 1, dir%n_bundles
          bnd2 => dir%bundle(jbn)
          if (.not. associated(bnd2%rays_fw))  allocate(bnd2%rays_fw (bnd2%n_rays))
          if (.not. associated(bnd2%rays_rev)) allocate(bnd2%rays_rev(bnd2%n_rays))

          do jpck = 1, bnd2%n_packets
            packet2 => bnd2%packet(jpck)

            do ray2=1, packet2%n_rays
              if (all(coords == packet2%coords(:,1,ray2))) then
                if (.not.associated(packet2%r_id_fw)) then            ! if chain connection data is not allocated
                  allocate(packet2%p_id_fw(packet2%n_rays))           ! packet ids forward
                  allocate(packet2%p_id_rev(packet2%n_rays))          ! packet ids reverse
                  allocate(packet2%r_id_fw(packet2%n_rays))           ! ray ids forward
                  allocate(packet2%r_id_rev(packet2%n_rays))          ! ray ids reverse
                  bnd2%b_id_fw    = 0; bnd2%b_id_rev    = 0           ! zero at first
                  packet2%p_id_fw = 0; packet2%p_id_rev = 0           ! zero at first
                  packet2%r_id_fw = 0; packet2%r_id_rev = 0           ! zero at first
                end if
                bnd1%b_id_fw           = bnd2%id                      ! get bundle id fw
                bnd2%b_id_rev          = bnd1%id                      ! get bundle id reverse
                bnd1%rays_fw(packet1%b_id(ray1))     = packet2%b_id(ray2)
                bnd2%rays_rev(packet2%b_id(ray2))    = packet1%b_id(ray1)
                packet1%p_id_fw(ray1)  = packet2%id                   ! get packet id fw
                packet2%p_id_rev(ray2) = packet1%id                   ! get packet id reverse
                packet1%r_id_fw(ray1)  = ray2                         ! get ray id fw
                packet2%r_id_rev(ray2) = ray1                         ! get ray id reverse
              end if
            end do
          end do
        end do
      end do
    end do
  end do

 ! if (io%verbose>3) then
  print *,'========== Direction',dir%id,'==============='
  print *, 'n_bundles:', dir%n_bundles
  do jbn =1, dir%n_bundles
    bnd2 => dir%bundle(jbn)
    print *,jbn,'BUNDLE', bnd2%id
    print *,'dir ID', bnd2%dir_ID
    print *,'axis, theta, phi',bnd2%axis, bnd2%theta, bnd2%phi
    print *,'b_id_fw', bnd2%b_id_fw
    print *,'b_id_rev', bnd2%b_id_rev
    print *,'n_packets',bnd2%n_packets
    print *,'n_rays', bnd2%n_rays
    print *,'lwall',bnd2%lwall
    print *,'uwall',bnd2%uwall
    print *, 'dstep', bnd2%dstep
    print *, 'bnd_rays_fw', bnd2%rays_fw
    print *, 'bnd_rays_rev', bnd2%rays_rev
    print *,'---------------------------------'
    do jpck =1, bnd2%n_packets
      packet2 => bnd2%packet(jpck)
      print *, jpck, 'packet ID', packet2%id
      print *, 'n_rays', packet2%n_rays
      print *, 'n_points',packet2%n_points
      do ray1 = 1, packet2%n_rays
        print *,'-- ray',ray1,'--'
        print *,'b_id', packet2%b_id(ray1)
        print *,'p_id_fw',packet2%p_id_fw(ray1)
        print *,'p_id_rev',packet2%p_id_rev(ray1)
        print *,'r_id_fw',packet2%r_id_fw(ray1)
        print *,'r_id_rev',packet2%r_id_rev(ray1)
        print *,'points'
        print *,'rev_coords', packet2%rev_coords(:,ray1)
        do ray2 = 1, packet2%n_points
          print *,ray2,packet2%coords(:,ray2,ray1)
        end do
        print *,'fwd_coords', packet2%fwd_coords(:,ray1)
      end do
    end do
  end do
 ! end if

  ! tests
  !neq_fw=0
  !neq_rv=0
  !diff_fw=0
  !diff_rv=0
  !do ibn =1, dir%n_bundles
  !  bnd1 => dir%bundle(ibn)
  !  d_bn = bnd1%packet(1)%b_id_fw
  !  do ipck =1, bnd1%n_packets
  !  packet1 => bnd1%packet(ipck)
  !    ray1=packet1%p_id_fw(1)
  !    ray2=packet1%p_id_rev(1)
  !    if (packet1%b_id_fw /= d_bn) print *, 'DIFFERENT BUNDLES FOR DIFFERENT PACKETS IN SAME BUNDLE!!'
  !    if (all(packet1%p_id_fw == ray1)) neq_fw = neq_fw+1
  !    if (all(packet1%p_id_rev == ray2))  neq_rv = neq_rv+1
  !    if (any(packet1%p_id_fw /= ray1)) diff_fw = diff_fw +1
  !    if (any(packet1%p_id_rev /= ray2))  diff_rv = diff_rv +1
  !  end do
  !end do

  !print *, 'same ,different, nbundles', neq_fw, diff_fw, dir%n_bundles
nullify (bnd1, bnd2, packet1, packet2)


END SUBROUTINE finish_direction


!===============================================================================
!> Print summary info for ray_t
!===============================================================================
SUBROUTINE print_ray (self, i)
  class(packet_t) :: self
  integer :: i, j, k
  !.............................................................................
  if (io%verbose < 3) return
  write (*,*) '      RAY_PACKET:', i,'   rays:', self%n_rays
  if (io%verbose < 4) return
  do j=1,self%n_rays
    print '("ray, start, stop:",i4,4(3x,3(i4)))', &
      j, self%fwd_coords(:,j), self%rev_coords(:,j) !, self%pt_wall(:,:,1,j)
    if (io%verbose > 4) then
      do k = 1, self%n_points
        print '("point,  position:",i4,4(3x,3(i4)))', &
      k, self%coords(:,k,j)
      end do
    end if
  end do
END SUBROUTINE print_ray

!===============================================================================
!> Print summary info for ray_bndl_t
!===============================================================================
SUBROUTINE print_bundle (self, i)
  class(ray_bndl_t) :: self
  integer :: j,i
  !.............................................................................
  if (io%verbose < 3) return

  write (*,'(1x,a,i3,2(3x,a,f5.0),2x,a,2(2x,3i3),2x,a,4i5)') &
    '    RAY_BUNDLE:', i, 'theta:', self%theta*radian, 'phi:', self%phi*radian, &
    'lwall, uwall:', self%lwall, self%uwall!, &
    !'packets, dstep:', self%n_packets, self%dstep
  do j=1,self%n_packets
    call self%packet(j)%print(j)

  end do
END SUBROUTINE print_bundle

!===============================================================================
!> Print summary info for ray_dirs_t
!===============================================================================
SUBROUTINE print_dirs (self, i)
  class(ray_dirs_t) :: self
  integer :: i, j
  !.............................................................................
  if (io%verbose < 2) return
  write (*,*) '  RAY_DIR:', i
  write (*,*) '  number of bundles:', self%n_bundles
  do j=1,self%n_bundles
    call self%bundle(j)%print (j)
  end do
END SUBROUTINE print_dirs

!===============================================================================
!> Print summary info for ray_geometry_t
!===============================================================================
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

FUNCTION round(a,digit) result (b)
  real :: a, b, rof
  integer, optional :: digit
! ----------------------------------------------------------------------
  if (present(digit)) then
    rof = 10**digit
  else
    rof = 1000.0
  end if

  b = nint(a*rof)
  b = b/rof
END FUNCTION round

END MODULE ray_geometry_v3_mod
