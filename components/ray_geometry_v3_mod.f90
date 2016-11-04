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
MODULE ray_geometry_v3_mod
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
    real    :: ds
    integer, dimension(:,:,:), pointer   :: coords => null()            ! x|y|z,        ray point, ray #
    integer, dimension(:,:  ), pointer   :: fwd_coords  => null()       ! x|y|z,                   ray #
    integer, dimension(:,:  ), pointer   :: rev_coords  => null()       ! x|y|z,                   ray #
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
    integer, dimension(:), pointer :: ray_chain_fw  => null()
    integer, dimension(:), pointer :: ray_chain_rev => null()
    integer :: ID
    integer :: dir_ID
    integer :: n_packets                                                ! number of packets in this bundle
    integer :: n_rays                                                   ! total number of rays in this bundle
    integer :: n_points
    integer :: l_chain         = 0                                      ! bundle-chain length
    integer :: next_chain_bndl_fw = -1
    integer :: next_chain_bndl_rev = -1
    integer :: lwall(3)
    integer :: uwall(3)
    integer :: dstep(3)                                                 ! traverse index increment
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
    logical :: space_diagonals = .true.
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
  type :: chain_data_t
    integer, dimension(:), pointer :: nr_ID => null()                   ! next ray ID
    integer, dimension(:), pointer :: n_ID => null()                    ! new ID in the final list
    integer, dimension(:), pointer :: p_ID => null()                    ! my packet ID
    integer, dimension(:), pointer :: b_ID => null()                    ! my bundle ID
    integer :: n_rays
  contains
    procedure :: init         => init_chain_data
    procedure :: add_ray      => chain_data_add_ray
    procedure :: destroy      => destroy_chain_data
    procedure :: add_next_ray => chain_data_add_next_ray
    procedure :: new_id       => chain_data_new_ID
    procedure :: packet_id    => chain_data_packet_ID
    procedure :: bundle_id    => chain_data_bundle_ID
    procedure :: finalize_direction
  end type
! ----------------------------------------------------------------------
  type(ray_geometry_t), public :: ray_geometry
  type(chain_data_t) :: chain_data
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
!>   dx, dy, dz     : ray directionality
!> =====================================================================
SUBROUTINE init_cartesian(self)
  class(ray_geometry_t) :: self
  integer                         :: i, d45
  real, dimension(:), allocatable :: xt, yt, zt                         ! FIXME: These seem to not be used
  real, dimension(5)              :: xti=0.0, yti=0.0, zti=0.0
  integer, save                   :: nx=1, ny=0, nz=0
  integer, save                   :: dx=0, dy=0, dz=0
  logical, save                   :: space_diagonals=.true.
! ......................................................................
  namelist /cartesian_ray_geometry_params/ &
    nx, ny, nz, xti, yti, zti, dx, dy, dz, space_diagonals
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
  self%space_diagonals = space_diagonals
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
  bndl%n_rays    = pmax(a)*pmax(b)                                       ! number of rays in this direction
  bndl%l_chain = 1
  bndl%lwall(a) =  0;  bndl%uwall(a) = 0
  bndl%lwall(b) =  0;  bndl%uwall(b) = 0
  bndl%lwall(c) = -1;  bndl%uwall(c) = 1                                ! lower/upper wall for this direction
  bndl%next_chain_bndl_fw  = 1
  bndl%next_chain_bndl_rev = 1

  ! --- Fill in rays info ---
  ray%head  => ray
  ray%tail  => ray
  ray%ID = 1                                                            ! it's packet ID
  ray%n_rays    = pmax(a)*pmax(b)                                       ! number of rays in this direction
  ray%n_points  = pmax(c)                                               ! number of points per ray corresponds to n_cells in that direction
  ray%ds = 1                                                            ! distance the ray travels between points

  allocate (ray%coords(3,pmax(c),ray%n_rays))                           ! create ray table, containing x|y|z loc, for all n_cells and n_rays
  allocate (ray%fwd_coords(3,ray%n_rays))                               ! create table, containing x|y|z coords for upstream position for particular ray, for all rays
  allocate (ray%rev_coords(3,ray%n_rays))                               ! create table, containing x|y|z coords for downstream position for particular ray, for all rays

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
        ray%coords(:,rp,ray_n) = p                                      ! cell position in particular ray
        rp = rp+1                                                       ! increment ray point
      end do
      ray%rev_coords(a,ray_n) = p(a)                                    ! downstream position in terms of first axis
      ray%rev_coords(b,ray_n) = p(b)                                    ! downstream position in terms of second axis
      ray%rev_coords(c,ray_n) = self%ui(c) + 1                          ! downstream position in terms of ray direction
      ray_n = ray_n + 1                                                 ! move to next ray
    end do
  end do
  bndl%dstep = 0;  bndl%dstep(dir) = 1                                  ! bundle dstep
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
  integer :: dir, nphi, a, b, c, pdir
  real :: dp(3), init_step(3), theta, t, phi
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
  real :: init_step(3), dp(3), theta, t, phi
  real, dimension(8) :: phis
  integer :: rep, dir, nphi, a, b, c, pdir
! ......................................................................
  call trace_begin('ray_geometry_t%add_diag',4)

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
      dp(a) = sin(t)*cos(phi)*self%cell_dims(a)                           ! step size along a direction
      dp(b) = sin(t)*sin(phi)*self%cell_dims(b)                           ! step size along b direction
      dp(c) = cos(t)         *self%cell_dims(c)                           ! step size along c direction

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
    !print *,self%theta*radian, self%phi*radian, self%axis, rep
      call self%traverse (dp, a, b, c, init_step)                         ! proceed with 3D traversal

    end do
  end if
  call trace_end
END SUBROUTINE add_diag

!> =====================================================================
!> Subroutine to traverse rays in 3D. Similar to "fast voxel traversal
!> algorithm" (Amanatides & Woo, 1987) and appendix A in Rijkhorst et al. (2006),
!> but without their serious drawbacks and is faster.
!> =====================================================================
SUBROUTINE traverse(self, dp, a, b, c, init_step)
  class(ray_geometry_t) :: self
  class(ray_point_t), pointer :: ray
  class(single_ray_bundle_t), pointer :: ray_bundle
  class(traverse_t), pointer :: t
  integer :: a, b, c, ia, ib, x, r_id, nr_id, or_id
  real :: init_step(3), dp(3)
  logical :: first_point
! ......................................................................
  call trace_begin('traverse_t%traverse')
  call chain_data%init(self%ui)                                         ! initialize chain data storage
  r_id = 0; nr_id = 0                                                   ! nullify ray and next ray IDs
  allocate(t)                                                           ! allocate traversal object
  allocate(ray_bundle)                                                  ! create ray tree to store rays
  call t%init(dp, self%ui, self%li, a, b, c)                            ! initialize traversal object
  ray_bundle%dstep = t%dstep                                            ! save for fast access later
  do ia = t%cell_min(1), t%cell_max(1)                                  ! loop over 1st direction
    do ib = t%cell_min(2), t%cell_max(2)                                ! loop over 2nd direction
      first_point = .true.
      allocate(ray)                                                     ! create new ray
      r_id = r_id + 1                                                   ! increment ray ID
      or_id = r_id
      call chain_data%add_ray                                           ! increment n_rays in chain data
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
          ray%head%rev = t%pos
          ray%head%fwd = t%upstream                                     ! save ray upstream cell

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
              call chain_data%add_next_ray(r_id, or_id)
          else                                                          ! otherwise
            allocate(ray)                                               ! get new ray
            call chain_data%add_ray
            call t%restart                                              ! and restart the 1st or 2nd direction
            nr_id = r_id + 1
            call chain_data%add_next_ray(r_id, nr_id)
            r_id = nr_id
            first_point = .true.
          end if
        end if

      end do
    end do
  end do
  call self%finish_single_ray_list(ray_bundle)                          ! save the ray tree to ray_geometry
  deallocate(t)
  !print *,'CHAIN DATA n_rays', chain_data%n_rays
  !print *,'nr_id', size(chain_data%nr_ID)
  !print *, chain_data%nr_ID
  !print *,'n_id', size(chain_data%n_ID)
  !print *, chain_data%n_ID
  !print *,'p_id', size(chain_data%p_ID)
  !print *, chain_data%p_ID
  !print *,'b_id',size(chain_data%b_ID)
  !print *, chain_data%b_ID
  call chain_data%finalize_direction(self%ray_list%tail%bundle)
  call chain_data%destroy
  !print *, '======'
  !print *, 'bundle ID, dir ID:',self%ray_list%tail%bundle%ID, self%ray_list%tail%bundle%dir_ID
  !print *, 'self%ray_list%tail%bundle%ray_chain_fw',self%ray_list%tail%bundle%ray_chain_fw
  !print *, 'next chain', self%ray_list%tail%bundle%next_chain_bndl_fw
  !print *,'---'
  !print *,'self%ray_list%tail%bundle%ray_chain_rev',self%ray_list%tail%bundle%ray_chain_rev
  !print *, 'prev chain', self%ray_list%tail%bundle%next_chain_bndl_rev
  
END SUBROUTINE traverse

!> =====================================================================
!> INPUT: self, step size along ray direction (dp), smallest cell coord (li)
!> highest cell coord (ui), nominal direction (a, b), ray direction (c)
!> =====================================================================
SUBROUTINE init_traverse(self, dp, ui, li, a, b, c)
  class(traverse_t) :: self
  integer :: a, b, c, li(3), ui(3)
  real :: dp(3)
! ......................................................................
  call trace_begin('traverse_t%init')
  self%p(1) = a; self%p(2) = b; self%p(3) = c                           ! convenience
  self%cell_max(1) = ui(a)                                              ! last cell point along a axis
  self%cell_max(2) = ui(b)                                              ! last cell point along b axis
  self%cell_max(3) = ui(c)                                              ! last cell point along c axis
  self%cell_min(1) = li(a)                                              ! first cell point along a axis
  self%cell_min(2) = li(b)                                              ! first cell point along b axis
  self%cell_min(3) = li(c)                                              ! first cell point along c axis
  self%ds = 1/dp(c)
  self%dstep = 0
  if (dp(a) /= 0) then
    self%dstep(a) = nint(dp(a)*self%ds)                                 ! distance between points along a axis
  end if
  if (dp(b) /= 0) then
    self%dstep(b) = nint(dp(b)*self%ds)                                 ! distance between points along b axis
  end if
  self%dstep(c) = 1                                                     ! distance between points along c axis
  call trace_end
END SUBROUTINE init_traverse

!> =====================================================================
!> INPUT: self, position in nominal directions (ia, ib),
!> distance to cell walls (init_step), upstream cell position for particular
!> ray (upstream)
!> =====================================================================
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

!> =====================================================================
!> restart the ray when it hits a side-wall (invoke pseudo-periodic
!> boundaries). This is for the general, 3D case.
!> =====================================================================
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
  call trace_begin('ray_geometry_t%add_ray_point', 4)
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
  integer :: pt, li(3), ui(3), w
! ......................................................................
  prev => ray%tail
  ray  => ray%head

  self%lwall = ray%lwall
  self%uwall = ray%uwall

  self%n_points = ray%n_points                                          ! store number of points in this ray

  allocate(self%coords(3, ray%n_points))                           ! create a stiff array for this ray with cell cords
  self%fwd_coords = ray%fwd                                             ! save upstream cell coords
  self%rev_coords = ray%rev                                             ! save downstream cell coords
  self%ds = ray%ds
  pt=1
  do while (associated(ray))                                            ! loop over ray points
    self%coords(:,pt)  = ray%coords                                     ! save ray point coords
    prev => ray                                                         ! save previous point
    ray => ray%next                                                     ! go to next point
    deallocate(prev)                                                    ! deallocate previous point
    pt = pt+1
  end do
END SUBROUTINE fill_single_ray_data

SUBROUTINE finish_ray(self, ray_point, id)
  class(single_ray_bundle_t) :: self
  class(single_ray_bundle_t), pointer :: bundle
  class(single_ray_packet_t), pointer :: packet
  class(single_ray_t), pointer  :: ray
  class(ray_point_t),  pointer  :: ray_point
  integer :: id
  logical :: fnd_p, fnd_b
! ......................................................................
  call trace_begin('ray_geometry_t%finish_ray')

  allocate (ray)                                                        ! allocate ray list
  ray%id = id
  call ray%fill_single_ray_data(ray_point)

  if (.not.associated(self%head)) then                                  ! if the ray tree is empty
    allocate(bundle)
    bundle%dstep = self%dstep
    self%head => bundle
    self%tail => bundle
    bundle%lwall = ray%lwall
    bundle%uwall = ray%uwall
    bundle%bundle_ID = 1
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
      bundle%dstep = self%dstep
      bundle%head => self%head
      bundle%prev => self%tail
      bundle%prev%next => bundle
      self%tail => bundle
      bundle%lwall = ray%lwall
      bundle%uwall = ray%uwall
      bundle%bundle_ID = bundle%prev%bundle_ID + 1
    end if
  end if

  if (.not.associated(bundle%packets)) then
    allocate(packet)
    bundle%packets => packet
    packet%head => packet
    packet%tail => packet
    packet%n_rays = 1
    packet%n_points = ray%n_points
    packet%packet_ID = 1
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
      packet%packet_ID = packet%prev%packet_ID + 1
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
  call chain_data%bundle_ID(ray%id,bundle%bundle_ID)
  call chain_data%packet_ID(ray%id,packet%packet_ID)
  nullify (ray, bundle, packet, ray_point)
  call trace_end
END SUBROUTINE finish_ray

SUBROUTINE fill_ray_list_data(self, sray,r_id)
  class(ray_t) :: self
  class(single_ray_t), pointer :: sray, prev_sray
  integer :: iray, r_id
! ......................................................................
  allocate (self%coords(3,self%n_points,self%n_rays))              ! cell coords per point per ray in the bundle
  allocate (self%fwd_coords(3,self%n_rays))                             ! upstream cell coords per ray
  allocate (self%rev_coords(3,self%n_rays))                             ! downstream cell coords per ray

  sray => sray%head
  iray = 0
  self%ds = sray%ds
  do while (associated(sray))                                           ! loop over all rays in the tree
    iray = iray + 1
    r_id = r_id +1
    self%coords(:,:,iray) = sray%coords                                 ! save cell cords
    self%fwd_coords(:,iray)    = sray%fwd_coords                        ! save upstream coords
    self%rev_coords(:,iray)    = sray%rev_coords                        ! save downstream coords
    call chain_data%new_id(sray%id,r_id)
    prev_sray   => sray                                                 ! hold pointer to this ray
    deallocate(sray%coords)
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
  integer :: r_id
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
    r_id = 0
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
    else
      bundle%head => dir%bundle
      bundle%prev => dir%bundle%tail
      bundle%prev%next => bundle
      dir%bundle%tail => bundle
    end if
    bundle%ID = sbundle%bundle_ID
    bundle%lwall = sbundle%lwall
    bundle%uwall = sbundle%uwall
    bundle%dstep = sbundle%dstep
    bundle%n_packets = 0
    bundle%n_rays = 0
    ! --- Go through single packets ---
    spacket => sbundle%packets
    do while (associated(spacket))
      bundle%n_packets = bundle%n_packets + 1
      allocate(ray)
      if (.not.associated(bundle%rays)) then
        bundle%rays => ray
        ray%head => ray
        ray%tail => ray
      else
        ray%head => bundle%rays
        ray%prev => bundle%rays%tail
        ray%prev%next => ray
        bundle%rays%tail => ray
      end if
      ray%ID = spacket%packet_ID
      ray%n_rays = spacket%n_rays
      bundle%n_rays = bundle%n_rays + spacket%n_rays
      ray%n_points = spacket%n_points
      call ray%fill_ray_list_data(spacket%rays, r_id)
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

!> =============================================================================
!> Print summary info for ray_bndl_t
!> =============================================================================
SUBROUTINE print_bundle (self, i)
  class(ray_bndl_t) :: self
  integer :: i, j, n
  class(ray_t), pointer :: ray
  !.............................................................................
  if (io%verbose < 3) return
  n = 0
  ray => self%rays
  do while (associated(ray))
    n = n+ray%n_points
    ray => ray%next
  end do
  self%n_points = n
  write (*,'(1x,a,i3,2(3x,a,f5.0),2x,a,2(2x,3i3),2x,a,4i5)') &
    '    RAY_BUNDLE:', i, 'theta:', self%theta*radian, 'phi:', self%phi*radian, &
    'lwall, uwall:', self%lwall, self%uwall, &
    'packets, dstep:', self%n_packets, self%dstep
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
    
SUBROUTINE init_chain_data(self, n)
  class(chain_data_t) :: self
  integer :: n(3)
! ----------------------------------------------------------------------
  allocate(self%nr_ID(n(1)*n(2)*3))
  allocate(self%n_ID(n(1)*n(2)*3))
  allocate(self%b_ID(n(1)*n(2)*3))
  allocate(self%p_ID(n(1)*n(2)*3))
  self%nr_id = 0; self%n_id = 0; self%b_id = 0; self%p_Id = 0
  self%n_rays = 0
END SUBROUTINE init_chain_data

SUBROUTINE destroy_chain_data(self)
  class(chain_data_t) :: self
! ----------------------------------------------------------------------
  deallocate(self%nr_ID)
  deallocate(self%n_ID)
  deallocate(self%b_ID)
END SUBROUTINE destroy_chain_data

SUBROUTINE chain_data_add_ray(self)
  class(chain_data_t) :: self
! ----------------------------------------------------------------------
  self%n_rays = self%n_rays+1
END SUBROUTINE chain_data_add_ray

SUBROUTINE chain_data_add_next_ray(self,m,n)
  class(chain_data_t) :: self
  integer :: m, n
! ----------------------------------------------------------------------
  self%nr_ID(m) = n
END SUBROUTINE chain_data_add_next_ray

SUBROUTINE chain_data_new_ID (self, m, n)
  class(chain_data_t) :: self
  integer :: m, n
! ----------------------------------------------------------------------
  self%n_ID(m) = n
END SUBROUTINE chain_data_new_ID

SUBROUTINE chain_data_packet_ID(self, m, n)
  class(chain_data_t) :: self
  integer :: m, n
! ----------------------------------------------------------------------
  self%p_ID(m) = n
END SUBROUTINE chain_data_packet_ID

SUBROUTINE chain_data_bundle_ID(self, m, n)
  class(chain_data_t) :: self
  integer :: m, n
! ----------------------------------------------------------------------
  self%b_ID(m) = n
END SUBROUTINE chain_data_bundle_ID

SUBROUTINE finalize_direction(self, bundle)
  class(chain_data_t) :: self
  class(ray_bndl_t), pointer   :: bundle, bnd, bndr
  integer :: ray
  logical :: fnd, fnd2
! ----------------------------------------------------------------------
  do ray = 1, self%n_rays
    ! find current bundle
    bnd => bundle
    fnd = .false.
    do while(associated(bnd).and..not.fnd)
      if (self%b_ID(ray) == bnd%id) then
        fnd = .true.
      else  
        bnd => bnd%next
      end if
    end do
    
    ! find forward bundle to store reverse data there
    bndr => bundle
    fnd2 = .false.
    do while(associated(bndr).and..not.fnd2)
      if (self%b_ID(self%nr_ID(ray)) == bndr%id) then
        fnd2 = .true.
      else  
        bndr => bndr%next
      end if
    end do
    
    if (fnd) then
      if (.not.(associated(bnd%ray_chain_fw))) allocate(bnd%ray_chain_fw(bnd%n_rays))
      if (self%nr_ID(ray) /= 0) then
        bnd%ray_chain_fw(self%n_ID(ray)) = self%n_ID(self%nr_ID(ray))   ! forward
        if (bnd%next_chain_bndl_fw < 0) then
          bnd%next_chain_bndl_fw = self%b_ID(self%nr_ID(ray))           ! forward
        else
          if (bnd%next_chain_bndl_fw /= self%b_ID(self%nr_ID(ray))) then
            error stop "finalize direction: b_ID mismatch"
          end if
        end if
        if (self%nr_id(ray) == ray) then
          bnd%l_chain =1
        else if (self%nr_id(self%nr_id(ray)) == ray) then 
          bnd%l_chain = 2
        else if (self%nr_id(self%nr_id(self%nr_id(ray))) == ray) then
          bnd%l_chain = 3
        else if (self%nr_id(self%nr_id(self%nr_id(self%nr_id(ray)))) == ray) then
          bnd%l_chain = 4
        else if (self%nr_id(self%nr_id(self%nr_id(self%nr_id(self%nr_id(ray))))) == ray) then
          bnd%l_chain = 5
        else
          error stop "finalize direction: when will it end???"
        end if 
      else  
        bnd%l_chain = 1
      end if
    else
      error stop "finalize direction: bundle not found by b_ID"
    end if
    
    if (fnd2) then
      if (.not.(associated(bndr%ray_chain_rev))) allocate(bndr%ray_chain_rev(bndr%n_rays))
      if (self%nr_ID(ray) /= 0) then
        bndr%ray_chain_rev(self%n_ID(self%nr_ID(ray))) = self%n_ID(ray) ! reverse
        if (bndr%next_chain_bndl_rev < 0) then
          bndr%next_chain_bndl_rev = self%b_ID(ray)                     ! reverse
        else
          if (bndr%next_chain_bndl_rev /= self%b_ID(ray)) then
            error stop "finalize direction: rev b_ID mismatch"
          end if
        end if
      end if
    else
      error stop "finalize direction: fw bundle not found by b_ID"
    end if
    
  end do
END SUBROUTINE finalize_direction

END MODULE ray_geometry_v3_mod
