MODULE ray_geometry_mod
  USE io_mod
  USE trace_mod
  implicit none
  type, public :: cell_t
    integer :: ray, point, lwall(3), uwall(3)
    class(cell_t), pointer :: head => null()
    class(cell_t), pointer :: next => null()
    class(cell_t), pointer :: tail => null()
  end type
    
  type, public                         :: ray_t
    integer                            :: n_rays
    integer                            :: lwall(3)
    integer                            :: uwall(3)
    integer, dimension(:    ), pointer :: n_points    => null()
    integer, dimension(:,:,:), pointer :: cell_coords => null()  ! x|y|z, ray point, ray #
    integer, dimension(:,:  ), pointer :: fwd_coords             ! x|y|z,            ray #
    integer, dimension(:,:  ), pointer :: rev_coords             ! x|y|z,            ray #
    real,    dimension(:,:  ), pointer :: ds          => null()  ! ds   , ray point, ray #
    integer, dimension(:,:,:,:), pointer :: pt_wall     => null()  ! lw|uw, ray point, ray #
    type(cell_t), dimension(:,:,:), pointer :: cell_domain => null()
  contains
    procedure :: store_cell_data
  end type
! ----------------------------------------------------------------------
  type, extends(ray_t), public:: ray_list_t
    class(ray_list_t), pointer :: head         => null()
    class(ray_list_t), pointer :: next_dir     => null()
    class(ray_list_t), pointer :: prev_dir     => null()
    class(ray_list_t), pointer :: last_dir     => null()
    class(ray_list_t), pointer :: first_bundle => null()
    class(ray_list_t), pointer :: next_bundle  => null()
    class(ray_list_t), pointer :: prev_bundle  => null()
    class(ray_list_t), pointer :: last_bundle  => null()
    integer                    :: max_n_points
    integer                    :: dir_ID
    integer                    :: bundle_ID
  contains
    procedure :: invert_my_rays
  end type
! ----------------------------------------------------------------------
  type, public :: ray_geometry_t
    integer                    :: n_bundles
    integer                    :: ray_dirs
    integer                    :: li(3), ui(3)
    integer                    :: cell_dims(3)
    integer                    :: geom = 1 !cartesian = 1, spherical = 2, cylindrical = 3
    logical                    :: init_done = .false.
    class(ray_list_t), pointer :: ray_list => null()
  contains
    procedure :: add_axis
    procedure :: add_diag
    procedure :: add_ray
    procedure :: get_my_rays
    procedure :: finish_single_ray_list
    procedure :: test_ray_geometry
  end type
! ----------------------------------------------------------------------
  type, extends(ray_geometry_t), public :: ray_geometry_list_t
    class(ray_geometry_list_t), pointer :: head => null()
    class(ray_geometry_list_t), pointer :: next => null()
    class(ray_geometry_list_t), pointer :: prev => null()
    class(ray_geometry_list_t), pointer :: tail => null()
  contains
    procedure :: init => init_cartesian
  end type
! ----------------------------------------------------------------------
  type, public                       :: single_ray_t
    integer                          :: lwall(3)
    integer                          :: uwall(3)
    integer                          :: n_points
    integer, dimension(:,:), pointer :: cell_coords => null()  ! x|y|z, ray point
    integer, dimension(3)            :: fwd_coords             ! x|y|z,           
    integer, dimension(3)            :: rev_coords             ! x|y|z,          
    real,    dimension(:  ), pointer :: ds          => null()  ! ds   , ray point
    integer, dimension(:,:,:), pointer :: cell_walls => null()  ! lw|uw(x,y,z), ray point
  end type
! ----------------------------------------------------------------------
  type, extends(single_ray_t), public:: single_ray_list_t
    class(single_ray_list_t), pointer :: head        => null()
    class(single_ray_list_t), pointer :: first_ray   => null()
    class(single_ray_list_t), pointer :: next_ray    => null()
    class(single_ray_list_t), pointer :: prev_ray    => null()
    class(single_ray_list_t), pointer :: last_ray    => null()
    class(single_ray_list_t), pointer :: next_bundle => null()
    class(single_ray_list_t), pointer :: prev_bundle => null()
    class(single_ray_list_t), pointer :: last_bundle => null()
    integer                           :: n_rays
    integer                           :: max_n_points
    integer                           :: dir_ID
    integer                           :: bundle_ID
  contains
    procedure :: finish_ray
  end type
! ----------------------------------------------------------------------
  type                         :: ray_point_t
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
  type(ray_geometry_list_t), public :: ray_geometry
  real(8), parameter:: pi = asin(1.0d0) * 2.0d0, radian = 180./pi
CONTAINS
!> =====================================================================
!> Initialisation subroutine
!> =====================================================================
SUBROUTINE init_cartesian(self,li,ui, cell_dims)
  class(ray_geometry_list_t)      :: self
  integer                         :: i, li(3), ui(3)
  integer, optional               :: cell_dims(3)
  integer                         :: n_dirs
  real, dimension(:), allocatable :: theta
  real, dimension(3)              :: thetas
! ......................................................................
  namelist /ray_geometry_params/ n_dirs, thetas
! ......................................................................
  call trace_begin('ray_geometry_t%init')
  !> TO DO: add an option for 2D simulation setup (i.e. any(cell_dims) == 0)
  if(.not.self%init_done) then
    print *,'INITIALISING RAY GEOMETRY'
    self%init_done = .true.
    allocate(self%ray_list)
    self%ray_dirs = 0
    self%n_bundles = 0
    rewind (io%input)
    read (io%input, ray_geometry_params)
    self%li = li
    self%ui = ui
    if(present(cell_dims)) then
      self%cell_dims = cell_dims
    else
      self%cell_dims(:) = 1               ! default - cubic cells
    end if
    allocate (theta(n_dirs))
    theta(1:n_dirs)=thetas(1:n_dirs)
    if (n_dirs == 0) then
      error stop "ERROR: no rays created. Increase theta to at least 1 angle."
    end if
    do i = 1, n_dirs
      if (theta(i) == 0) then    
        call self%add_axis   
      else if (theta(i) == 45.) then
        call self%add_diag
      else
        call self%add_ray(theta(i))
      end if
    end do
  end if 
!> TO DO: if init already done, but n_cells is different, create new geometry
  call trace_end
END SUBROUTINE init_cartesian
    

!> =====================================================================
!> Subroutine to add axisymmetric rays. It is done separately, because
!> it is much simpler and therefore faster than with rays with other
!> angles. Adds stiff ray arrays instantly.
!> =====================================================================
SUBROUTINE add_axis(self)
  class(ray_geometry_t)     :: self
  type(ray_list_t), pointer :: ray 
  integer                   :: ia, ib, ic, ray_n, dir, a, b, c, rp
  real                      :: ds
  integer     :: p(3), pmax(3)
! ......................................................................
  call trace_begin('ray_geometry_t%add_axis')
  pmax = (self%ui - self%li) + 1                         ! Total number of points
  do dir = 1, 3                                          ! Do to all 3 directions
    self%ray_dirs = self%ray_dirs + 1                    ! Adding new dir
    SELECT case(dir)
    case(1)
      a = 1; b = 2; c = 3                                ! along z axis
    case(2)
      a = 3; b = 1; c = 2                                ! along y axis
    case(3)
      a = 2; b = 3; c = 1                                ! along x axis
    end SELECT
    ! Create new ray bundle and add it to ray_geometry
    allocate (ray)                                       ! create new ray
    self%n_bundles   = self%n_bundles + 1                ! adding new bundle
    ray%dir_ID       = self%ray_dirs                     ! it's ID corresponds to it's direction
    ray%bundle_ID    = self%n_bundles                    ! it's bundle_ID
    ray%n_rays       = pmax(a)*pmax(b)                   ! number of rays in this direction
    ray%max_n_points = pmax(c)
    allocate(ray%cell_domain(self%ui(a),self%ui(b),self%ui(c)))
    if (associated(self%ray_list%head)) then             ! if ray list is not empty
      ray%prev_dir => self%ray_list%last_dir
      ray%prev_dir%next_dir  => ray
      self%ray_list%last_dir => ray
      ray%head => self%ray_list%head
    else
      self%ray_list => ray
      self%ray_list%head => ray                          ! if empty, point to itself
      self%ray_list%last_dir => ray
    end if
    ray%first_bundle => ray
    ray%last_bundle  => ray
    ! ---
    allocate (ray%cell_coords(3,pmax(c),ray%n_rays))     ! create ray table, containing x|y|z loc, for all n_cells and n_rays
    allocate (ray%ds(pmax(c),ray%n_rays))                ! create ds table, containing ds for all points in all rays
    allocate (ray%fwd_coords(3,ray%n_rays))              ! create table, containing x|y|z coords for upstream position for particular ray, for all rays
    allocate (ray%rev_coords(3,ray%n_rays))              ! create table, containing x|y|z coords for downstream position for particular ray, for all rays
    allocate (ray%n_points(ray%n_rays))                  ! number of points in the rays of this ray list
    allocate (ray%pt_wall(3,2,pmax(c), ray%n_rays))      ! lower and upper walls the ray crosses in the cell
    ray%n_points(:) = pmax(c)                            ! all n_points are equal
    ray%lwall(a) = 0
    ray%uwall(a) = 0
    ray%lwall(b) = 0
    ray%uwall(b) = 0
    ray%lwall(c) = -1                                    ! lower wall for this direction (in bi-directional aspect)
    ray%uwall(c) = 1                                     ! upper wall for this direction (in bi-directional aspect)
    ray_n = 1                                            ! ray number
    do ia = self%li(a), self%ui(a)                       ! iterate over cells in first axis
      do ib = self%li(b), self%ui(b)                     ! iterate over cells in second axis
        p(a) = ia; p(b) = ib                             ! point position in first two axes
        ray%fwd_coords(a,ray_n) = p(a)                   ! upstream position in terms of first axis
        ray%fwd_coords(b,ray_n) = p(b)                   ! upstream position in terms of second axis
        ray%fwd_coords(c,ray_n) = self%li(c) - 1         ! upstream position in terms of ray direcion
        rp = 1
        do ic = self%li(c), self%ui(c)                   ! move along ray direction
          p(c) = ic                                      ! current cell position in terms of third axis
          ray%cell_coords(:,rp,ray_n) = p                ! cell position in particular ray
          ray%ds(rp,ray_n)          = self%cell_dims(c)  ! ray patch length in particular cell
          call ray%store_cell_data(p,rp,ray_n, ray%lwall, ray%uwall)
          ray%pt_wall(:,1,rp,ray_n) = ray%lwall
          ray%pt_wall(:,2,rp,ray_n) = ray%uwall
          rp = rp+1
        end do
        ray%rev_coords(a,ray_n) = p(a)                   ! downstream position in terms of first axis
        ray%rev_coords(b,ray_n) = p(b)                   ! downstream position in terms of second axis
        ray%rev_coords(c,ray_n) = self%ui(c) + 1         ! downstream position in terms of ray direction
        ray_n = ray_n + 1                                ! move to next ray
      end do
    end do
  end do
  call trace_end
END SUBROUTINE add_axis

!> =====================================================================
!> Subroutine to add rays (except 0 and 45 degree angles). 
!> Done using modified "fast voxel traversal" algorithm 
!> (original, not modified, one is from Amanatides & Woo, 1987)
!> =====================================================================
SUBROUTINE add_ray(self, theta)
  class(ray_geometry_t)           :: self
  real                            :: theta, t, phi
  integer                         :: dir, nphi, a, b, c, pdir, round
  real, dimension(3)              :: dp, init_step
  real, dimension(2)              :: dp2d
  real, dimension(:), allocatable :: phis
  logical :: oob_dir
  
  ! ....................................................................
  call trace_begin('ray_geometry_t%add_ray')
  nphi=floor(360./theta)                                 ! equispaced solid angles for all rays
  allocate(phis(nphi))
  phis(1) = 0; t = 360./nphi                             ! phi = 0 direction

  do dir = 2, nphi-1                                     ! create phis array, avoiding phi = 360
    phis(dir) = phis(dir-1) + t
  end do
  
  phis(nphi) = phis(nphi-1) + t
  phis = phis/radian
  t = theta/radian
  
  do pdir = 1, nphi                                      ! loop over phis
    phi = phis(pdir)
    do dir = 1, 3                                        ! loop over x|y|z orientations
      self%ray_dirs = self%ray_dirs + 1                  ! increment ray dirs
      SELECT case(dir)
      case(1)
        a = 1; b = 2; c = 3                              ! along z
      case(2) 
        a = 3; b = 1; c = 2                              ! along y
      case(3)
        a = 2; b = 3; c = 1                              ! along x
      end SELECT
      dp(a) = sin(t)*cos(phi)*self%cell_dims(a)          ! step size along a direction
      dp(b) = sin(t)*sin(phi)*self%cell_dims(b)          ! step size along b direction
      dp(c) = cos(t)         *self%cell_dims(c)          ! step size along c direction
      
      round = anint(dp(a)*100.0)                         ! removing computer uncertainties
      dp(a) = round/100.0                                ! and giving angle error of around 1 degree
                                                         ! this is necessary in order to avoid infinite loops
      round = anint(dp(b)*100.0)                         ! and overabundance of ray points with machine epsilon 
      dp(b) = round/100.0                                ! ds
      
      round = anint(dp(c)*100.0)
      dp(c) = round/100.0

      if (dp(a)==0) then                                 ! if a direction == 0
        dp2d(1) = dp(b)                                  ! prepare b direction for 2D traversal
        dp2d(2) = dp(c)                                  ! prepare c direction for 2D traversal
        init_step(1) = 0.5                               ! cell center for b direction
        init_step(2) = 1                                 ! cell wall for c direction (ray propagation direction)
        init_step(3) = 0.5                               ! passive cell step
        call traverse_2d(self, dp2d, b, c, a, init_step) ! proceed with 2D traversal
        
      else if (dp(b) == 0) then                          ! if b direction == 0
        dp2d(1) = dp(a)                                  ! prepare a direction for 2D traversal
        dp2d(2) = dp(c)                                  ! prepare c direction for 2D traversal
        init_step(1) = 0.5                               ! cell center for a direction
        init_step(2) = 1                                 ! cell wall for c direction (ray propagation direction)
        init_step(3) = 0.5                               ! passive cell step
        call traverse_2d(self, dp2d, a, c, b, init_step) ! proceed with 2D traversal
        
      else                                               ! movement to all directions
        init_step(1) = 0.5                               ! cell center for a direction
        init_step(2) = 0.5                               ! cell center for b direction
        init_step(3) = 1                                 ! cell wall for c direction (ray propagation direction)
        call traverse_3d (self, dp, a, b, c, init_step)  ! proceed with 3D traversal
      end if
    end do
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> IMPORTANT: TO DO: implementing non-cubic cells - 45 degree rays will NOT
!> necessarily have constant ds. Implementation
!> should resemble add_ray, but with strict follow-up of ray directions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> =====================================================================
SUBROUTINE add_diag(self)
  class(ray_geometry_t)      :: self
  real                       :: theta, t, phi
  integer                    :: dir, nphi, a, b, c, pdir, round
  real, dimension(3)         :: init_step, dp
  real, dimension(8)         :: phis
  real, dimension(2)         :: dp2d
! ......................................................................
  call trace_begin('ray_geometry_t%add_diag')
  t = 45./radian
  do dir = 1, 2                         ! follow only 2 directions - the 3rd
    SELECT case(dir)                    ! one is overlapped by these 2
    case(1)
      a = 1; b = 2; c = 3
      nphi = 8
      phis(1:8) = (/0., 45., 90., 135., 180., 225., 270., 315./)
    case(2)
      a = 3; b = 1; c = 2
      nphi = 2
      phis(1:2) = (/45., 315./)                          ! the rest phis are overlapped by 
    end SELECT                                           ! 1st direction
    phis = phis/radian
    do pdir = 1, nphi                                    ! loop over phis
      self%ray_dirs = self%ray_dirs + 1                  ! increment number of ray directions
      phi=phis(dir)
      dp(a) = sin(t)*cos(phi)*self%cell_dims(a)          ! step size along a direction
      dp(b) = sin(t)*sin(phi)*self%cell_dims(b)          ! step size along b direction
      dp(c) = cos(t)         *self%cell_dims(c)          ! step size along c direction
      
      round = anint(dp(a)*100.0)                         ! removing computer uncertainties
      dp(a) = round/100.0                                ! and giving angle error of around 1 degree
      round = anint(dp(b)*100.0)                         ! this is necessary in order to avoid infinite loops
      dp(b) = round/100.0                                ! and overabundance of ray points with machine epsilon 
      round = anint(dp(c)*100.0)                         ! ds
      dp(c) = round/100.0
      
      if (dp(a)==0) then                                 ! if a direction == 0
        dp2d(1) = dp(b)                                  ! prepare b direction for 2D traversal
        dp2d(2) = dp(c)                                  ! prepare c direction for 2D traversal
        init_step(1) = 0.5                               ! cell center for a direction
        init_step(2) = 1                                 ! cell wall for c direction (ray propagation direction)
        init_step(3) = 0.5                               ! passive cell step
        call traverse_2d(self, dp2d, b, c, a, init_step) ! proceed with 2D traversal
        
      else if (dp(b) == 0) then                          ! if b direction == 0
        dp2d(1) = dp(a)                                  ! prepare a direction for 2D traversal
        dp2d(2) = dp(c)                                  ! prepare c direction for 2D traversal
        init_step(1) = 0.5                               ! cell center for a direction
        init_step(2) = 1                                 ! cell wall for c direction (ray propagation direction)
        init_step(3) = 0.5                               ! passive cell step
        call traverse_2d(self, dp2d, a, c, b, init_step) ! proceed with 2D traversal
        
      else                                               ! movement to all directions
        init_step(1) = 0.5                               ! cell center for a direction
        init_step(2) = 0.5                               ! cell center for b direction
        init_step(3) = 1                                 ! cell wall for c direction (ray propagation direction)
        call traverse_3d (self, dp, a, b, c, init_step)  ! proceed with 3D traversal
      end if
      
    end do
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
  class(ray_geometry_t)              :: self
  class(ray_point_t), pointer        :: ray
  class(single_ray_list_t) , pointer :: ray_bundle
  class(traverse_2d_t), pointer      :: t
  integer :: a, b, fc, ia, ib
  real    :: dp(2), init_step(3), I_up(3), lw(3)
  logical :: first_point
! ......................................................................
call trace_begin('traverse_2d_t%traverse_2d')
  allocate(t)                                            ! allocate traversal object
  allocate(ray_bundle)                                   ! create ray tree to store rays
  ray_bundle%dir_ID = self%ray_dirs                        ! identify ray_tree
  call t%init(dp, self%ui, self%li, a, b, fc)            ! initialize traversal object
  
  do ia = t%cell_min(1), t%cell_max(1), t%cell_step(1)   ! loop over nominal direction
    do ib = t%cell_min(3), t%cell_max(3), t%cell_step(3) ! loop over passive direction
      first_point = .true.
      allocate(ray)                                      ! create new ray
      call t%posit(ia, ib, init_step, I_up)              ! position traversal object                                
      do while (.not.t%oob)                              ! loop until ray reaches wall at ray direction
        call ray%add_ray_point(t%cell_pos)               ! add a ray point to the ray
        if (first_point) then
          ray%head%lwall = t%lw
          first_point = .false.
        else 
          ray%tail%lwall = t%lw
          ray%tail%prev%uwall = -ray%tail%lwall
        end if
        call t%hop                                       ! hop to the next cell
        ray%tail%ds = t%ds                               ! save distance traveled
        if (t%r_trace >= t%r_tot) then                   ! if the ray is traced
          ray%tail%uwall = -t%lw
          call t%get_downstream(ray%head%rev)            ! get ray downstream cell
          ray%head%fwd = I_up                            ! save ray upstream cell
          call ray_bundle%finish_ray(ray)                ! and save ray to the ray tree
          if (t%cell_pos(b) > t%cell_max(2)) then        ! if cell wall along the ray direction is reached
            t%oob = .true.                               ! this ray is done
          else                                           ! otherwise
            allocate(ray)                                ! get new ray
            first_point = .true.
            call t%restart(I_up)                         ! and restart the nominal direction
          end if
        end if
          
      end do
    end do
  end do
  call self%finish_single_ray_list(ray_bundle)           ! save the ray tree to ray_geometry
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
  class(single_ray_list_t), pointer :: ray_bundle
  class(traverse_3d_t),     pointer :: t
  integer            :: a, b, c, ia, ib
  real, dimension(3) :: init_step, dp, I_up
  logical :: first_point
! ......................................................................
  call trace_begin('traverse_3d_t%traverse_3d')
  allocate(t)                                            ! allocate traversal object
  allocate(ray_bundle)                                   ! create ray tree to store rays
  ray_bundle%dir_ID = self%ray_dirs                      ! identify ray_tree
  call t%init(dp, self%ui, self%li, a, b, c)             ! initialize traversal object
 
  do ia = t%cell_min(1), t%cell_max(1), t%cell_step(1)   ! loop over 1st direction
    do ib = t%cell_min(2), t%cell_max(2), t%cell_step(2) ! loop over 2nd direction
      first_point = .true.
      allocate(ray)                                      ! create new ray
      call t%posit(ia, ib, init_step, I_up)              ! position traversal object
      do while(.not.t%oob)                               ! loop until ray reaches wall at ray direction
        call ray%add_ray_point(t%cell_pos)               ! add a ray point to the ray
        if (first_point) then
          ray%head%lwall = t%lw
          first_point = .false.
        else 
          ray%tail%lwall = t%lw
          ray%tail%prev%uwall = -ray%tail%lwall
        end if
        call t%hop                                       ! hop to the next cell
        ray%tail%ds = t%ds                               ! save distance traveled
                
        if (t%r_trace >= t%r_tot) then                   ! if the ray is traced
          ray%tail%uwall = -t%lw
          call t%get_downstream(ray%head%rev)            ! get ray downstream cell
          ray%head%fwd = I_up                            ! save ray upstream cell
          call ray_bundle%finish_ray(ray)                ! and save ray to the ray tree
          if (t%cell_pos(c) >= t%cell_max(3)) then       ! if cell wall along the ray direction is reached
              t%oob = .true.                             ! this ray is done
          else                                           ! otherwise
            allocate(ray)                                ! get new ray
            call t%restart(I_up)                         ! and restart the 1st or 2nd direction
            first_point = .true.
          end if
        end if
        
      end do
    end do
  end do
  call self%finish_single_ray_list(ray_bundle)           ! save the ray tree to ray_geometry
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
  self%cell_max(1) = ui(a)                           ! last cell point along a axis
  self%cell_max(2) = ui(b)                           ! last cell point along b axis
  self%cell_max(3) = ui(fc)                          ! last cell point along fixed coord axis
  if (dp(1) < 0) then                                ! if ray goes from "top" to "bottom"
    self%cell_min (1) =  self%cell_max(1)            ! swap cell_max and cell_min
    self%cell_max (1) =  li(a)                       ! cell_max now is first cell point
    self%cell_step(1) = -1                           ! stepper now is negative
    dp(1) = -dp(1)                                   ! step size always has to be positive
  else
    self%cell_step(1) = 1                            ! step through one cell to another
    self%cell_min (1) = li(a)
  end if
  self%dstep = 1/dp                                  ! distance the ray travels to go one step
  self%cell_min(2)  = li(b);  self%cell_step(2) = 1  ! ray direction is always positive
  self%cell_min(3)  = li(fc); self%cell_step(3) = 1  ! fixed coordinate
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
  self%cell_pos(self%a ) = ia                                   ! position in a axis
  self%cell_pos(self%b ) = self%cell_min(2)                     ! position in ray direction axis
  self%cell_pos(self%fc) = ib                                   ! position in passive axis
  self%lw = 0; self%uw = 0
  
  if (init_step(1) == 1) then                                   ! initial ray position in a axis
    self%ray_p_init(1) = ia - self%cell_step(1)/2               ! it starts from the wall
  else
    self%ray_p_init(1) = ia                                     ! or at predetermined position
  end if
  
  if (init_step(2) == 1) then                                   ! initial ray position in ray direction
    self%ray_p_init(2) = self%cell_min(2) - self%cell_step(2)/2 ! it starts from the wall
  else
    self%ray_p_init(2) = self%cell_min(2)                       ! or at predetermined position
  end if
  
  self%ray_p_last(1) = self%cell_max(1) + self%cell_step(1)/2   ! ray ends at cell wall
  self%ray_p_last(2) = self%cell_max(2) + self%cell_step(2)/2   ! ray ends at cell wall
    
  if (self%cell_step(1) > 0) then                               ! distance ray travels in a direction
    self%ray_dist(1) = self%ray_p_last(1) - self%ray_p_init(1)  
  else
    self%ray_dist(1) = self%ray_p_init(1) - self%ray_p_last(1)
  end if
  
  self%ray_dist(2) = self%ray_p_last(2) - self%ray_p_init(2)    ! distance ray travels in ray direction axis
  self%ray_dist    = self%ray_dist * self%dstep                 ! total distance ray travels in each direction to reach the end of patch
  self%r_tot       = minval(self%ray_dist)                      ! determine the smallest distance
  self%ray_step(1) = init_step(1)
  self%ray_step(2) = init_step(2)
  self%ray_step    = self%ray_step * self%dstep                 ! distance along the ray to cross to another cell
  
  if (init_step(1) == init_step(2)) then                        ! find upstream cell position
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
  self%r_trace = 0.                                             ! distance, the ray has traveled
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
  if (ray_step(1) < ray_step(2)) then          ! if cell wall in direction 1 is closer than in
    cell_pos(a) = cell_pos(a) + cell_step(1)   ! direction 2, hop to cell towards direction 1
    ray_step(2) = ray_step(2) - ray_step(1)    ! reduce the distance to wall in 2 with distance traveled
             ds = ray_step(1)                  ! save distance traveled
    r_trace     = r_trace + ray_step(1)        ! save total distance traveled (for sanity check)
    ray_step(1) = dstep(1)                     ! reset distance to wall for 1
    self%lw(a)  = -1 
  else if (ray_step(1) == ray_step(2)) then    ! if both equal, move diagonally
    cell_pos(a) = cell_pos(a) + cell_step(1)   ! cell position is diagonal to the old one
    cell_pos(b) = cell_pos(b) + cell_step(2)
             ds = ray_step(1)                  ! distance of both 1 and 2 is the same
    r_trace     = r_trace + ray_step(1)        ! save total distance traveled
    ray_step    = dstep                        ! reset distance to wall for 1 and 2
    self%lw(a)  = -1
    self%lw(b)  = -1
  else
    cell_pos(b) = cell_pos(b) + cell_step(2)   ! if wall is closer in 2, move towards there
    ray_step(1) = ray_step(1) - ray_step(2)    ! reduce the distance to 1
             ds = ray_step(2)                  ! save distance traveled
    r_trace     = r_trace + ray_step(2)        ! save total distance traveled
    ray_step(2) = dstep(2)                     ! reset distance to wall for 2
    self%lw(b)  = -1
  end if
  end associate
END SUBROUTINE hop_2D

!> =====================================================================
!> =====================================================================
SUBROUTINE get_downstream_2D (self, downstream)
  class (traverse_2d_t) :: self
  real                  :: r_trace
  real, dimension(2)    :: ray_step
  integer, dimension(3) :: cell_pos, downstream
! ......................................................................
  ray_step = self%ray_step                 ! save current ray steps
  cell_pos = self%cell_pos                 ! save current cell position
  r_trace  = self%r_trace                  ! save current distance
  self%uw = 0
  call hop_2D(self)                        ! make a hop
  downstream = self%cell_pos               ! that is your downstream cell
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
  
  self%cell_pos = cell_pos                 ! reset cell position
  self%ray_step = ray_step                 ! reset ray steps
  self%r_trace  = r_trace                  ! reset total distance
END SUBROUTINE get_downstream_2D

!> =====================================================================
!> =====================================================================
SUBROUTINE init_traverse_3D(self, dp, ui, li, a, b, c)
  class(traverse_3d_t) :: self
  integer, dimension(3)   :: li, ui
  real, dimension(3)   :: dp
  integer              :: a, b, c
! ......................................................................
  call trace_begin('traverse_3d_t%init')
  self%a = a; self%b = b; self%c = c          ! convenience
  self%cell_max(1) = ui(a)                    ! last cell point along a axis
  self%cell_max(2) = ui(b)                    ! last cell point along b axis
  self%cell_max(3) = ui(c)                    ! last cell point along c axis
  if (dp(a) < 0) then                         ! if ray goes from "top" to "bottom"
    self%cell_min (1) =  self%cell_max(1)     ! swap cell_max and cell_min
    self%cell_max (1) =  li(a)                ! cell_max now is first cell point
    self%cell_step(1) = -1                    ! stepper now is negative
    dp(a) = -dp(a)                            ! step size always has to be positive
  else
    self%cell_step(1) = 1                     ! step through one cell to another
    self%cell_min (1) = li(a)
  end if
  if (dp(b) < 0) then                         ! if ray goes from "top" to "bottom"
    self%cell_min (2) = self%cell_max(2)      ! swap cell_max and cell_min
    self%cell_max (2) = li(b)                 ! cell_max now is first cell point
    self%cell_step(2) = -1                    ! stepper now is negative
    dp(b)        = -dp(b)                     ! step size always has to be positive
  else
    self%cell_step(2) = 1
    self%cell_min (2) = li(b)
  end if
  self%dstep(1) = 1/dp(a)                     ! distance the ray travels to go one step
  self%dstep(2) = 1/dp(b)                     ! distance the ray travels to go one step
  self%dstep(3) = 1/dp(c)                     ! distance the ray travels to go one step
  self%cell_min(3)  = li(c); self%cell_step(3) = 1 ! ray direction is always positive
  call trace_end
END SUBROUTINE init_traverse_3D

!> =====================================================================
!> =====================================================================
SUBROUTINE posit_3D(self, ia, ib, init_step, upstream)
  class (traverse_3d_t) :: self
  real, dimension(3) :: init_step
  real, dimension(3) :: upstream
  integer               :: ia, ib
! ......................................................................
  self%cell_pos(self%a) = ia                                  ! position in a axis
  self%cell_pos(self%b) = ib                                  ! position in b axis
  self%cell_pos(self%c) = self%cell_min(3)                    ! starting position for c axis
  self%lw = 0; self%uw = 0
  if (init_step(1) == 1) then                                 ! determine whether ray starts
    self%ray_p_init(1) = ia - self%cell_step(1)/2             ! at the wall or at the centre
  else                                                        ! in terms of thjs direction
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
  
  
  self%ray_p_last = self%cell_max + self%cell_step/2           ! ray finishes at the wall
  
  if (self%cell_step(1) > 0) then                              ! estimate ray distance for a dirrection
    self%ray_dist(1) = self%ray_p_last(1) - self%ray_p_init(1)
  else
    self%ray_dist(1) = self%ray_p_init(1) - self%ray_p_last(1)
  end if
  
  if (self%cell_step(2) > 0) then                              ! estimate ray distance for b dirrection
    self%ray_dist(2) = self%ray_p_last(2) - self%ray_p_init(2)
  else
    self%ray_dist(2) = self%ray_p_init(2) - self%ray_p_last(2)
  end if
  
  self%ray_dist(3) = self%ray_p_last(3) - self%ray_p_init(3)   ! estimate ray distance for c dirrection
  self%ray_dist    = self%ray_dist * self%dstep
  self%r_tot       = minval(self%ray_dist)                     ! determine the shortest path

  
  if (self%r_tot< 0) then                                      ! never happens, but just in case
   error stop "crap"
   end if
  self%ray_step    = init_step                                 ! set initial ray steps
  self%ray_step    = self%ray_step * self%dstep                ! normalise in axis point of view
  
  if (init_step(1) == init_step(2)) then                       ! determine upstream cell position
    if  (init_step(2) == init_step(3)) then
      upstream(self%a) = self%cell_pos(self%a) - self%cell_step(1)
      upstream(self%b) = self%cell_pos(self%b) - self%cell_step(2)
      upstream(self%c) = self%cell_pos(self%c) - self%cell_step(3)
      call self%wall(1); call self%wall(2); call self%wall(3) ! and the lower face
      
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
!> Convenience. Sets the lower wall index.
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
    error stop "So something else is wrong..."
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
    self%head%n_points = self%head%n_points + 1      ! increment number of points in this ray
  else
    self%head => point
    self%tail => point
    self%head%n_points = 1
  end if
  point%coords   = p                                 ! set coordinates
  call trace_end
END SUBROUTINE add_ray_point

!> =====================================================================
!> Subroutine to transform linked-list of ray points into a stiff (array)
!> ray. Arranging is done based on the walls the ray hits.
!> INPUT: self, linked list of ray points (ray_point)
!> =====================================================================
SUBROUTINE finish_ray(self, ray_point)
  class(single_ray_list_t)           :: self
  class(single_ray_list_t), pointer  :: ray, bundle
  class(ray_point_t),       pointer  :: ray_point, prev
  integer :: pt
  logical :: entry_found
! ......................................................................
  call trace_begin('ray_geometry_t%finish_ray')
  prev => ray_point%tail
  ray_point => ray_point%head
  allocate (ray)
  ray%lwall     = ray_point%lwall
  ray%uwall     = prev%uwall
  !print *,ray%uwall, ray_point%tail%uwall
  ray%n_points  = ray_point%n_points 
  if (.not.associated(self%head)) then                ! if the ray tree is empty
    self%head             => ray
    self%head%last_bundle => ray
    self%head%first_ray   => ray
    ray%first_ray         => ray
    ray%last_ray          => ray
    ray%bundle_ID    = 1
    ray%n_rays       = 1
    ray%max_n_points = ray%n_points
  else
    entry_found = .false.
    bundle => self%head
    ! Loop over already existing bundles to find similar rays
    do while ((associated(bundle)).and.(entry_found.eqv..false.))
      if (all(bundle%lwall == ray%lwall).and.all(bundle%uwall == ray%uwall)) then
        entry_found = .true.
      else
        bundle => bundle%next_bundle
      end if
    end do
    
    if (.not.entry_found) then                        ! If appropriate branch not found - create a new branch
      ray%head => self%head
      ray%prev_bundle => self%head%last_bundle
      ray%prev_bundle%next_bundle => ray
      self%head%last_bundle => ray
      ray%first_ray    => ray
      ray%last_ray     => ray
      ray%bundle_ID     = ray%prev_bundle%bundle_ID + 1
      ray%n_rays        = 1
      ray%dir_ID        = self%dir_ID
      ray%max_n_points  = ray%n_points
    else                                              ! If appropriate branch found - create a new ray
      bundle%n_rays = bundle%n_rays + 1
      ray%bundle_ID = bundle%bundle_ID
      ray%head      => self%head
      ray%first_ray => bundle
      ray%prev_ray  => bundle%last_ray
      ray%prev_ray%next_ray => ray
      bundle%last_ray => ray
      if (bundle%max_n_points < ray%n_points) then
        bundle%max_n_points = ray%n_points
      end if
    end if
  end if
  ! Fill in all the data
  allocate(ray%cell_coords(3, ray%n_points))          ! create a stiff array for this ray with cell cords
  allocate(ray%ds (ray%n_points))                     ! create a stiff array for this ray with ds
  allocate(ray%cell_walls (3,2,ray%n_points))
  ray%fwd_coords = ray_point%fwd                      ! save upstream cell coords
  ray%rev_coords = ray_point%rev                      ! save downstream cell coords

  pt=1
  do while (associated(ray_point))                    ! loop over ray points
    ray%cell_coords(:,pt) = ray_point%coords          ! save ray point coords
    ray%ds(pt)            = ray_point%ds              ! save ray point ds
    ray%cell_walls(:,1,pt)  = ray_point%lwall
    ray%cell_walls(:,2,pt)  = ray_point%uwall
    prev => ray_point                                 ! save previous point
    ray_point => ray_point%next                       ! go to next point
    deallocate(prev)                                  ! deallocate previous point
    pt = pt+1
  end do
  call trace_end
END SUBROUTINE finish_ray

!> =====================================================================
!> Subroutine to restructure ray tree into 2D sorted ray tables.
!> =====================================================================
SUBROUTINE finish_single_ray_list(self, bundle)
  class(ray_geometry_t)             :: self
  class(single_ray_list_t), pointer :: bundle, prev, next_bundle
  class(ray_list_t),        pointer :: ray 
  integer :: iray, npt, pt
  logical :: first_bundle = .true.
! ......................................................................
  bundle => bundle%head
  do while (associated(bundle))                                 ! loop over ray bundles to be saved
    allocate (ray)                                              ! alocate new ray bundle
    ray%dir_ID       = self%ray_dirs                            ! save direction ID
    ray%n_rays       = bundle%n_rays                            ! save number of rays in this bundle 
    ray%max_n_points = bundle%max_n_points                      ! save  number of points in the longest ray                
    self%n_bundles   = self%n_bundles + 1                       ! increment ray sets
    ray%bundle_ID    = self%n_bundles                           ! get ID for this bundle
    ray%lwall        = bundle%lwall                             ! save lower wall for this bundle
    ray%uwall        = bundle%uwall                             ! save upper wall for this bundle
    if (first_bundle) then
      first_bundle = .false.                                  ! more steps for first bundle are needed
      if (associated(self%ray_list%head)) then                ! if ray bundle list is already present
        ray%prev_dir => self%ray_list%last_dir                ! add this bundle to the end of the direction list
        ray%prev_dir%next_dir => ray                          ! connect to previous direction
        self%ray_list%last_dir => ray                         ! make it the last direction
        ray%head     => self%ray_list%head                    ! connect to the head of the ray list
      else                                                    ! if ray bundle list is empty
        self%ray_list%head => ray                             ! link to self
        self%ray_list%last_dir => ray
      end if 
      ray%first_bundle => ray                                 ! this is the first bundle for this direction
      ray%last_bundle  => ray                                 ! and also the last bundle for this direction
    else
      ray%head => self%ray_list%head                          ! connect to the head of the ray list
      ray%first_bundle => self%ray_list%last_dir%first_bundle ! connect to the first bundle in this direction
      ray%prev_bundle  => self%ray_list%last_dir%last_bundle  ! connect to previous bundle in this direction
      ray%prev_bundle%next_bundle => ray                      ! connect previous bundle to this ray
      self%ray_list%last_dir%last_bundle => ray               ! make it the last bundle in this direction
    end if 
  ! Save the data
    allocate (ray%n_points(ray%n_rays))                               ! number of points in all rays in the bundle
    allocate (ray%cell_coords(3,ray%max_n_points,ray%n_rays))         ! cell coords per point per ray in the bundle
    allocate (ray%ds(ray%max_n_points,ray%n_rays))                    ! ds per point per ray in the bundle
    allocate (ray%fwd_coords(3,ray%n_rays))                           ! upstream cell coords per ray
    allocate (ray%rev_coords(3,ray%n_rays))                           ! downstream cell coords per ray
    allocate (ray%cell_domain(self%ui(1),self%ui(2),self%ui(3)))
    allocate (ray%pt_wall(3,2,ray%max_n_points,ray%n_rays))
    next_bundle => bundle%next_bundle
    do iray = 1, ray%n_rays                                           ! loop over all rays in the tree
      npt = bundle%n_points                                           ! n_points of this ray
      ray%cell_coords(:,1:npt,iray) = bundle%cell_coords(:,1:npt)     ! save cell cords
      ray%ds(1:npt,iray)            = bundle%ds(:)                    ! save distances, traveled in those cells
      ray%n_points(iray)            = bundle%n_points
      ray%pt_wall(:,1,1:npt, iray) = bundle%cell_walls(:,1,1:npt)
      ray%pt_wall(:,2,1:npt, iray) = bundle%cell_walls(:,2,1:npt)
      do pt = 1, ray%n_points(iray)
        call ray%store_cell_data(ray%cell_coords(:,pt,iray),pt,iray,&
        ray%pt_wall(:,1,pt,iray),ray%pt_wall(:,2,pt,iray))
      end do
      ray%fwd_coords(:,iray)        = bundle%fwd_coords               ! save upstream coords
      ray%rev_coords(:,iray)        = bundle%rev_coords               ! save downstream coords
      prev   => bundle                                                ! hold pointer to this ray
      bundle => bundle%next_ray                                       ! go to next ray
      deallocate (prev)                                               ! deallocate previous ray
    end do
    bundle => next_bundle;                          ! move to next bundle and deallocate the previous one
  end do
END SUBROUTINE finish_single_ray_list

SUBROUTINE store_cell_data(self,c,pt,ray,lw,uw)
  class (ray_t)     :: self
  type (cell_t), pointer :: cell
  integer :: pt, ray, lw(3), uw(3),c(3)
!-----------------------------------------------------------------------
  allocate (cell)
  cell%ray   = ray
  cell%point = pt
  cell%lwall = lw
  cell%uwall = uw
  if (.not.associated(self%cell_domain(c(1),c(2),c(3))%head)) then
    cell%head => cell
    cell%tail => cell
    self%cell_domain(c(1),c(2),c(3))%head => cell
    self%cell_domain(c(1),c(2),c(3))%tail => cell
  else
    cell%head => self%cell_domain(c(1),c(2),c(3))%head
    self%cell_domain(c(1),c(2),c(3))%tail%next => cell
    self%cell_domain(c(1),c(2),c(3))%tail => cell
  end if
  
END SUBROUTINE store_cell_data


!> =====================================================================
!> Sort through the ray geometry list and select appropriate ray bundles
!> for particular direction
!> =====================================================================
SUBROUTINE get_my_rays (self, bundle, bundle_ID, inv)
  class(ray_geometry_t) :: self
  class(ray_list_t), pointer :: list_dir, list_bundle, bundle
  integer :: bundle_ID
  logical :: bundle_found
  logical :: inv
!-----------------------------------------------------------------------
  bundle_found = .false.
  list_dir => self%ray_list%head
  list_bundle => list_dir
  do while (associated(list_dir).and.(bundle_found.eqv..false.))
    do while (associated(list_bundle).and.(bundle_found.eqv..false.))
      if (bundle_ID == list_bundle%bundle_id) then
        bundle_found = .true.
        bundle => list_bundle
      else
        list_bundle => list_bundle%next_bundle
      end if
    end do
    if (bundle_found.eqv..true.) then 
      if (inv) then
        list_bundle => bundle
        call list_bundle%invert_my_rays(bundle)
      end if
      return
    else
      list_dir => list_dir%next_dir
      list_bundle => list_dir
    end if
  end do
  if (bundle_found.eqv..false.) then
  error stop "Appropriate direction not found. Very strange."
  end if
END SUBROUTINE get_my_rays

SUBROUTINE invert_my_rays (self, inv)
  class (ray_list_t) :: self
  class (ray_list_t), pointer :: inv
  integer :: ray, pt, ipt
! ----------------------------------------------------------------------

  allocate(inv)
  
  inv%n_rays = self%n_rays
  inv%lwall = self%uwall
  inv%uwall = self%lwall
  inv%max_n_points = self%max_n_points
  
  allocate(inv%n_points(inv%n_rays))
  allocate(inv%cell_coords(size(self%cell_coords,1)&
                          ,size(self%cell_coords,2)&
                          ,size(self%cell_coords,3)))
  allocate(inv%fwd_coords(3, inv%n_rays))
  allocate(inv%rev_coords(3, inv%n_rays))
  allocate(inv%ds(size(self%ds,1),size(self%ds,2)))
  allocate(inv%pt_wall(3,2,size(self%pt_wall,3),size(self%pt_wall,4)))
  allocate(inv%cell_domain(size(self%cell_domain,1),&
           size(self%cell_domain,2),size(self%cell_domain,3)))
  do ray = 1, self%n_rays
    inv%n_points(ray) = self%n_points(ray)
    inv%fwd_coords(:,ray) = self%rev_coords(:,ray)
    inv%rev_coords(:,ray) = self%fwd_coords(:,ray)
    do pt = 1, self%n_points(ray)
      ipt = self%n_points(ray) - pt + 1
      inv%ds(ipt,ray) = self%ds(pt,ray)
      inv%cell_coords(:,ipt, ray) = self%cell_coords(:,pt,ray)
      inv%pt_wall(:,1,ipt,ray) = self%pt_wall(:,2,pt,ray)
      inv%pt_wall(:,2,ipt,ray) = self%pt_wall(:,1,pt,ray)
      call inv%store_cell_data(inv%cell_coords(:,ipt, ray), ipt, ray, &
      inv%pt_wall(:,1,ipt, ray), inv%pt_wall(:,2,ipt, ray))
    end do
  end do

END SUBROUTINE invert_my_rays

SUBROUTINE test_ray_geometry(self, dir, bundle_ID, ray)
  class(ray_geometry_t) :: self
  class(ray_list_t), pointer :: list_dir, list_bundle, bundle
  integer :: dir, bundle_ID, ray, pt, p(3)
  
  call get_my_rays(self,bundle, bundle_ID, .false.)
  
  print *, 'Testing ray_geom for dir_ID, bundle_ID, ray:',bundle%dir_ID, bundle_ID, ray
  print *, 'n_rays in bundle:',bundle%n_rays
  print *, 'lwall, uwall', bundle%lwall, bundle%uwall
  print *, 'n_points in test ray:',bundle%n_points(ray)
  print *, 'cell coords in ray domain , ray info from cell domain'
  
  do pt=1, bundle%n_points(ray)
    p = bundle%cell_coords(:,pt,ray)
    print *, 'ray', ray,'pt', pt, bundle%cell_coords(:,pt,ray),&
    'is at ray,pt',bundle%cell_domain(p(1),p(2),p(3))%head%ray,&
     bundle%cell_domain(p(1),p(2),p(3))%head%point
    print *,'ray domain lwall, uwall', bundle%pt_wall(:,1,pt,ray),bundle%pt_wall(:,2,pt,ray)
    print *,'cell domain lwall, uwall',bundle%cell_domain(p(1),p(2),p(3))%head%lwall,&
    bundle%cell_domain(p(1),p(2),p(3))%head%uwall
  end do
  
END SUBROUTINE test_ray_geometry

END MODULE ray_geometry_mod
