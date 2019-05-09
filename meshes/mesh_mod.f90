!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for mesh
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mesh_mod
  USE trace_mod
  USE io_mod
  USE mpi_mod
  USE kinds_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! A generic orthogonal, curvilinear mesh type.
  !
  ! Rather than have one type for each cardinal direction, all metric factors
  ! for all directions are defined here, but only those for the current axis
  ! are allocated upon initialisation. This way, you can have one type for meshes
  ! in any direction.
  !
  ! Primary scale factors: h2c, h2f, h31c, h31f, h32c, h32f
  !   `c` means the quantity is located at cell-centre; `f` means the quantity
  !   is located at the cell-face. Regardless of whether you are using a
  !   staggered method or not, you will generally need these factors if using
  !   curvilinear coords.
  ! Secondary volume and area factors: dvol1c, dvol1f, dvol2c, dvol2f, dvol3c,
  !                                    dvol3f, dar1c, dar1f, dar2c, dar2f,
  !                                    dar31c, dar31f, dar32c, dar32f
  !
  ! Current requirements for curvilinear meshes:
  ! - The scale factors are either a function of one variable or are separable
  !   (e.g., h3 = r sin \theta in spherical coords. but can be expressed as
  !   h31(i) = r(i) and h32(j) = sin theta(j)).
  ! - The mesh spacing (self%d) is assumed to be uniform in all directions.
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  type, abstract, public :: mesh_t
    integer :: id=0          ! copy of task id, for info
    integer :: n = 16        ! number of "active" zones
    integer :: ng = 2        ! number of ghost zones
    integer :: nc            ! number of cells
    integer :: lb            ! lower boundary of mesh incl. ghost zones
    integer :: li            ! lower boundary of mesh excl. ghost zones
    integer :: lo
    integer :: ub            ! upper boundary of mesh incl. ghost zones
    integer :: ui            ! upper boundary of mesh excl. ghost zones
    integer :: uo
    integer :: gn = 20  ! total number of zones incl. ghost zones (= n + 2 * ng)
    integer :: idir          ! identify which direction the mesh applies to
    real(8) :: p = 0d0       ! position of patch centre in **Cartesian** coords.
    real(8) :: s = 0d0       ! spatial extent
    real(8) :: d = 0d0       ! grid spacing
    real(8) :: b = 1d0       ! box size
    real(8) :: o = 0.0       ! position offset in index space
    real(4) :: lf = 0.0, uf = 0.0  ! interpolation boundary indices (must be floating point)
    ! the variables `llc_cart` and `llc_nat` describe the same location, but in
    ! different reference frames.
    real(8) :: llc_cart = 0.0  ! position of the patch's lower left corner in Cartesian
                               ! coords.; replacement for `p`.
    real(8) :: llc_nat = 0.0 ! position of the patch's lower left corner in relative,
                             ! native coords.; always = 0 in Cartesian coords.
    real(8) :: centre_nat    ! volume centroid of the mesh in its *own* coords.
    real(8), dimension(:), pointer :: r => null()    ! coordinate relative to patch centre
    real(8), dimension(:), pointer :: h => null()    ! offset for staggered variables;
                                                     ! should be initialised in `mhd_t%init`
    real, dimension(:), pointer :: ddn   => null()     ! borrow for now (FIXME) 
    real, dimension(:), pointer :: dup   => null()     ! borrow for now (FIXME)
    real, dimension(:), pointer :: dsup  => null()     ! borrow for now (FIXME)
    real, dimension(:), pointer :: dsdn  => null()     ! borrow for now (FIXME)  
    real, dimension(:), pointer :: b_m   => null()     ! borrow for now (FIXME) -- bifrost xm/ym/zm
    real, dimension(:), pointer :: b_mdn => null()     ! borrow for now (FIXME) -- bifrost xmdn/ymdn/zmdn
    real, dimension(:), pointer :: b_1d  => null()     ! borrow for now (FIXME) -- bifrost dx1d/dy1d/dz1d
    real, dimension(:), pointer :: b_didup => null()   ! borrow for now (FIXME) -- bifrost dxidxup/...
    real, dimension(:), pointer :: b_diddn => null()   ! borrow for now (FIXME) -- bifrost dxidxdn/...
    real, dimension(:,:), pointer::  zlincoeff => null()! borrow for now (FIXME) -- bifrost zlincoeff
    real, dimension(:,:), pointer:: dzlincoeff => null()! borrow for now (FIXME) -- bifrost zlincoeff
    real, dimension(:,:), pointer::  zuincoeff => null()! borrow for now (FIXME) -- bifrost zlincoeff
    real, dimension(:,:), pointer:: dzuincoeff => null()! borrow for now (FIXME) -- bifrost zlincoeff
    real(8)   :: b_d  ! borrow for now (FIXME)  -- grid spacing in bifrost
    real      :: b_ad, b_bd, b_cd                      !borrow for now (FIXME)  -- bifrost cdx/bdx/adx...
    integer(4):: b_s  ! borrow for now (FIXME)  -- same as li, but in bifrost
    integer(4):: b_sb ! borrow for now (FIXME)  -- same as lb, but in bifrost
    integer(4):: b_e  ! borrow for now (FIXME)  -- same as ui, but in bifrost
    integer(4):: b_eb ! borrow for now (FIXME)  -- sane as ub, but in bifrost
    integer(kind=4) :: b_r, b_rb
    
    character(len=16):: type = ''
    real(8), allocatable :: zeta  ! position angle of the phi = 0 axis relative to the x-axis (cyl. and sph. coords.).
    real(8), allocatable :: xi    ! position angle of the theta = 0 relative to the neg. z-axis (sph. coords.).
    real(8):: origin=0d0
    !---------------------------------------------------------------------------
    ! Physical boundary condition info; indices beyond which BCs rule
    !---------------------------------------------------------------------------
    logical:: lower_boundary = .false.
    logical:: upper_boundary = .false.
    logical:: no_mans_land = .false.

    procedure (), pointer, nopass :: CurrentToCartesian => null()
    ! geometric factors as procedure pointers
    procedure(func_interface), pointer :: h2c => null()
    procedure(func_interface), pointer :: h2f => null()
    procedure(func_interface), pointer :: h31c => null()
    procedure(func_interface), pointer :: h31f => null()
    procedure(func_interface), pointer :: h32f => null()
    procedure(func_interface), pointer :: h32c => null()
    procedure(func_interface), pointer :: dvol1c => null()
    procedure(func_interface), pointer :: dvol1f => null()
    procedure(func_interface), pointer :: dvol2c => null()
    procedure(func_interface), pointer :: dvol2f => null()
    procedure(func_interface), pointer :: dvol3c => null()
    procedure(func_interface), pointer :: dvol3f => null()
    procedure(func_interface), pointer :: dar1c => null()
    procedure(func_interface), pointer :: dar1f => null()
    procedure(func_interface), pointer :: dar2c => null()
    procedure(func_interface), pointer :: dar2f => null()
    procedure(func_interface), pointer :: dar31c => null()
    procedure(func_interface), pointer :: dar31f => null()
    procedure(func_interface), pointer :: dar32c => null()
    procedure(func_interface), pointer :: dar32f => null()
  contains
    procedure :: init       => init_mesh
    procedure :: deallocate => deallocate_coord
    procedure :: print      => print_mesh
    procedure(init_geometryInterface), deferred :: init_geometry
  end type
  abstract interface
    subroutine init_geometryInterface(self, idir, origin)
      import mesh_t
      class(mesh_t) :: self
      integer, intent(in) :: idir
      real(8), intent(in) :: origin
    end subroutine

    function func_interface(self, i) result(x)
      import mesh_t
      class(mesh_t) :: self
      integer, intent(in) :: i
      real(8) :: x
    end function
  end interface

  !-----------------------------------------------------------------------------
  ! Utility list of coordinate systems available for meshes
  !-----------------------------------------------------------------------------
  type :: mesh_type
    integer :: Cartesian = 1
    integer :: spherical = 2
    integer :: cylindrical = 3
  end type
  type(mesh_type), public :: mesh_types

  !-----------------------------------------------------------------------------
  ! Cartesian: x,y,z
  !-----------------------------------------------------------------------------
  type, extends(mesh_t), public :: Cartesian_mesh
  contains
    procedure :: init_geometry => init_Cartesian_geometry
  end type Cartesian_mesh

  !-----------------------------------------------------------------------------
  ! cylindrical: z, R, Phi
  ! Assert that, when Phi is 0, the r-coord. and x-coord. overlie one another.
  ! The z-coords. should always overlap.
  ! Further, the angle `zeta` specifies the offset between the Phi = 0 and x axes.
  !-----------------------------------------------------------------------------
  type, extends(mesh_t), public :: cylindrical_mesh
  contains
    procedure :: init_geometry => init_cylindrical_geometry
  end type cylindrical_mesh

  !-----------------------------------------------------------------------------
  ! spherical: r, theta, phi
  ! Assert that, when theta = 0, the r-coord. and the *negative* z-coord overlie
  ! one another.
  ! Further, assert that, when phi = 0 and theta = pi/2, the r-coord. and x-coord.
  ! overlie one another.
  ! The angle `zeta` specifies the offset between the phi = 0 and x axes, while
  ! the angle `xi` specified the offset between the theta=0 and neg. z axes.
  !-----------------------------------------------------------------------------
  type, extends(mesh_t), public :: spherical_mesh
  contains
    procedure :: init_geometry => init_spherical_geometry
  end type spherical_mesh

  interface CylindricalToCartesian
    module procedure :: CylindricalToCartesian4, CylindricalToCartesian8
  end interface
  interface CylindricalToSpherical
    module procedure :: CylindricalToSpherical4, CylindricalToSpherical8
  end interface
  interface SphericalToCartesian
    module procedure :: SphericalToCartesian4, SphericalToCartesian8
  end interface
  interface CartesianToSpherical
    module procedure :: CartesianToSpherical4, CartesianToSpherical8
  end interface
  interface CartesianToCylindrical
    module procedure :: CartesianToCylindrical4, CartesianToCylindrical8
  end interface
  interface SphericalToCylindrical
    module procedure :: SphericalToCylindrical4, SphericalToCylindrical8
  end interface
  interface CurrentToCartesian
    module procedure :: CurrentToCartesian4, CurrentToCartesian8
  end interface
  interface CartesianToCurrent
    module procedure :: CartesianToCurrent4, CartesianToCurrent8
  end interface
  interface round
    module procedure :: round4, round8
  end interface

  real(8), parameter :: pi = 3.14159265358979323846

PUBLIC MeshFactory, MeshRecycler, &
       CartesianToCylindrical, CartesianToSpherical, CylindricalToCartesian, &
       CylindricalToSpherical, SphericalToCartesian, SphericalToCylindrical, &
       CurrentToCartesian, CartesianToCurrent, round
CONTAINS

!===============================================================================
!> Create a mesh collection.
!> This method creates three orthogonal meshes to form an 3-dimensional grid.
!> It is possible for up to 2 of the meshes to have extent one.
!===============================================================================
SUBROUTINE MeshFactory (meshtype, collection, n, ng, s, p, llc_native, nml)
  class(mesh_t), dimension(:), pointer :: collection
  integer, intent(in) :: meshtype
  integer, dimension(3), optional :: ng
  integer, dimension(3) :: ngi
  integer, dimension(3), optional :: n
  integer, dimension(3) :: ni
  real(8), dimension(3), optional :: s, p
  real(8), dimension(3) :: si, pi
  logical, optional :: nml
  integer :: i
  ! expected form of `origin` for Cartesian coords.: [0,0,0] (unused)
  !                           for cylindrical coords.: [rmin, 0, zeta]
  !                           for spherical coords.: [rmin, xi, zeta]
  real(8), dimension(3), optional :: llc_native
  real(8), dimension(3) :: llc_native_i
  !.............................................................................
  ! set reasonable defaults
  if (present(n)) then
    ni = n
  else
    ni(:) = 16
  end if
  if (present(s)) then
    si = s
  else
    si(:) = 1.0
  end if
  if (present(p)) then
    pi = p
  else
    pi(:) = 0.0
  end if
  if (present(llc_native)) then
    llc_native_i = llc_native
  else
    llc_native_i(:) = 0.0
  end if
  if (present(ng)) then
    ngi = ng
  else
    ngi(:) = 2
  end if


  if (associated(collection)) call MeshRecycler(collection)
  ! Allocate based on the desired mesh type.
  ! The grid is always assumed to be 3-D, even if one of the directions is
  ! collapsed/a symmetry direction (in which case, n = 1).
  if (meshtype .eq. mesh_types%Cartesian) then
    allocate(Cartesian_mesh :: collection(3))
    collection%type = 'Cartesian_mesh'
  else if ( meshtype .eq. mesh_types%spherical ) then
    allocate(spherical_mesh :: collection(3))
    collection%type = 'spherical_mesh'
  else if ( meshtype .eq. mesh_types%cylindrical ) then
    allocate(cylindrical_mesh :: collection(3))
    collection%type = 'cylindrical_mesh'
  else
    call mpi%abort("Invalid mesh type. Abort!")
  end if
  if (present(nml)) then
    collection%no_mans_land = nml
  end if

  do i=1,3
    ! initialise the mesh
    collection(i)%idir = i
    call collection(i)%init(ni(i),ngi(i),si(i),pi(i),llc_native_i(i))

    ! initialise geometric scale factors (if applicable).
    call collection(i)%init_geometry(i,llc_native_i(i))
  end do

END SUBROUTINE MeshFactory

!===============================================================================
!> A means to entirely deallocate a collection of meshes
!===============================================================================
SUBROUTINE MeshRecycler (collection)
  class(mesh_t), dimension(:), pointer :: collection
  integer :: i
  !.............................................................................
  do i=1,size(collection)
    call collection(i)%deallocate
  end do
  deallocate(collection)
  nullify (collection)
END SUBROUTINE MeshRecycler

!===============================================================================
!> The integer mesh is defined by specifying the number of points `n` and the
!> number of guard cells `ng`.  The real mesh is then fully specified by giving
!> the size `s`. The center of the mesh, `p` can also be supplied.
!===============================================================================
SUBROUTINE init_mesh (self, n, ng, s, p, llc_native)
  class(mesh_t) :: self
  integer :: i
  integer, optional :: n, ng
  real(8), optional :: s, p
  real(8), optional :: llc_native
  real(8):: h
  !.............................................................................
  call trace_begin ('mesh%init')

  if (present(n)) then
    self%n = n
  endif
  if (present(ng)) then
    self%ng = ng
  endif
  if (present(s)) then
    self%s = s
  endif
  if (present(p)) then
    self%p = p
  endif
  if (present(llc_native)) then
    self%llc_nat = llc_native  ! the lower left corner in the current patch coords.
  endif

  ! initialise mesh indices
  self%lb = 1
  self%li = self%lb + self%ng
  self%lo = self%li - 1
  self%ui = self%lo + self%n
  self%uo = self%ui + 1
  self%ub = self%ui + self%ng
  ! account for collapsed dimensions
  if (self%n == 1) then
    self%lo = self%li
    self%uo = self%ui
    self%ub = self%ui
  end if
  self%gn = self%ub
  allocate (self%r(self%gn))
  ! initialise mesh floating point attributes
  self%o = (self%li+self%ui)*0.5

  ! calculate the volume centroid of the current mesh
  self%centre_nat = VolumeCentroid(self)
  ! cell-centred coordinate relative to the patch *centre*; the absolute position
  ! of a cell can be recovered using `self%centre_nat + self%r(i)` in the coord. system
  ! of the patch, and `self%p + self%centre_nat + self%r(i)`.
  ! in Cartesian coords., self%centre_nat = self%p
  if (self%no_mans_land) then
    h = 0.5_8
    if (self%n == 1) h = 1.0_8 ! account for collapsed dimensions (otherwise `r` ends up non-zero).
    self%nc = self%ui - self%li + 1
  else
    h = 1.0_8
    self%nc = max(1, self%ui - self%li)
  end if
  self%d = self%s/max(1,self%nc)
  do i=1,self%gn
    self%r(i) = real(i - h - self%ng,kind=8) * self%d + self%llc_nat - self%centre_nat
  end do
  self%lf = -self%s*0.5_8
  self%uf = +self%s*0.5_8
  call trace_end
END SUBROUTINE init_mesh

!===============================================================================
!> deallocates the position coordinate of a mesh
!===============================================================================
SUBROUTINE deallocate_coord (self)
  class(mesh_t) :: self
  !.............................................................................
  if (associated(self%r))   deallocate (self%r)
  if (associated(self%h))   deallocate (self%h)
  if (associated(self%ddn)) deallocate (self%ddn)
  if (associated(self%dup)) deallocate (self%dup)
END SUBROUTINE deallocate_coord

!===============================================================================
!> print mesh information to the stdout
!===============================================================================
SUBROUTINE print_mesh (self)
  class(mesh_t):: self
  integer:: i
  !.............................................................................

  if (io%verbose <= 0) return
  print*,'::mesh%print::'
  select type(self)
  type is(Cartesian_mesh)
    print*,'Cartesian', self%idir
  type is(spherical_mesh)
    print*,'spherical', self%idir
  type is(cylindrical_mesh)
    print*,'cylindrical', self%idir
  end select
  print'(a,3i8)', "n, ng, gn: ", self%n, self%ng, self%gn
  print'(a,6i8)', "lb, lo, li, ub, uo, ui: ", self%lb, self%lo, self%li, self%ub, self%uo, self%ui
  print'(a,3x,3g16.8)', "lf, uf, p: ",self%lf, self%uf, self%p
  print'(a,3x,3g16.8)', "s, d, b: ", self%s, self%d, self%b
  print'(a,3x,4g16.8)', "llc_nat, centre_nat, llc_cart, o: ",self%llc_nat, self%centre_nat, self%llc_cart, self%o
  do i=1,self%gn
    print'(i4,f14.8)',i,self%r(i)
  end do

  if (allocated(self%zeta)) print*,"phi axis offset: ", self%zeta
  if (allocated(self%xi)) print*,"theta axis offset: ", self%xi

  print*,'----------'
END SUBROUTINE print_mesh

!===============================================================================
!> Dummy initialise for a Cartesian mesh: All scale factors are one and all
!> volume elements are equal to dx.
!===============================================================================
SUBROUTINE init_Cartesian_geometry (self, idir, origin)
  class(Cartesian_mesh) :: self
  integer, intent(in) :: idir
  real(8), intent(in) :: origin
  integer :: n
  !.............................................................................
  ! All the metric factors are 1 in Cartesian coords. and all the volume coords.
  ! are simply the cell sizes.

  select case (idir)
  case (1)
    self%h2c => func_one
    self%h2f => func_one
    self%h31c => func_one
    self%h31f => func_one
    self%dvol1c => func_dx
    self%dvol1f => func_dx

    self%dar1c => func_one
    self%dar1f => func_one
    self%dar2c => func_one
    self%dar2f => func_one
    self%dar31c => func_one
    self%dar31f => func_one
  case (2)
    self%h32c => func_one
    self%h32f => func_one
    self%dvol2c => func_dx
    self%dvol2f => func_dx

    self%dar32c => func_one
    self%dar32f => func_one
  case (3)
    self%dvol3c => func_dx
    self%dvol3f => func_dx
  end select

END SUBROUTINE init_Cartesian_geometry

!===============================================================================
!> Initalise the geometric factors for a cylindrical mesh.
!> A cylindrical mesh is defined as having the form (z,R,\phi).
!===============================================================================
SUBROUTINE init_cylindrical_geometry (self, idir, origin)
  class(cylindrical_mesh) :: self
  integer, intent(in) :: idir
  real(8), intent(in) :: origin
  !.............................................................................
  select case (idir)
  case (1)
    self%h2c => func_one
    self%h2f => func_one
    self%h31c => func_one
    self%h31f => func_one
    self%dvol1c => func_dx
    self%dvol1f => func_dx

    self%dar1c => func_one
    self%dar1f => func_one
    self%dar2c => func_one
    self%dar2f => func_one
    self%dar31c => func_one
    self%dar31f => func_one
  case (2)
    self%h32c => func_radius_c
    self%h32f => func_radius_f
    self%dvol2c => func_dr2_c
    self%dvol2f => func_dr2_f

    self%dar32c => func_dar32_c
    self%dar32f => func_dar32_f
  case (3)
    self%dvol3c => func_dx
    self%dvol3f => func_dx

    ! when phi = 0 and the radial vector overlies the x-axis, then zeta == 0.
    allocate(self%zeta)
    self%zeta = origin
  end select

END SUBROUTINE init_cylindrical_geometry

!===============================================================================
!> Initalise the geometry factors for a spherical mesh.
!> A spherical mesh is defined as having the form (r,\theta,\phi).
!===============================================================================
SUBROUTINE init_spherical_geometry (self, idir, origin)
  class(spherical_mesh) :: self
  integer, intent(in) :: idir
  real(8), intent(in) :: origin
  !.............................................................................

  select case (idir)
  case (1)
    self%h2c => func_radius_c
    self%h2f => func_radius_f
    self%h31c => func_radius_c
    self%h31f => func_radius_f
    self%dvol1c => func_dr3_c
    self%dvol1f => func_dr3_f

    self%dar1c => func_dar1_c
    self%dar1f => func_dar1_f
    self%dar2c => func_dar2_c
    self%dar2f => func_dar2_f
    self%dar31c => func_dar31_c
    self%dar31f => func_dar31_f
  case (2)
    self%h32c => func_sintheta_c
    self%h32f => func_sintheta_f
    self%dvol2c => func_dcostheta_c
    self%dvol2f => func_dcostheta_f

    self%dar32c => func_dar32_c
    self%dar32f => func_dar32_f

    ! when phi = theta = 0 and the radial vector overlies the x-axis, then xi == 0.
    allocate(self%xi)
    self%xi = origin
  case (3)
    self%dvol3c => func_dx
    self%dvol3f => func_dx

    ! when phi = theta = 0 and the radial vector overlies the x-axis, then zeta == 0.
    allocate(self%zeta)
    self%zeta = origin
  end select

END SUBROUTINE init_spherical_geometry

!===============================================================================
!> Collection of functions used to calculate geometric factors.
!===============================================================================
FUNCTION func_one (self, i) RESULT(one)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: one
  !.............................................................................
  one = 1.0

END FUNCTION func_one

FUNCTION func_dx (self, i) RESULT(dx)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dx
  !.............................................................................
  dx = self%d

END FUNCTION func_dx

FUNCTION func_position_c (self, i) RESULT(r)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: r
  !.............................................................................
  r = self%centre_nat + self%r(i)

END FUNCTION func_position_c

FUNCTION func_position_f (self, i) RESULT(r)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: r
  !.............................................................................
  if (i /= self%gn .or. self%gn == 1) then
    r = self%centre_nat + self%r(i) - 0.5 * self%d
  else
    r = self%centre_nat + self%r(i-1) + 0.5 * self%d
  end if

END FUNCTION func_position_f

FUNCTION func_radius_c (self, i) RESULT(r)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: r
  !.............................................................................
  r = abs(func_position_c(self, i))

END FUNCTION func_radius_c

FUNCTION func_radius_f (self, i) RESULT(r)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: r
  !.............................................................................
  r = abs(func_position_f(self, i))

END FUNCTION func_radius_f

FUNCTION func_sintheta_c (self, j) RESULT(sintc)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: sintc
  !.............................................................................
  sintc = abs(sin(func_position_c(self,j)))

END FUNCTION func_sintheta_c

FUNCTION func_sintheta_f (self, j) RESULT(sintf)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: sintf
  !.............................................................................
  sintf = abs(sin(func_position_f(self,j)))

END FUNCTION func_sintheta_f

!===============================================================================
!> Differential volume elements for non-Cartesian grids.
!****Note** that the 'f' and 'c' subscripts still denote the centring of a given
!> quantity, rather than the centring of the quantities that go into the
!> calculation. This is in contrast to ZEUS where the subscript ('a' or 'b') is
!> is determined by the centring of the inputs.
!===============================================================================
FUNCTION func_dr2_c (self, j) RESULT(dvol2c)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: dvol2c
  !.............................................................................
  dvol2c = 0.5 * abs( self%h32f(j+1) * func_position_f(self,j+1) &
                    - self%h32f(j  ) * func_position_f(self,j  ))

END FUNCTION func_dr2_c

FUNCTION func_dr2_f (self, j) RESULT(dvol2f)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: dvol2f
  !.............................................................................
  dvol2f = 0.5 * abs( self%h32c(j  ) * func_position_c(self,j  ) &
                    - self%h32c(j-1) * func_position_c(self,j-1))

END FUNCTION func_dr2_f

FUNCTION func_dr3_c (self, i) RESULT(dvol1c)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dvol1c
  !.............................................................................
  dvol1c = self%h2f(i+1) * self%h31f(i+1) * func_position_f(self,i+1) / 3.0d0 &
         - self%h2f(i  ) * self%h31f(i  ) * func_position_f(self,i  ) / 3.0d0

END FUNCTION func_dr3_c

FUNCTION func_dr3_f (self, i) RESULT(dvol1f)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dvol1f
  !.............................................................................
  dvol1f = self%h2c(i  ) * self%h31c(i  ) * func_position_c(self,i  ) / 3.0d0 &
         - self%h2c(i-1) * self%h31c(i-1) * func_position_c(self,i-1) / 3.0d0

END FUNCTION func_dr3_f

FUNCTION func_dcostheta_c (self, j) RESULT(dvol2c)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: dvol2c
  !.............................................................................
  dvol2c = -cos(func_position_f(self,j+1)) + cos(func_position_f(self,j))

END FUNCTION func_dcostheta_c

FUNCTION func_dcostheta_f (self, j) RESULT(dvol2f)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: dvol2f
  !.............................................................................
  dvol2f = -cos(func_position_c(self,j)) + cos(func_position_c(self,j-1))
END FUNCTION func_dcostheta_f

!===============================================================================
!> Secondary area factors; primarily used with operators.
!===============================================================================
FUNCTION func_dar1_c (self, i) RESULT(dar1c)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dar1c
  !.............................................................................
  dar1c = self%h2c(i) * self%h31c(i)

END FUNCTION func_dar1_c

FUNCTION func_dar1_f (self, i) RESULT(dar1f)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dar1f
  !.............................................................................
  dar1f = self%h2f(i) * self%h31f(i)

END FUNCTION func_dar1_f

FUNCTION func_dar2_c (self,i) RESULT(dar2c)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dar2c
  !.............................................................................
  dar2c = self%h31c(i) * self%d / self%dvol1c(i)

END FUNCTION func_dar2_c

FUNCTION func_dar2_f (self,i) RESULT(dar2f)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dar2f
  !.............................................................................
  dar2f = self%h31f(i) * self%d / self%dvol1f(i)

END FUNCTION func_dar2_f

FUNCTION func_dar31_c (self,i) RESULT(dar31c)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dar31c
  !.............................................................................
  dar31c = self%h2c(i) * self%d / self%dvol1c(i)

END FUNCTION func_dar31_c

FUNCTION func_dar31_f (self,i) RESULT(dar31f)
  class(mesh_t) :: self
  integer, intent(in) :: i
  real(8) :: dar31f
  !.............................................................................
  dar31f = self%h2f(i) * self%d / self%dvol1f(i)

END FUNCTION func_dar31_f

FUNCTION func_dar32_c (self, j) RESULT(dar32c)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: dar32c
  !.............................................................................
  dar32c = self%d / self%dvol2c(j)

END FUNCTION func_dar32_c

FUNCTION func_dar32_f (self, j) RESULT(dar32f)
  class(mesh_t) :: self
  integer, intent(in) :: j
  real(8) :: dar32f
  !.............................................................................
  dar32f = self%d / self%dvol2f(j)

END FUNCTION func_dar32_f

!===============================================================================
!> Mapping/transformation procedures.
!> `OCart` is the origin of the *target* coord. system in *Cartesian* coords.
!===============================================================================
FUNCTION CartesianToCylindrical4 (x, y, z, OCart) result (vec)
  real(4), intent(in) :: x, y, z
  real(4), dimension(3), optional:: OCart
  real(4), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0
  if (present(OCart)) lOCart = OCart

  vec(1) = z - lOCart(3)                              ! Z
  vec(2) = sqrt((x-lOCart(1))**2 + (y-lOCart(2))**2)  ! R
  if (lOCart(1) == 0.0 .and. lOCart(2) == 0.0) then
    vec(3) = atan2(y, x)                              ! Phi
  else
    vec(3) = atan2(y-lOCart(2), x-lOCart(1))          ! Phi
  end if
  if (vec(3) .lt. 0.0) vec(3) = vec(3) + 2.0 * real(pi,kind=4)

END FUNCTION CartesianToCylindrical4

FUNCTION CartesianToCylindrical8 (x, y, z, OCart) result (vec)
  real(8), intent(in) :: x, y, z
  real(8), dimension(3), optional:: OCart
  real(8), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0d0
  if (present(OCart)) lOCart = OCart

  vec(1) = z - lOCart(3)                              ! Z
  vec(2) = sqrt((x-lOCart(1))**2 + (y-lOCart(2))**2)  ! R
  if (lOCart(1) == 0.0 .and. lOCart(2) == 0.0) then
    vec(3) = atan2(y, x)                              ! Phi
  else
    vec(3) = atan2(y-lOCart(2), x-lOCart(1))          ! Phi
  end if
  if (vec(3) .lt. 0.0d0) vec(3) = vec(3) + 2.0d0 * pi

END FUNCTION CartesianToCylindrical8

FUNCTION CartesianToSpherical4 (x, y, z, OCart) result (vec)
  real(4), intent(in) :: x, y, z
  real(4), dimension(3), optional:: OCart
  real(4), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0
  if (present(OCart)) lOCart = OCart

  vec(1) = sqrt((x-lOCart(1))**2 + (y-lOCart(2))**2 + (z-lOCart(3))**2)  ! r
  vec(2) = pi - acos((z-lOCart(3)) / vec(1))                             ! theta
  if (lOCart(1) == 0.0 .and. lOCart(2) == 0.0) then
    vec(3) = -atan2(y, x)                                                ! phi
  else
    vec(3) = -atan2((y-lOCart(2)), (x-lOCart(1)))                        ! phi
  end if
  if (vec(3) .lt. 0.0) vec(3) = vec(3) + 2.0 * real(pi,kind=4)

END FUNCTION CartesianToSpherical4

FUNCTION CartesianToSpherical8 (x, y, z, OCart) result (vec)
  real(8), intent(in) :: x, y, z
  real(8), dimension(3), optional:: OCart
  real(8), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0
  if (present(OCart)) lOCart = OCart

  vec(1) = sqrt((x-lOCart(1))**2 + (y-lOCart(2))**2 + (z-lOCart(3))**2)  ! r
  vec(2) = pi - acos((z-lOCart(3)) / vec(1))                             ! theta
  if (lOCart(1) == 0.0 .and. lOCart(2) == 0.0) then
    vec(3) = -atan2(y, x)                                                ! phi
  else
    vec(3) = -atan2((y-lOCart(2)), (x-lOCart(1)))                        ! phi
  end if
  if (vec(3) .lt. 0.0d0) vec(3) = vec(3) + 2.0d0 * pi

END FUNCTION CartesianToSpherical8

FUNCTION CylindricalToCartesian4 (Z, R, Phi, OCart) result (vec)
  real(4), intent(in) :: Z, R, Phi
  real(4), dimension(3), optional :: OCart
  real(4), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0
  if (present(OCart)) lOCart = OCart

  vec(1) = R * cos(Phi) - lOCart(1)  ! x
  vec(2) = R * sin(Phi) - lOCart(2)  ! y
  vec(3) = Z - lOCart(3)             ! z

END FUNCTION CylindricalToCartesian4

FUNCTION CylindricalToCartesian8 (Z, R, Phi, OCart) result (vec)
  real(8), intent(in) :: Z, R, Phi
  real(8), dimension(3), optional :: OCart
  real(8), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0d0
  if (present(OCart)) lOCart = OCart

  vec(1) = R * cos(Phi) - lOCart(1)  ! x
  vec(2) = R * sin(Phi) - lOCart(2)  ! y
  vec(3) = Z - lOCart(3)             ! z

END FUNCTION CylindricalToCartesian8

FUNCTION SphericalToCartesian4 (r, theta, phi, OCart) result (vec)
  real(4), intent(in) :: r, theta, phi
  real(4), dimension(3), optional :: OCart
  real(4), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0
  if (present(OCart)) lOCart = OCart

  vec(1) =  r * sin(theta) * cos(phi) - lOCart(1)  ! x
  vec(2) = -r * sin(theta) * sin(phi) - lOCart(2)  ! y
  vec(3) = -r * cos(theta) - lOCart(3)             ! z

END FUNCTION SphericalToCartesian4

FUNCTION SphericalToCartesian8 (r, theta, phi, OCart) result (vec)
  real(8), intent(in) :: r, theta, phi
  real(8), dimension(3), optional :: OCart
  real(8), dimension(3) :: vec, lOCart
  !.............................................................................
  lOCart = 0.0d0
  if (present(OCart)) lOCart = OCart

  vec(1) =  r * sin(theta) * cos(phi) - lOCart(1)  ! x
  vec(2) = -r * sin(theta) * sin(phi) - lOCart(2)  ! y
  vec(3) = -r * cos(theta) - lOCart(3)             ! z

END FUNCTION SphericalToCartesian8

FUNCTION SphericalToCylindrical4 (r, theta, phi, OCart) result (vec)
  real(4), intent(in) :: r, theta, phi
  real(4), dimension(3), optional :: OCart
  real(4), dimension(3) :: vec, pCart
  !.............................................................................

  ! the most reliable way (that I've found so far) to convert between coord.
  ! systems that are *not guaranteed to have the same origins* is via Cartesian coords.
  pCart = SphericalToCartesian4(r, theta, phi, OCart)
  vec   = CartesianToCylindrical4(pCart(1), pCart(2), pCart(3), OCart)

END FUNCTION SphericalToCylindrical4

FUNCTION SphericalToCylindrical8 (r, theta, phi, OCart) result (vec)
  real(8), intent(in) :: r, theta, phi
  real(8), dimension(3), optional :: OCart
  real(8), dimension(3) :: vec, pCart
  !.............................................................................

  ! the most reliable way (that I've found so far) to convert between coord.
  ! systems that are *not guaranteed to have the same origins* is via Cartesian coords.
  pCart = SphericalToCartesian8(r, theta, phi, OCart)
  vec   = CartesianToCylindrical8(pCart(1), pCart(2), pCart(3), OCart)

END FUNCTION SphericalToCylindrical8

FUNCTION CylindricalToSpherical4 (Z, R, Phi, OCart) result (vec)
  real(4), intent(in) :: Z, R, Phi
  real(4), dimension(3), optional :: OCart
  real(4), dimension(3) :: vec, pCart
  !.............................................................................

  pCart = CylindricalToCartesian4(Z, R, Phi, OCart)
  vec   = CartesianToSpherical4(pCart(1), pCart(2), pCart(3), OCart)

END FUNCTION CylindricalToSpherical4

FUNCTION CylindricalToSpherical8 (Z, R, Phi, OCart) result (vec)
  real(8), intent(in) :: Z, R, Phi
  real(8), dimension(3), optional :: OCart
  real(8), dimension(3) :: vec, pCart
  !.............................................................................

  pCart = CylindricalToCartesian8(Z, R, Phi, OCart)
  vec   = CartesianToSpherical8(pCart(1), pCart(2), pCart(3), OCart)

END FUNCTION CylindricalToSpherical8

!===============================================================================
FUNCTION CurrentToCartesian4 (mtype, x1, x2, x3) result (vec)
  integer, intent(in) :: mtype
  real(4), intent(in) :: x1, x2, x3
  real(4), dimension(3) :: vec
  !.............................................................................
  if (mtype == mesh_types%Cartesian) then
    ! dummy
    vec(1) = x1
    vec(2) = x2
    vec(3) = x3
  else if (mtype == mesh_types%spherical) then
    vec = SphericalToCartesian4(x1, x2, x3)
  else if (mtype == mesh_types%cylindrical) then
    vec = CylindricalToCartesian4(x1, x2, x3)
  end if

END FUNCTION CurrentToCartesian4

FUNCTION CurrentToCartesian8 (mtype, x1, x2, x3) result (vec)
  integer, intent(in) :: mtype
  real(8), intent(in) :: x1, x2, x3
  real(8), dimension(3) :: vec
  !.............................................................................
  if (mtype == mesh_types%Cartesian) then
    ! dummy
    vec(1) = x1
    vec(2) = x2
    vec(3) = x3
  else if (mtype == mesh_types%spherical) then
    vec = SphericalToCartesian8(x1, x2, x3)
  else if (mtype == mesh_types%cylindrical) then
    vec = CylindricalToCartesian8(x1, x2, x3)
  end if

END FUNCTION CurrentToCartesian8

FUNCTION CartesianToCurrent4 (mtype, x1, x2, x3) result (vec)
  integer, intent(in) :: mtype
  real(4), intent(in) :: x1, x2, x3
  real(4), dimension(3) :: vec
  !.............................................................................
  if (mtype == mesh_types%Cartesian) then
    ! dummy
    vec(1) = x1
    vec(2) = x2
    vec(3) = x3
  else if (mtype == mesh_types%spherical) then
    vec = CartesianToSpherical4(x1, x2, x3)
  else if (mtype == mesh_types%cylindrical) then
    vec = CartesianToCylindrical4(x1, x2, x3)
  end if

END FUNCTION CartesianToCurrent4

FUNCTION CartesianToCurrent8 (mtype, x1, x2, x3) result (vec)
  integer, intent(in) :: mtype
  real(8), intent(in) :: x1, x2, x3
  real(8), dimension(3) :: vec
  !.............................................................................
  if (mtype == mesh_types%Cartesian) then
    ! dummy
    vec(1) = x1
    vec(2) = x2
    vec(3) = x3
  else if (mtype == mesh_types%spherical) then
    vec = CartesianToSpherical8(x1, x2, x3)
  else if (mtype == mesh_types%cylindrical) then
    vec = CartesianToCylindrical8(x1, x2, x3)
  end if

END FUNCTION CartesianToCurrent8

!===============================================================================
!> Calculate the volume centroid of `mesh`.
!> Assumes the mesh type, mesh%idir and mesh%s are already set.
!===============================================================================
FUNCTION VolumeCentroid (mesh) result (centre)
  class(mesh_t) :: mesh
  real(4) :: centre
  !.............................................................................
  call trace_begin('VolumeCentroid')

  associate(s => mesh%s, xmn => mesh%llc_nat)
  select type (mesh)
  type is (Cartesian_mesh)
    ! In Cartesian coords., this is simply half the size of the patch.
    centre = xmn + 0.5 * s
  type is (cylindrical_mesh)
    select case (mesh%idir)
    case (1,3)
      centre = xmn + 0.5 * s
    case (2)
      ! radial coord.
      ! = (2/3) delta(R**3) / delta(R**2)
      centre = 2.0 * (3.0 * xmn**2 + 3.0 * xmn * s + s**2) &
                 / (3.0 * (2.0 * xmn + s))
    end select
  type is (spherical_mesh)
    select case(mesh%idir)
    case (1)
      ! radial coord.
      ! = (3/4) delta(r**4) / delta(r**3)
      centre = 3.0 * (4.0 * xmn**3 + 6.0 * xmn**2 * s + 4.0 * xmn * s**2 + s**3) &
                 / (4.0 * (3.0 * xmn**2 + 3.0 * xmn * s + s**2))
    case (2)
      ! theta coord.
      ! = delta(sin(theta) - theta * cos(theta)) / delta(-cos(theta))
      centre = (sin(xmn + s) - (xmn + s) * cos(xmn + s) - sin(xmn) + xmn * cos(xmn)) &
                  / (cos(xmn) - cos(xmn + s))
    case (3)
      centre = mesh%llc_nat + 0.5 * s
    end select
  end select
  end associate

  call trace_end
END FUNCTION VolumeCentroid

!===============================================================================
!> Round a single-precision floating point number to `digit` sig. figs.
!===============================================================================
FUNCTION round4(a,digit) result (b)
  real(4) :: a, b, rof, q1, q2, q3
  integer :: digit
  real(4) :: tiny = 1.0e-32
  integer :: minexp = -6
  integer :: i1,i2
! ----------------------------------------------------------------------

  q1 = abs(a)
  if (q1 <= tiny .or. digit == 0.) then
    b = a
    return
  end if
  q1 = log10(q1)
  q2 = sign(0.5, q1)
  i1 = digit - int(q1 + q2 + 0.5)
  i2 = min(-3, digit + minexp)
  q1 = 10.0**i1
  q2 = 10.0**i2
  q3 = sign(q2, a)
  b = anint(a * q1 + q3) / q1

END FUNCTION round4

!===============================================================================
!> Round a double-precision floating point number to `digit` sig. figs.
!===============================================================================
FUNCTION round8(a,digit) result (b)
  real(8) :: a, b, rof, q1, q2,q3
  integer :: digit
  real(8) :: tiny = 1.0d-99
  integer :: minexp = -14
  integer :: i1,i2
! ----------------------------------------------------------------------

  q1 = abs(a)
  if (q1 <= tiny .or. digit == 0_8) then
    b = a
    return
  end if
  q1 = log10(q1)
  q2 = sign(0.5_8, q1)
  i1 = digit - int(q1 + q2 + 0.5_8)
  i2 = min(-3, digit + minexp)
  q1 = 10.0_8**i1
  q2 = 10.0_8**i2
  q3 = sign(q2, a)
  b = anint(a * q1 + q3) / q1

END FUNCTION round8

END MODULE mesh_mod
