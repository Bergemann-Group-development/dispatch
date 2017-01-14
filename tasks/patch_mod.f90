!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for patches, which adds pointers to memory and mesh, and
!> number of dimensions and variables.  The number of dimensions is hardwired
!> to 3, since we can just choose the actual dimensions to get 2-D and 1-D.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mesg_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE kinds_mod
  USE task_mod
  USE link_mod
  USE mesh_mod
  USE buffered_output_mod
  implicit none
  private
  integer, parameter:: ndim = 3
  type, extends(task_t), public:: patch_t
    real(kind=KindScalarVar), dimension(:,:,:,:,:,:), pointer:: mem      => null()
    real, dimension(:,:,:,:,:),   pointer:: vgas     => null()
    real, dimension(:,:,:,:),     pointer:: density  => null()
    real, dimension(:,:,:,:),     pointer:: pressure => null()
    class(mesh_t), dimension(:), pointer:: mesh => null()
    real:: courant                      ! courant number
    real:: gamma=1.0                    ! gas gamma
    real:: dt_fixed, u_fixed, u_max     ! for fixed time step Courant condition
    integer:: ndim = ndim               ! number of dimensions
    integer:: nv=5                      ! number of variables
    integer:: nw=2                      ! number of work indices
    integer:: msplit=2                  ! how to split a patch when refining
    integer:: max_files                 ! prevent disk overrun when NaN happens
    logical:: guard_zones, time_derivs  ! switches for I/O
    logical:: root_grid=.false.
    logical:: mem_allocated=.false.
    logical:: do_pebbles=.false.
    logical:: staggered=.true.
    class(link_t), pointer:: link       ! this could be useful
    !---------------------------------------------------------------------------
    ! HACKs below:  FIXME
    !---------------------------------------------------------------------------
    integer:: i_d=1, i_s=2, i_e=2, i_px=3, i_py=4, i_pz=5, i_bx=6, i_by=7, i_bz=8, i_tt=9, i_phi=10
    real:: grav                         ! gravitational consant (code units)
    real(8):: ds(3)                     ! convenience
    integer:: n(3)=0, ng(3)=2          ! convenient defaults
    integer:: li(3), ui(3), gn(3)       ! convenience -- avoid patch%mesg%li etc
    integer:: ipos(3)                   ! integer position
    real(8):: offset(3)                 ! this is used to compute indices, locations
    real(8):: llc_cart(3)=0d0           ! the lower left corner of the patch in *Cartesian* coords.
    real(8):: llc_nat(3)=0d0            ! the lower left corner of the patch in its *own* coords.
    real(8):: centre_cart(3)=0d0        ! the geometric centre of the patch in *Cartesian* coords.
    real(8):: centre_nat(3)=0d0         ! the geometric centre of the patch in its *own* coords.
    integer :: mesh_type = 1            ! default: Cartesian mesh; FIXME (there must be a better way)
    real(8):: grid_vel(3) = 0d0         ! velocity of the patch in its native coords. and
                                        ! units (e.g. rad/s for the phi coord.); for
                                        ! Cartesian coords., equal to `velocity`.
  contains
    procedure:: check_error
    procedure:: init
    procedure:: rotate
    procedure:: output                  ! legacy output
    procedure:: output_new              ! new output
    procedure:: output_legacy           ! legacy output
    procedure:: pack                    ! pack patch into buffer
    procedure:: unpack                  ! unpack buffer into patch
    procedure:: info                    ! print info to stdout
    procedure:: point_distances         ! distances btw two points
    procedure:: upload_range            ! returns the index ranges btw two points
    procedure:: dnload_range            ! returns the index ranges btw two points
    procedure:: index                   ! compute the index for a single point
    procedure:: nearest                 ! the nearest index that maps to patch point
    procedure:: lowest                  ! the lowest index that maps to patch point
    procedure:: highest                 ! the highest index that maps to patch point
    procedure:: myposition
    procedure:: myposition_staggered
    procedure:: periodic_grid
    procedure:: extrapolate_guards
    procedure, private :: fminval4
    procedure, private :: fminval8
    procedure, private :: fmaxval4
    procedure, private :: fmaxval8
    generic, public :: fminval => fminval4, fminval8  ! interior min value
    generic, public :: fmaxval => fmaxval4, fmaxval8  ! interior max value
    procedure:: courant_condition
    procedure:: check_density
    procedure:: patch_to_header
    procedure:: header_to_patch
    procedure:: local_iit
    procedure:: time_interpolate
    procedure:: stats
    procedure:: MapVectors
    procedure:: beyond_patch_edge
    procedure:: myfloatindex_scalar
    procedure:: log_interpolate
    procedure:: time_interval
  end type
  !-----------------------------------------------------------------------------
  ! Sequential header, for MPI and I/O
  !-----------------------------------------------------------------------------
  integer, save:: n_header = (14*4 + 7*8 + 5*8*3 + 2*8*8 + 1*4*8 + 2*4*3)/4
  type header_t
    sequence
    integer:: version=1
    integer:: n_header
    integer:: n_data
    integer:: id
    integer:: status
    integer:: rank
    integer:: istep
    integer:: iout
    integer:: level
    integer:: nt
    integer:: nv                        ! number of variables
    integer:: it
    integer:: nq=0
    integer:: npad=0
    real(8):: time
    real(8):: end_time
    real(8):: out_next
    real(8):: out_time
    real(8):: print_next
    real(8):: print_time
    real(8):: dtime
    real(8):: position(3)
    real(8):: velocity(3)
    real(8):: size(3)
    real(8):: box(3)                    ! the periodic wrapping size
    real(8):: ds(3)
    real(8):: send_time
    real(8), dimension(8):: t
    real(8), dimension(8):: dt
    integer, dimension(8):: iit
    integer:: n(3)=16                   ! convenience
    integer:: ng(3)=2                   ! convenience
  end type
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type box_t
    real:: size(3)=1d0
    real:: num(3)
    integer:: dim(3)
  end type
  public task2patch
CONTAINS

!> =============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!> =============================================================================
SUBROUTINE patch_to_header (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt
  !-----------------------------------------------------------------------------
  head%n_header   = n_header
  head%n_data     = product(self%gn)*self%nv
  head%id         = self%id
  head%status     = self%status
  head%rank       = self%rank
  head%istep      = self%istep
  head%iout       = self%iout
  head%level      = self%level
  head%time       = self%time
  head%end_time   = self%end_time
  head%out_next   = self%out_next
  head%out_time   = self%out_time
  head%print_next = self%print_next
  head%print_time = self%print_time
  head%dtime      = self%dtime
  head%position   = self%position
  head%velocity   = self%velocity
  head%size       = self%size
  head%it         = self%it
  head%nt         = self%nt
  nt              = self%nt
  head%t(1:nt)    = self%t
  head%dt(1:nt)   = self%dt
  head%iit(1:nt)  = self%iit
  head%nq         = self%nq
  head%send_time  = wallclock()
  if (io%verbose > 1) &
    write (io_unit%log,'(f12.6,3x,a,i7,i4,l3,2i4,f8.4)') &
      wallclock(), 'patch_to_header: id, nq =', self%id, &
      self%nq, self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n
END SUBROUTINE patch_to_header

!> =============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!> =============================================================================
SUBROUTINE header_to_patch (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt
  !-----------------------------------------------------------------------------
  self%id         = head%id
  self%status     = head%status
  self%rank       = head%rank
  self%istep      = head%istep
  self%iout       = head%iout
  self%level      = head%level
  self%time       = head%time
  self%end_time   = head%end_time
  self%out_next   = head%out_next
  self%out_time   = head%out_time
  self%print_next = head%print_next
  self%print_time = head%print_time
  self%dtime      = head%dtime
  self%position   = head%position
  self%velocity   = head%velocity
  self%size       = head%size
  self%it         = head%it
  self%nt         = head%nt
  nt              = head%nt
  self%t          = head%t(1:nt)
  self%dt         = head%dt(1:nt)
  self%iit        = head%iit(1:nt)
  self%nq         = head%nq
  self%latency    = wallclock() - head%send_time
  if (io%verbose > 1) &
    write (io_unit%log,'(f12.6,3x,a,i7,i4,l3,2i4,f8.4)') &
      wallclock(), 'header_to_patch: id, nq =', self%id, &
      self%nq, self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, self%latency
END SUBROUTINE header_to_patch

!> =============================================================================
!> Check for y-redundancy
!> =============================================================================
SUBROUTINE check_error (self, a, label, error)
  class(patch_t):: self
  character(len=*), optional:: label
  integer, optional:: error(3)
  integer:: loc(3), it, l(3), u(3), n(3), iy
  real, dimension(:,:,:), pointer:: a, b
  integer, save:: nerror=4
  !-----------------------------------------------------------------------------
  if (.not.io%do_debug) return
  if (self%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  n = u-l+1
  allocate (b(n(1),n(2),n(3)))
  do iy=l(2),u(2)
    b(:,iy,:)=a(:,iy,:)-a(:,1,:)
  end do
  if (any(b /= 0.0)) then
    loc = maxloc(abs(b)) + l-1
    if (present(label)) print*,'patch_t::error ', label
    print*,'patch_t::error loc        :', loc
    print*,'patch_t::error it,new,iit :', self%it, self%new, self%iit
    do iy=l(2),u(2)
      print*,'patch_t::error iy, a:', iy, a(loc(1),iy,loc(3))
    end do
    if (present(error)) then
      error=loc
      return
    else
      if (nerror<=0) stop
      nerror=nerror-1
    end if
  else
    if (present(label)) then
      print*,'patch_t::check_error maxval ', label, maxval(abs(a))
    else
      print*,'patch_t::check_error maxval', maxval(abs(a))
    end if
  end if
  deallocate(b)
END SUBROUTINE

!> =============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!> =============================================================================
SUBROUTINE init (self)
  class(patch_t):: self
  integer, save:: n(3)=16, ng(3)=2, nv=5, nt=5, nw, max_files=1000
  logical, save:: first_time = .true.
  logical, save:: guard_zones=.false., time_derivs=.false.
  real, save:: u_fixed=-1.0, dt_fixed=-1.0, grace=0.0, courant=0.2
  real, save:: end_time=1.0, out_time=0.1, print_time=0.1
  real(8):: position(3), size(3)
  integer:: iostat, i
  namelist /patch_params/ size, n, nv, nw, ng, nt, courant, u_fixed, dt_fixed, grace
  namelist /out_params/ out_time, end_time, print_time, max_files, &
    guard_zones, time_derivs
  !.............................................................................
  call trace_begin ('patch_mod::init', 2)
  call self%task_t%init                         ! sets unique id
  !-----------------------------------------------------------------------------
  ! Read and set patch parameters; these may have been set by the calling code
  !-----------------------------------------------------------------------------
  !$omp critical (patch_input_cr)
  if (first_time) then
    size = self%size
    n    = self%n
    ng   = self%ng
    nv   = self%nv
    nw   = self%nw
    nt   = self%nt
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (io%master) write (*, patch_params)
    rewind (io%input)
    read (io%input, out_params)
    if (io%master) write (*, out_params)
  end if
  !print *, 'MK', omp_mythread, self%id, self%nt, nt
  !$omp end critical (patch_input_cr)
  ! reasonable values of self%n should take precedence over a saved value of `n`.
  if (any(self%n == 0)) self%n = n              ! patch dimensions
  self%ng          = ng                         ! number of guard zones
  self%nv          = nv                         ! number of variables
  self%nw          = nw                         ! number of work indices
  self%nt          = nt                         ! number of time slices
  self%courant     = courant                    ! courant number
  self%grace       = grace                      ! grace fraction of dtime
  self%dt_fixed    = dt_fixed
  self%u_fixed     = u_fixed
  self%max_files   = max_files                  ! prevent run-away NaN
  !-----------------------------------------------------------------------------
  ! Output related parameters
  !-----------------------------------------------------------------------------
  self%guard_zones = guard_zones
  self%time_derivs = time_derivs
  self%end_time    = end_time
  self%out_time    = out_time
  self%print_time  = print_time
  self%print_next  = print_time
  self%ds = real(self%size,kind=8) / self%n     ! ds in patch
  if (all(self%size==self%box)) then
    write (io_unit%log,*) mpi%rank, 'task', self%id, ' is root grid'
    self%root_grid = .true.
    call self%set(bits%root_grid)
  end if
  !-----------------------------------------------------------------------------
  ! Special handling of collapsed dimensions
  !-----------------------------------------------------------------------------
  where (self%n==1)                             ! for collapsed dimensions
    self%ng = 0                                 ! no guard zones
    self%box = 1                                ! make sure it is non-zero
    self%ds = minval(self%ds)                   ! then set to reasonable val
    self%size = self%ds                         ! one cell thick
    self%llc_cart = self%position
  end where
  if (self%mesh_type == mesh_types%Cartesian) then
    self%llc_nat = self%llc_cart
  end if
  !-----------------------------------------------------------------------------
  ! Initialize and allocate the meshes, afterwards imposing the external values
  ! box size and position
  !-----------------------------------------------------------------------------
  call omp%set_stacksize (self%nv)
  call MeshFactory(self%mesh_type, self%mesh, n=self%n, ng=self%ng, s=self%size, &
                   p=self%position, llc_native=self%llc_nat)
  ! intialise staggering indices; it is the responsibility of the *solver* to
  ! assign non-zero values to this variable (if any).
  do i=1,3
    allocate(self%mesh(i)%h(self%nv))
    self%mesh(i)%h = 0.0
    if (self%kind == 'zeus_mhd_patch' .and. &
        self%staggered .and. self%n(i) /= 1) self%mesh(i)%gn = self%mesh(i)%gn + 1
  end do
  self%mesh%b = self%box
  self%mesh%llc_cart = self%llc_cart
  self%centre_nat = [self%mesh(1)%centre_nat,self%mesh(2)%centre_nat,self%mesh(3)%centre_nat]
  self%centre_cart = CurrentToCartesian(self%mesh_type, &
                       self%centre_nat(1), self%centre_nat(2), self%centre_nat(3))
  !-----------------------------------------------------------------------------
  ! Hack for backwards compatibility
  !-----------------------------------------------------------------------------
  self%gn = self%mesh%gn
  self%li = self%mesh%li
  self%ui = self%mesh%ui
  self%gsize = (self%size*self%gn/self%n)
  self%offset = (self%mesh%lb+self%mesh%ub)/2.0 ! mid point index offset
  !-----------------------------------------------------------------------------
  ! Measure levels in terms of the reinement factor
  !-----------------------------------------------------------------------------
  self%level = log(maxval(self%box)/minval(self%ds))/log(real(self%msplit)) + 0.1
  if (io%verbose>2) &
    write (io_unit%log,'(1x,a,i7,2(1p,3e12.4,2x))') 'patch_mod::init: id, size, pos =', &
      self%id, self%size, self%position
  if (io%verbose>3) &
    write (io_unit%log,*) 'id, LEVEL: ', self%id, self%level, maxval(self%box), minval(self%ds)
  !-----------------------------------------------------------------------------
  ! Allocate memory, if not already allocated
  !-----------------------------------------------------------------------------
  if (.not.self%mem_allocated) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
!    if (associated(self%mem)) deallocate (self%mem)
!    if (self%rank == mpi%rank .or. self%is_set (bits%virtual)) then
      allocate (self%mem(self%mesh(1)%gn, &
                         self%mesh(2)%gn, &
                         self%mesh(3)%gn, &
                         self%nv,self%nt,self%nw))
      self%mem_allocated = .true.               ! mark
!    end if
  end if
  self%dt = 0.0
  self%t =  0.0
  self%iit = 1                                  ! initial memory slots
  self%iit(nt) = 2                              ! first output slot
  self%it = self%iit(nt-1)                      ! short-hand (=1 initially)
  self%new = self%iit(nt)                       ! short-hand
  self%t(self%it) = 0.0                         ! the only valid slot (no. 1)
  self%wallclock = wallclock()                  ! starting time
  block
  type(header_t):: header
  n_header = storage_size(header)/32             ! 32 bits per word
  if (n_header*4 /= storage_size(header)/8) then
    print *, n_header*4, storage_size(header)/8
    error stop 'n_header is incorrect'
  end if
  end block
  if (all(self%size == self%box)) then
    call self%set (bits%root_grid)
    if (io%verbose > 1) print *, 'set bits%root_grid id =', self%id
  end if
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!> =============================================================================
SUBROUTINE rotate (self)
  class(patch_t):: self
  integer:: iv, nv, new
  real(8):: time
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin ('patch_t%rotate')
  !-----------------------------------------------------------------------------
  ! Periodic wrapping, time extrapolation
  !-----------------------------------------------------------------------------
  if (self%root_grid.or.self%is_set(bits%root_grid)) then
    call self%periodic_grid                     ! periodic root grid
  else
    call self%extrapolate_guards                ! time-extrapolate guard cells
  end if
  !-----------------------------------------------------------------------------
  ! One line summary
  !-----------------------------------------------------------------------------
  if (io%verbose>1) then
    write (io_unit%log,'(a,i6,2i3,2e13.6,2x,a,1p,6e10.2)') 'update: id,j,k,l,n,dt,t', &
      self%id, self%it, self%new, self%dtime, self%time, 'minmax:', &
      minval(self%mem(:,:,:,1,self%it ,self%nw)), &
      maxval(self%mem(:,:,:,1,self%it ,self%nw)), &
      minval(self%mem(:,:,:,1,self%new,1)), &
      maxval(self%mem(:,:,:,1,self%new,1))
  end if
  !-----------------------------------------------------------------------------
  ! Task time slot rotation
  !-----------------------------------------------------------------------------
  call self%task_t%rotate
  call trace_end
END SUBROUTINE rotate

!> =============================================================================
!> This function gives the physical position of an index triple (patch point)
!> in absolute Cartesian coords.
!> =============================================================================
FUNCTION myposition (self, ii) RESULT(out)
  class(patch_t):: self
  integer, intent(in):: ii(3)
  real(8):: out(3), cellpos(3)
  !.............................................................................
  !cellpos = (real(ii-self%li,kind=8) + merge(0.5,0.0,self%n>1)) * self%ds &
  !        + self%llc_nat ! in patch coords.
  !if (self%id == io%id_debug) &
  !  print *,'A', real(cellpos/self%size), real(self%llc_nat)
  cellpos = self%position + (real(ii)-self%offset)*self%ds
  !if (self%id == io%id_debug) &
  !  print *,'B',real(cellpos/self%size)
  ! convert to Cartesian coords.
  if (self%mesh_type == mesh_types%cylindrical) then
    cellpos = CylindricalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  else if (self%mesh_type == mesh_types%spherical) then
    cellpos = SphericalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  end if
  out = cellpos ! absolute position
END FUNCTION myposition

!> =============================================================================
!> Scalar version of the previous function, but **including staggering factor**.
!> For the conversion to Cartesian coords., all three indices are required!
!> =============================================================================
FUNCTION myposition_staggered (self, ii, idir, h) RESULT(out)
  class(patch_t):: self
  integer, intent(in):: ii(3), idir
  real(8):: out, cellpos(3)
  real(8), optional :: h
  real(4) :: hlocal
  !.............................................................................
  if (present(h)) then
    hlocal = h
  else
    hlocal = 0.0
  end if

  cellpos = (real(ii-self%li,kind=8) + merge(0.5 + hlocal,0.0,self%n(idir)>1)) &
          * self%ds + self%llc_nat ! in patch coords.
  ! convert to Cartesian coords.
  if (self%mesh_type == mesh_types%cylindrical) then
    cellpos = CylindricalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  else if (self%mesh_type == mesh_types%spherical) then
    cellpos = SphericalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  end if
  out = cellpos(idir)  ! absolute position

END FUNCTION myposition_staggered

!> =============================================================================
FUNCTION myfloatindex_scalar (self, cpos, idir, h) RESULT(out)
  class(patch_t) :: self
  real(8), intent(in) :: cpos(3)
  integer, intent(in) :: idir
  real(8), intent(in) :: h
  real(8) :: pp(3), ppCart(3), out, intermezzo
  !.............................................................................
  ! convert to native patch coords.
  if (self%mesh_type == mesh_types%Cartesian) then
    ! do nothing
    pp = cpos
  else if (self%mesh_type == mesh_types%spherical) then
    ppCart = cpos - self%position ! Cart. coords to centre of patch's coord. system
    pp = CartesianToSpherical(ppCart(1), ppCart(2), ppCart(3))
  else if (self%mesh_type == mesh_types%cylindrical) then
    ppCart = cpos - self%position ! Cart. coords to centre of patch's coord. system
    pp = CartesianToCylindrical(ppCart(1), ppCart(2), ppCart(3))
  end if

  intermezzo = (pp(idir) - self%llc_nat(idir)) / self%ds(idir) + self%li(idir) &
             - merge(0.5d0 + h,0.0d0,self%n(idir)>1)
  out = intermezzo

END FUNCTION myfloatindex_scalar

!> =============================================================================
!> Max absolute distance between two points in a periodic box
!> Assumes input is Cartesian.
!> =============================================================================
FUNCTION point_distances (self, p1, p2) RESULT (out)
  class(patch_t):: self
  real(8):: p1(3), p2(3), out(3)
  !.............................................................................
  out = abs(modulo(p1-p2+0.5_8*self%box,max(self%box,1d-30))-0.5_8*self%box)
END FUNCTION point_distances

!> =============================================================================
!> Periodic mapping of a grid -- typically the root grid
!> =============================================================================
SUBROUTINE periodic_grid (self)
  class(patch_t):: self
  integer:: ix, iy, iz, i(3), j(3), lb(3), ub(3), li(3), ui(3)
  !.............................................................................
  call trace_begin ('patch_mod::periodic_grid')
  lb = self%mesh%lb
  ub = self%mesh%ub
  li = self%mesh%li
  ui = self%mesh%ui
  where (self%mesh%lower_boundary)
    li = self%mesh%lb
  end where
  where (self%mesh%upper_boundary)
    ui = self%mesh%ub
  end where
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    i = [ix,iy,iz]
    if (all(i >= li .and. i <= ui)) cycle
    j = modulo(i-li,self%n)+li
    !---------------------------------------------------------------------------
    ! Make sure horizontal periodic wrapping continues over the whole vertical
    ! range when physical boundary conditions are present.  Note the use of
    ! mesh%ui and mesh%li below; this is because the physical boundaries are
    ! inside the normal patch boundaries.
    !---------------------------------------------------------------------------
    where (self%mesh%lower_boundary .and. i < self%mesh%li+4)
      j = i
    end where
    where (self%mesh%upper_boundary .and. i > self%mesh%ui-4)
      j = i
    end where
    self%mem(ix,iy,iz,:,self%it ,1) = self%mem(j(1),j(2),j(3),:,self%it ,1)
    self%mem(ix,iy,iz,:,self%new,1) = self%mem(j(1),j(2),j(3),:,self%new,1)
  end do
  end do
  end do
  if (io%verbose > 1) print *, &
    'periodic_grid: id, it, new =', self%id, self%it, self%new
  call trace_end
END SUBROUTINE periodic_grid

!> =============================================================================
!> Time extrapolate guard cell values, to allow using also guard zone values to
!> get boundary conditions for higher level patches.  This is done before time
!> slot rotation, so the slot that has just been updated is iit(nt).
!> =============================================================================
SUBROUTINE extrapolate_guards (self)
  class(patch_t):: self
  integer:: ix, iy, iz, i(3), it0, it1, it2, nt, iit(self%nt)
  integer:: li(3), ui(3)
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin ('patch_t%extrapolate_guards', 3, itimer=itimer)
  nt = self%nt
  iit = self%local_iit()
  it0 = iit(nt)                                         ! current time slot
  it1 = iit(nt-1)                                       ! previous time slot
  it2 = iit(nt-2)                                       ! the one before that
  li = self%mesh%li
  ui = self%mesh%ui
  where (self%mesh%lower_boundary)
    li = self%mesh%lb
  end where
  where (self%mesh%upper_boundary)
    ui = self%mesh%ub
  end where
  do iz=self%mesh(3)%lb,self%mesh(3)%ub                 ! loop over all cells
  do iy=self%mesh(2)%lb,self%mesh(2)%ub
  do ix=self%mesh(1)%lb,self%mesh(1)%ub
    i = [ix,iy,iz]
    if (all(i >= li .and. i <= ui)) cycle               ! skip interior
    self%mem(ix,iy,iz,:,it0,1)     = &
    self%mem(ix,iy,iz,:,it1,1)*2.0 - &
    self%mem(ix,iy,iz,:,it2,1)
  end do
  end do
  end do
  if (io%verbose>1) print *, &
    'extrapolate_guards: id, it =', self%id, it0, self%time
  call trace_end (itimer)
END SUBROUTINE extrapolate_guards

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(patch_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3)
  real:: pt, qt
  real, dimension(:,:,:), pointer:: buf
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  !$omp critical (out_cr)
  call trace_begin('patch_mod::output', itimer=itimer)
  if (self%out_time==0.0) then
    call output1
  else
    do while (self%time >= self%out_next)
      call output1
      self%out_next = (int(self%out_next/self%out_time)+1)*self%out_time
    end do
  end if
  call trace_end (itimer)
  !$omp end critical (out_cr)
  contains
  subroutine output1
    if (io%do_legacy) then
      if (io%verbose>1) print *, self%id, ' LEGACY OUTPUT'
      call self%output_legacy
    else
      call self%output_new
    end if
    self%iout = self%iout+1
  end subroutine output1
END SUBROUTINE output

!> =============================================================================
!> Write results to disk, using the buffered I/O in output_mod
!> =============================================================================
SUBROUTINE output_new (self)
  class(patch_t):: self
  !.............................................................................
  type(header_t):: header
  integer:: n_buf, n_data, ibuf, it1, it2, iv
  real(4), dimension(:), pointer:: buffer
  real, dimension(:,:,:), pointer:: var
  real:: pt
  !-----------------------------------------------------------------------------
  ! Compute size of buffer and allocate -- it will be freed by the writer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_mod::output_new')
  n_buf = product(self%gn)
  n_data = n_header + self%nv*n_buf
  allocate (buffer(n_data), var(self%gn(1),self%gn(2),self%gn(3)))
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info to sequenced header, and copy that to the buffer
  !-----------------------------------------------------------------------------
  call patch_to_header (self, header)
  ibuf = 1
  call anonymous_copy (n_header, header, buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer; these are first interpolated in time
  !-----------------------------------------------------------------------------
  it1 = self%iit(self%nt-2)
  it2 = self%iit(self%nt-1)
  pt = (self%out_next-self%t(it1))/max(self%t(it2)-self%t(it1),1d-30)
  do iv=1,self%nv
    var = self%mem(:,:,:,iv,it1,1)*(1.-pt) + self%mem(:,:,:,iv,it2,1)*pt
    call anonymous_copy (n_buf, var, buffer(ibuf))
    ibuf = ibuf + n_buf
  end do
  deallocate (var)
  call buffered_output%buffer (self%id, self%iout, self%out_next, buffer)
  call buffered_output%check
  call trace_end
END SUBROUTINE output_new

!> =============================================================================
!> Pack patch info into a buffer, for MPI communication
!> =============================================================================
SUBROUTINE pack (self, mesg)
  class(patch_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(header_t):: header
  integer:: n_buf, ibuf, iv, it
  real, dimension(:,:,:), pointer:: var
  real:: pt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer and allocate -- it will be freed by the writer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_mod::pack', itimer=itimer)
  allocate (mesg)
  n_buf = product(self%gn)      ! 4-bytes per element; `anonymous_copy` also uses 4-bytes.
  if (kind(var)==8) &
    n_buf = n_buf * 2           ! 8-bytes per element, but `anonymous_copy` uses 4-bytes!
  mesg%nbuf = n_header + self%nv*n_buf
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    mesg%nbuf = mesg%nbuf+ self%nv*n_buf
  end if
  allocate (mesg%buffer(mesg%nbuf))
  allocate (mesg%reqs(self%n_nbors))
  mesg%id = self%id                                             ! MPI message tag
  mesg%nreq = 0                                                 ! reset counter
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info to sequenced header, and copy that to the buffer
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,1) 'patch_mod::pack, id, nq, time, status:', &
    self%id, self%nq, self%time, &
    self%is_set (bits%internal), &
    self%is_set (bits%boundary), &
    self%is_set (bits%virtual), &
    self%is_set (bits%external), &
    self%is_set (bits%swap_request)
    1 format(a,2i9,f12.6,2x,5l1)
    flush (io_unit%log)
  end if
  call patch_to_header (self, header)
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer.  If the bdry+swap bit are set,
  ! copy as a minimum also the previous time slot
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    call anonymous_copy (n_buf, self%mem(:,:,:,iv,self%it,1), mesg%buffer(ibuf))
    ibuf = ibuf + n_buf
  end do
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    it = mod(self%it-1+self%nt-1,self%nt)+1
    if (io%verbose>1) then
      write (io_unit%log,*) self%id, 'extra time slot packed', it
      flush (io_unit%log)
    end if
    do iv=1,self%nv
      call anonymous_copy (n_buf, self%mem(:,:,:,iv,it,1), mesg%buffer(ibuf))
      ibuf = ibuf + n_buf
    end do
  end if
  if (io%verbose > 2) &
    write (io_unit%log,*) self%id, 'after pack', minval(self%mem(:,:,:,1,self%it,1))
  call trace_end (itimer)
END SUBROUTINE pack

!> =============================================================================
!> Pack patch info into a buffer, for MPI communication
!> =============================================================================
SUBROUTINE unpack (self, buffer)
  class(patch_t):: self
  !.............................................................................
  type(header_t):: header
  integer:: n_buf, ibuf, iv, it
  integer, dimension(:), pointer:: buffer
  real, dimension(:,:,:), pointer:: var
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_mod::unpack', itimer=itimer)
  n_buf = product(self%gn)      ! 4-bytes per element; `anonymous_copy` also uses 4-bytes.
  if (kind(var)==8) &
    n_buf = n_buf * 2           ! 8-bytes per element, but `anonymous_copy` uses 4-bytes!
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info from sequenced header
  !-----------------------------------------------------------------------------
  ibuf = 1
  call anonymous_copy (n_header, buffer(ibuf), header)
  ibuf = ibuf + n_header
  if (io%verbose > 2) write (io_unit%log,1) 'patch_mod::unpk id, nq, time before:', &
    self%id, self%nq, self%time
  call self%header_to_patch (header)
  if (io%verbose > 1) write (io_unit%log,1) 'patch_mod::unpk id, nq, time  after:', &
    self%id, self%nq, self%time, &
    self%is_set (bits%internal), &
    self%is_set (bits%boundary), &
    self%is_set (bits%virtual), &
    self%is_set (bits%external), &
    self%is_set (bits%swap_request)
    1 format(a,2i9,f12.6,2x,5l1)
  !-----------------------------------------------------------------------------
  ! Swap the roles of boundary and virtual on incoming patches.  Also clear the
  ! ready bit, which is an indicator of a patch being in the ready queue, which
  ! an incoming patch is not.
  !-----------------------------------------------------------------------------
  if (self%is_set(bits%boundary)) then
    call self%set   (bits%virtual)
    call self%clear (bits%boundary)
  else if (self%is_set(bits%virtual)) then
    call self%set   (bits%boundary)
    call self%clear (bits%virtual)
  end if
  call self%clear (bits%ready)
  !-----------------------------------------------------------------------------
  ! Copy the variables from the output buffer.  If this is swap of patch roles
  ! expect at least one more time slot.
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,self%it,1))
    ibuf = ibuf + n_buf
  end do
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%virtual)) then
    it = mod(self%it-1+self%nt-1,self%nt)+1
    if (io%verbose>0) then
      write (io_unit%log,*) self%id, 'extra time slot unpacked', it
      flush (io_unit%log)
    end if
    do iv=1,self%nv
      call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,it,1))
      ibuf = ibuf + n_buf
    end do
  end if
  if (io%verbose > 2) &
    write (io_unit%log,*) self%id, 'after unpack', minval(self%mem(:,:,:,1,self%it,1))
  call trace_end (itimer)
END SUBROUTINE unpack

!> =============================================================================
!> Write results to disk
!> =============================================================================
SUBROUTINE output_legacy (self, experiment_name)
  USE iso_fortran_env, only : real32
  class(patch_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3)
  real:: pt, qt
  real, dimension(:,:,:), pointer:: buf
  !-----------------------------------------------------------------------------
  call trace_begin('patch_mod::output_legacy')
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (self%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (self%time_derivs) then
    no = min(2,self%nw)
  else
    no = 1
  end if

  !-----------------------------------------------------------------------------
  ! Write a log file, with files written by this rank
  !-----------------------------------------------------------------------------
  write (logname,'(a)') trim(io%outputname)//'log'
  if (first_time) then
    first_time = .false.
    if (io%master) print*, 'output data file log written to '//trim(logname)
    open (io%data_unit,file=trim(logname),status='unknown')
  else
    open (io%data_unit,file=trim(logname),status='unknown',position='append')
  end if
  write (filename,'(a,i4.4,"_",i4.4)') trim(io%outputname), self%id, self%iout
  write (io%data_unit,'(a)') trim(filename)
  ds = self%box/2.0**io%levelmax
  do i=1,3
    if (self%size(i) > 0.0) then
      li(i) = (self%position(i)+0.5*(self%box(i)-self%size(i)))/ds(i)
      ui(i) = (self%position(i)+0.5*(self%box(i)+self%size(i)))/ds(i) - 1.0
    end if
  end do
  write (io%data_unit,'(1p,e14.6,3(2x,3e14.6),1x,i3.1,1x,e10.3)') self%time, self%position, self%ds, &
        self%velocity, self%level, self%quality
  write (io%data_unit,'(4(3i7,2x))') int(self%box/self%ds+0.5), u-l+1, li, ui
  close (io%data_unit)

  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  if (io%master .and. io%verbose>=0) then
    block
    real(8):: time
    if (self%out_time==0.0) then
      time = self%time
    else
      time = self%out_next
    end if
    print'(a,3x,"time:",2f10.6)', 'patch_mod::output: file= '//trim(filename), time, self%time
    end block
  end if
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  write (io%data_unit,'(1p,e18.10,6(2x,3e18.10)1x,i3.1,1x,e10.3)') self%time, self%position, self%ds, &
        self%velocity, self%llc_nat, self%llc_cart, self%mesh%centre_nat, self%level, self%quality
  write (io%data_unit,'(6(3i10.1,2x),3i2.1,1x,1i2.1)') int(self%box/self%ds+0.5), u-l+1, &
        li, ui, self%n, self%gn, self%ng, self%nv
  write (io%data_unit,'(3(1i8.1,1x),i2.1,1x,a)') self%id, self%iout, self%istep, &
        self%mesh_type, trim(self%kind)
  close (io%data_unit)

  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube.
  !-----------------------------------------------------------------------------
  n = u-l+1
  open (io%data_unit, file=trim(filename)//'.dat', form='unformatted', access='direct', &
    status='unknown', recl=real32*product(n))
  do jo=1,no
    do iv=1,self%nv
      allocate (buf(n(1),n(2),n(3)))
      !-------------------------------------------------------------------------
      ! If we write out time derivatives we want a pair of consistent values
      ! and time derivatives, which we have at iit(nt-2).
      !-------------------------------------------------------------------------
      it1 = self%iit(self%nt-2)
      if (self%out_time == 0.0) then
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,jo)
      !-------------------------------------------------------------------------
      ! Normally, we write out only variable values, which can then be interpolated
      ! to the exact output time, self%out_next, the update of which must come
      ! after this point
      !-------------------------------------------------------------------------
      else if (self%time_derivs) then
        it = it1
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)
        call check_error(self, buf, "output 1")
      else
        it2 = self%iit(self%nt-1)
        it = it2
        pt = (self%out_next-self%t(it1))/max(self%t(it2)-self%t(it1),1d-30)
        if (self%t(it2)==self%t(it1)) pt=0.0
        qt = 1.0-pt
        buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it1,jo)*qt &
            + self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*pt
      end if
      call check_error(self, buf, "output 2")
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),iv,it,jo,minval(buf),maxval(buf),pt
      write (io%data_unit,rec=iv+(jo-1)*self%nv) buf
      deallocate (buf)
    end do
  end do
  close (io%data_unit)

  call trace_end
END SUBROUTINE output_legacy

!> =============================================================================
!> Print info to stdout.  If
!>   verbose==0:        print only id=1 to stdout and io_unit%log
!>   verbose==1:        print only id=1 to stdout but all to io_unit%log
!>   verbose==2:        print iv=1 to stdout and io_unit%log for all id
!>   verbose==3:        print iv=1,nv to stdout and io_unit%log for all id
!> =============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(patch_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
  integer:: iv, nv
  integer, save:: itimer=0
  integer, parameter:: max_lines=25
  integer, save:: n_lines=0
  integer, save:: counts(7)=0
  !-----------------------------------------------------------------------------
  1 format(a,i8,2i4,i3,f10.6,1p,e12.4,0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e10.2,3i5,8i4,l3,l1,o12)
  if (io%verbose>=0) then
    if (self%time>=self%print_next) then
      if (self%print_time > 0) &
        self%print_next = (nint(self%time/self%print_time)+1)*self%print_time
      call trace%begin ('patch_t%info', itimer=itimer)
      nv = 1
      if (io%verbose>2) nv = self%nv
      if (io%verbose>1) then
        do iv=1,nv
          write (io_unit%log,1) &
            'upd',self%id, self%level, iv, omp_mythread, self%time, self%dtime, self%u_max, &
            self%fminval(self%mem(:,:,:,iv,self%it ,self%nw)), &
            self%fmaxval(self%mem(:,:,:,iv,self%it ,self%nw)), &
            self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
            self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
            real(n_update), ntask, nq, &
            mpi_mesg%n_check, mpi_mesg%n_ready, mpi_mesg%n_update, &
            mpi_mesg%n_send, mpi_mesg%n_recv, mpi_mesg%n_unpk, mpi_mesg%n_delay, &
            mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, &
            self%is_set(bits%boundary), self%is_set(bits%virtual)
          if (io%verbose>0 .and. io_unit%mpi/=io_unit%log) then
            !$omp critical (info_cr)
            write (io_unit%mpi,1) &
              'upd',self%id, self%level, iv, omp_mythread, self%time, self%dtime, self%u_max, &
              self%fminval(self%mem(:,:,:,iv,self%it ,self%nw)), &
              self%fmaxval(self%mem(:,:,:,iv,self%it ,self%nw)), &
              self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
              self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
              real(n_update), ntask, nq, &
              mpi_mesg%n_check, mpi_mesg%n_ready, mpi_mesg%n_update, &
              mpi_mesg%n_send, mpi_mesg%n_recv, mpi_mesg%n_unpk, mpi_mesg%n_delay, &
              mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, &
              self%is_set(bits%boundary), self%is_set(bits%virtual)
            flush (io_unit%mpi)
            !$omp end critical (info_cr)
          end if
        end do
      end if
      !-------------------------------------------------------------------------
      ! Print a one-liner if verbose==1 and id==1, or verbose>1, while
      ! accumulating counters also when not printing
      !-------------------------------------------------------------------------
      if (io%verbose>0 .or. self%id==1) then
        if (io%master) &
          counts = counts + &
            [mpi_mesg%n_check, mpi_mesg%n_ready, mpi_mesg%n_update, &
             mpi_mesg%n_send , mpi_mesg%n_recv , mpi_mesg%n_unpk  , &
             mpi_mesg%n_delay]
        if (io%master .or. self%id==1) then
          if (n_lines<=0) then
            !$omp critical (info_cr2)
            print '(a)', '0....+....1....+....2....+....3....+....4....+....5....+....' &
                         //'6....+....7....+....8....+....9....+....0....+....ntsk   nq ' &
                       //' chk rdy upd snt prb unp dly snq rcq  BV'
            n_lines = max_lines
            !$omp end critical (info_cr2)
          end if
          !$omp atomic
          n_lines = n_lines-1
          do iv=1,nv
            print 1, &
            'upd',self%id, self%level, iv, omp_mythread, self%time-self%dtime, self%dtime, self%u_max, &
            self%fminval(self%mem(:,:,:,iv,self%it ,self%nw)), &
            self%fmaxval(self%mem(:,:,:,iv,self%it ,self%nw)), &
            self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
            self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
            real(n_update), ntask, nq, counts, &
            mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, &
            self%is_set(bits%boundary), self%is_set(bits%virtual)
          end do
          counts(:) = 0
        end if
        !-----------------------------------------------------------------------
        ! Reset the message counters after every timestep, if relevant
        !-----------------------------------------------------------------------
        mpi_mesg%n_check=0; mpi_mesg%n_ready=0; mpi_mesg%n_update=0
        mpi_mesg%n_send=0 ; mpi_mesg%n_recv=0 ; mpi_mesg%n_unpk=0
        mpi_mesg%n_delay=0
      end if
      call trace%end (itimer)
    end if
  end if
END SUBROUTINE info

!> =============================================================================
!> Find the p1 target and p2 source index ranges that correspond to the overlap
!> region, for upload from p2 to p1 (p1 has lower level than p2)
!> =============================================================================
SUBROUTINE upload_range (p1, p2, l1, u1, l2, u2)
  class(patch_t):: p1
  class(patch_t):: p2
  integer, dimension(3):: l1, u1, l2, u2
  real, dimension(3):: f, pos1, pos2
  !-----------------------------------------------------------------------------
  ! Lower index limits
  !-----------------------------------------------------------------------------
  !print*, p1%id, ' target position, size:', p1%position, p1%size
  !print*, p2%id, ' source position, size:', p2%position, p2%size
  l2 = max(p2%mesh%li, p2%lowest  (p1, p1%mesh%lb))
  u2 = min(p2%mesh%ui, p2%highest (p1, p1%mesh%ub))
  l1 = p1%nearest (p2, l2)
  u1 = p1%nearest (p2, u2)
  !print*,'l1,l2', l1, l2
  !print*,'u1,u2', u1, u2
END SUBROUTINE upload_range

!> =============================================================================
!> Find the p1 target index ranges that correspond to the overlap region, for
!> download from p2 to p1 (p1 has higher or the same level as p2)
!> =============================================================================
SUBROUTINE dnload_range (target, source, l, u, o)
  class(patch_t):: target
  class(patch_t):: source
  integer, dimension(3):: l, u, o
  real, dimension(3):: f, pos1, pos2, dist
  !-----------------------------------------------------------------------------
  if (target%level == source%level) then
    l = target%mesh%li
    u = target%mesh%ui
    o = 0
    dist = source%distance(target)  ! positive if source to the right of target
    where (dist > 0.0)
      l = target%mesh%uo
      u = target%mesh%ub
      o = -target%mesh%n
    else where (dist < 0)
      l = target%mesh%lb
      u = target%mesh%lo
      o = target%mesh%n
    end where
  else
    print*,'dnload_range: ERROR, not implemented yet'
    stop
  end if
END SUBROUTINE dnload_range
!> =============================================================================
!> The nearest smaller index position in the self patch coordinate system.  The
!> point with index self%mesh%lb is half a cell inside the patch guard size, so the
!> floating point position is counted relative to that point, and is split into
!> an integer part and a remainder.  For patches smaller than the box, the index
!> may become both smaller than self%mesh%lb and larger than self%mesh%ub.  For the root
!> grid patch, the index is in the range (li-1,ui) inside the interior.
!> =============================================================================
FUNCTION index (self, other, oindex, f)
  class(patch_t):: self, other
  real, dimension(3):: position, f, pos, pos_lb
  integer, dimension(3):: index, oindex
  !.............................................................................
  pos = other%position + (oindex-other%offset)*other%ds         ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  index = f                                                     ! smaller integer
  f = f - index                                                 ! remainder
END FUNCTION index
!> =============================================================================
!> The nearest index position in the self patch of the point in another patch.
!> The cost of this operation is 4 plus and 2 multiply per direction = 18 flops
!> =============================================================================
FUNCTION nearest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset)*other%ds          ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! nearest integer
END FUNCTION nearest
!> =============================================================================
!> The lowest index position in the self patch of the point in another patch
!> =============================================================================
FUNCTION lowest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset-0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! lowest integer
END FUNCTION lowest
!> =============================================================================
!> The highest index position in the self patch of the point in another patch
!> =============================================================================
FUNCTION highest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset+0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! highest integer
END FUNCTION highest

!> =============================================================================
!> Courant condition, one way or another
!> =============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(patch_t):: self
  logical, optional:: detailed_timer
  real(8):: dtime
  integer, save:: itimer=0
  !.............................................................................
  if (present(detailed_timer)) then
    if (detailed_timer) then
      call trace_begin('patch_mod::courant', itimer=itimer)
    else
      call trace_begin('patch_mod::courant')
    end if
  else
    call trace_begin('patch_mod::courant')
  end if
  if (io%verbose>1) &
    write (io_unit%log,*) self%id, ' courant: u_max, ds', self%u_max, minval(self%ds)
  if (self%dt_fixed > 0.0) then
    dtime = self%dt_fixed
    self%courant = dtime*self%u_max/minval(self%ds)
  else if (self%u_fixed > 0.0) then
    dtime = self%courant*minval(self%ds)/self%u_fixed
  else
    dtime = self%courant*minval(self%ds)/self%u_max
  endif
  !$omp atomic write
  self%dtime = dtime
  !-----------------------------------------------------------------------------
  ! Prevent large (cheap) patches from taking time steps longer than the out_time
  !-----------------------------------------------------------------------------
  if (self%out_time > 0.0) &
    self%dtime = min (self%dtime, self%out_time)
  self%dt(self%it) = self%dtime
  if (io%verbose>3) &
    print'(a,i7,f12.6,f8.4,1p,e10.2)', 'courant:', &
      self%id, self%dtime, self%courant, self%u_max
  call trace_end (itimer)
END SUBROUTINE courant_condition

!> =============================================================================
FUNCTION fmaxval4 (self, f) RESULT (fmaxval)
  real:: fmaxval
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION fmaxval4
!> =============================================================================
FUNCTION fminval4 (self, f) RESULT(fminval)
  real:: fminval
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION fminval4
!> =============================================================================
FUNCTION fmaxval8 (self, f) RESULT (fmaxval)
  real(8):: fmaxval
  class(patch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION fmaxval8
!> =============================================================================
FUNCTION fminval8 (self, f) RESULT(fminval)
  real(8):: fminval
  class(patch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION fminval8

!> =============================================================================
!> Check if the density is negative somewhere, counting points in the interior
!> and guard zones separately.  Attempt to patch the guard zones, if applicable.
!> =============================================================================
SUBROUTINE check_density (self, label)
  class(patch_t):: self
  character(len=*), optional:: label
  real(kind=KindScalarVar):: mind, posd
  real:: tiny = 1.0e-32
  integer:: ix, iy, iz, ii(3), ntot, nins, l(3), u(3)
  logical:: inside
  !-----------------------------------------------------------------------------
  if (io%verbose < 1) return
  if (.not. associated(self%mem)) return
  l = self%mesh%lb
  u = self%mesh%ub
  mind = minval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,1))
  if (present(label) .and. io%verbose>1) write (io_unit%log,*) label, self%id, mind
  if (mind <= 0.0) then
    !---------------------------------------------------------------------------
    ! Negative densities appear to develop at a few points in the guard zones,
    ! so this tries to patch up, sampling densitites in the guard zone.
    !---------------------------------------------------------------------------
    posd = tiny
    ntot = 0
    nins = 0
    do iz=1,self%gn(3)
    do iy=1,self%gn(2)
    do ix=1,self%gn(1)
      ii = [ix,iy,iz]
      inside = all([ix,iy,iz]>self%mesh%lo.and.ii<self%mesh%uo)
      posd = merge (posd, max(posd,self%mem(ix,iy,iz,1,self%it,1)), inside)
      ntot = ntot + merge(1,0,self%mem(ix,iy,iz,1,self%it,1)<=0.0)
      nins = nins + merge(1,0,self%mem(ix,iy,iz,1,self%it,1)<=0.0 .and. inside)
      !-------------------------------------------------------------------------
      ! Pick the sampled positive value if a negative one is found in the
      ! guard zone
      !-------------------------------------------------------------------------
      !self%mem(ix,iy,iz,1,self%it,1) = merge (self%mem(ix,iy,iz,1,self%it,1), &
      !  posd, self%mem(ix,iy,iz,1,self%it,1) > 0.0 .or. inside)
    end do
    end do
    end do
    1 format(a,i8,3x,a,i8,a,i8,a)
    if (present(label)) then
      if (io%verbose>1 .or. nins>0) &
        write (io_unit%log,1) 'ID ',self%id, &
         'WARNING: density non-positive in',ntot,' points,',nins,' inside, '// &
         trim(label)
    else
      if (io%verbose>1 .or. nins>0) &
        write (io_unit%log,1) 'ID ',self%id, &
         'WARNING: density non-positive in',ntot,' points,',nins,' inside'
    end if
  end if
END SUBROUTINE check_density

!> =============================================================================
!> =============================================================================
FUNCTION local_iit (self) RESULT (iit)
  class(patch_t):: self
  integer:: i, it, iit(self%nt)
  !-----------------------------------------------------------------------------
  ! Make a consistent local copy of self%iit, based on self%it, which is
  ! written and read in atomic operations.  The array iit is such that its last
  ! value, iit(self%nt) is equal to self%new.
  !-----------------------------------------------------------------------------
  !!$omp atomic read
  !it = self%it
  !iit = [(1+mod(it+i,self%nt), i=1,self%nt)]
  iit = self%iit
END FUNCTION local_iit

!> =============================================================================
!> Find the right time interval in source%t(source%iit), which are the source
!> times in increasing order.  To avoid having to lock access to source%iit while
!> reading it, we read instead (with an atomic read) source%new, from which a
!> consisten local copy of itt can be constructed.
!> =============================================================================
SUBROUTINE time_interpolate (self, source, jt, pt)
  class(patch_t):: self, source
  integer:: jt(2)
  integer:: i, id, it, iit(self%nt), verbose
  real(8):: pt(2), dt1, t(self%nt), dt(self%nt)
  !.............................................................................
  call trace_begin ('task_t::time_interpolate',3)
  verbose = io_unit%verbose
  !-----------------------------------------------------------------------------
  ! In case of e.g. an error in initial conditions, choose a valid default
  !-----------------------------------------------------------------------------
  !!omp critical (time_cr)
  iit = source%iit
  t   = source%t
  dt  = source%dt
  id  = source%id
  !!omp end critical (time_cr)
  jt = iit(1:2)
  !-----------------------------------------------------------------------------
  ! As long as the first source time is earlier than, or equal to the target
  ! time, renew the interval.  Initially, when all times are zero, we are left
  ! with the interval nt-2:nt-1 and dt=0d0, which means slot nt-1 is used.
  !-----------------------------------------------------------------------------
  do i=2,self%nt-2
    if (t(iit(i)) <= self%time) then
      jt = iit(i:i+1)
    end if
  end do
  !-----------------------------------------------------------------------------
  ! A vanishing source time interval happens only at the initial condition ...
  !-----------------------------------------------------------------------------
  dt1 = t(jt(2))-t(jt(1))
  if (dt1==0d0) then
    pt(2) = 1d0
  !-----------------------------------------------------------------------------
  ! If a task lags behind, prevent backwards time interpolation, with a warning
  !-----------------------------------------------------------------------------
  else
    pt(2) = (self%time - t(jt(1)))/dt1
    if (pt(2) < -self%grace) then
      write (io_unit%log,1) &
        'patch_mod::time_interpolate WARNING: backward extrapolation, trg,src,pt:', &
        self%id, id, pt(2)
      1 format(a,2i9,f7.2,1p,3g16.7)
      pt(2) = -self%grace
      verbose = 3
    else if (pt(2) > 1.001+self%grace) then
      write (io_unit%log,1) &
        'patch_mod::time_interpolate WARNING: forward extrapolation, trg,src,pt:', &
        self%id, id, pt(2), dt1, dt(jt(1)), dt(jt(2))
      pt(2) = 1.0+self%grace
      verbose = 3
    end if
  end if
  pt(1) = 1d0-pt(2)
  if (verbose > 2) then
    write (io_unit%log,'(a,2i5,2f7.3,1p,g16.6)') &
      'patch_mod::time_interpolation jt, pt, target', &
      jt, pt, self%time
    do i=1,self%nt
      write (io_unit%log,'(a,3i4,1p,2g16.6)') &
        'patch_mod::time_interpolation i,iit,src%iit,src%t,src%dt:', &
        i, iit(i), iit(i), t(i), dt(i)
    end do
    flush (io_unit%log)
  end if
  !print'(a,5i3,3x,2i3,2x,2f7.3)', 'MK:', iit-source%iit, jt, pt
  !-----------------------------------------------------------------------------
  ! Facilitate exact time evolution when dt_fixed is being used
  !-----------------------------------------------------------------------------
  if (self%dt_fixed > 0.0) then
    if (abs(pt(1))< 1d-5) then
      if (io_unit%verbose>2) write (io_unit%log,*) 'EXACT'
      pt(1)=0d0; pt(2)=1d0
    else if (abs(pt(2)) < 1d-5) then
      if (io_unit%verbose>2) write (io_unit%log,*) 'EXACT'
      pt(1)=1d0; pt(2)=0d0
    end if
  end if
  call trace_end
END SUBROUTINE time_interpolate

!> =============================================================================
FUNCTION task2patch (task) RESULT (patch)
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  select type (task)
  class is (patch_t)
  patch => task
  class default
  nullify(patch)
  stop
  end select
END FUNCTION

!> =============================================================================
SUBROUTINE stats (self, label)
  class (patch_t):: self
  character(len=*):: label
  integer:: iv
  real(8), dimension(4):: si, sb
  do iv=1,self%nv
    call scalar_stats (self%mem(:,:,:,iv,self%it,1), self%mesh%li, self%mesh%ui, si)
    call scalar_stats (self%mem(:,:,:,iv,self%it,1), self%mesh%lb, self%mesh%ub, sb)
    write (io_unit%log,'(a,i3,1p,2(3x,4g12.4))') &
      trim(label)//' stats: iv, min, aver, rms, max =', iv, si, sb
  end do
contains
  subroutine scalar_stats (f, l, u, s)
    real(kind=KindScalarVar), dimension(:,:,:):: f
    integer:: l(3), u(3)
    real(8):: s(4)
    integer:: ix, iy, iz
    s(2) = 0d0
    s(3) = 0d0
    s(1) = f(1,1,1)
    s(4) = f(1,1,1)
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      s(2) = s(2) + f(ix,iy,iz)
      s(3) = s(3) + f(ix,iy,iz)**2
      s(1) = min(s(1),f(ix,iy,iz))
      s(4) = max(s(4),f(ix,iy,iz))
    end do
    end do
    end do
    s(2) = s(2)/product(u-l+1)
    s(3) = sqrt(max(0d0,s(3)/product(u-l+1)-s(2)**2))
  end subroutine
END SUBROUTINE

!> =============================================================================
FUNCTION MapVectors (self, in, target, ii, jj) result (out)
  class(patch_t) :: self
  class(patch_t), pointer :: target
  real(kind=KindScalarVar), dimension(:), intent(in):: in
  real(kind=KindScalarVar):: out(size(in))
  integer, intent(in) :: ii(3), jj(3)
  !.............................................................................
  out = in
  error stop 'Not implemented error! patch_t%MapVectors'
END FUNCTION MapVectors

!> =============================================================================
FUNCTION beyond_patch_edge (self, idir, position) result (outside)
  class(patch_t) :: self
  integer, intent(in) :: idir
  real(8), intent(in) :: position(3)
  logical :: outside
  real(8) :: edge, curvpos(3)
  !.............................................................................

  if (self%n(idir) == 1) then
    outside = .false.
    return
  endif

  select case (self%mesh_type)
  case (1) ! Cartesian
    edge = 0.5 * self%gsize(idir)
    outside = abs(position(idir)-self%position(idir)) > edge
  case (2) ! spherical
    curvpos = position - self%centre_cart
    curvpos = CartesianToSpherical(curvpos(1), curvpos(2), curvpos(3))
    select case (idir)
    case (1)
      edge = self%gsize(idir) + self%llc_nat(idir)
      outside = abs(curvpos(idir)) > edge .or. abs(curvpos(idir)) < self%llc_nat(idir)
    case (2)
      if (associated(self%mesh)) then
        edge = self%gsize(idir) + self%mesh(2)%xi
      else
        edge = self%gsize(idir)
      end if
      outside = abs(curvpos(idir)) > edge
    case (3)
      if (associated(self%mesh)) then
        edge = self%gsize(idir) + self%mesh(3)%zeta
      else
        edge = self%gsize(idir)
      end if
      outside = abs(curvpos(idir)) > edge
    end select
    outside = abs(curvpos(idir)) > edge
  case (3) ! cylindrical
    curvpos = position - self%centre_cart
    curvpos = CartesianToCylindrical(curvpos(1), curvpos(2), curvpos(3))
    select case (idir)
    case (1)
      edge = 0.5 * self%gsize(idir)
      outside = abs(curvpos(idir)) > edge
    case (2)
      edge = self%gsize(idir) + self%llc_nat(idir)
      outside = abs(curvpos(idir)) > edge .or. abs(curvpos(idir)) < self%llc_nat(idir)
    case (3)
      if (associated(self%mesh)) then
        edge = self%gsize(idir) + self%mesh(3)%zeta
      else
        edge = self%gsize(idir)
      end if
      outside = abs(curvpos(idir)) > edge
    end select
  end select
END FUNCTION beyond_patch_edge

!> =====================================================================
!> Interpolate -- only in time for now -- from the MHD time slots save in
!> this module to the time requested from the solve procedure
!> =====================================================================
SUBROUTINE log_interpolate (self, time, i, a)
  class (patch_t)                 :: self
  real(8)                         :: time
  integer                         :: i
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: a, b
  !.....................................................................
  integer                         :: jt(2)
  real                            :: pt(2)
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%log_interpolate', itimer=itimer)
  call self%time_interval (time, jt, pt)
  if (io%verbose > 1) then
    b => self%mem(:,:,:,i,jt(1),1)
    print *, self%id, 'log_interp(1)', jt(1), minval(b), maxval(b)
    b => self%mem(:,:,:,i,jt(2),1)
    print *, self%id, 'log_interp(2)', jt(2), minval(b), maxval(b)
  end if
  a = log(self%mem(:,:,:,i,jt(1),1))*pt(1) &
    + log(self%mem(:,:,:,i,jt(2),1))*pt(2)
  call trace%end (itimer)
END SUBROUTINE log_interpolate

!> =====================================================================
!> Search for a time interval to interpolate in, and return the time
!> slot indices and weight factors
!> =====================================================================
SUBROUTINE time_interval (self, time, jt, pt)
  class (patch_t)                 :: self
  real(8)                         :: time
  integer                         :: jt(2)
  real                            :: pt(2)
  !.....................................................................
  integer                         :: i, j
  real(8)                         :: t(self%nt), dt
  ! --------------------------------------------------------------------
  call trace%begin ('rmhd_t%time_interval')
  t = self%t(self%iit)
  j = 1
  do i=1,self%nt-2
    if (t(i) < time) j=i
  end do
  jt = self%iit(j:j+1)
  dt = t(j+1)-t(j)
  if (dt==0.0) then
    pt(2) = 0.0
  else
    pt(2) = (time-t(j))/dt
  end if
  pt(1) = 1.-pt(2)
  if (io%verbose>1 .and. pt(2)>0.0) then
    print '(a,1p,3e12.4,2i5,0p,2f8.4)', 'time_interval: times  =', t(j), time, t(j+1), jt, pt
  end if
  call trace%end
END SUBROUTINE time_interval

END MODULE
