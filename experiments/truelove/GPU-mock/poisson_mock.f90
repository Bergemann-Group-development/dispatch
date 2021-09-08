!===============================================================================
! Function to do simple timing
!===============================================================================
FUNCTION wallclock() result (time)
  real(kind=8)           :: time
  real(kind=8), save     :: offset
  real(kind=8), external :: omp_get_wtime
  !.............................................................................
  time = omp_get_wtime()
  if (offset == 0d0) then
    offset=time
  end if
  time = time-offset
END FUNCTION wallclock 
!===============================================================================
! General parameters needed to run experiment 
!===============================================================================
MODULE params_mod
  implicit none
  integer, parameter, public :: realkind=8 ! = 4 (8) for  single (double) precision calculation
  type :: io_t
    integer :: io_input=23
    character(len=200) :: run_dir='./'
  end type
  type(io_t) :: io
  !
  integer, parameter, public :: ndim = 3
  integer, dimension(ndim) :: n=16  ! number of interior points in a patch
  integer, dimension(ndim) :: ng=2  ! number of cells in a ghostzone of a patch
  !
  real(kind=realkind) :: Gnewton=1.0_realkind ! G for gravity
END MODULE params_mod
!===============================================================================
! Cut-down version of the list, patch, and mesh modules to the minimum necessary
!===============================================================================
MODULE patch_mod
  USE params_mod
  implicit none
  private
  ! Mesh variables.
  !    ghost         interior    ghost
  !   |+++++|+++++++++++++++++++|+++++|
  !    |   | |                 | |   |
  !   lb  lo li               ui uo  ub
  !   Total grid size == gn = 2xghost-zone + interior = 2 ng + n
  !
  type, public:: patch_t
    integer, dimension(ndim) :: lb, ub, li, ui, n, gn, lo, uo
    integer :: id
    real(kind=realkind), dimension(ndim) :: s, ds
    real(kind=8),        dimension(ndim) :: p
    logical :: is_periodic
    logical, dimension(4*ndim) :: ready ! if all true, patch is ready for update
    !
    real(kind=realkind), dimension(:,:,:), allocatable :: density ! density inside patch
  contains
    procedure :: init
    procedure :: fsum
    procedure :: fmaxval
  end type 
  !
  type, public:: link_t
    type(patch_t), pointer :: p
    type(link_t),  pointer :: next, prev
  end type
  !
  type(link_t), pointer :: head, tail
  !
  integer :: id=0                   ! running id counter, used for unique patch ids
  !
CONTAINS
!===============================================================================
! Setup data for a patch
!===============================================================================
SUBROUTINE init (self, s, p)
  implicit none
  class(patch_t) :: self
  real(kind=realkind), dimension(ndim), intent(in) :: s
  real(kind=8),        dimension(ndim), intent(in) :: p
  !
  ! Setup patch indices on number of interior points and size of ghost zones
  !
  self%n  = n               ! number of interior grid points
  self%ng = ng              ! number of ghost zones
  self%lb = 1               ! lowest grid point
  self%li = self%lb + ng    ! lower grid point of interior mesh
  self%ui = self%li + n - 1 ! upper grid point of interior mesh
  self%lo = self%li - 1     ! upper grid point in lower ghostzone
  self%uo = self%ui + 1     ! lower grid point in upper ghostzone
  self%ub = self%ui + ng    ! upper grid point in upper ghostzone
  self%gn = self%ub         ! total number of grid points in patch
  !
  self%s  = s               ! size of patch in code units
  self%ds = s / n           ! size of a grid cell in code units
  self%p  = p               ! position of center of patch in code units
  !
  self%is_periodic = .false. ! by default patch is non-periodic. Initialising function can make root-patch periodic by hand
  !
  !Use global counter to obtain unique id and increment it in a threadsafe manner
  !
  !$omp atomic capture
  self%id = id
  id = id + 1
  !$omp end atomic
END SUBROUTINE init
!===============================================================================
! Find maxval in the inside of a patch 
!===============================================================================
REAL(kind=realkind) FUNCTION fmaxval(self,f)
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  fmaxval = maxval(f(self%li(1):self%ui(1),self%li(2):self%ui(2),self%li(3):self%ui(3)))
END FUNCTION fmaxval 
!===============================================================================
! Find sum in the inside of a patch 
!===============================================================================
REAL(kind=realkind) FUNCTION fsum(self,f)
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  real(kind=8) :: s1, s2, s3
  integer, i, j, k
  !-----------------------------------------------------------------------------
  ! sum by cells, pencils, and planes to minimise round-off error
  s3 = 0.0_8
  do k = self%li(3),self%ui(3)
    s2 = 0.0_8
    do j = self%li(2),self%ui(2)
      s1 = 0.0_8
      do i = self%li(1),self%ui(1)
        s1 = s1 + f(i,j,k)
      enddo
      s2 = s2 + s1
    enddo
    s3 = s3 + s2
  enddo
  fsum = real(s3,kind=realkind)
END FUNCTION fsum 
!
END MODULE
!===============================================================================
!> Methods for solving the Poisson equation; partially derived from routines
!> developed by Troels Haugbølle for the course Computational Astrophysics.
!>
!> Boundary conditions are either periodic, in which case the mean source term
!> must vanish, or else the potential is considered given in the guard zones.
!>
!> The error is measured in terms of the equivalent density of the residual,
!> written as residual = (Laplace(phi)- 4 pi G rho) / (4 pi G).  The error
!> is then normalized with the maximum of the actual density, and a floor
!> value.  Note that the floor value must be expressed in code density units.
!===============================================================================
MODULE poisson_mod
  USE patch_mod
  implicit none
  private
  type, public:: poisson_t
    integer:: id=0
    real(8):: omega_sor ! Overrelaxation parameter
    real(8):: rjac2     ! Jacobi spectral radius squared
    procedure(sor), pointer:: update => null()
  contains
    procedure:: init
    procedure:: sor
    procedure, nopass:: dump
    procedure, nopass:: zero_boundaries
  end type
  real,    save:: floor=1e-2, tolerance=1e-4
  real(kind=realkind), save :: fourPiG=4.0_realkind*acos(-1.0_realkind)
  integer, save:: verbose=0
  integer, save:: max_iter=500
  integer, save:: dump_unit=0
  logical, save:: do_dump=.false.
  real(kind=realkind), parameter:: pi = asin(1.0) * 2.0
  logical, save:: chebyshev = .true. ! Activate Chebyshev acceleration of SOR?
  type(poisson_t), public:: poisson
CONTAINS

!===============================================================================
!> Initialize poisson_params
!===============================================================================
SUBROUTINE init (self)
  implicit none
  class(poisson_t)  :: self
  !.............................................................................
  character(len=12), save :: solver='sor'
  logical,           save :: first_time = .true.
  namelist /poisson_params/ tolerance, floor, max_iter, fourPiG, &
                            verbose, solver, chebyshev, do_dump
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input,poisson_params)
    write (*, poisson_params)
  end if
  select case(trim(solver))
  case ('sor')
    self%update => sor
  case default
    write(*,*) 'Poisson solver '//trim(solver)//' is unknown'
    write(*,*) 'Valid solvers are: sor'
    stop
  end select
END SUBROUTINE init
!===============================================================================
!> 2nd order successive over relaxation with red-black ordering
!===============================================================================
SUBROUTINE sor (self, patch, source, phi1, d, lb, ub)
  class(poisson_t)  :: self
  class(patch_t)    :: patch
  real(var_k), dimension(:,:,:)             :: phi1
  real(var_k), dimension(:,:,:), pointer, contiguous :: source
  real(var_k), dimension(:,:,:), intent(in) :: d
  integer, optional:: lb(3), ub(3)
  !.............................................................................
  integer                :: i, j, k, sweep, ss, ione, jone, kone, m(3), iter, &
                            niter, unit
  integer, dimension(3)  :: l, u
  real                   :: error, error0
  real(8)                :: a, b, b1, b2, b3, res, numer, denom, wc, usppt
  integer, save          :: itimer=0
  real(8), dimension(:,:,:), pointer, contiguous:: phi
  !
  !-----------------------------------------------------------------------------
  unit = merge (stdout, stdlog, verbose > 1)
  associate(n=>patch%n, ds=>patch%ds, gn=>patch%gn)
  allocate (phi(gn(1),gn(2),gn(3)))
  if (present(lb)) then
    l = lb
    u = ub
  else
    l = patch%li
    u = patch%ui
  end if
  !
  !-----------------------------------------------------------------------------
  ! Calculate optimal guess for the SOR prefactor depending on
  ! dimensionality, ds, and grid points n:
  !    delta = (1/sum(1/ds^2)) sum(1/ds^2 * cos(pi/n))
  !    omega_sor = 2 / (1 + (1 - delta^2)^0.5)
  ! For an even dimensionality grid
  !    delta = cos(pi / n)
  ! which for n > 10 is well approximated by
  !    delta ~ 1 - 1/2 (pi/n)**2
  ! and then
  !    omega_sor = 2 / (1 + pi/n)
  !-----------------------------------------------------------------------------
  numer = 0.0
  denom = 0.0
  do i=1,3
    if (patch%n(i) > 1) then
      numer = numer + cos(pi / patch%n(i))
      denom = denom + 1d0
    end if
  end do
  self%rjac2 = (numer / denom)**2
  if (chebyshev) then
    ! initial value for Chebyshev acceleration
    self%omega_sor = 1.0
  else
    ! optimal value for the relaxation parameter
    self%omega_sor = 2d0 / (1d0 + sqrt(1d0 - self%rjac2))
  end if
  !
  wc = wallclock()
  phi = phi1
  !-----------------------------------------------------------------------------
  ! If the root grid patch is periodic the mean density MUST be zero and
  ! the guard zones values must be consistent with the periodicity. Refined
  ! patches are then not periodic.
  !-----------------------------------------------------------------------------
  ione = 1
  jone = 1
  kone = 1
  if (n(1) <= 1) ione = 0
  if (n(2) <= 1) jone = 0
  if (n(3) <= 1) kone = 0
  !
  ! Calculate 1 / [(2/dx + 2/dy + 2/dz)] prefactor according to dimensionality
  a = 0.0
  if (n(1) > 1) a = a + 2. / ds(1)**2
  if (n(2) > 1) a = a + 2. / ds(2)**2
  if (n(3) > 1) a = a + 2. / ds(3)**2
  a = 1. / a
  b1 = a / ds(1)**2 * ione
  b2 = a / ds(2)**2 * jone
  b3 = a / ds(3)**2 * kone
  !call self%zero_boundaries (patch, source)
  !.............................................................................
  call dump_init (gn(1:3), l, u, a, b1, b2, b3)
  call self%dump (source, 0)
  !
  do iter=1,max_iter
    error = 0.0
    do sweep=0,1
      do k=l(3),u(3)
        do j=l(2),u(2)
          ss = modulo(j + k + sweep,2)
          do i=l(1)+ss,u(1),2
            !-------------------------------------------------------------------
            ! This form of the residual has no factor in front of phi, which
            ! makes it clear how to apply the over-relaxation factor
            !-------------------------------------------------------------------
            res = b1 * (phi(i+ione,j,k) + phi(i-ione,j,k)) &
                + b2 * (phi(i,j+jone,k) + phi(i,j-jone,k)) &
                + b3 * (phi(i,j,k+kone) + phi(i,j,k-kone)) &
                - phi(i,j,k) - a * source(i,j,k)
            phi(i,j,k) = phi(i,j,k) + self%omega_sor * res
            error = max(abs(res / (d(i,j,k) + floor)), error)
          end do
        end do
      end do
      if (chebyshev) then
        if (iter == 1) then
          self%omega_sor = 1.0 / (1.0 - self%rjac2 * 0.5)
        else
          self%omega_sor = 1.0 / (1.0 - self%rjac2 * 0.25 * self%omega_sor)
        end if
      end if
    end do
    if (verbose > 1) &
      call self%dump (phi, iter)
    !---------------------------------------------------------------------------
    ! We are supposed to return the residual as (Laplace(phi)-source)/d, so we
    ! need to divide with the factor that density was multiplied by
    !---------------------------------------------------------------------------
    error = error / (a*fourPiG)
    if (iter==1) error0 = error
    !...........................................................................
    if (verbose > 2) then
      write (unit,'(f12.6,a,i6,f12.6)') &
        wallclock(), pad_trace(1:indent)// &
        'poisson_t%sor: iter, error =', &
        iter, error
    end if
    !...........................................................................
    if (error < tolerance) exit
  end do
  !
  if (iter==max_iter+1) then
    write (stdout,'(a,i6,i5,1p,2e11.2)') &
      'poisson_t%sor: no convergence, id,niter,error0,error=', &
      patch%id, max_iter, error0, error
  else if (verbose > 0) then
    usppt = 1e6 * (wallclock() - wc) / product(shape(phi))
    !...........................................................................
    write (unit,'(f12.6,a,i7,i5,3x,2e12.3,f8.2,1p,g12.2)') &
      wallclock(), pad_trace(1:indent)// &
      'poisson_t%sor: id, niter, initial/final error, us/pt, 4piG=', &
      patch%id, iter, error0, error, usppt, fourPiG
    !...........................................................................
  end if
  niter = iter
  phi1 = phi
  call self%dump (phi1, -2)
  deallocate (phi)
  end associate
contains
  subroutine dump_init (m, l, u, a, b1, b2, b3)
    integer, dimension(3):: m, l, u
    real(8):: a, b1, b2, b3
    integer:: sentinel
    !-----------------------------------------------------------------------------
    if (do_dump) then
      if (dump_unit==0) then
        open (newunit=dump_unit, file=trim(io%run_dir)//'sor.dump', &
              form='unformatted', status='unknown')
      end if
      sentinel = 0
      write (dump_unit) sentinel, sentinel
      write (dump_unit) self%id, m, l, u, a, b1, b2, b3, fourPiG
    end if
  end subroutine dump_init
END SUBROUTINE sor

!=============================================================================
!> Set scratch array outer boundary to zero
!=============================================================================
subroutine zero_boundaries (patch, f)
  class(patch_t):: patch
  real(kind=realkind):: f(:,:,:)
  class(patch_t), pointer:: m(:)
  !...........................................................................
  if (patch%is_periodic) return
  m => patch
  if (size(f,1) > 1) then
    f(m%lb(1):m%lo(1),:,:) = 0.0_realkind
    f(m%uo(1):m%ub(1),:,:) = 0.0_realkind
  end if
  if (size(f,2) > 1) then
    f(:,m%lb(2):m%lo(2),:) = 0.0_realkind
    f(:,m%uo(2):m%ub(2),:) = 0.0_realkind
  end if
  if (size(f,3) > 1) then
    f(:,:,m%lb(3):m%lo(3)) = 0.0_realkind
    f(:,:,m%uo(3):m%ub(3)) = 0.0_realkind
  end if
end subroutine zero_boundaries

!===============================================================================
!> Dump values for analysis
!===============================================================================
SUBROUTINE dump (f, iter)
  real(kind=realkind):: f(:,:,:)
  integer:: iter, sentinel
  !-----------------------------------------------------------------------------
  if (do_dump) then
    sentinel = 1
    write (dump_unit) sentinel, iter
    write (dump_unit) real(f,kind=4)
  end if
END SUBROUTINE dump

END MODULE poisson_mod

MODULE experiment_mod
  USE patch_mod
  USE poisson_mod
  implicit none
  private
  public: initial_conditions
CONTAINS
!===============================================================================
!> Setup initial conditions 
!===============================================================================
SUBROUTINE initial_conditions
  real(kind=realkind) :: amplitude, width, global_center(ndim)

  ! parameters related to patch geometry
  integer :: dims(ndim) = (/1, 1, 1/) ! geometry of grid of patches
  real(kind=8) :: grace=0.1 /         ! parameter deciding amount of extrapolation 
  logical,           save :: first_time = .true.
  integer :: i,j,k
  real(kind=realkind), dimension(ndim) :: patch_size
  real(kind=8) :: position
  class(patch_t), pointer :: p
 
  namelist /cartesian_params/ n, ng, dims,  
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (cartesian_params)
    write (*, cartesian_params)
  end if
  !-----------------------------------------------------------------------------
  ! Allocate patch data and compute initial condition
  !-----------------------------------------------------------------------------
  patch_size = 1.0_realkind / real(dims, kind=realkind)
  do k=1,dims(3)
  do j=1,dims(2)
  do i=1,dims(1)
    allocate(p)
    position = (real((/i,j,k/),kind=8) - 0.5_8) * patch_size ! center of patch
    call p%init (patch_size, position)                       ! set up counters
    call init_patch(p)                                       ! initial condition
    
  enddo
  enddo
  enddo
  !-----------------------------------------------------------------------------
  ! Establish neighbourhood links 
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Update ghostzones 
  !-----------------------------------------------------------------------------
CONTAINS
!
SUBROUTINE initial_conditions(p)
  implicit none
  class(patch_t) :: p
  integer :: i,j,k
  real(kind=8), dimension(ndim) :: pos,center, ids
  !
  if (.not. allocated(p%density)) allocate(p%density(p%gn,p%gn,p%gn))
  ids = 1.0_8 / p%ds                        ! inverse of cell size
  center = p%p - (p%n / 2 + p%ng + 0.5)*ids ! offset to compute cell center position
  do k=p%li(3),p%ui(3) 
  do j=p%li(2),p%ui(2) 
  do i=p%li(1),p%ui(1) 
    pos = (/i, j, k/)*ids + center 
    ! make density distribution a Gaussian
    p%density(i,j,k) =  amplitude / sqrt(2.0_realkind * pi * width) * exp( - sum((pos - global_center)**2) / (2.0_realkind * width**2) )
  enddo
  enddo
  enddo
END SUBROUTINE initial_conditions

END MODULE experiment_mod

!===============================================================================
!> Program that sets up an experiment and solves the Poisson equation 
!===============================================================================
PROGRAM poisson-mockup
  USE experiment_mod
  USE patch_mod
  USE poisson_mod
  implicit none

END

