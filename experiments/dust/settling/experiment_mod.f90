!===============================================================================
!> Experiment setup for experimenting with dust settling, against a background
!> of an isothermal, hydrostatic disk stratification.  The experiment may be
!> linked and run with any SOLVER, and the intent is to keep the dust module
!> separate from the HD/MHD solver, either as an "extras" module, or a separate
!> task.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE mesh_mod
  USE bits_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: initial_conditions
    procedure:: update
    procedure:: boundary_conditions
  end type
  integer:: verbose = 1
  real:: scale_height = 1.0
  logical:: first_time = .true.
CONTAINS

!===============================================================================
!> Read experiment parameters and setup the initial state
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iostat
  real:: z_lower, z_upper
  namelist /experiment_params/ verbose, scale_height
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  ! Read the namelist file with parameters, making sure to onlya do it once.
  ! The OpenMP "critical region" only allows one thread at a time, and the first
  ! one sets first_time = .false., which prevents the nexts ones from reading.
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%output, experiment_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Initialize the patch, allocating %mem and %mesh etc, and then call
  ! initial_conditions to set initial values of the variables
  !-----------------------------------------------------------------------------
  call self%solver_t%init
  call self%initial_conditions
  !-----------------------------------------------------------------------------
  ! Set bits indicating upper and lower boundaries.  The box size is self%box,
  ! and the patch position, size, and cell size are self%position, self%size,
  ! and self%ds, respectively.
  !-----------------------------------------------------------------------------
  if (self%position(3) - self%size(3)/2.0 - self%ds(3) < -self%box(3)/2.0) then
    call self%boundaries%set (bits%zl)
  end if
  if (self%position(3) + self%size(3)/2.0 + self%ds(3) > +self%box(3)/2.0) then
    call self%boundaries%set (bits%zu)
  end if
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    write (io%output,'(a,i4,3(1x,3f7.3),1x,2l2)') &
      'init: id, box, pos, size, lower, upper =', &
      self%id, self%box, self%position, self%size, &
      self%boundaries%is_set (bits%zl), self%boundaries%is_set (bits%zu)
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Setup initial conditions
!===============================================================================
SUBROUTINE initial_conditions (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iz, m(3)
  real:: z, d
  namelist /forces_params/ verbose, scale_height
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%initial_conditions')
  !-----------------------------------------------------------------------------
  ! A patch has one mesh for each direction -- choose the last index as vertical
  !-----------------------------------------------------------------------------
  mesh => self%mesh(3)
  self%gamma = 1d0
  !-----------------------------------------------------------------------------
  ! Loop over the internal part -- the external parts are handled by the frame-
  ! framework
  !-----------------------------------------------------------------------------
  do iz=mesh%lb,mesh%ub
    !---------------------------------------------------------------------------
    ! A mesh has attributes:
    !
    !   mesh%li   : lower internal index
    !   mesh%ui   : upper internal index
    !   mesh%p    : patch mid-point
    !   mesh%d    : cell size
    !   mesh%r(i) : index i distance from midpoint
    !   mesh%h(v) : offset of variable v in grid units
    !
    ! The vertical distance from the mid plane is thus:
    !---------------------------------------------------------------------------
    z = mesh%p + mesh%r(iz)
    !---------------------------------------------------------------------------
    ! An iso-thermal pressure balance, where the pressure is = the density, and
    ! the linear force of gravity is balanced by a pressure gradient has an
    ! analytical solution where the density is a Gaussian.
    !---------------------------------------------------------------------------
    d = exp(-z**2/scale_height)
    !---------------------------------------------------------------------------
    ! self%mem(ix,iy,iz,iv,it,iw) is a 6-dimensional array, where the indices are
    !   ix = x-index
    !   iy = y-index
    !   iz = z-index
    !   iv = variable index (often 1:density, 2-4:velocity/momentum, 5:energy)
    !   it = time slot in circular time slice buffer
    !   iw = work space (often just 1)
    !---------------------------------------------------------------------------
    self%mem(:,:,iz,self%idx%d ,self%it,1) = d
    self%mem(:,:,iz,self%idx%px,self%it,1) = 0.0
    self%mem(:,:,iz,self%idx%py,self%it,1) = 0.0
    self%mem(:,:,iz,self%idx%pz,self%it,1) = 0.0
    self%mem(:,:,iz,self%idx%e ,self%it,1) = d
  end do
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE initial_conditions

!===============================================================================
!> Update the stat, by computing the self%force_per_unit_mass() that the solver
!> needs, and then calling the solver to do the rest
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iz, m(3)
  real:: z
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  !----------------------------------------------------------------------------
  ! Allocate an (nx,ny,nz,3) array and increment the memory count
  !----------------------------------------------------------------------------
  m = self%mesh%gn
  allocate (self%force_per_unit_mass(m(1),m(2),m(3),3))
  call io%bits_mem (storage_size(self%force_per_unit_mass), &
                    product(shape(self%force_per_unit_mass)), 'fpm')
  !----------------------------------------------------------------------------
  ! The expression for the force is linear in z; cf. the intitial condition.
  ! The mesh%h(iv) (in cell units) is the offset for staggered variables.
  !----------------------------------------------------------------------------
  mesh => self%mesh(3)
  do iz=mesh%li,mesh%ui
    z = mesh%p + mesh%r(iz) + mesh%h(self%idx%pz)*mesh%d
    self%force_per_unit_mass(:,:,iz,1) = 0.0
    self%force_per_unit_mass(:,:,iz,2) = 0.0
    self%force_per_unit_mass(:,:,iz,3) = -2.0*z/scale_height
  end do
  !-----------------------------------------------------------------------------
  ! Enforce boundary conditions, conditionally write output,  before calling
  ! the MHD solver
  !-----------------------------------------------------------------------------
  call self%boundary_conditions
  call self%output
  call self%solver_t%update
  !-----------------------------------------------------------------------------
  ! Deallocate the force field, to save memory
  !-----------------------------------------------------------------------------
  call io%bits_mem (-storage_size(self%force_per_unit_mass), &
                    product(shape(self%force_per_unit_mass)), '-fpm')
  deallocate (self%force_per_unit_mass)
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Boundary conditions, for patches that contain a physical boundary.  Note that
!> the code here is anlogous to the one in the initial conditions.
!===============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iz, m(3)
  real:: z, d
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%boundary_conditions')
  if (verbose > 1) then
    write (io%output,*) &
      'experiment_t%boundary_conditions: id, lower, upper =', self%id, &
      self%boundaries%is_set (bits%zl), self%boundaries%is_set (bits%zu)
  end if
  !-----------------------------------------------------------------------------
  ! Lower boundary
  !-----------------------------------------------------------------------------
  mesh => self%mesh(3)
  if (self%boundaries%is_set (bits%zl)) then
    do iz=mesh%lb,mesh%lo
      z = mesh%p + mesh%r(iz)
      d = exp(-z**2/scale_height)
      self%mem(:,:,iz,self%idx%d ,self%it,1) = d
      self%mem(:,:,iz,self%idx%px,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%py,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%pz,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%e ,self%it,1) = d
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Upper boundary
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%zu)) then
    do iz=mesh%uo,mesh%ub
      z = mesh%p + mesh%r(iz)
      d = exp(-z**2/scale_height)
      self%mem(:,:,iz,self%idx%d ,self%it,1) = d
      self%mem(:,:,iz,self%idx%px,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%py,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%pz,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%e ,self%it,1) = d
    end do
  end if
  call trace%end ()
END SUBROUTINE boundary_conditions

END MODULE experiment_mod
