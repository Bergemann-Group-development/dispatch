!===============================================================================
!> Experiment setup for experimenting with dust settling, against a background
!> of an isothermal, hydrostatic disk stratification.  The experiment may be
!> linked and run with any SOLVER, and the intent is to keep the dust module
!> separate from the HD/MHD solver, either as an "extras" module, or a separate
!> task.
!>
!> Calling sequence (shows only data type and procedure):
!>   experiment_t%update
!>     mhd_t%update
!>     extras_t%post_update
!>       particle_solver_t%update
!>       while (particle):
!>         particle_solver%kdk (particle)
!>           dust_t%drag (particle)
!>           particle_solver_t%force (particle)
!>             experiment_t%pforce (position)
!>         if (particle_solver_t%is_outside):
!>           particle_solver_t%particles%remove
!>           particle_solver_t%exports%append
!>       dust_t%gather
!>       dust_t%renormalize
!>       dust_t%evolve
!        dust_t%scatter
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE mesh_mod
  USE bits_mod
  USE task_mod
  USE link_mod
  USE list_mod
  USE particles_mod
  USE scaling_mod
  USE units_mod
  USE omp_lock_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    real:: dz
  contains
    procedure:: init
    procedure:: initial_conditions
    procedure:: init_task_list
    procedure:: update
    procedure:: boundary_conditions
  end type
  integer:: verbose = 1
CONTAINS

!===============================================================================
!> Read experiment parameters and setup the initial state
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iostat, l(3), u(3)
  real:: z_lower, z_upper
  namelist /experiment_params/ verbose
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  ! Read the namelist file with parameters, making sure to onlya do it once.
  ! The OpenMP "critical region" only allows one thread at a time, and the first
  ! one sets first_time = .false., which prevents the nexts ones from reading.
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%output, experiment_params)
    first_time = .false.
  end if
  call io%lock%unset
  call scaling%init
  self%gamma = scaling%gamma
  self%csound = 1.0
  !-----------------------------------------------------------------------------
  ! Initialize the patch, allocating %mem and %mesh etc
  !-----------------------------------------------------------------------------
  call self%solver_t%init
  call self%initial_conditions
  call self%particle_solver%initial_condition
  self%particle_solver%pforce => pforce
  self%periodic = [.true., .true., .false.]
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
  real:: z, d, g2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%initial_conditions', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! A patch has one mesh for each direction -- choose the last index as vertical
  !-----------------------------------------------------------------------------
  mesh => self%mesh(3)
  !-----------------------------------------------------------------------------
  ! Gamma factors for relating thermal energy to temperature and pressure
  !-----------------------------------------------------------------------------
  if (self%gamma == 1d0) then
    g2 = 1.0
  else
    g2 = self%gamma*(self%gamma-1.0)
  end if
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
    d = exp(-0.5*z**2)
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
    self%mem(:,:,iz,self%idx%e ,self%it,1) = d*g2*scaling%u**2
  end do
  !-----------------------------------------------------------------------------
  ! Make sure force_per_unit_mass is allocated
  !-----------------------------------------------------------------------------
  if (.not.associated(self%force_per_unit_mass)) then
    m = self%mesh%gn
    allocate (self%force_per_unit_mass(m(1),m(2),m(3),3))
    self%force_per_unit_mass = 0.0
  end if
  call trace%end (itimer)
END SUBROUTINE initial_conditions

!===============================================================================
!> Initialize particle solver nbor list, after task list is created
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(experiment_t):: self
  class(list_t), pointer:: task_list
  class(nbor_t), pointer:: nbor, nbor2
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init_task_list')
  call self%particle_solver%init_task_list (task_list)                  ! part
  call self%gpatch_t%init_task_list (task_list)                         ! keep
  nbor => self%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (experiment_t)
      !-------------------------------------------------------------------------
      ! Allocate an nbor, point to the particle_solver_t part of it, and prepend
      !-------------------------------------------------------------------------
      allocate (nbor2)
      nbor2%task => task%particle_solver
      nbor2%next => self%particle_solver%nbor
      self%particle_solver%nbor => nbor2
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update the stat, by computing the self%force_per_unit_mass() that the solver
!> needs, and then calling the solver to do the rest
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iz
  real:: z
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  if (self%id == 1 .and. self%particle_solver%verbose > 0) &
    call total_dust_mass (self)
  !----------------------------------------------------------------------------
  ! Allocate an (nx,ny,nz,3) array and increment the memory count
  !----------------------------------------------------------------------------
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
    self%force_per_unit_mass(:,:,iz,3) = -z
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
  call trace%end(itimer)
END SUBROUTINE update

!===============================================================================
!> Boundary conditions, for patches that contain a physical boundary.  Note that
!> the code here is anlogous to the one in the initial conditions.
!===============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  !.............................................................................
  class(mesh_t), pointer:: mesh
  integer:: iz, m(3), l(3), u(3)
  real:: z, d
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%boundary_conditions')
  if (verbose > 2) then
    write (io%output,*) &
      'experiment_t%boundary_conditions: id, lower, upper =', self%id, &
      self%particle_solver%lower(3), self%particle_solver%upper(3)
  end if
  !-----------------------------------------------------------------------------
  ! Lower boundary
  !-----------------------------------------------------------------------------
  mesh => self%mesh(3)
  if (self%particle_solver%lower(3)) then
    do iz=mesh%lb,mesh%lo
      z = mesh%p + mesh%r(iz)
      d = exp(-0.5*z**2)
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
  if (self%particle_solver%upper(3)) then
    do iz=mesh%uo,mesh%ub
      z = mesh%p + mesh%r(iz)
      d = exp(-0.5*z**2)
      self%mem(:,:,iz,self%idx%d ,self%it,1) = d
      self%mem(:,:,iz,self%idx%px,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%py,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%pz,self%it,1) = 0.0
      self%mem(:,:,iz,self%idx%e ,self%it,1) = d
    end do
  end if
  call trace%end ()
END SUBROUTINE boundary_conditions

!===============================================================================
!> Point force, assuming lengh is measured in vertical scale heights
!===============================================================================
FUNCTION pforce (r) RESULT (f)
  real(8):: r(3), f(3), z
  !-----------------------------------------------------------------------------
  f(1) = 0.0
  f(2) = 0.0
  f(3) = -r(3)
END FUNCTION pforce

!===============================================================================
!> For validation purposes, compute the total dust mass per bin
!===============================================================================
SUBROUTINE total_dust_mass (self)
  class(experiment_t):: self
  !.............................................................................
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  class(particle_t), pointer:: p
  real(8):: mass(self%particle_solver%dust%n_bin)
  integer(8):: n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%total_dust_mass', itimer=itimer)
  n = 0
  mass = 0d0
  link => self%task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
      call task%particle_solver%particles%lock%set ('total_dust_mass')
      p => task%particle_solver%particles%head
      do while (associated(p))
        mass(p%bin) = mass(p%bin) + p%w
        n = n+1
        p => p%next
      end do
      call task%particle_solver%particles%lock%unset ('total_dust_mass')
    end select
    link => link%next
  end do
  !.............................................................................
  if (self%particle_solver%dust%verbose > 0) then
    write(io_unit%log,'(a,f10.3,i12,2x,1p,30e11.2)') &
      'total_dust_mass:', self%time, n, mass
    if (self%particle_solver%dust%verbose > 1) then
      write(io%output,'(a,f10.3,i10,2x,1p,30e11.2)') &
        'total_dust_mass:', self%time, n, mass
    end if
  end if
  !.............................................................................
  if (any(abs(mass) > 1e20)) then
    link => self%task_list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (experiment_t)
        call task%particle_solver%particles%lock%set ('total_dust_mass')
        p => task%particle_solver%particles%head
        do while (associated(p))
          if (abs(p%w) > 1e20) then
            write(io_unit%log,'("WARNING:",i6,f10.3,4i4,1p,e11.2)') &
              task%id, task%time, p%q, p%bin, p%w
          end if
          p => p%next
        end do
        call task%particle_solver%particles%lock%unset ('total_dust_mass')
      end select
      link => link%next
    end do
  end if
  call trace%end (itimer)
END SUBROUTINE total_dust_mass

END MODULE experiment_mod
