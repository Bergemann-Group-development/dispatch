!> =============================================================================
!> This module contains all experiment specific information.  The experiment_t
!> object extends the generic MHD object mhd_t, which is again an extension of
!> the basic patch_t object for tasks (object task_t) with meshes.
!>
!> The MHD object is chosen with the SOLVER macro in the Makefile in the curren
!> directory, while the generic task_t and patch_t objects are defined in the
!> $(TOP)/tasks/ directory, in task_mod.f90 and patch_mod.f90.
!>
!> The init procedure here (experiment_t::init) is called from the cartesian_mod
!> module, which is specified as the patch distributor in the main program
!> dispatch.f90.   The update procedure here is called by the task_list_mod
!> update procedure in $(TOP)/lists/task_list_mod.f90.
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE eos_mod
  USE trace_mod
  USE omp_timer_mod
  USE task_mod
  USE gpatch_mod
  USE link_mod
  USE rt_mod
  USE rt_solver_mod
  USE solver_mod
  USE scaling_mod
  USE mesh_mod
  USE download_mod
  USE bits_mod
  USE spitzer_mod
  USE scratch_mod
  implicit none
  private
  type, public, extends(rt_solver_t):: experiment_t
    real, dimension(:,:,:,:), pointer:: lower => null()
    real, dimension(:,:,:,:), pointer:: upper => null()
    type(spitzer_t):: spitzer
    logical :: first_bc = .true.
  contains
    procedure:: init
    procedure:: update
    procedure:: dnload
    procedure:: boundary_condition
    procedure, nopass:: is_ready
  end type
  integer, save:: verbose=1
  real, save:: grav=2.75, d_bot=0.0, e_bot=0.0, top_temperature=4000.
  real, save:: bottom_bc_rate=0.01
  type(experiment_t), public:: experiment
CONTAINS

!> =============================================================================
!> Experiment setup
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /experiment_params/ grav, d_bot, e_bot, bottom_bc_rate, &
                               top_temperature, verbose
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, experiment_params, iostat=iostat)
      if (io%master) write (*, experiment_params)
      top_temperature = top_temperature / scaling%temp
    end if
    !$omp end critical (input_cr)
  end if
  if (.not. io%master) verbose=0
  !-----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !-----------------------------------------------------------------------------
  call self%rt_solver_t%init
  call self%spitzer%init(self%mesh, self%idx)
  scratch%force_per_unit_mass => force
  self%grav = grav
  self%unsigned = .false.
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  self%periodic(3) = .false.                     ! non-periodic vertically
  call self%init_bdries
  if (verbose > 0) then
    if (self%mesh(3)%lower_boundary) write(io%output,*) self%id, 'at lower B'
    if (self%mesh(3)%upper_boundary) write(io%output,*) self%id, 'at upper B'
  end if
  call trace%end()
END SUBROUTINE init

!> =============================================================================
!> Experiment update
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  ! -- apply BCs only to the parent MHD task
  if (self%mesh(3)%lower_boundary .or. &
      self%mesh(3)%upper_boundary) then
    if (self%i_omega==0) then
      call self%boundary_condition
    end if
  end if
  call self%rt_solver_t%update
  call self%spitzer%update(self%mem,self%it,self%dtime)
  call trace%end()
END SUBROUTINE update

!> =============================================================================
!> =============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  if (self%i_omega > 0) then
    call download%download_link (self%link, only=only, all_times=.true.)
  else
    call download%download_link (self%link, only=only)
  end if
END SUBROUTINE dnload

!> =============================================================================
!> =============================================================================
SUBROUTINE boundary_condition (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, tt
  real, dimension(:,:,:), allocatable:: ux, uy, uz, ex, ey, ez
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz
  real:: d_ratio, h_scale, f_scale, dlnddz, dlnedz
  logical, save:: first_lower=.true., first_upper=.true.
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%boundary_condition')
  call self%check_density(' before BC ')
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  m3 => self%mesh(3)
  z  => m3%r
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  tt => self%mem(:,:,:,self%idx%tt,self%it,1)
  ! --------------------------------------------------------------------
  ! Lookup temperature
  ! --------------------------------------------------------------------
  if (self%first_bc) then
    call eos%lookup_table (shape(d), d=d, e=e,tt=tt)
    self%first_bc = .false.
  end if
  if (m3%lower_boundary) then
    !---------------------------------------------------------------------------
    ! Clamp the energy per unit mass at the table limit
    !---------------------------------------------------------------------------
    e = max(e,3.5*d)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    do iz=lb(3),lo(3)
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      ! -- adjust the scale height so it becomes larger at high temperatures
      f_scale = (tt(ix,iy,li(3))/top_temperature)**0.05
      h_scale = tt(ix,iy,li(3))/(scaling%temp*grav)*f_scale
      d_ratio = exp((z(iz)-z(li(3)))/h_scale)
      if (first_lower) then
        first_lower = .false.
        write(io%output,*) self%id, 'experiment_t%boundary_condition (lower)', &
          tt(ix,iy,li(3)), scaling%temp, h_scale, d_ratio
      end if
      ! -- velocity / momentum boundary conditions
      d (ix,iy,iz) = d (ix,iy,li(3))*d_ratio
      e (ix,iy,iz) = e (ix,iy,li(3))*d_ratio
      e (ix,iy,iz) = e (ix,iy,iz)/f_scale
      px(ix,iy,iz) = px(ix,iy,li(3))*d_ratio
      py(ix,iy,iz) = py(ix,iy,li(3))*d_ratio
      pz(ix,iy,iz) = pz(ix,iy,li(3))*d_ratio
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Let the magnetic field take care of itself, regularizing only the layer
    ! that gets contaminated by the vertical derivatives.  The time derivative
    ! of the horizontal components involve a vertical derivative, and the horizontal
    ! components of the electric field controlling the vertical component of the
    ! magnetic field depend on vertical velocities coming out of vertical staggering
    !---------------------------------------------------------------------------
    if (self%mhd) then
      do iz=lb(3),lo(3)-1
        bx(:,:,iz) = bx(:,:,li(3)-1)
        by(:,:,iz) = by(:,:,li(3)-1)
        bz(:,:,iz) = bz(:,:,li(3)-1)
      end do
    end if
  else if (m3%upper_boundary) then
    do iz=uo(3),ub(3)
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      if (pz(ix,iy,ui(3)) < 0.0) then
        d (ix,iy,ui(3)) = d_bot*bottom_bc_rate + d(ix,iy,ui(3))*(1.-bottom_bc_rate)
        e (ix,iy,ui(3)) = e_bot*bottom_bc_rate + e(ix,iy,ui(3))*(1.-bottom_bc_rate)
        px(ix,iy,ui(3)) = px(ix,iy,ui(3))*(1.-bottom_bc_rate)
        py(ix,iy,ui(3)) = py(ix,iy,ui(3))*(1.-bottom_bc_rate)
      end if
      dlnddz = log(d(ix,iy,ui(3))/d(ix,iy,ui(3)-2))/(z(ui(3))-z(ui(3)-2))
      dlnedz = log(e(ix,iy,ui(3))/e(ix,iy,ui(3)-2))/(z(ui(3))-z(ui(3)-2))
      d_ratio = exp((z(iz)-z(ui(3)))*dlnddz)
      d (ix,iy,iz) = d (ix,iy,ui(3))*d_ratio
      px(ix,iy,iz) = px(ix,iy,ui(3))*d_ratio
      py(ix,iy,iz) = py(ix,iy,ui(3))*d_ratio
      pz(ix,iy,iz) = pz(ix,iy,ui(3))*d_ratio
      e(ix,iy,iz) = e(ix,iy,ui(3))*exp((z(iz)-z(ui(3)))*dlnedz)
      if (first_upper) then
        first_upper = .false.
        write(io%output,*) self%id, 'experiment_t%boundary_condition (upper)', &
          dlnddz, dlnedz, z(ui(3)-1), z(ui(3)), d_ratio
      end if
    end do
    end do
    end do
    if (self%mhd) then
      do iz=uo(3)+1,ub(3)
        bx(:,:,iz) = bx(:,:,uo(3))
        by(:,:,iz) = by(:,:,uo(3))
        bz(:,:,iz) = bz(:,:,uo(3))
      end do
    end if
  end if
  call self%check_density(' after BC ')
  call trace%end()
END SUBROUTINE boundary_condition

!> =============================================================================
!> Stellar atmosphere with gravity pointing in the positive y-direction
!> =============================================================================
FUNCTION force (self) RESULT (ff)
  class(gpatch_t)   :: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3) :: ff
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%force')
  ff(:,:,:,1) = 0.0
  ff(:,:,:,2) = 0.0
  ff(:,:,:,3) = grav
  call trace%end()
END FUNCTION force

!> =============================================================================
!> Check if a task is ready to update, by changing if all nbor task are "ahead
!> of" the task in question
!> =============================================================================
LOGICAL FUNCTION is_ready (link)
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor
  class(rt_t), pointer:: target, source
  integer:: state
  logical:: debug1, debug2
  !-----------------------------------------------------------------------------
  target => rt%task2rt(link%task)
  !!$omp atomic
  !n_check = n_check+1
  !$omp atomic read
  state = link%task%state
  debug2 = .false.
  if (target%is_set (bits%virtual)) then
    is_ready = .false.
  else
    is_ready = .true.
    nbor => link%nbor
    do while (associated(nbor))
     source => rt%task2rt(nbor%task)
     debug1 = source%id==io%id_debug
     debug2 = debug2 .or. debug1
     if (debug1) write(io_unit%log,*) 'dbg check_ready: id, nbor, needed, ahead', &
       link%task%id, source%id, nbor%needed, source%is_ahead_of(target), &
       source%is_set (bits%virtual)
     if (nbor%needed) then
        !$omp atomic read
        state = source%state
        if (state==3) then
          is_ready = .false.
        else
          is_ready = is_ready .and. source%is_ahead_of (target)
        end if
        !if (source%id == io%id_debug) then
         ! write(io%output,*) '    dispatcher3_t%is_ready: id, nbid, ready =', &
          !  target%id, source%id, is_ready
        !end if
        if (.not.is_ready) then
          if (verbose > 2) &
            write(io_unit%mpi,*) &
              target%id, 'is_ready: failed  on', source%id, &
              source%time, target%time
          exit
        end if
      end if
      nbor => nbor%next
    end do
  end if
  if (is_ready) then
    !!$omp atomic
    !n_ready = n_ready+1
    if (target%id == io%id_debug) &
      write(io%output,*) 'DBG dispatcher3_t%is_ready: id =', target%id
    if (debug2) &
      write(io_unit%log,*) 'DBG dispatcher3_t%is_ready: id =', target%id
!    write(io_unit%mpi,'(f12.6,2x,a,i6,f12.6)') &
!      wallclock(), 'is_ready:', target%id, target%time
  end if
contains
  !-----------------------------------------------------------------------------
  function task2rt(task) result(rt)
    class(task_t), pointer:: task
    class(rt_t), pointer:: rt
    !...........................................................................
    select type (task)
    class is (rt_t)
    rt => task
    end select
  end function task2rt
  !-----------------------------------------------------------------------------
END FUNCTION is_ready

END MODULE
