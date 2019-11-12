!===============================================================================
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
!===============================================================================
MODULE experiment_mod
  USE rt_nbors_mod
  USE solver_mod
  USE list_mod
  USE download_mod
  USE patch_mod
  USE link_mod
  USE rt_mod
  USE mesh_mod
  USE task_mod
  USE eos_mod
  USE data_io_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE scaling_mod
  USE bits_mod
  USE omp_timer_mod
  implicit none
  private
  !.............................................................................
  type, public, extends(solver_t):: experiment_t
    real, dimension(:,:,:,:), pointer:: lower => null()
    real, dimension(:,:,:,:), pointer:: upper => null()
  contains
    procedure:: init
    procedure:: init_task_list
    procedure:: update
    procedure:: boundary_condition
  end type
  !.............................................................................
  integer, save:: verbose=1
  logical, save:: smudge_boundary = .false.
  real, save:: d_bot=0.0, e_bot=0.0
  real, save:: bottom_bc_rate=0.01
  type(experiment_t), public:: experiment
CONTAINS

!===============================================================================
!> Experiment setup
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /experiment_params/ d_bot, e_bot, bottom_bc_rate, verbose
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (io%master) write (*, experiment_params)
    !d_bot = d_bot/scaling%d
    !e_bot = e_bot/scaling%e
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  if (.not. io%master) verbose=0
  !-----------------------------------------------------------------------------
  ! Initialize the MHD task and the RT sub-tasks, and setup their nbor lists
  !-----------------------------------------------------------------------------
  self%mhd = .true.
  call self%solver_t%init
  self%unsigned = .false.
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  self%periodic(3) = .false.                     ! non-periodic vertically
  call self%init_bdries
  call self%boundaries%init(self%mesh)
  self%boundaries%id = self%id
  if (verbose > 0) then
    if (self%mesh(3)%lower_boundary) write(io%output,*) self%id, 'at lower B'
    if (self%mesh(3)%upper_boundary) write(io%output,*) self%id, 'at upper B'
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Intercept the init_task_list call to finalize RT setup
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(experiment_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('experiment_t%init_task_list')
  ! -- punt to the next level; picked up by extras_t ---------------------------
  call self%solver_t%init_task_list (task_list)
  ! -- for RT tasks, modify the nbor lists and add tasks to the task list-------
  if (associated (self%rt)) then
    if (self%rt%on) then
      call rt_nbors%init (self)
      call self%rt%add_to_task_list
    end if
  end if
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Experiment update
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  if (self%mesh(3)%lower_boundary .or. &
      self%mesh(3)%upper_boundary) then
      call self%boundary_condition
  end if
  call self%solver_t%update
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Top boundary conditions
!===============================================================================
SUBROUTINE boundary_condition (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss, d1, e1, p1
  real, dimension(:,:,:), allocatable:: ux, uy, uz, ex, ey, ez, pg
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz, iz1
  real:: d_ratio, e_ratio, h_scale, p, ee, dlnddz, dlnedz
  logical, save:: first_lower=.true., first_upper=.true.
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_condition')
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
  if (self%is_set(bits%zl+bits%zu)) print *, 'doing BC according to bits'
  if (m3%lower_boundary) then
    !---------------------------------------------------------------------------
    ! Clamp the energy per unit mass at the table limit
    !---------------------------------------------------------------------------
    e = max(e,3.5*d)
    !---------------------------------------------------------------------------
    ! Boundary conditions need the gas pressure in one layer
    !---------------------------------------------------------------------------
    allocate (pg(size(d,1),size(d,2),1))
    d1 => d(:,:,li(3):li(3))
    e1 => e(:,:,li(3):li(3))
    call eos%lookup_table (shape(pg), pg=pg, d=d1, e=e1)
    !---------------------------------------------------------------------------
    ! Note that the loop for simplicity includes iz=li(3), but d_ratio=1.0 for
    ! that index, except for the vertical momentum
    !---------------------------------------------------------------------------
    do iz=lb(3),li(3)
      do iy=lb(2),ub(2)
      do ix=lb(1),ub(1)
        !-----------------------------------------------------------------------
        ! 1) The BC for energy per unit mass is that it is = the value at li,
        ! while the density is assumed to drop with a scale height that is to
        ! lowest order proportional to temperature.
        !-----------------------------------------------------------------------
        iz1 = li(3)
        h_scale = pg(ix,iy,1)/(d(ix,iy,iz1)*self%gravity%const_grav)
        ee = e(ix,iy,iz1)/d(ix,iy,iz1)
        d_ratio = exp((z(iz)-z(iz1))/h_scale)
        d(ix,iy,iz) = d(ix,iy,iz1)*d_ratio
        e(ix,iy,iz) = d(ix,iy,iz)*ee
        !-----------------------------------------------------------------------
        ! 2) The horizontal components of the momenta are assumed to taper off
        ! with density, consistent with assuming a no-shear boundary condition
        ! for the horizontal velocity components
        !-----------------------------------------------------------------------
        px(ix,iy,iz) = px(ix,iy,iz1)*d_ratio
        py(ix,iy,iz) = py(ix,iy,iz1)*d_ratio
        !-----------------------------------------------------------------------
        ! 3) The vertical momenta for z-index 1-li are locked to the vertical
        ! momentum at the point li+1, which is the 1st point inside the physical
        ! domain.  The vertical velocity is tapered of with the d_ratio factor.
        !-----------------------------------------------------------------------
        iz1 = li(3)+1
        d_ratio = exp((z(iz)-z(iz1))/h_scale)
        pz(ix,iy,iz) = pz(ix,iy,iz1)*d_ratio**2
      end do
      end do
    end do
    deallocate (pg)
    if (self%mhd) then
      do iz=lb(3),lo(3)-1
        bx(:,:,iz) = bx(:,:,li(3)-1)
        by(:,:,iz) = by(:,:,li(3)-1)
        bz(:,:,iz) = bz(:,:,li(3)-1)
      end do
    end if
  else if (m3%upper_boundary) then
    !---------------------------------------------------------------------------
    ! Compute the pressure pg at the bottom BC, and what we would like it to be
    !---------------------------------------------------------------------------
    d1 => d(:,:,ui(3):ui(3))
    e1 => e(:,:,ui(3):ui(3))
    allocate (pg(size(d,1),size(d,2),1))
    call eos%lookup_table (shape(pg), pg=pg, d=d1, e=e1)
    allocate (d1(size(d,1),size(d,2),1))
    allocate (e1(size(d,1),size(d,2),1))
    allocate (p1(size(d,1),size(d,2),1))
    d1 = d_bot
    e1 = e_bot
    call eos%lookup_table (shape(p1), pg=p1, d=d1, e=e1)
    !---------------------------------------------------------------------------
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      dlnddz = log(d(ix,iy,ui(3)-1)/d(ix,iy,ui(3)-3))/(z(ui(3)-1)-z(ui(3)-3))
      dlnedz = log(e(ix,iy,ui(3)-1)/e(ix,iy,ui(3)-3))/(z(ui(3)-1)-z(ui(3)-3))
      !-------------------------------------------------------------------------
      ! We assume that the iz index of the last free values is ui(3)-1, except
      ! for pz, for which it us ui(3).  Hence, if there are 2 guard zones, we
      ! set 2 values for pz, and 3 values for the other variables
      !-------------------------------------------------------------------------
      do iz=ui(3),ub(3)
        !-----------------------------------------------------------------------
        ! For inflow, set d,e to d_bot,e_bot at li, and set the momenta constant,
        ! corresponding to velocities that taper off with depth
        !-----------------------------------------------------------------------
        if (pz(ix,iy,ui(3)-1) < 0.0) then
          d(ix,iy,iz) = bottom_bc_rate*d_bot + (1.0-bottom_bc_rate)*d(ix,iy,iz)
          e(ix,iy,iz) = bottom_bc_rate*e_bot + (1.0-bottom_bc_rate)*e(ix,iy,iz)
          d(ix,iy,iz) = d_bot*exp((z(iz)-z(ui(3)))*dlnddz)
          e(ix,iy,iz) = e_bot*exp((z(iz)-z(ui(3)))*dlnedz)
          px(ix,iy,iz) = px(ix,iy,ui(3)-1)
          py(ix,iy,iz) = py(ix,iy,ui(3)-1)
          pz(ix,iy,iz) = pz(ix,iy,ui(3)-1)
        !-----------------------------------------------------------------------
        ! For outflow, extrapolate d and e, and assume velocity is constant.
        ! Leave values untouched at iz=ui(3), to affect outflow as little as
        ! possible
        !-----------------------------------------------------------------------
        else
          !---------------------------------------------------------------------
          ! At iz=ui(3), nudge the density towards what would give a constant
          ! pressure = pg(d_bot,e_bot), while keeping ee=e/d unchanged
          !---------------------------------------------------------------------
          if (iz == ui(3)) then
            d1(ix,iy,1) = d(iz,iy,iz)*(p1(ix,iy,1)/pg(ix,iy,1))
            e1(ix,iy,1) = e(iz,iy,iz)*(p1(ix,iy,1)/pg(ix,iy,1))
            associate (w => bottom_bc_rate)
            d(ix,iy,iz) = w*d1(ix,iy,1) + (1.0-w)*d(ix,iy,iz)
            e(ix,iy,iz) = w*e1(ix,iy,1) + (1.0-w)*e(ix,iy,iz)
            end associate
          !---------------------------------------------------------------------
          ! Note that pz(ix,iy,iz) is centered half a cell above pressure, so
          ! by imposing pz at ui(3)+1=uo(3), we are in principle imposing the
          ! density evolution at uo(3), but that is overruled in any case here,
          ! hence it is not necessary to stagger the BC for pz
          !---------------------------------------------------------------------
          else
            d_ratio = exp((z(iz)-z(ui(3)))*dlnddz)
            e_ratio = exp((z(iz)-z(ui(3)))*dlnedz)
            d(ix,iy,iz)  =  d(ix,iy,ui(3))*d_ratio
            e(ix,iy,iz)  =  e(ix,iy,ui(3))*e_ratio
            px(ix,iy,iz) = px(ix,iy,ui(3))*d_ratio
            py(ix,iy,iz) = py(ix,iy,ui(3))*d_ratio
            pz(ix,iy,iz) = pz(ix,iy,ui(3))*d_ratio
          end if
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
    deallocate (d1, e1, p1, pg)
  end if
  call trace_end
END SUBROUTINE boundary_condition
 
END MODULE
