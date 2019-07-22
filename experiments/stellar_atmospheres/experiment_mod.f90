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
  USE nbors_mod
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
    procedure(boundary_condition1), pointer:: &
      boundary_condition => null()
  contains
    procedure:: init
    procedure:: update
  end type
  !.............................................................................
  integer, save:: verbose=1
  logical, save:: smudge_boundary = .false.
  real, save:: d_bot=0.0, e_bot=0.0, ee_top=5.2, p_top=0.5, tt_top=4000.
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
  namelist /experiment_params/ d_bot, e_bot, bottom_bc_rate, ee_top, &
                               tt_top, verbose, smudge_boundary
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
      !d_bot = d_bot/scaling%d
      !e_bot = e_bot/scaling%e
    end if
    !$omp end critical (input_cr)
  end if
  if (.not. io%master) verbose=0
  !-----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !-----------------------------------------------------------------------------
  self%mhd = .true.
  call self%solver_t%init
  self%unsigned = .false.
  self%boundary_condition => boundary_condition3
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
!> Save a pointer to the task list, for use in underlings that need it
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(experiment_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  self%task_list => task_list
  call nbors%init (self)
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
!===============================================================================
SUBROUTINE boundary_condition1 (self)
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
  tt => self%mem(:,:,:,self%idx%tt,self%it,1)
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
      f_scale = (tt(ix,iy,li(3))/tt_top)**0.05
      h_scale = tt(ix,iy,li(3))/(scaling%temp*self%gravity%const_grav)*f_scale
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
    call bottom_condition1 (self)
  end if
  call trace_end
END SUBROUTINE boundary_condition1

!===============================================================================
!===============================================================================
SUBROUTINE boundary_condition2 (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, tt, ss
  real, dimension(:,:,:), allocatable:: ux, uy, uz, ex, ey, ez
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz
  real(8):: d_ratio, h_scale, f_scale, dlnddz, dlnedz, grav, mean1, mean2
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
  tt => self%mem(:,:,:,self%idx%tt,self%it,1)
  grav = self%gravity%const_grav
   if (self%is_set(bits%zl+bits%zu)) print *, 'doing BC according to bits'
  if (m3%lower_boundary) then
    !---------------------------------------------------------------------------
    ! Clamp the energy per unit mass at the table limit
    !---------------------------------------------------------------------------
    e = max(e,3.5*d)
    call eos%lookup_table (shape(d), d=d, e=e, tt=tt)
  allocate(ss,mold=d)
  !ss=e/d
  !call self%e2ss(ss,self%it)
  !call self%boundaries%log_extrapolate(ss,self%idx%e)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! -- adjust the scale height so it becomes larger at high temperatures
    if (smudge_boundary) then
      mean1 = sum(tt(li(1):ui(1),li(2):ui(2),li(3)+1))/(self%n(1)*self%n(2))
      mean2 = sum(tt(li(1):ui(1),li(2):ui(2),li(3)))/(self%n(1)*self%n(2))
      f_scale = (mean2/mean1)**0.05
      h_scale = mean2/grav*f_scale
    end if
    do iz=lb(3),lo(3)
    do iy=lb(2),ub(2)
    do ix=lb(1),ub(1)
      ! -- adjust the scale height so it becomes larger at high temperatures
      if (.not.smudge_boundary) then
        f_scale = (tt(ix,iy,li(3))/tt_top)**0.05
        h_scale = tt(ix,iy,li(3))/grav*f_scale
        ! -- option: more localized, not constant clamping:
        !f_scale = (tt(ix,iy,li(3)+1)/tt(ix,iy,li(3)))**0.05
        !h_scale = tt(ix,iy,li(3))/grav*f_scale
      end if
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
    !e=ss*d
    !call self%ss2e(ss,self%it)
    !deallocate(ss)
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
    ! -------------------------------------------------
    ! avoid going outside the EOS table.
    ! -------------------------------------------------
    !where (e < exp(eos%lne_min)) e = exp(eos%lne_min)
  else if (m3%upper_boundary) then
    call bottom_condition2 (self)
  end if
  call trace_end
END SUBROUTINE boundary_condition2

!===============================================================================
!===============================================================================
SUBROUTINE boundary_condition3 (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz, ss
  real, dimension(:,:,:), allocatable:: ux, uy, uz, ex, ey, ez
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz
  real:: d_ratio, h_scale, f, p, ee, tt
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
    ! Estimate scale height
    !---------------------------------------------------------------------------
    do iz=lo(3),lb(3),-1
      ! -- weight from 0 to 1 in inteval lb to li
      p = (z(iz)-z(lb(3)))/(z(li(3))-z(lb(3)))
      ! -- remap to interval p_top - 1.0
      p = p_top*(1.-p) + 1.*p
      do iy=lb(2),ub(2)
      do ix=lb(1),ub(1)
        if (associated(self%tt)) then
          tt = self%tt(ix,iy,li(3))
        else
          tt = tt_top
        end if
        f = (tt/tt_top)**0.03
        if (iz==lo(3)) then
          h_scale =     tt/(scaling%temp*self%gravity%const_grav)
        else
          h_scale = tt_top/(scaling%temp*self%gravity%const_grav)
        end if
        d_ratio = exp((z(iz)-z(iz+1))/h_scale)
        ee = e(ix,iy,li(3))/d(ix,iy,li(3))
        d (ix,iy,iz) = d (ix,iy,iz+1)*d_ratio
        e (ix,iy,iz) = d (ix,iy,iz  )*ee
        px(ix,iy,iz) = px(ix,iy,iz+1)*d_ratio*p
        py(ix,iy,iz) = py(ix,iy,iz+1)*d_ratio*p
        pz(ix,iy,iz) = pz(ix,iy,iz+1)*d_ratio*P
      end do
      end do
    end do
    if (self%mhd) then
      do iz=lb(3),lo(3)-1
        bx(:,:,iz) = bx(:,:,li(3)-1)
        by(:,:,iz) = by(:,:,li(3)-1)
        bz(:,:,iz) = bz(:,:,li(3)-1)
      end do
    end if
  else if (m3%upper_boundary) then
    call bottom_condition2 (self)
  end if
  call trace_end
END SUBROUTINE boundary_condition3

!===============================================================================
!===============================================================================
SUBROUTINE bottom_condition1 (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz
  real(8):: d_ratio, h_scale, dlnddz, dlnedz, taper
  logical, save:: first_lower=.true., first_upper=.true.
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%bottom_condition')
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
  call trace_end
END SUBROUTINE bottom_condition1

!===============================================================================
!===============================================================================
SUBROUTINE bottom_condition2 (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, bx, by, bz
  real(8), dimension(:), pointer:: z
  class(mesh_t), pointer:: m3
  integer:: ix, iy, iz
  real(8):: d_ratio, h_scale, dlnddz, dlnedz, taper
  logical, save:: first_lower=.true., first_upper=.true.
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%bottom_condition')
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
  call trace_end
END SUBROUTINE bottom_condition2
 
END MODULE
