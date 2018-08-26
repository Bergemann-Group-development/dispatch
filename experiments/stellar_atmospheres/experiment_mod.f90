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
  USE trace_mod
  USE gpatch_mod
  USE link_mod
  USE rt_solver_mod
  USE solver_mod
  USE scaling_mod
  USE mesh_mod
  implicit none
  private
  type, public, extends(rt_solver_t):: experiment_t
    real, dimension(:,:,:,:), pointer:: lower => null()
    real, dimension(:,:,:,:), pointer:: upper => null()
  contains
    procedure:: init
    procedure:: update
    procedure:: boundary_condition
  end type
  integer, save:: verbose=1
  real, save:: grav=2.75, d_bot=0.0, e_bot=0.0, top_temperature=4000.
  real, save:: bottom_bc_rate=0.01
CONTAINS

!> =============================================================================
!> Experiment setup
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  namelist /experiment_params/ grav, d_bot, e_bot, bottom_bc_rate, &
                               top_temperature, verbose
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, experiment_params, iostat=iostat)
      if (io%master) write (*, experiment_params)
    end if
    !$omp end critical (input_cr)
  end if
  if (.not. io%master) verbose=0
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call self%rt_solver_t%init
  self%force_per_unit_mass => force
  self%grav = grav
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  self%periodic(3) = .false.                     ! non-periodic vertically
  call self%init_bdries
  if (verbose > 0) then
    if (self%mesh(3)%lower_boundary) print *, self%id, 'at lower B'
    if (self%mesh(3)%upper_boundary) print *, self%id, 'at upper B'
  end if
  call trace_end (itimer=itimer)
END SUBROUTINE init

!> =============================================================================
!> Experiment update
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  ! -- apply BCs only to the parent MHD task
  if (self%mesh(3)%lower_boundary .or. &
      self%mesh(3)%upper_boundary) then
    if (self%i_omega==0) then
      call self%boundary_condition
    end if
  end if
  call self%rt_solver_t%update
  call trace_end
END SUBROUTINE update

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
      f_scale = (tt(ix,iy,li(3))/top_temperature)**0.05
      h_scale = tt(ix,iy,li(3))/(scaling%temp*grav)*f_scale
      d_ratio = exp((z(iz)-z(li(3)))/h_scale)
      if (first_lower) then
        first_lower = .false.
        print *, self%id, 'experiment_t%boundary_condition (lower)', &
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
    !-----------------------------------------------------------------------------
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
      dlnddz = log(d(ix,iy,ui(3))/d(ix,iy,ui(3)-2))/(z(ui(3)-z(ui(3)-2)))
      dlnedz = log(e(ix,iy,ui(3))/e(ix,iy,ui(3)-2))/(z(ui(3)-z(ui(3)-2)))
      d_ratio = exp((z(iz)-z(ui(3)))*dlnddz)
      d (ix,iy,iz) = d (ix,iy,ui(3))*d_ratio
      px(ix,iy,iz) = px(ix,iy,ui(3))*d_ratio
      py(ix,iy,iz) = py(ix,iy,ui(3))*d_ratio
      pz(ix,iy,iz) = pz(ix,iy,ui(3))*d_ratio
      e(ix,iy,iz) = e(ix,iy,ui(3))*exp((z(iz)-z(ui(3)))*dlnedz)
      if (first_upper) then
        first_upper = .false.
        print *, self%id, 'experiment_t%boundary_condition (upper)', &
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
  call trace_end
END SUBROUTINE boundary_condition

!> =============================================================================
!> Stellar atmosphere with gravity pointing in the positive y-direction
!> =============================================================================
FUNCTION force (self) RESULT (ff)
  class(gpatch_t)   :: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3) :: ff
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%force')
  ff(:,:,:,1) = 0.0
  ff(:,:,:,2) = 0.0
  ff(:,:,:,3) = grav
  call trace_end
END FUNCTION force

END MODULE
