MODULE rt_data_mod
  USE bits_mod
  USE io_mod
  USE microphysics_mod
  USE ray_geometry_v2_mod
  USE rt_mod
  USE trace_mod
  implicit none
  private
!-----------------------------------------------------------------------
  type, extends(rt_t), public :: rt_data_t
    integer :: n_cells(3)
    integer :: n_tau_dat = 0                                            ! used for mpi comm
    real, dimension(:,:,:,:), pointer :: opacity => null()              ! (x,y,z,wl_bin)
    real, dimension(:,:,:,:), pointer :: source  => null()              ! (x,y,z,wl_bin)
  contains
    procedure :: init   => init_rt_data
    procedure :: update => update_rt_data
    procedure :: check  => check_rt_data
    procedure :: add_cr_heating
  end type
!-----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialise RT_data tables, which takes required MHD variables and gets
!> the source function and opacity per cell, mimicking the MHD patch cell
!> structure. It also stores intensity per direction per cell to communi-
!> cate it if patches overlap.
!> =====================================================================
SUBROUTINE init_rt_data(self)
  class(rt_data_t) :: self
  integer :: i, j, k, d
  !---------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_data')
  call self%rt_t%init
  self%type = 'rt_data_t'
  associate (mhd_patch => self%mhd_patch, n => self%n_cells, &
                   nwb => microphysics%n_wl_bins)
  !---------------------------------------------------------------------
  ! Copy relevant data from mhd_task and it's patch to perfectly shadow it
  !---------------------------------------------------------------------
  self%li       = mhd_patch%li
  self%ui       = mhd_patch%ui
  self%n_cells  = mhd_patch%gn
  self%position = mhd_patch%position                ! task position
  self%offset   = mhd_patch%offset
  self%level    = mhd_patch%level                   ! task level
  self%velocity = mhd_patch%velocity                ! task velocity
  self%size     = mhd_patch%size                    ! task size
  self%cell_ds  = mhd_patch%ds                      ! cell size
  self%it       = mhd_patch%it                      ! task time
  self%end_time = self%mhd_patch%end_time
  self%out_time = self%mhd_patch%out_time
  self%time     = -3                                ! initial time
  self%dtime    = 1e-6                              ! initial dtime
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set(bits%root_grid)
  !---------------------------------------------------------------------
  ! The rest of the setup
  !---------------------------------------------------------------------
  allocate(self%name)
  self%name       = 'rt_data'
  self%n_ray_dirs = ray_geometry%head%n_ray_dirs
  allocate(self%q      (n(1), n(2), n(3), self%nt)) ! net heating
  allocate(self%opacity(n(1), n(2), n(3), nwb))     ! opacity table
  allocate(self%source (n(1), n(2), n(3), nwb))     ! Source function table
  self%q = 0.0; self%opacity = 0;   self%source = 0
  end associate
  !
  if (io%verbose>1) then
    write (*,'(a,i5,3x,a,i5,3x,a,i3)') &
      '    rt_data_t%id:', self%id, &
      'mhd_patch%id:', self%mhd_patch%id
  end if
  if (io%verbose>2) call self%check(1)             ! check if rt_data initialized ok
  call trace_end
END SUBROUTINE init_rt_data

!> =====================================================================
!> Update RT data task values
!> =====================================================================
SUBROUTINE update_rt_data (self)
  class(rt_data_t) :: self
  character(len=16):: stime
  real, dimension(:,:,:), pointer :: rho, temp
  integer, save:: itimer=0
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_data', itimer=itimer)

  self%it  = self%mhd_patch%it

  allocate (temp(size(self%source,1),size(self%source,2),size(self%source,3)))
  associate  (l => self%li, u=>self%ui)
  rho => self%mhd_patch%mem(:,:,:,1,self%it,1)

  temp(l(1):u(1),l(2):u(2),l(3):u(3)) = self%mhd_patch%gas_temperature_ngz() ! estimate temperature

  self%opacity(l(1):u(1),l(2):u(2),l(3):u(3),:) = &
  microphysics%opacity%solve (rho(l(1):u(1),l(2):u(2),l(3):u(3)),&
  temp(l(1):u(1),l(2):u(2),l(3):u(3)), microphysics%n_wl_bins)  !FIXME : should pass wl_bins with info ! estimate opacity
  self%source (l(1):u(1),l(2):u(2),l(3):u(3),:) = &
  microphysics%source_fn%solve (temp(l(1):u(1),l(2):u(2),l(3):u(3)), microphysics%wl_bins) ! estimate the source function

  ! add cosmic ray heating
  call self%add_cr_heating(self%mhd_patch%time,temp(l(1):u(1),l(2):u(2),l(3):u(3)))  ! add cr heating
  call self%update_time

  end associate
  deallocate (temp)

  if (io%verbose>2) call self%check(2)

  call trace_end (itimer)
END SUBROUTINE update_rt_data

!> =============================================================================
!> Optionally add cosmic ray heating -- generic if applied before 'apply_heating'
!> =============================================================================
SUBROUTINE add_cr_heating (self, time, temp)
  class(rt_data_t) :: self
  real, dimension(:,:,:) :: temp
  real(kind=8) :: time
  integer :: l(3), u(3)
  ! ----------------------------------------------------------------------------
  call trace%begin ('mhd_rt_t_%add_cr_heating')
  l = self%mhd_patch%mesh%li
  u = self%mhd_patch%mesh%ui
  associate (rho => self%mhd_patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,self%it,1), &
               q => self%q(l(1):u(1),l(2):u(2),l(3):u(3),self%it))
  q = q + rho * microphysics%cr_heating
  end associate
  call trace%end
END SUBROUTINE add_cr_heating

!> =====================================================================
!> Test if rt_data module has any errors
!> =====================================================================
SUBROUTINE check_rt_data(self, plc)
  class(rt_data_t) :: self
  integer :: plc
!-----------------------------------------------------------------------
  if (any(isnan(self%opacity))) then
    print *,'=================== opacity ================'
    print *,self%opacity
    if (plc ==1) error stop "rt_data%init: self%opacity has NaN"
    if (plc ==2) error stop "rt_data%update: self%opacity has NaN"
  end if

  if (any(isnan(self%source))) then
    print *,'=================== source ================'
    print *,self%source
    if (plc ==1) error stop "rt_data%init: self%source has NaN"
    if (plc ==2) error stop "rt_data%update: self%source has NaN"
  end if
END SUBROUTINE check_rt_data

END MODULE rt_data_mod
