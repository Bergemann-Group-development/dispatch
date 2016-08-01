MODULE rt_data_mod
  USE bits_mod
  USE microphysics_mod
  USE rt_mod
  USE trace_mod
  USE io_mod
  implicit none
!----------------------------------------------------------------------- 
  type, extends(rt_t), public :: rt_data_t
    integer                             :: n_cells(3)
    real(kind=8), dimension(:,:,:,:),   pointer :: opacity => null()  ! (x,y,z,wl_bin)
    real(kind=8), dimension(:,:,:,:),   pointer :: source  => null()  ! (x,y,z,wl_bin)
  contains
    procedure :: init   => init_rt_data
    procedure :: update => update_rt_data
    procedure :: check  => check_rt_data
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
  !---------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_data')
  call self%task_t%init
  associate (mhd_patch => self%mhd_patch, nc => self%n_cells, &
                   nwb => microphysics%n_wl_bins)
  !---------------------------------------------------------------------
  ! Copy relevant data from mhd_task and it's patch to perfectly shadow it
  !---------------------------------------------------------------------
  self%li         = mhd_patch%li
  self%ui         = mhd_patch%ui
  self%n_cells(:) = mhd_patch%gn(:)
  self%position   = mhd_patch%position                ! task position
  self%offset     = mhd_patch%offset
  self%level      = mhd_patch%level                   ! task level
  self%velocity   = mhd_patch%velocity                ! task velocity
  self%size       = mhd_patch%size                    ! task size
  self%cell_ds    = mhd_patch%ds                      ! cell size
  self%it         = mhd_patch%it                      ! task time
  self%end_time   = self%mhd_patch%end_time
  self%out_time   = self%mhd_patch%out_time
  self%time       = -3
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set(bits%root_grid)
  !---------------------------------------------------------------------
  ! The rest of the setup
  !---------------------------------------------------------------------
  allocate(self%name)
  self%name       = 'rt_data'
  self%n_ray_dirs = ray_geometry%ray_dirs
  allocate(self%opacity(nc(1), nc(2), nc(3), nwb)) ! opacity table
  allocate(self%source (nc(1), nc(2), nc(3), nwb)) ! Source function table
  self%opacity = 0; self%source = 0
  end associate
  
  if (io%verbose>2) call self%check       
  call trace_end  
END SUBROUTINE init_rt_data

!> =====================================================================
!> Update RT data task values
!> =====================================================================
SUBROUTINE update_rt_data (self)
  class(rt_data_t) :: self
  character(len=16):: stime
  real, dimension(:,:,:), allocatable :: temp
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_data '//stime)
  
  self%it  = self%mhd_patch%it
  
  allocate (temp(size(self%source,1),size(self%source,2),size(self%source,3)))
  associate  (rho => self%mhd_patch%mem(:,:,:,1,self%it,1), &
                l => self%li, u=>self%ui) 
  
  temp = self%mhd_patch%gas_temperature()         
  self%opacity = microphysics%opacity%solve   (rho, temp, microphysics%n_wl_bins)  !FIXME : should pass wl_bins with info
  self%source  = microphysics%source_fn%solve (     temp, microphysics%wl_bins)
  ! add cosmic ray heating
  call self%mhd_patch%add_cr_heating(self%time,temp)
  call self%update_time
  
  end associate
  deallocate (temp)
  
  if (io%verbose>2) call self%check
    
  call trace_end
END SUBROUTINE update_rt_data

SUBROUTINE check_rt_data(self)
  class(rt_data_t) :: self
!-----------------------------------------------------------------------
  if (any(isnan(self%opacity))) then
    print *,'=================== opacity ================'
    print *,self%opacity
    error stop "rt_data%init: self%opacity has NaN"
  end if
  
  if (any(isnan(self%source))) then
    print *,'=================== source ================'
    print *,self%source
    error stop "rt_data%init: self%source has NaN"
  end if
END SUBROUTINE check_rt_data

END MODULE rt_data_mod
