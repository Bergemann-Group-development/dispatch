MODULE rt_corrections_mod
  USE bits_mod
  USE link_mod
  USE microphysics_mod
  USE rt_mod
  USE rt_data_mod
  USE rt_internal_mod
  USE ray_geometry_mod
  USE trace_mod
  implicit none
  private
!-----------------------------------------------------------------------
  type :: rt_stream_t
    type (rt_stream_t), pointer :: head => null()
    type (rt_stream_t), pointer :: next => null()
    type (rt_stream_t), pointer :: tail => null()
    type (link_t)     , pointer :: link => null()
    integer, dimension(:), pointer :: ray, pt
  end type

!-----------------------------------------------------------------------
  type, extends(rt_t), public :: rt_correction_t
    type (rt_data_t),     pointer :: my_rt_data     => null()
    type (rt_internal_t), pointer :: my_rt_internal => null()
    type (rt_stream_t),   pointer :: send_list => null()
    type (rt_stream_t),   pointer :: recv_list => null()
    integer :: solver
  contains
    procedure :: init => init_rt_correction
    procedure :: update => update_rt_correction
    procedure :: check => check_rt_corrections
    procedure :: send_downstream
  end type
!-----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialise boundary table for upstream values
!> =====================================================================
SUBROUTINE init_rt_correction (self)
  class(rt_correction_t) :: self
  class(ray_list_t), pointer :: my_rays
  integer :: ray
  real(kind=8) :: p(3)
!-----------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_correction')
  call self%task_t%init
  associate (mhd_patch => self%mhd_patch)
  self%position  = mhd_patch%position
  self%offset    = mhd_patch%offset
  self%level     = mhd_patch%level
  self%velocity  = mhd_patch%velocity
  self%size      = mhd_patch%size
  self%cell_ds   = mhd_patch%ds
  self%it        = mhd_patch%it
  self%end_time  = self%mhd_patch%end_time
  self%out_time  = self%mhd_patch%out_time
  self%time       = -1
  self%li        = mhd_patch%li
  self%ui        = mhd_patch%ui
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set (bits%root_grid)
  
  allocate(self%name)
  
  ! FIXME: add the right corrections for the right solvers.
  
  SELECT case(self%rt_type)
    case(2)
      call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .false.)                ! retrieve ray list for this direction
      self%name = 'rt_forward'
    case(3)
      call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .true. )                ! retrieve ray list for this direction
      self%name = 'rt_reverse'
  end SELECT
  self%my_rays => my_rays

  allocate(self%intensity(microphysics%n_wl_bins, self%my_rays%max_n_points, self%my_rays%n_rays))
  allocate(self%i_in(microphysics%n_wl_bins, self%my_rays%n_rays))
  allocate(self%i_in_time(self%my_rays%n_rays))
  self%intensity = 0  
  self%i_in(:,:) = 0
  self%i_in_time = self%time
  end associate
  call trace_end
END SUBROUTINE init_rt_correction

!> =====================================================================
!> Forward/reverse RT correction
!> =====================================================================
SUBROUTINE update_rt_correction (self)
  class(rt_correction_t) :: self
  real, dimension(microphysics%n_wl_bins) :: I_in
  integer :: pt, ray, ipt, wl, p(3)
  real(kind=8), dimension(:,:), pointer :: q
  real(kind=8) :: dir_wgt
  real :: temp
  character(len=16) :: stime
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_correction '//stime)
  dir_wgt = 2.*pi/self%n_ray_dirs
  self%it = self%mhd_patch%it
  associate (npt   => self%my_rays%n_points,    nwb => microphysics%n_wl_bins,  &
             nrays => self%my_rays%n_rays,      ex0 => self%my_rt_internal%ex0, &
                 s => self%my_rt_data%source,   op  => self%my_rt_data%opacity, &
              cell => self%my_rays%cell_coords, Ifw => self%intensity,          &
            q_cool => self%mhd_patch%rt_q)
  
  
  do ray = 1, nrays                                                     ! loop over rays
  
    SELECT case(self%rt_type)                                           ! select appropriate internal direction
    case(2)
      q => self%my_rt_internal%qfw(:,:,ray)
    case(3)
      q => self%my_rt_internal%qrv(:,:,ray)
    end SELECT
    
    I_in = self%i_in(:,ray)                                             ! get incoming intensity
   
    p = cell(:,1,ray)                                                   ! rt data cell coords
    q(:,1) = I_in - s(p(1),p(2),p(3),:)                                 ! correction for first (or single) point
    
    temp = self%mhd_patch%cell_gas_temperature(p)
    call self%mhd_patch%update_q( self%time, q(:,1),p,nwb,&
    op(p(1),p(2),p(3),:),temp, dir_wgt)                                 ! update q
    Ifw(:,1,ray) = q(:,1) * ex0(:,1,ray) + s(p(1),p(2),p(3),:)          ! store forward intensity

    if (npt(ray) > 1) then                                              ! for non-single-point rays only
      do pt = 2, npt(ray)                                               ! loop over all points in ray
      
        q(:,pt) = q(:,pt) + q(:,pt-1) * ex0(:,pt,ray)                   ! apply corrections
        p = cell(:,pt,ray)
        temp = self%mhd_patch%cell_gas_temperature(p)                   ! rt data and mhd mesh cell coords
        call self%mhd_patch%update_q( self%time, q(:,pt),p,nwb,&
        op(p(1),p(2),p(3),:),temp, dir_wgt)                             ! update q
        Ifw(:,pt,ray) = q(:,pt) + s(p(1),p(2),p(3),:)                   ! store forward intensities

      end do
      
    end if
    
    
  end do
  end associate
  
  if (io%verbose>2) call self%check
  call self%update_time
  
  call trace_end
END SUBROUTINE update_rt_correction

SUBROUTINE send_downstream(self)
  class(rt_correction_t) :: self

END SUBROUTINE send_downstream

SUBROUTINE check_rt_corrections(self)
  class(rt_correction_t) :: self
!-----------------------------------------------------------------------

  if(.not.associated(self%my_rt_data)) then
    error stop "rt_corrections%check: my_rt_data is not associated"
  end if
  
  if(.not.associated(self%my_rt_internal)) then
    error stop "rt_corrections%check: my_rt_internal is not associated"
  end if
  
  !if(.not.associated(self%my_rt_downstream)) then
  !  error stop "rt_corrections%check: my_rt_downstream is not associated"
  !end if
  
  if(.not.associated(self%mhd_patch)) then
    error stop "rt_corrections%check: rt_t%mhd_patch is not associated"
  end if
  
  if(.not.associated(self%mhd_link)) then
    error stop "rt_corrections%check: rt_t%mhd_link is not associated"
  end if
  
  if(.not.associated(self%my_nb_link)) then
    error stop "rt_corrections%check: rt_t%my_nb_link is not associated"
  end if
  
  if(.not.associated(self%my_rays)) then
    error stop "rt_corrections%check: rt_t%my_rays is not associated"
  end if
  
  if(.not.associated(self%intensity)) then
    error stop "rt_corrections%check: rt_t%intensity is not associated"
  else
    if (any(isnan(self%intensity))) then
      print *,'=================== ds ================'
      print *,self%intensity
      error stop "rt_corrections%check: self%intensity has NaN"
    end if
  end if
  
  if(.not.associated(self%i_in)) then
    error stop "rt_corrections%check: rt_t%i_in is not associated"
  else
    if (any(isnan(self%i_in))) then
      print *,'=================== ds ================'
      print *,self%i_in
      error stop "rt_corrections%check: self%i_in has NaN"
    end if
  end if
  
  if(.not.associated(self%i_in_time)) then
    error stop "rt_corrections%check: rt_t%i_in_time is not associated"
  else
    if (any(isnan(self%i_in_time))) then
      print *,'=================== ds ================'
      print *,self%i_in_time
      error stop "rt_corrections%check: self%i_in_time has NaN"
    end if
  end if
  
END SUBROUTINE check_rt_corrections

END MODULE rt_corrections_mod
