MODULE rt_corrections_mod
  USE bits_mod
  USE link_mod
  USE microphysics_mod
  USE rt_mod
  USE rt_data_mod
  USE rt_internal_mod
  USE rt_boundary_mod
  USE ray_geometry_mod
  USE trace_mod
  implicit none
  private
!-----------------------------------------------------------------------
  type :: rt_stream_t
    type (link_t)     , pointer :: nbor => null()
    type (rt_stream_t), pointer :: head => null()
    type (rt_stream_t), pointer :: next => null()
    type (rt_stream_t), pointer :: tail => null()
    integer :: method = -1
    integer :: n_entries = 0
    real(kind=8), dimension(:), pointer :: t   => null()
    integer,      dimension(:), pointer :: iit => null()
    integer,      dimension(:), pointer :: ray => null()
    integer,      dimension(:), pointer :: pt  => null()
  contains
    procedure :: init    => init_rt_stream
    procedure :: add_src => add_rt_stream_src
  end type

!-----------------------------------------------------------------------
  type, extends(rt_t), public     :: rt_correction_t
    type (rt_stream_t),   pointer :: recv
    type (rt_data_t),     pointer :: my_rt_data     => null()
    type (rt_internal_t), pointer :: my_rt_internal => null()
    real (kind=8), dimension(:,:,:), pointer :: intensity => null()  ! (wl_bin, point, ray), to communicate the corrections
    real (kind=8), dimension(:,:),   pointer :: I_in => null()       ! (wl_bin, ray) incoming intensity. 
    real (kind=8), dimension(:),     pointer :: I_in_time  => null() ! (ray, time)
    integer :: solver
  contains
    procedure :: init => init_rt_correction
    procedure :: update => update_rt_correction
    procedure :: check => check_rt_corrections
    procedure :: get_upstream
    procedure :: update_BC
  end type
!-----------------------------------------------------------------------
CONTAINS

SUBROUTINE init_rt_stream(self)
  class(rt_stream_t) :: self
  
END SUBROUTINE init_rt_stream

SUBROUTINE add_rt_stream_src(self, src)
  class(rt_stream_t) :: self
  class(link_t), optional :: src
  type(rt_stream_t), pointer :: new
!-----------------------------------------------------------------------

  if (.not.present(src)) then
    self%method = 0                    ! boundary
    
  else
    self%method = 1 ! determine the method. 
    allocate(new)
    if (.not.associated(self%head)) then
      self%head => new
      self%tail => new
    else
      self%tail%next => new
      self%tail => new
      new%head => self%head
    end if
    new%nbor  => src
    
  end if
  
END SUBROUTINE add_rt_stream_src

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
  call self%get_upstream
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

SUBROUTINE update_BC (self)
  class(rt_correction_t) :: self
  integer :: ray
  real(kind=8) :: ppos(3)
  !---------------------------------------------------------------------
  call trace_begin('rt_t%set_BC')
  SELECT case(self%rt_type)
  case(0)
  case(1)
  case(2,3)  
    do ray = 1, self%my_rays%n_rays
      ppos = self%absolute_point_position(ray, 1, self%rt_type)
      if (self%obct == 0) then
        call rt_boundary%set_obc(ppos,self%cell_ds,self%I_in(:,ray))
      end if
    end do
  end SELECT
  call trace_end
END SUBROUTINE update_BC

SUBROUTINE get_upstream(self)
  class(rt_correction_t) :: self
  !---------------------------------------------------------------------
  select CASE(self%recv%method)
    case (0)
      call self%set_BC
    case (1)
      call self%get_upstream_direct(self%recv%nbor) 

SUBROUTINE get_upstream_corr(self, nbor, my_ray, upd)
  class(rt_correction_t)  :: self
  class(rt_correction_t)  :: nbor
  integer      :: my_ray, nb_c(3), nb_ray, nb_pt, try
  logical      :: fnd, upd
  real(kind=8) :: I_in, p(3), mean(microphysics%n_wl_bins), n_mean ,my_c_pos(3)
!-----------------------------------------------------------------------
  upd = .false.
  if (self%level >= nbor%level) then
  
    p = self%cell2pos(self%my_rays%fwd_coords(:,my_ray))
    fnd = .false.
    call nbor%pos2ray(p,self%my_rays%lwall, nb_ray, nb_pt, fnd)  !here?
    if (fnd.and.(self%time <= nbor%time)) then
      upd = .true.
      self%I_in(:,my_ray) =  nbor%intensity(:,nb_pt,nb_ray)
      self%I_in_time(my_ray) = nbor%time
    end if
  else
    my_c_pos = self%cell2pos(self%my_rays%fwd_coords(:,my_ray))  !here
    n_mean = 0
    mean = 0
    do nb_ray = 1, nbor%my_rays%n_rays
      try = -1
      do nb_pt = 1, nbor%my_rays%n_points(nb_ray)    ! find the last point in the ray before the wall in question
        if (point_in_cell(nbor%ray2pos(nb_ray,nb_pt),my_c_pos, self%cell_ds)) try = nb_pt
      end do
      if (try > 0) then
        mean = mean + nbor%intensity(:, try, nb_ray)
        n_mean = n_mean + 1.
      end if
    end do
    if ((n_mean > 0).and.(self%time <= nbor%time)) then
      self%I_in(:,my_ray) = mean/n_mean
      self%I_in_time(my_ray) = nbor%time
      upd = .true.
    end if
  end if
END SUBROUTINE get_upstream_corr

SUBROUTINE get_downstream(self)
  class(rt_correction_t) :: self

END SUBROUTINE get_downstream

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
