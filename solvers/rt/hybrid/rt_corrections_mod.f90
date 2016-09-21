MODULE rt_corrections_mod
  USE io_mod
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
  type, extends(rt_t), public     :: rt_correction_t
    type (rt_data_t),     pointer :: my_rt_data     => null()       ! corresponding rt data task
    type (rt_internal_t), pointer :: my_rt_internal => null()       ! corresponding internal RT task
    real (kind=4), dimension(:,:),  pointer :: I_in => null()       ! (wl_bin, ray) incoming intensity. 
    real (kind=8), dimension(:),    pointer :: I_in_time  => null() ! (ray, time)
    integer :: solver
  contains
    procedure :: init => init_rt_correction
    procedure :: update => update_rt_correction
    procedure :: check => check_rt_corrections
    procedure :: check_intensities => check_rt_get_intensities
    procedure :: get_upstream
    procedure :: get_local_static
    procedure :: get_local_dynamic
    procedure :: get_remote_static
    procedure :: get_remote_dynamic
    procedure :: get_correction_from_BC
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
  self%time      = -1
  self%dtime     = 3e-6
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

  allocate(self%i_in(microphysics%n_wl_bins, self%my_rays%n_rays))
  allocate(self%i_in_time(self%my_rays%n_rays))
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
  integer :: pt, ray, ipt, wl, p(3), lw(3), uw(3)
  real(kind=4), dimension(:,:), pointer :: q
  real(kind=4) :: dir_wgt
  real(kind=4), dimension(:), allocatable :: Ifw
  real :: temp
  character(len=16) :: stime
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_correction '//stime)
  dir_wgt = 2.*pi/self%n_ray_dirs
  self%it = self%mhd_patch%it
  call self%get_upstream                                                ! get upstream data, whatever it means
  associate (npt   => self%my_rays%n_points,    nwb => microphysics%n_wl_bins,  &
             nrays => self%my_rays%n_rays,      ex0 => self%my_rt_internal%ex0, &
                 s => self%my_rt_data%source,   op  => self%my_rt_data%opacity, &
              cell => self%my_rays%cell_coords, q_cool => self%mhd_patch%rt_q)
  
  allocate(Ifw(microphysics%n_wl_bins))
  
  
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
    
    Ifw = q(:,1) * ex0(:,1,ray) + s(p(1),p(2),p(3),:) 
    lw = self%my_rays%pt_wall(:,1,1,ray)
    uw = self%my_rays%pt_wall(:,2,1,ray) 
    call self%my_rt_data%I_store(Ifw,p,1,ray,self%dir_ID,lw,uw)         ! store forward intensity for 1st point

    if (npt(ray) > 1) then                                              ! for non-single-point rays only
      do pt = 2, npt(ray)                                               ! loop over all points in ray
      
        q(:,pt) = q(:,pt) + q(:,pt-1) * ex0(:,pt,ray)                   ! apply corrections
        p = cell(:,pt,ray)
        temp = self%mhd_patch%cell_gas_temperature(p)                   ! rt data and mhd mesh cell coords
        call self%mhd_patch%update_q( self%time, q(:,pt),p,nwb,&
        op(p(1),p(2),p(3),:),temp, dir_wgt)                             ! update q
        
        Ifw = q(:,pt) + s(p(1),p(2),p(3),:)
        lw = self%my_rays%pt_wall(:,1,pt,ray)
        uw = self%my_rays%pt_wall(:,2,pt,ray)
        call self%my_rt_data%I_store(Ifw,p,pt,ray,self%dir_ID,lw,uw)    ! store forward intensities                  

      end do  
    end if
    
    
  end do
  deallocate(Ifw)
  end associate
  
  if (io%verbose>2) call self%check
  call self%update_time                                                 ! update time
  
  call trace_end
END SUBROUTINE update_rt_correction

!> =====================================================================
!> Select appropriate method to get the upstream values
!> =====================================================================
SUBROUTINE get_upstream (self)
  class(rt_correction_t) :: self
  logical:: static, local
  class(link_t), pointer:: nbor
  !---------------------------------------------------------------------
  call trace_begin('rt_corrections_mod%get_upstream')
  
  if (self%is_set(bits%has_BC)) then  ! FIXME - make inner/outer BCs separate. Currently only outer BC -9.9000002145767212
    if (self%obc == 0) then
      error stop "huh?"
    end if
    call self%get_correction_from_BC                                    ! get boundary value
  else
    static = .true.; local = .true.                                     ! not moving? on this rank?
    if (self%is_clear(bits%static)) static = .false.
    nbor => self%link%nbor
    do while (associated(nbor))                                         ! check nbors for these logicals
      associate (this => nbor%task)
      SELECT type(this)
        class is (rt_correction_t)
          if (this%is_clear(bits%static)) static = .false.
          if (this%is_set(bits%virtual)) local = .false.
      end SELECT
      end associate
      nbor => nbor%next
    end do
    
    if (local) then                                                     ! select appropriate mode
      if (static) then
        call self%get_local_static
      else 
        call self%get_local_dynamic
      end if
    else
      if (static) then
        call self%get_remote_static
      else 
        call self%get_remote_dynamic
      end if  
    end if
  end if
  call trace_end
END SUBROUTINE get_upstream

!> =====================================================================
!> Upstream data is on same MPI rank and tasks are not moving.
!> =====================================================================
SUBROUTINE get_local_static (self)
  class(rt_correction_t) :: self
  class(rt_correction_t), pointer :: rt_nb
  class(link_t), pointer :: nbor
  integer:: ray, nb_cell(3), lw(3), uw(3), nb_ray, nb_pt, try
  logical, dimension(:), allocatable :: updated
  real(kind=4), dimension(microphysics%n_wl_bins) :: mean, Iin
  real(kind=8) :: n_mean ,my_c_pos(3)
  !---------------------------------------------------------------------  
  call trace_begin('rt_corrections_mod::get_local_static')
  
  allocate (updated(self%my_rays%n_rays))                               ! make sure all intensities are updated
  updated = .false.
  
  nbor => self%link%nbor
  do while(associated(nbor))
    nullify(rt_nb)
      
    associate (rt => nbor%task)
    SELECT type (rt)
    class is (rt_correction_t)
    
      if (self%the_right_wall(self%link, nbor,self%my_rays%lwall)) then                   ! make sure nbor is on the right wall
        if (rt%rt_type == self%rt_type) then                                              ! make sure nbor is of the right RT type
          rt_nb => rt
          do ray = 1, self%my_rays%n_rays                                                 ! loop  over my rays
            if (self%level >= rt_nb%level) then                                           ! if i am not larger than my nbor
              nb_cell = rt_nb%pos2cell(self%cell2pos(self%my_rays%fwd_coords(:,ray)))     ! corresponding cell at the nbor mesh
              if (all(nb_cell>=rt_nb%li).and.all(nb_cell<=rt_nb%ui)) then                 ! double check that we pointed to the right spot
                lw      = -self%my_rays%pt_wall(:,2,1,ray)                                ! select right walls for the ray
                uw      = -self%my_rays%pt_wall(:,1,1,ray)
                self%I_in   (:,ray) = rt_nb%my_rt_data%I_retrv(nb_cell,self%dir_ID,lw,uw) ! get the right intensity
                self%I_in_time(ray) = rt_nb%time                                          ! store the time of the update
                updated(ray) = .true.                                                     ! the ray was updated
              end if
            else                                                                          ! i am larger than my nbor
            
              my_c_pos = self%cell2pos(self%my_rays%fwd_coords(:,ray))                    ! get the absolute cell position in the comp. domain
              n_mean = 0                                                                  ! number of rays to be averaged
              mean = 0                                                                    ! mean intensity
              do nb_ray = 1, rt_nb%my_rays%n_rays                                         ! loop over nbor rays
                try = -1                                                                  ! indicator whether this ray meets the criterion
                do nb_pt = 1, rt_nb%my_rays%n_points(nb_ray)                              ! find the last point in the ray before the wall in question
                  if (self%point_in_cell(rt_nb%ray2pos(nb_ray,nb_pt),my_c_pos, self%cell_ds)) try = nb_pt
                end do
                if (try > 0) then
                  Iin = rt_nb%my_rt_data%I_retrv_ray(nb_ray,nb_pt,self%dir_ID)            ! retrieve this intensity
                  mean = mean + Iin
                  n_mean = n_mean + 1.
                end if
              end do
              if ((n_mean > 0).and.(self%time <= rt_nb%time)) then
                self%I_in(:,ray) = mean/n_mean
                self%I_in_time(ray) = rt_nb%time                                          ! store the time of the update
                updated(ray) = .true.                                                     ! the ray was updated
              end if
            end if
          end do
        end if
      end if

    end SELECT
    end associate
                
  nbor => nbor%next
  end do
  
  if(.not.all(updated)) then                                                              ! TO DO: old test, to be moved to self%check
    print *, 'task ID, level, time, n_nbors', self%id, self%level, self%time, self%n_nbors
    print *, updated
    print *, self%i_in_time
    print *, self%i_in
    error stop "not all rays updated for in local static rt"
  end if
  
  if (io%verbose>2) call self%check_intensities(1)
  call trace_end  
END SUBROUTINE get_local_static

!> =====================================================================
!> Upstream data is on different MPI rank, but tasks are not moving.
!> !!! TEMPORARY REMOVED !!!
!> =====================================================================
SUBROUTINE get_remote_static (self)
  class(rt_correction_t) :: self
  !---------------------------------------------------------------------  
  error stop "ERROR: remote same static not implemented"
END SUBROUTINE get_remote_static

SUBROUTINE get_local_dynamic (self)
  class(rt_correction_t) :: self
  !---------------------------------------------------------------------  
  error stop "ERROR: local different dynamic not implemented"
END SUBROUTINE get_local_dynamic

SUBROUTINE get_remote_dynamic (self)
  class(rt_correction_t) :: self
  !---------------------------------------------------------------------  
  error stop "ERROR: remote same dynamic not implemented"
END SUBROUTINE get_remote_dynamic


!> =====================================================================
!> Upstream data is from boundary conditions.
!> =====================================================================
SUBROUTINE get_correction_from_BC (self)
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
END SUBROUTINE get_correction_from_BC

!SUBROUTINE get_upstream_corr(self, nbor, my_ray, upd)
!  class(rt_correction_t)  :: self
!  class(rt_correction_t)  :: nbor
!  integer      :: my_ray, nb_c(3), nb_ray, nb_pt, try
!  logical      :: fnd, upd
!  real(kind=8) :: I_in, p(3), mean(microphysics%n_wl_bins), n_mean ,my_c_pos(3)
!-----------------------------------------------------------------------
!  upd = .false.
!  if (self%level >= nbor%level) then
!  
!    p = self%cell2pos(self%my_rays%fwd_coords(:,my_ray))
!    fnd = .false.
 !   call nbor%pos2ray(p,self%my_rays%lwall, nb_ray, nb_pt, fnd)  !here?
 !   if (fnd.and.(self%time <= nbor%time)) then
!@      upd = .true.
!      self%I_in(:,my_ray) =  nbor%intensity(:,nb_pt,nb_ray)
!      self%I_in_time(my_ray) = nbor%time
!    end if
!  else
!    my_c_pos = self%cell2pos(self%my_rays%fwd_coords(:,my_ray))  !here
!    n_mean = 0
!    mean = 0
!    do nb_ray = 1, nbor%my_rays%n_rays
!      try = -1
!      do nb_pt = 1, nbor%my_rays%n_points(nb_ray)    ! find the last point in the ray before the wall in question
!        if (self%point_in_cell(nbor%ray2pos(nb_ray,nb_pt),my_c_pos, self%cell_ds)) try = nb_pt
!      end do
!      if (try > 0) then
!        mean = mean + nbor%intensity(:, try, nb_ray)
!        n_mean = n_mean + 1.
!      end if
!    end do
!    if ((n_mean > 0).and.(self%time <= nbor%time)) then
!      self%I_in(:,my_ray) = mean/n_mean
!      self%I_in_time(my_ray) = nbor%time
!      upd = .true.
!    end if
!  end if
!END SUBROUTINE get_upstream_corr

SUBROUTINE get_downstream(self)
  class(rt_correction_t) :: self

END SUBROUTINE get_downstream

SUBROUTINE check_rt_get_intensities(self,method)
  class(rt_correction_t) :: self
  integer :: method
  !---------------------------------------------------------------------
   
  if (any(isnan(self%I_in))) then
    select case(method)
    case(1)
      print *,'ERROR: rt_corrections, get intensities: local static has NaN in I_in'
    case(2)
      print *,'ERROR: rt_corrections, get intensities: local dynamic has NaN in I_in'
    end select
    
    print *,'rt_type, task ID, BC bit:', self%rt_type, self%id, self%is_set(bits%has_bc)
    error stop "STOPPING"
  end if
  
END SUBROUTINE check_rt_get_intensities

SUBROUTINE check_rt_corrections(self)
  class(rt_correction_t) :: self
!-----------------------------------------------------------------------

  if(.not.associated(self%my_rt_data)) then
    error stop "rt_corrections%check: my_rt_data is not associated"
  end if
  
  if(.not.associated(self%my_rt_internal)) then
    error stop "rt_corrections%check: my_rt_internal is not associated"
  end if
  
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
