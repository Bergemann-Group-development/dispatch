MODULE rt_internal_mod
  USE bits_mod
  USE link_mod
  USE microphysics_mod
  USE rt_mod
  USE rt_data_mod
  USE rt_boundary_mod
  USE ray_geometry_mod
  USE trace_mod
  implicit none
!-----------------------------------------------------------------------
  type,extends(rt_t), public :: rt_internal_t
    type (rt_data_t),     pointer :: my_rt_data => null()
    real, dimension(:,:), pointer :: ds         => null()               ! (point,ray)
    real, dimension(:,:), pointer :: I_in => null()                     ! (wl_bin, ray) incoming intensity. 
    real, dimension(:),   pointer :: I_in_time  => null()               ! (ray, time)
    integer :: solver
    procedure (solve_internal_integral), pointer :: solve => null()
  contains
    procedure :: init   => init_rt_internal
    procedure :: update => update_rt_internal
    procedure :: check  => check_rt_internal
    procedure :: check_intensities => check_rt_get_intensities
    procedure :: get_upstream
    procedure :: get_local_static
    procedure :: get_local_dynamic
    procedure :: get_remote_static
    procedure :: get_remote_dynamic
    procedure :: get_correction_from_BC
    procedure :: solve_internal_integral
    procedure :: solve_internal_feautrier
  end type
!-----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialise tables for internal RT solution
!> =====================================================================
SUBROUTINE init_rt_internal (self)
  class(rt_internal_t) :: self
  class(ray_list_t),    pointer :: my_rays
  integer :: point, ray
  real(kind=4) :: p(3)
!-----------------------------------------------------------------------
  call trace_begin('rt_t%init_rt_internal')
  call self%task_t%init
  associate (mhd_patch => self%mhd_patch)
  self%position  = mhd_patch%position
  self%offset    = mhd_patch%offset
  self%level     = mhd_patch%level
  self%velocity  = mhd_patch%velocity
  self%size      = mhd_patch%size
  self%cell_ds   = mhd_patch%ds
  self%it        = mhd_patch%it
  self%end_time  = mhd_patch%end_time
  self%out_time  = mhd_patch%out_time
  self%li        = mhd_patch%li
  self%ui        = mhd_patch%ui
  self%time      = -2
  self%dtime     = 2e-6
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set (bits%root_grid)
  !if (mhd_patch%is_set(bits%has_BC))    call self%set (bits%has_BC)

  allocate(self%name)
  SELECT case(self%rt_type)
    case(2)
      call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .false.)                ! retrieve ray list for this direction
      self%name = 'rt_forward'
    case(3)
      call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .false. )                ! retrieve ray list for this direction
      self%name = 'rt_reverse'
  end SELECT
  self%my_rays => my_rays

  if (self%solver ==  0) then                                        ! set the right solver
    self%solve => solve_internal_feautrier
  else
    self%solve => solve_internal_integral
  end if

  !---------------------------------------------------------------------
  ! Allocate memory slots.
  !---------------------------------------------------------------------
  allocate(self%i_in (microphysics%n_wl_bins, self%my_rays%n_rays))
  allocate(self%i_in_time(self%my_rays%n_rays))
  allocate(self%ds    (my_rays%max_n_points, my_rays%n_rays))
  self%i_in(:,:) = 0
  self%i_in_time = self%time
  
  do ray = 1, self%my_rays%n_rays                                       ! loop over rays to copy coords and ds
    do point = 1, self%my_rays%n_points(ray)
      self%ds(point, ray) = my_rays%ds(point, ray) * self%cell_ds(1)
    end do
  end do

  if (io%verbose>2) call self%check
  end associate
  call trace_end
END SUBROUTINE init_rt_internal

!> =====================================================================
!> Update internal RT solution
!> =====================================================================
SUBROUTINE update_rt_internal (self)
  class(rt_internal_t) :: self
  integer :: wl, ray, point, p(3)
  character(len=16) :: stime
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_internal '//stime)

  call self%solve
  call self%update_time

  if (io%verbose>2) call self%check

  call trace_end

END SUBROUTINE update_rt_internal

!> =====================================================================
!> Select appropriate method to get the upstream values
!> =====================================================================
SUBROUTINE get_upstream (self)
  class(rt_internal_t) :: self
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
        class is (rt_internal_t)
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
  class(rt_internal_t) :: self
  class(rt_internal_t), pointer :: rt_nb
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
    class is (rt_internal_t)
    
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
              n_mean = 0d0                                                                ! number of rays to be averaged
              mean = 0                                                                    ! mean intensity
              do nb_ray = 1, rt_nb%my_rays%n_rays                                         ! loop over nbor rays
                try = -1                                                                  ! indicator whether this ray meets the criterion
                do nb_pt = 1, rt_nb%my_rays%n_points(nb_ray)                              ! find the last point in the ray before the wall in question
                  if (point_in_cell(rt_nb%ray2pos(nb_ray,nb_pt),my_c_pos, self%cell_ds)) try = nb_pt
                end do
                if (try > 0) then
                  Iin = rt_nb%my_rt_data%I_retrv_ray(nb_ray,nb_pt,self%dir_ID)            ! retrieve this intensity
                  mean = mean + Iin
                  n_mean = n_mean + 1d0
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
  class(rt_internal_t) :: self
  !---------------------------------------------------------------------  
  error stop "ERROR: remote same static not implemented"
END SUBROUTINE get_remote_static

SUBROUTINE get_local_dynamic (self)
  class(rt_internal_t) :: self
  !---------------------------------------------------------------------  
  error stop "ERROR: local different dynamic not implemented"
END SUBROUTINE get_local_dynamic

SUBROUTINE get_remote_dynamic (self)
  class(rt_internal_t) :: self
  !---------------------------------------------------------------------  
  error stop "ERROR: remote same dynamic not implemented"
END SUBROUTINE get_remote_dynamic


!> =====================================================================
!> Upstream data is from boundary conditions.
!> =====================================================================
SUBROUTINE get_correction_from_BC (self)
  class(rt_internal_t) :: self
  integer :: ray
  real(kind=8) :: ppos(3)
  !---------------------------------------------------------------------
  call trace_begin('rt_internal_t%get_correction_from_BC')
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

!> =====================================================================
!> Solution of RTE using the Integral formulation. Iwf and Irv are 
!> stored in case of overlapping patches.
!> =====================================================================

SUBROUTINE solve_internal_integral (self)
  class(rt_internal_t) :: self
  real(kind=4), dimension(:,:), allocatable :: ex0, ex1, ex2, ds1, d1s, d2s
  real(kind=4), dimension(:,:), allocatable :: s, dtau, idtau, ctau, q
  real(kind=4), dimension(microphysics%n_wl_bins)   :: Ifw, I_in
  integer :: pt, ipt, it, ray, wl, p(3), lw(3), uw(3)
  real    :: temp
  real(kind=4) :: dir_wgt
  !---------------------------------------------------------------------
  call trace_begin('rt_t%update_internal%solve_internal_integral')
  call self%get_upstream 
  
  dir_wgt = 2.*pi/self%n_ray_dirs
  associate (npt   => self%my_rays%n_points,nwb => microphysics%n_wl_bins, &
             nrays => self%my_rays%n_rays,  rtd => self%my_rt_data, &
               cell=> self%my_rays%cell_coords, &
               mnp => self%my_rays%max_n_points)

  allocate( dtau(nwb, mnp))
  allocate(idtau(nwb, mnp))
  allocate(ex0  (nwb, mnp))
  allocate(ex1  (nwb, mnp))
  allocate(ex2  (nwb, mnp))
  allocate(ds1  (nwb, mnp))
  allocate(d1s  (nwb, mnp))
  allocate(d2s  (nwb, mnp))
  allocate(s    (nwb, mnp))
  allocate(q    (nwb, mnp))
  allocate(ctau (nwb, mnp))



  do ray = 1, nrays

    do pt = 1, npt(ray)
      p = cell(:,pt,ray)
      do wl = 1, nwb
        dtau (wl,pt) = rtd%opacity(p(1),p(2),p(3),wl)*self%ds(pt,ray)
           s (wl,pt) = rtd%source (p(1),p(2),p(3),wl)
      end do
    end do
    
    I_in = self%i_in(:,ray)

    pt = npt(ray)
    idtau(:,1:pt) = 1_8 / dtau(:,1:pt)
      ex0(:,2:pt) = exp(-dtau(:,2:pt))                                  ! first exponential coeff
      ex1(:,2:pt) = 1_8-ex0(:,2:pt)                                     ! second exponential coeff
      ex2(:,2:pt) = ex1(:,2:pt)-dtau(:,2:pt) *   ex0(:,2:pt)            ! third exponential coeff
      ds1(:,2:pt) = (s(:,2:pt) -s(:,1:pt-1)) * idtau(:,2:pt)            ! half-point centered dS/dtau
      !-----------------------------------------------------------------
      !  dS/dtau and d2S/dtau2
      !-----------------------------------------------------------------
     ctau(:,2:pt-1) = 1_8 / (dtau(:,3:pt) + dtau(:,2:pt-1))
      d1s(:,2:pt-1) = &
     (ds1(:,3:pt  ) * dtau(:,2:pt-1) + ds1(:,2:pt-1) * &
     dtau(:,3:pt) ) * ctau(:,2:pt-1)                                    ! dS/dtau, 2nd. order accurate
      d2s(:,2:pt-1) = (ds1(:,3:pt  ) - ds1(:,2:pt-1))* &
               2_8*ctau(:,2:pt-1)                                       ! d2S/dtau2, 2nd. order accurate


      d2s(:,1)  = d2s(:,2)                                              ! assume constant d2S/dtau2
      d1s(:,1)  = d1s(:,2)    - dtau(:,2)  * d2s(:,2)                   ! extrapolate to 1
      d2s(:,pt) = d2s(:,pt-1)                                           ! assume constant d2S/dtau2
      d1s(:,pt) = d1s(:,pt-1) + dtau(:,pt) * d2s(:,pt)                  ! extrapolate to n
      !-----------------------------------------------------------------
      !  Internal forward or reverse solutions
      !-----------------------------------------------------------------
      
      select case(self%rt_type)
        case(2)                                                         ! forward
        q(:,1) = I_in(:) - s(:,1)                                       ! initial forward Q
        Ifw = I_in * ex0(:,1) 
        lw = self%my_rays%pt_wall(:,1,1,ray)
        uw = self%my_rays%pt_wall(:,2,1,ray)
        p = cell(:,1,ray)
        call self%my_rt_data%I_store(Ifw,p,1,ray,self%dir_ID,lw,uw)     ! store forward intensities 
        temp = self%mhd_patch%cell_gas_temperature(p)
        call self%mhd_patch%update_q( self%time, q(:,1),p,nwb,&
        rtd%opacity(p(1),p(2),p(3),:),temp, dir_wgt)                    ! update q
        do pt = 2, npt(ray)
          q(:,pt) = q(:,pt-1) * ex0(:,pt) - d1s(:,pt) * ex1(:,pt) &
                  + d2s(:,pt) * ex2(:,pt)
          Ifw = q(:,pt) + s(:,pt)
          lw = self%my_rays%pt_wall(:,1,pt,ray)
          uw = self%my_rays%pt_wall(:,2,pt,ray)
          p = cell(:,pt,ray)
          call self%my_rt_data%I_store(Ifw,p,pt,ray,self%dir_ID,lw,uw)  ! store forward intensities 
          temp = self%mhd_patch%cell_gas_temperature(p)
          call self%mhd_patch%update_q( self%time, q(:,pt),p,nwb,&
          rtd%opacity(p(1),p(2),p(3),:),temp, dir_wgt)                  ! update q
        end do
        
        case(3)                                                         ! reverse
        q(:,pt) = I_in(:) - s(:,pt)                                       ! initial forward Q
        Ifw = I_in * ex0(:,pt) 
        lw = self%my_rays%pt_wall(:,2,pt,ray)
        uw = self%my_rays%pt_wall(:,1,pt,ray)
        p = cell(:,pt,ray)
        call self%my_rt_data%I_store(Ifw,p,pt,ray,self%dir_ID,lw,uw)     ! store forward intensities 
        temp = self%mhd_patch%cell_gas_temperature(p)
        call self%mhd_patch%update_q( self%time, q(:,pt),p,nwb,&
        rtd%opacity(p(1),p(2),p(3),:),temp, dir_wgt)                    ! update q
        do pt = npt(ray)-1, 1, -1
          q(:,pt) = q(:,pt+1) * ex0(:,pt+1) - d1s(:,pt) * ex1(:,pt+1) &
                  + d2s(:,pt) * ex2(:,pt+1)
          Ifw = q(:,pt) + s(:,pt)
          lw = self%my_rays%pt_wall(:,2,pt,ray)
          uw = self%my_rays%pt_wall(:,1,pt,ray)
          p = cell(:,pt,ray)
          call self%my_rt_data%I_store(Ifw,p,pt,ray,self%dir_ID,lw,uw)  ! store forward intensities 
          temp = self%mhd_patch%cell_gas_temperature(p)
          call self%mhd_patch%update_q( self%time, q(:,pt),p,nwb,&
          rtd%opacity(p(1),p(2),p(3),:),temp, dir_wgt)                  ! update q
        end do
      end select
  end do

  deallocate(ex0, ex1, ex2, ds1, d1s, d2s, ctau, s, dtau, idtau, q)
  end associate
  call trace_end
END SUBROUTINE solve_internal_integral

!> =====================================================================
!> Solution of RTE using the Feautrier formulation. Forward and reverse
!> incoming radiation is ignored.
!> =====================================================================
SUBROUTINE solve_internal_feautrier (self)
  class(rt_Internal_t) :: self
  real, dimension(:,:,:), pointer :: a1, a2, a3, d
  integer                         :: j, k, l
  real                            :: a, dsdtau_1, dsdtau_2, d2sdtau2

END SUBROUTINE solve_internal_feautrier

!> =====================================================================
!> Test if rt_internal module has any errors
!> =====================================================================
SUBROUTINE check_rt_internal(self)
  class(rt_internal_t) :: self
!-----------------------------------------------------------------------
  if(.not.associated(self%my_rt_data)) then
    error stop "rt_internal%check: my_rt_data is not associated"
  end if

  if(.not.associated(self%solve)) then
    error stop "rt_internal%check: solver is not set"
  end if

  if (any(isnan(self%ds))) then
    print *,'=================== ds ================'
    print *,self%ds
    error stop "rt_internal%check: self%ds has NaN"
  end if

END SUBROUTINE check_rt_internal

SUBROUTINE check_rt_get_intensities(self,method)
  class(rt_internal_t) :: self
  integer :: method
  !---------------------------------------------------------------------
   
  if (any(isnan(self%I_in))) then
    select case(method)
    case(1)
      print *,'ERROR: rt_internal, get intensities: local static has NaN in I_in'
    case(2)
      print *,'ERROR: rt_internal, get intensities: local dynamic has NaN in I_in'
    end select
    
    print *,'rt_type, task ID, BC bit:', self%rt_type, self%id, self%is_set(bits%has_bc)
    error stop "STOPPING"
  end if
  
END SUBROUTINE check_rt_get_intensities

SUBROUTINE check_rt_corrections(self)
  class(rt_internal_t) :: self
!-----------------------------------------------------------------------

  if(.not.associated(self%my_rt_data)) then
    error stop "rt_corrections%check: my_rt_data is not associated"
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

END MODULE rt_internal_mod
