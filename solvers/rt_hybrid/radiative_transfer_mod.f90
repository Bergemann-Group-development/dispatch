MODULE radiative_transfer_mod
  USE rt_mod
  USE experiment_mod
  USE bits_mod
  USE io_mod
  USE list_mod
  USE task_list_mod
  USE omp_mod
  USE patch_mod
  USE mhd_mod
  USE ray_geometry_mod
  USE microphysics_mod
  USE scaling_mod
  USE task_mod
  USE trace_mod
  USE units_mod
  USE experiment_mod
  USE rt_boundary_mod
!----------------------------------------------------------------------- 
  type, extends(rt_t), public :: rt_data_t
    integer                             :: n_cells(3)
    real(kind=8), dimension(:,:,:,:),   pointer :: opacity => null()  ! (x,y,z,wl_bin)
    real(kind=8), dimension(:,:,:,:),   pointer :: source  => null()  ! (x,y,z,wl_bin)
  contains
    procedure :: init   => init_rt_data
    procedure :: update => update_rt_data
  end type
!-----------------------------------------------------------------------
 type,extends(rt_t), public :: rt_internal_t
    type (rt_data_t),               pointer :: my_rt_data => null()
    real,         dimension(:,:),   pointer :: ds         => null()  ! (point,ray)
    real(kind=8), dimension(:,:,:), pointer :: source     => null()  ! (wl_bin, point, ray) source function in ray domain
    real(kind=8), dimension(:,:,:), pointer :: dtau       => null()  ! (wl_bin, point, ray)
    real(kind=8), dimension(:,:,:), pointer :: qfw        => null()  ! (wl_bin, point, ray)
    real(kind=8), dimension(:,:,:), pointer :: qrv        => null()  ! (wl_bin, point, ray)
    real(kind=8), dimension(:,:,:), pointer :: ex0        => null()  ! (wl_bin, point, ray), first exponential coeff, for Integral solver
    procedure (solve_internal_integral), pointer :: solve_internal => null()
  contains
    procedure :: init   => init_internal
    procedure :: update => update_internal
  end type
!-----------------------------------------------------------------------
  type, extends(rt_t), public :: rt_correction_t
    type (rt_data_t),     pointer :: my_rt_data     => null()
    type (rt_internal_t), pointer :: my_rt_internal => null()
  contains
    procedure :: init   => init_correction
    procedure :: update => update_correction
  end type
!-----------------------------------------------------------------------

  PUBLIC rt_setup
CONTAINS

SUBROUTINE rt_setup (self)
  class(list_t)                :: self
  class(task_list_t),  pointer :: task_list                  ! local task list
  class(list_t),       pointer :: list
  class(link_t),       pointer :: t_list, link, rt_data_l, internal_l
  class(experiment_t), pointer :: mhd_patch    
  class(rt_data_t),    pointer :: rt_data
  class(rt_internal_t),pointer :: internal
  class(rt_correction_t), pointer :: correction
  class(ray_list_t),   pointer :: ray_bundle, ray_dir
  logical, save                :: first_time=.true.
  character(len=16)            :: RT_solver ='integral'      ! choice of RT solver
  integer:: solver, rt_type, d, b, r
  real(kind=8) :: cr_heating=4.4e-27
  namelist /rt_params/ RT_solver, cr_heating
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer::rt_setup')
  associate (mhd => self%head%task)
  SELECT TYPE(mhd)
    type is (experiment_t)
    mhd_patch => mhd
  end SELECT
  end associate
  !$omp critical (rt_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, rt_params)
    if (io%master) write (*, rt_params)
    call ray_geometry%init(mhd_patch%li,mhd_patch%ui)  ! FIXME: send cell wall sizes and dimensions as well 

    print *, 'Ray geometry: n_bundles, n_ray_dirs:',ray_geometry%n_bundles, ray_geometry%ray_dirs
    call microphysics%init
    call rt_boundary%init
    rt_boundary%ubox = mhd_patch%box
    rt_boundary%lbox = mhd_patch%origin
    SELECT case(RT_solver)
    case('feautrier')
      solver = 0
    case('integral')
      solver = 1
    end SELECT
  end if
  allocate(task_list)
  call task_list%list_t%init ('rt task list')
  t_list => self%head
  do while (associated(t_list))

    !-------------------------------------------------------------------
    ! Get the right type of the task (mhd_t)
    !-------------------------------------------------------------------
    associate (mhd => t_list%task)
    SELECT TYPE(mhd)
    type is (experiment_t)
      mhd_patch => mhd
      call mhd_patch%init_rt_heating
      !-----------------------------------------------------------------------------
      ! Construct the RT task for the RT data storage and update and add it to the 
      ! task list, giving it a nbor list consisting of only the MHD task
      !-----------------------------------------------------------------------------
      allocate (rt_data, link)             ! the internal RT task and its link
      rt_data%mhd_patch => mhd_patch       ! Pointer to mhd_patch
      rt_data%mhd_link  => t_list          ! Pointer to mhd_pathc's link
      link%task         => rt_data         ! connect to RT task
      rt_data_l => link
      rt_data%rt_type = 0
      call rt_data%init                    ! initialise
      call connect_nbors(t_list,link)
      call task_list%append_link(link)     ! insert into the task list
      !---------------------------------------------------------------------
      ! Loop over the ray bundles to construct the solutions and corrections
      !---------------------------------------------------------------------
      ray_dir    => ray_geometry%ray_list%head
      do while (associated(ray_dir))                     ! loop over all directions
        ray_bundle => ray_dir                            
        do while (associated(ray_bundle))                ! loop over all bundles in this direction
          allocate (internal, link)                      ! the RT task and its link
          internal%bundle_ID  = ray_bundle%bundle_ID        ! ray bundle ID
          internal%dir_ID     = ray_bundle%dir_ID           ! ray direction ID
          internal%my_rt_data => rt_data                    ! Pointer to rt_data task
          internal%mhd_patch  => mhd_patch                  ! Pointer to mhd_patch
          internal%rt_type = 1
          if (solver ==  0) then
            internal%solve_internal => solve_internal_feautrier
          else
            internal%solve_internal => solve_internal_integral
          end if
          link%task => internal
          internal_l => link                               ! save pointer to internal task link
          call internal%init                               ! initialize in ray direction
          call connect_nbors(rt_data_l,link)               ! data checks if internal is ready, internal does not do that
          call task_list%append_link(link)                 ! insert into the task list
          ! -------------------------------
          do rt_type = 2,3
            allocate (correction, link)                      ! the RT task and its link
            correction%bundle_ID      = ray_bundle%bundle_ID ! ray bundle ID
            correction%dir_ID         = ray_bundle%dir_ID    ! ray direction ID
            correction%my_rt_data     => rt_data             ! Pointer to rt_data task
            correction%my_rt_internal => internal            ! Pointer to internal task
            correction%mhd_patch      => mhd_patch           ! Pointer to mhd_patch
            correction%rt_type = rt_type
            correction%mhd_link       => t_list              ! Pointer to mhd_pathc's link
            link%task => correction
            call correction%init                             ! initialize in ray direction
            call connect_nbors(link,t_list)                  ! correction checks if MHD is ready, MHD does not do that
            call connect_corr_nbors(internal_l,link, &
                                    correction%my_nb_link)   ! internal checks if correction is ready, corr does not do that
            call task_list%append_link(link)                 ! insert into the task list
          end do
          ! -------------------------------
          
          ray_bundle => ray_bundle%next_bundle     
        end do
        ray_dir => ray_dir%next_dir
      end do
    end SELECT 
    end associate
    t_list => t_list%next                                  ! go to next task
  end do
  call add_correction_nbors(task_list)
  !call add_internal_nbors(task_list)
  list => task_list
  call self%append_list(list)
  call trace_end 
  print *, '=================================='
  print *, 'INIT RT DONE'
  print *, '=================================='
  !$omp end critical (rt_cr)
END SUBROUTINE rt_setup

!> =====================================================================
!> My separate nbor connectors. The original one in link_mod does not do
!> the job right for rt. Second one is for connecting rt_internal and 
!> corresponding corrections. Both checks are deactivated by default. If
!> task becomes a BC one, it activates itself's link.
!> =====================================================================
SUBROUTINE connect_nbors(self, nbor)
  class(link_t), pointer :: self
  class(link_t), pointer :: nbor, link1, link2
! -------------------------------------------------------------------
  allocate(link1, link2)
  link1%task   => self%task
  link1%link   => self
  link1%parent => self%parent
  link1%check = .false.
  call nbor%add_nbor_by_rank(link1)
  link2%task   => nbor%task
  link2%link   => nbor
  link2%parent => nbor%parent     
  call self%add_nbor_by_rank(link2)
END SUBROUTINE connect_nbors

SUBROUTINE connect_corr_nbors(self, nbor, nb_l)
  class(link_t), pointer           :: self
  class(link_t), pointer           :: nbor, link1, link2
  class(link_t), pointer, optional :: nb_l
! -------------------------------------------------------------------
  allocate(link1, link2)
  link1%task   => self%task
  link1%link   => self
  link1%parent => self%parent
  link1%check = .false.
  call nbor%add_nbor_by_rank(link1)
  link2%task   => nbor%task
  link2%link   => nbor
  link2%parent => nbor%parent
  link2%check =.false.     
  call self%add_nbor_by_rank(link2)
  if (present(nb_l)) nb_l => link2
  
END SUBROUTINE connect_corr_nbors

SUBROUTINE add_rt_nbor(self, nbor)
  class(link_t), pointer :: self
  class(link_t), pointer :: nbor, link
! -------------------------------------------------------------------
  allocate(link)
  link%task   => nbor%task
  link%link   => nbor
  link%parent => nbor%parent    
  call self%add_nbor_by_rank(link)
END SUBROUTINE add_rt_nbor

!> =============================================================================
!> Finalize the neighbor lists for upstream correction tasks by connecting 
!> them to appropriate up/down nbors.
!> =============================================================================
SUBROUTINE add_correction_nbors (self)
  class(list_t):: self
  class(link_t),       pointer :: rt_list, my_mhd_nbor, nb_mhd
  class(link_t),       pointer :: rt, my_mhd_l, nb_mhd_nbor, nb_rt_l, nb_mhd_l
  class(rt_correction_t), pointer :: rtc, nb_rt
  integer:: wall(3), n_connected, n_found, n_non_child
  !.....................................................................
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  !$omp critical (list_cr)
  rt_list => self%head
  
  do while (associated(rt_list))                                    ! loop through task_list
    associate (this => rt_list%task)
    SELECT TYPE (this)
    type is (rt_correction_t)
      rt  => rt_list
      rtc => this
      n_connected = 0
      call rtc%find_upstream(nb_mhd, n_found)
      if (n_found < 1) then
        call rtc%set_BC
      else
        do while (associated(nb_mhd))
          nb_mhd_nbor => nb_mhd%link%nbor
          do while (associated(nb_mhd_nbor))
            associate (nb_rt_nb => nb_mhd_nbor%task)
            SELECT TYPE (nb_rt_nb)
            type is (rt_correction_t)
              nb_rt   => nb_rt_nb
              nb_rt_l => nb_mhd_nbor
              
              if ((rtc%dir_ID  == nb_rt%dir_ID).and.&
                  (rtc%rt_type == nb_rt%rt_type)) then
                if (rtc%id/=nb_rt%id) then
                  if (abs(rtc%level-nb_rt%level)<2) then
                    call connect_nbors (nb_rt_l%link, rt)
                    n_connected = n_connected+1
                  end if
                end if
              end if
              
            end SELECT
            end associate
            nb_mhd_nbor => nb_mhd_nbor%next
          end do
          nb_mhd => nb_mhd%next
        end do
      end if
    end SELECT
    end associate
    rt_list => rt_list%next
  end do
  !$omp end critical (list_cr)
  call trace_end
  
END SUBROUTINE add_correction_nbors


!> =============================================================================
!> Finalize the neighbor lists for internal solution tasks by connecting 
!> them to appropriate higher/lower level nbors.
!> =============================================================================
SUBROUTINE add_internal_nbors (self)
  class(list_t):: self
  class(link_t),        pointer :: rt_tasks, rt_tasks2
  real:: me_li(3), me_ui(3), nb_li(3), nb_ui(3)
  !.....................................................................
  call trace_begin('radiative_transfer_mod%add_rt_nbors ')
  !$omp critical (list_cr)
  rt_tasks => self%head
  
  do while (associated(rt_tasks))                                    ! loop through task_list
  
    associate (me => rt_tasks%task)
    SELECT TYPE (me)
    
      type is (rt_internal_t)
      me_li = me%position - 0.5*me%size
      me_ui = me%position + 0.5*me%size
     ! slow, but clean
    
      rt_tasks2 => self%head
    
      do while (associated(rt_tasks2))
      
        associate (nb => rt_tasks2%task)
        SELECT TYPE (nb)
          type is (rt_internal_t)
          if (.not.associated(rt_tasks,rt_tasks2)) then
            if (me%level /= nb%level) then
              nb_li = me%position - 0.5*me%size
              nb_ui = me%position + 0.5*me%size
              
              SELECT case(me%level - nb%level)
              case (1)
                if ((all(me_li<nb_li)).and.(all(me_ui>nb_ui))) &
                     call connect_nbors (rt_tasks, rt_tasks2)
              case (-1)
                if ((all(me_li>nb_li)).and.(all(me_ui<nb_ui))) &
                     call connect_nbors (rt_tasks, rt_tasks2)
                     
              end SELECT
            end if
          end if
          
        end SELECT
        end associate
        rt_tasks2 => rt_tasks2%next
      end do
    
    end SELECT
    end associate    
    rt_tasks => rt_tasks%next
  end do
  !$omp end critical (list_cr)
  call trace_end
  
END SUBROUTINE add_internal_nbors

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
                   nwb => microphysics%n_wl_bins, nr => self%n_ray_dirs)
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
  self%time       = -3!self%mhd_patch%time
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set(bits%root_grid)
  !---------------------------------------------------------------------
  ! The rest of the setup
  !---------------------------------------------------------------------
  allocate(self%name)
  self%name       = 'rt_data'
  self%n_ray_dirs = ray_geometry%ray_dirs
  allocate(self%opacity(nc(1),nc(2),nc(3),   nwb)) ! opacity table
  allocate(self%source (nc(1),nc(2),nc(3),   nwb)) ! Source function table
  self%opacity = 0; self%source = 0
  end associate         
  call trace_end  
END SUBROUTINE init_rt_data

!> =====================================================================
!> Initialise tables for internal RT solution
!> =====================================================================
SUBROUTINE init_internal (self)
  class(rt_internal_t) :: self
  class(ray_list_t),    pointer :: my_rays
  integer :: point, ray
  logical :: fnd
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
  self%time       = -2!elf%mhd_patch%time 
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set (bits%root_grid)
  !if (mhd_patch%is_set(bits%has_BC))    call self%set (bits%has_BC)
  call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .false.)                ! retrieve ray list for this direction
  self%my_rays => my_rays
  allocate(self%name)
  self%name       = 'rt_internal'
  
  !---------------------------------------------------------------------
  ! Allocate memory slots. 
  !---------------------------------------------------------------------
  allocate(self%ds    (my_rays%max_n_points, my_rays%n_rays))
  allocate(self%dtau  (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))
  allocate(self%ex0   (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))

  allocate(self%source(microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))
  allocate(self%qfw   (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))
  allocate(self%qrv   (microphysics%n_wl_bins, my_rays%max_n_points, my_rays%n_rays))

  do ray = 1, self%my_rays%n_rays                                       ! loop over rays to copy coords and ds
    do point = 1, self%my_rays%n_points(ray)
      self%ds(point, ray) = my_rays%ds(point, ray) * self%cell_ds(1) 
    end do
  end do
  self%qfw = 0; self%qrv = 0; self%source = 0; self%dtau = 0
  self%ex0 = 0
  
  end associate
  call trace_end
END SUBROUTINE init_internal

!> =====================================================================
!> Initialise boundary table for upstream values
!> =====================================================================
SUBROUTINE init_correction (self)
  class(rt_correction_t) :: self
  class(ray_list_t),       pointer :: my_rays
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
  self%time       = -1 !self%mhd_patch%time
  self%li        = mhd_patch%li
  self%ui        = mhd_patch%ui
  if (mhd_patch%is_set(bits%static))    call self%set (bits%static)
  if (mhd_patch%is_set(bits%root_grid)) call self%set (bits%root_grid)
 ! allocate (my_rays)
  allocate(self%name)
  
  SELECT case(self%rt_type)
    case(2)
      call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .false.)                ! retrieve ray list for this direction
      self%name       = 'rt_forward'
    case(3)
      call ray_geometry%get_my_rays(my_rays, self%bundle_ID, .true. )                ! retrieve ray list for this direction
      self%name       = 'rt_reverse'
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
END SUBROUTINE init_correction

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
  
  end associate
  deallocate (temp)
  
  self%dtime = self%mhd_patch%dtime
  self%time = self%mhd_patch%time+ max(1e-7,self%dtime*1e-3)
  
  call trace_end
END SUBROUTINE update_rt_data

!> =====================================================================
!> Update internal RT solution
!> =====================================================================
SUBROUTINE update_internal (self)
  class(rt_internal_t) :: self
  integer :: wl, ray, point, p(3)
  character(len=16) :: stime
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_internal '//stime)
  associate (mhd => self%mhd_patch,   cell => self%my_rays%cell_coords, &
             rtd => self%my_rt_data)
  do ray = 1, self%my_rays%n_rays
    do point = 1, self%my_rays%n_points(ray)
      do wl = 1, microphysics%n_wl_bins
        p = cell(:,point,ray)
        self%dtau  (wl,point,ray) = rtd%opacity(p(1),p(2),p(3),wl)*self%ds(point,ray)
        self%source(wl,point,ray) = rtd%source(p(1),p(2),p(3),wl)
      end do
    end do
  end do
  call self%solve_internal
  end associate
  self%dtime = self%mhd_patch%dtime
  self%time = self%mhd_patch%time+ max(2e-7,self%dtime*2e-3)
  call trace_end
END SUBROUTINE update_internal

!> =====================================================================
!> Forward/reverse RT correction
!> =====================================================================
SUBROUTINE update_correction (self)
  class(rt_correction_t) :: self
  real, dimension(microphysics%n_wl_bins) :: I_in
  integer :: pt, ray, ipt, wl, p(3)
  real(kind=8), dimension(:,:,:), pointer :: q
  real(kind=8) :: dir_wgt
  real :: temp
  character(len=16) :: stime
!-----------------------------------------------------------------------
  write (stime,'(1p,g16.6)') self%time
  call trace_begin('rt_t%update_rt_correction '//stime)
  dir_wgt = 2.*pi/self%n_ray_dirs
  self%it = self%mhd_patch%it
  associate (npt   => self%my_rays%n_points,      nwb => microphysics%n_wl_bins,  &
             nrays => self%my_rays%n_rays,        ex0 => self%my_rt_internal%ex0, &
                 s => self%my_rt_internal%source, &
              cell => self%my_rays%cell_coords,   Ifw => self%intensity,  &
               dir => self%dir_ID, q_cool => self%mhd_patch%rt_q, &
               op  => self%my_rt_data%opacity)
  SELECT case(self%rt_type)
    case(2)
      q => self%my_rt_internal%qfw
    case(3)
      q => self%my_rt_internal%qrv
  end SELECT
  
  do ray = 1, self%my_rays%n_rays
    
    I_in = self%i_in(:,ray)                                                          ! get incoming intensity
   ! print *,'i_in:', self%i_in
    q(:,1,ray) = I_in - s(:,1,ray)                                                 ! correction for first (or single) point
    p = cell(:,1,ray)
    temp = self%mhd_patch%cell_gas_temperature(p)
    call self%mhd_patch%update_q( self%time, q(:,1,ray),p,nwb,&
    op(p(1),p(2),p(3),:),temp, dir_wgt)                                      ! update q
    Ifw(:,1,ray) = q(:,1,ray) * ex0(:,1,ray) + s(:,1,ray)                    ! store forward intensity

    if (npt(ray) > 1) then                                                           ! for non-single-point rays only
      do pt = 2, npt(ray)                                                            ! loop over all points in ray
      
        q(:,pt,ray) = q(:,pt,ray) + q(:,pt-1,ray) * ex0(:,pt,ray)              ! apply corrections
        p = cell(:,pt,ray)
        temp = self%mhd_patch%cell_gas_temperature(p)
        call self%mhd_patch%update_q( self%time, q(:,pt,ray),p,nwb,&
        op(p(1),p(2),p(3),:),temp, dir_wgt)                                        ! update q
        Ifw(:,pt,ray) = q(:,pt,ray) + s(:,pt,ray)                                  ! store forward intensities

      end do
      
    end if
    
    
  end do
  if (any(isnan(Ifw)).or.any(isnan(ex0))) then 
  !print *,'=================== Ifw ================'
  !print *,Ifw
  !print *,'=================== ex0 ================'
  !print *,exo
  !print *,'=================== qfw ================'
  !print *, qfw
  !print *,'=================== s ================'
  !print *, s

  error stop "update corr is a bitch"
  end if
  end associate
  
  self%dtime = self%mhd_patch%dtime
  self%time = self%mhd_patch%time+ max(3e-7,self%dtime*3e-3)
  call trace_end
END SUBROUTINE update_correction

!> =====================================================================
!> Solution of RTE using the Integral formulation. Forward and reverse
!> incoming radiation is ignored. dtau is stored for every point to later
!> apply the corrections. Iwf and Irv are stored in case of overlapping 
!> patches.
!> =====================================================================

! FIXME: storing Ifw and Irv should be ooptional, depending whether the
! patches are moving in the simulation, i.e. they might overlap each other

SUBROUTINE solve_internal_integral (self)
  class(rt_internal_t) :: self
  real(kind=8), dimension(:,:,:), allocatable :: ex1, ex2, ds1, d1s, d2s
  real(kind=8), dimension(:)    , allocatable :: ctau
  integer :: pt, ipt, it, ray
  !---------------------------------------------------------------------
  call trace_begin('rt_t%update_internal%solve_internal_integral')

  associate (npt   => self%my_rays%n_points,nwb => microphysics%n_wl_bins, &
             nrays => self%my_rays%n_rays, dtau => self%dtau,   &
               ex0 => self%ex0,               s => self%source, & 
               qfw => self%qfw,             qrv => self%qrv)
               
  allocate(ex1 (nwb, self%my_rays%max_n_points, nrays))
  allocate(ex2 (nwb, self%my_rays%max_n_points, nrays))
  allocate(ds1 (nwb, self%my_rays%max_n_points, nrays))
  allocate(d1s (nwb, self%my_rays%max_n_points, nrays))
  allocate(d2s (nwb, self%my_rays%max_n_points, nrays))
  allocate(ctau(nwb))

  do ray = 1, nrays
    if (npt(ray) > 1) then
      do pt=2,npt(ray)
        ex0(:,pt,ray) = dexp(-dtau(:,pt,ray))                                 ! first exponential coeff
        ex1(:,pt,ray) = 1.-ex0(:,pt,ray)                                     ! second exponential coeff
        ex2(:,pt,ray) = ex1(:,pt,ray)-dtau(:,pt,ray)*ex0(:,pt,ray)           ! third exponential coeff
        ds1(:,pt,ray) = (s(:,pt,ray)-s(:,pt-1,ray))/dtau(:,pt,ray)           ! half-point centered dS/dtau
      end do
      !---------------------------------------------------------------------
      !  dS/dtau and d2S/dtau2 
      !---------------------------------------------------------------------
      do pt=2,npt(ray)-1
        ctau          = 1. / (dtau(:,pt+1,ray) + dtau(:,pt,ray))          
        d1s(:,pt,ray) = &
        (ds1(:,pt+1,ray)*dtau(:,pt,ray)+ds1(:,pt,ray)*dtau(:,pt+1,ray))*ctau ! dS/dtau, 2nd. order accurate
        d2s(:,pt,ray) = (ds1(:,pt+1,ray)-ds1(:,pt,ray))*2.*ctau              ! d2S/dtau2, 2nd. order accurate
      end do

      d2s(:,1,ray) = d2s(:,2,ray)                                            ! assume constant d2S/dtau2
      d1s(:,1,ray) = d1s(:,2,ray) - dtau(:,2,ray) * d2s(:,2,ray)             ! extrapolate to 1
      d2s(:,npt(ray),ray) = d2s(:,npt(ray)-1,ray)                            ! assume constant d2S/dtau2
      d1s(:,npt(ray),ray) = d1s(:,npt(ray)-1,ray) + &
                            dtau(:,npt(ray) ,ray) * d2s(:,npt(ray),ray)      ! extrapolate to n
      !---------------------------------------------------------------------
      !  Internal forward and reverse formal solutions, skipping first and last points
      !---------------------------------------------------------------------
      do pt=2,npt(ray)
        ipt = npt(ray)-pt+1
        it  = ipt+1
        qfw(:,pt,ray) =-d1s(:, pt,ray)*ex1(:,pt,ray)+d2s(:, pt,ray)*ex2(:,pt,ray)
        qrv(:,ipt,ray)= d1s(:,ipt,ray)*ex1(:,it,ray)+d2s(:,ipt,ray)*ex2(:,it,ray)
      end do

    else
      ex0(:,npt(ray),ray) = exp(-dtau(:,npt(ray),ray))                  ! exp_dtau is still needed
    end if
  end do
  if (any(isnan(ex0)).or.any(isnan(qfw)).or.any(isnan(qrv))) then
  !print *,'=================== ex0 ================'
  !print *,ex0
  !print *,'=================== qfw ================'
  !print *,qfw
  !print *,'=================== qrv ================'
  !print *, qrv
  !print *,'=================== s ================'
  !print *, s

  error stop "internal is a bitch"
  end if
  deallocate(ex1, ex2, ds1, d1s, d2s, ctau)
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

END MODULE radiative_transfer_mod
