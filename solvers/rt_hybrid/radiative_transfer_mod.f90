MODULE radiative_transfer_mod
  USE bits_mod
  USE experiment_mod
  USE io_mod
  USE list_mod
  USE mhd_mod
  USE microphysics_mod
  USE omp_mod
  USE patch_mod
  USE ray_geometry_mod
  USE rt_mod
  USE rt_boundary_mod
  USE rt_corrections_mod
  USE rt_data_mod
  USE rt_internal_mod
  USE scaling_mod
  USE task_mod
  USE task_list_mod
  USE trace_mod
  USE units_mod
  implicit none

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
  integer:: solver, rt_type
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
    call ray_geometry%init(mhd_patch%li,mhd_patch%ui)                   ! FIXME: send cell wall sizes and dimensions as well 

    print *, 'Ray geometry: n_bundles, n_ray_dirs:',&
              ray_geometry%n_bundles, ray_geometry%ray_dirs
              
    call microphysics%init
    call rt_boundary%init
    
    rt_boundary%ubox = mhd_patch%box                                    ! FIXME : add bndry patches
    rt_boundary%lbox = 0.0
    
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
      !-----------------------------------------------------------------
      ! Construct the RT task for the RT data storage and update and add
      ! it to the task list, giving it a nbor list consisting of only
      ! the MHD task.
      !-----------------------------------------------------------------
      allocate (rt_data, link)                                          ! the internal RT task and its link
      rt_data%mhd_patch => mhd_patch                                    ! Pointer to mhd_patch
      rt_data%mhd_link  => t_list                                       ! Pointer to mhd_patch's link
      link%task         => rt_data                                      ! connect to RT task
      rt_data_l => link                                                 ! save the link
      rt_data%rt_type = 0                                               ! set the rt type
      call rt_data%init                                                 ! initialise
      call connect_nbors(t_list,link)                                   ! connect rt_data and mhd task
      call task_list%append_link(link)                                  ! insert rt_data into the task list
      !-----------------------------------------------------------------
      ! Loop over the ray bundles to construct the solutions and
      ! corrections
      !-----------------------------------------------------------------
      ray_dir    => ray_geometry%ray_list%head
      do while (associated(ray_dir))                                    ! loop over all directions
        ray_bundle => ray_dir                            
        do while (associated(ray_bundle))                               ! loop over all bundles in this direction
          allocate (internal, link)                                     ! the RT task and its link
          internal%bundle_ID  = ray_bundle%bundle_ID                    ! ray bundle ID
          internal%dir_ID     = ray_bundle%dir_ID                       ! ray direction ID
          internal%my_rt_data => rt_data                                ! Pointer to rt_data task
          internal%mhd_patch  => mhd_patch                              ! Pointer to mhd_patch
          internal%rt_type = 1                                          ! set the rt type
          internal%solver  = solver                                     ! set the rt solver type
          link%task => internal                                         ! connect to RT task
          internal_l => link                                            ! save pointer to internal task link
          call internal%init                                            ! initialize in ray direction
          call connect_nbors(rt_data_l,link)                            ! data checks if internal is ready, internal does not do that
          call task_list%append_link(link)                              ! insert into the task list
          ! ------------------------------------------------------------
          ! Loop over forward and reverse corrections
          ! ------------------------------------------------------------
          do rt_type = 2,3
            allocate (correction, link)                                 ! the RT task and its link
            correction%bundle_ID      = ray_bundle%bundle_ID            ! ray bundle ID
            correction%dir_ID         = ray_bundle%dir_ID               ! ray direction ID
            correction%my_rt_data     => rt_data                        ! Pointer to rt_data task
            correction%my_rt_internal => internal                       ! Pointer to internal task
            correction%mhd_patch      => mhd_patch                      ! Pointer to mhd_patch
            correction%rt_type = rt_type                                ! set the rt type
            correction%solver  = solver                                 ! set the rt solver type
            correction%mhd_link       => t_list                         ! Pointer to mhd_patch's link
            link%task => correction                                     ! connect to RT task
            call correction%init                                        ! initialize in ray direction
            call connect_nbors(link,t_list)                             ! correction checks if MHD is ready, MHD does not do that
            call connect_corr_nbors(internal_l,link, &
                                    correction%my_nb_link)              ! internal checks if correction is ready, corr does not do that
            call task_list%append_link(link)                            ! insert into the task list
          end do
          ! -------------------------------
          
          ray_bundle => ray_bundle%next_bundle                          ! take the next ray bundle
        end do
        ray_dir => ray_dir%next_dir                                     ! go to the next ray direction
      end do
    end SELECT 
    end associate
    t_list => t_list%next                                               ! go to the next mhd task
  end do
  call add_correction_nbors(task_list)                                  ! connect the corrections
  !call add_internal_nbors(task_list)
  list => task_list
  call self%append_list(list)                                           ! add rt task list to the main task list
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
      n_found = 0
      call rtc%find_upstream(nb_mhd, n_found)
      if (n_found < 1) then
        call rtc%set_BC
        call rtc%add_src()
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
                    call rtc%add_src(nb_rt_l%link)
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

END MODULE radiative_transfer_mod
