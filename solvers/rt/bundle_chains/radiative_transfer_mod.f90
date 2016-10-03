!> =============================================================================
!> Module for setting up the radiative transfer, and the relations between the
!> RT tasks and ray systems.
!> =============================================================================
MODULE radiative_transfer_mod
  USE bits_mod
  USE experiment_mod
  USE io_mod
  USE list_mod
  USE mhd_mod
  USE mhd_rt_mod
  USE microphysics_mod
  USE omp_mod
  USE patch_mod
  USE ray_geometry_v3_mod
  USE rt_mod
  USE rt_boundary_mod
  USE rt_data_mod
  USE rt_internal_mod
  USE scaling_mod
  USE task_mod
  USE task_list_mod
  USE trace_mod
  USE units_mod
  implicit none
  private
  PUBLIC rt_setup
CONTAINS

!> =============================================================================
!> The only procedure callable from the outside
!> =============================================================================
SUBROUTINE rt_setup (self)
  class(list_t)                :: self
  class(task_list_t),  pointer :: task_list                  ! local task list
  class(list_t),       pointer :: list
  class(link_t),       pointer :: t_list, link, rt_data_l, internal_l
  class(mhd_t),        pointer :: mhd_patch
  class(mhd_rt_t),     pointer :: mhd_rt_patch
  class(rt_data_t),    pointer :: rt_data
  class(rt_internal_t),pointer :: internal
  class(ray_dirs_t),   pointer :: ray_dir
  class(ray_bndl_t),   pointer :: ray_bundle
  logical, save                :: first_time=.true.
  character(len=16)            :: RT_solver ='integral'      ! choice of RT solver
  integer :: n_rt
  real,    save :: cr_heating=4.4e-27
  integer, save :: solver, rt_type, a, b
  integer, save :: rt_freq=1
  namelist /rt_params/ RT_solver, cr_heating, rt_freq
  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer::rt_setup')
  associate (mhd => self%head%task)
  SELECT TYPE(mhd)
  class is (mhd_t)
    mhd_patch => mhd
  class default
    error stop '(1) in rt_setup'
  end SELECT
  end associate
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, rt_params)
    if (io%master) write (*, rt_params)
    SELECT case(RT_solver)
    case('feautrier')
      solver = 0
    case('integral')
      solver = 1
    end SELECT
  end if
  call ray_geometry%init(1, mhd_patch%li,mhd_patch%ui)    ! FIXME: send cell wall sizes and dimensions as well
  ! IMPORTANT: change first constant into mhd_patch%mesh%mesh_type after merge with curvilinear coords
  write (io_unit%log,*) 'Ray geometry: n_ray_dirs:', ray_geometry%head%n_ray_dirs
  call microphysics%init
  call rt_boundary%init
  rt_boundary%ubox = mhd_patch%box                                    ! FIXME : add bndry patches
  rt_boundary%lbox = 0.0
  allocate(task_list)
  call task_list%list_t%init ('rt task list')
  t_list => self%head
  do while (associated(t_list))
    !-------------------------------------------------------------------
    ! Get the right type of the task for the init call (mhd_rt_t)
    !-------------------------------------------------------------------
    associate (mhd => t_list%task)
    SELECT TYPE(mhd)
    class is (mhd_rt_t)
      mhd_rt_patch => mhd
      allocate (rt_data, link)                                          ! the internal RT task and its link
      mhd_rt_patch%rt_data => rt_data
      call mhd_rt_patch%init_rt
    END SELECT
    !-------------------------------------------------------------------
    ! Get the right type of the task to pick up other properties (mhd_t)
    !-------------------------------------------------------------------
    SELECT TYPE(mhd)
    class is (mhd_t)
      mhd_patch => mhd
      n_rt = 0
      !-----------------------------------------------------------------
      ! Construct the RT task for the RT data storage and update and add
      ! it to the task list, giving it a nbor list consisting of only
      ! the MHD task.
      !-----------------------------------------------------------------
      rt_data%mhd_patch => mhd_patch                                    ! Pointer to mhd_patch
      rt_data%mhd_link  => t_list                                       ! Pointer to mhd_patch's link
      link%task         => rt_data                                      ! connect to RT task
      rt_data_l => link                                                 ! save the link
      rt_data%link => link
      rt_data%rt_type = 0                                               ! set the rt type
      rt_data%n_ray_dirs = ray_geometry%n_ray_dirs
      rt_data%rt_freq    =  rt_freq                                     ! integer frequency, relative to MHD task
      call rt_data%init                                                 ! initialise
      n_rt = n_rt+1
      call connect_nbors(t_list,link)                                   ! connect rt_data and mhd task
      call task_list%append_link(link)                                  ! insert rt_data into the task list
      !-----------------------------------------------------------------
      ! Loop over the ray bundles to construct the solutions
      !-----------------------------------------------------------------
      ray_dir    => ray_geometry%head%ray_list
      do while (associated(ray_dir))                                    ! loop over all directions
        ray_bundle => ray_dir%bundle
        do while (associated(ray_bundle))                               ! loop over all bundles in this direction
          ! ------------------------------------------------------------
          ! Loop over forward and/or reverse directions
          ! ------------------------------------------------------------
          select case(ray_dir%directionality(ray_bundle%axis))
          case(1)
          a = 2; b = 3
          case(2)
          a = 2; b = 2
          case(3)
          a = 3; b = 3
          end select
          do rt_type = a,b
            allocate (internal, link)                                   ! the RT task and its link
            internal%bundle_ID  =  ray_bundle%ID                        ! ray bundle ID
            internal%dir_ID     =  ray_dir%ID                           ! ray direction ID
            internal%my_rt_data => rt_data                              ! Pointer to rt_data task
            internal%mhd_patch  => mhd_patch                            ! Pointer to mhd_patch
            internal%rt_type    =  rt_type                              ! set the rt type
            internal%n_ray_dirs =  ray_geometry%head%n_ray_dirs
            internal%solver     =  solver                               ! set the rt solver type
            internal%rt_freq    =  rt_freq                              ! integer frequency, relative to MHD task
            internal%mhd_link   => t_list                               ! Pointer to mhd_patch's link
            internal%link       => link                                 ! convenience
            internal%BC_direction = mhd_rt_patch%BC_direction           ! pass on
            link%task => internal                                       ! connect to RT task
            internal%nt = mhd_patch%nt
            call internal%init                                          ! initialize in ray direction
            n_rt = n_rt+1
            call connect_nbors(link,t_list)                             ! correction checks if MHD is ready, MHD does not do that
            call connect_corr_nbors(rt_data_l,link, &
                                    internal%my_nb_link)                ! internal checks if correction is ready, corr does not do that
            call internal%print
            call task_list%append_link(link)                            ! insert into the task list
          end do
          ! -------------------------------
          ray_bundle => ray_bundle%next                                 ! take the next ray bundle
        end do
        ray_dir => ray_dir%next                                         ! go to the next ray direction
      end do
    end SELECT
    end associate
    write (io_unit%log,*) 'mhd_rt_t%id, n_nbors:', mhd_patch%id, mhd_patch%n_nbors, n_rt
    t_list => t_list%next                                               ! go to the next mhd task
  end do
  call add_upstream_nbors(task_list)                                    ! connect the corrections
  !call add_internal_nbors(task_list)
  call task_list%print
  list => task_list
  call self%append_list(list)                                           ! add rt task list to the main task list
  call trace_end
  if (io%master) then
    print *, '=================================='
    print *, 'INIT RT DONE'
    print *, '=================================='
  end if
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

!> =====================================================================
!> My separate nbor connectors. This version deactivates both links and
!> saves the link of the correction task for it to be altered upon need
!> =====================================================================
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

!> =============================================================================
!> Finalize the neighbor lists for upstream correction tasks by connecting
!> them to appropriate up/down nbors.
!> =============================================================================
SUBROUTINE add_upstream_nbors (self)
  class(list_t):: self
  class(link_t),       pointer :: rt_list, my_mhd_nbor, nb_mhd, nb_mhd_prev
  class(link_t),       pointer :: rt, my_mhd_l, nb_mhd_nbor, nb_rt_l, nb_mhd_l
  class(rt_internal_t),pointer :: rtc, nb_rt
  integer:: wall(3), n_connected, n_found, n_non_child
  !.....................................................................
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  !$omp critical (list_cr)
  rt_list => self%head
  do while (associated(rt_list))                                        ! loop through task_list
    associate (this => rt_list%task)
    SELECT TYPE (this)
    type is (rt_internal_t)                                             ! select correction tasks
      rt  => rt_list
      rtc => this
      n_connected = 0
      n_found = 0
      nullify(nb_mhd)
      call rtc%find_upstream(nb_mhd, n_found)                           ! find upstream mhd nbors for this task
      if (n_found < 1) then
        call rtc%init_obc                                               ! if no upstream found - initialize it as an outer boundary
      else
        do while (associated(nb_mhd))                                   ! loop over upstream mhd patches
          nb_mhd_nbor => nb_mhd%link%nbor
          do while (associated(nb_mhd_nbor))
            associate (nb_rt_nb => nb_mhd_nbor%task)
            SELECT TYPE (nb_rt_nb)
            type is (rt_internal_t)                                     ! select upstream rt correction tasks
              nb_rt   => nb_rt_nb
              nb_rt_l => nb_mhd_nbor
              if ((rtc%dir_ID  == nb_rt%dir_ID).and.&                   ! check if they are the needed ones
                  (rtc%rt_type == nb_rt%rt_type)) then
                if (rtc%id/=nb_rt%id) then
                  if (abs(rtc%level-nb_rt%level)<2) then                ! only +/- one level difference is concerned
                    call connect_nbors (nb_rt_l%link, rt)               ! connect these correction tasks
                    n_connected = n_connected+1
                  end if
                end if
              end if
            end SELECT
            end associate
            nb_mhd_nbor => nb_mhd_nbor%next                             ! move to next upstream mhd nbor's nbor
          end do
          nb_mhd_prev => nb_mhd
          nb_mhd => nb_mhd%next                                         ! move to next upstream nbor
          deallocate(nb_mhd_prev)
        end do
      end if
    end SELECT
    end associate
    rt_list => rt_list%next                                             ! move to next rt task
  end do
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE add_upstream_nbors

!> =============================================================================
!> Finalize the neighbor lists for internal solution tasks by connecting
!> them to appropriate higher/lower level nbors.
!> TO DO: NEEDS LOGICS, to be done: internal alteration, ray cutting, etc
!> =============================================================================
SUBROUTINE add_internal_nbors (self)
  class(list_t):: self
  class(link_t),        pointer :: rt_tasks, rt_tasks2
  real:: me_li(3), me_ui(3), nb_li(3), nb_ui(3)
  !.....................................................................
  call trace_begin('radiative_transfer_mod%add_rt_nbors ')
  !$omp critical (list_cr)
  rt_tasks => self%head
  do while (associated(rt_tasks))                                      ! loop through task_list
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
