MODULE radiative_transfer_mod
  USE bits_mod
  USE eos_mod
  USE io_mod
  USE list_mod
  USE mhd_mod
  USE microphysics_mod
  USE omp_mod
  USE opacity_mod
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
  implicit none
  private
  PUBLIC rt_setup, connect_nbors, add_upstream_nbors
CONTAINS

SUBROUTINE rt_setup (mhd_link, mp_link, rt_task_list, &
                             solver, rt_freq, store_rt)
  class(link_t),       pointer :: mhd_link
  class(task_list_t),  pointer :: rt_task_list                          ! rt task list
  class(list_t),       pointer :: list
  class(link_t),       pointer :: mp_link, link, rt_data_l, internal_l
  class(mhd_t),        pointer :: mhd_patch
  class(microphysics_t), pointer :: mp
  class(rt_data_t),    pointer :: rt_data
  class(rt_internal_t),pointer :: internal
  class(ray_dirs_t),   pointer :: ray_dir
  class(ray_bndl_t),   pointer :: ray_bundle 
  integer :: n_rt, rt_type, a, b, solver, rt_freq
  logical :: store_rt

  !---------------------------------------------------------------------
  call trace_begin('radiative_transfer::rt_setup')

  !---------------------------------------------------------------------
  ! Get the right type of the task (mhd_t)
  !---------------------------------------------------------------------
  associate (mhd => mhd_link%task)
  SELECT TYPE(mhd)
  class is (mhd_t)
    mhd_patch => mhd
    n_rt = 0
    !-------------------------------------------------------------------
    ! Construct the RT task for the RT data storage and update and add
    ! it to the task list, giving it a nbor list consisting of only
    ! the MHD task.
    !-------------------------------------------------------------------
    allocate (rt_data, link)                                            ! the internal RT task and its link
    rt_data%mhd_patch => mhd_patch                                      ! Pointer to mhd_patch
    rt_data%mhd_link  => mhd_link                                       ! Pointer to mhd_patch's link
    link%task         => rt_data                                        ! connect to RT task
    rt_data_l => link                                                   ! save the link
    rt_data%link => link                                                ! link to oneself link in the task list
    rt_data%rt_type = 0                                                 ! set the rt type
    if (store_rt) rt_data%do_output = .true.
    rt_data%n_ray_dirs = ray_geometry%n_ray_dirs                        ! self explanatory
    rt_data%rt_freq    =  rt_freq                                       ! execution frequency, relative to MHD task
    call rt_data%init                                                   ! initialise
    n_rt = n_rt+1
    call connect_nbors(mhd_link,link)                                   ! connect rt_data and mhd task
    call rt_task_list%append_link(link)                                 ! insert rt_data into the task list
    !-----------------------------------------------------------------
    ! Loop over the ray bundles to construct the solutions
    !-----------------------------------------------------------------
    ray_dir    => ray_geometry%head%ray_list
    associate (mph => mp_link%task)
    select type (mph)
    class is (microphysics_t)
      mp => mph
    end select
    end associate
    do while (associated(ray_dir))                                      ! loop over all directions
      ray_bundle => ray_dir%bundle
      do while (associated(ray_bundle))                                 ! loop over all bundles in this direction
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
          allocate (internal, link)                                     ! the RT task and its link
          internal%bundle_ID  =  ray_bundle%ID                          ! ray bundle ID
          internal%dir_ID     =  ray_dir%ID                             ! ray direction ID
          internal%my_rt_data => rt_data                                ! Pointer to rt_data task
          internal%mhd_patch  => mhd_patch                              ! Pointer to mhd_patch
          internal%rt_type    =  rt_type                                ! set the rt type
          if (store_rt) internal%do_output = .true.
          internal%n_ray_dirs =  ray_geometry%head%n_ray_dirs
          internal%solver     =  solver                                 ! set the rt solver type
          internal%rt_freq    =  rt_freq                                ! integer frequency, relative to MHD task
          internal%mhd_link   => mhd_link                               ! Pointer to mhd_patch's link
          internal%link       => link                                   ! convenience
          link%task => internal                                         ! connect to RT task
          internal%nt = mhd_patch%nt
          internal%mp => mp
          call internal%init                                            ! initialize in ray direction
          n_rt = n_rt+1
          call connect_nbors(link,mp_link)                              ! correction checks if forge is ready, forge does not do that
          call connect_corr_nbors(rt_data_l,link, &
                                  internal%my_nb_link)                  ! internal checks if correction is ready, corr does not do that
          call internal%print
          call rt_task_list%append_link(link)                            ! insert into the task list
        end do
        ! -------------------------------
       ray_bundle => ray_bundle%next                                    ! take the next ray bundle
      end do
      ray_dir => ray_dir%next                                           ! go to the next ray direction
    end do
  end SELECT
  end associate

END SUBROUTINE rt_setup

!=======================================================================
!> My separate nbor connectors. The original one in link_mod does not do
!> the job right for rt. Second one is for connecting rt_internal and
!> corresponding corrections. Both checks are deactivated by default. If
!> task becomes a BC one, it activates itself's link.
!=======================================================================
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

!=======================================================================
!> My separate nbor connectors. This version deactivates both links and
!> saves the link of the correction task for it to be altered upon need
!=======================================================================
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

!===============================================================================
!> Finalize the neighbor lists for upstream correction tasks by connecting
!> them to appropriate up/down nbors.
!===============================================================================
SUBROUTINE add_upstream_nbors (self)
  class(list_t):: self
  class(link_t), pointer :: rt_list     !, my_mhd_nbor, nb_mhd, nb_mhd_prev
  class(link_t), pointer :: rt, try     !my_mhd_l, nb_mhd_nbor, nb_rt_l, nb_mhd_l
  class(task_t), pointer :: this, that
  class(rt_internal_t), pointer :: ithis, ithat
  integer:: wall(3)
  !.....................................................................
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  !$omp critical (list_cr)
  rt_list => self%head

  do while (associated(rt_list))                                        ! loop through task_list
    this => rt_list%task
    SELECT TYPE (this) 
    type is (rt_internal_t)                                             ! select internal tasks
      ithis => this
      if (ithis%has_OB()) then
        call ithis%init_obc
        go to 9
      end if
      try => self%head
      do while (associated(try))
        that => try%task
        if (that%id /= this%id) then
          SELECT TYPE (that)
          class is (rt_internal_t)
            ithat => that
            if ((ithat%dir_id == ithis%dir_id).and. &
                (ithat%rt_type==ithis%rt_type)) then
              select case(this%rt_type)                                     ! select the right upstream wall
              case(2)
                wall = this%bundle%lwall
              case(3)
                wall = this%bundle%uwall
              end select
              if (abs(ithis%level-ithat%level)<2) then
                if (ithis%the_right_wall(try,wall)) then
                  if (ithis%bundle_connects_to(ithat)) then
                    call connect_nbors (try, rt_list)                         ! connect these correction tasks
                    call ithis%add_upstream(try)
                    go to 9
                  end if
                end if
              else
                error stop "radiative_transfer_mod::add_upstream_nbors:: (1)"
              end if
            end if
          END SELECT
        end if
        try => try%next
      end do
    end SELECT
9 rt_list => rt_list%next                                             ! move to next rt task
  end do
  !$omp end critical (list_cr)
  call trace_end
 
END SUBROUTINE add_upstream_nbors


!===============================================================================
!> Finalize the neighbor lists for internal solution tasks by connecting
!> them to appropriate higher/lower level nbors.
!> TO DO: NEEDS LOGICS, to be done: internal alteration, ray cutting, etc
!===============================================================================
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
