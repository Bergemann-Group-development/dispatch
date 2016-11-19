!> =============================================================================
!> Module for setting up the radiative transfer, and the relations between the
!> RT tasks and ray systems.
!> =============================================================================
MODULE radiative_transfer_mod
  USE bits_mod
  USE experiment_mod
  USE io_mod
  USE list_mod
  USE mesh_mod
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
  class(mesh_t),       pointer :: mesh
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
  logical, save :: store_rt = .false.
  namelist /rt_params/ RT_solver, cr_heating, rt_freq, store_rt
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
  
  mesh => mhd_patch%mesh(1) 
  select type(mesh)
    type is(Cartesian_mesh)
      call ray_geometry%init(1, mhd_patch%li,mhd_patch%ui)    ! FIXME: send cell wall sizes and dimensions as well
    type is(spherical_mesh)
      call ray_geometry%init(2, mhd_patch%li,mhd_patch%ui)    ! FIXME: send cell wall sizes and dimensions as well
    type is(cylindrical_mesh)
      call ray_geometry%init(3, mhd_patch%li,mhd_patch%ui)    ! FIXME: send cell wall sizes and dimensions as well
  end select
  nullify(mesh)
  
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
      call mhd_rt_patch%init_rt                                         ! NOTE: isn't it already initialised in experiment_mod?
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
      if (store_rt) rt_data%do_output = .true.
      rt_data%n_ray_dirs = ray_geometry%n_ray_dirs
      rt_data%rt_freq    =  rt_freq                                     ! integer frequency, relative to MHD task
      rt_data%position = mhd_patch%position
      call rt_data%init                                                 ! initialise
      n_rt = n_rt+1
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
            internal%bundle     => ray_bundle
            internal%my_rt_data => rt_data                              ! Pointer to rt_data task
            internal%mhd_patch  => mhd_patch                            ! Pointer to mhd_patch
            internal%rt_type    =  rt_type                              ! set the rt type
            internal%n_ray_dirs =  ray_geometry%head%n_ray_dirs
            internal%solver     =  solver                               ! set the rt solver type
            if (store_rt) internal%do_output = .true.
            internal%rt_freq    =  rt_freq                              ! integer frequency, relative to MHD task
            internal%mhd_link   => t_list                               ! Pointer to mhd_patch's link
            internal%link       => link                                 ! convenience
            internal%BC_direction = mhd_rt_patch%BC_direction           ! pass on
            link%task => internal                                       ! connect to RT task
            internal%nt = mhd_patch%nt
            internal%position = mhd_patch%position
            call internal%init                                          ! initialize in ray direction
            n_rt = n_rt+1
            if (internal%chain_start()) then                            ! bundle chain start?
              if (associated(mhd_rt_patch%first_internal)) then
                internal%next => mhd_rt_patch%first_internal%next       ! point to the next internal
              end if
              mhd_rt_patch%first_internal => internal                   ! prepend the new one
              write (io_unit%log,*) &
                'bundle chain: id, axis =', internal%id, internal%bundle%axis, &
                real(internal%bundle%theta*radian), real(internal%bundle%phi*radian), internal%bundle%lwall, internal%bundle%uwall
            end if
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
  call add_downstream (task_list)                                       ! connect the downstrem tasks
  !call add_internal_nbors(task_list)
  call task_list%print
  call trace_end
  if (io%master) then
    print *, '=================================='
    print *, 'INIT RT DONE'
    print *, '=================================='
  end if
END SUBROUTINE rt_setup

!> =============================================================================
!> Finalize the neighbor lists for upstream correction tasks by connecting
!> them to appropriate up/down nbors.
!> An rt_internal_t has integer dir_IDs that must match
!> A task has an double precision position
!> =============================================================================
SUBROUTINE add_downstream (self)
  class(list_t):: self
  class(link_t),       pointer :: rt_list
  class(link_t),       pointer :: rt, try
  class(rt_internal_t),pointer :: rtc, nb_rt
  integer:: wall(3), n_connected, n_found, n_non_child
  class(task_t), pointer :: this, that
  class(rt_internal_t), pointer :: ithis, ithat
  !.....................................................................
  call trace_begin('radiative_transfer_mod%add_correction_nbors ')
  !$omp critical (list_cr)
  rt_list => self%head
  do while (associated(rt_list))                                        ! loop through task_list
    this => rt_list%task
    SELECT TYPE (this)
    class is (rt_internal_t)                                             ! select correction tasks
      try => self%head
      do while (associated(try))
        that => try%task
        SELECT TYPE (that)
        class is (rt_internal_t)
          ithis => this
          ithat => that
          if (ithis%bundle_connects_to(ithat)) then
            this%downstream => that
            that%upstream => this
            if (io%verbose > 1) &
              write(io_unit%log,'(1x,a,i5,a,i5,3x,a,2(3f6.1,3x),a,3i3)') &
                'downstream nbor of',this%id,' is',this%downstream%id, &
                'at positions:', this%mhd_patch%position/this%mhd_patch%ds, &
                                 that%mhd_patch%position/this%mhd_patch%ds, &
                'dstep:', this%bundle%dstep
            go to 9
          end if
        END SELECT
        try => try%next
      end do
      if (io%verbose > 2) &
        write(io_unit%log,'(1x,a,i5,3x,a,3f6.1,3x,a,3i3,3x,a)') &
         'downstream nbor of', this%id, &
         'at pos:', this%mhd_patch%position/this%mhd_patch%ds, &
         'dstep:', this%bundle%dstep, 'not found'
    end SELECT
9   rt_list => rt_list%next                                             ! move to next rt task
  end do
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE add_downstream

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
                !if ((all(me_li<nb_li)).and.(all(me_ui>nb_ui))) &
                     !call connect_nbors (rt_tasks, rt_tasks2)
              case (-1)
                !if ((all(me_li>nb_li)).and.(all(me_ui<nb_ui))) &
                     !call connect_nbors (rt_tasks, rt_tasks2)
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
