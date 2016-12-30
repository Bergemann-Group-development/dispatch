MODULE microphysics_factory_mod
  USE bits_mod
  USE eos_mod
  USE experiment_mod
  USE io_mod
  USE list_mod
  USE mesh_mod
  USE mhd_mod
  USE microphysics_mod
  USE omp_mod
  USE patch_mod
  USE radiative_transfer_mod
  USE ray_geometry_v3_mod
  USE rt_boundary_mod
  USE scaling_mod
!  USE source_function_mod
  USE task_mod
  USE task_list_mod
  USE trace_mod
  USE units_mod  
  implicit none
  private  
  PUBLIC Microphysics_Factory
  
CONTAINS

!> =============================================================================
!> Microphysics factory, which spawns microphysics objects. Such objects can 
!> contain radiative transfer, chemistry, etc. tasks. 
!> =============================================================================
SUBROUTINE Microphysics_Factory (self)
  class(list_t)               :: self
  class(microphysics_t), pointer :: mp
  class(task_list_t), pointer :: rt_task_list, task_list                ! task list to hold rt_internal tasks and mp tasks
  class(list_t),      pointer :: list
  class(link_t),      pointer :: t_list, link, mp_l
  class(mhd_t),       pointer :: mhd_patch
  class(mesh_t),      pointer :: mesh
  logical, save :: first_time=.true.
  logical, save :: do_rt = .false.
  logical, save :: do_chemistry = .false.
  real (kind=8), save :: cr_heating=4.4e-27
  integer, save :: solver
  integer, save :: rt_freq=1
  logical, save :: store_rt = .false.
  character(len=16), save :: RT_solver ='integral'                      ! choice of RT solver
  namelist /microphysics_params/ do_rt, do_chemistry
  namelist /rt_params/ RT_solver, cr_heating, rt_freq, store_rt
  !---------------------------------------------------------------------
  call trace_begin('microphysics_t::microphysics_factory')

  !---------------------------------------------------------------------
  ! Read input parameters
  !---------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, microphysics_params)
    if (io%master) write (*, microphysics_params)
  end if
  
  !---------------------------------------------------------------------
  ! Initialize what is needed
  !---------------------------------------------------------------------
  allocate(task_list)
  call task_list%list_t%init ('rt task list')
  !---------------------------------------------------------------------
  ! If do_rt then initialize ray geometry, allocate rt task list for
  ! internal rt tasks, initialize rt_boundary, read rt parameters from
  ! the input file.
  !---------------------------------------------------------------------
  if (do_rt) then
    associate (mhd => self%head%task)
    SELECT TYPE(mhd)
      class is (mhd_t)
      mhd_patch => mhd
    end SELECT
    end associate
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
  write (io_unit%log,*) 'Ray geometry: n_ray_dirs:', &
  ray_geometry%head%n_ray_dirs
  allocate(rt_task_list)
  call rt_task_list%list_t%init ('rt task list')
  call rt_boundary%init
  rt_boundary%ubox = mhd_patch%box                                    ! FIXME : add bndry patches
  rt_boundary%lbox = 0.0
  rewind (io%input); read (io%input, rt_params)
  if (io%master) write (*, rt_params)
  SELECT case(RT_solver)
  case('feautrier')
    solver = 0
  case('integral')
    solver = 1
  end SELECT
  end if
  !---------------------------------------------------------------------
  ! Loop through the task list, pick out the mhd tasks and create forge
  ! and other tasks around them
  !---------------------------------------------------------------------
  t_list => self%head
  do while (associated(t_list))
    !-------------------------------------------------------------------
    ! Get the right type of the task for the init call (mhd_t)
    !-------------------------------------------------------------------
    associate (mhd => t_list%task)
    SELECT TYPE(mhd)
    class is (mhd_t)
      mhd_patch => mhd
      allocate (mp, link)                                               ! the forge patch and its link
      mp%rt        = do_rt
      mp%chemistry = do_chemistry
      mp%mhd_patch => mhd_patch                                         ! pointer to mhd_patch
      mp%mhd_link  => t_list                                            ! pointer to mhd_patch's link
      link%task    => mp                                                ! connect link to forge task
      mp%link      => link
      mp%store_rt  = store_rt
      call mp%init                                                      ! initialize the forge patch
      call connect_nbors(t_list,link)                                   ! connect forge and mhd task
      call task_list%append_link(link)                                  ! insert rt_data into the task list
      !-----------------------------------------------------------------
      ! Set up whatever needed for the forge: RT, chemistry, anything
      !-----------------------------------------------------------------
      if (do_rt) call rt_setup(t_list, link, rt_task_list, solver,&
                               rt_freq, store_rt)                       ! initialize rt tasks for this mhd patch  
      !if (do_chemistry) call chemistry%whatever
    END SELECT
    end associate
    t_list => t_list%next                                               ! go to the next mhd task
  end do
  
  if (do_rt) then
    call add_upstream_nbors(rt_task_list)                               ! connect the rt internal tasks to appropriate upstreams
    call rt_task_list%print
    list => rt_task_list
    call task_list%append_list(list)                                    ! append rt task list to the forge task list
  end if
  list => task_list
  call self%append_list(list)                                           ! add forge task list to the main task list
  call trace_end
  if (io%master) then
    print *, '=================================='
    print *, 'MICROPHYSICS FACTORY COMPLETE'
    print *, '=================================='
  end if
END SUBROUTINE Microphysics_factory

END MODULE microphysics_factory_mod
