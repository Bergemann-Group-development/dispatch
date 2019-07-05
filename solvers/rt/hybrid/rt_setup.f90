!===============================================================================
!> Add RT subtasks into the general task_list.  Each subtask is added to the
!> nbor list of the parent task, and the parent task should be added to the nbor
!> list of the subtask.
!>
!> Reminder about links and nbors:   A task should be represented on the task_list
!> by a link, which has %next and %prev pointers that are not used for any other
!> purpose.   The link must contain an %nbor pointer, which is the head of a chain
!> of links whose tasks are to provide guard zones.  The main and nbor links have
!> the relations
!>
!> link%nbor%task => nbtask
!> nbtask%link%task => nbtask
!>
!> The nbor lists have special behavior, in that even though the parent RT task
!> needs to check that the sub-tasks are up-to-date, it does not need to, and
!> should not, use a standard download_link call visavis the RT sub-tasks
!===============================================================================
SUBROUTINE rt_setup (task_list)
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  USE task_list_mod
  USE list_mod
  USE link_mod
  USE task_mod
  USE solver_mod
  USE patch_mod
  USE trace_mod
  USE rt_mod
  USE bits_mod
  USE data_io_mod
  USE download_mod
  USE ray_geometry_v3_mod
  USE mesh_mod
  implicit none
  type(task_list_t):: task_list
  !.............................................................................
  class(link_t), pointer:: link, nbor, new
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  class(mesh_t),  pointer:: mesh
  class(rt_t), pointer:: omega=>null(), rt_patch
  class(solver_t), pointer:: mhd
  class(task_list_t), pointer:: rt_task_list
  integer:: i_omega, ntask, iter, iostat, id_save, n_mhd, n_rt, id, ng(3)
  integer, save:: verbose=0, n_warmup=5, n_omega=0
  real(8):: boundary_distance(3)
  logical:: on = .false.
  logical, save:: unique_ids=.false.
  logical, save:: only_vertical=.false.
  logical, save:: omp_warmup = .false.
  namelist /rt_setup_params/ n_warmup, verbose, unique_ids, only_vertical, omp_warmup
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup')
  call io%header ('begin rt_setup')
  !-----------------------------------------------------------------------------
  rewind (io_unit%input)
  read (io_unit%input, rt_setup_params, iostat=iostat)
  if (mpi%size > 1 .and. .not.unique_ids) then
    write(io%output,*) &
      'WARNING: unique_ids set to .true. when using MPI size =', mpi%size
    unique_ids = .true.
  end if
  if (io%master) write (io_unit%output, rt_setup_params)
  !-----------------------------------------------------------------------------
  allocate (rt_task_list)
  call task_list%initialize
  call rt_task_list%init ('RT_list')
  if (verbose > 1) &
    call task_list%print ('task_list')
  n_rt = 0
  n_mhd = 0
  link => task_list%head
  task => link%task
  SELECT TYPE(task)
  class is (solver_t)
    mhd => task
  class default
    error stop '(1) in rt_setup: unknown task type'
  end SELECT
  !> -------------------------------------------------------------------
  !> Init the ray_geometry objet.
  !> -------------------------------------------------------------------
  ng = 1
  mesh => mhd%mesh(1)
  select type(mesh)
    type is(Cartesian_mesh)
      call ray_geometry%init(geometry=1, &
      n=mhd%n, ng=ng, cell_dims=mhd%ds, no_mans_land=mhd%no_mans_land)
    type is(spherical_mesh)
      call ray_geometry%init(geometry=2, &
      n=mhd%n,  ng=ng, cell_dims=mhd%ds, no_mans_land=mhd%no_mans_land)
    type is(cylindrical_mesh)
      call ray_geometry%init(geometry=3, &
      n=mhd%n,  ng=ng, cell_dims=mhd%ds, no_mans_land=mhd%no_mans_land)
    class default
      error stop '(2) in rt_setup: unknown mesh type'
  end select
  print *, 'ray geometry initialized'
  
  ! ----------------------------------------------------------------------------
  ! Init rt tasks
  ! ----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    call io%assert (associated(link%task), 'link%task missing')
    task => link%task
    select type (task)
    class is (solver_t)
      mhd => task
      call io%assert (associated(mhd%link%task), 'mhd%link%task missing')
      
      call mhd%rt%init
      link => link%next
    end select
  end do
  !
  link => task_list%head
  do while (associated(link))
    call io%assert (associated(link%task), 'link%task missing')
    task => link%task
    select type (task)
    class is (solver_t)
      mhd => task
      !-------------------------------------------------------------------------
      ! Add nbor lists to the sub-tasks, if RT is active in the MHD patch
      !-------------------------------------------------------------------------
      if (mhd%is_clear(bits%no_rt)) then
        call mhd%rt%init_nbors
        !-----------------------------------------------------------------------
        ! Append the main task links to the task list
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        ! Add a link to the rt task to the MHD task nbor list
        !-----------------------------------------------------------------------
        omega => mhd%rt
        if (omega%on) then
          on = .true.
          allocate (new)
          new%link => omega%link
          new%task => omega
          !call link%add_nbor_by_rank(link%nbors_by_level, new)
          call link%add_nbor_by_rank(link%nbor, new)
          !-----------------------------------------------------------------------
          ! Add a link to the MHD task to the rt task nbor list
          !-----------------------------------------------------------------------
          allocate (new)
          new%link => link
          new%task => task
          !call omega%link%add_nbor_by_rank(omega%link%nbors_by_level,new)
          call omega%link%add_nbor_by_rank(omega%link%nbor,new)
          !-----------------------------------------------------------------------
          ! Add the rt task to the RT task list
          !-----------------------------------------------------------------------
          call omega%clear (bits%ready)
          call rt_task_list%append_link (omega%link)
          n_rt = n_rt+1
        end if
      else
        n_mhd = n_mhd+1
      end if
    end select
    if (verbose > 1) &
      call link%print_nbors
    link => link%next
  end do
  write(io%output,*) 'rt_setup: n_rt, n_mhd =', n_rt, n_mhd
  if (associated(omega)) then
    patch => omega
    call data_io%register_kind (patch)
  end if
  write(io%output,'(1x,a,i5)') 'rt_setup: io%nwrite =', io%nwrite
  !-----------------------------------------------------------------------------
  ! EOS for warm up
  !-----------------------------------------------------------------------------
  write(io%output,*) 'EOS for RT'
  !!omp parallel default (shared) private(link, task, mhd)
  !!!mp master
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (solver_t)
      mhd => task
      if (task%is_clear(bits%no_rt)) then
        !-----------------------------------------------------------------------
        ! Spawn OMP tasks to compute EOS
        !-----------------------------------------------------------------------
        if (omp%nthreads > 1) then
          !!omp task default(shared) firstprivate(mhd)
            mhd%mem(:,:,:,:,mhd%new,1) = mhd%mem(:,:,:,:,mhd%it,1)
            call mhd%rt%rt_eos (mhd%time)
          !!omp end task
        else
            mhd%mem(:,:,:,:,mhd%new,1) = mhd%mem(:,:,:,:,mhd%it,1)
          call mhd%rt%rt_eos (mhd%time)
        end if
      end if
    end select
    link => link%next
  end do
  !!omp end master
  !!omp end parallel
  !-----------------------------------------------------------------------------
  ! Merge task lists, if RT is present
  !-----------------------------------------------------------------------------
  if (verbose > 1) &
    call rt_task_list%print ('rt_task_list')
  call task_list%append_task_list (rt_task_list)
  link => task_list%head
  ! Sort out the nbors list...
  do while (associated(link))
    call link%sort_nbors_by_level(link%nbors_by_level,link%nbors_by_level)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Reset task_list status to establish nbor relations
  !-----------------------------------------------------------------------------
  call task_list%reset_status()
  !
  if (on) then
    !---------------------------------------------------------------------------
    ! Warm up the RT
    !---------------------------------------------------------------------------
    do iter=1,n_warmup
      write(io%output,*) 'warm-up, iteration', iter
      if (omp_warmup) then
        !$omp parallel default(shared)
        link => task_list%head
        id =0
        do while (associated(link))
          id = id+1
          if (mod(id,omp%nthreads)==omp%thread) then
            task => link%task
            select type (task)
            class is (rt_t)
            if (task%on) then
              call task%dnload
              call task%warm_up()
            end if
            end select
          end if
          link => link%next
        end do
        !$omp end parallel
        if (omp_warmup .and. omp%nthreads>1) then
          task%mem_thread = omp%thread
        else
          task%mem_thread = -1
        end if
      else
        link => task_list%head
        do while (associated(link))
          task => link%task
          select type (task)
          class is (rt_t)
          if (task%on) then
            call task%dnload
            call task%warm_up()
          end if
          end select
          link => link%next
        end do
      end if
      
    end do
    !---------------------------------------------------------------------------
    ! Collect the first heating rates
    !---------------------------------------------------------------------------
    write(io%output,*) 'collecting the heat'
    link => task_list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (rt_t)
      if (task%on) then
        call task%warm_up
      end if
      end select
      link => link%next
    end do
  end if
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    call task_list%print ('new task_list')
  call trace%end
END SUBROUTINE
