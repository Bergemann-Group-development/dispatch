MODULE restart_mod
  USE io_mod
  USE os_mod
  USE omp_mod
  USE trace_mod
  USE timer_mod
  USE mpi_mod
  USE data_io_mod
  USE kinds_mod
  USE task_factory_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE experiment_mod
  USE link_mod
  USE halo_mod
  USE list_mod
  USE task_list_mod
  USE read_bifrost_mod
  USE read_stagger_mod
  USE load_balance_mod
  USE bits_mod
  USE shared_mod
  USE divb_cleaner_mod
  implicit none
  private
  type, public:: restart_t
  contains
    procedure:: init
    procedure:: read
  end type
  character(120) :: data        = 'data/'
  character(120) :: run         = ''
  character(16)  :: method      = 'amr'
  integer        :: snapshot    = -1
  integer        :: verbose     = 0
  integer        :: nml_version = 4
  type(restart_t), public:: restart
CONTAINS

!===============================================================================
!> Select I/O method
!===============================================================================
SUBROUTINE init (self, task_list)
  class(restart_t)             :: self
  class(task_list_t),  pointer :: task_list
  integer                      :: iostat
  character(120)               :: filename
  namelist /restart_params/ snapshot, run, data, method, verbose
  !.............................................................................
  call trace%begin ('restart_t%init')
  !-----------------------------------------------------------------------------
  ! Namelist
  !-----------------------------------------------------------------------------
  run = io%run_name
  data = io%data_dir
  rewind (io%input)
  read (io%input, restart_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, restart_params)
  end if
  call io%namelist_warning ('restart_params', iostat)
  call io%ensure_dirname (data)
  write (io%namelist, restart_params)
  !-----------------------------------------------------------------------------
  ! Only allocate a task list if actually restarting (snapshot non-negative)
  !-----------------------------------------------------------------------------
  io%restart = snapshot
  io%inputdir = trim(data)//trim(run)
  call io%ensure_dirname (io%inputdir)
  if (snapshot >= 0) then
    !---------------------------------------------------------------------------
    ! Select method
    !---------------------------------------------------------------------------
    if (.not.associated(task_list)) then
      allocate (task_list)
      call task_list%init ('restart')
    end if
    select case (method)
    case ('bifrost')
      call read_bifrost%init (task_list)
    case ('stagger')
      call read_stagger%init (task_list)
    case default
      call self%read (task_list)
    end select
    !---------------------------------------------------------------------------
    ! Prepare task list for execution
    !---------------------------------------------------------------------------
    !call task_list%init_all_nbors ()
    call task_list%init_bdries
    call task_list%print_tasks
    call task_list%check_hash_table
  end if
  !.............................................................................
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Generic DISPATCH formats
!===============================================================================
SUBROUTINE read (self, task_list)
  class(restart_t)             :: self
  class(task_list_t),  pointer :: task_list
  !.............................................................................
  class(experiment_t), pointer :: experiment
  class(patch_t),      pointer :: patch
  class(task_t),       pointer :: task
  class(list_t),       pointer :: list
  class(link_t),       pointer :: link, next
  class(nbor_t),       pointer :: nbor
  integer                      :: n, unit, level, levelmin, levelmax, rank, &
                                  it, iout, i
  logical                      :: ok
  integer,                save :: itimer(3)=0
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read(1)', itimer=itimer(1))
              call io%header ('-restart_t%read begin', left=25)
              call timer%tic ()
  !-----------------------------------------------------------------------------
  ! Populate a temporart task list, so far with only light-weight tasks, based
  ! on metadata, and establish the level range
  !-----------------------------------------------------------------------------
  n = 0
  unit = 0
  levelmin = 9999
  levelmax = 0
  !-----------------------------------------------------------------------------
  ! Read in and count all tasks, and determine the largest ID
  !-----------------------------------------------------------------------------
  do while (.true.)
    allocate (experiment)
    patch => experiment
    call data_io%input_nml (patch, unit, verbose)
    call patch%set (bits%internal)
    iout = patch%iout
    if (data_io%iostat == unit) then
      call experiment%dealloc
      exit
    end if
    allocate (halo_t:: link)
    call link%init ()
    link%task => experiment
    experiment%link => link
    call task_list%append (link)
    n = n+1
    if (mod(n,1000)==0) &
      write (stdout,*) n, 'tasks'
  end do
  close (unit)
  task_list%next_free = task_list%max_id+1
  write (stdout,*) n, 'tasks, next_free =', task_list%next_free
  !-----------------------------------------------------------------------------
  ! Loop through list and pre_init() tasks, which may prepend more tasks
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    next => link%next
    level = link%task%level
    task => link%task
    select type (task)
    class is (experiment_t)
      task%task_list => task_list
      call task%pre_init ()
      task%level = level
      task%key = task%hash_key()
    class default
      call task%pre_init ()
      task%level = level
      task%key = task%hash_key()
    end select
    if (verbose > 1) &
      write (stdout,'(a,i6,":",i2.2,3i4,3f10.4,i4)') &
        'restart_t%read: id,key,pos,rank=', &
        link%task%id, link%task%key, link%task%position, link%task%rank
    call task_list%hash_table%add (link%task%key, link)
    io%start_time = max(io%start_time, experiment%time)
    levelmin = min(link%task%level,levelmin)
    levelmax = max(link%task%level,levelmax)
    link => next
  end do
  call os%mkdir_next_output (iout+1)
  levelmin = max(1,levelmin)
  shared%levelmin = levelmin
  shared%levelmax = levelmax
  if (verbose > 0) &
    write (stdout,*) 'task_list%n,levelmin,levelmax:', &
      task_list%n, levelmin, levelmax
  if (task_list%n > 10000) verbose = verbose+1
            if (verbose > 0) then
              call timer%toc ()
              call io%header ('-sort task list', left=25)
              call timer%tic ()
            end if
  call task_list%sort_by_level ()
  call task_list%info ("after sort")
  !.............................................................................
  if (verbose > 1) then
    link => task_list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (patch_t)
        write (stdlog,'(a,i6,2i3,f10.3)') &
          'restart_t: id,level,idx%d,z=', &
          task%id, task%level, task%idx%d, task%position(3)
      end select
      link => link%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Decide which tasks to keep, and append them to a new task list
  !-----------------------------------------------------------------------------
            if (verbose > 0) then
              call timer%toc ()
              call io%header ('-creating nbor lists', left=25)
              call timer%tic ()
            end if
  call trace%end (itimer(1))
  call trace%begin ('restart_t%read(2)', itimer=itimer(3))
  call task_list%init_all_nbors ()
  call init_all_parents (task_list)
  call mpi%barrier ('1')
  list => task_list
            if (verbose > 0) then
              call timer%toc ()
              call io%header ('-initial load balancing', left=25)
              call timer%tic ()
            end if
  call load_balance%init (list)
  call load_balance%initial (list)
            if (verbose > 0) then
              call timer%toc ()
              call io%header ('-creating local task list', left=25)
              call timer%tic ()
            end if
  !call task_list%init_all_nbors ()
  !call init_all_parents (task_list)
  call task_list%set_io_counts()
  call mpi%barrier ('restart_mod%read: task%init')
  !-----------------------------------------------------------------------------
  ! Allocate task memory and read in data for remaining tasks
  !-----------------------------------------------------------------------------
            if (verbose > 0) then
              call timer%toc ()
              call io%header ('-reading in patch data', left=25)
              call timer%tic ()
            end if
  !$omp parallel default(shared)
  call init_parallel (task_list)
  !$omp end parallel
  call mpi%barrier ('restart_mod%read')
  !-----------------------------------------------------------------------------
  ! Clean div(B) from patch interiors (if divb_cleaner_params on=t)
  !-----------------------------------------------------------------------------
  call divb_cleaner%init 
  !$omp parallel default(shared)
  call divb_cleaner%update (task_list)
  !$omp end parallel
  call mpi%barrier ('restart_mod%read')
  !.............................................................................
  if (io%in%do_generic) then
    write (stdout,*) task_list%n, ' tasks read with generic I/O'
  else
    write (stdout,*) task_list%n, ' tasks read with solver-specific I/O'
  end if
  !.............................................................................
              call timer%toc ()
              call io%header ('-restart_t%read end', left=25)
  call trace%end (itimer(3))
END SUBROUTINE read

!===============================================================================
!> Initialize tasks in parallel, to ensure memory affinity spread
!===============================================================================
SUBROUTINE init_parallel (task_list)
  class(task_list_t),  pointer :: task_list
  class(link_t), pointer:: link  
  class(task_t), pointer:: task  
  integer:: i, n, it, l(3), u(3)
  logical:: ok
  !-----------------------------------------------------------------------------
  i = 0
  n = 0
  link => task_list%head
  do while (associated(link))
    if (mod(i,omp%nthreads)==omp%thread) then           ! thread selected
      !$omp critical (in_cr)
      n = n+1
      task => link%task
      !PRINT *,omp%thread, task%id
      select type (task)
      class is (experiment_t)
        call task%init ()
        task%t = task%time
        task%dt = 0d0
        call task%input (ok)
        !-----------------------------------------------------------------------
        ! Copy to other mem slots
        !-----------------------------------------------------------------------
        do it=1,task%nt
          if (it /= task%it) then
            task%mem(:,:,:,:,it,:) = task%mem(:,:,:,:,task%it,:)
          end if
        end do
        if (verbose > 1) then
          l = task%li
          u = task%ui
          write (stdlog,'(a,i6,1p,2g12.3)')  &
            'init_parallel: id,min/maxval(d):', task%id,&
            minval(task%mem(l(1):u(1),l(2):u(2),l(3):u(3),task%idx%d,task%it,1)), &
            maxval(task%mem(l(1):u(1),l(2):u(2),l(3):u(3),task%idx%d,task%it,1))
          flush (stdlog)
        end if
      class is (gpatch_t)
        call task%init ()
        task%t = task%time
      end select
      !$omp end critical (in_cr)
    end if
    i = i+1
    link => link%next
  end do
END SUBROUTINE init_parallel

!===============================================================================
!> Add nbor lists to all task in the task_list, based only on overlap
!===============================================================================
SUBROUTINE init_all_nbors (task_list)
  class(list_t)          :: task_list
  !.............................................................................
  class(link_t), pointer :: link, scan
  class(nbor_t), pointer :: nbor
  integer                :: level, levelmin, levelmax, n
  integer                :: i, i0, i1, i2
  real(8)                :: limit(3), distance(3)
  logical                :: rejected
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%init_all_nbors', itimer=itimer)
  link => task_list%head
  if (.not.associated(link)) then
    call io%abort ('no task_list in list_t%init_all_nbors')
  end if
  levelmin = link%task%level
  levelmax = link%task%level
  do while (associated(link))
    levelmin = min(link%task%level,levelmin)
    levelmax = max(link%task%level,levelmax)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  i2 = task_list%n/omp%nthreads + 1                       ! tasks per thread
  !!omp parallel private(i,i0,i1,link,scan,distance) &
  !!omp & shared(task_list,i2,levelmin,levelmax,nbor,limit) default(none)
  i0 = omp%thread*i2
  i1 = i0 + i2
  link => task_list%head                                  ! start checking links
  i = 0                                                   ! task counter
  do while (associated(link))
    do level=max(link%task%level-2,levelmin),min(link%task%level+2,levelmax)
      if (i >= i0 .and. i<i1) then                        ! thread interval
      if (associated(task_list%first(level)%pointer)) then
        scan => task_list%first(level)%pointer
        do while (associated(scan))
          if (scan%task%level == level) then
            if (scan%task%id /= link%task%id) then
              limit = 0.501d0*(link%task%size + scan%task%size)
              distance = link%task%distance_to(scan%task) 
              if (all(abs(distance) < limit)) then
                call link%add_nbor (scan, nbor=nbor, rejected=rejected)
                if (rejected) then
                  call task_list%check_duplicates (scan)
                end if
                if (associated(nbor)) then
                  nbor%task => scan%task
                end if
              end if
            end if
          else
            exit
          end if
          scan => scan%next
        end do
        end if
      end if
    end do
    i = i + 1
    link => link%next
  end do
  !!omp end parallel
  call task_list%reset_status ()
  call trace%end (itimer)
END SUBROUTINE init_all_nbors

!===============================================================================
!> For multi-level snapshots, find the parent tasks and set pointers
!===============================================================================
SUBROUTINE init_all_parents (task_list)
  class(list_t)          :: task_list
  class(link_t), pointer :: link, parent
  class(task_t), pointer :: task
  integer                :: key(4), id_max, id_master
  !.............................................................................
  call trace%begin ('restart_t%init_all_nbors')
  id_max = 0
  link => task_list%head
  do while (associated(link))
    task => link%task
    id_max = max(task%id,id_max)
    if (task%level > shared%levelmin .and. &
        task%is_clear (bits%virtual)) then
      key = task%hash_key (level=task%level-1, &
        position=task%position, size=task%size*shared%refine_ratio)             ! call to hash_key overloadable ove rall levels
      call task_list%hash_table%get_link (key, parent)
      if (associated(parent)) then
        link%parent => parent
        select type (task)
        class is (patch_t)
        task%parent => parent%task
        end select
      else
        write (stdout,*) &
          'WARNING: orphan task:', task%id, mpi%rank
      end if
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  id_master = mpi%id%update (0)
  if (id_max > id_master) then
    id_master = mpi%id%update (id_max-id_master)
  end if
  call trace%end ()
END SUBROUTINE init_all_parents

END MODULE restart_mod
