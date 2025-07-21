!===============================================================================
!> Reads a tabulated model atmosphere from a text file
!===============================================================================
MODULE read_external_mod
  USE io_mod
  USE os_mod
  USE mpi_coords_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE units_mod
  USE scaling_mod
  USE trace_mod
  USE task_mod
  USE link_mod
  USE mesh_mod
  USE bits_mod
  USE halo_mod
  USE list_mod
  USE shared_mod
  USE task_list_mod
  USE experiment_mod, only: experiment_t
  USE composition_mod
  USE composition_patch_mod, only: FeH, save_atmos, save_eos, external_abund, read_composition_params
  USE patch_mod
  USE const_mod, only:cc
  USE m3d_patch_mod, only: maxiter, ilambd, m3d_patch_t, atom_update_counter, &
  global_maxerr, global_pos, nz_total, use_rho, use_ne, atmos_format, amr, vmic, &
  opac_table_vmic, atmos_file, use_vel
  USE m3d_halo_mod, only: m3d_halo_t
  USE read_atmos_mod, only: ia, read_stagger, read_co5bold, read_stagger2, &
  legacy_stagger_interp, read_marcs, read_multi_snapshot, read_text_1d, read_multi_1d, &
  read_siu_bin, timestamp
  USE handyfuncs_mod
  USE ray_mod, only: one_dimensional
  implicit none

  type, public:: read_external_t
    class(task_list_t), pointer:: task_list
    integer :: dims(3)=[1,1,1], n(3), verbose
    real(8) :: ds(3), size(3)
  contains
    procedure:: init
    procedure:: init_exp
    procedure:: setup_patch
    procedure:: read_saved_snapshot
  end type read_external_t

  type, public, extends(task_t):: rank_t
  end type

  type(read_external_t), public:: read_external
  type(rank_t) :: rank
  real(4) :: origin(3) = [0., 0., 0.]

contains

!===============================================================================
SUBROUTINE init(self, task_list, name)
  !-----------------------------------------------------------------------------
  class(read_external_t):: self
  class(task_list_t),  pointer:: task_list
  character(len=*),   optional:: name
  !-----------------------------------------------------------------------------
  class(experiment_t), pointer:: patch
  class(link_t), pointer:: link, next
  integer :: lunit, iostat, k, ip
  integer :: patches_per_mpi=1
  integer :: dims=1, nx=-1, ny=-1, nz=128, snap=0
  integer :: verbose=0
  logical :: use_abnd=.false.
  character(len=512) :: command
  character(len=512) :: save_file
  namelist /atmos_params/ atmos_file, dims, verbose, patches_per_mpi, FeH, &
  save_atmos, atmos_format, vmic, use_rho, use_ne, nx, ny, nz, snap, &
  legacy_stagger_interp, save_eos, amr, use_abnd, timestamp, use_vel
  !-----------------------------------------------------------------------------
  if (associated(task_list)) return
  call trace%begin('read_external_t%init')
  call io%header ('read_external_t%init')
  ! default values
  self%dims = 1
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and saves them
  !-----------------------------------------------------------------------------
  call io%lock%set
  rewind (io%input)
  read (io%input, nml=atmos_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, atmos_params)
  end if
  call io%namelist_warning ('atmos_params', iostat)

  if (amr) then
    ! making sure that dims is not 1, otherwise amr has no effect
    dims = max(dims, 8)
  endif

  if (mpi%master) then
    call os%mkdir(trim(io%run_dir)//'save')
    write (command, *) 'rm -f ', trim(io%run_dir), '*.txt ', trim(io%run_dir), '*.bin'
    call execute_command_line (trim(command), exitstat=iostat)
    write (command, *) 'rm -f ', trim(io%run_dir), 'save/*.txt ', trim(io%run_dir), 'save/*.bin'
    call execute_command_line (trim(command), exitstat=iostat)
    call os%mkdir(trim(io%run_dir)//'atoms')
    write (command, *) 'rm -f ', trim(io%run_dir), 'atoms/*'
    call execute_command_line (trim(command), exitstat=iostat)
  endif

  if (patches_per_mpi > 1) then
    if (dims > 1) then
      print *, 'Cannot give PATCHES_PER_MPI if DIMS is also specified in ATMOS_PARAMS. Stopping...'
      stop
    else
      self%dims(3) = mpi%size * patches_per_mpi
    endif
  else
    if (dims > 1) then
      self%dims(3) = dims
    else
      self%dims(3) = mpi%size
    endif
  endif

  if (mpi%size > self%dims(3)) then
    print *, 'More mpi ranks requested than spatial atmosphere splits. Stopping...'
    stop
  endif

  self%verbose = verbose
  call io%lock%unset
  
  atmos_format = lowercase(atmos_format)
  ia%kind = atmos_format

  if (vmic > 0) opac_table_vmic = vmic

  if (atmos_format == 'save') then
    save_file = atmos_file
    open(file=trim(save_file)//'/params.nml', newunit=lunit, form='FORMATTED')
    read (lunit, nml=atmos_params, iostat=iostat)
    close(lunit)
    atmos_file = save_file
    atmos_format = 'save'
    call self%read_saved_snapshot(atmos_file, dims)
  else
    select case (atmos_format)
    case('stagger')
      call read_stagger(atmos_file, snap, use_ne)
    case('stagger2')
      call read_stagger2(atmos_file)
    case('co5bold')
      call read_co5bold(atmos_file)
    case ('marcs')
      call read_marcs(atmos_file, FeH, vmic, external_abund)
      if (.not. use_abnd) deallocate(external_abund)
    case ('multi', 'must')
      call read_multi_snapshot(atmos_file)
    case ('text')
      call read_text_1d(atmos_file, vmic, FeH)
    case ('multi1d')
      call read_multi_1d(atmos_file, vmic)
    case ('siu_bin')
      call read_siu_bin(atmos_file, vmic)
    case default
      write(stdout, *) 'unknown atmos_format ', trim(atmos_format), &
      ". Options are: 'Stagger', 'Co5bold', 'MUST', 'Multi', &
      &'Marcs', 'SIU-bin', 'Multi1D' and 'Text'."
      stop
    end select
  endif

  ! this is read here, in case atmosphere has optical depth scale
  call read_composition_params()

  if (atmos_format /= 'save') then
    call ia%thick2thin()
    call ia%get_box(self%size)
    self%size = self%size / scaling%l

    if (.not. ia%is_3d) then
      nx = 1; ny = 1
    elseif (nx < 0 .or. ny < 0) then
      nx = ia%nx; ny = ia%ny
    endif
    self%n = [nx, ny, nz]

    if (amr) then
      call ia%setup_amr_mesh(dims, use_ne)
    else
      self%ds = self%size / self%n
      if (save_atmos) then
        open(file=trim(io%run_dir)//'atmos_mesh.txt', newunit=lunit, form='FORMATTED')
        write(lunit, *) self%n
        write(lunit, *) self%ds
        close(lunit)
      endif
    endif
  endif

  one_dimensional = self%n(1)*self%n(2) == 1

  write (io%namelist, atmos_params)

  !-----------------------------------------------------------------------------
  ! Generate a (periodic) MPI geometry, with dims specified in the input file,
  ! or assigned automatically by mpi_coords%init.
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) then
     mpi_coords%dims = 1
  else
    call mpi_coords%init (mpi_dims=[1,1,mpi%size])
  end if
  patches_per_mpi = self%dims(3)/mpi%size
  mpi_coords%npatch(3) = patches_per_mpi
  io%mpi_dims = mpi_coords%dims
  !-----------------------------------------------------------------------------
  ! Set attritbutes of the rank_t data type -- just enough to use distance().
  !-----------------------------------------------------------------------------
  rank%periodic = .true.
  rank%box = self%size
  rank%size = rank%box/mpi_coords%dims
  rank%ds(1:2) = rank%size(1:2)
  rank%ds(3) = rank%size(3)/patches_per_mpi
  rank%position = (mpi_coords%rank_to_coords(mpi%rank) + 0.5_8)*rank%size + origin
  !-----------------------------------------------------------------------------
  ! Begin on task list
  !-----------------------------------------------------------------------------
  allocate (task_list)
  self%task_list => task_list
  call self%task_list%init ('read_external_mod')
  self%task_list%size = self%size
  self%task_list%dims = self%dims
  self%task_list%position = 0.5*self%size + origin
  self%task_list%n_tasks = product(self%dims)

  !-----------------------------------------------------------------------------
  ! Allocate an array of links, pointing to patches, and give the patches each
  ! id, rank, size, position, integer position, and status bits.  Discard
  ! external patches and add the rest to the task_list.
  !-----------------------------------------------------------------------------

  ip = 0
  do k=0, self%dims(3)-1
    allocate (patch)
    patch%kind = 'mhd_patch'
    ! --- patch 3D and 1D coordinates ---
    patch%ipos = [0,0,k]
    patch%rank = mpi_coords%coords_to_rank ([0,0,k]/[1,1,patches_per_mpi])
    if (patch%rank == mpi%rank) then
      ip = ip + 1
      patch%ip = ip
    end if
    if (all(self%dims==1)) call patch%task_t%set(bits%root_grid)
    call self%setup_patch(patch)

    if (amr .and. save_atmos) then
      open(file=trim(io%run_dir)//'atmos_mesh.txt', newunit=lunit, form='FORMATTED', access='append')
      write(lunit, *) patch%n
      write(lunit, *) patch%ds
      close(lunit)
    endif

    if (verbose > 0) then
      if (k==0) write(stdout,'(I2, A28)') self%dims(3), ' subdomains in z-direction:'
      write(stdout,'(I4, I6,A3,I4,4X,3I4, 2X, F7.4, E14.4)') k+1, patch%zs+1, ' - ', patch%ze, &
      patch%n, (patch%position(3)-patch%origin(3))/patch%box(3), patch%ds(3)
    endif
  end do
  
  nz_total = patch%ze
  
  if (verbose > 0) then
    write (stdout, '(A, 3(F10.4))') 'box = ', patch%box
    write (stdout,*) ip, 'tasks generated'
  endif

  io%end_time = (maxiter + ilambd + 1) * patch%dtime
  io%print_every = huge(1)
  timer%dead_mans_hand = huge(1.0)

  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_levels'
  call task_list%init_levels()
  shared%levelroot = shared%levelmin
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_ids'
  call init_ids(task_list)
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_all_nbors'
  call init_all_nbors(task_list)
  if (verbose > 0) write(stdout,*) 'read_external_t%init: remove_external'
  call remove_external(task_list)
  !if (verbose > 1) then
  !  write(stdout,*) 'read_external_t%init: info'
  !  call task_list%info
  !endif
  !-----------------------------------------------------------------------------
  ! This counts are not fully correct until the patch meshes have been allocated
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(stdout,*) 'read_external_t%init: count status'
  call self%task_list%count_status
  !-----------------------------------------------------------------------------
  ! The I/O conunters must be set before initializing patches
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  if (verbose > 0) write (stdout,'(a,i8)') ' number of tasks in task list:', self%task_list%n
  io%nwrite = product(self%dims)
  write (stdlog,*) 'read_external_mod: ntask,nwrite=', io%ntask, io%nwrite

  !-----------------------------------------------------------------------------
  ! Initialize mpi counters
  !-----------------------------------------------------------------------------
  call atom_update_counter%init(0)
  call global_maxerr%init(0)
  call global_pos%init(-1)

  !-----------------------------------------------------------------------------
  ! Initialize tasks
  !-----------------------------------------------------------------------------
  !$omp parallel default(shared)
  call self%init_exp()
  !$omp end parallel
  
  call print_list(task_list)
  !-----------------------------------------------------------------------------
  ! Refresh the counts, now with bits%has_output set on tasks
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdlog,*) 'read_external_mod: ntask,nwrite=', io%ntask, io%nwrite
  if (verbose > 0) write(stdout,*) 'read_external_t%init: print tasks'
  call self%task_list%print_tasks
  !call self%task_list%check_hash_table
  !call self%task_list%hash_table%stats
  if (verbose > 0) write(stdout,*) 'read_external_t%init: done'
  flush (stdout)
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE read_saved_snapshot (self, atmos_file, dims)
  class(read_external_t):: self
  character(len=*), intent(in):: atmos_file
  integer, intent(in):: dims
  integer:: lunit, iostat, i, nz
  !-----------------------------------------------------------------------------
  save_atmos = .false.
  self%dims(3) = dims
  allocate(ia%sub_size(dims))
  allocate(ia%sub_pos(dims))

  call os%check_file(trim(atmos_file)//'/atmos_mesh.txt')
  open(file=trim(atmos_file)//'/atmos_mesh.txt', newunit=lunit, form='FORMATTED', iostat=iostat)
  nz = 0
  self%size = 0
  do i=1, dims
    read(lunit, *) self%n
    read(lunit, *) self%ds
    nz = nz + self%n(3)
    ia%sub_size(i) = self%n(3) * self%ds(3)
    ia%sub_pos(i) = 0.5*ia%sub_size(i) + self%size(3)
    self%size(3) = self%size(3) + self%n(3)*self%ds(3)
  end do
  self%n(3) = nz
  self%size(1:2) = self%n(1:2)*self%ds(1:2)
  close(lunit)
END SUBROUTINE read_saved_snapshot

!===============================================================================
!> Initialize a minimum amount of patch information in a Cartesian arrangement
!===============================================================================
SUBROUTINE setup_patch (self, patch)
  class (read_external_t):: self
  class (experiment_t), intent(inout), target :: patch
  integer :: md(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%setup_patch')
  patch%task_list => self%task_list
  call patch%patch_t%pre_init() ! this overwrites patch%n
  patch%pre_initialized = .true.

  patch%box    = self%size
  patch%nbox   = self%n
  patch%dims   = self%dims
  patch%origin = origin
  
  call patch%pre_init()
  call trace%end ()
END SUBROUTINE setup_patch

!===============================================================================
SUBROUTINE init_exp (self)
  class(read_external_t):: self
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  logical:: periodic(3)
  integer:: id, thread
  !-----------------------------------------------------------------------------
  id = 0
  link => self%task_list%head
  do while (associated(link))
    id = id+1
    thread = mod(id,omp%nthreads)
    if (thread==omp%thread) then
      task => link%task
      select type (task)
      class is (experiment_t)
        call task%init()
      end select
    endif
    link => link%next
  end do
END SUBROUTINE init_exp

!===============================================================================
SUBROUTINE init_all_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: i, i0, i1, i2, m, n=0, nv=0, nb=0
  real(8):: time
  integer, save:: itimer=0
  !.............................................................................
  if (self%hands_off_nbor_lists) return
  call trace%begin('read_external_mod%init_all_nbors', itimer=itimer)
  if (self%verbose > 0) &
    call timer%tic (time)
  i2 = self%n/omp%nthreads + 1                         ! tasks per thread
  n = 0
  !$omp parallel private(i,m,i0,i1,link) shared(self,i2,n,stdout) default(none) num_threads(omp%nthreads)
  i0 = omp%thread*i2
  i1 = i0 + i2
  if (omp%master .and. self%verbose > 0) &
    write (stdout,*) self%n, 'task to generate nbor lists for'
  link => self%head                                     ! start checking links
  i = 0                                                 ! task counter
  do while (associated(link))
    if (i >= i0 .and. i<i1) then                        ! thread interval
      call init_nbors (self, link)
    end if
    i = i+1                                             ! increment counter
    link => link%next                                   ! continue list
  end do
  !$omp end parallel
  link => self%head
  do while (associated(link))
    call link%task%sort_nbors()
    link => link%next
  end do
  if (self%verbose > 0) &
    call timer%toc ('time to initialize nbor lists', 1, time)
  call trace%end (itimer)
END SUBROUTINE init_all_nbors

!===============================================================================
SUBROUTINE init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(nbor_t), pointer:: new_nbors, old_nbors
  integer, save:: itimer=0, n_max=10000
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_mod%init_nbors', itimer=itimer)
  call link%task%log ('init_nbors')
  nullify (new_nbors)
  if (self%use_hash) then
    !$omp atomic write
    link%task%n_nbors = 0
    select type (task => link%task)
    class is (m3d_patch_t)
      call task%init_nbors_hash(new_nbors)
    class is (experiment_t)
      call task%init_nbors_hash(new_nbors)
    end select
  else
    call init_nbors_by_level (self, link, new_nbors)
  end if
  !-----------------------------------------------------------------------------
  ! Swap the link%nbor pointer under lock
  !-----------------------------------------------------------------------------
  call                                              link%lock%set ('init_nbors')
  old_nbors => link%nbor
  link%nbor => new_nbors
  !-----------------------------------------------------------------------------
  ! Update the status bits based on the new nbor list and adjust counts
  !-----------------------------------------------------------------------------
  call link%reset_status (link%nbor)
  link%task%previous_status = link%task%status
  !call self%update_status (link)
  call                                            link%lock%unset ('init_nbors')
  !-----------------------------------------------------------------------------
  ! Clear the init_nbors bit and remove the old nbor list
  !-----------------------------------------------------------------------------
  call link%task%clear (bits%init_nbors)
  call link%remove_nbor_list (old_nbors)
  !.............................................................................
  if (self%verbose > 0) &
    call link%info ('init_nbors')
  call trace%end (itimer)
END SUBROUTINE init_nbors

!===============================================================================
SUBROUTINE init_nbors_by_level (self, link, nbors)
  class(list_t)          :: self
  class(link_t)          :: link
  class(nbor_t), pointer :: nbors
  !.............................................................................
  class(link_t), pointer :: scan
  class(task_t), pointer :: task, nbtask
  class(nbor_t), pointer :: nbor, nbor2
  real(8)                :: distance(3)
  logical                :: added
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  ! Save a pointer to the current nbor list, to remove after adding a new
  ! (to ensure that task%n_needed does not touch zero prematurely)
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_mod%init_nbors_by_level', itimer=itimer)
  nullify (nbors)
  !-----------------------------------------------------------------------------
  ! Run through (only) tasks at the relevant levels
  !-----------------------------------------------------------------------------
  task => link%task
  scan => self%head
  do while (associated(scan))
    nbtask => scan%task
    if (task%rank == mpi%rank .or. nbtask%rank == mpi%rank) then
      !--- does it have a chance to overlap? ---
      distance = task%distance_to (nbtask)
      if (all(distance < huge(1.0))) then
        select type (task)
        class is (patch_t)
          if (task%interest%overlaps (nbtask%authority, distance)) then
            if (nbtask%id /= task%id) then
              call add_nbor (link, scan, nbors, nbor, added)
              if (added) nbor%distance = distance
            end if
          end if
        end select
      endif
    endif
    scan => scan%next
  end do
  call trace%end (itimer)
END SUBROUTINE init_nbors_by_level

!===============================================================================
SUBROUTINE add_nbor (self, link, nbors, nbor, added, by_quality, by_rank)
  class(link_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer, optional:: nbors
  class(nbor_t), pointer, optional:: nbor
  logical,                optional:: added
  logical,                optional:: by_quality, by_rank
  !.............................................................................
  logical:: ok
  class(nbor_t), pointer:: next, prev, nbor_l
  !integer, save:: itimer=0
  !.............................................................................
  !call trace%begin ('read_external_mod%add_nbor', itimer=itimer)
  call link%lock%set ('add_nbor')
  !-----------------------------------------------------------------------------
  ! Use a class method to overload the allocation
  !-----------------------------------------------------------------------------
  call self%alloc_nbor (link, nbor_l)
  !-----------------------------------------------------------------------------
  ! Consider this carefully!
  ! From Andrius: please consult me before changing this.
  !-----------------------------------------------------------------------------
  !call self%task%nbor_relations (nbor_l%task, nbor_l%needed, nbor_l%needs_me, &
  !                               nbor_l%download)
  call self%check_nbor_relations(nbor_l%task, nbor_l%needed, nbor_l%needs_me, & 
                                 nbor_l%download)
  if (.not. nbor_l%needed   .and. &
      .not. nbor_l%needs_me .and. &
      .not. nbor_l%download) then
    deallocate(nbor_l)
    call link%lock%unset('add_nbor')
    !call trace%end(itimer)
    if (present(added)) added=.false.
    return
  end if
  ok = .true.
  nullify (prev)
  if (present(nbors)) then
    next => nbors
  else
    next => self%nbor
  end if
  do while (associated(next))
    if (next%task%id == link%task%id) then
      !-------------------------------------------------------------------------
      ! Return if link task is already in the nbor list
      !-------------------------------------------------------------------------
      ok = .false.
      exit
    else
      if (next%task%level > link%task%level) then
        !-----------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-----------------------------------------------------------------------
        nbor_l%next => next
        exit
      end if
    end if
    prev => next
    next => next%next
  end do
  if (ok) then
    !---------------------------------------------------------------------------
    ! Found a place to insert, it link it in
    !---------------------------------------------------------------------------
    if (associated(prev)) then
      prev%next => nbor_l
    else
      !-------------------------------------------------------------------------
      ! If prev is not associated, this is the first nbor
      !-------------------------------------------------------------------------
      if (present(nbors)) then
        nbors => nbor_l
      else
        self%nbor => nbor_l
      end if
    end if
    if (present(added)) added = .true.
    link%task%n_nbors = link%task%n_nbors + 1
    !$omp atomic update
    link%task%n_needed = link%task%n_needed + 1
  else
    if (present(added)) added = .false.
    write (stderr,*) 'WARNING: trying to add existing nbor', next%task%id
  end if
  if (present(nbor)) then
    nbor => nbor_l
  end if
  call link%lock%unset ('add_nbor')
  !call trace%end (itimer)
END SUBROUTINE add_nbor

!===============================================================================
SUBROUTINE init_ids (list)
  class(list_t):: list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer, save :: level=0
  !-----------------------------------------------------------------------------
  call list%lock%set ('init_ids')
  link => list%head
  do while (associated(link))
    task => link%task
    call task%init_unique (same=.true.) ! obtain a unique ID
    !if (task%level > level) then
    !  level = task%level
    !  list%first(level)%pointer => link
    !endif
    link => link%next
  end do
  call list%lock%unset ('init_ids')
END SUBROUTINE init_ids

!===============================================================================
SUBROUTINE remove_external (list)
  class(list_t):: list
  class(link_t), pointer:: link, next
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  ! Remove tasks that are external, and make sure to also remove them from the
  ! nbor lists of virtual tasks
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_mod%remove_external', itimer=itimer)
  link => list%head
  do while (associated(link))
    next => link%next
    if (link%task%is_set (bits%external)) then
      !call link%remove_nbor_list2 (link%nbor)
      call list%remove (link)
      call link%task%dealloc()
      !deallocate (link%task)
      !deallocate (link)
    end if
    link => next
  end do
  call trace%end(itimer)
END SUBROUTINE remove_external

!===============================================================================
SUBROUTINE print_list (list)
  class(list_t):: list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_mod%print_list', itimer=itimer)
  call list%lock%set ('print_list')
  write (io_unit%queue, *) 'print_list'
  link => list%head
  do while (associated(link))
    task => link%task
    select type(task)
    class is (m3d_patch_t)
      write (io_unit%queue, *) task%is_set(bits%external), &
                               task%is_set(bits%internal), &
                               task%is_set(bits%boundary), &
                               task%is_set(bits%virtual), &
                               task%m3d_patch_info()
                              
    class is (experiment_t)
      write(io_unit%queue, '(4(L2), 1x, I3, 1x, I3, 1x, F8.4, A12, I8, I4)') &
                    task%is_set(bits%external), &
                    task%is_set(bits%internal), &
                    task%is_set(bits%boundary), &
                    task%is_set(bits%virtual), &
                    omp%place, task%place, task%time, &
                    trim(task%kind), task%id, task%ipos(3)+1
    end select
    link => link%next
  end do
  call list%lock%unset ('print_list')
  flush(io_unit%queue)
  call trace%end(itimer)
END SUBROUTINE print_list

END MODULE