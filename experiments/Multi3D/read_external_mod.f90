!===============================================================================
!> Reads a tabulated model atmosphere from a text file
!===============================================================================
MODULE read_external_mod
  USE io_mod
  USE os_mod
  USE m3io_mod, only: NextEntry
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
  USE experiment_mod
  USE ion_mod, only: vmic
  USE composition_mod
  USE composition_patch_mod, only: FeH
  USE patch_mod
  USE mconst, only:cc
  USE atom_mod, only: m1d_legacy_mode
  USE ray_patch_mod, only: maxvel
  USE m3d_patch_mod, only: maxiter, ilambd, m3d_patch_t, atom_update_counter, &
  global_maxerr, global_pos, nz_total, use_density, use_ne
  USE m3d_halo_mod, only: m3d_halo_t
  USE read_stagger_snapshot_mod
  USE read_marcs_mod, only: read_marcs, zfac
  implicit none

  type, public:: read_external_t
    class(task_list_t), pointer:: task_list
    integer :: dims(3)=[1,1,1], n(3), verbose
    real(4), dimension(:,:,:), pointer, contiguous :: temp, ne, pe, rho, px, py, pz, ee
    real(8) :: ds(3), size(3), position(3), origin(3)
  contains
    procedure:: init
    procedure:: init_exp
    procedure:: read_1d_snapshot
    procedure:: read_multi_snapshot
    procedure:: setup_patch
    procedure:: load_authority
    procedure:: save_atmos
    procedure:: save_multi_atmos
  end type read_external_t

  type, public, extends(task_t):: rank_t
  end type

  type(read_external_t), public:: read_external
  type(rank_t) :: rank
  character(len=126), public :: atmos_file='./input_multi3d/atmos/atmos.sun_MARCS_50'

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
  real(4) :: origin(3)=[0.0,0.0,0.0]
  integer :: iostat, k, ip
  integer :: patches_per_mpi=1
  integer :: dims=1, nx=5, ny=5, nz=256
  integer :: verbose=1
  character(len=126) :: atmos_format='Multi'
  real(4), dimension(:,:,:), pointer, contiguous :: vel
  namelist /atmos_params/ atmos_file, dims, verbose, patches_per_mpi, FeH, &
  atmos_format, vmic, use_density, use_ne, nx, ny, nz
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
  write (io%namelist, atmos_params)

  if (atmos_format == 'Stagger') then
    call read_stagger_params(nx, ny, nz)
    if (irregular) then
      dims = size(stagger_levels)
    endif
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
  self%origin = origin
  call io%lock%unset

  allocate (task_list)
  self%task_list => task_list

  if (irregular .and. atmos_format == 'Stagger') then
    call read_stagger_size(self%n, self%size)
    self%size = self%size / scaling%l

  else
    select case (atmos_format)
    case ('Stagger')
      call read_stagger_snapshot(self%ds, self%rho, self%px, self%py, &
      self%pz, self%ee, self%temp, self%ne)
      atmos_file = stagger_file
    case ('Marcs')
      call read_marcs(atmos_file, nz, self%ds, self%rho, self%temp, self%pe, FeH)
    case ('Multi')
      call self%read_multi_snapshot(atmos_file, self%ds, self%rho, self%temp, &
                                    self%ne, self%px, self%py, self%pz)
    case ('Text')
      call self%read_1d_snapshot(atmos_file, self%ds, self%rho, self%temp, self%pe)
    case ('MUST')
      call self%read_multi_snapshot(atmos_file, self%ds, self%rho, self%temp, &
                                    self%ne, self%px, self%py, self%pz)
      self%pz = - self%pz
    case default
      write(stdout, *) 'unknown atmos_format ', atmos_format, &
      '. Options are: "Stagger", "Multi", "Marcs" and "Text".'
      stop
    end select

    self%n = shape(self%rho)
    if (.not. associated(self%pe)) then
      call os%allocate(self%pe, self%n)
      self%pe = self%ne * cgs%k_b * self%temp
    elseif (.not. associated(self%ne)) then
      call os%allocate(self%ne, self%n)
      self%ne = self%pe / (cgs%k_b * self%temp)
    endif

    self%ds  = self%ds  / scaling%l
    self%pe  = self%pe  / scaling%p
    self%rho = self%rho / scaling%d

    if (.not. associated(self%pz)) then
      call os%allocate(self%px, self%n, zero=.true.)
      call os%allocate(self%py, self%n, zero=.true.)
      call os%allocate(self%pz, self%n, zero=.true.)
    else
      self%px = self%px / (scaling%d * scaling%u)
      self%py = self%py / (scaling%d * scaling%u)
      self%pz = self%pz / (scaling%d * scaling%u)
    endif

    if (self%rho(1,1,1) < self%rho(1,1,self%n(3))) then
      self%temp =  self%temp(:,:,self%n(3):1:-1)
      self%rho  =  self%rho (:,:,self%n(3):1:-1)
      self%ne   =  self%ne  (:,:,self%n(3):1:-1)
      self%pe   =  self%pe  (:,:,self%n(3):1:-1)
      self%px   =  self%px  (:,:,self%n(3):1:-1)
      self%py   =  self%py  (:,:,self%n(3):1:-1)
      self%pz   = -self%pz  (:,:,self%n(3):1:-1)
    endif

    call os%allocate(vel, self%n)
    vel = sqrt(self%px**2 + self%py**2 + self%pz**2) / self%rho
    maxvel = scaling%u / cc * maxval(vel)
    call os%deallocate(vel)

    maxvel = max(1e-8, maxvel)

    self%size = self%n * self%ds
  endif

  self%position = 0.5*self%size + self%origin

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
  rank%position = (mpi_coords%rank_to_coords(mpi%rank) + 0.5_8)*rank%size + self%origin
  !-----------------------------------------------------------------------------
  ! Begin on task list
  !-----------------------------------------------------------------------------
  call self%task_list%init ('read_external_mod')
  self%task_list%size = self%size
  self%task_list%dims = self%dims
  self%task_list%position = self%position
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

    if (k==0) write(stdout,'(I2, A28)') self%dims(3), ' subdomains in z-direction:'
    write(stdout,'(I4, I6,A3,I4,4X,3I4, 2X, F7.4, E14.4)') k+1, patch%zs+1, ' - ', patch%ze, &
    patch%n, (patch%position(3)-patch%origin(3))/patch%box(3), patch%ds(3)
  end do
  
  nz_total = patch%ze
  
  write (stdout, '(A, 3(F10.4))') 'box = ', patch%box
  write (stdout,*) ip, 'tasks generated'

  if (.not. irregular) then
    call self%save_atmos(patch)
    call self%save_multi_atmos()
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
  write (stdout,'(a,i8)') ' number of tasks in task list:', self%task_list%n
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

  call os%deallocate(self%temp)
  call os%deallocate(self%rho)
  call os%deallocate(self%ne)
  call os%deallocate(self%pe)
  call os%deallocate(self%px)
  call os%deallocate(self%py)
  call os%deallocate(self%pz)
  
  call print_list(task_list)
  !-----------------------------------------------------------------------------
  ! Refresh the counts, now with bits%has_output set on tasks
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdlog,*) 'read_external_mod: ntask,nwrite=', io%ntask, io%nwrite
  if (verbose > 0) write(stdout,*) 'read_external_t%init: print tasks'
  call self%task_list%print_tasks
  call self%task_list%check_hash_table
  !call self%task_list%hash_table%stats
  if (verbose > 0) write(stdout,*) 'read_external_t%init: done'
  flush (stdout)
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE init

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
  patch%origin = self%origin

  if (.not. irregular) then
    patch%ds = self%ds
  endif
  
  call patch%pre_init()
  call trace%end ()
END SUBROUTINE setup_patch

!===============================================================================
SUBROUTINE init_exp (self)
  class(read_external_t):: self
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  logical:: periodic(3)
  integer:: id
  !-----------------------------------------------------------------------------
  id = 0
  link => self%task_list%head
  do while (associated(link))
    id = id+1
    if (mod(id,omp%nthreads)==omp%thread) then
      task => link%task
      select type (task)
      class is (experiment_t)
        call task%init()
        if (.not. irregular) then
          call self%load_authority(task)
        endif
      end select
    endif
    link => link%next
  end do
END SUBROUTINE init_exp

!===============================================================================
!> Reads the actual text files and stores variables under self
!===============================================================================
SUBROUTINE read_1d_snapshot(self, atmos_file, ds, rho, temp, pe)
  class(read_external_t):: self
  real(4), intent(out), dimension(:,:,:), pointer, contiguous :: rho, temp, pe
  real(8), intent(out), dimension(:) :: ds
  integer :: nz, lunit, i, iostat
  character(len=126) :: atmos_file
  real(4), dimension(:), allocatable :: zi
  real(4) :: vturb
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  call os%check_file(atmos_file)
  open(file=atmos_file, newunit=lunit)
  call nextentry(lunit)
  read(lunit, *) ! skip name
  call nextentry(lunit)
  read(lunit, *) nz

  allocate(zi(nz))
  call os%allocate(temp, [1,1,nz])
  call os%allocate(rho , [1,1,nz])
  call os%allocate(pe  , [1,1,nz])

  do i=1, nz
    call nextentry(lunit)
    read(lunit, *) zi(i), temp(1,1,i), pe(1,1,i), rho(1,1,i), vturb
  end do

  close(lunit)

  ds(1) = 1d10
  ds(2) = 1d10 
  ds(3) = zi(2) - zi(1)

  if (abs(ds(3) - (zi(nz) - zi(nz-1))) > 1d-4) then
    write(stdout,*) 'read_external_t%read_1d_snapshot: Warning: Non-equidistant atmosphere grid!'
  endif
END SUBROUTINE read_1d_snapshot

!===============================================================================
SUBROUTINE read_multi_snapshot(self, atmos_file, ds, rho, temp, ne, px, py, pz)
  class(read_external_t):: self
  real(4), intent(out), dimension(:,:,:), pointer, contiguous :: rho, temp, ne, px, py, pz
  real(8), intent(out), dimension(:) :: ds
  integer :: lunit, n(3)
  character(len=126) :: atmos_file
  real(4) :: test
  real(4), dimension(:), allocatable :: xi, yi, zi
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  call os%check_file(trim(atmos_file) // '.mesh')
  open(file=trim(atmos_file) // '.mesh', newunit=lunit, form='FORMATTED')
  read(lunit, *) n
  allocate(xi(n(1)))
  allocate(yi(n(2)))
  allocate(zi(n(3)))
  read(lunit, *) xi
  read(lunit, *) yi
  read(lunit, *) zi
  close(lunit)

  ds(1) = xi(2) - xi(1)
  ds(2) = yi(2) - yi(1)
  ds(3) = zi(2) - zi(1)

  if (abs(ds(3) / (zi(n(3)) - zi(n(3)-1)) - 1) > 1d-4) then
    write(stdout,*) 'read_external_t%read_1d_snapshot: Warning: Non-equidistant atmosphere grid! Stopping!'
    stop
  endif

  ds = abs(ds)

  call os%allocate(temp, n)
  call os%allocate(rho, n)
  call os%allocate(ne, n)
  call os%allocate(px, n)
  call os%allocate(py, n)
  call os%allocate(pz, n)

  call os%check_file(trim(atmos_file) // '.atmos')
  open(file=trim(atmos_file) // '.atmos', newunit=lunit, form='UNFORMATTED', access='stream')
  read(lunit) ne 
  read(lunit) temp
  read(lunit) px
  read(lunit) py
  read(lunit) pz
  read(lunit) rho
  close(lunit)

  px =   px * rho
  py =   py * rho
  pz = - pz * rho ! z-velocity is flipped in Multi format
END SUBROUTINE read_multi_snapshot

!===============================================================================
SUBROUTINE save_multi_atmos (self)
  class (read_external_t) :: self
  real(4), dimension(:,:,:), allocatable :: ne, temp, vx, vy, vz, rho
  real(8) :: ds(3)
  integer :: lunit, i, j, recl
  !-----------------------------------------------------------------------------
  open(file=trim(io%run_dir)//'mesh.multi', newunit=lunit, form='FORMATTED')
  ds = self%ds * scaling%l
  do i=1, 3
    write(lunit, *) self%n(i)
    write(lunit, '(E14.8)') [(ds(i)*j, j=0, self%n(i)-1)]
  end do
  close(lunit)

  allocate(ne, temp, vx, vy, vz, rho, mold=self%temp)

  ne   = self%ne  (:,:,self%n(3):1:-1)
  temp = self%temp(:,:,self%n(3):1:-1)
  vx   = self%px  (:,:,self%n(3):1:-1) / scaling%u
  vy   = self%py  (:,:,self%n(3):1:-1) / scaling%u
  vz   = self%pz  (:,:,self%n(3):1:-1) / scaling%u
  rho  = self%rho (:,:,self%n(3):1:-1) * scaling%d

  vx = vx / rho
  vy = vy / rho
  vz = vz / rho

  recl = kind(1.0)*product(self%n)

  open(file=trim(io%run_dir)//'atm3d.multi', newunit=lunit, form='UNFORMATTED', access='direct', recl=recl)
  write(lunit, rec=1) ne  
  write(lunit, rec=2) temp
  write(lunit, rec=3) vx  
  write(lunit, rec=4) vy  
  write(lunit, rec=5) vz  
  write(lunit, rec=6) rho 
  close(lunit)
END SUBROUTINE save_multi_atmos

!===============================================================================
!> Takes variables from self and loads them into patch memory
!===============================================================================
SUBROUTINE load_authority(self, patch)
  class (read_external_t) :: self
  class (experiment_t) :: patch
  integer :: i,j,k, it, iz, xs, xe, ys, ye, zs, ze
  !-----------------------------------------------------------------------------
  it = 1 ! time index, it=1 is the initial condition
  xs = patch%li(1); xe = patch%ui(1)
  ys = patch%li(2); ye = patch%ui(2)
  zs = patch%li(3); ze = patch%ui(3)

  do k=zs,ze
    iz = patch%zs+1 + k - zs

    patch%mem(xs:xe,ys:ye,k,patch%idx%d ,it,1) = self%rho(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%px,it,1) = self%px(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%py,it,1) = self%py(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%pz,it,1) = self%pz(:,:,iz)

    if (associated(self%temp)) then
      patch%mem(xs:xe,ys:ye,k,patch%idx%e,it,1) = self%temp(:,:,iz)
    elseif (associated(self%ee)) then
      patch%mem(xs:xe,ys:ye,k,patch%idx%e,it,1) = self%ee(:,:,iz)
    endif
    
    patch%mem(xs:xe,ys:ye,k,6,it,1) = self%pe(:,:,iz)
  end do

  if (associated(self%temp)) then
    patch%v5 = 'tt'
  elseif (associated(self%ee)) then
    patch%v5 = 'ee'
  endif

  call os%allocate(patch%filled, patch%gn)
  patch%filled = -1
END SUBROUTINE

!===============================================================================
SUBROUTINE init_all_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: i, i0, i1, i2, m, n=0, nv=0, nb=0
  real(8):: time
  integer, save:: itimer=0
  integer :: fixed_threads
  !.............................................................................
  if (self%hands_off_nbor_lists) return
  call trace%begin('read_external_mod%init_all_nbors', itimer=itimer)
  if (self%verbose > 0) &
    call timer%tic (time)
  fixed_threads = min(omp%nthreads, 4)
  i2 = self%n/fixed_threads + 1                         ! tasks per thread
  n = 0
  !$omp parallel private(i,m,i0,i1,link) shared(self,i2,n,stdout) default(none) num_threads(fixed_threads)
  i0 = omp%thread*i2
  i1 = i0 + i2
  if (omp%master) &
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
SUBROUTINE save_atmos (self, patch)
  class (read_external_t) :: self
  class(experiment_t), pointer:: patch
  class(composition_type), pointer:: composition
  real(4), dimension(:,:,:,:), allocatable :: chi, abslte, sumsca, eta
  real(4), dimension(:,:,:), allocatable :: tau
  real(8), dimension(1) :: nu
  real(8) :: dz
  character(len=82) :: out_atmos_file='atmos_data.bin'
  character(len=82) :: out_mesh_file='atmos_mesh.txt'
  integer :: lunit, n,step, m(3), i, j, k
  !-----------------------------------------------------------------------------
  open(file=trim(io%run_dir)//out_mesh_file, newunit=lunit, form='FORMATTED')
  write(lunit, *) self%n
  write(lunit, *) self%ds
  close(lunit)
  
  n = self%n(3)
  m = shape(self%temp)
  allocate(chi(m(1), m(2), m(3),1))
  allocate(abslte, sumsca, eta, mold=chi)
  allocate(tau, mold=self%temp)
  nu = 1e8 * cc / 5000.0
  dz = self%ds(3) * scaling%l

  composition => patch%composition%composition
  call composition%load_state_variables(dble(self%temp), dble(self%ne), &
      dble(self%rho)*scaling%d, use_rho=use_density, use_ne=use_ne)
  call composition%absorption(nu, chi, abslte, sumsca, sorted=.true.)
  call composition%scattering(nu, abslte, eta, sumsca)
  call compute_tau_scale(chi(:,:,n:1:-1,1), dz, tau(:,:,n:1:-1))

  open(file=trim(io%run_dir)//out_atmos_file, newunit=lunit, form='UNFORMATTED')
  write(lunit) self%temp
  write(lunit) self%rho 
  write(lunit) self%px  
  write(lunit) self%py  
  write(lunit) self%pz  
  write(lunit) self%pe  
  write(lunit) tau      
  write(lunit) chi      
  write(lunit) eta   
  close(lunit)
END SUBROUTINE save_atmos

!===============================================================================
SUBROUTINE compute_tau_scale (chi, ds, tau)
  real(4), dimension(:,:,:), intent(in) :: chi
  real(4), dimension(:,:,:), intent(out):: tau
  real(8) :: ds
  integer :: k
  !-----------------------------------------------------------------------------
  tau(:,:,2) = 0.5*(chi(:,:,2) + chi(:,:,1)) * ds
  tau(:,:,3) = 0.5*(chi(:,:,3) + chi(:,:,2)) * ds + tau(:,:,2)
  tau(:,:,1) = exp(2*log(tau(:,:,2))-log(tau(:,:,3)))

  do k=2, size(chi,3)
    tau(:,:,k) = tau(:,:,k-1) + 0.5*(chi(:,:,k) + chi(:,:,k-1)) * ds
  end do
END SUBROUTINE compute_tau_scale

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
  call trace%end(itimer)
END SUBROUTINE print_list

END MODULE