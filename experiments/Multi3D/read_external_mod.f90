!===============================================================================
!> Reads a tabulated model atmosphere from a text file
!===============================================================================
MODULE read_external_mod
  USE io_mod
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
  USE task_list_mod
  USE experiment_mod
  USE patch_mod
  USE ray_patch_mod
  USE m3d_patch_mod, only: maxiter, ilambd
  USE read_stagger_snapshot_mod
  implicit none

  type, public:: read_external_t
    class(task_list_t), pointer:: task_list
    integer :: dims(3), n(3), ng(3), verbose
    real(4), dimension(:,:,:), pointer, contiguous :: temp, ne, pe, rho, vturb, px, py, pz, ee
    real(4), dimension(:), pointer, contiguous :: xi, yi, zi
    real(4) :: ds(3), size(3), position(3)!, origin(3)
    character(len=126) :: atmos_name
  contains
    procedure:: init
    procedure:: init_exp
    procedure:: read_1d_snapshot
    procedure:: read_multi_snapshot
    procedure:: setup_patch
    procedure:: load_authority
  end type read_external_t

  type, public, extends(task_t):: rank_t
  end type

  type(read_external_t), public:: read_external
  type(rank_t) :: rank
  logical, save :: read_Stagger=.false.
  logical, save :: read_Multi=.false.

CONTAINS

SUBROUTINE init(self, task_list, name)
  !-----------------------------------------------------------------------------
  class(read_external_t):: self
  class(task_list_t),  pointer:: task_list
  character(len=*),   optional:: name
  !-----------------------------------------------------------------------------
  class(experiment_t), pointer:: patch
  class(link_t), pointer:: link, next
  !real(4) :: origin(3)=[0.0,0.0,0.0]
  integer :: iostat, k, ip
  integer :: mpi_dims(3)=[1,1,1]
  integer :: patches_per_mpi(3)=[1,1,1]
  integer :: dims(3)=[1,1,1]
  integer :: verbose=1
  logical :: face_nbors=.false.
  character(len=126) :: atmos_file='./input_multi3d/atmos/atmos.sun_MARCS_50'
  namelist /atmos_params/ atmos_file, dims, verbose, read_Stagger, read_Multi, mpi_dims
  !-----------------------------------------------------------------------------
  if (associated(task_list)) return
  call trace%begin('read_external_t%init')
  call io%header ('read_external_t%init')
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

  self%dims = dims
  self%verbose = verbose
  !self%origin = origin
  call io%lock%unset

  allocate (task_list)
  self%task_list => task_list
  call scaling%init()
  
  if (read_Stagger) then
    call read_stagger_snapshot(self%ds, self%rho, self%px, self%py, &
    self%pz, self%ee, tt=self%temp, ne=self%ne)
    self%n = shape(self%rho)
  elseif (read_Multi) then
    call self%read_multi_snapshot(atmos_file)
  else
    call self%read_1d_snapshot(atmos_file)
  endif
  self%size = self%n * self%ds
  self%position = 0.5*self%size! + self%origin

  !===============================================================================
  if (mpi%size==1) &
    mpi_dims = 1
  !-----------------------------------------------------------------------------
  ! Generate a (periodic) MPI geometry, with dims specified in the input file,
  ! or assigned automatically by mpi_coords%init.
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) then
     mpi_coords%dims = 1
  else
    if (product(mpi_dims)==1) then
      call mpi_coords%init (dims=self%dims)
    else
      call mpi_coords%init (mpi_dims=mpi_dims)
    end if
  end if
  mpi_dims = mpi_coords%dims
  patches_per_mpi = self%dims/mpi_dims
  mpi_coords%npatch = patches_per_mpi
  io%mpi_dims = mpi_dims
  !-----------------------------------------------------------------------------
  ! Set attritbutes of the rank_t data type -- just enough to use distance().
  !-----------------------------------------------------------------------------
  rank%periodic = .true.
  rank%box = self%size
  rank%size = rank%box/mpi_dims
  rank%ds = rank%size/patches_per_mpi
  rank%position = (mpi_coords%rank_to_coords(mpi%rank) + 0.5_8)*rank%size! + self%origin
  !-----------------------------------------------------------------------------
  ! Begin on task list
  !-----------------------------------------------------------------------------
  call self%task_list%init ('read_external_mod')
  self%task_list%size = self%size
  self%task_list%dims = self%dims
  self%task_list%position = self%position
  self%task_list%n_tasks = product(self%dims)
  self%task_list%face_nbors = face_nbors

  !-----------------------------------------------------------------------------
  ! Allocate an array of links, pointing to patches, and give the patches each
  ! id, rank, size, position, integer position, and status bits.  Discard
  ! external patches and add the rest to the task_list.
  !-----------------------------------------------------------------------------
  ip = 0
  do k=0, self%dims(3)-1
    allocate (patch)
    ! --- patch 3D and 1D coordinates ---
    patch%ipos = [0,0,k]
    patch%rank = mpi_coords%coords_to_rank ([0,0,k]/patches_per_mpi)
    if (patch%rank == mpi%rank) then
      ip = ip + 1
      patch%ip = ip
    end if
    if (all(self%dims==1)) call patch%task_t%set(bits%root_grid)
    call self%setup_patch(patch)
    call patch%set (bits%virtual)
    call patch%set (bits%static)

    if (verbose > 0) then 
      if (k==0) write(stdout,'(I2, A28)') dims(3), ' subdomains in z-direction:'
      write(stdout,'(I4, I6,A3,I4)') k+1, patch%zs, ' - ', patch%ze
    endif

  end do
  write (stdout,*) ip, 'task generated'

  io%end_time = (maxiter + ilambd) * patch%dtime
  io%print_every = huge(1)
  timer%dead_mans_hand = huge(1.0)

  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_levels'
  call self%task_list%init_levels
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init_all_nbors'
  call self%task_list%init_all_nbors
  if (verbose > 0) write(stdout,*) 'read_external_t%init: reset_status'
  call self%task_list%reset_status
  if (verbose > 0) write(stdout,*) 'read_external_t%init: info'
  call task_list%info

  !-----------------------------------------------------------------------------
  ! Remove tasks that are external, and make sure to also remove them from the
  ! nbor lists of virtual tasks
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    next => link%next
    if (link%task%is_set (bits%external)) then
      call link%remove_nbor_list2 (link%nbor)
      call self%task_list%remove (link)
      call link%task%dealloc()
      deallocate (link%task)
      deallocate (link)
    end if
    link => next
  end do
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
  ! Initialize tasks
  !-----------------------------------------------------------------------------
  !$omp parallel default(shared)
  call self%init_exp()
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! Refresh the counts, now with bits%has_output set on tasks
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdlog,*) 'read_external_mod: ntask,nwrite=', io%ntask, io%nwrite
  !-----------------------------------------------------------------------------
  ! This call is not valid until the patch mem and meshes have been allocated
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(stdout,*) 'read_external_t%init: init boundaries'
  call self%task_list%init_bdries
  !-----------------------------------------------------------------------------
  !if (verbose > 0) write(stdout,*) 'read_external_t%init: diagnostics'
  !call self%diagnostics()
  if (verbose > 0) write(stdout,*) 'read_external_t%init: print tasks'
  call self%task_list%print_tasks
  call self%task_list%check_hash_table
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
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%setup_patch')
  patch%task_list => self%task_list

  call patch%init_unique (same=.true.) ! obtain an unique ID
  call patch%gpatch_t%pre_init() ! this overwrites patch%n
  patch%pre_initialized = .true.

  patch%mhd = .false.
  patch%kind = 'mhd_patch'

  patch%n = self%n / self%dims
  patch%ng = 1 ! needed for nbor_relation
  patch%zs = patch%ipos(3) * patch%n(3) + 1
  if (patch%ipos(3)+1 == self%dims(3)) then
    patch%n(3) = patch%n(3) + mod(self%n(3), self%dims(3))
  endif
  patch%ze = patch%zs + patch%n(3) - 1

  patch%box = self%size
  patch%ds  = self%ds

  patch%size = patch%ds * patch%n
  patch%periodic = [.true., .true., .false.]

  patch%position = 0.5d0*patch%size
  patch%position(3) = patch%position(3) + (patch%zs - 1) * patch%ds(3)
  patch%origin = patch%position - 0.5d0*patch%size

  !-----------------------------------------------------------------------------
  ! Preliminary region information, relative to the task position.
  ! The settings below allow overlap() calls to work
  !-----------------------------------------------------------------------------
  patch%authority%l = - patch%size*0.5_8
  patch%authority%u = + patch%size*0.5_8
  patch%interest%l  = patch%authority%l - patch%ng*patch%size/patch%n
  patch%interest%u  = patch%authority%u + patch%ng*patch%size/patch%n

  allocate (patch%boundaries)
  call patch%boundaries%pre_init (patch)
  if (patch%ipos(3)==0)              call patch%boundaries%set(bits%zl)
  if (patch%ipos(3)==self%dims(3)-1) call patch%boundaries%set(bits%zu)
  
  call patch%pre_init()
  
  allocate(patch%composition)
  patch%composition%patch => patch
  call patch%composition%pre_init()
  
  patch%nv = 6
  patch%nt = 4 ! Download returns wrong stuff if nt <= 2
  patch%nw = 1

  call patch%idx%init (energy=5, mhd=patch%mhd)
  !---------------------------------------------------------------------------
  ! Add the patch to the task list -- this creates the patch%link
  !---------------------------------------------------------------------------
  call patch%append_to(patch%task_list)

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
        call self%load_authority(task)
      end select
    endif
    link => link%next
  end do
END SUBROUTINE init_exp

!===============================================================================
!> Reads the actual text files and stores variables under self
!===============================================================================
SUBROUTINE read_1d_snapshot(self, atmos_1d_file)
  !-----------------------------------------------------------------------------
  USE m3io_mod, only: NextEntry
  class(read_external_t):: self
  integer :: nz, lunit, i, iostat
  character(len=126) :: atmos_1d_file
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  open(file=atmos_1d_file, newunit=lunit)
  call nextentry(lunit)
  read(lunit, *) self%atmos_name
  call nextentry(lunit)
  read(lunit, *) nz

  self%n(3) = nz - 2*self%ng(3)

  allocate(self%zi(1-self%ng(3):self%n(3)+self%ng(3)))
  allocate(self%temp(1,1, 1-self%ng(3):self%n(3)+self%ng(3)))
  allocate(self%px, self%py, self%pz, self%pe, self%rho, self%vturb, mold=self%temp)

  do i=lbound(self%zi,1), ubound(self%zi,1)
    call nextentry(lunit)
    read(lunit, *) self%zi(i), self%temp(1,1,i), self%pe(1,1,i), self%rho(1,1,i), self%vturb(1,1,i)
  end do

  close(lunit)

  self%zi  = self%zi  / scaling%l
  self%rho = self%rho / scaling%d
  self%pe  = self%pe  / scaling%p
  self%px  = 0
  self%py  = 0
  self%pz  = 0

  self%n(:2) = 1
  self%ds(1) = 1e-2 ! 10km spacing as for 1d atmospheres in Multi
  self%ds(2) = 1e-2
  self%ds(3) = self%zi(2) - self%zi(1)

  if (abs(self%ds(3) - (self%zi(self%n(3)) - self%zi(self%n(3)-1))) > 1d-4) then
    write(stdout,*) 'read_external_t%init: Warning: Non-equidistant atmosphere grid!'
  endif
END SUBROUTINE read_1d_snapshot

!===============================================================================
SUBROUTINE read_multi_snapshot(self, stagger_file)
  class(read_external_t):: self
  integer :: lunit
  character(len=126) :: stagger_file
  real(4) :: test
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  open(file=trim(stagger_file) // '.mesh', newunit=lunit, form='FORMATTED')
  read(lunit, *) self%n
  allocate(self%xi(self%n(1)))
  allocate(self%yi(self%n(2)))
  allocate(self%zi(self%n(3)))
  read(lunit, *) self%xi
  read(lunit, *) self%yi
  read(lunit, *) self%zi
  close(lunit)

  allocate(self%temp(self%n(1), self%n(2), self%n(3)))
  allocate(self%rho, self%ne, self%px, self%py, self%pz, self%pe, mold=self%temp)

  open(file=trim(stagger_file) // '.atmos', newunit=lunit, form='UNFORMATTED', access='stream')
  read(lunit) self%ne 
  read(lunit) self%temp
  read(lunit) self%px
  read(lunit) self%py
  read(lunit) self%pz
  read(lunit) self%rho
  close(lunit)

  self%xi  = self%xi  / scaling%l
  self%yi  = self%yi  / scaling%l
  self%zi  = self%zi  / scaling%l
  self%rho = self%rho / scaling%d
  self%px  = self%px * self%rho / scaling%u
  self%py  = self%py * self%rho / scaling%u
  self%pz  = self%pz * self%rho / scaling%u
  self%pe  = self%ne * cgs%k_b * self%temp  / scaling%p

  self%ds(1) = self%xi(2) - self%xi(1)
  self%ds(2) = self%yi(2) - self%yi(1)
  self%ds(3) = self%zi(2) - self%zi(1)
END SUBROUTINE read_multi_snapshot

!===============================================================================
!> Takes variables from self and loads them into patch memory
!===============================================================================
SUBROUTINE load_authority(self, patch)
  !-----------------------------------------------------------------------------
  class (read_external_t) :: self
  class (experiment_t) :: patch
  integer :: i,j,k, it, iz, xs, xe, ys, ye, zs, ze

  it = 1 ! time index, it=1 is the initial condition
  xs = patch%li(1); xe = patch%ui(1)
  ys = patch%li(2); ye = patch%ui(2)
  zs = patch%li(3); ze = patch%ui(3)

  do k=zs,ze
    if (read_Stagger) then
      iz = patch%zs + k - zs
    else
      ! flipping z-axis, so that the first patch is the bottom of the atmosphere
      iz = self%n(3) - (patch%zs - 1 + k - zs) 
    endif

    patch%mem(xs:xe,ys:ye,k,patch%idx%d ,it,1) = self%rho(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%px,it,1) = self%px(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%py,it,1) = self%py(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%pz,it,1) = self%pz(:,:,iz)

    if (associated(self%temp)) then
      patch%mem(xs:xe,ys:ye,k,patch%idx%e ,it,1) = self%temp(:,:,iz)
    elseif (associated(self%ee)) then
      patch%mem(xs:xe,ys:ye,k,patch%idx%e ,it,1) = self%ee(:,:,iz)
    endif
    
    if (associated(self%pe)) then
      patch%mem(xs:xe,ys:ye,k,6,it,1) = self%pe(:,:,iz)
    elseif (associated(self%ne) .and. associated(self%temp)) then
      patch%mem(xs:xe,ys:ye,k,6,it,1) = self%ne(:,:,iz) * cgs%k_b * self%temp(:,:,iz)  / scaling%p
    endif
  end do

  if (associated(self%temp)) then
    patch%v5 = 'tt'
  elseif (associated(self%ee)) then
    patch%v5 = 'ee'
  endif

  allocate(patch%filled(patch%gn(1), patch%gn(2), patch%gn(3)))
  patch%filled = -1

END SUBROUTINE

END MODULE