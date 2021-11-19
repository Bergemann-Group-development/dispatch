!===============================================================================
!> Reads a tabulated model atmosphere from a text file
!===============================================================================
MODULE read_external_mod
  USE io_mod
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
  USE download_mod
  USE m3d_rt_mod

  type, public:: read_external_t
    character(len=126) :: atmos_file
    class(task_list_t), pointer:: task_list
    integer :: nz, npatches(3), n(3), ng(3), verbose
    real(8), dimension(:), allocatable :: zi, temp, pe, rho, vturb
    real(8) :: size, ds
    character(len=126) :: atmos_name
  contains
    procedure:: init
    procedure:: read_snapshot
    procedure:: create_patch
    procedure:: load_authority
    procedure:: fill_guard_zones
  end type read_external_t

  type(read_external_t), public:: read_external

CONTAINS

SUBROUTINE init(self, task_list, name)
  !-----------------------------------------------------------------------------
  class(read_external_t):: self
  class(task_list_t),  pointer:: task_list
  character(len=*),   optional:: name
  !-----------------------------------------------------------------------------
  class(task_t),       pointer:: task
  class(experiment_t), pointer:: patch
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  logical :: searching
  integer :: iostat, npatches, iz, i, ng(3), verbose, itaskl
  character(len=126) :: atmos_file
  namelist /atmos_params/ atmos_file, npatches, ng, verbose
  !-----------------------------------------------------------------------------
  if (associated(task_list)) return
  call trace%begin('read_external_t%init')
  call io%init ('read_external_t%init')
  call io%header ('read_external_t%init')
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and saves them
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    ! default values
    atmos_file = './input_multi3d/atmos/atmos.sun_MARCS'
    npatches = 1
    ng = [4,4,1]
    verbose = 1
    
    rewind (io%input)
    read (io%input, nml=atmos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, atmos_params)
    end if
    call io%namelist_warning ('atmos_params', iostat)
    write (io%namelist, atmos_params)

    self%atmos_file = atmos_file
    self%npatches(1) = 1
    self%npatches(2) = 1
    self%npatches(3) = npatches
    self%ng = ng
    self%verbose = verbose
  endif
  call io%lock%unset

  allocate (task_list)
  self%task_list => task_list
  call self%task_list%init ('multi3d_mod')

  call scaling%init()
  call self%read_snapshot()

  self%size = self%zi(self%nz) - self%zi(1)

  if (mod(self%nz, npatches) /= 0) then
    write(stdout,*) 'read_external_t%init: Error: npatches must be a factor of the atmosphere size!'
    write(stdout,*) 'npatches = ', self%npatches(3), 'nz = ', self%nz
    write(stdout,*) 'Stopping.'
    stop
  endif

  do i=1, npatches
    call self%create_patch(i, patch) ! patch is allocated here
    call self%load_authority(patch, i)

    if (verbose > 0) then 
      iz = (i-1)*patch%n(3) + 1
      if (i==1) write(stdout,'(I2, A28)') npatches, ' subdomains in z-direction:'
      write(stdout,'(I4, I6,A3,I4)') i, iz, ' - ', iz + patch%n(3) - 1
    endif
  end do

  !-----------------------------------------------------------------------------
  ! Add all tasks to the ready queue, preventing a slow default add_by_time
  !-----------------------------------------------------------------------------
  itaskl = 0
  link => self%task_list%head
  do while(associated(link))
    itaskl = itaskl + 1
    link => link%next
  enddo
  
  self%task_list%n_tasks = itaskl

  call self%task_list%queue_all
  print *, self%task_list%n_tasks
  !call self%task_list%print_queue

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
  if (verbose > 0) write(stdout,*) 'read_external_t%init: fill_guard_zones'
  call self%fill_guard_zones

  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize a minimum amount of patch information in a Cartesian arrangement
!===============================================================================
SUBROUTINE create_patch (self, i, patch)
  class (read_external_t):: self
  class (experiment_t), intent(out), pointer:: patch
  class(link_t), pointer:: link
  !.............................................................................
  integer :: i
  class (task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%create_patch')
  allocate (patch)

  call patch%init_unique (same=.true.) ! obtain an unique ID

  patch%task_list => self%task_list

  patch%mhd = .false.
  patch%kind = 'hlld_mhd_patch'

  call patch%gpatch_t%pre_init ()

  patch%n = self%n / self%npatches
  patch%ng = self%ng
  patch%nv = 6
  patch%nt = 2
  patch%nw = 1

  patch%box(1) = self%ds * 1
  patch%box(2) = self%ds * 1
  patch%box(3) = self%ds * self%nz

  patch%ds = self%ds

  patch%size = patch%ds * patch%n
  patch%periodic = [.true., .true., .false.]

  patch%position(1) = 0.5d0*patch%size(1)
  patch%position(2) = 0.5d0*patch%size(2)
  patch%position(3) = (i-0.5d0)*patch%size(3)

  call patch%idx%init (energy=5, mhd=patch%mhd)
  call patch%extras_t%pre_init ()
  patch%pre_initialized = .true.

  call patch%init

  patch%level  = 1 ! this is set manually

  !---------------------------------------------------------------------------
  ! Make sure interpolations are done optimally
  !---------------------------------------------------------------------------
  patch%unsigned(patch%idx%d) = .true.
  patch%unsigned(patch%idx%e) = .true.
  patch%pervolume(patch%idx%px) = .true.
  patch%pervolume(patch%idx%py) = .true.
  patch%pervolume(patch%idx%pz) = .true.

  call patch%set (bits%has_output)
  !---------------------------------------------------------------------------
  allocate (patch%filled(patch%gn(1),patch%gn(2),patch%gn(3)))
  call io%bits_mem (storage_size(patch%filled), &
                    product(shape(patch%filled)), 'fill')
  patch%filled = -1
  !---------------------------------------------------------------------------
  ! Add the patch to the task list -- this creates the patch%link
  !---------------------------------------------------------------------------
  allocate (halo_m3d_t:: link)
  link%task => patch
  patch%link => link
  call self%task_list%prepend_link (link)
  call trace%end ()
END SUBROUTINE create_patch


!===============================================================================
!> Reads the actual text files and stores variables under self
!===============================================================================
SUBROUTINE read_snapshot(self)
  !-----------------------------------------------------------------------------
  USE m3io_mod, only: NextEntry
  class(read_external_t):: self
  !-----------------------------------------------------------------------------
  integer :: lunit, i, iostat
  !-----------------------------------------------------------------------------
  ! Open atmosphere file and read MHD variables
  !-----------------------------------------------------------------------------
  open(file=self%atmos_file, newunit=lunit)
  call nextentry(lunit)
  read(lunit, *) self%atmos_name
  call nextentry(lunit)
  read(lunit, *) self%nz

  self%nz = self%nz - 2*self%ng(3)

  allocate(self%zi(1-self%ng(3):self%nz+self%ng(3)))
  allocate(self%temp, self%pe, self%rho, self%vturb, mold=self%zi)

  do i=lbound(self%zi,1), ubound(self%zi,1)
    call nextentry(lunit)
    read(lunit, *) self%zi(i), self%temp(i), self%pe(i), self%rho(i), self%vturb(i)
  end do

  self%zi  = self%zi  / scaling%l
  self%rho = self%rho / scaling%d
  self%pe  = self%pe  / scaling%p

  self%n = 1
  self%n(3) = self%nz
  self%ds = self%zi(2) - self%zi(1)

  if (self%ds /= self%zi(self%nz) - self%zi(self%nz-1)) then
    write(stdout,*) 'read_external_t%read_snapshot: Warning: Non-equidistant atmosphere grid!'
  endif
END SUBROUTINE read_snapshot

!===============================================================================
!> Takes variables read from text file and loads them into patch memory
!===============================================================================
SUBROUTINE load_authority(self, patch, ii)
  !-----------------------------------------------------------------------------
  class (read_external_t):: self
  class (experiment_t), pointer:: patch
  integer :: i,j,k, it, iz, ii

  it = 1 ! time index, it=1 is the initial condition

  do k=patch%li(3), patch%ui(3)
    ! flipping z-axis, so that the first patch is the bottom of the atmosphere
    iz = self%nz - ((ii-1)*patch%n(3) + k - patch%li(3))

    do j=patch%li(2), patch%ui(2)
      do i=patch%li(1), patch%ui(1)
        patch%mem(i,j,k,patch%idx%d ,it,1) = self%rho(iz)
        patch%mem(i,j,k,patch%idx%px,it,1) = 0
        patch%mem(i,j,k,patch%idx%py,it,1) = 0
        patch%mem(i,j,k,patch%idx%pz,it,1) = 0
        patch%mem(i,j,k,patch%idx%e ,it,1) = self%temp(iz)
        patch%mem(i,j,k,          6 ,it,1) = self%pe(iz)
        patch%filled(i,j,k) = patch%level
      end do
    end do
  end do

  if (ii == 1) then
    patch%mesh(3)%lower_boundary = .true.
    patch%boundaries%bottop(1) = .true.
    patch%boundaries%temp(1) = self%temp(ubound(self%zi,1))
    patch%boundaries%rho (1) = self%rho (ubound(self%zi,1))
    patch%boundaries%pe  (1) = self%pe  (ubound(self%zi,1))
  elseif (ii == self%npatches(3)) then
    patch%mesh(3)%upper_boundary = .true.
    patch%boundaries%bottop(2) = .true.
    patch%boundaries%temp(2) = self%temp(lbound(self%zi,1))
    patch%boundaries%rho (2) = self%rho (lbound(self%zi,1))
    patch%boundaries%pe  (2) = self%pe  (lbound(self%zi,1))
  endif

  call patch%rt%load_initial_conditions()
END SUBROUTINE

!===============================================================================
!> Download guard zones from neighbouring patches and boundary conditions
!===============================================================================
SUBROUTINE fill_guard_zones (self)
  class(read_external_t):: self
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin('read_external_t%dnload')
  call io%header ('-dnload: begin', left=25)
  call timer%tic ()
  link => self%task_list%head
  do while (associated(link))
    call link%task%dnload
    call link%task%bndries
    link => link%next
  end do
  call timer%toc ('time to dnload', 1)
  call io%header ('-dnload: end', left=25)
  call trace%end ()
END SUBROUTINE fill_guard_zones

END MODULE