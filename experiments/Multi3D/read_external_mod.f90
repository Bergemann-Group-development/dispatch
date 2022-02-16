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
  USE m3d_patch_mod

  type, public:: read_external_t
    class(task_list_t), pointer:: task_list
    integer :: npatches(3), n(3), ng(3), verbose
    real(4), dimension(:,:,:), allocatable :: temp, ne, pe, rho, vturb, px, py, pz
    real(4), dimension(:,:,:,:), allocatable :: nh
    real(8), dimension(:), allocatable :: xi, yi, zi
    real(8) :: ds(3)
    character(len=126) :: atmos_name
  contains
    procedure:: init
    procedure:: read_1d_snapshot
    procedure:: read_stagger_snapshot
    procedure:: create_patch
    procedure:: load_authority
  end type read_external_t

  type(read_external_t), public:: read_external

CONTAINS

SUBROUTINE init(self, task_list, name)
  !-----------------------------------------------------------------------------
  class(read_external_t):: self
  class(task_list_t),  pointer:: task_list
  character(len=*),   optional:: name
  !-----------------------------------------------------------------------------
  class(experiment_t), pointer:: patch
  logical, save:: first_time=.true.
  logical, save :: read_Stagger=.false.
  integer :: iostat, iz, i
  integer, save :: npatches=1, ng(3)=1, verbose=1
  character(len=126) :: atmos_file='./input_multi3d/atmos/atmos.sun_MARCS_50'
  namelist /atmos_params/ atmos_file, npatches, ng, verbose, read_Stagger
  !-----------------------------------------------------------------------------
  if (associated(task_list)) return
  call trace%begin('read_external_t%init')
  call io%header ('read_external_t%init')
  !-----------------------------------------------------------------------------
  ! The first thread reads the namelist parameters and saves them
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, nml=atmos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, atmos_params)
    end if
    call io%namelist_warning ('atmos_params', iostat)
    write (io%namelist, atmos_params)

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
  
  if (read_Stagger) then
    call self%read_stagger_snapshot(atmos_file)
  else
    call self%read_1d_snapshot(atmos_file)
  endif

  if (abs(self%ds(3) - (self%zi(self%n(3)) - self%zi(self%n(3)-1))) > 1d-4) then
    write(stdout,*) 'read_external_t%init: Warning: Non-equidistant atmosphere grid!'
  endif

  if (mod(self%n(3), npatches) /= 0) then
    write(stdout,*) 'read_external_t%init: Error: npatches must be a factor of the atmosphere size!'
    write(stdout,*) 'npatches = ', self%npatches(3), 'nz = ', self%n(3)
    write(stdout,*) 'Stopping.'
    stop
  endif

  do i=1, npatches
    call self%create_patch(i, patch) ! patch is allocated here

    if (verbose > 0) then 
      iz = (i-1)*self%n(3)/npatches + 1
      if (i==1) write(stdout,'(I2, A28)') npatches, ' subdomains in z-direction:'
      write(stdout,'(I4, I6,A3,I4)') i, iz, ' - ', iz + self%n(3)/npatches - 1
    endif
  end do

  io%end_time = (patch%maxiter + patch%ilambd) * patch%dtime
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

  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize a minimum amount of patch information in a Cartesian arrangement
!===============================================================================
SUBROUTINE create_patch (self, i, patch)
  class (read_external_t):: self
  class (experiment_t), intent(out), pointer:: patch
  integer :: i
  !-----------------------------------------------------------------------------
  call trace%begin ('read_external_t%create_patch')
  allocate (patch)
  patch%task_list => self%task_list

  call patch%init_unique (same=.true.) ! obtain an unique ID

  patch%mhd = .false.
  patch%kind = 'hlld_mhd_patch'

  call patch%pre_init()
  patch%pre_initialized = .true.

  patch%n = self%n / self%npatches
  patch%ng = self%ng
  patch%nv = 6
  patch%nt = 4 ! Download returns wrong stuff if nt <= 2
  patch%nw = 1

  patch%box = self%ds * self%n
  patch%ds  = self%ds

  patch%size = patch%ds * patch%n
  patch%periodic = [.true., .true., .false.]

  patch%position(1) = 0.5d0*patch%size(1)
  patch%position(2) = 0.5d0*patch%size(2)
  patch%position(3) = (i - 0.5d0)*patch%size(3)
  patch%origin = patch%position - 0.5d0*patch%size

  call patch%init()
  call patch%idx%init (energy=5, mhd=patch%mhd)

  call self%load_authority(patch, i)
  
  !---------------------------------------------------------------------------
  ! Add the patch to the task list -- this creates the patch%link
  !---------------------------------------------------------------------------
  call patch%append_to(patch%task_list)

  call patch%init_sub_patches()

  call trace%end ()
END SUBROUTINE create_patch

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
END SUBROUTINE read_1d_snapshot

!===============================================================================
SUBROUTINE read_stagger_snapshot(self, stagger_file)
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
  !allocate(self%nh(self%n(1), self%n(2), self%n(3), 6))
  allocate(self%rho, self%ne, self%px, self%py, self%pz, self%pe, mold=self%temp)

  open(file=trim(stagger_file) // '.atmos', newunit=lunit, form='UNFORMATTED', access='stream')
  read(lunit) self%ne 
  read(lunit) self%temp
  read(lunit) self%px
  read(lunit) self%py
  read(lunit) self%pz
  read(lunit) self%rho
  !read(lunit) self%nh
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
END SUBROUTINE read_stagger_snapshot


!===============================================================================
!> Takes variables read from text file and loads them into patch memory
!===============================================================================
SUBROUTINE load_authority(self, patch, ii)
  !-----------------------------------------------------------------------------
  class (read_external_t):: self
  class (experiment_t), pointer:: patch
  integer :: i,j,k, it, iz, ii, xs, xe, ys, ye, zs, ze

  it = 1 ! time index, it=1 is the initial condition
  xs = patch%li(1); xe = patch%ui(1)
  ys = patch%li(2); ye = patch%ui(2)
  zs = patch%li(3); ze = patch%ui(3)

  do k=zs,ze
    ! flipping z-axis, so that the first patch is the bottom of the atmosphere
    iz = self%n(3) - ((ii-1)*patch%n(3) + k - zs)

    patch%mem(xs:xe,ys:ye,k,patch%idx%d ,it,1) = self%rho(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%px,it,1) = self%px(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%py,it,1) = self%py(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%pz,it,1) = self%pz(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx%e ,it,1) = self%temp(:,:,iz)
    patch%mem(xs:xe,ys:ye,k,patch%idx_pe,it,1) = self%pe(:,:,iz)

    if (allocated(self%nh)) then
      patch%nh(:,:,k-zs+1,:) = self%nh(:,:,iz,:)
    endif
  end do

  if (allocated(self%nh)) patch%nh_missing = .false.

  patch%rho  = dble(patch%mem(xs:xe,ys:ye,zs:ze,patch%idx%d , it, 1)) * scaling%d
  patch%temp = dble(patch%mem(xs:xe,ys:ye,zs:ze,patch%idx%e , it, 1))
  patch%ne   = dble(patch%mem(xs:xe,ys:ye,zs:ze,patch%idx_pe, it, 1)) * scaling%p / (cgs%k_b * patch%temp)
  patch%vx   = dble(patch%mem(xs:xe,ys:ye,zs:ze,patch%idx%px, it, 1)) * scaling%d * scaling%u / patch%rho
  patch%vy   = dble(patch%mem(xs:xe,ys:ye,zs:ze,patch%idx%py, it, 1)) * scaling%d * scaling%u / patch%rho
  patch%vz   = dble(patch%mem(xs:xe,ys:ye,zs:ze,patch%idx%pz, it, 1)) * scaling%d * scaling%u / patch%rho

  ! Top of the atmosphere is bottom of index space!!!
  !    -> Top of atmosphere is top in index space (swiched if conditions)
  if (ii == self%npatches(3)) then
    patch%mesh(3)%upper_boundary = .true.
    patch%boundaries%topbot(1) = .true.
  endif

  if (ii == 1) then
    patch%mesh(3)%lower_boundary = .true.
    patch%boundaries%topbot(2) = .true.
  endif

  allocate(patch%filled(patch%gn(1), patch%gn(2), patch%gn(3)))
  patch%filled = -1

END SUBROUTINE

END MODULE