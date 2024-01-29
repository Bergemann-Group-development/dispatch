!===============================================================================
!> Top level MHD task definition -- can either rely entirely on lower levels,
!> including extras_t%initial_conditions, or can optionally do such things here
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE os_mod
  USE data_io_mod
  USE bits_mod
  USE trace_mod
  USE shared_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t, link_t
  USE list_mod, only: list_t
  USE halo_mod, only: halo_t
  USE solver_mod, only: solver_t
  USE units_mod, only: cgs
  USE scaling_mod, only: scaling
  USE handyfuncs_mod
  USE m3d_halo_mod, only: m3d_halo_t
  USE m3d_patch_mod, only: verbose, print_execution_history, use_density, use_ne, &
  long_vec, pivot_z, rotate_atmos, atmos_format
  USE composition_patch_mod, only: composition_patch_t
  USE composition_mod, only: composition_type
  USE read_stagger_snapshot_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    class(composition_patch_t), pointer :: composition => null()
    real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
    integer :: iter=0, maxiter, zs, ze
    integer :: dims(3), nbox(3)
    character(len=2) :: v5
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: info
    procedure:: task_info
    procedure:: pre_update
    procedure:: update
    procedure:: post_update
    procedure:: rotate
    procedure:: finish
    procedure:: init_level
    procedure:: bndries
    procedure:: append_to
    procedure:: statistics
    procedure:: hash_key
    procedure:: init_nbors_hash
    procedure:: check_nbor
    procedure:: nbor_relations
  end type

  logical, public :: interpolated_3d_atmos=.false.
  
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self)
  class(experiment_t), target :: self
  integer, save:: itimer=0
  real(8), save:: timestep=1.d0
  integer :: md(3), thread
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%pre_init', itimer=itimer)

  self%periodic = [.true., .true., .false.]

  self%n = self%nbox / self%dims
  md = mod(self%nbox, self%dims)

  self%zs = self%ipos(3) * self%n(3) + min(md(3), self%ipos(3))
  if (self%ipos(3)+1 <= md(3)) then
    self%n(3) = self%n(3) + 1
  endif
  self%ze = self%zs + self%n(3)

  self%size = self%ds * self%n
  self%position = 0.5d0*self%size
  self%position(3) = self%position(3) + self%zs*self%ds(3)

  !-----------------------------------------------------------------------------
  ! Preliminary region information, relative to the task position.
  ! The settings below allow overlap() calls to work
  !-----------------------------------------------------------------------------
  self%authority%l = - self%size*0.5_8
  self%authority%u = + self%size*0.5_8
  self%interest%l  = self%authority%l - self%ng*self%size/self%n
  self%interest%u  = self%authority%u + self%ng*self%size/self%n

  allocate (self%boundaries)
  if (self%ipos(3)==0)              call self%boundaries%set(bits%zl)
  if (self%ipos(3)==self%dims(3)-1) call self%boundaries%set(bits%zu)
  call self%boundaries%pre_init (self)

  self%time     = 0.0
  self%dtime    = timestep
  self%out_next = huge(1.0)

  !self%mhd = .false.
  !call self%solver_t%pre_init ()
  call self%idx%init (energy=5, mhd=.false.)

  call self%append_to(self%task_list)

  allocate(self%composition)
  self%composition%patch => self
  self%composition%task_list => self%task_list
  self%composition%zs = self%zs
  call self%composition%pre_init()

  self%nv = 6
  self%nt = 4 ! Download returns wrong stuff if nt <= 2
  self%nw = 1 + size(long_vec, 2)

  call trace%end (itimer=itimer)
END SUBROUTINE pre_init

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  class(composition_type), pointer:: composition
  real(4), dimension(:,:,:,:), pointer, contiguous :: authority
  real(4), dimension(:,:,:), pointer, contiguous :: d, px, py, pz, tt, pe
  real(4), dimension(:,:,:), allocatable :: ne, vel
  real(8), dimension(:,:,:), allocatable :: ne8, rho8, temp8
  real(4), dimension(:,:,:,:), pointer, contiguous :: chi
  real(8) :: nu(1), sub(3), pos(3), bottom, bot2pivot
  integer :: i, n(3), li(3), ui(3), thread, loaded=0
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%init', itimer=itimer)
  call self%patch_t%init()

  li = self%li
  ui = self%ui
  n  = self%n

  self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)

  if (interpolated_3d_atmos) then
    sub = scaling%l * self%size
    pos = scaling%l * self%position

    authority => self%authority_mem(:,:,:,:,self%it,1)
    d  => authority(:,:,:,self%idx%d)
    px => authority(:,:,:,self%idx%px)
    py => authority(:,:,:,self%idx%py)
    pz => authority(:,:,:,self%idx%pz)
    tt => authority(:,:,:,self%idx%e)
    pe => authority(:,:,:,6)

    call load_atmos_subregion(n, sub, pos, d, px, py, pz, tt, pe)

    if (rotate_atmos) then
      do i=1, self%nw-1
        authority => self%authority_mem(:,:,:,:,self%it,1+i)
        d  => authority(:,:,:,self%idx%d)
        px => authority(:,:,:,self%idx%px)
        py => authority(:,:,:,self%idx%py)
        pz => authority(:,:,:,self%idx%pz)
        tt => authority(:,:,:,self%idx%e)
        pe => authority(:,:,:,6)
        bottom = self%position(3) - 0.5*self%size(3)
        bot2pivot = scaling%l * (pivot_z * self%ds(3) - bottom)
        call load_atmos_subregion(n, sub, pos, d, px, py, pz, tt, pe, bot2pivot, long_vec(:,i))
      end do
    endif

    loaded = loaded + 1
    if (loaded == self%dims(3)) then
      call deallocate_stagger()
    endif
  endif

  call os%allocate(self%filled, self%gn)
  self%filled = -1

  if (associated(self%boundaries)) then
    call self%init_bdries ()
  endif
  call data_io%init (self)
  !call self%init_level()
  call trace%end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%bndries', itimer=itimer)
  call trace%end (itimer)
END SUBROUTINE bndries

!===============================================================================
SUBROUTINE pre_update (self)
  class(experiment_t), target:: self
  integer, save:: itimer=0
  logical :: ready, busy
  !----------------------------------------------------------------------------
  call self%task_t%pre_update()
  call self%set(bits%no_download)
  call print_execution_history(self, 'PRE')
END SUBROUTINE pre_update

!===============================================================================
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  call self%post_update()
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
SUBROUTINE rotate (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call self%task_t%rotate
  if (self%iter == 0) then
    call self%finish()
  endif
  self%iter = self%iter + 1
END SUBROUTINE rotate

!===============================================================================
SUBROUTINE post_update (self)
  class(experiment_t):: self
  integer :: p(3)
  !----------------------------------------------------------------------------
  !if (self%verbose > 1) then
  !    p = [1,1,1]
  !    p = p + self%li - 1
  !    write(*,'(F5.1, A12, I4, 12(ES15.7))') self%time, trim(self%kind), self%id, &
  !    self%mem(p(1),p(2),p(3),:,1,1)
  !endif
  call print_execution_history(self, 'END')
END SUBROUTINE post_update

!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(experiment_t)    :: self
  integer            :: nq, ntask
  character(len=64)  :: experiment_name
  optional           :: nq, ntask, experiment_name
  !-----------------------------------------------------------------------------
  !print *, 'Experiment update called for patch', self%id
END SUBROUTINE info

!===============================================================================
!> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
!===============================================================================
SUBROUTINE task_info (self, unit)
  class(experiment_t):: self
  integer, optional:: unit
  integer:: unit_L
  !-----------------------------------------------------------------------------
  unit_L = shared%option (unit, io_unit%mpi)
  write(unit_L,'(1x,a,2x,a10,i6,4i4,7x,2l1,2x,3f9.3)') &
    'experiment_t%task_info: kind, id, rank, level, nu, mu, BV, pos', &
    self%kind, self%id, self%rank, self%level, 0, 0, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
END SUBROUTINE task_info

!===============================================================================
SUBROUTINE finish (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  !self%finished = .true.
  self%time = io%end_time
  self%t = io%end_time
  call self%set(bits%frozen)
END SUBROUTINE finish

!===============================================================================
SUBROUTINE init_level (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  !if (irregular) then
  !  self%level = stagger_levels(self%ipos(3)+1)
  !elseif (any(self%n > 1)) then
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  !end if
END SUBROUTINE init_level

!===============================================================================
SUBROUTINE append_to(self, list)
  class(experiment_t), target :: self
  class(list_t),  pointer :: list
  class(m3d_halo_t),  pointer :: link
  !-----------------------------------------------------------------------------
  allocate (link)
  call link%init()
  link%task => self
  self%link => link
  call self%init_level()
  !call link%append_to(list)
  call list%append_link(self%link)
  self%key = self%hash_key()
  call list%hash_table%add (self%key, self%link)
END SUBROUTINE append_to

!!===============================================================================
!SUBROUTINE experiment_nbor_relations(self, another, needed, needs_me, download)
!  class(experiment_t) :: self
!  class(task_t), pointer :: another         ! another task in the list
!  logical :: needed, needs_me, download    ! Is the other task: needed, 
!                                           !    does it need me, should I download?
!  !-----------------------------------------------------------------------------
!  call trace%begin ('experiment_t%experiment_nbor_relations')
!  call self%extras_t%nbor_relations(another, needed, needs_me, download)
!  call trace%end()
!END SUBROUTINE experiment_nbor_relations

!===============================================================================
SUBROUTINE statistics (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  ! do nothing!
END SUBROUTINE statistics

!===============================================================================
!> Return a unique hash key for the patch
!===============================================================================
FUNCTION hash_key (self, position, size, level, box, origin) RESULT (key)
  class(experiment_t):: self
  real(8), dimension(3), optional:: position, origin, size, box
  integer, optional:: level
  integer:: key(4)
  !-----------------------------------------------------------------------------
  key(1) = self%ipos(3) + 1 
  key(2) = 0
  key(3) = 0
  key(4) = hash_string(self%kind)
END FUNCTION hash_key
  
!===============================================================================
!> Set default nbor relations, and allow underlings to modify them.  This is
!> called from list_t%init_nbors and from list_t%reset_status()
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(experiment_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%nbor_relations')
  if (another%kind == self%kind .or. another%kind == 'composition') then
    call self%task_t%nbor_relations(another, needed, needs_me, download) ! keep
    select type (another)
    class is (composition_patch_t)
      if (another%patch%id == self%id) then
        needed=.false.; download=.false.
      else
        needed=.false.; needs_me=.false.; download=.false.
      endif
    class default
      needed=.false.; needs_me=.false.; download=.false.
    end select
  else
    needed=.false.; needs_me=.false.; download=.false.
  endif
  call trace%end()
END SUBROUTINE nbor_relations

!===============================================================================
!> Fast method to generate nbor lists in cases with very long task lists
!===============================================================================
SUBROUTINE init_nbors_hash (self, nbors)
  class(experiment_t):: self
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  integer :: i, key(4)
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init_nbors_hash', itimer=itimer)
  nullify (nbors)
  !$omp atomic write
  self%n_nbors = 0
  key(1:3) = self%key(1:3)
  key(4) = hash_string('composition')
  call self%check_nbor(key, nbors)
  call trace%end (itimer)
END SUBROUTINE init_nbors_hash

!===============================================================================
SUBROUTINE check_nbor (self, key, nbors)
  class(experiment_t):: self
  class(nbor_t), pointer:: nbors
  class(link_t), pointer:: nblink
  integer:: key(4)
  logical:: added
  !-----------------------------------------------------------------------------
  call self%task_list%hash_table%get (key, nblink)
  if (associated(nblink)) then
    select type (link => self%link)
    class is (halo_t)
      call link%conditionally_add_nbor (nblink, nbors, added)
      if (added) then
        !$omp atomic update
        self%n_nbors = self%n_nbors + 1
      end if
    end select
  end if
END SUBROUTINE check_nbor

END MODULE
