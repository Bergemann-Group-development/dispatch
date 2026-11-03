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
  USE gpatch_mod, only: gpatch_t
  USE handyfuncs_mod
  USE m3d_halo_mod, only: m3d_halo_t, add_nbor_by_key
  USE m3d_patch_mod, only: verbose, print_execution_history, amr
  USE rt_params_mod, only: rt_out
  USE composition_patch_mod, only: composition_patch_t
  USE read_atmos_mod, only: ia
  implicit none
  private
  type, public, extends(gpatch_t):: experiment_t
    class(composition_patch_t), pointer :: composition => null()
    integer :: iter=0, zs, ze
    integer :: dims(3), nbox(3)
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

CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
  SUBROUTINE pre_init (self)
    class(experiment_t), target :: self
    integer, save:: itimer=0
    real(8), save:: timestep=1.d0
    !-----------------------------------------------------------------------------
    call trace%begin('experiment_t%pre_init', itimer=itimer)

    self%periodic = [.true., .true., .false.]

    self%n = self%nbox / self%dims
    call unequal_intervals(self%nbox(3), self%ipos(3)+1, self%dims(3), self%zs, self%ze)
    self%n(3) = self%ze - self%zs + 1

    self%size(1:2) = self%box(1:2)
    self%ds(1:2) = self%size(1:2) / self%n(1:2)
    self%position(1:2) = 0.5d0*self%size(1:2)

    if (amr=='') then
      self%ds(3)       = self%box(3) / self%nbox(3)
      self%size(3)     = self%ds(3) * self%n(3)
      self%position(3) = 0.5*self%size(3) + (self%zs - 1)*self%ds(3)
    else
      self%size(3)     = ia%sub_size(self%ipos(3)+1)
      self%position(3) = ia%sub_pos (self%ipos(3)+1)
      self%ds(3)       = self%size(3) / self%n(3)
    endif

    !-----------------------------------------------------------------------------
    ! Preliminary region information, relative to the task position.
    ! The settings below allow overlap() calls to work
    !-----------------------------------------------------------------------------
    self%authority%l = - self%size*0.5_8
    self%authority%u = + self%size*0.5_8
    self%interest%l  = self%authority%l - self%ng*self%ds
    self%interest%u  = self%authority%u + self%ng*self%ds

    allocate (self%boundaries)
    if (self%ipos(3)==0)              call self%boundaries%set(bits%zl)
    if (self%ipos(3)==self%dims(3)-1) call self%boundaries%set(bits%zu)
    call self%boundaries%pre_init (self)

    self%time     = 0.0
    self%dtime    = timestep
    self%out_next = huge(1.0)

    call self%append_to(self%task_list)

    allocate(self%composition)
    self%composition%patch => self
    self%composition%task_list => self%task_list
    self%composition%zs = self%zs
    call self%composition%pre_init()

    ! no mem: this patch carries only geometry and the task-list lifecycle; the
    ! atmosphere is interpolated straight into its composition patch
    self%nv = 0
    self%nt = 4 ! Download returns wrong stuff if nt <= 2
    self%nw = 1 + size(rt_out%vec, 2)

    ! being memory-less (nv=0), this patch has no density slot: mark it so the
    ! task list's pre-update validate%check does not index a null/absent mem
    call self%set (bits%no_density)
    call self%set (bits%no_filled)

    call trace%end (itimer=itimer)
  END SUBROUTINE pre_init


!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
  SUBROUTINE init (self)
    class(experiment_t):: self
    integer, save:: itimer=0
    !----------------------------------------------------------------------------
    call trace%begin('experiment_t%init', itimer=itimer)
    call self%patch_t%init()

    call os%allocate(self%filled, self%gn)
    self%filled = -1

    if (associated(self%boundaries)) then
      call self%init_bdries ()
    endif
    call data_io%init (self)
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
    !----------------------------------------------------------------------------
    call print_execution_history(self, 'END')
  END SUBROUTINE post_update

  !===============================================================================
  SUBROUTINE info (self, nq, ntask, experiment_name)
    class(experiment_t)    :: self
    integer            :: nq, ntask
    character(len=64)  :: experiment_name
    optional           :: nq, ntask, experiment_name
    !-----------------------------------------------------------------------------
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
    self%time = io%end_time
    self%t = io%end_time
    call self%set(bits%frozen)
  END SUBROUTINE finish

  !===============================================================================
  SUBROUTINE init_level (self)
    class(experiment_t):: self
    !-----------------------------------------------------------------------------
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
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
    call list%append_link(self%link)
    self%key = self%hash_key()
    call list%hash_table%add (self%key, self%link)
  END SUBROUTINE append_to

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
    integer :: key(4)
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
    integer:: key(4)
    !-----------------------------------------------------------------------------
    call add_nbor_by_key (self, key, nbors)
  END SUBROUTINE check_nbor

END MODULE
