!===============================================================================
!> Depth patches of the atmosphere: geometry and task-list placeholders that own
!> the composition patches.
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
  !===============================================================================
  !> One depth slab of the atmosphere, spanning the whole horizontal box; it holds
  !> no data and owns the slab's composition patch.
  !===============================================================================
  type, public, extends(gpatch_t):: experiment_t
    class(composition_patch_t), pointer :: composition => null()
    integer :: iter=0  !< number of rotations
    integer :: maxiter !< unused
    integer :: zs      !< first depth point of the slab in the global mesh
    integer :: ze      !< last depth point of the slab in the global mesh
    integer :: dims(3) !< number of patches along each axis
    integer :: nbox(3) !< number of mesh points of the whole box
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
  !> Sets the slab's geometry: periodic in x and y, the whole horizontal box, and
  !> its share zs:ze of the nbox(3) depth points (or, with amr set, the sub-domain
  !> of read_atmos). Then adds the patch to the task list and creates its
  !> composition patch. Holding no data (nv=0), the patch is marked no_density and
  !> no_filled.
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
  !> Initialises the patch, its filled map, its boundaries and its data I/O.
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
  !> Does nothing: the patch holds no data to set boundary values on.
  !===============================================================================
  SUBROUTINE bndries (self)
    class(experiment_t):: self
    integer, save:: itimer=0
    !----------------------------------------------------------------------------
    call trace%begin('experiment_t%bndries', itimer=itimer)
    call trace%end (itimer)
  END SUBROUTINE bndries

  !===============================================================================
  !> Prepares the update like any task, without downloading neighbour data, and
  !> logs it in the execution history.
  !===============================================================================
  SUBROUTINE pre_update (self)
    class(experiment_t), target:: self
    !----------------------------------------------------------------------------
    call self%task_t%pre_update()
    call self%set(bits%no_download)
    call print_execution_history(self, 'PRE')
  END SUBROUTINE pre_update

  !===============================================================================
  !> Only calls post_update; the composition, NLTE and RT patches do the work.
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
  !> Counts the rotation; the first one finishes the patch.
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
  !> Logs the end of the update in the execution history.
  !===============================================================================
  SUBROUTINE post_update (self)
    class(experiment_t):: self
    !----------------------------------------------------------------------------
    call print_execution_history(self, 'END')
  END SUBROUTINE post_update

  !===============================================================================
  !> Does nothing.
  !===============================================================================
  SUBROUTINE info (self, nq, ntask, experiment_name)
    class(experiment_t)    :: self
    integer            :: nq              !< unused
    integer            :: ntask           !< unused
    character(len=64)  :: experiment_name !< unused
    optional           :: nq, ntask, experiment_name
    !-----------------------------------------------------------------------------
  END SUBROUTINE info

  !===============================================================================
  !> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
  !===============================================================================
  SUBROUTINE task_info (self, unit)
    class(experiment_t):: self
    integer, optional:: unit !< output unit, default io_unit%mpi
    integer:: unit_L
    !-----------------------------------------------------------------------------
    unit_L = shared%option (unit, io_unit%mpi)
    write(unit_L,'(1x,a,2x,a10,i6,4i4,7x,2l1,2x,3f9.3)') &
      'experiment_t%task_info: kind, id, rank, level, nu, mu, BV, pos', &
      self%kind, self%id, self%rank, self%level, 0, 0, &
      self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
  END SUBROUTINE task_info

  !===============================================================================
  !> Freezes the patch at the end time, so the dispatcher stops scheduling it.
  !===============================================================================
  SUBROUTINE finish (self)
    class(experiment_t):: self
    !-----------------------------------------------------------------------------
    self%time = io%end_time
    self%t = io%end_time
    call self%set(bits%frozen)
  END SUBROUTINE finish

  !===============================================================================
  !> Sets level 1 on every patch.
  !===============================================================================
  SUBROUTINE init_level (self)
    class(experiment_t):: self
    !-----------------------------------------------------------------------------
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  END SUBROUTINE init_level

  !===============================================================================
  !> Appends the patch to list through an m3d_halo_t link and adds it to the list's
  !> hash table.
  !===============================================================================
  SUBROUTINE append_to(self, list)
    class(experiment_t), target :: self
    class(list_t),  pointer :: list !< task list
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
  !> Does nothing.
  !===============================================================================
  SUBROUTINE statistics (self)
    class(experiment_t):: self
    !-----------------------------------------------------------------------------
    ! do nothing!
  END SUBROUTINE statistics

  !===============================================================================
  !> The unique hash key of the patch: its slab number and the hash of its kind.
  !===============================================================================
  FUNCTION hash_key (self, position, size, level, box, origin) RESULT (key)
    class(experiment_t):: self
    real(8), dimension(3), optional:: position !< unused
    real(8), dimension(3), optional:: origin   !< unused
    real(8), dimension(3), optional:: size     !< unused
    real(8), dimension(3), optional:: box      !< unused
    integer, optional:: level                  !< unused
    integer:: key(4)
    !-----------------------------------------------------------------------------
    key(1) = self%ipos(3) + 1
    key(2) = 0
    key(3) = 0
    key(4) = hash_string(self%kind)
  END FUNCTION hash_key

  !===============================================================================
  !> Sets the neighbour relations: the patch needs no neighbour, and only its own
  !> composition patch keeps the default needs_me. Called from list_t%init_nbors and
  !> list_t%reset_status.
  !===============================================================================
  SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
    class(experiment_t):: self
    class(task_t), pointer:: another !< the other task
    logical:: needed                 !< self needs another, set
    logical:: needs_me               !< another needs self, set
    logical:: download               !< self downloads data from another, set
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
  !> Fast method to generate nbor lists in cases with very long task lists: the
  !> patch's only neighbour is its composition patch, found by its hash key.
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(experiment_t):: self
    class(nbor_t), pointer:: nbors !< neighbour list, set
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
  !> Adds the task with the hash key to nbors.
  !===============================================================================
  SUBROUTINE check_nbor (self, key, nbors)
    class(experiment_t):: self
    class(nbor_t), pointer:: nbors !< neighbour list
    integer:: key(4)               !< hash key of the task to add
    !-----------------------------------------------------------------------------
    call add_nbor_by_key (self, key, nbors)
  END SUBROUTINE check_nbor

END MODULE
