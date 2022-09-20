!===============================================================================
!> Top level MHD task definition -- can either rely entirely on lower levels,
!> including extras_t%initial_conditions, or can optionally do such things here
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE units_mod, only: cgs
  USE scaling_mod, only: scaling
  USE kinds_mod
  USE m3d_patch_mod, only: hash_string
  USE m3d_halo_mod, only: m3d_halo_t
  USE shared_mod
  USE halo_mod
  USE list_mod
  USE link_mod
  USE task_mod
  USE math_mod
  USE timer_mod
  USE os_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    integer :: iter=0, maxiter, zs, ze
    character(len=2) :: v5
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: info
    procedure:: task_info
    procedure:: pre_update
    procedure:: update
    procedure:: post_update
    procedure:: finish
    procedure:: init_level
    procedure:: bndries
    procedure:: append_to
    procedure:: statistics
    procedure:: hash_key
    procedure:: init_nbors_hash
    procedure:: check_nbor
    procedure:: nbor_relations => experiment_nbor_relations
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

  self%time     = 0.0
  self%dtime    = timestep
  self%out_next = huge(1.0)
  self%place    = -1
  !call self%set(bits%remove)

  call self%solver_t%pre_init ()
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
  select type (link => self%link)
  class is (m3d_halo_t)
    call self%solver_t%init
  class default
    call self%task_t%init
  end select
  call self%init_level()
  call trace%end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%bndries', itimer=itimer)
  mem => self%mem(:,:,:,:,self%it,:)
  if (self%boundaries%lower .or. self%boundaries%upper) call self%boundaries%boundary(mem)
  call trace%end (itimer)
END SUBROUTINE bndries

!===============================================================================
SUBROUTINE pre_update (self)
  class(experiment_t), target:: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call self%task_t%pre_update()
  call self%set(bits%no_download)
END SUBROUTINE pre_update

!===============================================================================
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  !call self%solver_t%update
  call self%post_update()
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
SUBROUTINE post_update (self)
  class(experiment_t):: self
  integer :: p(3)
  !----------------------------------------------------------------------------
  if (self%verbose > 1) then
      p = [1,1,1]
      p = p + self%li - 1
      write(*,'(F5.1, A12, I4, 12(ES15.7))') self%time, trim(self%kind), self%id, &
      self%mem(p(1),p(2),p(3),:,1,1)
  endif

  if (self%iter == 0) then
    call self%finish()
  endif
  
  self%iter = self%iter + 1
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
    'experiment_t%task_info: kind, id, rank, level, mu, nu, BV, pos', &
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
  if (any(self%n > 1)) then
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  end if
END SUBROUTINE init_level

!===============================================================================
SUBROUTINE append_to(self, list)
  class (experiment_t), target :: self
  class(list_t),  pointer :: list
  class(m3d_halo_t),  pointer :: link
  !-----------------------------------------------------------------------------
  allocate (link)
  call link%init()
  link%task => self
  self%link => link
  call link%append_to(list)
  self%key = self%hash_key()
  call list%hash_table%add (self%key, self%link)
END SUBROUTINE append_to

!===============================================================================
SUBROUTINE experiment_nbor_relations(self, another, needed, needs_me, download)
  class(experiment_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%experiment_nbor_relations')
  call self%extras_t%nbor_relations(another, needed, needs_me, download)
  call trace%end()
END SUBROUTINE experiment_nbor_relations

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
