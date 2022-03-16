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
  USE m3d_patch_mod, only: halo_m3d_t
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
    procedure:: update
    procedure:: post_update
    procedure:: finish
    procedure:: init_level
    procedure:: bndries
    procedure:: append_to
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
  call trace%begin('experiment_mod::pre_init', itimer=itimer)

  self%time     = 0.0
  self%dtime    = timestep
  self%out_next = huge(1.0)

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
  call trace%begin('experiment_mod::init', itimer=itimer)
  call self%solver_t%init
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
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: it_l, in_l, in_lte_l
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
  integer :: pcheck(3)
  !----------------------------------------------------------------------------
  if (self%verbose > 1) then
      pcheck = [1,1,1]
      write(*,'(F5.1, A12, I4, 12(ES15.7))') self%time, trim(self%kind), self%id, &
      self%mem(self%li(1)+pcheck(1)-1,self%li(2)+pcheck(2)-1,self%li(3)+pcheck(3)-1,:,1,1)
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
  print *, 'Experiment update called for patch', self%id
END SUBROUTINE info

!===============================================================================
SUBROUTINE finish (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  self%finished = .true.
  self%time = io%end_time
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
SUBROUTINE append_to(self, task_list)
  class (experiment_t), target :: self
  class(list_t),  pointer :: task_list
  class(link_t),  pointer :: link
  class (task_t), pointer :: task
  !-----------------------------------------------------------------------------
  allocate (halo_m3d_t:: link)
  call link%init()
  link%task => self
  self%link => link
  call task_list%append (link)
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

END MODULE
