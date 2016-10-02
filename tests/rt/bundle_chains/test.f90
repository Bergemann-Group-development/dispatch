!> =============================================================================
!> This test demonstrates a module hierarchy, where the rt_*_t tasks are 
!> able to access data in the mhd_t task, while at the same time the mhd_rt_t
!> task can access data in both mhd_t and rt_*_t tasks.  Each of the modules
!> below use only the module definitions that have already occured.  If the
!> modules were put in separate files, the Makefile dependencies would not have
!> any circular dependencies. Graphically, one can draw the hierarchy as
!>
!>            task_t
!>           /  |   \
!>          /   |    link_t
!>         /    |   /   \
!>        /    patch_t   list_t
!>       /       |          |
!>      /       mhd_t       |
!>     /       /   |        |
!>    rt_data_t    |        |
!>    |       \    |        |
!> rt_bundle_t \   |        |
!>           \  \  |        |
!>           mhd_rt_t       |
!>              |           |
!>         experiment_t     |
!>                    \     |
!>                    task_list_t
!>
!> =============================================================================

!> -----------------------------------------------------------------------------
!> The task_t object holds generic task information, such as task id and time.
!> -----------------------------------------------------------------------------
MODULE task_mod
  implicit none
  private
  type, public :: task_t
    integer :: id=0
    real(8) :: time=0.0
  end type
END MODULE task_mod

!> -----------------------------------------------------------------------------
!> The link_t object allows building a linked list of tasks
!> -----------------------------------------------------------------------------
MODULE link_mod
  USE task_mod
  implicit none
  private
  type, public :: link_t
    class(link_t), pointer :: next => null()
    class(task_t), pointer :: task => null()
  end type
END MODULE link_mod

!> -----------------------------------------------------------------------------
!> The list_t object defins a list of links
!> -----------------------------------------------------------------------------
MODULE list_mod
  USE link_mod
  implicit none
  private
  type, public :: list_t
    class(link_t), pointer :: head => null()
  end type
END MODULE list_mod

!> -----------------------------------------------------------------------------
!> The patch_t object extends the task_t with mesh and memory
!> -----------------------------------------------------------------------------
MODULE patch_mod
  USE task_mod
  USE link_mod
  implicit none
  private
  type, public, extends(task_t) :: patch_t
    class (link_t), pointer :: nbor => null()
    class (task_t), pointer :: task => null()
    real, dimension(:,:,:,:,:,:), pointer :: mem => null()
  end type
END MODULE patch_mod

!> -----------------------------------------------------------------------------
!> The mhd_t object holds data and procedures for MHD
!> -----------------------------------------------------------------------------
MODULE mhd_mod
  USE patch_mod
  implicit none
  private
  type, public, extends(patch_t) :: mhd_t
  contains
    procedure :: update
    procedure :: apply_heating
  end type
CONTAINS
!> -----------------------------------------------------------------------------
!> Update the MHD time derivatives
!> -----------------------------------------------------------------------------
SUBROUTINE update (self)
  class (mhd_t) :: self
  real, dimension(:,:,:), pointer:: q
  !.............................................................................
  ! ... compute the MHD update ...
END SUBROUTINE update
!> -----------------------------------------------------------------------------
!> Apply heating
!> -----------------------------------------------------------------------------
SUBROUTINE apply_heating (self, q)
  class (mhd_t) :: self
  real, dimension(:,:,:), pointer:: q
  !.............................................................................
  self%mem(:,:,:,2,1,2) = self%mem(:,:,:,2,1,2) + q
END SUBROUTINE apply_heating
END MODULE mhd_mod

!> -----------------------------------------------------------------------------
!> The rt_data_t object type contains a pointer to the associated MHD data, and
!> adds storage for source functions, opacities, and net heating rate q.
!> -----------------------------------------------------------------------------
MODULE rt_data_mod
  USE task_mod
  USE mhd_mod
  implicit none
  private
  type, public, extends(task_t) :: rt_data_t
    class (mhd_t),          pointer :: mhd
    real, dimension(:,:,:), pointer :: source, opacity
    real, dimension(:,:,:), pointer :: q
  contains
    procedure :: update
  end type
CONTAINS
!> -----------------------------------------------------------------------------
!> The update procedure computes source function and opacity
!> -----------------------------------------------------------------------------
SUBROUTINE update (self)
  class (rt_data_t) :: self
  !.............................................................................
END SUBROUTINE update
END MODULE rt_data_mod

!> -----------------------------------------------------------------------------
!> The rt_bundle_t object is an extension of the rt_data_t object, and contains
!> pointers to related object. The next pointer points to the next bundle
!> (if any) at a boundary (real or virtual) in this patch. The downstream
!> pointer points to the downstream object.  Input and output are pointers to
!> initial and final values of the net intensities q along the bundle chain.
!> -----------------------------------------------------------------------------
MODULE rt_bundle_mod
  USE rt_data_mod
  implicit none
  private
  type, public, extends(rt_data_t) :: rt_bundle_t
    class (rt_bundle_t),  pointer :: next       => null()
    class (rt_bundle_t),  pointer :: downstream => null()
    real, dimension(:,:), pointer :: input      => null()
    real, dimension(:,:), pointer :: output     => null()
  contains
    procedure :: update
  end type
CONTAINS
!> -----------------------------------------------------------------------------
!> Updating consist of computing the local solution, pointing the input net
!> intensity if the downstream task to its output, and calling that task update
!> -----------------------------------------------------------------------------
SUBROUTINE update (self)
  class (rt_bundle_t) :: self
  !.............................................................................
  ! ... solve the RT for the bundle
  self%downstream%input => self%output
  call self%downstream%update
END SUBROUTINE update
END MODULE rt_bundle_mod

!> -----------------------------------------------------------------------------
!> The mhd_rt_t object serves to coordinate and combine the MHD and RT tasks.
!> It has a pointer to the mhd object, the rt_data object, and to a list of
!> budles to start the chains from.  Each instance has a pointer to the next
!> downstream instance.
!> -----------------------------------------------------------------------------
MODULE mhd_rt_mod
  USE mhd_mod
  USE rt_data_mod
  USE rt_bundle_mod
  implicit none
  private
  type, public, extends(mhd_t) :: mhd_rt_t
    class (mhd_t),        pointer :: mhd
    class (rt_data_t),    pointer :: rt_data
    class (rt_bundle_t),  pointer :: bundle_chain
    class (mhd_rt_t),     pointer :: downstream
  contains
    procedure :: init
    procedure :: update
  end type
CONTAINS
!> -----------------------------------------------------------------------------
!> Initialize an mhd_rt_t, setting pointers and loading initial values
!> -----------------------------------------------------------------------------
SUBROUTINE init (self)
  class (mhd_rt_t) :: self
  !.............................................................................
END SUBROUTINE init
!> -----------------------------------------------------------------------------
!> Updating an mhd_rt_t task consists of computing the necessary RT data, and
!> launching a number of linked rt_bundle tasks, which are RT tasks for all the
!> angles involved.  In patches that are NOT boundaries (real or virtual), the
!> link to rt_bundle is null, so only data are computed.  In boundary patches,
!> the linked list of bundles consist of only the bundles that start out on the
!> boundary -- not all of the bundles the belong to the patch.
!> -----------------------------------------------------------------------------
SUBROUTINE update (self)
  class (mhd_rt_t) :: self
  class (rt_bundle_t), pointer :: bundle
  !.............................................................................
  call self%rt_data%update
  bundle => self%bundle_chain
  do while (associated(bundle))
    !$omp task firstprivate(bundle)
    call bundle%update
    !$omp end task
    bundle => bundle%next
  end do
  !$omp taskwait
  !----------------------------------------------------------------------------
  ! At this point all the bundle chains have executed, and hence the rt_data%q
  ! has the net heating.  For different solvers, applying this differs in 
  ! details, so the appply_heating procedure should be part of the mhd_t task.
  ! The net heating is therefore supplied as an (optional) argument.
  !----------------------------------------------------------------------------
  call self%mhd_t%update
  call self%mhd_t%apply_heating (self%rt_data%q)
END SUBROUTINE update
END MODULE mhd_rt_mod

!> -----------------------------------------------------------------------------
!> The main purpose of the experiment_t object is to present a generic name to
!> the task_list procedures.
!> -----------------------------------------------------------------------------
MODULE experiment_mod
  USE mhd_rt_mod
  implicit none
  private
  type, public, extends(mhd_rt_t) :: experiment_t
  end type
END MODULE experiment_mod

!> -----------------------------------------------------------------------------
!> In the simplest form, the task_list_t object just contains a pointer to the
!> first element in a list of links that contain pointes to tasks.
!> -----------------------------------------------------------------------------
MODULE task_list_mod
  USE list_mod
  implicit none
  private
  type, public, extends(list_t) :: task_list_t
  contains
    procedure :: execute
  end type
CONTAINS
!> -----------------------------------------------------------------------------
!> -----------------------------------------------------------------------------
SUBROUTINE execute (self)
  class (task_list_t) :: self
  !.............................................................................
END SUBROUTINE execute
END MODULE task_list_mod

!> -----------------------------------------------------------------------------
!> MOdule for setting up the RT part of the data and procedures
!> -----------------------------------------------------------------------------
MODULE rt_mod
  USE link_mod
  USE task_list_mod
  USE experiment_mod
  USE rt_bundle_mod
  implicit none
  private
  type, public :: rt_t
  contains
    procedure :: init
  end type
CONTAINS
!> -----------------------------------------------------------------------------
!> Run through the task list.  We have access at this point to the links to the
!> MHD tasks, and their nbor lists.  It is assumed that this is used to define
!> the downstream links in the rt_internal tasks.
!> -----------------------------------------------------------------------------
SUBROUTINE init (self, task_list)
  class (rt_t):: self
  type (task_list_t), pointer :: task_list
  type (link_t),       pointer :: link
  !.............................................................................
  link => task_list%head
  do while (associated(link))
    ! ... use the link to figure out downstream bundle connections ...
    associate (task=>link%task)
    select type (task)
    class is (experiment_t)
      call task%init
    end select
    end associate
    link => link%next
  end do
END SUBROUTINE init
END MODULE rt_mod

!> -----------------------------------------------------------------------------
!> Setting up the experiment consists of first setting up the MHD experiment,
!> and then just initializing the extra rt_*_t objects
!> -----------------------------------------------------------------------------
PROGRAM test
  USE link_mod
  USE rt_mod
  USE task_list_mod
  USE experiment_mod
  implicit none
  type (task_list_t), pointer :: task_list
  type (rt_t) :: rt
  !.............................................................................
  ! ... use some standard MHD setup to generate a task list ...
  call rt%init (task_list)
  call task_list%execute
END PROGRAM
