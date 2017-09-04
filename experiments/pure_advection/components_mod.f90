!> =============================================================================
!> Set up the pulse advection benchmark
!> =============================================================================
MODULE components_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE link_mod
  USE patch_mod
  USE experiment_mod
  USE cartesian_mod
  implicit none
  private

  type, public:: components_t
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure:: init
  end type

CONTAINS

!> =============================================================================
!> Initialisation
!> =============================================================================
SUBROUTINE init (self)
  class(components_t):: self
  type(cartesian_t):: cartesian
  class(link_t), pointer:: head, new_link
  class(patch_t), pointer:: patch, floating_patch
  !----------------------------------------------------------------------------
  call trace_begin('components_t%init')

  ! initialise base domain
  call cartesian%init('pulse_advection')
  self%task_list => cartesian%task_list
  head => self%task_list%head
  patch => task2patch(head%task)

  ! initialise a single patch, with the same resolution as the base domain patches,
  ! but centred and encompassing the initial pulse.
  select type (patch)
  type is (experiment_t)
    if (patch%add_patch) then
      allocate (experiment_t:: floating_patch)
      floating_patch%size = patch%psize * (patch%pbufferfrac + 1.0)
      floating_patch%position = patch%pcentre
      floating_patch%llc_cart = floating_patch%position - 0.5 * floating_patch%size
      floating_patch%llc_nat = floating_patch%llc_cart !FIXME
      floating_patch%quality = patch%quality + 1.0
      call floating_patch%init()
      floating_patch%mesh(:)%origin = patch%mesh(:)%origin

      ! move the patch at the bulk velocity
      if (patch%pmoving) floating_patch%velocity = [patch%pvx, patch%pvy, patch%pvz]
      
      ! clear status bits
      call floating_patch%clear(bits%root_grid)
      !call floating_patch%clear(bits%static)

      ! now add the patch to the task list
      allocate (new_link)
      new_link%task => floating_patch
      new_link%parent => head
      call self%task_list%append_link (new_link)
      call self%task_list%queue_by_time (new_link)
    end if
  end select

  ! ensure the neighbour list is correct
  call self%task_list%init_all_nbors

  call trace_end
END SUBROUTINE init

END MODULE components_mod

