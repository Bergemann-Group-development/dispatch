!> =============================================================================
!> Add the RT tasks to the task_list and queue them
!> =============================================================================
SUBROUTINE rt_setup (task_list)
  USE link_mod
  USE task_list_mod
  USE rt_mod
  USE rmhd_mod
  USE task_mod
  USE patch_mod
  USE task_list_mod
  implicit none
  type(task_list_t)      :: rt_task_list
  type(task_list_t)      :: task_list
  class(link_t), pointer :: mhd_link, rt_link, nbor_link, next
  class(task_t), pointer :: task, rt_task
  !-----------------------------------------------------------------------------
  mhd_link => task_list%head
  call rt_task_list%init
  do while (associated(mhd_link))
    !---------------------------------------------------------------------------
    ! For tasks in the task list of type rmht_t
    !---------------------------------------------------------------------------
    task => mhd_link%task
    select type (task)
    class is (rmhd_t)
      call task%dnload
      rt_task => task%rt
      !-------------------------------------------------------------------------
      ! Create a task list element for the RT task, and append to rt_task_list
      !-------------------------------------------------------------------------
      allocate (rt_link)
      rt_link%task      => rt_task
      select type (rt_task)
      class is (patch_t)
        rt_task%link => rt_link
      end select
      call rt_task_list%append_link (rt_link)
      !-------------------------------------------------------------------------
      ! Make the MHD task an nbor of the RT task
      !-------------------------------------------------------------------------
      allocate (nbor_link)
      nbor_link%task    => task
      nbor_link%link    => mhd_link
      nbor_link%needed   = .true.
      nbor_link%needs_me = .true.
      rt_link%nbor      => nbor_link
      !-------------------------------------------------------------------------
      ! Make an nbor link pointing to the RT task and append it to the nbor list
      ! of the mhd_task
      !-------------------------------------------------------------------------
      allocate (nbor_link)
      nbor_link%task      => rt_task
      nbor_link%link      => rt_link
      call mhd_link%nbor%add_nbor_by_rank (nbor_link)
      nbor_link%needed   = .true.
      nbor_link%needs_me = .true.
      !-------------------------------------------------------------------------
      ! Adding upstream and downstream links (not needed in the current setup
      ! for Marconi since we use a single patch in the vertical)
      !-------------------------------------------------------------------------
      !allocate (upstream)
      !upstream      => some_link
      !upstream%task => some_link%task
      !rt_bundle%nbor%add_nbor_by_rank (upstream)
      !upstream%needed   = .true.
      !upstream%needs_me = .false.
      !allocate (dnstream)
      !dnstream      => some_link
      !dnstream%task => some_link%task
      !rt_bundle%nbor%add_nbor_by_rank (dnstream)
      !dnstream%needed   = .false.
      !dnstream%needs_me = .true.
    end select
    mhd_link => mhd_link%next
  end do
  call task_list%print ('task_list before RT')
  call rt_task_list%print ('rt_task_list')
  rt_link => rt_task_list%head
  rt_link => rt_task_list%head
  do while (associated(rt_link))
    task => rt_link%task
    next => rt_link%next
    !---------------------------------------------------------------------------
    ! Append the new task list link to the task_list and queue it
    !---------------------------------------------------------------------------
    select type (task)
    class is (rt_t)
      call task_list%append_link (rt_link)
      rt_link%task%time = -999d0
    end select
    rt_link => next
  end do
  call task_list%print ('task_list after adding RT')
END SUBROUTINE rt_setup

PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  USE eos_mod
  USE opacity_mod
  USE scaling_mod
  type (cartesian_t):: cartesian
  type (task_list_t), pointer:: task_list
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call mpi%init
  call io%init('turbulence')            ! I/O system
  call scaling%init                     ! Scaling to code units
  call eos%init                         ! Equation-of-state
  call opacity%init                     ! Opacity
  call cartesian%init                   ! Initialize the patch list and patches
  call rt_setup (cartesian%task_list)   ! add RT tasks to task_list
  call cartesian%task_list%execute      ! execute the list = update member patches
  call mpi%end
!> =============================================================================
END PROGRAM dispatch
