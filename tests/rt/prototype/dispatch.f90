!> =============================================================================
!> Add the RT tasks to the task_list and queue them
!> =============================================================================
SUBROUTINE rt_setup (task_list)
  USE link_mod
  USE task_list_mod
  USE rt_mod
  USE rmhd_mod
  USE task_mod
  USE task_list_mod
  implicit none
  type (task_list_t)      :: task_list
  class (link_t), pointer :: mhd_link, rt_link
  class (task_t), pointer :: task
  type (task_list_t)      :: rt_task_list
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
      !-------------------------------------------------------------------------
      ! Create a new task list element for the RT task
      !-------------------------------------------------------------------------
      allocate (rt_link)
      rt_link%task => task%rt
      !-------------------------------------------------------------------------
      ! Make the MHD task an nbor of the RT task, and vice versa
      !-------------------------------------------------------------------------
      rt_link%nbor => mhd_link
      call mhd_link%nbor%add_nbor_by_rank (rt_link)
      !-------------------------------------------------------------------------
      ! Append to rt_task_list
      !-------------------------------------------------------------------------
      call rt_task_list%append_link (rt_link)
      print *, 'added RT task', task%rt%id, rt_link%task%id
    end select
    mhd_link => mhd_link%next
  end do
  rt_link => rt_task_list%head
  do while (associated(rt_link))
    task => rt_link%task
    print *,task%id
    !---------------------------------------------------------------------------
    ! Append the new task list link to the task_list and queue it
    !---------------------------------------------------------------------------
    select type (task)
    class is (rt_t)
      call task_list%append_link (rt_link)
      rt_link%task%time = 0d0
      call task_list%queue_by_time (rt_link)
      print *, 'added RT task', task%id
    end select
    rt_link => rt_link%next
  end do
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
