!> =============================================================================
!> This test demonstrates a module hierarchy, where the rt_*_t tasks are 
!> able to access data in the mhd_t task, while at the same time the mhd_rt_t
!> task can access data in both mhd_t and rt_*_t tasks.
!> =============================================================================
MODULE task_mod
  implicit none
  private
  type, public:: task_t
    integer:: id
  end type
END MODULE task_mod

MODULE link_mod
  USE task_mod
  implicit none
  private
  type, public:: link_t
    type(link_t), pointer:: next => null()
    type(task_t), pointer:: task => null()
  end type
END MODULE link_mod

MODULE patch_mod
  USE task_mod
  USE link_mod
  implicit none
  private
  type, public, extends(task_t):: patch_t
    type(link_t), pointer:: nbor => null()
    type(task_t), pointer:: task => null()
  end type
END MODULE patch_mod

MODULE mhd_mod
  USE patch_mod
  implicit none
  private
  type, public, extends(patch_t):: mhd_t
  end type
END MODULE mhd_mod

MODULE rt_mod
  USE task_mod
  implicit none
  private
  type, public, extends(task_t):: rt_t
  end type
END MODULE rt_mod

MODULE rt_data_mod
  USE rt_mod
  USE mhd_mod
  implicit none
  private
  type, public, extends(rt_t):: rt_data_t
    type(mhd_t),         pointer:: mhd
    type(rt_t),          pointer:: rt
  end type
END MODULE rt_data_mod

MODULE rt_internal_mod
  USE rt_data_mod
  implicit none
  private
  type, public, extends(rt_data_t):: rt_internal_t
    type(rt_data_t),     pointer:: rt_data
  end type
END MODULE rt_internal_mod

MODULE mhd_rt_mod
  USE mhd_mod
  USE rt_mod
  USE rt_data_mod
  USE rt_internal_mod
  implicit none
  private
  type, public, extends(mhd_t):: mhd_rt_t
    type(mhd_t),         pointer:: mhd
    type(rt_t),          pointer:: rt
    type(rt_data_t),     pointer:: rt_data
    type(rt_internal_t), pointer:: rt_internal
  end type
END MODULE mhd_rt_mod

PROGRAM test
  USE mhd_rt_mod
  implicit none
END
