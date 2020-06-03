!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE task_list_mod
  USE solver_mod
  implicit none
  private
  type, public:: link_t
    class(solver_t), pointer:: task=>null()
    class(link_t), pointer:: next=>null(), prev=>null()
  end type
  type, public:: task_list_t
    integer:: n=0
    class(link_t), pointer:: head=>null(), tail=>null()
  contains
    procedure:: append
    procedure:: pop
  end type
CONTAINS
!===============================================================================
SUBROUTINE append (self, link)
  !!$omp declare target
  class(task_list_t):: self
  class(link_t), pointer:: link
  !$omp critical (tail_cr)
  if (associated(self%tail)) then
    link%prev => self%tail
    self%tail%next => link
    self%tail => link
  else
    self%tail => link
    self%head => link
  end if
  self%n = self%n + 1
  !$omp end critical (tail_cr)
END SUBROUTINE append
!===============================================================================
FUNCTION pop (self) RESULT (link)
  !!$omp declare target
  class(task_list_t):: self
  class(link_t), pointer:: link
  !$omp critical (head_cr)
  link => self%head
  self%head => self%head%next
  if (associated(self%head)) nullify(self%head%prev)
  self%n = self%n - 1
  !$omp end critical (head_cr)
END FUNCTION pop
END MODULE task_list_mod
!===============================================================================
PROGRAM DISPATCH_mockup
  USE solver_mod
  USE task_list_mod
  class(solver_t), pointer:: task
  class(link_t), pointer:: link
  type(task_list_t):: task_list
  integer:: i
#ifdef __OPENMP__
  double precision, external :: omp_get_wtime
#endif
  real(8):: used
  !-----------------------------------------------------------------------------
  ! Main parameters
  !-----------------------------------------------------------------------------
  integer:: n=16
  integer:: ntask=30
  real(8):: end_time=10d0
  !-----------------------------------------------------------------------------
  do i=1,ntask
    allocate (task, link)
    link%task => task
    task%n = n
    call task%init
    call task_list%append (link)
  end do
  used = omp_get_wtime()
  !$omp single
  do while (task_list%n > 0)
    link => task_list%pop ()
    call link%task%update
    if (link%task%time < end_time) call task_list%append (link)
  end do
  !$omp end single
  used = omp_get_wtime() - used
  task => link%task
  print "('FINISHED: the answer is ...',f10.6,9x,'used:',f7.3,'s')", &
    task%mem(5,5,5,1,task%it,1), used
#ifndef __OPENMP__
contains
  real(8) function omp_get_wtime()
  integer:: count, count_rate
  call system_clock (count, count_rate)
  omp_get_wtime = real(count)/real(count_rate)
  end function
#endif
END PROGRAM DISPATCH_mockup
