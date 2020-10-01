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
  class(task_list_t):: self
  class(link_t), pointer:: link
  !$omp critical (head_cr)
  link => self%head
  self%head => self%head%next

  if (associated(self%head)) then
   if (associated(self%head%prev)) then
    nullify(self%head%prev)
   end if
  end if

  self%n = self%n - 1
  !$omp end critical (head_cr)
END FUNCTION pop
END MODULE task_list_mod
!===============================================================================
PROGRAM DISPATCH_mockup
  USE solver_mod
  USE omp_lib
  USE task_list_mod
  class(solver_t), pointer:: task
  class(link_t), pointer:: link
  type(task_list_t):: task_list
  integer:: i
  !-----------------------------------------------------------------------------
  ! Main parameters
  !-----------------------------------------------------------------------------
  integer:: n=NSIZE
  integer:: ntask=30
  integer(kind=omp_lock_kind):: LCK
  real(8):: end_time=100d0
  real(8):: time
  !-----------------------------------------------------------------------------
  do i=1,ntask
    allocate (task, link)
    link%task => task
    task%n = n
    call task%init
    call task_list%append (link)
  end do

  time = omp_get_wtime()

  call omp_init_lock(lck)
  
  !$omp parallel shared (lck, task_list, end_time) private(link) default (none)
  
  call omp_set_lock(lck)

  do while (task_list%n > 0)
    link => task_list%pop ()
    call omp_unset_lock(lck)
    
    call link%task%update(end_time)
    ! print *, 'step: ', omp_get_thread_num(),  link%task%time
    
    call omp_set_lock(lck)
    !if (link%task%time < end_time) call task_list%append (link)
  end do
  call omp_unset_lock(lck)

  !$omp barrier
  
  print *, 'FINISHED: the answer is ...', link%task%mem(1,1,1,1,1,1)
  
  !$omp end parallel
  
  time = omp_get_wtime() - time
 
  print *, 'FINISHED: the answer is ...', link%task%mem(1,1,1,1,1,1), time
END PROGRAM DISPATCH_mockup
