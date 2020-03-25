!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE solver_mod
  implicit none
  private
  type, public:: solver_t
    real(8):: time=0d0, dtime=1d0
    real, pointer, dimension(:,:,:,:):: uin, uout, scr
  contains
    procedure:: init
    procedure:: update
  end type
  integer:: n=20, nv=8, nscr=50
CONTAINS
!===============================================================================
SUBROUTINE init (self)
  class(solver_t):: self
  allocate (self%uin(n,n,n,nv), self%uout(n,n,n,nv))
  self%uin = 0d0
END SUBROUTINE init
!===============================================================================
SUBROUTINE update (self)
  class(solver_t):: self
  integer:: i
  allocate (self%scr(n,n,n,nscr))
  self%scr(:,:,:,1:nv) = self%uin
  do i=2,nscr
    self%scr(:,:,:,i) = self%scr(:,:,:,i-1) + 1.0
  end do
  self%uout = self%scr(:,:,:,nscr-nv+1:nscr)
  self%time = self%time + self%dtime
END SUBROUTINE update
END MODULE solver_mod
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
  nullify(self%head%prev)
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
  integer:: i, ntask=30
  real(8):: end_time=100d0
  do i=1,ntask
    allocate (task, link)
    link%task => task
    call task%init
    call task_list%append (link)
  end do
  !$omp parallel private(link)
  do while (task_list%n > 0)
    link => task_list%pop ()
    call link%task%update
    if (link%task%time < end_time) call task_list%append (link)
  end do
  !$omp end parallel
  print *, 'FINISHED: the answer is ...', link%task%uout(1,1,1,1)
END PROGRAM DISPATCH_mockup