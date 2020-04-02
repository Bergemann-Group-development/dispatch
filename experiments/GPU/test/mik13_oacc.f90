!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE solver_mod
  implicit none
  public
  type, public:: solver_t
    real(8):: time=0d0, dtime=1d0
    real, pointer, contiguous, dimension(:,:,:,:):: uin, uout, scr
  contains
    procedure:: init
  end type
  integer, parameter:: n_threads=1, tasks_per_block=64
  integer:: n=48, nv=8, nscr=50, ntask=256
  real(8):: end_time=10d0
CONTAINS
!===============================================================================
SUBROUTINE init (self)
  class(solver_t):: self
  allocate (self%uin(n,n,n,nv), self%uout(n,n,n,nv))
  self%uin = 1d0
END SUBROUTINE init
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
  type, public:: links_t
    class(link_t), pointer:: link
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
  if (associated(self%tail)) then
    link%prev => self%tail
    self%tail%next => link
    self%tail => link
  else
    self%tail => link
    self%head => link
  end if
  self%n = self%n + 1
END SUBROUTINE append
!===============================================================================
FUNCTION pop (self) RESULT (link)
  class(task_list_t):: self
  class(link_t), pointer:: link
  link => self%head
  if (associated(link)) then
    self%head => self%head%next
    if (associated(self%head)) then
      nullify(self%head%prev)
    end if
  end if
  self%n = self%n - 1
END FUNCTION pop
END MODULE task_list_mod
!===============================================================================
PROGRAM DISPATCH_mockup
  USE omp_lib
  USE solver_mod
  USE task_list_mod
  implicit none
  class(solver_t), pointer:: task
  class(link_t), pointer:: link
  type(task_list_t):: task_list
  integer:: i, j, k, l, m
  real(8):: ntime=0d0, ndtime=1d0
  real(8):: used, t0, t1, time(3)
  real:: out
  real, pointer, contiguous, dimension(:,:,:,:):: scr, p
  real, pointer, contiguous, dimension(:,:,:,:,:):: uin, uout
  integer:: id, n_tasks
  integer(8):: n3
  type(links_t), pointer:: links(:)
  !-----------------------------------------------------------------------------
  ! Initialize tasks and the task_list
  !-----------------------------------------------------------------------------
  do i=1,ntask
    allocate (task, link)
    link%task => task
    call task%init
    call task_list%append (link)
  end do
  allocate (links(tasks_per_block))
  !-----------------------------------------------------------------------------
  ! Loop over task updates
  !-----------------------------------------------------------------------------
  m = 0
  allocate (uin (n,n,n,nv,tasks_per_block))
  time = 0d0
  used = omp_get_wtime()
  do while (task_list%n > 0)
    t0 = omp_get_wtime()
    !---------------------------------------------------------------------------
    ! Prepare a block of tasks
    !---------------------------------------------------------------------------
    do i=1,tasks_per_block
      link => task_list%pop()
      if (associated(link)) then
        links(i)%link => link
      else
        exit
      end if
    end do
    n_tasks = i-1
    do i=1,n_tasks
      uin(:,:,:,:,i) = links(i)%link%task%uin(:,:,:,:)
      !print *,i, uin(1,1,1,1,i)
    end do
    !print *, m, n_tasks
    !---------------------------------------------------------------------------
    ! Distribute tasks over GPUs, with e.g. 64 tasks per block, leaving 12
    ! threads per block, and overloading the 6 SMs with about 10 tasks
    !---------------------------------------------------------------------------
    t1 = omp_get_wtime()
    time(1) = time(1) + (t1-t0)
    t0 = t1
    !test = 11.
    allocate (scr(1:n,1:n,1:n,1:nscr))
    !$acc data create(scr) copy(n,nv,nscr,uin) copyout(out)
    !$acc parallel loop gang    
    do id=1,n_tasks
      !$acc loop collapse(3)  
      do k=1,8
        do j=1,n
          do i=1,n
            !$acc loop vector  
            do l=1,n
              scr(l,i,j,k) = uin(l,i,j,k,id)
            end do
          end do
        end do
      end do
     
    !$acc loop collapse(2)  
      do j=1,n
        do i=1,n
          do k=2,nscr
            !$acc loop vector  
            do l=1,n
              scr(l,i,j,k) = scr(l,i,j,k-1) + 1.0
            end do
          end do
        end do
      end do
     !$acc  loop collapse(3)   
      do k=1,8
        do j=1,n
          do i=1,n
            !$acc loop vector  
            do l=1,n
              uin(l,i,j,k,id) = scr(l,i,j,nscr-nv+k)
            end do
          end do
        end do
      end do
      out = uin(1,1,1,1,1)
    end do
    !$acc end data
    deallocate (scr)
    m = m + n_tasks
    t1 = omp_get_wtime()
    time(2) = time(2) + (t1-t0)
    t0 = t1
    !---------------------------------------------------------------------------
    ! Add tasks back onto the task list
    !---------------------------------------------------------------------------
    do id=1,n_tasks
      link => links(id)%link
      link%task%time = link%task%time + link%task%dtime
      link%task%uout = uin(:,:,:,:,id)
      !print *, maxval(link%task%uout)
      if (link%task%time < end_time) call task_list%append (link)
    end do
    t1 = omp_get_wtime()
    time(3) = time(3) + (t1-t0)
  end do
  used = omp_get_wtime() - used
  print *, 'mik13.f90: Using task blocks and collapse(2-3)'
  print *, 'FINISHED: the answer is ...', link%task%uout(1,1,1,1), real(used), ' sec'
  n3 = n**3
  n3 = n3*m*nv
  print '(a,f8.1,i8,a,i12,a)', &
    ' Mz/s =', real(n3)/used*1e-6, m, ' updates', n3, ' zones'
END PROGRAM DISPATCH_mockup
