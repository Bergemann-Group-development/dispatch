!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE solver_mod
  implicit none
  public
  type, public:: solver_t
    real(8):: time=0d0, dtime=1d0
    real, pointer, dimension(:,:,:,:):: uin, uout, scr
  contains
    procedure:: init
  end type
  integer, parameter:: tasks_per_block=64
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
  !$omp critical (app_cr)
  if (associated(self%tail)) then
    link%prev => self%tail
    self%tail%next => link
    self%tail => link
  else
    self%tail => link
    self%head => link
  end if
  !$omp end critical (app_cr)
  self%n = self%n + 1
END SUBROUTINE append
!===============================================================================
FUNCTION pop (self) RESULT (link)
  class(task_list_t):: self
  class(link_t), pointer:: link
  !$omp critical (pop_cr)
  link => self%head
  if (associated(link)) then
    self%head => self%head%next
    if (associated(self%head)) then
      nullify(self%head%prev)
    end if
  end if
  self%n = self%n - 1
  !$omp end critical (pop_cr)
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
  real, pointer, contiguous, dimension(:,:,:,:):: scr
  integer:: id, n_tasks, n_threads
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
  !-----------------------------------------------------------------------------
  ! Loop over task updates
  !-----------------------------------------------------------------------------
  m = 0
  used = omp_get_wtime()
  !$omp parallel private(i,j,k,l,link) shared(m,task_list)
  do while (task_list%n > 0)
    !$omp atomic update
    m = m+1
    link => task_list%pop()
    !---------------------------------------------------------------------------
    ! Distribute tasks over GPUs
    !---------------------------------------------------------------------------
    allocate (scr(1:n,1:n,1:n,1:nscr))
    do k=1,8
      do j=1,n
        do i=1,n
          !$omp simd
          do l=1,n
            scr(l,i,j,k) = link%task%uin(l,i,j,k)
          end do
        end do
      end do
    end do
    do k=2,nscr
      do j=1,n
        do i=1,n
          !$omp simd
          do l=1,n
            scr(l,i,j,k) = scr(l,i,j,k-1) + 1.0
          end do
        end do
      end do
    end do
    do k=1,8
      do j=1,n
        do i=1,n
         !$omp simd
          do l=1,n
            link%task%uout(l,i,j,k) = scr(l,i,j,nscr-nv+k)
          end do
        end do
      end do
    end do
    deallocate (scr)
    out = link%task%uout(1,1,1,1)
    link%task%time = link%task%time + link%task%dtime
    if (link%task%time < end_time) then
      call task_list%append(link)
    end if
  end do
  !$omp master
  n_threads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel
  used = omp_get_wtime() - used
  print *, 'mik13c.f90: Using CPUs, n_threads = ', n_threads
  print *, 'FINISHED: the answer is ...', link%task%uout(1,1,1,1), real(used), ' sec'
  n3 = n**3
  n3 = n3*m*nv
  print '(a,f8.1,i8,a,i12,a)', &
    ' Mz/s =', real(n3)/used*1e-6, m, ' updates', n3, ' zones'
END PROGRAM DISPATCH_mockup
