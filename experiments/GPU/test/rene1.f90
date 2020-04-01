!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE solver_mod
  implicit none
  public
  type, public:: solver_t
    real(8):: time=0d0, dtime=1d0
    real, pointer, dimension(:,:,:,:):: uin, uout, scr
  contains
    procedure:: init
    procedure:: update
  end type
  integer:: n=48, nv=8, nscr=50
CONTAINS
!===============================================================================
SUBROUTINE init (self)
  class(solver_t):: self
  allocate (self%uin(n,n,n,nv), self%uout(n,n,n,nv))
  self%uin = 1d0
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
  deallocate (self%scr)
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
  class(solver_t), pointer:: task
  class(link_t), pointer:: link
  type(task_list_t):: task_list
  integer:: i, j, k, l, ntask=256, m
  integer(8):: n3
  real(8):: ntime=0d0, ndtime=1d0
  real(8):: end_time=10d0
  real(8):: used
  real, pointer, dimension(:,:,:,:):: uin, uout, scr
  do i=1,ntask
    allocate (task, link)
    link%task => task
    call task%init
    call task_list%append (link)
  end do
  m = 0
  used = omp_get_wtime()
  do while (task_list%n > 0)
    link => task_list%pop ()
      ! unpack solver task
      uin  => link%task%uin
      uout => link%task%uout
      ntime = link%task%time
      ndtime = link%task%dtime
      allocate (scr(n,n,n,nscr))
      !$OMP TARGET DATA map(alloc:scr)
      do k=1,nv
        scr(:,:,:,k) = uin(:,:,:,k)
      end do
      !$OMP TARGET map(from:uout)
        !$omp teams distribute num_teams(64) thread_limit(2)
        do j=1,n
          do k=2,nscr
            !$omp parallel do
            do i=1,n
              !$omp simd
              do l=1,n
                scr(l,i,j,k) = scr(l,i,j,k-1) + 1.0
              end do
            end do
          end do
        enddo
        !$omp end teams distribute
        !$omp teams distribute parallel do collapse(3)
        do k=1,nv
          do j=1,n
            do i=1,n
              !$omp simd
              do l=1,n
                uout(l,i,j,k) = scr(l,i,j,nscr-nv+k)
              end do
            end do
          end do
        end do
        !$omp end teams distribute parallel do
      !$OMP END TARGET
      !$OMP END TARGET DATA
      m = m+1
      deallocate (scr)
      ntime = ntime + ndtime
      link%task%time = ntime
      !if (mod(m,100)==0) &
        !print *, m, link%task%time
    if (link%task%time < end_time) call task_list%append (link)
  end do
  used = omp_get_wtime() - used
  print *, 'rene1.f90: Using thread_limit(2) with omp parallel, and collapse(3)'
  print *, 'FINISHED: the answer is ...', link%task%uout(1,1,1,1), real(used), ' sec'
  n3 = n**3
  n3 = n3*nv*m
  print '(a,f8.1,i8,a,i12,a)', &
    ' Mz/s =', real(n3)/used*1e-6, m, ' updates', n3, ' zones'
END PROGRAM DISPATCH_mockup
