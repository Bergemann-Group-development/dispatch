MODULE p_arrays_mod
  USE omp_lib
  implicit none
  private
  type, public:: p_arrays_t
    integer:: n
    real:: cs
    real:: pc
    real:: vg(3)
    real:: pa(3)
    integer(2), allocatable:: q(:,:) ! 4*16 =  64
    real      , allocatable:: w(:)   ! 1*32 =  32
    real      , allocatable:: r(:,:) ! 3*32 =  96
    real      , allocatable:: v(:,:) ! 3*32 =  96
  contains                           ! total: 288
    procedure:: init
    procedure:: update
    procedure:: dealloc
    procedure:: test
  end type
CONTAINS

SUBROUTINE init (self, n)
  class(p_arrays_t):: self
  integer:: n
  !-----------------------------------------------------------------------------
  self%n = n
  allocate (self%w(n), self%q(n,4), self%r(n,3), self%v(n,3))
END SUBROUTINE init

SUBROUTINE dealloc (self)
  class(p_arrays_t):: self
  !-----------------------------------------------------------------------------
  deallocate (self%w, self%q, self%r, self%v)
END SUBROUTINE dealloc

SUBROUTINE update (self, dt)
  class(p_arrays_t):: self
  real:: dt
  integer:: j, q(3)
  real:: pa(3), dv(3), pd(3), pc
  !-----------------------------------------------------------------------------
  !$omp simd
  do j=1,self%n
    !- kick 1 ------------------------------------------------------------
    dv(1) = self%vg(1)-self%v(j,1)
    dv(2) = self%vg(2)-self%v(j,2)
    dv(3) = self%vg(3)-self%v(j,3)
    pc = self%pc*sqrt(1.0 &
       + (dv(1)**2+dv(2)**2+dv(3)**2)/self%cs**2)
    pd(1) = pc*dv(1)
    pd(2) = pc*dv(2)
    pd(3) = pc*dv(3)
    self%v(j,1) = self%v(j,1) + dt*(self%pa(1)+pd(1))/(1.+pc*dt)
    self%v(j,2) = self%v(j,2) + dt*(self%pa(2)+pd(2))/(1.+pc*dt)
    self%v(j,3) = self%v(j,3) + dt*(self%pa(3)+pd(3))/(1.+pc*dt)
    !- drift -------------------------------------------------------------
    self%r(j,1) = self%r(j,1) + 2.*dt*self%v(j,1)
    self%r(j,2) = self%r(j,2) + 2.*dt*self%v(j,2)
    self%r(j,3) = self%r(j,3) + 2.*dt*self%v(j,3)
    pd(1) = pc*(self%vg(1)-self%v(j,1))
    pd(2) = pc*(self%vg(2)-self%v(j,2))
    pd(3) = pc*(self%vg(3)-self%v(j,3))
    !- kick 2 ------------------------------------------------------------
    self%v(j,1) = self%v(j,1) + dt*(self%pa(1)+pd(1))/(1.+pc*dt)
    self%v(j,2) = self%v(j,2) + dt*(self%pa(2)+pd(2))/(1.+pc*dt)
    self%v(j,3) = self%v(j,3) + dt*(self%pa(3)+pd(3))/(1.+pc*dt)
    !- update q+r representation -----------------------------------------
    q(1) = floor(self%r(j,1))
    q(2) = floor(self%r(j,2))
    q(3) = floor(self%r(j,3))
    self%r(j,1) = self%r(j,1) - q(1)
    self%r(j,2) = self%r(j,2) - q(2)
    self%r(j,3) = self%r(j,3) - q(3)
    self%q(j,1) = self%q(j,1) + q(1)
    self%q(j,2) = self%q(j,2) + q(2)
    self%q(j,3) = self%q(j,3) + q(3)
  end do
END SUBROUTINE update

SUBROUTINE test (self, np, n, niter)
  class(p_arrays_t):: self
  integer:: n, np, niter
  type(p_arrays_t), allocatable:: pp(:,:,:), pp1(:,:,:)
  real:: dt
  real(8):: start, used
  integer:: i1, i2, i3, j, iter, q(3)
  !-----------------------------------------------------------------------------
  dt = 1.0
  do iter=1,niter
    !---------------------------------------------------------------------------
    start = omp_get_wtime()
    allocate (pp(n,n,n))
    do i3=1,n
    do i2=1,n
    do i1=1,n
      call pp(i1,i2,i3)%init(np)
      pp(i1,i2,i3)%q(:,1) = i1
      pp(i1,i2,i3)%q(:,2) = i2
      pp(i1,i2,i3)%q(:,3) = i3
      pp(i1,i2,i3)%pc    = 1.
      pp(i1,i2,i3)%cs    = 1.
      pp(i1,i2,i3)%pa(:) = 1.
      pp(i1,i2,i3)%vg(:) = 1.
    end do
    end do
    end do
    print *, iter, np, n, (omp_get_wtime()-start)*1e9/(n**3*np), 'init'
    !---------------------------------------------------------------------------
    start = omp_get_wtime()
    allocate (pp1(n,n,n))
    do i3=1,n
    do i2=1,n
    do i1=1,n
      do j=1,pp(i1,i2,i3)%n
        q = pp(i1,i2,i3)%q(j,1:3)
        pp1(q(1),q(2),q(3)) = pp(i1,i2,i3)
      end do
    end do
    end do
    end do
    print *, iter, np, n, (omp_get_wtime()-start)*1e9/(n**3*np), 'copy'
    !---------------------------------------------------------------------------
    start = omp_get_wtime()
    do i3=1,n
    do i2=1,n
    do i1=1,n
      do j=1,pp(i1,i2,i3)%n
        q = pp(i1,i2,i3)%q(j,1:3)
        pp1(q(1),q(2),q(3))%r = pp(i1,i2,i3)%r
        pp1(q(1),q(2),q(3))%v = pp(i1,i2,i3)%v
        pp1(q(1),q(2),q(3))%w = pp(i1,i2,i3)%w
        pp1(q(1),q(2),q(3))%vg = pp(i1,i2,i3)%vg
      end do
    end do
    end do
    end do
    print *, iter, np, n, (omp_get_wtime()-start)*1e9/(n**3*np), 'copy'
    !---------------------------------------------------------------------------
    start = omp_get_wtime()
    do i3=1,n
    do i2=1,n
    do i1=1,n
      call pp1(i1,i2,i3)%update (dt)
    end do
    end do
    end do
    print *, iter, np, n, (omp_get_wtime()-start)*1e9/(n**3*np), 'update'
    deallocate (pp, pp1)
  end do
END SUBROUTINE test

END MODULE p_arrays_mod
