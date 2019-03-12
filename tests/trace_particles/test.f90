!===============================================================================
!> Program to test functionality and speed of particles, when using a list of
!> structures, or temporary arrays.
!>
!> The final, optimal code is to use the vectorized loops as long as n remains
!> the same, to do a direct loop, which caching the positions and velocities
!> if/when n decreases, and to deallocate / allocate if the number of points is
!> larger than the size of the cache arrays.
!===============================================================================
MODULE particles_mod
  USE dll_mod
  USE omp_timer_mod
  integer, parameter:: nt=4
  type, public:: particle_t
    real(8):: r(3,nt)=0.0, v(3,nt)=0.0, w=1.0
  end type
  type, public, extends(dll_t):: particles_t
    integer:: it=1
    real(8), allocatable:: r(:,:,:), v(:,:,:), w(:)
    real(8):: wt(nt)=1.0
  contains
    procedure:: update
    procedure:: interpolate
  end type
CONTAINS

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self, type)
  class(particles_t):: self
  integer:: type
  integer:: new, o0, o1, o2, o3, i ,j
  class(dll_node_t), pointer:: p
  class(*), pointer:: car
  real(8), pointer:: r(:,:), v(:,:)
  !-----------------------------------------------------------------------------
  new = 1+mod(self%it,nt)
  o0 = self%it
  o1 = 1+mod(o0+nt-2,nt)
  o2 = 1+mod(o0+nt-3,nt)
  o3 = 1+mod(o0+nt-4,nt)
  select case (type)
  !-----------------------------------------------------------------------------
  ! Implicit loops -- slow
  !-----------------------------------------------------------------------------
  case(1)
    self%r(:,:,new) = &
    self%r(:,:,o0) + &
    self%v(:,:,o0)*self%wt(1) + &
    self%v(:,:,o1)*self%wt(2) + &
    self%v(:,:,o2)*self%wt(3)
  !-----------------------------------------------------------------------------
  ! Explicit loop over j,i -- this is the fastest option
  !-----------------------------------------------------------------------------
  case(2)
    do j=1,3
    do i=1,self%n
      self%r(i,j,new) = &
      self%r(i,j,o0) + &
      self%v(i,j,o0)*self%wt(1) + &
      self%v(i,j,o1)*self%wt(2) + &
      self%v(i,j,o2)*self%wt(3)
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Implicit loop over i -- not the fastest
  !-----------------------------------------------------------------------------
  case(3)
    do j=1,3
      self%r(:,j,new) = &
      self%r(:,j,o0) + &
      self%v(:,j,o0)*self%wt(1) + &
      self%v(:,j,o1)*self%wt(2) + &
      self%v(:,j,o2)*self%wt(3)
    end do
  !-----------------------------------------------------------------------------
  ! Direct, non-vectorized update
  !-----------------------------------------------------------------------------
  case(4)
    p => self%head
    do while (associated(p))
      car => p%car
      select type (car)
      type is (particle_t)
        r=>car%r
        v=>car%v
        r(:,new) = &
        r(:,o0) + &
        v(:,o0)*self%wt(1) + &
        v(:,o1)*self%wt(2) + &
        v(:,o2)*self%wt(3)
      end select
      p => p%next
    end do
  end select
END SUBROUTINE update

!===============================================================================
!> Tri-linear interpolate in a velocity array
!===============================================================================
SUBROUTINE interpolate (self, v, type)
  class(particles_t):: self
  integer:: type
  real:: v(:,:,:,:)
  real, allocatable:: f(:,:), p(:,:)
  integer:: sz(4)
  integer:: ir(3)
  real:: r(3)
  !-----------------------------------------------------------------------------
  sz = size(v)
  select case (type)
  !-----------------------------------------------------------------------------
  ! Straight loops
  !-----------------------------------------------------------------------------
  case (1)
    do i=1,self%n
      r = 1.0 + self%r(i,:,self%it)
      ir = floor(r)
      ir = max(1,min(ir,sz(1:3)))
      r = r-ir
      do j=1,3
        self%v(i,j,self%it) = (1.-r(3))*((1.-r(2))*(v(ir(1),ir(2)  ,ir(3)  ,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)  ,ir(3)  ,j)*r(1))  &
                                            +r(2) *(v(ir(1),ir(2)+1,ir(3)  ,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)+1,ir(3)  ,j)*r(1))) &
                            +     r(3) *((1.-r(2))*(v(ir(1),ir(2)  ,ir(3)  ,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)  ,ir(3)  ,j)*r(1))  &
                                            +r(2) *(v(ir(1),ir(2)+1,ir(3)+1,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)+1,ir(3)+1,j)*r(1)))
      end do
    end do
  !-----------------------------------------------------------------------------
  ! Swap loops
  !-----------------------------------------------------------------------------
  case (2)
    do j=1,3
      do i=1,self%n
        r = 1.0 + self%r(i,:,self%it)
        ir = floor(r)
        ir = max(1,min(ir,sz(1:3)))
        r = r-ir
        self%v(i,j,self%it) = (1.-r(3))*((1.-r(2))*(v(ir(1),ir(2)  ,ir(3)  ,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)  ,ir(3)  ,j)*r(1))  &
                                            +r(2) *(v(ir(1),ir(2)+1,ir(3)  ,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)+1,ir(3)  ,j)*r(1))) &
                            +     r(3) *((1.-r(2))*(v(ir(1),ir(2)  ,ir(3)  ,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)  ,ir(3)  ,j)*r(1))  &
                                            +r(2) *(v(ir(1),ir(2)+1,ir(3)+1,j)*(1.0-r(1)) + v(ir(1)+1,ir(2)+1,ir(3)+1,j)*r(1)))
      end do
    end do
  !-----------------------------------------------------------------------------
  ! Cache corner values and use them in a vectorized loop
  !-----------------------------------------------------------------------------
  case (3)
    allocate (p(self%n,3), f(self%n,8))
    do j=1,3
    do i=1,self%n
      r = 1.0 + self%r(i,:,self%it)
      ir = floor(r)
      ir = max(1,min(ir,sz(1:3)))
      p(i,j) = r(j)-ir(j)
      f(i,1) = v(ir(1)  ,ir(2)  ,ir(3)  ,j) 
      f(i,2) = v(ir(1)+1,ir(2)  ,ir(3)  ,j) 
      f(i,3) = v(ir(1)  ,ir(2)+1,ir(3)  ,j) 
      f(i,4) = v(ir(1)+1,ir(2)+1,ir(3)  ,j) 
      f(i,5) = v(ir(1)  ,ir(2)  ,ir(3)+1,j) 
      f(i,6) = v(ir(1)+1,ir(2)  ,ir(3)+1,j) 
      f(i,7) = v(ir(1)  ,ir(2)+1,ir(3)+1,j) 
      f(i,8) = v(ir(1)+1,ir(2)+1,ir(3)+1,j) 
    end do
    do i=1,self%n
      self%v(i,j,self%it) = (1.-p(i,3))*((1.-p(i,2))*(f(i,1)*(1.0-p(i,1)) + f(i,2)*p(i,1))  &
                                            +p(i,2) *(f(i,3)*(1.0-p(i,1)) + f(i,4)*p(i,1))) &
                          +     p(i,3) *((1.-p(i,2))*(f(i,5)*(1.0-p(i,1)) + f(i,6)*p(i,1))  &
                                            +p(i,2) *(f(i,7)*(1.0-p(i,1)) + f(i,8)*p(i,1)))
    end do
    end do
    deallocate (f, p)
  end select
END SUBROUTINE interpolate

END MODULE particles_mod

!===============================================================================
PROGRAM test
  USE random_mod
  USE particles_mod
  integer:: n=1000000, m=16, type
  class(*), pointer:: particle
  class(dll_node_t), pointer:: node
  type(particles_t):: particles
  real(8):: used
  real, allocatable:: v(:,:,:,:)
  type(random_t):: random
  !-----------------------------------------------------------------------------
  call particles%init
  particles%n = n
  used = wallclock()
  allocate (particles%r(n,3,nt))
  allocate (particles%v(n,3,nt))
  do i=1,n
    allocate (node)
    allocate (particle_t:: particle)
    node%car => particle
    call particles%append (node)
    particles%r(i,:,:)=0.0
    particles%v(i,:,:)=0.0
  end do
  used = wallclock()-used
  print *, n, used, used/n
  !-----------------------------------------------------------------------------
  print *, n
  particles%n = n
  do type=1,4
    used = wallclock()
    call particles%update (type)
    used = wallclock()-used
    print *, type, used, used/n
  end do
  !-----------------------------------------------------------------------------
  allocate (v(m,m,m,3))
  call random%init()
  v = 0.0
  do j=1,3
  do i=1,particles%n
    particles%r(i,j,:) = 1.0 + random%ran1()*(m-1)
  end do
  end do
  do type=1,3
    used = wallclock()
    call particles%interpolate (v, type)
    used = wallclock()-used
  print *, type, used, used/n
  end do
END PROGRAM
