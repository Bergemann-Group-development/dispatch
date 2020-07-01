!===============================================================================
!> Program to test functionality and speed of particles
!===============================================================================

MODULE pic_mod
  implicit none
  private
  type, public:: particles_t
    integer:: np=100
    real, pointer, contiguous:: r(:,:), v(:,:)
    integer(kind=2), pointer, contiguous:: q(:,:)
  end type
  type, public:: cells_t
    integer:: n=32, n_bin=5, mp=100
    real(8):: time=0d0, dtime=0d0
    type(particles_t), pointer:: particles(:,:,:,:)
    real, pointer, contiguous:: s(:)
    real, pointer, contiguous:: c(:,:,:)
    real, pointer, contiguous, dimension(:,:,:,:):: a, v
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

SUBROUTINE init (self)
  class(cells_t):: self
  integer:: n, mp, np, nb
  integer:: i1, i2, i3, ib
  type(particles_t), pointer:: p
  n = self%n
  nb = self%n_bin
  allocate (self%particles(n,n,n,nb))
  allocate (self%v(n,n,n,3))
  allocate (self%a(n,n,n,3))
  allocate (self%c(n,n,n))
  allocate (self%s(nb))
  mp = self%mp
  do ib=1,nb
    self%s(ib) = 1.0
    do i3=1,n
    do i2=1,n
    do i1=1,n
      p => self%particles(i1,i2,i3,ib)
      allocate (p%r(mp,3),p%v(mp,3), p%q(mp,3))
      p%r(:,:) = 0.0
      p%v(:,:) = 0.0
      self%c(i1,i2,i3) = 1.0
    end do
    end do
    end do
  end do
END SUBROUTINE

SUBROUTINE update (self)
  class(cells_t):: self
  integer:: n, nb, np
  integer:: i1, i2, i3, ib, i, j, q
  integer(8):: nu
  real:: dt, cs, a, v, c
  real(8):: wc, omp_get_wtime
  type(particles_t), pointer:: p
  real, pointer, contiguous:: pv(:,:), pr(:,:)
  n  = self%n
  nb = self%n_bin
  dt = self%dtime
  wc = omp_get_wtime()
  nu = 0
  do i3=1,self%n
  do i2=1,self%n
  do i1=1,self%n
    cs = self%c(i1,i2,i3)
    do ib=1,nb
      c = cs/self%s(ib)
      p => self%particles(i1,i2,i3,ib)
      np = p%np
      do i=1,3
        a = self%a(i,i1,i2,i3)
        v = self%v(i,i1,i2,i3)
        !$omp simd
        do j=1,np
          p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
          p%r(j,i) = p%r(j,i) + dt*p%v(j,i)*2.0
          p%v(j,i) = p%v(j,i) + dt*(a+c*(v-p%v(j,i)))/(1.+c*dt)
          q = floor(p%r(j,i))
          p%r(j,i) = p%r(j,i) - q
          p%q(j,i) = p%q(j,i) + q
        end do
      end do
      nu = nu + np
    end do
  end do
  end do
  end do
  wc = omp_get_wtime() - wc
  print '(f6.3," s",6p,f8.3," mu/p",i10," updates")', wc, wc/nu, nu
END SUBROUTINE
END MODULE

PROGRAM test
  USE pic_mod
  implicit none
  type(cells_t):: cells
  call cells%init
  call cells%update
END PROGRAM
