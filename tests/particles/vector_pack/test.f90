!===============================================================================
!> Program to test functionality and speed of particles
!===============================================================================

MODULE pic_mod
  implicit none
  private
  type, public:: particles_t
    integer:: np=30
    real, pointer, contiguous:: r(:,:), v(:,:)
    integer(kind=2), pointer, contiguous:: q(:,:)
  end type
  type, public:: cells_t
    integer:: n=32, n_bin=5, mp=60
    real(8):: time=0d0, dtime=1d0
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
  integer:: n, mp, np, nb, nmax
  integer:: i1, i2, i3, ib, j
  type(particles_t), pointer:: p
  real:: v
  n = self%n
  nb = self%n_bin
  allocate (self%particles(n,n,n,nb))
  allocate (self%v(3,n,n,n))
  allocate (self%a(3,n,n,n))
  allocate (self%c(n,n,n))
  allocate (self%s(nb))
  mp = self%mp
  nmax = 0.7*n
  print *,(real(nmax)/n)**3,' fraction'
  do ib=1,nb
    self%s(ib) = 1.0
    do i3=1,n
    do i2=1,n
    do i1=1,n
      p => self%particles(i1,i2,i3,ib)
      allocate (p%r(mp,3),p%v(mp,3), p%q(mp,3))
      if (i1 < nmax .and. i2 < nmax .and. i3 < nmax) then
        v = 0.5
      else
        v = 0.0
      end if
      do j=1,mp
        p%r(j,:) = 0.5
        p%v(j,:) = v
        p%q(j,:) = [i1,i2,i3]
      end do
      self%c(i1,i2,i3) = 0.0
      self%a(:,i1,i2,i3) = 0.0
      self%v(:,i1,i2,i3) = 0.0
    end do
    end do
    end do
  end do
END SUBROUTINE

SUBROUTINE update (self)
  class(cells_t):: self
  integer:: n, nb, np
  integer:: i1, i2, i3, ib, i, j, q, jp
  integer:: q1, q2, q3
  integer(8):: nu
  real:: dt, cs, a, v, c
  real(8):: wc, wu, omp_get_wtime
  type(particles_t), pointer:: p, pn
  real, pointer, contiguous:: pv(:,:), pr(:,:)
  n  = self%n
  nb = self%n_bin
  dt = self%dtime
  wc = omp_get_wtime()
  nu = 0
  !-----------------------------------------------------------------------------
  ! Kick-drift-kick update
  !-----------------------------------------------------------------------------
  do i3=1,self%n
  do i2=1,self%n
  do i1=1,self%n
    cs = self%c(i1,i2,i3)
    do ib=1,nb
      c = cs/self%s(ib)
      c = 0.0
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
  wu = omp_get_wtime() - wc
  print '(f6.3," s",6p,f8.3," mu/p",i10," KDK")', wu, wu/nu, nu
  !-----------------------------------------------------------------------------
  ! Re-assign particles to other cells
  !-----------------------------------------------------------------------------
  wc = omp_get_wtime()
  n = 0
  do i3=1,self%n
  do i2=1,self%n
  do i1=1,self%n
    do ib=1,nb
      p => self%particles(i1,i2,i3,ib)
      np = p%np
      do j=1,np
        if (any(p%q(j,:) /= [i1,i2,i3])) then
          n = n+1
          q1 = p%q(j,1)
          q2 = p%q(j,2)
          q3 = p%q(j,3)
          pn => self%particles(q1,q2,q3,ib)
          !pn => self%particles(p%q(j,1),p%q(j,2),p%q(j,3),ib)
          pn%np = pn%np + 1
          jp = pn%np
          do i=1,3
            pn%r(jp,i) = p%r(j,i)
            pn%v(jp,i) = p%v(j,i)
            pn%q(jp,i) = p%q(j,i)
          end do
          p%q(j,1) = 0
        end if
      end do
    end do
  end do
  end do
  end do
  wu = wu + omp_get_wtime() - wc
  print '(f6.3," s",6p,f8.3," mu/p",i10," move",i10)', wu, wu/nu, nu, n
  !-----------------------------------------------------------------------------
  ! Remove the holes in the arrays
  !-----------------------------------------------------------------------------
  wc = omp_get_wtime()
  n = 0
  do i3=1,self%n
  do i2=1,self%n
  do i1=1,self%n
    do ib=1,nb
      p => self%particles(i1,i2,i3,ib)
      np = p%np
      j = 1
      do while (j <= p%np)
        jp = p%np
        if (p%q(j,1) == 0) then
          p%q(j,:) = p%q(jp,:)
          p%r(j,:) = p%r(jp,:)
          p%v(j,:) = p%v(jp,:)
          p%np = p%np-1
          n = n+1
        end if
        j = j+1
      end do
    end do
  end do
  end do
  end do
  wu = wu + omp_get_wtime() - wc
  print '(f6.3," s",6p,f8.3," mu/p",i10," pack",i10)', wu, wu/nu, nu, n

END SUBROUTINE
END MODULE

PROGRAM test
  USE pic_mod
  implicit none
  type(cells_t):: cells
  call cells%init
  call cells%update
END PROGRAM
