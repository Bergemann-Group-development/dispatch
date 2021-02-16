!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Alfven wave problem to test the hlld solver implementation
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
  USE index_mod
  USE solver_mod
  USE eos_mod
  implicit none
  private
  !
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
  end type
CONTAINS

!===============================================================================
!> Setup a patch with an Alfven wave (see Masson et al. 2012)
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d, Bx, By, Bz, px, py, pz, s, e, pg
  real, dimension(:,:,:,:), pointer:: b, p
  integer:: j, k
  real:: x, y, z, cc, cs, sc, ss, p0
  real(8):: k_w, lambda, gm1, ds2
  real, save:: ampl=0.1, d0=1.0, B0=1.0, csound=1.0
  integer, save:: waves=1, nj=1
  character(len=16):: type='transverse'
  character(len=1):: direction='z'
  logical, save:: first_time=.true., mhd=.false.
  namelist /experiment_params/ type, mhd, direction, ampl, d0, B0, csound, waves, nj
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  !-----------------------------------------------------------------------------
  ! Mandatory namelist
  !-----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  self%nv = 0
  self%mhd = mhd
  call self%mhd_t%init
  print *, 'nv, mhd =', self%nv, self%mhd
  !-----------------------------------------------------------------------------
  ! Variables via pointers
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  s => self%mem(:,:,:,self%idx%s,self%it,1)
  e => self%mem(:,:,:,self%idx%e,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  if  (self%mhd) then
    B => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
  end if
  allocate(pg(size(d,1),size(d,2),size(d,3)))
  !-----------------------------------------------------------------------------
  ! Allow isothermal case (gamma==1)
  !-----------------------------------------------------------------------------
  if (eos%gamma==1.0) then
    gm1 = 1d0
  else
    gm1 = eos%gamma-1d0
  end if
  p0 = d0*csound**2/eos%gamma
  pg = p0
  !-----------------------------------------------------------------------------
  ! Waves in the x-direction
  !-----------------------------------------------------------------------------
  if (direction=='x') then
    !-----------------------------------------------------------------------------
    ! Wave number
    !-----------------------------------------------------------------------------
    lambda=self%box(1)/waves
    k_w=2.*acos(-1.0d0)/lambda
    print *, lambda, k_w
    if (self%kind(1:6) == 'ramses') then
      ds2 = 0.0
    else
      ds2 = self%ds(1)/2d0
    end if
    do k=self%mesh(1)%lb,self%mesh(1)%ub
      x = self%mesh(1)%p + self%mesh(1)%r(k)
      cc=0.0; cs=0.0; sc=0.0; ss=0.0
      do j=1,nj
        cc = cc + ampl*cos((j*2-1)*k_w*x)/(j*2-1)
        sc = sc + ampl*sin((j*2-1)*k_w*x)/(j*2-1)
        cs = cs + ampl*cos((j*2-1)*k_w*(x-ds2))/(j*2-1)
        ss = ss + ampl*sin((j*2-1)*k_w*(x-ds2))/(j*2-1)
      end do
      d (k,:,:)= d0
      if (trim(type)=='transverse') then
        py(k,:,:)= d(k,:,:)*cc
        pz(k,:,:)=-d(k,:,:)*sc
        px(k,:,:)= 0.
        if (self%mhd) then
        Bx(k,:,:)= B0
        By(k,:,:)= 0.
        Bz(k,:,:)= 0.
        end if
      else if (trim(type)=='fast') then
        py(k,:,:)= 0.
        pz(k,:,:)= 0.
        px(k,:,:)= ampl*cos(k_w*(x-ds2))
        if (self%mhd) then
        Bx(k,:,:)= 0.0
        By(k,:,:)= B0
        Bz(k,:,:)= 0.0
        end if
      else if (trim(type)=='sound') then
        py(k,:,:)= 0.
        pz(k,:,:)= 0.
        px(k,:,:)= ampl*cos(k_w*(x-ds2))
        if (self%mhd) then
        Bx(k,:,:)= 0.0
        By(k,:,:)= 0.0
        Bz(k,:,:)= 0.0
        end if
      else if (trim(type)=='sound2') then
        cc=0.0; cs=0.0
        do j=1,nj
          cc = cc + ampl*sin((j*2-1)*k_w*x)/(j*2-1)
          cs = cs + ampl*sin((j*2-1)*k_w*(x-ds2))/(j*2-1)
        end do
        d (k,:,:)= d0*(1.+ampl*cc)
        pg(k,:,:)= csound*d(k,:,:)/eos%gamma
        px(k,:,:)= ampl*d0*csound*cs
        py(k,:,:)= 0.
        pz(k,:,:)= 0.
        if (self%mhd) then
          Bx(k,:,:)= 0.0
          By(k,:,:)= 0.0
          Bz(k,:,:)= 0.0
        end if
      else if (trim(type)=='advect') then
        d (k,:,:)= d0*(1.+1e-3*cc)
        py(k,:,:)= 0.
        pz(k,:,:)= 0.
        px(k,:,:)= ampl*d0*(1.+1e-3*cs)
        if (self%mhd) then
          Bx(k,:,:)= 0.0
          By(k,:,:)= 0.0
          Bz(k,:,:)= 0.0
        end if
      else
        call io%abort ('experiment_t%init: unknown experiment type')
      end if
      print '(i4,5f12.6)',k,x,cc,cs,d(k,1,1),px(k,1,1)
    end do
    if (self%mhd) then
      e = 0.5*(px**2+py**2+pz**2)/d + 0.5*(bx**2+by**2+bz**2)
    else
      e = 0.5*(px**2+py**2+pz**2)/d
    end if
  !-----------------------------------------------------------------------------
  ! Waves in the z-direction
  !-----------------------------------------------------------------------------
  else if (direction=='z') then
    lambda=self%box(3)/waves
    k_w=2.*acos(-1.0d0)/lambda
    print *, lambda, k_w
    ds2 = self%ds(3)/2d0
    do k=self%mesh(3)%lb,self%mesh(3)%ub
      z = self%mesh(3)%p + self%mesh(3)%r(k)
      cc = ampl*cos(k_w*z)
      sc = ampl*sin(k_w*z)
      print '(i4,3f12.6)',k,z,sc,cc
      d (:,:,k)= d0
      if (trim(type)=='transverse') then
        px(:,:,k)= d(:,:,k)*cc
        py(:,:,k)=-d(:,:,k)*sc
        pz(:,:,k)= 0.
        Bx(:,:,k)= 0.
        By(:,:,k)= 0.
        Bz(:,:,k)= B0
      else if (trim(type)=='fast') then
        px(:,:,k)= 0.
        py(:,:,k)= 0.
        pz(:,:,k)= ampl*cos(k_w*(z-ds2))
        Bx(:,:,k)= B0
        By(:,:,k)= 0.0
        Bz(:,:,k)= 0.0
      else if (trim(type)=='sound') then
        px(:,:,k)= 0.
        py(:,:,k)= 0.
        pz(:,:,k)= ampl*cos(k_w*(z-ds2))
        Bx(:,:,k)= 0.0
        By(:,:,k)= 0.0
        Bz(:,:,k)= 0.0
      else
        call io%abort ('experiment_t%init: unknown experiment type')
      end if
    end do
    if (self%mhd) then
      e = 0.5*(px**2+py**2+pz**2)/d + 0.5*(bx**2+by**2+bz**2)
    else
      e = 0.5*(px**2+py**2+pz**2)/d
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Energy variable
  !-----------------------------------------------------------------------------
  if (self%kind(1:6) == 'ramses') then
    e = pg/gm1 + e
  else if (self%kind(1:9) == 'stagger2e') then
    e = pg/gm1
  else if (self%kind(1:8) == 'stagger2') then
    s = (log(pg)-log(d)*eos%gamma)/gm1 * d
  end if
  print *,minval(d),minval(e)
  !
  self%periodic = .true.
  print*, 'done!'
  call trace_end
END SUBROUTINE init

END MODULE
