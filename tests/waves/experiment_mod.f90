!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Alfven wave problem to test the hlld solver implementation
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE patch_mod
  USE index_mod
  USE solver_mod
  implicit none
  private
  !
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
  end type
  !
  type(index_t):: idx
  !
CONTAINS

!> =============================================================================
!> Setup a patch with an Alfven wave (see Masson et al. 2012)
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d, Bx, By, Bz, px, py, pz, s, e
  real, dimension(:,:,:,:), pointer:: b, p
  integer:: k
  real:: x, y, z, u1, u2, p0
  real(8):: k_w, lambda, gm1, ds2
  real, save:: ampl=0.1, d0=1.0, B0=1.0, csound=1.0
  character(len=16):: type='transverse'
  character(len=1):: direction='z'
  logical, save:: first_time = .true.
  namelist /experiment_params/ type, direction, ampl, d0, B0, csound
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  self%nv = 8
  self%mhd = .true.
  call self%mhd_t%init
  !-----------------------------------------------------------------------------
  ! Mandatory namelist
  !-----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Variables via pointers
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,idx%d,self%it,1)
  s => self%mem(:,:,:,idx%s,self%it,1)
  e => self%mem(:,:,:,idx%e,self%it,1)
  p => self%mem(:,:,:,idx%px:idx%pz,self%it,1)
  B => self%mem(:,:,:,idx%bx:idx%bz,self%it,1)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
  !-----------------------------------------------------------------------------
  ! Wave number
  !-----------------------------------------------------------------------------
  lambda=self%box(3)
  k_w=2.*acos(-1.0d0)/lambda
  !-----------------------------------------------------------------------------
  ! Allow isothermal case (gamma==1)
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    gm1 = 1d0
  else
    gm1 = self%gamma-1d0
  end if
  p0 = d0*csound**2/self%gamma
  !-----------------------------------------------------------------------------
  ! Waves in the x-direction
  !-----------------------------------------------------------------------------
  if (direction=='x') then
    ds2 = self%ds(1)/2d0
    do k=self%mesh(1)%lb,self%mesh(1)%ub
      x = self%mesh(1)%p + self%mesh(1)%r(k)
      u1 = ampl*cos(k_w*x)
      u2 = ampl*sin(k_w*x)
      print '(i4,3f12.6)',k,x,u1,u2
      d (k,:,:)= d0
      if (trim(type)=='transverse') then
        py(k,:,:)= d(k,:,:)*u1
        pz(k,:,:)=-d(k,:,:)*u2
        px(k,:,:)= 0.
        By(k,:,:)= 0.
        Bz(k,:,:)= 0.
        Bx(k,:,:)= B0
      else if (trim(type)=='fast') then
        py(k,:,:)= 0.
        pz(k,:,:)= 0.
        px(k,:,:)= ampl*cos(k_w*(x-ds2))
        By(k,:,:)= B0
        Bz(k,:,:)= 0.0
        Bx(k,:,:)= 0.0
      else if (trim(type)=='sound') then
        py(k,:,:)= 0.
        pz(k,:,:)= 0.
        px(k,:,:)= ampl*cos(k_w*(x-ds2))
        By(k,:,:)= 0.0
        Bz(k,:,:)= 0.0
        Bx(k,:,:)= 0.0
      else
        call io%abort ('experiment_t%init: unknown experiment type')
      end if
      e(k,:,:) = 0.5*d(k,:,:)*(u1**2 + u2**2) + 0.5*(Bx(k,:,:)**2 + By(k,:,:)**2 + Bz(k,:,:)**2)
    end do
  !-----------------------------------------------------------------------------
  ! Waves in the z-direction
  !-----------------------------------------------------------------------------
  else if (direction=='z') then
    ds2 = self%ds(3)/2d0
    do k=self%mesh(3)%lb,self%mesh(3)%ub
      z = self%mesh(3)%p + self%mesh(3)%r(k)
      u1 = ampl*cos(k_w*z)
      u2 = ampl*sin(k_w*z)
      print '(i4,3f12.6)',k,z,u1,u2
      d (:,:,k)= d0
      if (trim(type)=='transverse') then
        px(:,:,k)= d(:,:,k)*u1
        py(:,:,k)=-d(:,:,k)*u2
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
      e(:,:,k) = 0.5*d(:,:,k)*(u1**2 + u2**2) + 0.5*(Bx(:,:,k)**2 + By(:,:,k)**2 + Bz(:,:,k)**2)
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Energy variable
  !-----------------------------------------------------------------------------
  if (self%kind(1:6) == 'ramses') then
    e = p0/gm1 + e
  else if (self%kind(1:9) == 'stagger2e') then
    e = p0/gm1
  else if (self%kind(1:8) == 'stagger2') then
    s = (log(p0)-log(d)*self%gamma)/gm1 * d
  end if
  !
  self%periodic = .true.
  self%mhd = .true.
  print*, 'done!'
  call trace_end
END SUBROUTINE init

END MODULE
