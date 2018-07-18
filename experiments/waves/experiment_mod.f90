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
  real, dimension(:,:,:,:), pointer:: b, p
  real, dimension(:,:,:), pointer:: d, Bx, By, Bz, px, py, pz, s, e
  integer:: k
  real:: z, ux, uy, p0
  real(8):: k_w, lambda, gm1
  real, save:: ampl=0.1, d0=1.0, B0=1.0, csound=1.0
  character(len=16):: type='transverse'
  logical, save:: first_time = .true.
  namelist /experiment_params/ ampl, d0, B0, type

  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  print*, 'start init'
  call trace_begin ('experiment_t%init')
  self%nv = 8
  self%mhd = .true.
  call self%mhd_t%init
 
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if

  d => self%mem(:,:,:,idx%d,self%it,1)
  s => self%mem(:,:,:,idx%s,self%it,1)
  e => self%mem(:,:,:,idx%e,self%it,1)
  p => self%mem(:,:,:,idx%px:idx%pz,self%it,1)
  B => self%mem(:,:,:,idx%bx:idx%bz,self%it,1)
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);

  print*, 'gn is:', self%gn(1)
  print*,'self%mesh%lo!!!',self%mesh%lo
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

  do k=self%mesh(3)%lb,self%mesh(3)%ub
    z = self%mesh(3)%p + self%mesh(3)%r(k)
    ux = ampl*cos(k_w*z)
    uy = ampl*sin(k_w*z)
    d (:,:,k)= d0
    if (trim(type)=='transverse') then
      px(:,:,k)= d(:,:,k)*ux
      py(:,:,k)=-d(:,:,k)*uy
      pz(:,:,k)= 0.
      Bx(:,:,k)= B0*ampl*cos(k_w*z) 
      By(:,:,k)=-B0*ampl*sin(k_w*z)
      Bz(:,:,k)= B0
    else if (trim(type)=='fast') then
      px(:,:,k)= 0.
      py(:,:,k)= 0.
      pz(:,:,k)= ampl*cos(k_w*(z-0.5*self%ds(3)))
      Bx(:,:,k)= B0
      By(:,:,k)= 0.0
      Bz(:,:,k)= 0.0
    else if (trim(type)=='sound') then
      px(:,:,k)= 0.
      py(:,:,k)= 0.
      pz(:,:,k)= ampl*cos(k_w*(z-0.5*self%ds(3)))
      Bx(:,:,k)= 0.0
      By(:,:,k)= 0.0
      Bz(:,:,k)= 0.0
    else
      call io%abort ('experiment_t%init: unknown experiment type')
    end if
  end do

  if (self%kind(1:6) == 'ramses') then
    e= p0/gm1 + 0.5*d*(ux**2 + uy**2) + 0.5 * (Bx**2 + By**2 + Bz**2)  
  else if (self%kind(1:9) == 'stagger2e') then
    e= p0/gm1
  else if (self%kind(1:8) == 'stagger2') then
    s= (log(p0)-log(d)*self%gamma)/gm1 * d
  end if

  self%periodic = .true.
  self%mhd = .true.
  print*, 'done!'

  print*,'self%mem(1,1,1,:,self%it,1) in experiment setup',self%mem(1,1,1,:,self%it,1)
  if (self%initial%time > 0d0) then
    self%time = self%initial%time
    self%t(self%it) = self%initial%time
  end if

  call trace_end
END SUBROUTINE init

END MODULE
