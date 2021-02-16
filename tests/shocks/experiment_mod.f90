!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE index_mod
  USE solver_mod
  USE stagger_mod
  USE mesh_mod
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
!> Setup patches with shock conditions
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d, Bx, By, Bz, px, py, pz, s, e, p, v
  class(mesh_t), pointer:: m
  integer:: k, iv
  real(8):: x, ds2, gm1
  character(len=16):: type='standoff'
  character(len=1):: direction='x'
  real(4), save:: P_l=1.0 , P_r=1.0 , d_l=1.0 , d_r=1.0 , &
                  vx_l=5.0, vx_r=-5., vy_l=0.0, vy_r=0.0, vz_l=0.0, vz_r=0.0, &
                  Bx_l=1.0, Bx_r=1.0, By_l=1.0, By_r=0.0, Bz_l=0.0, Bz_r=0.0, &
                  gamma=1.0, x_interface=0.5
  namelist /experiment_params/ type, direction, &
                               P_l, P_r, d_l, d_r, &
                               vx_l, vx_r, vy_l, vy_r, vz_l, vz_r, &
                               Bx_l, Bx_r, By_l, By_r, Bz_l, Bz_r, &
                               gamma, x_interface
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call trace_begin ('experiment_t%init')
  self%nv = 8
  self%mhd = .true.
  call self%mhd_t%init
  self%periodic = .true.
  self%mhd = .true.
  !-----------------------------------------------------------------------------
  ! Mandatory namelist
  !-----------------------------------------------------------------------------
  if (first_time) then
    gamma = eos%gamma
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
    first_time = .false.
  end if
  eos%gamma = gamma
  !-----------------------------------------------------------------------------
  ! Fool proof
  !-----------------------------------------------------------------------------
  call io%assert(Bx_L==Bx_r, "Shock tubes in the x-direction must have a constant Bx!")
  !-----------------------------------------------------------------------------
  ! Variables via pointers
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,idx%d,self%it,1)
  s => self%mem(:,:,:,idx%s,self%it,1)
  e => self%mem(:,:,:,idx%e,self%it,1)
  px=> self%mem(:,:,:,idx%px,self%it,1)
  py=> self%mem(:,:,:,idx%py,self%it,1)
  pz=> self%mem(:,:,:,idx%pz,self%it,1)
  bx=> self%mem(:,:,:,idx%bx,self%it,1)
  by=> self%mem(:,:,:,idx%by,self%it,1)
  bz=> self%mem(:,:,:,idx%bz,self%it,1)
  allocate (p(self%gn(1),self%gn(2),self%gn(3)))
  !-----------------------------------------------------------------------------
  ! Allow isothermal case (gamma==1)
  !-----------------------------------------------------------------------------
  if (eos%gamma==1.0) then
    gm1 = 1d0
  else
    gm1 = eos%gamma-1d0
  end if
  !-----------------------------------------------------------------------------
  ! Waves in the x-direction
  !-----------------------------------------------------------------------------
  if (direction=='x') then
    ds2 = self%ds(1)/2d0
    do k=self%mesh(1)%lb,self%mesh(1)%ub
      x = self%mesh(1)%p + self%mesh(1)%r(k)
      if (trim(type)=='standoff') then
        if (x < x_interface) then
          d (k,:,:)= d_l
          p (k,:,:)= p_l
          px(k,:,:)= vx_l
          py(k,:,:)= vy_l
          pz(k,:,:)= vz_l
          By(k,:,:)= Bx_l
          Bz(k,:,:)= By_l
          Bx(k,:,:)= Bz_l
        else
          d (k,:,:)= d_r
          p (k,:,:)= p_r
          px(k,:,:)= vx_r
          py(k,:,:)= vy_r
          pz(k,:,:)= vz_r
          By(k,:,:)= Bx_r
          Bz(k,:,:)= By_r
          Bx(k,:,:)= Bz_r
        end if
      else
        call io%abort ('experiment_t%init: unknown experiment type')
      end if
    end do
  else
    call io%abort ('experiment_t%init: only x-direction supported')
  end if
  !-----------------------------------------------------------------------------
  ! Energy variable
  !-----------------------------------------------------------------------------
  if (self%kind(1:6) == 'ramses') then
    e = p/gm1 + 0.5*(px**2+py**2+pz**2) + 0.5*(bx**2+by**2+bz**2)
  else if (self%kind(1:9) == 'stagger2e') then
    e = p/gm1
  else if (self%kind(1:8) == 'stagger2') then
    s = (log(p)-log(d)*eos%gamma)/gm1 * d
  end if
  deallocate (p)
  !-----------------------------------------------------------------------------
  ! Stagger variables
  !-----------------------------------------------------------------------------
  call self%periodic_grid
  do iv=1,self%nv
    v => self%mem(:,:,:,iv,self%it,1)
    m => self%mesh(1)
    if (m%h(iv) /= 0.0) then
      v= xdn(v)
      if (iv == idx%px) then
        px = xdn(px)*exp(xdn(log(d)))
      end if
    end if
    m => self%mesh(2)
    if (m%h(iv) /= 0.0) then
      v = ydn(v)
      if (iv == idx%py) then
        py = py*exp(ydn(log(d)))
      end if
    end if
    m => self%mesh(3)
    if (m%h(iv) /= 0.0) then
      v = zdn(v)
      if (iv == idx%pz) then
        pz = pz*exp(zdn(log(d)))
      end if
    end if
  end do
  !-----------------------------------------------------------------------------
  call trace_end
END SUBROUTINE init

END MODULE
