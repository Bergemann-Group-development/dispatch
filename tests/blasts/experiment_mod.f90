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
  real, dimension(:,:,:), pointer:: d, px, py, pz, s, e, p, v
  class(mesh_t), pointer:: m
  integer:: ix, iy, iz, iv
  real(8):: gm1
  real, save:: d0=1., d1=0., p0=1., p1=10., w=4., power=4., gamma=1.666667
  namelist /experiment_params/ d0, d1, p0, p1, w, power, gamma
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call trace_begin ('experiment_t%init')
  self%nv = 5
  self%mhd = .false.
  call self%mhd_t%init
  print *,'nv =',self%nv
  self%periodic = .true.
  self%mhd = .false.
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
  ! Variables via pointers
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,idx%d,self%it,1)
  s => self%mem(:,:,:,idx%s,self%it,1)
  e => self%mem(:,:,:,idx%e,self%it,1)
  px=> self%mem(:,:,:,idx%px,self%it,1)
  py=> self%mem(:,:,:,idx%py,self%it,1)
  pz=> self%mem(:,:,:,idx%pz,self%it,1)
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
  ! Spherical initial state
  !-----------------------------------------------------------------------------
  allocate (p(self%gn(1),self%gn(2),self%gn(3)))
  block
    real(8):: x, y, z, r2, w2
    associate (m1=>self%mesh(1), m2=>self%mesh(2), m3=>self%mesh(3))
    w2 = (self%ds(1)*w)**2
    do iz=1,self%gn(3)
      z = m3%p + m3%r(iz) - 0.5d0*self%box(3)
      do iy=1,self%gn(2)
        y = m2%p + m2%r(iy) - 0.5d0*self%box(2)
        do ix=1,self%gn(1)
          x = m1%p + m1%r(ix) - 0.5d0*self%box(1)
          r2 = x**2+y**2+z**2
          d(ix,iy,iz) = d0 + d1*exp(-(r2/w2)**(power/2.))
          p(ix,iy,iz) = p0 + p1*exp(-(r2/w2)**(power/2.))
          px(ix,iy,iz) = 0.0
          py(ix,iy,iz) = 0.0
          pz(ix,iy,iz) = 0.0
        end do
      end do
    end do
    end associate
  end block
  !-----------------------------------------------------------------------------
  ! Energy variable
  !-----------------------------------------------------------------------------
  if (self%kind(1:6) == 'ramses') then
    e = p/gm1 + 0.5*(px**2+py**2+pz**2)
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
