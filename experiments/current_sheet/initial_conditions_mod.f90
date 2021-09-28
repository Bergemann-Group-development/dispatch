!===============================================================================
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  USE eos_mod
  USE gpatch_mod
  implicit none
  private
  type, public :: initial_conditions_t
    integer :: idir
    logical :: mhd=.true.
    real    :: d0, beta, ampl
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure :: pre_init
    procedure :: init
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t) :: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  integer, save :: idir = 3
  real, save :: d0, beta, ampl
  logical, save :: first_time=.true.
  namelist /IC_params/ d0, beta, ampl, idir
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%init')
  self%patch => patch

  call io%lock%set ('initial_t%init')
  if (first_time) then
      d0   = 1.0
      beta = 0.1
      ampl = 0.1
    rewind (io%input); read (io%input,IC_params)
    if (io%master) write(*,IC_params)
    first_time = .false.
  end if
  call io%lock%unset ('initial_t%init')

  self%d0   = d0
  self%beta = beta
  self%ampl = ampl
  self%idir = idir

  call trace_end
END SUBROUTINE pre_init

!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t) :: self
  ! ............................................................................
  type(index_t), pointer:: idx
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(:,:,:,:), pointer :: ff
  !.............................................................................
  integer:: ix, iy, iz                                       ! loop indices
  real:: entropy, pg, pb, bb, rhov2, b2
  real(8):: xs, ys, zs, xp, yp, zp
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  if (idx%bx > 0) then
    bb=1.0/sqrt(math%pi*4)
  else
    bb=0.0
  end if


  select case(self%idir)
  case(1)  !YZ
    do iz=m(3)%lb,m(3)%ub
      zp = m(3)%p + r3(iz)
      zs = zp - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xp = m(1)%p + r1(ix)
          xs = xp - 0.5*m(1)%d
          ff(ix,iy,iz,idx%d) = self%d0
          ff(ix,iy,iz,idx%px) = 0.0
          ff(ix,iy,iz,idx%py) = self%ampl*sin(2.*math%pi*zp)*self%d0
          ff(ix,iy,iz,idx%pz) = 0.0
          ff(ix,iy,iz,idx%bx) = 0.0
          ff(ix,iy,iz,idx%by) = 0.0
          if (abs(yp)>0.25) then
            ff(ix,iy,iz,idx%bz) = bb
          else
            ff(ix,iy,iz,idx%bz) = -bb
          end if
          !pb = ff(ix,iy,iz,idx%by)**2/4/math%pi
          pg = self%beta/2.*bb
          if (eos%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else if (idx%e>0) then
            ff(ix,iy,iz,idx%e) = pg / (eos%gamma - 1.0)
          else if (idx%s>0) then
            ff(ix,iy,iz,idx%s) = (eos%s0+ log(pg)-log(self%d0)*eos%gamma)/(eos%gamma-1.0)*self%d0
          end if
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          if (idx%e>0) then
            if (self%patch%kind(1:6)=='ramses') then 
              rhov2 = (ff(ix,iy,iz,idx%px)**2 + &
                       ff(ix,iy,iz,idx%py)**2 + &
                       ff(ix,iy,iz,idx%pz)**2)/ff(ix,iy,iz,idx%d) 
              ff(ix,iy,iz,idx%e) = ff(ix,iy,iz,idx%e) + 0.5*rhov2 + 0.5*bb*bb
            end if
          end if
        end do
      end do
    end do

  case(2) !XZ
    do iz=m(3)%lb,m(3)%ub
      zp = m(3)%p + r3(iz)
      zs = zp - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xp = m(1)%p + r1(ix)
          xs = xp - 0.5*m(1)%d
          ff(ix,iy,iz,idx%d) = self%d0
          ff(ix,iy,iz,idx%px) = self%ampl*sin(2.*math%pi*zp)*self%d0
          ff(ix,iy,iz,idx%py) = 0.0
          ff(ix,iy,iz,idx%pz) = 0.0
          ff(ix,iy,iz,idx%bx) = 0.0
          ff(ix,iy,iz,idx%by) = 0.0
          if (abs(xp)>0.25) then
            ff(ix,iy,iz,idx%bz) = bb
          else
            ff(ix,iy,iz,idx%bz) = -bb
          end if
          !pb = ff(ix,iy,iz,idx%by)**2/4/math%pi
          pg = self%beta/2.*bb
          if (eos%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else if (idx%s>0) then
            ff(ix,iy,iz,idx%s) = (eos%s0+log(pg)-log(self%d0)*eos%gamma)/(eos%gamma-1.0)*self%d0
          else if (idx%e>0) then
            ff(ix,iy,iz,idx%e) = pg / (eos%gamma - 1.0)
          end if
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          if (idx%e>0) then
            if (self%patch%kind(1:6)=='ramses') then 
              rhov2 = (ff(ix,iy,iz,idx%px)**2 + &
                       ff(ix,iy,iz,idx%py)**2 + &
                       ff(ix,iy,iz,idx%pz)**2)/ff(ix,iy,iz,idx%d) 
              ff(ix,iy,iz,idx%e) = ff(ix,iy,iz,idx%e) + 0.5*rhov2 + 0.5*bb*bb
            end if
          end if
        end do
      end do
    end do

  case(3) !XY
    do iz=m(3)%lb,m(3)%ub
      zp = m(3)%p + r3(iz)
      zs = zp - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xp = m(1)%p + r1(ix)
          xs = xp - 0.5*m(1)%d
          ff(ix,iy,iz,idx%d) = self%d0
          ff(ix,iy,iz,idx%px) = self%ampl*sin(2.*math%pi*yp)*self%d0
          ff(ix,iy,iz,idx%py) = 0.0
          ff(ix,iy,iz,idx%pz) = 0.0
          ff(ix,iy,iz,idx%bx) = 0.0
          ff(ix,iy,iz,idx%bz) = 0.0
          if (abs(xp)>0.25) then
            ff(ix,iy,iz,idx%by) = bb
          else
            ff(ix,iy,iz,idx%by) = -bb
          end if
          !pb = ff(ix,iy,iz,idx%by)**2/4/math%pi
          pg = self%beta/2.*bb
          if (eos%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else if (idx%s>0) then
            ff(ix,iy,iz,idx%s) = (eos%s0+log(pg)-log(self%d0)*eos%gamma)/(eos%gamma-1.0)*self%d0
          else if (idx%e>0) then
            ff(ix,iy,iz,idx%e) = pg / (eos%gamma - 1.0)
          end if
        end do
      end do
    end do
    
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          if (idx%e>0) then
            if (self%patch%kind(1:6)=='ramses') then 
              rhov2 = (ff(ix,iy,iz,idx%px)**2 + &
                       ff(ix,iy,iz,idx%py)**2 + &
                       ff(ix,iy,iz,idx%pz)**2)/ff(ix,iy,iz,idx%d) 
              ff(ix,iy,iz,idx%e) = ff(ix,iy,iz,idx%e) + 0.5*rhov2 + 0.5*bb*bb
            end if
          end if
        end do
      end do
    end do

  end select

  call trace_end
  end associate
END SUBROUTINE init

END MODULE initial_conditions_mod
