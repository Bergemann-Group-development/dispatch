!===============================================================================
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  implicit none
  private
  type, public :: initial_t
    integer :: idir
    real    :: gamma=1.66667
    logical :: mhd=.true.
    real    :: d0, beta, ampl
    character(len=16):: solver
    procedure(initial_conditions), pointer :: condition => null()
  contains
    procedure :: init
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t) :: self
  integer, save :: idir = 3
  real, optional:: gamma
  !.............................................................................
  character(len=*):: solver
  real, save :: d0, beta, ampl
  logical, save :: first_time=.true.
  namelist /IC_params/ d0, beta, ampl, idir
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%init')
  if (present(gamma)) self%gamma = gamma

  d0   = 1.0
  beta = 0.1
  ampl = 0.1

  call io%lock%set ('initial_t%init')
  if (first_time) then
    rewind (io%input); read (io%input,initial_params)
    if (io%master) write(*,initial_params)
    first_time = .false.
  end if
  call io%lock%unset ('initial_t%init')

  self%d0   = d0
  self%beta = beta
  self%ampl = ampl
  self%solver = trim(solver)
  self%condition => initial_conditions
  self%idir = idir

  call trace_end
END SUBROUTINE init

!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE initial_conditions (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  !.............................................................................
  integer:: ix, iy, iz                                       ! loop indices
  real:: entropy, pg, pb, bb
  real(8):: xs, ys, zs, xp, yp, zp
  real,dimension(:,:,:),allocatable :: vecpot
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')
  
  bb=1.0/sqrt(math%pi*4)


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
          if (self%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:10) == 'stagger2e_') then
            ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
          else if (self%solver(1:9) == 'stagger2_') then
            ff(ix,iy,iz,idx%s) = (log(pg)-log(self%d0)*self%gamma)/(self%gamma-1.0)*self%d0
          else if (self%solver(1:6) == 'ramses')  then
            ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
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
          if (self%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:10) == 'stagger2e_') then
            ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
          else if (self%solver(1:9) == 'stagger2_') then
            ff(ix,iy,iz,idx%s) = (log(pg)-log(self%d0)*self%gamma)/(self%gamma-1.0)*self%d0
          else if (self%solver(1:6) == 'ramses')  then
            ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
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
          if (self%gamma==1.0) then
            ff(ix,iy,iz,idx%s) = 1.0
          else if (self%solver(1:14) == 'zeus_mhd_patch' .or. self%solver(1:10) == 'stagger2e_') then
            ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
          else if (self%solver(1:9) == 'stagger2_') then
            ff(ix,iy,iz,idx%s) = (log(pg)-log(self%d0)*self%gamma)/(self%gamma-1.0)*self%d0
          else if (self%solver(1:6) == 'ramses')  then
            ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
          end if
        end do
      end do
    end do
  end select

  call trace_end
  end associate
END SUBROUTINE initial_conditions

END MODULE initial_mod
