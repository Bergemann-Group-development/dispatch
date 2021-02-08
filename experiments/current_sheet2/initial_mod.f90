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
    integer :: idirect
    real    :: gamma=1.66667
    logical :: mhd=.true.
    real    :: d0=1.0, beta=0.1, alpha=0.2
    real(8) :: llc(3)
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
  integer, save :: idirect, idir
  real(8), optional:: gamma
  !.............................................................................
  character(len=*):: solver
  real, save :: d0, beta, alpha
  logical, save :: first_time=.true.
  namelist /initial_params/ d0, beta, alpha, idir
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%init')
  if (present(gamma)) self%gamma = gamma

  !$omp critical (initial_t)
  if (first_time) then
    rewind (io%input); read (io%input,initial_params)
    if (io%master) write(*,initial_params)
    first_time = .false.
  end if
  !$omp end critical (initial_t)

  self%d0   = d0
  self%beta = beta
  self%alpha = alpha
  self%solver = trim(solver)
  self%condition => initial_conditions
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
  integer:: ix, iy, iz, n                                  ! loop indices
  real:: entropy, pg, pb, bb, bx_sum, d, bx, by, bz, bxp, byp, bzp
  real(8):: xs, ys, zs, xp, yp, zp
  real,dimension(:,:,:),allocatable :: vecpot
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')
  
  !pg = 0.01
  !bb = sqrt(2*pg/self%beta/4/math%pi) !mu_0 possibly
  bb = 3.!1./sqrt(4./math%pi) !mu_0 possibly
  d = 20.

  do iz=m(3)%lb,m(3)%ub
    zp = m(3)%p + r3(iz)
    zs = zp - 0.5*m(3)%d
    !print*, 'mp, r3', m(3)%p, r3(iz)
    do iy=m(2)%lb,m(2)%ub
      yp = m(2)%p + r2(iy)
      ys = yp - 0.5*m(2)%d
      do ix=m(1)%lb,m(1)%ub
        xp = m(1)%p + r1(ix)
        xs = xp - 0.5*m(1)%d
!        bx_sum = 0.
!        do n=1,25
!          bx_sum=bx_sum+(((1-(yp**2+zp**2)/(d**2*n**2))*((8*zp)/(d**2*n**2)- &
!                        (4*zp*(1-(yp**2+zp**2)/(d**2*n**2)))/(d**2*n**2)))/  &
!                        (n**2*((4*zp**2)/(d**2*n**2)+(1-(yp**2+zp**2)/       &
!                        (d**2*n**2))**2)**2))+(2*zp)/(d**2*n**4*((4*zp**2)/  &
!                        (d**2*n**2)+(1-(yp**2+zp**2)/(d**2*n**2))**2))
!        end do 
        bx = 0.0 !self%alpha*bb*yp*(2/d)**2*bx_sum
        by = 4*self%alpha*bb*ys*zp/(ys**2+zp**2)**2
        bz = -bb*(1+2*self%alpha*(yp**2-zs**2)/(yp**2+zs**2)**2)
        bxp = 0.0 !self%alpha*bb*yp*(2/d)**2*bx_sum
        byp = 4*self%alpha*bb*yp*zp/(yp**2+zp**2)**2
        bzp = -bb*(1+2*self%alpha*(yp**2-zp**2)/(yp**2+zp**2)**2)
        ff(ix,iy,iz,idx%d)=self%d0
        ff(ix,iy,iz,idx%px)=0.0
        ff(ix,iy,iz,idx%py)=0.0
        ff(ix,iy,iz,idx%pz)=0.0
        ff(ix,iy,iz,idx%bx)=bx !0.0!self%alpha*bb*yp*(2/d)**2*bx_sum 
        ff(ix,iy,iz,idx%by)=by !4*self%alpha*bb*ys*zp/(ys**2+zp**2)**2
        ff(ix,iy,iz,idx%bz)=bz !-bb*(1+2*self%alpha*(yp**2-zs**2)/(yp**2+zs**2)**2)
        !pb = 0.5*(bx**2+by**2+bz**2)
        !pb = 0.5*(bxp**2+byp**2+bzp**2)
        pg = self%beta*bb
        if (self%gamma==1.0) then
          ff(ix,iy,iz,idx%s) = 1.0
        else if (self%solver(1:14) == 'zeus_mhd_patch' .or. &
                 self%solver(1:10) == 'stagger2e_') then 
          ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
        else if (self%solver(1:9) == 'stagger2_') then
          ff(ix,iy,iz,idx%s) = &
              (log(pg)-log(self%d0)*self%gamma)/(self%gamma-1.0)*self%d0
        else if (self%solver(1:6) == 'ramses')  then
          ff(ix,iy,iz,idx%e) = pg / (self%gamma - 1.0)
        end if
      end do
    end do
  end do

  call trace_end
  end associate
END SUBROUTINE initial_conditions

END MODULE initial_mod

