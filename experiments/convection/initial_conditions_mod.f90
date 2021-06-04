!===============================================================================
!> Initial conditions for Whole Sun convection benchmark
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  USE eos_mod
  USE gpatch_mod
  USE random_mod
  USE math_mod
  implicit none
  private
  type, public :: initial_conditions_t
    real:: theta        = 10.0
    real:: mpoly        = 1.0
    real:: amplitude(2) = [1e-3, 0.0]
    integer:: seed      = 11
    logical:: mhd       = .false.
    character(len=64):: solver
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
  end type
  type(random_t):: random
CONTAINS

!===============================================================================
!> Pre-initialize the initial conditions.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_conditions_t%pre_init')
  self%patch => patch
  self%solver = patch%kind
  call random%init (self%seed)
  call trace%end()
END subroutine pre_init

!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t) :: self
  !.............................................................................
  type(index_t), pointer:: idx
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(:,:,:,:), pointer :: ff
  integer:: ix, iy, iz                                       ! loop indices
  real:: z,x,y,zz, pg, p_pert, rhov2, B2, d
  real,dimension(:,:,:),allocatable :: vecpot, emag
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%initial_conditions')
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  if (self%solver(1:6)=='ramses') then
    allocate(emag(m(1)%gn,m(2)%gn,m(3)%gn))
    emag=0.0
  end if
!
  do iz=m(3)%lb,m(3)%ub
    z = max (m(3)%p+m(3)%r(iz)-0.5*m(3)%d, 0.0)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + m(2)%r(iy)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + m(1)%r(ix)
        p_pert = random%ran1()*self%amplitude(1)
        pg = (1.0+p_pert)*(1.0+self%theta*z)**(self%mpoly+1.0)
        d = (1.0+self%theta*z)**self%mpoly
        ff(ix,iy,iz,idx%d) = d
        if (idx%s>0) ff(ix,iy,iz,idx%s) = &
        log(pg*d**-eos%gamma)*d / (eos%gamma-1d0)
        if (idx%e>0) ff(ix,iy,iz,idx%e) = pg/(eos%gamma-1.0)
        ff(ix,iy,iz,idx%px) = 0.0
        ff(ix,iy,iz,idx%py) = 0.0
        ff(ix,iy,iz,idx%pz) = d*self%amplitude(2) &
                            * sin(math%pi*x/2)
        if (idx%bx > 0) then
          ff(ix,iy,iz,idx%bx) = 0.0
          ff(ix,iy,iz,idx%by) = 0.0
          ff(ix,iy,iz,idx%bz) = 0.0
        end if
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  if (self%solver(1:6)=='ramses') then
    if (idx%e>0) then
      do iz=m(3)%lb,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            rhov2 = (ff(ix,iy,iz,idx%px)**2 + &
                     ff(ix,iy,iz,idx%py)**2 + &
                     ff(ix,iy,iz,idx%pz)**2)/ff(ix,iy,iz,idx%d)
            B2 = 0.0
            ff(ix,iy,iz,idx%e) = ff(ix,iy,iz,idx%e) + 0.5*rhov2 + 0.5*B2
          end do
        end do
      end do
    end if
    deallocate(emag)
  end if
  !if (idx%s>0) then
  !do iz=m(3)%lb,m(3)%ub
  !  print *, iz, ff(12,12,iz,idx%d), ff(12,12,iz,idx%s)
  !end do
  !end if
  !-----------------------------------------------------------------------------
  self%patch%boundaries%theta = self%theta
  self%patch%boundaries%mpoly = self%mpoly
  !-----------------------------------------------------------------------------
  call trace_end
  end associate
END SUBROUTINE init

END MODULE initial_conditions_mod
