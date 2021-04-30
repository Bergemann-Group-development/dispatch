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
  USE random_mod
  implicit none
  private
  type, public :: initial_conditions_t
    logical :: mhd=.true.
    character(len=64):: solver
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
  end type
  real, save:: theta=10.0
  real, save:: sigma=1.0
  type(random_t):: random
CONTAINS

!===============================================================================
!> Pre-initialize the initial conditions.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  integer, save:: idirect
  real(8), save :: p0, v0, b0
  logical, save :: first_time=.true.
  namelist /IC_params/ theta, sigma
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%pre_init')
  self%patch => patch
  self%solver = patch%kind
  self%mhd = .true.
  call io%lock%set('initial_conditions_t%pre_init')
  if (first_time) then
    ! default values for namelist parameters
    rewind (io%input)
    read (io%input, IC_params)
    write (io%namelist, IC_params)
    first_time = .false.
  end if
  call io%lock%unset('initial_conditions_t%pre_init')
  !
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
  real:: mpoly=1.0
  real:: amplitude=1e-3
  real:: z,zz, pg, p_pert, rhov2, B2
  real,dimension(:,:,:),allocatable :: vecpot, emag
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%initial_conditions')
  call random%init
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  if (self%solver(1:6)=='ramses') then
    allocate(emag(m(1)%gn,m(2)%gn,m(3)%gn))
    emag=0.0
  end if
! d0 = eos%gamma * self%p0
! "pre"-initialisation of all variables, independent of `idirect`.
! if (eos%gamma==1.0) then
!   entropy = 1.0
! else
!   entropy = (log(self%p0)-log(d0)*eos%gamma)/(eos%gamma-1.0)
! end if
!
!
  do iz=m(3)%lb,m(3)%ub
    z = max(m(3)%p + m(3)%r(iz),0.0)
    !z=1.0-zz
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        p_pert = random%ran1()*amplitude
        pg = (1.0+p_pert)*(1.0+theta*z)**(mpoly+1.0)
        ff(ix,iy,iz,idx%d) = (1.0+theta*z)**mpoly
        if (idx%e>0) ff(ix,iy,iz,idx%e) = pg/(eos%gamma-1.0)
        ff(ix,iy,iz,idx%px) = 0.0
        ff(ix,iy,iz,idx%py) = 0.0
        ff(ix,iy,iz,idx%pz) = 0.0
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
  !-----------------------------------------------------------------------------
  self%patch%boundaries%theta=theta
  self%patch%boundaries%mpoly=mpoly
  !-----------------------------------------------------------------------------
  call trace_end
  end associate
END SUBROUTINE init

END MODULE initial_conditions_mod
