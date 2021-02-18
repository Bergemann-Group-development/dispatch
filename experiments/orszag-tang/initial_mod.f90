!===============================================================================
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  USE eos_mod
  implicit none
  private
  type, public :: initial_t
    integer :: idirect
    logical :: mhd=.true.
    real(8) :: p0, v0, b0
    character(len=64):: solver
    procedure(initial_conditions), pointer :: condition => null()
  contains
    procedure :: init
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, solver)
  class(initial_t) :: self
  character(len=*):: solver
  integer, save :: idirect
  !.............................................................................
  real(8), save :: p0, v0, b0
  logical, save :: first_time=.true.
  namelist /IC_params/ idirect, b0, p0, v0
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%init')
  call io%lock%set('initial_t%init')
  if (first_time) then
    ! default values for namelist parameters
    idirect = 3
    b0 = 0.5 / sqrt(math%pi)
    p0 = 5.0 / (12.0 * math%pi)
    v0 = 1.0
    rewind (io%input)
    read (io%input, IC_params)
    write (io%namelist, IC_params)
    first_time = .false.
  end if
  call io%lock%unset('initial_t%init')
  !
  self%idirect = idirect
  self%p0 = p0
  self%v0 = v0
  self%b0 = b0
  self%solver=solver
  !
  self%condition => initial_conditions
  !
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
  real:: entropy, twopi, xa, ya, za, xb, yb, zb, d0 , rhov2, b2         ! local variables
  real,dimension(:,:,:),allocatable :: vecpot, emag
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')

  if (self%solver(1:6)=='ramses') then
    allocate(emag(m(1)%gn,m(2)%gn,m(3)%gn))
    emag=0.0
  end if

  twopi = 2.0 * math%pi
  d0 = eos%gamma * self%p0

  ! "pre"-initialisation of all variables, independent of `idirect`.
  if (eos%gamma==1.0) then
    entropy = 1.0
  else
    entropy = (log(self%p0)-log(d0)*eos%gamma)/(eos%gamma-1.0)
  end if


  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ff(ix,iy,iz,idx%d) = d0
        if (idx%s>0) ff(ix,iy,iz,idx%s) = entropy*d0
        if (idx%e>0) ff(ix,iy,iz,idx%e) = self%p0/(eos%gamma-1.0)
        ff(ix,iy,iz,idx%px) = 0.0
        ff(ix,iy,iz,idx%py) = 0.0
        ff(ix,iy,iz,idx%pz) = 0.0
        ff(ix,iy,iz,idx%bx) = 0.0
        ff(ix,iy,iz,idx%by) = 0.0
        ff(ix,iy,iz,idx%bz) = 0.0
      end do
    end do
  end do
 

  ! initialisation based on value of `idirect`
  ! ---- x1-x2 plane ----
  if (self%idirect == 3) then

    ! initialise vector potential
    allocate(vecpot(m(1)%gn+1,m(2)%gn+1,m(3)%gn))
    
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        ya = m(2)%p + r2(iy) - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(ix) - 0.5*m(1)%d
          vecpot(ix,iy,iz) = (self%b0 / twopi) * ( 0.5 * cos(2.0 * twopi * xa) + cos(twopi  * ya) )
        end do
      end do
    end do
    
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%ub,m(2)%ub+1
        ya = m(2)%p + r2(iy-1) + 0.5*m(2)%d
        do ix=m(1)%ub,m(1)%ub+1
          xa = m(1)%p + r1(ix-1) + 0.5*m(1)%d
          vecpot(ix,iy,iz) = (self%b0 / twopi) * ( 0.5 * cos(2.0 * twopi * xa) + cos(twopi  * ya) )
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        yb = m(2)%p + r2(iy)
        do ix=m(1)%lb,m(1)%ub
          xb = m(1)%p + r1(ix)
          ff(ix,iy,iz,idx%px) = -d0 * self%v0 * sin(twopi * yb)
          ff(ix,iy,iz,idx%py) =  d0 * self%v0 * sin(twopi * xb)
          ff(ix,iy,iz,idx%bx) =  (vecpot(ix,iy+1,iz) - vecpot(ix,iy,iz) ) / m(2)%d
          ff(ix,iy,iz,idx%by) = -(vecpot(ix+1,iy,iz) - vecpot(ix,iy,iz) ) / m(1)%d
        end do
      end do
    end do
    deallocate(vecpot)

  ! ---- x2-x3 plane ----
  else if (self%idirect == 1) then

    ! initialise vector potential
    allocate(vecpot(m(1)%gn,m(2)%gn+1,m(3)%gn+1))
    !
    do iz=m(3)%lb,m(3)%ub
      za = m(3)%p + r3(iz) - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        ya = m(2)%p + r2(iy) - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          vecpot(ix,iy,iz) = self%b0 / twopi * ( 0.5d0 * cos (2.0 * twopi * ya) + cos (twopi  * za) )
        end do
      end do
    end do

    do iz=m(3)%ub,m(3)%ub+1
      za = m(3)%p + r3(iz-1) + 0.5*m(3)%d
      do iy=m(2)%ub,m(2)%ub+1
        ya = m(2)%p + r2(iy-1) + 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          vecpot(ix,iy,iz) = self%b0 / twopi * ( 0.5d0 * cos (2.0 * twopi * ya) + cos (twopi  * za) )
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      zb = m(3)%p + r3(iz)
      do iy=m(2)%lb,m(2)%ub
        yb = m(2)%p + r2(iy)
        do ix=m(1)%lb,m(1)%ub
          ff(ix,iy,iz,idx%py) = -d0 * self%v0 * sin(twopi * zb)
          ff(ix,iy,iz,idx%pz) =  d0 * self%v0 * sin(twopi * yb)
          ff(ix,iy,iz,idx%by) =  (vecpot(ix,iy,iz+1) - vecpot(ix,iy,iz  ) ) / m(3)%d
          ff(ix,iy,iz,idx%bz) = -(vecpot(ix,iy+1,iz) - vecpot(ix,iy  ,iz) ) / m(2)%d
        end do
      end do
    end do

    deallocate(vecpot)
  
  ! ---- x3-x1 plane ----
  else if (self%idirect == 2) then

    ! initialise vector potential
    allocate(vecpot(m(1)%gn+1,m(2)%gn,m(3)%gn+1))
    !
    do iz=m(3)%lb,m(3)%ub
      za = m(3)%p + r3(iz) - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(ix) - 0.5*m(1)%d
          vecpot(ix,iy,iz) = self%b0 / twopi * ( 0.5d0 * cos (2.0 * twopi * za) + cos (twopi  * xa) )
        end do
      end do
    end do

    do iz=m(3)%ub,m(3)%ub+1
      za = m(3)%p + r3(iz-1) + 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%ub,m(1)%ub+1
          xa = m(1)%p + r1(ix-1) + 0.5*m(1)%d
          vecpot(ix,iy,iz) = self%b0 / twopi * ( 0.5d0 * cos (2.0 * twopi * za) + cos (twopi  * xa) )
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      zb = m(3)%p + r3(iz)
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          xb = m(1)%p + r1(ix)
          ff(ix,iy,iz,idx%pz) = -d0 * self%v0 * sin(twopi * xb)
          ff(ix,iy,iz,idx%px) =  d0 * self%v0 * sin(twopi * zb)
          ff(ix,iy,iz,idx%bz) =  (vecpot(ix+1,iy,iz) - vecpot(ix,iy,iz) ) / m(1)%d
          ff(ix,iy,iz,idx%bx) = -(vecpot(ix,iy,iz+1) - vecpot(ix,iy,iz) ) / m(3)%d
        end do
      end do
    end do

    deallocate(vecpot)
  else
    call io%abort('initial%condition: idirect not understood.')
  end if
  
  do iz=m(3)%lb,m(3)%ub
    ff(:,:,iz,:) = ff(:,:,m(3)%li,:)
  end do

  if (self%solver(1:6)=='ramses') then 
    if (idx%e>0) then
      do iz=m(3)%lb,m(3)%ub-1
        do iy=m(2)%lb,m(2)%ub-1
          do ix=m(1)%lb,m(1)%ub-1
            rhov2 = (ff(ix,iy,iz,idx%px)**2 + &
                    ff(ix,iy,iz,idx%py)**2 + &
                    ff(ix,iy,iz,idx%pz)**2)/ff(ix,iy,iz,idx%d) 
            if (idx%bx > 0) then
              B2 = 0.25*(ff(ix,iy,iz,idx%bx) + ff(ix+1,iy,iz,idx%bx))**2 + &
                  0.25*(ff(ix,iy,iz,idx%by) + ff(ix,iy+1,iz,idx%by))**2 + &
                  0.25*(ff(ix,iy,iz,idx%bz) + ff(ix,iy,iz+1,idx%bz))**2
            else
              B2 = 0.0
            end if
            ff(ix,iy,iz,idx%e) = ff(ix,iy,iz,idx%e) + 0.5*rhov2 + 0.5*B2
          end do
        end do
      end do
    end if
    deallocate(emag)
  end if

  call trace_end
  end associate
END SUBROUTINE initial_conditions

END MODULE initial_mod
