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
    integer :: idirect
    logical :: mhd=.true.
    real(8) :: y0, d1, d2, deltay, deltavy, dymin, dymax
    character(len=64):: solver
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
  end type
CONTAINS

!===============================================================================
!> Pre-initialize the initial conditions.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  integer, save:: idirect
  real(8), save :: y0, d1, d2, deltay, deltavy, dymin, dymax
  logical, save :: first_time=.true.
  namelist /IC_params/ y0, d1, d2, deltay, deltavy, dymin, dymax
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%pre_init')
  self%patch => patch
  self%solver = patch%kind
  self%mhd = .true.
  call io%lock%set('initial_conditions_t%pre_init')
  if (first_time) then
    ! default values for namelist parameters
    y0 = 0.5
    d1 = 1.0
    d2 = 2.0
    deltay = 0.025
    deltavy = 0.025
    dymin = 0.3
    dymax = 0.7    

    rewind (io%input)
    read (io%input, IC_params)
    write (io%namelist, IC_params)
    first_time = .false.
  end if
  call io%lock%unset('initial_conditions_t%pre_init')
  self%y0 = y0
  self%d1 = d1
  self%d2 = d2
  self%deltay = deltay
  self%deltavy = deltavy
  self%dymin = dymin
  self%dymax = dymax
  !
  call trace%end()
END subroutine pre_init


!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t) :: self
  ! ----------------------------------------------------------------------------
  type(index_t), pointer:: idx
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(:,:,:,:), pointer :: ff
  !.............................................................................
  integer:: ix, iy, iz                                       ! loop indices
  real:: entropy, twopi, xa, ya, za, xb, yb, zb, d0, p0, rhov2, b2         ! local variables
  real(8):: y_grav_acc=-0.5
  real,dimension(:,:,:),allocatable :: emag
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')

  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)

  if (self%solver(1:6)=='ramses') then
    allocate(emag(m(1)%gn,m(2)%gn,m(3)%gn))
    emag=0.0
  end if

  twopi = 2.0 * math%pi
  !d0 = eos%gamma * self%p0


  ! "pre"-initialisation of all variables, independent of `idirect`.
  

  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ff(ix,iy,iz,idx%d) = 0.0 !d0
        !if (idx%s>0) ff(ix,iy,iz,idx%s) = entropy*d0
        !if (idx%e>0) ff(ix,iy,iz,idx%e) = p0/(eos%gamma-1.0)
        ff(ix,iy,iz,idx%px) = 0.0
        ff(ix,iy,iz,idx%py) = 0.0
        ff(ix,iy,iz,idx%pz) = 0.0
        ff(ix,iy,iz,idx%bx) = 0.0
        ff(ix,iy,iz,idx%by) = 0.0
        ff(ix,iy,iz,idx%bz) = 0.0
      end do
    end do
  end do
 

    
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        ya = m(2)%p + r2(iy) - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(ix) - 0.5*m(1)%d
          
          d0 = self%d1 + (self%d2 - self%d1) / (1.0 + exp((self%y0 -ya)/self%deltay))
          ff(ix,iy,iz,idx%d) = d0

          if (ya >= self%dymin .and. ya <= self%dymax) then
            ff(ix,iy,iz,idx%py) = self%deltavy * (1.0 + cos(8.0 * math%pi * (xa +  0.25))) &
                                  * (1.0 + cos(5.0 * math%pi * (ya - 0.50)))* d0
          endif                     
         
          p0 = self%d2 / eos%gamma + y_grav_acc *d0 * (ya - self%y0)        
 
          if (eos%gamma==1.0) then
            entropy = 1.0
          else
            entropy = (eos%s0+log(p0)-log(d0)*eos%gamma)/(eos%gamma-1.0)
          end if

          if (idx%s>0) ff(ix,iy,iz,idx%s) = entropy*d0
          if (idx%e>0) ff(ix,iy,iz,idx%e) = p0/(eos%gamma-1.0)

        end do
      end do
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
END SUBROUTINE init

END MODULE initial_conditions_mod
