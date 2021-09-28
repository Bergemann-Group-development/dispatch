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
    real(8) :: b0, ux, uy, uz
    class(gpatch_t), pointer:: patch
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
  real(8), save :: b0, ux, uy, uz
  integer, save:: idirect
  logical, save :: first_time=.true.
  namelist /IC_params/ idirect, b0, ux, uy, uz 
  !-----------------------------------------------------------------------------

  call trace_begin('initial_conditions_t%pre_init')
  call io%lock%set('initial_conditions_t%pre_init')
  if (first_time) then
    ! default values for namelist parameters
    idirect = 3
    b0 = 0.001
    ux = 1.0
    uy = 1.0
    uz = 0.0
    rewind (io%input); read (io%input,IC_params)
    if (io%master) write(*,IC_params)
    first_time = .false.
  end if
  call io%lock%unset('initial_conditions_t%pre_init')

  self%patch => patch
  self%mhd = .true.

  self%idirect = idirect
  self%b0 = b0
  self%ux = ux
  self%uy = uy
  self%uz = uz

  call trace_end
END SUBROUTINE pre_init

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
  real:: entropy, xa, ya, za, d0, p0, radius     ! local variables
  real,dimension(:,:,:),allocatable :: vecpot, emag
  !.............................................................................
  call trace_begin ('initial_t%initial_conditions')

  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)

  if (self%patch%kind(1:6)=='ramses') then
    allocate(emag(m(1)%gn,m(2)%gn,m(3)%gn))
    emag=0.0
  end if

  d0 = 1.0
  p0 = 1.0

  ! "pre"-initialisation of all variables, independent of `idirect`.
  if (eos%gamma==1.0) then
    entropy = 1.0
  else
    entropy = (eos%s0+log(p0)-log(d0)*eos%gamma)/(eos%gamma-1.0)
  end if


  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ff(ix,iy,iz,idx%d) = d0
        if (idx%s>0) ff(ix,iy,iz,idx%s) = entropy*d0
        if (idx%e>0) ff(ix,iy,iz,idx%e) = p0/(eos%gamma-1.0)
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
    
    allocate(vecpot(m(1)%gn+1, m(2)%gn+1, m(3)%gn))
    
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        ya = m(2)%p + r2(iy) - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(ix) - 0.5*m(1)%d

          radius = sqrt(ya**2 + xa**2) 
          if ( radius .le. 0.3 ) then 
              vecpot(ix,iy,iz) = self%b0*(0.3-radius)
          else
              vecpot(ix,iy,iz) = 0.0
          endif

        end do
      end do
    end do
    
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%ub, m(2)%ub+1
        ya = m(2)%p + r2(iy-1) + 0.5*m(2)%d 
        do ix=m(1)%ub, m(1)%ub+1
          xa = m(1)%p + r1(ix-1) + 0.5*m(1)%d 

          radius = sqrt(ya**2 + xa**2) 
          if ( radius .le. 0.3 ) then 
              vecpot(ix,iy,iz) = self%b0*(0.3-radius)
          else
              vecpot(ix,iy,iz) = 0.0
          endif

        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          ff(ix,iy,iz,idx%px) =  d0 * self%ux
          ff(ix,iy,iz,idx%py) =  d0 * self%uy
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
          
          radius = sqrt(ya**2 + za**2) 
          if ( radius .le. 0.3 ) then 
              vecpot(ix,iy,iz) = self%b0*(0.3-radius)
          else
              vecpot(ix,iy,iz) = 0.0
          endif
           
        end do
      end do
    end do

    do iz=m(3)%ub,m(3)%ub+1
      za = m(3)%p + r3(iz-1) + 0.5*m(3)%d
      do iy=m(2)%ub,m(2)%ub+1
        ya = m(2)%p + r2(iy-1) + 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          
          radius = sqrt(ya**2 + za**2) 
          if ( radius .le. 0.3 ) then 
              vecpot(ix,iy,iz) = self%b0*(0.3-radius)
          else
              vecpot(ix,iy,iz) = 0.0
          endif
         
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          ff(ix,iy,iz,idx%py) =  d0 * self%uy
          ff(ix,iy,iz,idx%pz) =  d0 * self%uz
          ff(ix,iy,iz,idx%by) =  (vecpot(ix,iy,iz+1) - vecpot(ix,iy,iz) ) /m(3)%d
          ff(ix,iy,iz,idx%bz) = -(vecpot(ix,iy+1,iz) - vecpot(ix,iy,iz) ) /m(2)%d
        end do
      end do
    end do

    deallocate(vecpot)

  ! ---- x3-x1 plane ----
  else if (self%idirect == 2) then

    ! initialise vector potential
    allocate(vecpot(m(1)%gn+1,m(2)%gn,m(3)%gn+1))
    
    do iz=m(3)%lb,m(3)%ub
      za = m(3)%p + r3(iz) - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(ix) - 0.5*m(1)%d
         
          radius = sqrt(xa**2 + za**2) 
          if ( radius .le. 0.3 ) then 
              vecpot(ix,iy,iz) = self%b0*(0.3-radius)
          else
              vecpot(ix,iy,iz) = 0.0
          endif 

        end do
      end do
    end do

    do iz=m(3)%ub,m(3)%ub+1
      za = m(3)%p + r3(iz-1) + 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%ub,m(1)%ub+1
          xa = m(1)%p + r1(ix-1) + 0.5*m(1)%d
          
          radius = sqrt(xa**2 + za**2) 
          if ( radius .le. 0.3 ) then 
              vecpot(ix,iy,iz) = self%b0*(0.3-radius)
          else
              vecpot(ix,iy,iz) = 0.0
          endif
           
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          ff(ix,iy,iz,idx%pz) = d0 * self%uz
          ff(ix,iy,iz,idx%px) = d0 * self%ux
          ff(ix,iy,iz,idx%bz) =  (vecpot(ix+1,iy,iz) - vecpot(ix,iy,iz)) /m(1)%d
          ff(ix,iy,iz,idx%bx) = -(vecpot(ix,iy,iz+1) - vecpot(ix,iy,iz)) /m(3)%d
        end do
      end do
    end do

    deallocate(vecpot)

  endif


  if (self%patch%kind(1:6)=='ramses') deallocate(emag)

  call trace_end
  end associate
END SUBROUTINE init

END MODULE initial_conditions_mod

