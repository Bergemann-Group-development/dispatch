!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  USE stagger_mod
  USE index_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: id=0
    integer:: seed
    logical:: mhd
    real(8):: gamma=1d0, alpha
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
    procedure, nopass:: cleanup
  end type
  logical, save:: first_time=.true.
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  type(index_t):: idx
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self, solver, gamma)
   class(initial_t):: self
   character(len=64):: solver
   real, optional:: gamma
   !............................................................................
   logical, save:: mhd=.true.
   real, save:: alpha=0.8
   character(len=32), save:: type='politano+95'
   namelist /initial_params/ alpha
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, initial_params)
     if (io%master) write (*, initial_params)
   end if
   self%solver = trim(solver)
   !$omp end critical (input_cr)
   self%mhd = .true.
   self%alpha = alpha
   call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!> =============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!> =============================================================================
SUBROUTINE condition (self, m, ff)
  class(initial_t)      :: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t:single_solenoidal')
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b
        ff(ix,iy,iz,idx%d ) = 1.0
        ff(ix,iy,iz,idx%e ) = 1.0
        ff(ix,iy,iz,idx%px) = -2.0*sin(ky)
        ff(ix,iy,iz,idx%py) = +2.0*sin(kx)
        ff(ix,iy,iz,idx%pz) = 0.0
        ff(ix,iy,iz,idx%bx) = self%alpha*(-2.*sin(2.*ky)+sin(kz))
        ff(ix,iy,iz,idx%by) = self%alpha*(+2.*sin(2.*kx)+sin(kz))
        ff(ix,iy,iz,idx%bz) = self%alpha*(    sin(   kx)+sin(ky))
      end do
    end do
  end do
  end associate
  call trace_end
END SUBROUTINE condition

END MODULE
