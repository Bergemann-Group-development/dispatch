!===============================================================================
!> Simple initial condition module for tests.  The variables used in this
!> solver are: density, momentum, entropy per unit volume, and magnetic vector 
!> potential.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE omp_mod
  USE index_mod
  USE eos_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  integer, parameter:: dp=8
  logical, save:: isothermal=.false.
  type, public:: initial_t
    integer:: seed
    real:: d0
    real(8):: time=0d0
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
    procedure, nopass:: cleanup
  end type
  integer, save:: verbose=0
  real(8), save:: g1, g2
  character(len=20):: type='simple_advection'
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver)
  class(initial_t):: self
  character(len=64):: solver
  character(len=32), save:: type='void'
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /IC_params/ type, verbose
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, IC_params, iostat=iostat)
    call io%namelist_warning ('muscl_A%initial_mod.f90: IC_params', iostat)
    write (io%output, IC_params)
    if (eos%gamma == 1d0) then
      g1 = 1d0
    else
      g1 = (eos%gamma-1d0)
    end if
    g2 = eos%gamma*g1
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  call trace_end
END SUBROUTINE init

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!===============================================================================
!> Initialize a single 3-D wavenumber, or constant values.
!===============================================================================
SUBROUTINE condition (self, m, ff, idx)
  class(initial_t) :: self
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:) :: ff
  type(index_t):: idx
  !.............................................................................
  integer        :: ix, iy, iz, iv, iostat
  real           :: d
  real(8)        :: x, y, z, kx, ky, kz, pb, bx, by, bz
  class(mesh_t), pointer:: m1, m2, m3
  real(8), pointer  :: r1(:), r2(:), r3(:)
  !.............................................................................
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, s0=0.0, b0(3)=0.0, u0(3)=0.0, &
                    p0=1.0, phi(3)=0.0, c=5.0
  namelist /IC_simple_params/ k, p0, a0, u0, b0, d0, s0, u0, c
  logical, save  :: first_time=.true.
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%condition', itimer=itimer)
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, IC_simple_params, iostat=iostat)
    call io%namelist_warning ('IC_simple_params', iostat)
    write (io%output, IC_simple_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  r1=>m(1)%r
  r2=>m(2)%r
  r3=>m(3)%r
  do iz=m(3)%lb,m(3)%ub
    z = r3(iz)+m(3)%p-0.5d0*m(3)%s
    kz = pi2*k(3)*z/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      y = r2(iy)+m(2)%p-0.5d0*m(2)%s
      ky = pi2*k(2)*y/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        x = r1(ix)+m(1)%p-0.5d0*m(1)%s
        kx = pi2*k(1)*x/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%s ) = s0
        select case (trim(type))
        case ('single_solenoidal')
          ff(ix,iy,iz,idx%d ) = d0
          ff(ix,iy,iz,idx%px) = d0*(a0(1)*(sin(ky)+sin(kz)) + u0(1))
          ff(ix,iy,iz,idx%py) = d0*(a0(2)*(sin(kz)+sin(kx)) + u0(2))
          ff(ix,iy,iz,idx%pz) = d0*(a0(3)*(sin(kx)+sin(ky)) + u0(3))
          if (size(ff,4) >= 8) then
            ff(ix,iy,iz,idx%ax) = b0(2)*z - b0(3)*y
            ff(ix,iy,iz,idx%ay) = b0(3)*x - b0(1)*z
            ff(ix,iy,iz,idx%az) = b0(1)*y - b0(2)*x
          end if
        case ('single_compressive')
          ff(ix,iy,iz,idx%d ) = d0
          ff(ix,iy,iz,idx%px) = d0*u0(1)*sin(kx)
          ff(ix,iy,iz,idx%py) = d0*u0(2)*sin(ky)
          ff(ix,iy,iz,idx%pz) = d0*u0(3)*sin(kz)
          if (size(ff,4) >= 8) then
            ff(ix,iy,iz,idx%ax) = b0(2)*z - b0(3)*y
            ff(ix,iy,iz,idx%ay) = b0(3)*x - b0(1)*z
            ff(ix,iy,iz,idx%az) = b0(1)*y - b0(2)*x
          end if
        case ('simple_advection')
          d = d0*(1.+a0(1)*atan(c*sin(kx)) &
                    +a0(2)*atan(c*sin(ky)) &
                    +a0(3)*atan(c*sin(kz)))
          ff(ix,iy,iz,idx%d ) = d
          ff(ix,iy,iz,idx%px) = u0(1)*d
          ff(ix,iy,iz,idx%py) = u0(2)*d
          ff(ix,iy,iz,idx%pz) = u0(3)*d
          if (size(ff,4) >= 8) then
            ff(ix,iy,iz,idx%ax) = b0(2)*sin(kz) - b0(3)*sin(ky)
            ff(ix,iy,iz,idx%ay) = b0(3)*sin(kx) - b0(1)*sin(kz)
            ff(ix,iy,iz,idx%az) = b0(1)*sin(ky) - b0(2)*sin(kx)
          end if
          bx = b0(1)*(cos(ky)-cos(kz))
          by = b0(2)*(cos(kz)-cos(kx))
          bz = b0(3)*(cos(kx)-cos(ky))
          pb = 0.5*(bx**2+by**2+bz**2)
          ff(ix,iy,iz,idx%s ) = (log(p0-pb)-eos%gamma*log(d))/g1
        case ('single_advection')
          ff(ix,iy,iz,idx%d ) = d0*(1.+a0(1)*sin(kx)+a0(2)*sin(ky)+a0(3)*sin(kz))
          ff(ix,iy,iz,idx%px) = d0*u0(1)
          ff(ix,iy,iz,idx%py) = d0*u0(2)
          ff(ix,iy,iz,idx%pz) = d0*u0(3)
          if (size(ff,4) >= 8) then
            ff(ix,iy,iz,idx%ax) = b0(2)*sin(kz) - b0(3)*sin(ky)
            ff(ix,iy,iz,idx%ay) = b0(3)*sin(kx) - b0(1)*sin(kz)
            ff(ix,iy,iz,idx%az) = b0(1)*sin(ky) - b0(2)*sin(kx)
          end if
        case default
          call io%abort ('UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':')
        end select
      end do
    end do
  end do
  call trace_end (itimer)
END SUBROUTINE condition

END MODULE
