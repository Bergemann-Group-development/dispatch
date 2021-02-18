!===============================================================================
!> Simple initical condition module for tests
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE patch_mod
  USE gpatch_mod
  USE omp_mod
  USE index_mod
  USE eos_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_conditions_t
    integer:: seed
    real:: d0
    logical:: mhd
    real(8):: time=0d0
    character(len=64):: solver
    class(gpatch_t), pointer:: patch => null()
    procedure(single_solenoidal), pointer:: condition=>null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure, nopass:: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  character(len=32), save:: type='void'
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /IC_params/ type
  !----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_params, iostat=iostat)
    call io%namelist_warning ('ramses/hydro/initial_conditions_mod.f90: IC_params', iostat)
    write (io%namelist, IC_params)
 end if
  call io%lock%unset
  select case (trim(type))
  case ('single_solenoidal')
    self%condition => single_solenoidal
  case ('single_compressive')
    self%condition => single_compressive
  case ('single_advection')
    self%condition => single_advection
  case ('void')
    self%condition => void
  case default
    print*,'UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':'
  end select
  call trace_end
END SUBROUTINE pre_init

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%init')
  call self%condition
  call trace_end
END SUBROUTINE init

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!===============================================================================
!> Void initializer
!===============================================================================
SUBROUTINE void (self)
  class(initial_conditions_t) :: self
END SUBROUTINE void

!===============================================================================
!> Initialize a single 3-D wavenumber
!===============================================================================
SUBROUTINE single_solenoidal (self)
  class(initial_conditions_t)       :: self
  class(mesh_t),            pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  type(index_t),            pointer :: idx
  !.............................................................................
  integer        :: ix, iy, iz, iostat
  real(8)        :: kx, ky, kz
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1d0, gamma=1.4d0
  logical, save  :: first_time=.true.
  namelist /IC_solenoidal_params/ k, a0, u0, b0, d0, e0, mhd, gamma, t_turn
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%single_solenoidal', itimer=itimer)
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !.............................................................................
  call io%lock%set
  if (first_time) then
    first_time = .false.
    gamma = eos%gamma
    rewind (io%input)
    read (io%input, IC_solenoidal_params, iostat=iostat)
    call io%namelist_warning ('IC_solenoidal_params', iostat)
    write (io%namelist, IC_solenoidal_params)
  end if
  eos%gamma = gamma
  call io%lock%unset
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0
        ff(ix,iy,iz,idx%e ) = e0
        ff(ix,iy,iz,idx%px) = d0*a0(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,idx%py) = d0*a0(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,idx%pz) = d0*a0(3)*(sin(kx)+sin(ky))
        if (idx%bx > 0) then
          ff(ix,iy,iz,idx%bx) = b0(1)
          ff(ix,iy,iz,idx%by) = b0(2)
          ff(ix,iy,iz,idx%bz) = b0(3)
        end if
      end do
    end do
  end do
  call energy (self, m, idx, ff)
  end associate
  call trace_end (itimer)
END SUBROUTINE single_solenoidal

!===============================================================================
!> Initialize a single 3-D wavenumber, with non-random amplitude and zero phase
!===============================================================================
SUBROUTINE single_compressive (self)
  class(initial_conditions_t)       :: self
  class(mesh_t),            pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  type(index_t),            pointer :: idx
  !.............................................................................
  integer        :: ix, iy, iz, iostat
  real(8)        :: kx, ky, kz
  logical, save  :: mhd=.false.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  real(8), save  :: t_turn=0.1d0, gamma=1.4d0
  logical, save  :: first_time=.true.
  namelist /IC_compressive_params/ k, a0, u0, b0, d0, e0, mhd, gamma, t_turn
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%single_compressive', itimer=itimer)
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !.............................................................................
  call io%lock%set
  if (first_time) then
    first_time = .false.
    gamma = eos%gamma
    rewind (io%input)
    read (io%input, IC_compressive_params,iostat=iostat)
    mhd = mhd .and. any(b0 /= 0.0)
    call io%namelist_warning ('IC_compressive_params', iostat)
    write (io%namelist, IC_compressive_params)
  end if
  eos%gamma = gamma
  call io%lock%unset
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0*(1.0+a0(1)*sin(kx)+a0(2)*sin(ky)+a0(3)*sin(kz))
        ff(ix,iy,iz,idx%e ) = e0
        ff(ix,iy,iz,idx%px) = d0*u0(1)*sin(kx)
        ff(ix,iy,iz,idx%py) = d0*u0(2)*sin(ky)
        ff(ix,iy,iz,idx%pz) = d0*u0(3)*sin(kz)
        if (iy==1.and.iz==1) &  
          print'(i6,5f10.4)',ix,r1(ix),m(1)%p,kx,ff(ix,iy,iz,idx%px)
        if (idx%bx > 0) then
          ff(ix,iy,iz,idx%bx) = b0(1)
          ff(ix,iy,iz,idx%by) = b0(2)
          ff(ix,iy,iz,idx%bz) = b0(3)
        end if
      end do
    end do
  end do
  call energy (self, m, idx, ff)
  end associate
  call trace_end (itimer)
END SUBROUTINE single_compressive

!===============================================================================
!> Initialize a single 3-D wavenumber, with non-random amplitude and phase
!===============================================================================
SUBROUTINE single_advection (self)
  class(initial_conditions_t)       :: self
  class(mesh_t),            pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  type(index_t),            pointer :: idx
  !.............................................................................
  integer        :: ix, iy, iz
  real(8)        :: kx, ky, kz, gamma=1.4d0
  logical, save  :: mhd=.false.
  logical, save  :: first_time=.true.
  real, save     :: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0, phi(3)=0.0
  namelist /IC_advection_params/ k, a0, u0, b0, d0, e0, mhd, gamma
  integer, save  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%single_advection', itimer=itimer)
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !.............................................................................
  call io%lock%set
  if (first_time) then
    first_time = .false.
    gamma = eos%gamma
    rewind (io%input)
    read (io%input, IC_advection_params)
    mhd = mhd .and. any(b0 /= 0.0)
    if (io%master) write (*, IC_advection_params)
  end if
  eos%gamma = gamma
  call io%lock%unset
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz   =  pi2*k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky   =  pi2*k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx   =  pi2*k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + phi(1)
        ff(ix,iy,iz,idx%d ) = d0*(1.+a0(1)*sin(kx)+a0(2)*sin(ky)+a0(3)*sin(kz))
        ff(ix,iy,iz,idx%e ) = e0
        ff(ix,iy,iz,idx%px) = d0*u0(1)
        ff(ix,iy,iz,idx%py) = d0*u0(2)
        ff(ix,iy,iz,idx%pz) = d0*u0(3)
        if (idx%bx > 0) then
          ff(ix,iy,iz,idx%bx) = b0(3)*sin(kz)
          ff(ix,iy,iz,idx%by) = b0(1)*sin(kx)
          ff(ix,iy,iz,idx%bz) = b0(2)*sin(ky)
        end if
      end do
    end do
  end do
  call energy (self, m, idx, ff)
  end associate
  call trace_end (itimer)
END SUBROUTINE single_advection

!===============================================================================
!> Calculate total energy
!===============================================================================
SUBROUTINE energy (self, m, idx, ff)
  class(initial_conditions_t) :: self
  class(mesh_t):: m(3)
  type(index_t):: idx
  real, dimension(:,:,:,:):: ff
  !.............................................................................
  real(8):: fgamma
  !-----------------------------------------------------------------------------
  call trace%begin ('initial_conditions_t%energy')
  if (eos%gamma==1d0) then
    fgamma = 1d0
  else
    fgamma = eos%gamma*(eos%gamma-1d0)
  end if
  ff(:,:,:,idx%e) = 0.5*(ff(:,:,:,idx%px)**2+ &
                         ff(:,:,:,idx%py)**2+ &
                         ff(:,:,:,idx%pz)**2)/ff(:,:,:,idx%d) &
                  + ff(:,:,:,idx%e)
  if (idx%bx > 0) then
    ff(:,:,:,idx%e) = ff(:,:,:,idx%e) + &
      0.5*(ff(:,:,:,idx%bx)**2 + &
           ff(:,:,:,idx%by)**2 + &
           ff(:,:,:,idx%bz)**2)
  end if
  !print *, 'IIII', idx%e, fgamma, &
    !minval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%d)), &
    !maxval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%d)), &
    !minval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%e)), &
    !maxval(ff(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub,idx%e))
  call trace%end()
END SUBROUTINE energy

END MODULE
