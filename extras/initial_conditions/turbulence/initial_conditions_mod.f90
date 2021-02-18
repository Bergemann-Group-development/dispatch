!===============================================================================
!> Initical condition module with a generic turbulence driver
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
  type, public:: initial_conditions_t
    integer:: seed
    logical:: mhd
    real(8):: time=0d0
    class(patch_t), pointer:: patch => null()
    procedure(turbulence), pointer:: condition=>null()
  contains
    procedure:: pre_init
    procedure:: init
  end type
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
    call io%namelist_warning ('IC_params', iostat)
    write (io%namelist, IC_params)
 end if
  call io%lock%unset
  select case (trim(type))
  case ('turbulence')
    self%condition => turbulence
  case ('void')
    self%condition => void
  case default
    print*,'UNKNOWN INITIAL CONDITION TYPE :'//trim(type)//':'
  end select
  call trace_end
END SUBROUTINE pre_init

!===============================================================================
!> Generic entry point
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%init')
  call self%condition
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Void initializer
!===============================================================================
SUBROUTINE void (self)
  class(initial_conditions_t) :: self
END SUBROUTINE void

!===============================================================================
!> Initialize to zero, and let the forces drive from there
!===============================================================================
SUBROUTINE turbulence (self)
  class(initial_conditions_t):: self
  !.............................................................................
  type(index_t),            pointer:: idx
  class(mesh_t),            pointer:: m(:)
  real, dimension(:,:,:,:), pointer:: ff
  real, save                       :: b0(3)=[0.0,0.0,0.1]
  integer                          :: iostat
  logical, save                    :: first_time=.true.
  integer, save                    :: itimer=0
  namelist /IC_turbulence_params/ b0
  !-----------------------------------------------------------------------------
  call trace_begin('initial_conditions_t%turbulence', itimer=itimer)
  m => self%patch%mesh
  ff => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !..............................................................
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, IC_turbulence_params, iostat=iostat)
    call io%namelist_warning ('IC_turbulence_params', iostat)
    write (io%namelist, IC_turbulence_params)
  end if
  call io%lock%unset
  ff(:,:,:,idx%d        ) = 1.0
  ff(:,:,:,idx%e        ) = 1.0
  ff(:,:,:,idx%px:idx%pz) = 0.0
  if (idx%bx > 0) then
    ff(:,:,:,idx%bx       ) = b0(1)
    ff(:,:,:,idx%by       ) = b0(2)
    ff(:,:,:,idx%bz       ) = b0(3)
  end if
  call energy (self, m, idx, ff)
  call trace_end (itimer)
END SUBROUTINE turbulence

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
