!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    integer:: id
    real(8):: a(3), phase(3)
    real(8):: t_turn, t_save=0.0
    type(random_t):: random
    character(len=64):: solver
    procedure(single_solenoidal), pointer:: selected=>null()
  contains
    procedure:: init
    procedure:: dealloc
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, id, mesh)
   class(force_t):: self
   character(len=64):: solver
   integer:: id, iostat
   class(mesh_t), pointer:: mesh(:)
   character(len=32), save:: type='void'
   integer, save  :: seed=22
   namelist /force_params/ type, seed
   !----------------------------------------------------------------------------
   call trace_begin('force_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params, iostat=iostat)
     if (io%master) write (*, force_params)
   end if
   !$omp end critical (input_cr)
   self%id     = id
   self%solver = trim(solver)
   call self%random%init (seed)
   select case (type)
   case ('single_solenoidal')
     self%selected => single_solenoidal
   case ('void')
     nullify(self%selected)
   case default
     print*,'UNKNOWN FORCING TYPE'
   end select
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Empty procedure; nothing to deallocate
!===============================================================================
SUBROUTINE dealloc (self)
  class(force_t):: self
END SUBROUTINE dealloc

!===============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!===============================================================================
FUNCTION single_solenoidal (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                             :: self
  real(8)                                    :: time
  real, dimension(:,:,:), pointer            :: d, Ux, Uy, Uz
  real, dimension(:,:,:,:), pointer          :: p
  class(mesh_t), dimension(:), pointer       :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3) :: ff
  !.............................................................................
  integer        :: ix, iy, iz, iostat
  real(8)        :: kx, ky, kz
  logical, save  :: first_time=.true.
  real, save     :: t_turn=0.3, a0(3)=0.1, k(3)=1.0
  namelist /force_solenoidal_params/ k, a0, t_turn
  !-----------------------------------------------------------------------------
  call trace_begin('force_t%single_solenoidal')
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, force_solenoidal_params, iostat=iostat)
      if (io%master) write (*, force_solenoidal_params)
    end if
    !$omp end critical (input_cr)
  end if
  if (time >= self%t_save) then
    self%t_save = self%t_save + t_turn
    self%a = a0*2./3.*(1.0+0.1*self%random%ran3())/t_turn
    self%phase = pi2*0.5*self%random%ran3()
    if (self%id==1) print *,'new force at', real(time), self%a, real(self%phase/pi2)
  end if
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r, l=>m%li, u=>m%ui)
  do iz=m(3)%lb,m(3)%ub
    kz = pi2*k(3)*(r3(iz)+m(3)%p)/m(3)%b + self%phase(3)
    do iy=m(2)%lb,m(2)%ub
      ky = pi2*k(2)*(r2(iy)+m(2)%p)/m(2)%b + self%phase(2)
      do ix=m(1)%lb,m(1)%ub
        kx = pi2*k(1)*(r1(ix)+m(1)%p)/m(1)%b + self%phase(1)
        ff(ix,iy,iz,1) = self%a(1)*(sin(ky)+sin(kz))
        ff(ix,iy,iz,2) = self%a(2)*(sin(kz)+sin(kx))
        ff(ix,iy,iz,3) = self%a(3)*(sin(kx)+sin(ky))
      end do
    end do
  end do
  end associate
  call trace_end
END FUNCTION single_solenoidal

END MODULE
