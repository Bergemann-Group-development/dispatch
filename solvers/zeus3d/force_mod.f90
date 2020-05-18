!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE random_mod
  USE kinds_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    integer:: seed=-11
    real(8):: a(3), phi(3)
    real(8):: t_turn, t_save=0.0
    real(8):: a0(3), k(3)
    character(len=64):: solver
    procedure(single_solenoidal), pointer:: selected=>null()
  contains
    procedure init
  end type
 logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver)
   class(force_t):: self
   character(len=64):: solver
   integer,save:: seed=-11, iostat
   character(len=32), save:: type='zero'
   real, save:: t_turn=0.3, a0(3)=0.1, k(3)=1.0
   namelist /force_params/ k, a0, t_turn, type
   !----------------------------------------------------------------------------
   call trace_begin('force_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params)
     if (io%master) write (*, force_params)
   end if
   call io%lock%unset
   self%seed   = seed
   self%t_turn = t_turn
   self%a0     = a0
   self%k      = k
   select case (type)
   case ('single_solenoidal')
     self%selected => single_solenoidal
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN FORCING TYPE'
   end select
   self%solver = trim(solver)
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Forcing on a single 3-D wavenumber, changing amplitude and phase periodically
!===============================================================================
FUNCTION single_solenoidal (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer             :: p
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  !.............................................................................
  type(random_t)        :: random
  integer               :: ix, iy, iz
  real(8)               :: kx, ky, kz
  !-----------------------------------------------------------------------------
  call trace_begin('force_t:single_solenoidal')
  if (time >= self%t_save) then
    self%a = self%a0*random%ran3()
    self%phi = pi2*random%ran3()
    self%t_save = self%t_save + self%t_turn
!     if (io%verbose>3) print *,'forcing', self%a, self%phi
  end if

  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    kz = pi2*self%k(3)*(r3(iz)+m(3)%p-0.5d0*m(3)%s)/m(3)%b + self%phi(3)
    do iy=m(2)%lb,m(2)%ub
      ky = pi2*self%k(2)*(r2(iy)+m(2)%p-0.5d0*m(2)%s)/m(2)%b + self%phi(2)
      do ix=m(1)%lb,m(1)%ub
        kx = pi2*self%k(1)*(r1(ix)+m(1)%p-0.5d0*m(1)%s)/m(1)%b + self%phi(1)
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
