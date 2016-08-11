!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: seed
    logical:: mhd
    real(8):: a(3), b(3), phi(3), d0, e0
    real(8):: a0(3), k(3), b0(3), u0(3)
    real(8):: time
    procedure(raw), pointer:: selected=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
  logical, save:: first_time=.true.
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
   class(initial_t):: self
   integer, save:: seed=-33
   logical, save:: mhd=.true.
   real, save:: a0(3)=0.1, k(3)=1.0, d0=1.0, e0=1.0, b0(3)=0.0, u0(3)=0.0
   character(len=32), save:: type='single_solenoidal'
   namelist /initial_params/ seed, k, a0, u0, b0, d0, e0, mhd, type
   !----------------------------------------------------------------------------
   call trace_begin('initial_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, initial_params)
     mhd = any(b0 /= 0.0)
     if (io%master) write (*, initial_params)
   end if
   !$omp end critical (input_cr)
   self%seed   = seed
   self%k      = k
   self%a0     = a0
   self%b0     = b0
   self%u0     = u0
   self%d0     = d0
   self%e0     = e0
   self%mhd    = mhd
   self%time   = 0d0
   select case (type)
   case ('raw')
     self%selected => raw
   case default
     print*,'UNKNOWN INITIAL CONDITION TYPE'
     stop
   end select
   call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!> =============================================================================
!> Initialize from a raw data file
!> =============================================================================
FUNCTION raw (self, m) RESULT (ff)
  USE io_unit_mod
  class(initial_t):: self
  type(mesh_t):: m(3)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8):: ff
  !.............................................................................
  integer:: n(3), i(3), iv
  integer, save:: mv=4
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t::raw')
  self%mhd = .false.
  self%time = 0.02
  !-----------------------------------------------------------------------------
  ! At first call, allocate a buffer and read the entire file -- this is OK since
  ! we only use one thread per rank to initialize, computing also how many calls
  ! will be made, so the buffer can be deallocated when all patches are done
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    n = nint(m%b/m%d)
    allocate (buffer(n(1),n(2),n(3),mv))
    open (io_unit%data, file='data/ramses256/t0.02.raw', form='unformatted', &
      access='direct', recl=4*product(n))
    do iv=1,mv
      read (io_unit%data, rec=iv) buffer(:,:,:,iv)
    end do
    close (io_unit%data)
  end if
  i = 1 + nint(m%p/m%d)-m%n/2                            ! index offsets
  ff = 0.0
  do iv=1,mv
    ff(m(1)%li:m(1)%ui, &
       m(2)%li:m(2)%ui, &
       m(3)%li:m(3)%ui,iv) = buffer(i(1):i(1)+m(1)%n-1, &
                                    i(2):i(2)+m(2)%n-1, &
                                    i(3):i(3)+m(3)%n-1,iv)
  end do
  do iv=2,4
    ff(:,:,:,iv) = ff(:,:,:,1)*ff(:,:,:,iv)
  end do
  if (io%verbose > 2) &
    write (io_unit%log,'(a,3i4,1p,2g12.3)') 'ix,iy,iz,dmin,dmax =', i, &
      minval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1)), &
      maxval(ff(m(1)%li:m(1)%ui,m(2)%li:m(2)%ui,m(3)%li:m(3)%ui,1))
  call trace_end
END FUNCTION raw

END MODULE
