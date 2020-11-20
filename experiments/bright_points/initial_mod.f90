!===============================================================================
!> Read a snapshot given by Klaus Galsgaard and convert to DISPATCH data structures.
!> 
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private

  type, public :: initial_t
    integer :: expe
    logical :: mhd
    real(8) :: time
    logical :: pic_needed
    logical :: first_pic
    procedure(initial_conditions), pointer :: init_mhd_vars => null()
  contains
    procedure :: init
    procedure, nopass :: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()

CONTAINS
!===============================================================================
!> Facilitate call to `ramses%read`.  All heavy lifting is done there
!===============================================================================
SUBROUTINE init (self)
  class(initial_t):: self
  logical:: file_exist
  !-----------------------------------------------------------------------------
  inquire (file='data.fe', exist=file_exist)
  if (file_exist) then
    self%init_mhd_vars => read_conditions
  else
    self%init_mhd_vars => initial_conditions
  end if
END SUBROUTINE init

!===============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!===============================================================================
!> This routine should give initial conditions when starting from scratch
!===============================================================================
SUBROUTINE initial_conditions (self, m, ff, idx)
  class(initial_t):: self
  type(index_t):: idx
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  !-----------------------------------------------------------------------------
  error stop 'bright_points/initial_mod::initial_conditions, no IC data or IC code'
END SUBROUTINE initial_conditions

!===============================================================================
!> Reads a snapshot from Klaus and store them in Dispatch variables.  No need
!> to be fancy here; we just hardwire the file dimensions and record size unit.
!===============================================================================
SUBROUTINE read_conditions (self, m, ff)
  class(initial_t):: self
  class(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  !---- variables for reading the snapshot ----
  integer(kind=4)    :: ilun=10
  logical            :: first_time = .true.
  integer            :: i(3), j(3), l(3), u(3), iv, iostat
  integer, save      :: itimer=0
  integer, parameter :: rrlen=4 ! This could also be 1, depending on the compiler
  integer, save      :: nx=1170, ny=351, nz=351
  character(len=64), save :: filename='data.fe'
  namelist /IC_params/ nx, ny, nz, filename
  !.............................................................................
  call trace_begin ('initial_t%init_mhd_vars', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! First time, read the entire file into a buffer to select patch data from
  !-----------------------------------------------------------------------------
  !$omp critical (io_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, initial_params, iostat=iostat)
    if (io%master) write (*, initial_params)
    print*, "000000000000 mesh", m%gn 
    allocate (buffer(nx,ny,nz,8))
    ! Here we read d,px,py,pz,e,Bx,By,Bz from the file
    open(ilun, file=trim(filename), status='old', form='unformatted', access='direct', &
      recl=rrlen*nx*ny*nz)
    do iv=1,8
      read(ilun,rec=iv) buffer(:,:,:,iv)
      print *,'init_mhd_vars: iv,min,max =', &
        iv,minval(buffer(:,:,:,iv)),maxval(buffer(:,:,:,iv))
    end do
    close (ilun)
  end if
  !$omp end critical (io_cr)
  !-----------------------------------------------------------------------------
  ! Select patch data; beware of convention for position zero point !
  !-----------------------------------------------------------------------------
  i = 1 + nint((m%p-m%origin)/m%d - m%o) + m%li
  j = i + m%n - 1 
  print *, m%p/m%d, m%o, i, j
  l = m%li
  u = m%ui
  ff(l(1):u(1),l(2):u(2),l(3):u(3),1) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),1)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),2) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),5)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),3) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),2)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),4) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),3)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),5) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),4)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),6) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),6)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),7) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),7)
  ff(l(1):u(1),l(2):u(2),l(3):u(3),8) = buffer(i(1):j(1),i(2):j(2),i(3):j(3),8)
  print 1,'init', i, m%p, minval(ff(l(1):u(1),l(2):u(2),l(3):u(3),1)), &
                          maxval(ff(l(1):u(1),l(2):u(2),l(3):u(3),1)), &
                          minval(ff(l(1):u(1),l(2):u(2),l(3):u(3),6)), &
                          maxval(ff(l(1):u(1),l(2):u(2),l(3):u(3),6))
  1 format(a,5x,"offset:",3i5,5x,"position:",3f10.2,5x,"dmin/max:",1p,2(2e12.3,2x))
  !-----------------------------------------------------------------------------
  call trace_end
END SUBROUTINE read_conditions
 
END MODULE initial_mod
