!===============================================================================
!> Module to handle allocation and deallocation of memory, with accounting.
!> Allows generic calls for 1-D to 6-D for reals, and for 1-D to 4-D for double
!> and integers.
!>
!> FIXME: The accounting could/should be moved from os_mod to here, and the
!>        interface io%bits_mem should be deprecated.
!===============================================================================
MODULE memory_mod
  USE io_mod
  USE os_mod
  implicit none
  private
  type, public:: memory_t
  contains
    procedure, nopass:: allocate_1d_r
    procedure, nopass:: allocate_2d_r
    procedure, nopass:: allocate_3d_r
    procedure, nopass:: allocate_4d_r
    procedure, nopass:: allocate_5d_r
    procedure, nopass:: allocate_6d_r
    procedure, nopass:: allocate_1d_d
    procedure, nopass:: allocate_2d_d
    procedure, nopass:: allocate_3d_d
    procedure, nopass:: allocate_4d_d
    procedure, nopass:: allocate_1d_i
    procedure, nopass:: allocate_2d_i
    procedure, nopass:: allocate_3d_i
    procedure, nopass:: allocate_4d_i
    generic:: allocate => &
      allocate_1d_r, allocate_2d_r, allocate_3d_r, allocate_4d_r, &
      allocate_1d_d, allocate_2d_d, allocate_3d_d, allocate_4d_d, &
      allocate_1d_i, allocate_2d_i, allocate_3d_i, allocate_4d_i, &
      allocate_5d_r, allocate_6d_r
    procedure, nopass:: deallocate_1d_r
    procedure, nopass:: deallocate_2d_r
    procedure, nopass:: deallocate_3d_r
    procedure, nopass:: deallocate_4d_r
    procedure, nopass:: deallocate_5d_r
    procedure, nopass:: deallocate_6d_r
    procedure, nopass:: deallocate_1d_d
    procedure, nopass:: deallocate_2d_d
    procedure, nopass:: deallocate_3d_d
    procedure, nopass:: deallocate_4d_d
    procedure, nopass:: deallocate_1d_i
    procedure, nopass:: deallocate_2d_i
    procedure, nopass:: deallocate_3d_i
    procedure, nopass:: deallocate_4d_i
    generic:: deallocate => &
      deallocate_1d_r, deallocate_2d_r, deallocate_3d_r, deallocate_4d_r, &
      deallocate_1d_d, deallocate_2d_d, deallocate_3d_d, deallocate_4d_d, &
      deallocate_1d_i, deallocate_2d_i, deallocate_3d_i, deallocate_4d_i, &
      deallocate_5d_r, deallocate_6d_r
    procedure, nopass:: info
  end type
  integer:: verbose=0
  type(memory_t), public:: memory
CONTAINS

!===============================================================================
!> Allocate a 1-D real array
!===============================================================================
SUBROUTINE allocate_1d_r (a, m, zero)
  real(4), pointer, contiguous:: a(:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_r

!===============================================================================
!> Allocate a 2-D real array
!===============================================================================
SUBROUTINE allocate_2d_r (a, m, zero)
  real(4), pointer, contiguous:: a(:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_r

!===============================================================================
!> Allocate a 3-D real array
!===============================================================================
SUBROUTINE allocate_3d_r (a, m, zero)
  real(4), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_r

!===============================================================================
!> Allocate a 4-D real array
!===============================================================================
SUBROUTINE allocate_4d_r (a, m, zero)
  real(4), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_r

!===============================================================================
!> Allocate a 5-D real array
!===============================================================================
SUBROUTINE allocate_5d_r (a, m, zero)
  real(4), pointer, contiguous:: a(:,:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_5d_r

!===============================================================================
!> Allocate a 6-D real array
!===============================================================================
SUBROUTINE allocate_6d_r (a, m, zero)
  real(4), pointer, contiguous:: a(:,:,:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4),m(5),m(6)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_6d_r

!===============================================================================
!> Allocate a 1-D integer array
!===============================================================================
SUBROUTINE allocate_1d_i (a, m, zero)
  integer(4), pointer, contiguous:: a(:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_i

!===============================================================================
!> Allocate a 2-D integer array
!===============================================================================
SUBROUTINE allocate_2d_i (a, m, zero)
  integer(4), pointer, contiguous:: a(:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_i

!===============================================================================
!> Allocate a 3-D integer array
!===============================================================================
SUBROUTINE allocate_3d_i (a, m, zero)
  integer(4), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_i

!===============================================================================
!> Allocate a 4-D integer array
!===============================================================================
SUBROUTINE allocate_4d_i (a, m, zero)
  integer(4), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_i

!===============================================================================
!> Allocate a 1-D double array
!===============================================================================
SUBROUTINE allocate_1d_d (a, m, zero)
  real(8), pointer, contiguous:: a(:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_1d_d

!===============================================================================
!> Allocate a 2-D double array
!===============================================================================
SUBROUTINE allocate_2d_d (a, m, zero)
  real(8), pointer, contiguous:: a(:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_2d_d

!===============================================================================
!> Allocate a 3-D double array
!===============================================================================
SUBROUTINE allocate_3d_d (a, m, zero)
  real(8), pointer, contiguous:: a(:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_3d_d

!===============================================================================
!> Allocate a 4-D double array
!===============================================================================
SUBROUTINE allocate_4d_d (a, m, zero)
  real(8), pointer, contiguous:: a(:,:,:,:)
  integer:: m(:)
  logical, optional:: zero
  !-----------------------------------------------------------------------------
  allocate (a(m(1),m(2),m(3),m(4)))
  call io%bits_mem (storage_size(a), product(shape(a)), 'mem')
  if (option(zero)) a = 0.0
END SUBROUTINE allocate_4d_d

!===============================================================================
!> Dealloccate a 1-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_r (a)
  real(4), pointer, contiguous:: a(:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_r

!===============================================================================
!> Dealloccate a 2-D integer array
!===============================================================================
SUBROUTINE deallocate_2d_r (a)
  real(4), pointer, contiguous:: a(:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_r

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_r (a)
  real(4), pointer, contiguous:: a(:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_r

!===============================================================================
!> Dealloccate a 4-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_r (a)
  real(4), pointer, contiguous:: a(:,:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_r

!===============================================================================
!> Dealloccate a 5-D integer array
!===============================================================================
SUBROUTINE deallocate_5d_r (a)
  real(4), pointer, contiguous:: a(:,:,:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_5d_r

!===============================================================================
!> Dealloccate a 6-D integer array
!===============================================================================
SUBROUTINE deallocate_6d_r (a)
  real(4), pointer, contiguous:: a(:,:,:,:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_6d_r

!===============================================================================
!> Dealloccate a 1-D integer array
!===============================================================================
SUBROUTINE deallocate_1d_i (a)
  integer(4), pointer, contiguous:: a(:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_i

!===============================================================================
!> Dealloccate a 2-D integer array
!===============================================================================
SUBROUTINE deallocate_2d_i (a)
  integer(4), pointer, contiguous:: a(:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_i

!===============================================================================
!> Dealloccate a 3-D integer array
!===============================================================================
SUBROUTINE deallocate_3d_i (a)
  integer(4), pointer, contiguous:: a(:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_i

!===============================================================================
!> Dealloccate a 4-D integer array
!===============================================================================
SUBROUTINE deallocate_4d_i (a)
  integer(4), pointer, contiguous:: a(:,:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_i

!===============================================================================
!> Dealloccate a 1-D double array
!===============================================================================
SUBROUTINE deallocate_1d_d (a)
  real(8), pointer, contiguous:: a(:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_1d_d

!===============================================================================
!> Dealloccate a 2-D double array
!===============================================================================
SUBROUTINE deallocate_2d_d (a)
  real(8), pointer, contiguous:: a(:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_2d_d

!===============================================================================
!> Dealloccate a 3-D double array
!===============================================================================
SUBROUTINE deallocate_3d_d (a)
  real(8), pointer, contiguous:: a(:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_3d_d

!===============================================================================
!> Dealloccate a 4-D double array
!===============================================================================
SUBROUTINE deallocate_4d_d (a)
  real(8), pointer, contiguous:: a(:,:,:,:)
  !-----------------------------------------------------------------------------
  if (associated(a)) then
    call io%bits_mem (storage_size(a), -product(shape(a)), 'mem')
    deallocate (a)
  end if
END SUBROUTINE deallocate_4d_d

!===============================================================================
!> Print memory use, gigabytes.  Expected use:
!>
!>   call memory%info ('some label', verbose)
!>
!> For verbose==0 no checking is done
!> For verbose==1 or missing a log message is printed if memory has increased.
!> For verbose==2 or larger a log message is always printed
!===============================================================================
SUBROUTINE info (label, verbose)
  real(8):: dgb
  real(8), save:: gb=0d0
  integer, optional:: verbose
  character(*), optional:: label
  integer:: verbose_L
  !.............................................................................
  if (present(verbose)) then
    verbose_L = verbose
  else
    verbose_L = 1
  end if
  if (verbose_L > 0) then
    call io%lock%set
    if (gb > 0d0) then
      dgb = os%gb-gb
      if (verbose_L > 1 .or. dgb > 0d0) then
        if (present(label)) then
          write (stdout,'(1x,a,f15.9,a,f15.9,a)') &
            'memory allocated:', os%gb,' GB, change:', dgb, ' GB, at '//trim(label)
        else
          write (stdout,'(1x,a,f15.9,a,f15.9,a)') &
            'memory allocated:', os%gb,' GB, change:', dgb, ' GB'
        end if
      end if
    end if
    gb = os%gb
  end if
  call io%lock%unset
END SUBROUTINE info

!===============================================================================
!> Local utility
!===============================================================================
LOGICAL FUNCTION option (o)
  logical, optional:: o
  if (present(o)) then
    option = o
  else
    option = .false.
  end if
END FUNCTION option

END MODULE memory_mod

