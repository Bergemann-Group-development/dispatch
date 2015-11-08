!> =============================================================================
MODULE io_m
  implicit none
  private
  integer             :: data_unit = 52
  integer, parameter  :: flen = 200
  character(len=flen) :: ext='.dat'
  !
  ! Default values 
  !
  character(len=flen), public :: fname='snapshot'

  public :: write_field, open_file, close_file
CONTAINS
!> =============================================================================
SUBROUTINE write_field(f)
  implicit none
  real, dimension(:,:,:) :: f
  integer                :: n(3)
  n = shape(f)
  write(data_unit) n
  write(data_unit) f
END SUBROUTINE write_field

!> =============================================================================
SUBROUTINE open_file(cnt)
  implicit none
  integer :: cnt
  character(len=flen) :: file
  !
  write(file,'(i5.5)') cnt
  file = trim(fname)//'_'//trim(file)//trim(ext)
  !
  open(data_unit, file=trim(file), form='unformatted', status='unknown')
END SUBROUTINE

!> =============================================================================
SUBROUTINE close_file
  implicit none
  !
  close(data_unit)
END SUBROUTINE
    
END MODULE
!> =============================================================================
