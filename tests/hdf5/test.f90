!> =============================================================================
!> Demo code that creates a file hdf5.dat with a 'dispatch' dataset, or opens it
!> if it exists, and writes an array to it.
!> =============================================================================
PROGRAM test
  USE hdf5
  implicit none
  integer:: error
  integer(HID_T):: file_id, dset_id, dspc_id
  real, dimension(16,16,16):: a
  character(len=16):: file
  logical:: exist
  !-----------------------------------------------------------------------------
  a = 1.0
  call h5open_f (error)
  file='hdf5.dat'
  inquire (file=file, exist=exist)
  if (exist) then
    print *,'updating '//trim(file)
    call h5fopen_f (file, H5F_ACC_RDWR_F, file_id, error)
    call h5dopen_f (file_id, 'dispatch', dset_id, error)
  else
    print *,'creating '//trim(file)
    call h5fcreate_f (file, H5F_ACC_EXCL_F, file_id, error)
    call h5screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), dspc_id, error)
    call h5dcreate_f (file_id, 'dispatch', H5T_NATIVE_INTEGER, dspc_id, dset_id, error)
  end if
  call h5dwrite_f (dset_id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), error)
  call h5dclose_f (dset_id, error)
  call h5close_f (error)
END PROGRAM test
