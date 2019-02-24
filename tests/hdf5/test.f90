!> =============================================================================
!> Demo code that creates a file hdf5.dat with a 'dispatch' dataset, or opens it
!> if it exists, and writes an array to it.
!> =============================================================================
PROGRAM test
  USE hdf5
  implicit none
  integer:: error
  integer(HID_T):: file_id, set_id, dspc_id, grp_id, seq_id
  real, dimension(16,16,16):: a
  character(len=16):: file
  logical:: new
  integer:: verbose=1
  !-----------------------------------------------------------------------------
  file='hdf5.dat'
  file_id = h5_open (file, new)
  if (verbose >= 0) &
    print *, merge('creating ','updating ',new)//trim(file)
  grp_id = h5_group_open (file_id, 'sequence')
  seq_id = h5_group_open (grp_id, '000000')
  set_id = h5_set_open (seq_id, 'value', a)
  a = 1.0
  call h5_set_write (set_id, a)
  call h5_set_close (set_id)
  call h5_close (file_id)
contains
  !-----------------------------------------------------------------------------
  ! Interface function to open or create an HDF5 file
  !-----------------------------------------------------------------------------
  integer(HID_T) function h5_open (name, new) result (fid)
  character(len=*):: name
  integer:: error
  logical, optional:: new
  logical:: exists
  !.............................................................................
  call H5open_f (error)
  inquire (file=name, exist=exists)
  if (exists) then
    call H5Fopen_f (name, H5F_ACC_RDWR_F, fid, error)
  else
    call H5Fcreate_f (name, H5F_ACC_EXCL_F, fid, error)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (verbose > 0) &
    print *, 'h5_open', exists, error
  end function h5_open

  !-----------------------------------------------------------------------------
  ! Close an HDF5 file
  !-----------------------------------------------------------------------------
  subroutine h5_close (id)
  integer(HID_T):: id
  integer:: error
  !.............................................................................
  call H5Fclose_f (id, error)
  call H5close_f (error)
  if (verbose > 0) &
    print *, 'h5_close', error
  end subroutine h5_close

  !-----------------------------------------------------------------------------
  ! Check for group existence
  !-----------------------------------------------------------------------------
  logical function h5_group_exists (id, name)
  integer(HID_T):: id
  character(len=*):: name
  integer:: type, links, corder, error
  !.............................................................................
  call H5Gget_info_by_name_f (id, name, type, links, corder, error)
  h5_group_exists = (error == 0)
  if (verbose > 0) &
    print *, 'h5_group_exists ', name, h5_group_exists, error
  end function h5_group_exists

  !-----------------------------------------------------------------------------
  ! Open group
  !-----------------------------------------------------------------------------
  integer(HID_T) function h5_group_open (id, name, new) result (gid)
  integer(HID_T):: id
  character(len=*):: name
  logical, optional:: new
  integer:: type, links, corder, error
  logical:: exists
  !.............................................................................
  call H5Lexists_f (id, name, exists, error)
  if (exists) then
    call H5Gopen_f (id, name, gid, error)
  else
    call H5Gcreate_f (id, name, gid, error)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (verbose > 0) &
    print *, 'h5_group_open ', name, new, error
  end function h5_group_open

  !-----------------------------------------------------------------------------
  ! Open or create dataset
  !-----------------------------------------------------------------------------
  integer(HID_T) function h5_set_open (id, name, a, new) result (did)
  integer(HID_T):: id, sid
  character(len=*):: name
  real:: a(:,:,:)
  integer:: type, links, corder, error
  logical, optional:: new
  logical:: exists
  !.............................................................................
  call H5Lexists_f (id, name, exists, error)
  if (exists) then
    call H5Dopen_f (id, name, did, error)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, error)
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, error)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (verbose > 0) &
    print *, 'h5_set_open ', name, exists, error
  end function h5_set_open

  !-----------------------------------------------------------------------------
  ! Write a data set
  !-----------------------------------------------------------------------------
  subroutine h5_set_write (id, a)
  integer(HID_T):: id
  real:: a(:,:,:)
  integer:: error
  !.............................................................................
  call H5Dwrite_f (id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), error)
  if (verbose > 0) &
    print *, 'h5_set_write', SHAPE(a), error
  end subroutine h5_set_write

  !-----------------------------------------------------------------------------
  ! Close a data set
  !-----------------------------------------------------------------------------
  subroutine h5_set_close (id)
  integer(HID_T):: id
  integer:: error
  !.............................................................................
  call H5Dclose_f (id, error)
  if (verbose > 0) &
    print *, 'h5_set_close', error
  end subroutine h5_set_close

END PROGRAM test
