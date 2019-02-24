!> =============================================================================
!> Demo code that creates a file hdf5.dat with a 'dispatch' dataset, or opens it
!> if it exists, and writes an array to it.
!> =============================================================================

MODULE h5_mod
  USE hdf5
  implicit none
  private
  type, public:: h5_t
    integer:: err
  contains
    procedure:: open
    procedure:: close
    procedure:: group_open
    procedure:: set_open
    procedure:: set_write
    procedure:: set_close
  end type
  integer:: verbose=0
  type(h5_t), public:: h5
CONTAINS

  !-----------------------------------------------------------------------------
  ! Interface FUNCTION to open or create an HDF5 file
  !-----------------------------------------------------------------------------
  integer(HID_T) FUNCTION open (self, name, new) result (fid)
  class (h5_t):: self
  character(len=*):: name
  logical, optional:: new
  logical:: exists
  !.............................................................................
  call H5open_f (self%err)
  inquire (file=name, exist=exists)
  if (exists) then
    call H5Fopen_f (name, H5F_ACC_RDWR_F, fid, self%err)
  else
    call H5Fcreate_f (name, H5F_ACC_EXCL_F, fid, self%err)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (verbose > 0) &
    print *, 'open', exists, self%err
  END FUNCTION open

  !-----------------------------------------------------------------------------
  ! Close an HDF5 file
  !-----------------------------------------------------------------------------
  SUBROUTINE close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Fclose_f (id, self%err)
  call H5close_f (self%err)
  if (verbose > 0) &
    print *, 'close', self%err
  END SUBROUTINE close

  !-----------------------------------------------------------------------------
  ! Open group
  !-----------------------------------------------------------------------------
  integer(HID_T) FUNCTION group_open (self, id, name, new) result (gid)
  class (h5_t):: self
  integer(HID_T):: id
  character(len=*):: name
  logical, optional:: new
  integer:: type, links, corder
  logical:: exists
  !.............................................................................
  call H5Lexists_f (id, name, exists, self%err)
  if (exists) then
    call H5Gopen_f (id, name, gid, self%err)
  else
    call H5Gcreate_f (id, name, gid, self%err)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (verbose > 0) &
    print *, 'group_open ', name, new, self%err
  END FUNCTION group_open

  !-----------------------------------------------------------------------------
  ! Open or create dataset
  !-----------------------------------------------------------------------------
  integer(HID_T) FUNCTION set_open (self, id, name, a, new) result (did)
  class (h5_t):: self
  integer(HID_T):: id, sid
  character(len=*):: name
  real:: a(:,:,:)
  integer:: type, links, corder
  logical, optional:: new
  logical:: exists
  !.............................................................................
  call H5Lexists_f (id, name, exists, self%err)
  if (exists) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (verbose > 0) &
    print *, 'set_open ', name, exists, self%err
  END FUNCTION set_open

  !-----------------------------------------------------------------------------
  ! Write a data set
  !-----------------------------------------------------------------------------
  SUBROUTINE set_write (self, id, a)
  class (h5_t):: self
  integer(HID_T):: id
  real:: a(:,:,:)
  !.............................................................................
  call H5Dwrite_f (id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  if (verbose > 0) &
    print *, 'set_write', SHAPE(a), self%err
  END SUBROUTINE set_write

  !-----------------------------------------------------------------------------
  ! Close a data set
  !-----------------------------------------------------------------------------
  SUBROUTINE set_close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Dclose_f (id, self%err)
  if (verbose > 0) &
    print *, 'set_close', self%err
  END SUBROUTINE set_close

END MODULE h5_mod

PROGRAM test
  USE hdf5
  USE h5_mod
  implicit none
  integer(HID_T):: file_id, set_id, dspc_id, grp_id, seq_id
  real, dimension(16,16,16):: a
  character(len=16):: file
  logical:: new
  integer:: verbose=1
  !-----------------------------------------------------------------------------
  file='hdf5.dat'
  file_id = h5%open (file, new)
  if (verbose >= 0) &
    print *, merge('creating ','updating ',new)//trim(file)
  grp_id = h5%group_open (file_id, 'sequence')
  seq_id = h5%group_open (grp_id, '000000')
  set_id = h5%set_open (seq_id, 'value', a)
  a = 1.0
  call h5%set_write (set_id, a)
  call h5%set_close (set_id)
  call h5%close (file_id)

END PROGRAM test
