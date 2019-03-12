!===============================================================================
!> Module interface to the HDF5 library
!===============================================================================

MODULE h5_mod
  USE hdf5
  USE io_unit_mod
  USE omp_lock_mod
  implicit none
  private
  type, public:: h5_t
    integer:: err
    integer:: verbose=0
    integer(HID_T):: fid=0, gid=0
    type(lock_t):: lock
  contains
    procedure:: open
    procedure:: close
    procedure:: group_open
    procedure:: set_open
    procedure:: set_write
    procedure:: set_close
    procedure:: att_open
    procedure:: att_write
    procedure:: att_close
    procedure:: exists
    procedure, private:: ints
    procedure, private:: ints_1d
    procedure, private:: ints_2d
    procedure, private:: ints_3d
    procedure, private:: ints_4d
    procedure, private:: real
    procedure, private:: real_1d
    procedure, private:: real_2d
    procedure, private:: real_3d
    procedure, private:: real_4d
    generic, public:: output => ints, ints_1d, ints_2d, ints_3d, ints_4d, &
                                real, real_1d, real_2d, real_3d, real_4d
  end type
  type(h5_t), public:: h5
CONTAINS

!> -----------------------------------------------------------------------------
!> Interface FUNCTION to open or create an HDF5 file
!> -----------------------------------------------------------------------------
INTEGER(HID_T) FUNCTION open (self, name, new, verbose) result (fid)
  class (h5_t):: self
  character(len=*):: name
  logical, optional:: new
  integer, optional:: verbose
  logical:: exists
  !.............................................................................
  call H5open_f (self%err)
  inquire (file=name, exist=exists)
  if (exists) then
    call H5Fopen_f (name, H5F_ACC_RDWR_F, fid, self%err)
  else
    call H5Fcreate_f (name, H5F_ACC_EXCL_F, fid, self%err)
  end if
  self%fid = fid
  if (present(new)) then
    new = .not.exists
  end if
  if (present(verbose)) &
    self%verbose = verbose
  if (self%verbose > 1) &
    print *, 'open ', trim(name), exists, self%err
END FUNCTION open

!> -----------------------------------------------------------------------------
!>  Close an HDF5 file
!> -----------------------------------------------------------------------------
SUBROUTINE close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Fclose_f (id, self%err)
  call H5close_f (self%err)
  if (self%verbose > 1) &
    print *, 'close', self%err
END SUBROUTINE close

!> -----------------------------------------------------------------------------
!>  Open group
!> -----------------------------------------------------------------------------
INTEGER(HID_T) FUNCTION group_open (self, id, name, new) result (gid)
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
  self%gid = gid
  if (present(new)) then
    new = .not.exists
  end if
  if (self%verbose > 1) &
    print *, 'group_open ', name, new, self%err
END FUNCTION group_open

!> -----------------------------------------------------------------------------
!>  Open or create dataset
!> -----------------------------------------------------------------------------
INTEGER(HID_T) FUNCTION set_open (self, id, name, a, new) result (did)
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
  if (self%verbose > 1) &
    print *, 'set_open ', name, exists, self%err
END FUNCTION set_open

!> -----------------------------------------------------------------------------
!>  Write a data set
!> -----------------------------------------------------------------------------
SUBROUTINE set_write (self, id, a)
  class (h5_t):: self
  integer(HID_T):: id
  real:: a(:,:,:)
  !.............................................................................
  call H5Dwrite_f (id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  if (self%verbose > 1) &
    print *, 'set_write', SHAPE(a), self%err
END SUBROUTINE set_write

!> -----------------------------------------------------------------------------
!>  Close a data set
!> -----------------------------------------------------------------------------
SUBROUTINE set_close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Dclose_f (id, self%err)
  if (self%verbose > 1) &
    print *, 'set_close', self%err
END SUBROUTINE set_close

!> -----------------------------------------------------------------------------
!>  Open or create dataset
!> -----------------------------------------------------------------------------
INTEGER(HID_T) FUNCTION att_open (self, id, name, a, new) result (did)
  class (h5_t):: self
  integer(HID_T):: id, sid
  character(len=*):: name
  real:: a(:,:,:)
  integer:: type, links, corder
  logical, optional:: new
  logical:: exists
  !.............................................................................
  call H5Aexists_f (id, name, exists, self%err)
  if (exists) then
    call H5Aopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (1, SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Acreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (self%verbose > 1) &
    print *, 'att_open ', name, exists, self%err
END FUNCTION att_open

!> -----------------------------------------------------------------------------
!>  Write a data set
!> -----------------------------------------------------------------------------
SUBROUTINE att_write (self, id, a)
  class (h5_t):: self
  integer(HID_T):: id
  real:: a(:,:,:)
  !.............................................................................
  call H5Awrite_f (id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  if (self%verbose > 1) &
    print *, 'att_write', SHAPE(a), self%err
END SUBROUTINE att_write

!> -----------------------------------------------------------------------------
!>  Close a data set
!> -----------------------------------------------------------------------------
SUBROUTINE att_close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Aclose_f (id, self%err)
  if (self%verbose > 1) &
    print *, 'att_close', self%err
END SUBROUTINE att_close

!> -----------------------------------------------------------------------------
!>  Check if /sequence/RRRRRR/variable
!> -----------------------------------------------------------------------------
LOGICAL FUNCTION exists (self, sequence, record, name, id, attribute)
  class (h5_t)     :: self
  character(len=*) :: sequence, name
  integer          :: record
  character(len=6) :: srecord
  integer(HID_T)   :: id
  logical, optional:: attribute
  !.............................................................................
  if (self%fid == 0) then
    self%fid = self%open (trim(io_unit%outputname)//'/hdf5.dat')
  end if
  write (srecord,'(i6.6)') record
  id = self%group_open (self%group_open (self%fid, sequence), srecord)
  if (present(attribute)) then
    call H5Aexists_f (id, name, exists, self%err)
  else
    call H5Lexists_f (id, name, exists, self%err)
  end if
END FUNCTION exists

!> -----------------------------------------------------------------------------
SUBROUTINE ints (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  integer          :: a
  integer(HID_T)   :: id, sid, did, rank(1)
  !.............................................................................
  if (self%exists(seq, record, name, id, attribute=.true.)) then
    call H5Aopen_f (id, name, did, self%err)
  else
    rank = 1
    call H5Screate_simple_f (size(rank), rank, sid, self%err)
    call H5Acreate_f (id, name, H5T_NATIVE_INTEGER, sid, did, self%err)
  end if
  call H5Awrite_f (did, H5T_NATIVE_INTEGER, a, rank, self%err)
  if (self%verbose > 1) &
    print *, 'h5_t%ints: attribute =', name, a, self%err
  call H5Aclose_f (did, self%err)
END SUBROUTINE ints
!> -----------------------------------------------------------------------------
SUBROUTINE ints_1d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  integer          :: a(:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_INTEGER, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_INTEGER, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE ints_1d
!> -----------------------------------------------------------------------------
SUBROUTINE ints_2d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  integer          :: a(:,:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_INTEGER, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_INTEGER, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE ints_2d
!> -----------------------------------------------------------------------------
SUBROUTINE ints_3d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  integer          :: a(:,:,:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_INTEGER, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_INTEGER, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE ints_3d
!> -----------------------------------------------------------------------------
SUBROUTINE ints_4d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  integer          :: a(:,:,:,:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_INTEGER, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_INTEGER, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE ints_4d

!> -----------------------------------------------------------------------------
SUBROUTINE real (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  real             :: a
  integer(HID_T)   :: id, sid, did, rank(1)
  !.............................................................................
  if (self%exists(seq, record, name, id, attribute=.true.)) then
    call H5Aopen_f (id, name, did, self%err)
  else
    rank = 1
    call H5Screate_simple_f (size(rank), rank, sid, self%err)
    call H5Acreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  call H5Awrite_f (did, H5T_NATIVE_REAL, a, rank, self%err)
  if (self%verbose > 1) &
    print *, 'h5_t%real: attribute =', name, a, self%err
  call H5Aclose_f (did, self%err)
END SUBROUTINE real
!> -----------------------------------------------------------------------------
SUBROUTINE real_1d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  real             :: a(:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE real_1d
!> -----------------------------------------------------------------------------
SUBROUTINE real_2d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  real             :: a(:,:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE real_2d
!> -----------------------------------------------------------------------------
SUBROUTINE real_3d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  real             :: a(:,:,:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE real_3d
!> -----------------------------------------------------------------------------
SUBROUTINE real_4d (self, seq, record, name, a)
  class (h5_t)     :: self
  character(len=*) :: seq, name
  integer          :: record
  real             :: a(:,:,:,:)
  integer(HID_T)   :: id, sid, did
  !.............................................................................
  if (self%exists(seq, record, name, id)) then
    call H5Dopen_f (id, name, did, self%err)
  else
    call H5Screate_simple_f (SIZE(SHAPE(a)), SHAPE(a,KIND=HSIZE_T), sid, self%err)
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err)
  end if
  call H5Dwrite_f (did, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  call H5Dclose_f (did, self%err)
END SUBROUTINE real_4d

END MODULE h5_mod
