!===============================================================================
!> Module interface to the HDF5 library
!===============================================================================
MODULE h5_mod
#ifdef HDF5
  USE hdf5
  USE io_unit_mod
  USE io_mod
  USE omp_lock_mod
  USE trace_mod
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
    procedure:: flush
    procedure:: group_open
    procedure:: set_open
    procedure:: set_open_ext_1d
    procedure:: set_write
    procedure:: set_write_ext_1d
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
  integer(kind=HSIZE_T), parameter :: chunk_size=16 ! number of records to extend dataspace with in an extendible dataset
  logical, private :: hdf5_initialized=.false.
  type(h5_t), public:: h5
CONTAINS

!===============================================================================
!> Interface FUNCTION to open or create an HDF5 file
!===============================================================================
INTEGER(HID_T) FUNCTION open (self, name, new, verbose) result (fid)
  class (h5_t):: self
  character(len=*):: name
  logical, optional:: new
  integer, optional:: verbose
  logical:: exists
  !.............................................................................
  if (.not. hdf5_initialized) then
    !$omp critical
    if (.not. hdf5_initialized) then
      call H5open_f (self%err)
      hdf5_initialized=.true.
    endif
    !$omp end critical
  endif
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
    write (stdout,*) 'open ', trim(name), exists, self%err
END FUNCTION open

!===============================================================================
!>  Close an HDF5 file
!===============================================================================
SUBROUTINE close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Fclose_f (id, self%err)
  !call H5close_f (self%err) ! never close hdf5 interface in program
  if (self%verbose > 1) &
    write (stdout,*) 'close', self%err
END SUBROUTINE close

!===============================================================================
!>  Flush an HDF5 file
!===============================================================================
SUBROUTINE flush (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  ! scope can be either local or global depending on if the file
  ! has virtual links to other files. Assume it is not the case.
  call H5Fflush_f (id, H5F_SCOPE_LOCAL_F, self%err)
  if (self%verbose > 1) &
    write (stdout,*) 'flush', self%err
END SUBROUTINE flush

!===============================================================================
!>  Open group
!===============================================================================
INTEGER(HID_T) FUNCTION group_open (self, id, name, new) result (gid)
  class (h5_t):: self
  integer(HID_T):: id
  character(len=*):: name
  logical, optional:: new
  integer:: type, links
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
    write (stdout,*) 'group_open ', name, new, self%err
END FUNCTION group_open

!===============================================================================
!> Open or create extensible 1d dataset. The number of records used is stored 
!> in an attribute.  When opening, modify dataset creation properties to enable 
!> chunking, which is required for variable length dataset
!===============================================================================
INTEGER(HID_T) FUNCTION set_open_ext_1d (self, id, name, a, new) result (did)
  USE ISO_C_BINDING
  class (h5_t):: self
  integer(HID_T):: id ! group/file id
  character(len=*):: name
  real(kind=8), dimension(:), intent(in) :: a
  logical, optional:: new
  !.............................................................................
  integer(HID_T):: sid, pid           ! dataspace id, property identifier id
  integer(HSIZE_T) :: maxdims(1), maxdims2(2)
  integer:: type, links, nrec
  integer(HSIZE_T):: asize=1
  logical:: exists
  !-----------------------------------------------------------------------------
  call H5Lexists_f (id, name, exists, self%err)
  if (exists) then
    call H5Dopen_f (id, name, did, self%err)
  else
    !---------------------------------------------------------------------------
    ! if a has one element dataspace is a vector ...
    !---------------------------------------------------------------------------
    if (size(a)==1) then
      maxdims  =(/H5S_UNLIMITED_F/)
      call H5Screate_simple_f (1, int((/chunk_size/),kind=HSIZE_T), sid, &
                               self%err, maxdims)
      call h5pcreate_f (H5P_DATASET_CREATE_F, pid, self%err)
      call h5pset_chunk_f (pid, 1, (/chunk_size/), self%err)
    !---------------------------------------------------------------------------
    ! ... otherwise it is 2-dimensional
    !---------------------------------------------------------------------------
    else
      maxdims2 =(/2_HSIZE_T, H5S_UNLIMITED_F/)
      call H5Screate_simple_f (2, (/int(size(a),kind=HSIZE_T), chunk_size/), sid, &
                               self%err, maxdims2)
      call h5Pcreate_f(H5P_DATASET_CREATE_F, pid, self%err)
      call h5Pset_chunk_f (pid, 2, (/int(size(a),kind=HSIZE_T), chunk_size/), &
                           self%err)
    endif
    call H5Dcreate_f (id, name, H5T_NATIVE_REAL, sid, did, self%err, pid)
    call h5Pclose_f (pid, self%err)
    call h5Sclose_f (sid, self%err)
    !---------------------------------------------------------------------------
    ! write number of records currently stored in dataset
    nrec=0
    !call H5LTset_attribute_int_f (id, name, 'length', nrec, asize, self%err)
    !call H5LTset_attribute_int_f (id, name, 'nchunk', chunk_size, asize, self%err)
  end if
  if (present(new)) then
    new = .not.exists
  end if
  if (self%verbose > 1) &
    write (stdout,*) 'set_open_ext_1d ', name, exists, self%err
END FUNCTION set_open_ext_1d

!===============================================================================
!>  Open or create dataset
!===============================================================================
INTEGER(HID_T) FUNCTION set_open (self, id, name, a, new) result (did)
  class (h5_t):: self
  integer(HID_T):: id, sid
  character(len=*):: name
  real:: a(:,:,:)
  integer:: type, links
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
    write (stdout,*) 'set_open ', name, exists, self%err
END FUNCTION set_open

!===============================================================================
!>  Write to an extensible 1d data set
!===============================================================================
SUBROUTINE set_write_ext_1d (self, id, a)
  class (h5_t):: self
  integer(HID_T):: id, sid
  real(kind=8), dimension(:), intent(in) :: a
  !.............................................................................
  integer(kind=HSIZE_T) :: sza
  integer(HID_T):: aid
  integer(kind=HSIZE_T), dimension(1) :: offset1, count1
  integer(kind=HSIZE_T), dimension(2) :: offset2, count2
  integer, target::nrec,nchunk
  type(c_ptr) :: f_ptr
  !-----------------------------------------------------------------------------
  sza = size(a)
  ! read atrtibute to know how many records have been written until now
  call H5Aopen_f(id, 'length', aid, self%err)
  f_ptr = c_loc(nrec)
  call H5Aread_f(aid, H5T_NATIVE_INTEGER, f_ptr, self%err)  ! read number of records
  nrec = nrec+1                                             ! increment counter
  call H5Awrite_f(aid, H5T_NATIVE_INTEGER, f_ptr, self%err) ! update number of records
  call H5Aclose_f(aid, self%err)
  !-----------------------------------------------------------------------------
  ! Read atttibute to know how big the the data set buffer is 
  !-----------------------------------------------------------------------------
  call H5Aopen_f(id, 'nchunk', aid, self%err)
  f_ptr = c_loc(nchunk)
  call H5Aread_f(aid, H5T_NATIVE_INTEGER, f_ptr, self%err) ! read number of records
  if (nrec > nchunk) then                                  ! increment extend if needed
    nchunk = nchunk + chunk_size
    call H5Awrite_f(aid, H5T_NATIVE_INTEGER, f_ptr, self%err) ! update number of records
    if (sza==1) then
      call h5dset_extent_f(id, int((/nchunk/),kind=HSIZE_T), self%err)
    else
      call h5dset_extent_f(id, int((/sza, int(nchunk,kind=HSIZE_T)/),kind=HSIZE_T), self%err)
    endif
  endif
  call H5Aclose_f(aid, self%err)
  !-----------------------------------------------------------------------------
  ! select hyperslab at the correct position and create corresponding dataspace
  !-----------------------------------------------------------------------------
  call H5Dget_space_f(id, sid, self%err) ! get a copy of the dataspace of the dataset
  if (sza==1) then
    count1  = (/1/) 
    offset1 = int((/nrec/),kind=HSIZE_T)
    call H5Sselect_hyperslab_f(sid, H5S_SELECT_SET_F, offset1, count1, self%err) 
    ! write data at the right position using hyperslabbed dataspace
    call H5Dwrite_f (id, H5T_NATIVE_DOUBLE, a, count1, self%err, sid)
  else
    count2  = (/sza,1_HSIZE_T/) 
    offset2 = int((/0,nrec/),kind=HSIZE_T)
    call H5Sselect_hyperslab_f(sid, H5S_SELECT_SET_F, offset2, count2, self%err) 
    ! write data at the right position using hyperslabbed dataspace
    call H5Dwrite_f (id, H5T_NATIVE_DOUBLE, a, count2, self%err, sid)
  endif
  if (self%verbose > 1) &
    write (stdout,*) 'set_write_ext_1d', SHAPE(a), self%err
END SUBROUTINE set_write_ext_1d

!===============================================================================
!>  Write a data set
!===============================================================================
SUBROUTINE set_write (self, id, a)
  class (h5_t):: self
  integer(HID_T):: id
  real:: a(:,:,:)
  !.............................................................................
  call H5Dwrite_f (id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  if (self%verbose > 1) &
    write (stdout,*) 'set_write', SHAPE(a), self%err
END SUBROUTINE set_write

!===============================================================================
!>  Close a data set
!===============================================================================
SUBROUTINE set_close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Dclose_f (id, self%err)
  if (self%verbose > 1) &
    write (stdout,*) 'set_close', self%err
END SUBROUTINE set_close

!===============================================================================
!>  Open or create dataset
!===============================================================================
INTEGER(HID_T) FUNCTION att_open (self, id, name, a, new) result (did)
  class (h5_t):: self
  integer(HID_T):: id, sid
  character(len=*):: name
  real:: a(:,:,:)
  integer:: type, links
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
    write (stdout,*) 'att_open ', name, exists, self%err
END FUNCTION att_open

!===============================================================================
!>  Write a data set
!===============================================================================
SUBROUTINE att_write (self, id, a)
  class (h5_t):: self
  integer(HID_T):: id
  real:: a(:,:,:)
  !.............................................................................
  call H5Awrite_f (id, H5T_NATIVE_REAL, a, SHAPE(a,KIND=HSIZE_T), self%err)
  if (self%verbose > 1) &
    write (stdout,*) 'att_write', SHAPE(a), self%err
END SUBROUTINE att_write

!===============================================================================
!>  Close a data set
!===============================================================================
SUBROUTINE att_close (self, id)
  class (h5_t):: self
  integer(HID_T):: id
  !.............................................................................
  call H5Aclose_f (id, self%err)
  if (self%verbose > 1) &
    write (stdout,*) 'att_close', self%err
END SUBROUTINE att_close

!===============================================================================
!>  Check if /sequence/RRRRRR/variable
!===============================================================================
LOGICAL FUNCTION exists (self, sequence, record, name, id, attribute)
  class (h5_t)     :: self
  character(len=*) :: sequence, name
  integer          :: record
  character(len=6) :: srecord
  integer(HID_T)   :: id
  logical, optional:: attribute
  !.............................................................................
  if (self%fid == 0) then
    self%fid = self%open (trim(io_unit%run_dir)//'/hdf5.dat')
  end if
  write (srecord,'(i6.6)') record
  id = self%group_open (self%group_open (self%fid, sequence), srecord)
  if (present(attribute)) then
    call H5Aexists_f (id, name, exists, self%err)
  else
    call H5Lexists_f (id, name, exists, self%err)
  end if
END FUNCTION exists

!===============================================================================
!> Generic procedures for scalar, 1D, 2D, 3D, and 4D integer and real output
!===============================================================================
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
    write (stdout,*) 'h5_t%ints: attribute =', name, a, self%err
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
    write (stdout,*) 'h5_t%real: attribute =', name, a, self%err
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

#endif
END MODULE h5_mod
