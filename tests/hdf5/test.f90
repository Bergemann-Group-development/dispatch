!===============================================================================
!> Demo code that creates a file hdf5.dat with a 'dispatch' dataset, or opens it
!> if it exists, and writes an array to it.
!===============================================================================

PROGRAM test
  USE hdf5
  USE h5_mod
  implicit none
  integer(HID_T):: file_id, set_id, dspc_id, grp_id, seq_id
  real:: a0
  integer:: i0
  real, dimension(16):: a1
  real, dimension(16,16):: a2
  real, dimension(16,16,16):: a3
  real, dimension(16,16,16,5):: a4
  integer, dimension(16,16,16):: i3
  character(len=16):: file, arg
  logical:: new
  integer:: verbose=1
  !-----------------------------------------------------------------------------
  a0 = 1.0
  a1 = 2.0
  a2 = 3.0
  a3 = 4.0
  i0 = 1
  i3 = 4
  call getarg(1,arg)
  if (trim(arg) /= '') & 
    read (arg,*) verbose
  file='hdf5.dat'
  file_id = h5%open (file, new, verbose=verbose)
  if (verbose > 0) &
    print *, merge('creating ','updating ',new)//trim(file)
  grp_id = h5%group_open (file_id, 'sequence')
  seq_id = h5%group_open (grp_id, '000000')
  set_id = h5%set_open (seq_id, 'value', a3)
  call h5%set_write (set_id, a3)
  call h5%output ('sequence', 10, 'a0', a0)
  call h5%output ('sequence', 10, 'a1', a1)
  call h5%output ('sequence', 10, 'a2', a2)
  call h5%output ('sequence', 10, 'a3', a3)
  call h5%output ('sequence', 10, 'a4', a4)
  call h5%output ('sequence', 10, 'i0', i0)
  call h5%output ('sequence', 10, 'i3', i3)
  call h5%set_close (set_id)
  call h5%close (file_id)
END PROGRAM test
