!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Source of authoritative I/O unit number information
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE io_unit_mod
  implicit none
  private
  type io_unit_t
    logical:: master=.true., do_validate=.false.
    integer:: verbose=0, input=1, output=6, data=3, log=4, debug=7, &
      flag=8, index=9, dump=10, trace=11, os=12, mpi=13, queue=14, direct=15, &
      validate=16, dispatcher=17, nml=18, copy1=19, copy2=20, datain=21, &
      nml1=22, nml2=23, dbg=24, hash=25, sinks=26, task=27, tmp=50, sent=150
    integer:: iodir = -1
    character(len=64):: inputname, outputname, datadir='data', top
    character(len=64):: rundir, inputdir
    character(len=64):: rankbase, threadbase
  end type
  integer, parameter, public:: mch=120
  type(io_unit_t), public:: io_unit
  integer, public:: stdout=6, stderr=6, stdin=1
  !$omp threadprivate (io_unit)
END MODULE io_unit_mod
