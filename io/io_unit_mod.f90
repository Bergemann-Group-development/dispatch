!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Source of authoritative I/O unit number information
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE io_unit_mod
  implicit none
  private
  type io_unit_t
    logical:: master=.true., do_validate=.false., &
      do_trace=.false., omp_trace=.false., flush_logs=.false.
    integer:: verbose=0, input=1, output=6, data=3, log=4, debug=7, &
      flag=8, index=9, dump=10, trace=11, os=12, mpi=13, queue=14, direct=15, &
      validate=16, dispatcher=17, nml=18, copy1=19, copy2=20, datain=21, &
      nml1=22, nml2=23, dbg=24, hash=25, sinks=26, task=27, lock1=28, lock2=29, &
      aux=30, particle_d=31, particle_i=32, tmp=50, sent=150, &
      mesa=200 ! Mesa has reserved numbers >= 200
    integer:: iodir = -1
    character(len=64):: inputname='input.nml', data_dir='data/', top='../../'
    character(len=64):: run_dir='data/', inputdir='data/'
    character(len=64):: rankbase='rank', threadbase='thread'
  contains
    procedure, nopass:: print_id
  end type
  integer, parameter, public:: mch=120
  type(io_unit_t), public:: io_unit
  integer, public:: stdout=6, stderr=6, stdin=1, stdlog=4
  !$omp threadprivate (io_unit,stdlog)
CONTAINS

!===============================================================================
!> Write information from Id string to data/[run/]hash.log
!===============================================================================
SUBROUTINE print_id (id)
  character(len=120) id
  !.............................................................................
  if (id .ne. '') then
    !$omp critical (print_id_cr)
    if (len_trim(id) > 8) then
      write (io_unit%hash,'(a)') id(6:len_trim(id)-2)
      flush (io_unit%hash)
      id = ''
    end if
    !$omp end critical (print_id_cr)
  end if
END SUBROUTINE print_id

END MODULE io_unit_mod
