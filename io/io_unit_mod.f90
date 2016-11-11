!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Source of authoritative I/O unit number information
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE io_unit_mod
  implicit none
  private
  type io_unit_t
    logical:: master=.true., do_validate=.false.
    integer:: verbose=0, input=1, output=6, data=3, log=4, debug=7, &
      flag=8, index=9, dump=10, trace=11, os=12, mpi=13
  end type
  integer, parameter, public:: mch=120
  type(io_unit_t), public:: io_unit
  integer, public:: log_unit=6
  !$omp threadprivate (io_unit, log_unit)
END MODULE io_unit_mod
