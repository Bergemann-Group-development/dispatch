!> =============================================================================
!> Skeleton external code.
!> =============================================================================
MODULE code
  USE interface_mod
  implicit none
  real, pointer:: d(:,:,:)
  real, pointer:: e(:,:,:)
  real, pointer:: p(:,:,:,:)
  logical:: master=.true.
  integer, target:: id
  integer, target:: log_unit
  integer, target:: mpi_rank
  real(8), target:: t, dt
  real(8), target:: pos(3)
  real,    target:: umax
  !$omp threadprivate (d,e,p,t,dt,pos)
CONTAINS

!> =============================================================================
!> Start up code
!> =============================================================================
SUBROUTINE code_begin
  if (master) print *, 'starting up the external code'
END SUBROUTINE code_begin

!> =============================================================================
!> Connect the code variables to the interface arrays and scalars
!> =============================================================================
SUBROUTINE code_connect
  d       => interface%mem(:,:,:,1)
  e       => interface%mem(:,:,:,2)
  p       => interface%mem(:,:,:,3:5)
  interface%id       => id
  interface%time     => t
  interface%dtime    => dt
  interface%u_max    => umax            ! for Courant condition!
  interface%position => pos
  interface%mpi_rank => mpi_rank
  interface%log_unit => log_unit
  master   = (mpi_rank==0)
END SUBROUTINE code_connect

!> =============================================================================
!> Start up code
!> =============================================================================
SUBROUTINE code_init
  call code_connect
  !call whatever_init
END SUBROUTINE code_init

!> =============================================================================
!> Take a time step, after first connecting to the patch data
!> =============================================================================
SUBROUTINE code_update
  call code_connect
  !call whatever_update
  t = t+dt
END SUBROUTINE code_update

!> =============================================================================
!> A chance to do I/O.  May be left empty if generic patch info is enough.
!> =============================================================================
SUBROUTINE code_output
END SUBROUTINE code_output

!> =============================================================================
!> Print a line to standard out and a thread specific log file.  If generic
!> DISPATCH info is enough, remove the experiment_t%output procedure
!> =============================================================================
SUBROUTINE code_info
  print *, id, 'code time: ', t
  write (log_unit, *) id, 'code time: ', t
END SUBROUTINE code_info

!> =============================================================================
!> Pack code info to send to virtual neighbors. May be left empty if generic patch
!> data is sufficient
!> =============================================================================
SUBROUTINE code_pack
  call code_connect
END SUBROUTINE code_pack

!> =============================================================================
!> Unpack code infor from virtual neighbors.  May be left empty if generic patch
!> data is sufficient
!> =============================================================================
SUBROUTINE code_unpack
  call code_connect
END SUBROUTINE code_unpack

!> =============================================================================
!> Close down code
!> =============================================================================
SUBROUTINE code_end
END SUBROUTINE code_end

END MODULE code
