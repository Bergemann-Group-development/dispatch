!> =============================================================================
!> Skeleton external code.  To make the external code thread safe (able to run
!> simultaneously on any number of OMP threads), all dynamic variables -- that
!> change during a time step -- must be declared "!$omp threadprivate" here.
!>
!> Any variables and parameters that need to be communicated between DISPATCH and
!> the external code (note that most variable can and should be kept hidden from
!> DISPATCH) should be entered in the interface_mod.f90 module, and corresponding
!> pointer assignements should be added in experiment_mod::connect and in 
!> code_connect here (for more info see the procedure headers).
!> =============================================================================
MODULE dispatch_mod
  USE interface_mod
  implicit none
  !---------------- static variables -------------------------------------------
  integer, target:: id                          ! patch ID - only for info
  integer, target:: log_unit                    ! thread specific log unit
  integer, target:: mpi_rank                    ! MPI rank
  logical:: master=.true.
  logical:: add_to_buffer=.false.               ! do we want to add to MPI buffer?
  !---------------- dynamic scalar ---------------------------------------------
  real(8), target:: t, dt                       ! time and time step
  real(8), target:: pos(3)                      ! patch mid point position
  real,    target:: umax                        ! max wave speed - for Courant cond.
  !---------------- dynamicarray pointers --------------------------------------
  real, pointer:: d(:,:,:)
  real, pointer:: e(:,:,:)
  real, pointer:: p(:,:,:,:)
  !$omp threadprivate (d,e,p,t,dt,pos,umax)     ! IMPORTANT: add dynamic vars here
CONTAINS

!> =============================================================================
!> Start up code
!> =============================================================================
SUBROUTINE code_begin
  if (master) print *, 'starting up the external code'
  !call whatever_begin
END SUBROUTINE code_begin

!> =============================================================================
!> Connect the code variables to the interface arrays and scalars.  The code
!> variables that need to be connected must have the 'target' attribute.
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
!> Initialize a patch; this includes generating initial conditions (possible from
!> previous snapshots), reading namelist parameter files (but only once!), etc
!> =============================================================================
SUBROUTINE code_init
  call code_connect
  !call whatever_init
END SUBROUTINE code_init

!> =============================================================================
!> Take a time step, after first connecting to the patch data.  This must also
!> include computing the max wave speed in the patch, storing it in a variable
!> that is connected to the DISPATCH patch_t%u_max variable via pointers.
!> =============================================================================
SUBROUTINE code_update
  call code_connect
  !call whatever_update
  umax = 1.0
  t = t+dt
END SUBROUTINE code_update

!> =============================================================================
!> A chance to do I/O.  May be left empty if generic patch info is enough.
!> =============================================================================
SUBROUTINE code_output (out_time)
  real(8):: out_time
  if (t > out_time) then
    if (master) print *, 'external code may want to write additional output now'
  end if
END SUBROUTINE code_output

!> =============================================================================
!> Print a line to standard out and a thread specific log file.
!> =============================================================================
SUBROUTINE code_info
  print *, id, 'code time: ', t
  write (log_unit, *) id, 'code time: ', t
END SUBROUTINE code_info

!> =============================================================================
!> Pack code info to send to virtual neighbors. May be left empty if generic patch
!> data is sufficient.
!> =============================================================================
SUBROUTINE code_pack (buffer, nbuffer, ibuffer)
  integer, pointer:: buffer(:)
  integer:: nbuffer, ibuffer, i
  !-----------------------------------------------------------------------------
  if (add_to_buffer) then
    nbuffer = 100
    allocate (buffer(nbuffer))
    buffer = [(real(i), i=1,nbuffer)]
  else
    nbuffer = 0
    nullify (buffer)
  end if
END SUBROUTINE code_pack

!> =============================================================================
!> Unpack code info from virtual neighbors.  May be left empty if generic patch
!> data is sufficient.
!> =============================================================================
SUBROUTINE code_unpack (buffer, nbuffer, ibuffer)
  integer, pointer:: buffer(:)
  integer:: nbuffer, ibuffer, i
  !-----------------------------------------------------------------------------
  if (add_to_buffer) then
    if (master) then
      do i=1,nbuffer
        if (master) print *,'buffer:', i, buffer(i)
      end do
    end if
    deallocate (buffer)
  end if
END SUBROUTINE code_unpack

!> =============================================================================
!> Close down code
!> =============================================================================
SUBROUTINE code_end
  !call whatever_shut_down
END SUBROUTINE code_end

END MODULE dispatch_mod
