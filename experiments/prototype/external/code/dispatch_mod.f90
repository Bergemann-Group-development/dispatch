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
  private
  !---------------- code local variables ---------------------------------------
  real,    save:: whatever=1.0
  logical, save:: master=.true.
  logical, save:: add_to_buffer=.false.         ! do we want to add to MPI buffer?
  integer, parameter:: mbuffer=10               ! size of extra buffer content
  !---------------- static shared variables ------------------------------------
  integer, target:: id                          ! patch ID - only for info
  integer, target:: log_unit                    ! thread specific log unit
  integer, target:: input_unit                  ! namelist input unit
  integer, target:: mpi_rank                    ! MPI rank
  public:: id, log_unit, input_unit, mpi_rank
  !$omp threadprivate (id,log_unit,input_unit,mpi_rank)
  !---------------- dynamic scalar ---------------------------------------------
  real(8), target:: t, dt                       ! time and time step
  real(8), target:: pos(3)                      ! patch mid point position
!  real(8), target:: size(3)                     ! patch size
  real,    target:: umax                        ! max wave speed - for Courant cond.
  public:: t, dt, pos, umax
  !$omp threadprivate (t,dt,pos,umax)           ! IMPORTANT: add dynamic vars here
  !---------------- dynamic array pointers -------------------------------------
  real, pointer:: d(:,:,:)
  real, pointer:: e(:,:,:)
  real, pointer:: p(:,:,:,:)
  real, pointer:: b(:,:,:,:)
  real, pointer:: dd_dt(:,:,:)
  real, pointer:: de_dt(:,:,:)
  real, pointer:: dp_dt(:,:,:,:)
  real, pointer:: db_dt(:,:,:,:)
  public:: d, e, p, b, dd_dt, de_dt, dp_dt, db_dt
  !$omp threadprivate (d,e,p,b,dd_dt,de_dt,dp_dt,db_dt)
  public:: code_begin, code_init, code_connect, code_update, code_output, &
           code_pack, code_unpack, code_info, code_end
CONTAINS

!> =============================================================================
!> Start up code
!> =============================================================================
SUBROUTINE code_begin
  master = (mpi_rank==0)
  if (master) print *, 'starting up the external code'
  !call whatever_begin
END SUBROUTINE code_begin

!> =============================================================================
!> Connect the code variables to the interface arrays and scalars.  The code
!> scalar variables that need to be connected must have the 'target' attribute.
!> =============================================================================
SUBROUTINE code_connect
  d       => interface%mem(:,:,:,1  ,1)
  e       => interface%mem(:,:,:,2  ,1)
  p       => interface%mem(:,:,:,3:5,1)
  b       => interface%mem(:,:,:,6:8,1)
  dd_dt   => interface%mem(:,:,:,1  ,2)
  de_dt   => interface%mem(:,:,:,2  ,2)
  dp_dt   => interface%mem(:,:,:,3:5,2)
  db_dt   => interface%mem(:,:,:,6:8,2)
  interface%id         => id
  interface%time       => t
  interface%dtime      => dt
  interface%u_max      => umax
  interface%position   => pos
!  interface%size       => size
  interface%mpi_rank   => mpi_rank
  interface%log_unit   => log_unit
  interface%input_unit => input_unit
END SUBROUTINE code_connect

!> =============================================================================
!> Initialize a patch; this includes generating initial conditions (possible from
!> previous snapshots), reading namelist parameter files (but only once!), etc
!> =============================================================================
SUBROUTINE code_init
  logical, save:: first_time = .true.
  namelist /code_params/ whatever
  !-----------------------------------------------------------------------------
  !$omp critical (code_cr)
  if (first_time) then
    first_time = .false.
    print *, input_unit
    rewind (input_unit)
    read (input_unit, code_params)
    if (master) write (*, code_params)
  end if
  !$omp end critical (code_cr)
  !call whatever_init
END SUBROUTINE code_init

!> =============================================================================
!> Take a time step, after first connecting to the patch data.  This may also
!> include computing the max wave speed in the patch, storing it in a variable
!> that is connected to the DISPATCH patch_t%u_max variable via pointers (only
!> relevant when DISPATCH takes care of time stepping.
!> =============================================================================
SUBROUTINE code_update
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
  if (id==1) print *, id, 'code time: ', t
  write (log_unit, *) id, 'code time: ', t
END SUBROUTINE code_info

!> =============================================================================
!> Pack code info to send to virtual neighbors. May be left empty if generic patch
!> data is sufficient.
!> =============================================================================
SUBROUTINE code_pack (buffer, nbuffer, ibuffer)
  integer, pointer:: buffer(:)
  integer, pointer:: tmp(:)
  integer:: nbuffer, ibuffer, i, mbuffer
  !-----------------------------------------------------------------------------
  if (add_to_buffer) then
    mbuffer = 100
    allocate (tmp(nbuffer+mbuffer))
    tmp(1:nbuffer) = buffer
    do i=1,mbuffer
      tmp(nbuffer+i) = i
    end do
    deallocate (buffer)
    buffer => tmp
    nbuffer = nbuffer + mbuffer
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
      do i=1,mbuffer
        if (master) print *,'buffer:', i, buffer(nbuffer-(mbuffer-i))
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
