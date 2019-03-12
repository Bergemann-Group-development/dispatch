!===============================================================================
!
! Simple algorithm to do dynamic work stealing between threads
!
!===============================================================================
MODULE mpi_mod
  USE mpi
  implicit none
  logical :: master
  integer :: mpi_mode, err, rank, ncpu
  integer, dimension(MPI_STATUS_SIZE) :: status
END MODULE mpi_mod
!
MODULE work_stealing
  USE mpi_mod
  implicit none
  private
  type thief_t
    integer :: nworker                           ! active workers
    integer :: nresult                           ! outstanding results to receive
    logical, allocatable, dimension(:) :: worker ! is rank still working?
    integer, allocatable, dimension(:) :: req    ! arrays of MPI requests
    integer, allocatable, dimension(:) :: data   ! data package to receive
  end type
  !
  type(thief_t) :: t
  !
  public :: t, thief_present, init_stealing, steal_work
CONTAINS
  SUBROUTINE init_stealing()
    implicit none
    integer :: i
    !
    ! initialize thief structure
    !
    allocate(t%worker(0:ncpu-1),t%req(0:ncpu-1),t%data(0:ncpu-1))
    t%nworker = ncpu-1
    t%nresult = 0
    t%worker  = .true.
    t%worker(rank) = .false.
    t%req = MPI_REQUEST_NULL
    !
    ! Post receives for receiving data packages
    !
    do i=0,ncpu-1
      if (i==rank) cycle
      call MPI_IRECV(t%data(i),1,MPI_INTEGER,i,i+1,mpi_comm_world,t%req(i),err)
    enddo
    !
  END SUBROUTINE init_stealing
  !
  LOGICAL FUNCTION thief_present(cycles)
    implicit none
    integer, intent(in) :: cycles
    integer :: source, tag, dummy
    logical :: retry
    !
    integer :: result ! This is the return payload
    !
    ! Probe to see if we are being approached by a thief
    !
    retry = .true.
    do while (retry)
      call MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,mpi_comm_world,&
                      thief_present,status,err)
      !
      retry = .false.      ! Default is to only probe once
      if (thief_present) then
        source = status(MPI_SOURCE)
        tag    = status(MPI_TAG)
        if (tag==0) then ! A zero tag means that thief needs work
          call MPI_RECV(dummy,0,MPI_INTEGER,source,tag,mpi_comm_world,MPI_STATUS_IGNORE,err) ! drop message from queue
          call MPI_RSEND(cycles,1,MPI_INTEGER,source,rank+1,mpi_comm_world,err)
          t%nresult = t%nresult + 1 ! One more outstanding result to receive
        else              ! A positive tag means that thief is returning result
          call MPI_RECV(result,1,MPI_INTEGER,source,tag,mpi_comm_world,status,err)
          t%nresult = t%nresult - 1
          ! tag should contain the info to insert result in global data structure 
          retry = .true.
        endif
      endif
    enddo
  END FUNCTION thief_present
  !
  SUBROUTINE steal_work(do_work)
    implicit none
    interface
      subroutine do_work(data)
        integer,intent(in) :: data
      end subroutine
    end interface
    !
    integer :: irank, cycles, result, rtag, source, tag, nthief
    integer, dimension(:), allocatable :: preq, rreq
    logical, dimension(:), allocatable :: reqsend
    logical :: flag, retry
    !
    allocate(preq(0:ncpu-1),rreq(0:ncpu-1),reqsend(0:ncpu))
    preq = MPI_REQUEST_NULL
    rreq = MPI_REQUEST_NULL
    reqsend = .false.
    cycles=0
    nthief=1
    !
    ! Send everybody a message telling them that we solicit work
    ! FIXME -- this should be only a subset and gradually increased
    do irank=0,ncpu-1
      if (irank==rank) cycle
      reqsend(irank) = .true.
      call MPI_ISEND(cycles,0,MPI_INTEGER,irank,0,mpi_comm_world,preq(irank),err)
    enddo
    !
    ! Check if somebody wants to send us work or results
    do while(t%nworker > 0 .or. t%nresult > 0 .or. nthief < ncpu)

      ! Test if a data package is being send by a worker:
      flag = .false.
      if (t%nworker > 0) call MPI_TESTANY(ncpu,t%req,irank,flag,status,err)

      retry = flag
      if (flag) then
        irank = irank - 1 ! offset because arrays start at 0, but MPI makes irank 1-based
        if (t%data(irank)<0) then ! sender responded that it has no work left
          t%worker(irank) = .false.
          t%nworker = t%nworker - 1
        else
          call do_work(t%data(irank)) ! process work

          result=0    ! result payload
          rtag=ncpu+1 ! tag telling where to insert result

          call MPI_WAIT(rreq(irank),MPI_STATUS_IGNORE,err) ! Make sure that there are no outstanding results send to this rank
          call MPI_ISEND(result,1,MPI_INTEGER,irank,rtag,mpi_comm_world,rreq(irank),err)

          ! Post new mpi_recv to receive more work
          call MPI_IRECV(t%data(irank),1,MPI_INTEGER,irank,irank+1,mpi_comm_world,t%req(irank),err)

          ! Post new send to signal we want to steal
          call MPI_ISEND(cycles,0,MPI_INTEGER,irank,0,mpi_comm_world,preq(irank),err)
        endif
      endif

      ! Test if a thief is present, or we are receiving a result
      call MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,mpi_comm_world,flag,status,err)
      retry = retry .or. flag

      if (flag) then
        source = status(MPI_SOURCE)
        tag    = status(MPI_TAG)
        if (tag==0) then      ! A zero tag means that thief needs work
          nthief = nthief + 1 ! inc thief counter
          call MPI_RECV(result,0,MPI_INTEGER,source,tag,mpi_comm_world,MPI_STATUS_IGNORE,err)
          ! Tell thief there is no work here
          call MPI_RSEND(-1,1,MPI_INTEGER,source,rank+1,mpi_comm_world,err)
        else if (tag > ncpu) then           ! A positive tag means that worker is returning result
          call MPI_RECV(result,1,MPI_INTEGER,source,tag,mpi_comm_world,MPI_STATUS_IGNORE,err)
          t%nresult = t%nresult - 1
          ! tag should contain the info to insert result in global data structure 
        endif
      endif
      
      if (.not. retry) call do_work(10) ! if nothing in incoming queues, cycle for 10 ms
    enddo

    ! Wait for outstanding mpi_sends and mpi_recvs
    call MPI_WAITALL(ncpu,rreq,MPI_STATUSES_IGNORE, err)
    call MPI_WAITALL(ncpu,preq,MPI_STATUSES_IGNORE, err)
    call MPI_WAITALL(ncpu,t%req,MPI_STATUSES_IGNORE, err)

  END SUBROUTINE steal_work
END MODULE
!===============================================================================
!> Ancient Numerical Recipes code.  Will do for now.
!===============================================================================
FUNCTION ran1(seed)
  implicit none
  integer :: seed
  integer :: idum,IA,IM,IQ,IR,NTAB,NDIV
  real    :: ran1,AM,EPS,RNMX
  parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
    NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  integer :: k
  !.............................................................................
  idum = seed
  k =idum/IQ
  idum = IA*(idum-k*IQ)-IR*k
  if (idum<0) idum = idum+IM
  ran1 = min(AM*idum,RNMX)
  seed = idum
END FUNCTION
!
SUBROUTINE do_work(cycles)
  USE mpi_mod, only : master
  implicit none
  integer, intent(in) :: cycles
  !
  integer(kind=8) :: calibrate=0
  integer(kind=8) :: i,k
  integer         :: count1, count2, count_rate
  real            :: time
  !
  if (calibrate==0) then
    k = 0
    call system_clock(count=count1,count_rate=count_rate)
    do calibrate=1,100000_8
      k = calibrate - floor(exp(sin(real(k))*10.))
    enddo
    call system_clock(count=count2)
    time = 1000.*real(count2-count1) / count_rate ! time in milliseconds
    calibrate = ceiling(100000 / time)           ! nr of iterations for 1ms
    if (master) print *, 'Time in milliseconds, calibrate :', time, calibrate
  endif
  !
  k=0*k
  do i=1_8, calibrate*cycles
    k = i - floor(exp(sin(real(k))*10.))
  enddo
  !
END SUBROUTINE do_work
!
PROGRAM main
  USE mpi_mod
  USE work_stealing
  implicit none
  !
  integer :: seed,i,cycles,nwork
  real    :: x
  real, external :: ran1
  interface
    subroutine do_work(data)
      integer,intent(in) :: data
    end subroutine
  end interface
  !
  !  Start up MPI and get number of nodes and our node rank
  call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, mpi_mode, err)
  call MPI_Comm_size (mpi_comm_world, ncpu, err)
  call MPI_Comm_rank (mpi_comm_world, rank, err)
  master = (rank == 0)
  !
  if (master) then
    print 1, '--------------------------------------------------------------------------------'
    print 1, ' MPI initialized, ncpu =', ncpu
    select case (mpi_mode)
    case(MPI_THREAD_MULTIPLE)
      print 1, ' MPI implementation provides support for simultaneous MPI calls by different threads'
    case(MPI_THREAD_SINGLE)
      print 1, ' MPI implementation provides only support for MPI calls outside openmp regions'
    case(MPI_THREAD_FUNNELED)
      print 1, ' MPI implementation provides only support for MPI calls in master regions'
    case(MPI_THREAD_SERIALIZED)
      print 1, ' MPI implementation provides only support for serial MPI calls by any thread'
    case default
      print 1, ' unknown MPI implementation'
    end select
    print 1, '--------------------------------------------------------------------------------'
    1 format(a,i7)
  end if
  !
  ! Setup random generator by burn-in
  !
  seed = -1*rank
  do i=1,10000
    x = ran1(seed)
  enddo
  !
  ! do work
  !
  nwork = 20 + floor(50 * ran1(seed))
  call init_stealing()
  do i=1,nwork
    !
    cycles = 100 + floor(50 * ran1(seed))
    if (.not. thief_present(cycles)) then
      call do_work(cycles)
    endif
  enddo
  !
  ! The steal work takes as argument a function that does work.
  ! It returns when there is no more work to do.
  call steal_work(do_work)
  !  
END PROGRAM
