!> =============================================================================
!> Test program for MPI communication
!> =============================================================================
MODULE mpi_mod
  USE mpi
  USE omp_lib
  implicit none
  type mpi_t
    integer:: err, comm, size, rank, provided
    logical:: wait=.false.
    integer:: n, max_probe
    real:: del
  contains
    procedure, nopass:: init
    procedure, nopass:: send
    procedure, nopass:: check
    procedure, nopass:: delay
    procedure, nopass:: end
  end type
  type(mpi_t):: mp
  type mesg_t
    integer:: req, msg=-1
    real(4), dimension(:), pointer:: buffer => null()
  end type
  integer, save:: nbuf=8*8**3
  integer, save:: nrecv=0, nsent=0
  real(8):: total=0d0
  real(8):: wc
  !-----------------------------------------------------------------------------
  ! Threadprivate send/recv data
  !-----------------------------------------------------------------------------
  type (mesg_t), dimension(:), pointer:: recv => null()
  real(4), dimension(:), pointer:: buffer => null()
  !$omp threadprivate (buffer,recv)
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init
  integer:: rank
  integer:: ierr
  logical:: wait=.true.
  integer, save:: n=8, max_probe=2
  real(8):: delay=1e-2
  namelist /input/ wait, n, delay, max_probe
  !-----------------------------------------------------------------------------
  !print*,'init'
  open (1, file='input.nml', form='formatted', status='old')
  read (1,input)
  close (1)
  mp%n = n
  mp%wait = wait
  mp%del = delay
  mp%max_probe = max_probe
  call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, mp%provided, ierr)
  mp%comm = mpi_comm_world
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mp%size, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, mp%rank, ierr)
  nbuf = 8*mp%n**3
  if (mp%rank==0) then
    write (*,input)
    select case (mp%provided)
      case (MPI_THREAD_MULTIPLE)
        print *, ' MPI implementation provides support for simultaneous MPI calls by different threads'
      case(MPI_THREAD_SERIALIZED)
        print *, ' MPI implementation provides only support for serial MPI calls by any thread'
      case (MPI_THREAD_SINGLE)
        print *, ' MPI implementation provides only support for MPI calls outside openmp regions'
      case default
        print *, mp%provided, MPI_THREAD_MULTIPLE, MPI_THREAD_SERIALIZED
    end select
    print *, 'buffer size:', nbuf*sizeof(buffer(1))/1024., ' kb'
  end if
  !$omp parallel
  allocate (buffer(nbuf))
  allocate (recv(0:mp%size-1))
  do rank=0,mp%size-1
    allocate (recv(rank)%buffer(nbuf))
  end do
  !$omp end parallel
END SUBROUTINE init

!> =============================================================================
!> Delay for a certain time, emulating the patch update time
!> =============================================================================
SUBROUTINE delay (s)
  real:: s
  real(8):: wc
  !-----------------------------------------------------------------------------
  wc = omp_get_wtime()
  do while (omp_get_wtime()-wc < s)
  end do
END SUBROUTINE delay

!> =============================================================================
!> Wrap up MPI, finalizing first all pending MPI receives
!> =============================================================================
SUBROUTINE end
  integer:: rank, ierr, i
  !-----------------------------------------------------------------------------
  mp%max_probe=1000
  do i=1,10
    call mp%check
    call mp%delay (mp%del)
  end do
  print*,'average latency', total/max(1,nrecv), nsent, nrecv
  do rank=0,mp%size-1
    deallocate (recv(rank)%buffer)
  end do
  deallocate (recv)
  call MPI_FINALIZE (ierr)
END SUBROUTINE end

!> =============================================================================
!> Send a buffer to a random other rank, freeing the request, since we do not
!> bother about waiting on the request in this test
!> =============================================================================
SUBROUTINE send
  integer:: rank, req, ierr
  integer:: stat(MPI_STATUS_SIZE)
  integer, save:: seed=1
  real(8):: r
  !-----------------------------------------------------------------------------
  rank = mp%rank
  !do while (rank==mp%rank)
    !r = ran()
    !rank = floor((mp%size-1e-7)*r)
    !if (rank==mp%rank) cycle
  !end do
  rank = mod(mp%rank+1,mp%size)
  buffer(1) = omp_get_wtime()-wc
  if (mp%provided==MPI_THREAD_MULTIPLE) then
    call doit
  else
    !$omp critical (mpi_cr)
    call doit
    !$omp end critical (mpi_cr)
  end if
contains
  subroutine doit
    !print *, 'sending to', rank
    call MPI_ISEND (buffer, nbuf, MPI_REAL, rank, mp%rank, mp%comm, req, ierr)
    call MPI_REQUEST_FREE (req, ierr)
    !$omp atomic
    nsent = nsent+1
  end subroutine doit
END SUBROUTINE send

!> =============================================================================
!> Wait for recv completion, either with MPI_WAIT or repeated MPI_TEST
!> =============================================================================
SUBROUTINE test_wait (from)
  integer:: from
  integer:: req
  integer:: stat(MPI_STATUS_SIZE)
  integer:: ierr, nr
  logical:: flag
  real:: latency
  type(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  mesg => recv(from)
  !print *, 'test_wait begin', from, mesg%req
  call doit
  if (flag) then
    latency = omp_get_wtime()-wc-recv(from)%buffer(1)
    !print*, mp%rank, 'recv from', from, ' latency', latency, mesg%buffer(1)
    !flush (6)
    !$omp atomic
    total = total + latency
    !$omp atomic
    nrecv = nrecv + 1
    mesg%msg = -1
  end if
  !print *, 'test_wait end'
contains
  subroutine doit
  if (mp%wait) then
    call MPI_WAIT (mesg%req, stat, ierr)
    flag = .true.
  else
    flag = .false.
    nr = 0
    do while (.not.flag)
      call MPI_TEST (mesg%req, flag, stat, ierr)
      nr = nr+1
    end do
    !print*, 'test',mesg%req, nr
  end if
  end subroutine doit
END SUBROUTINE test_wait

!> =============================================================================
!> Probe for incoming MPI messages
!> =============================================================================
SUBROUTINE check
  logical:: flag
  integer:: msg, ierr, rank, req, from, n
  integer:: stat(MPI_STATUS_SIZE)
  real, dimension(:), pointer:: buffer => null()
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (mp%provided==MPI_THREAD_MULTIPLE) then
    call doit
  else
    !$omp critical (mpi_cr)
    call doit
    !$omp end critical (mpi_cr)
  end if
contains
  subroutine doit
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mp%comm, flag, msg, stat, ierr)
    n = 0
    !print *, 'flag', flag
    do while (flag)
      n = n+1
      call MPI_GET_COUNT (stat, MPI_REAL, nbuf, ierr)
      from = stat(MPI_SOURCE)
      recv(from)%msg = msg
      call MPI_IMRECV (recv(from)%buffer, nbuf, MPI_REAL, msg, req, ierr)
      recv(from)%req = req
      !print *, 'from', from, req
      call test_wait (from)
      if (n > mp%max_probe) exit
      call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mp%comm, flag, msg, stat, ierr)
      !print *, 'flag', flag
    end do
  end subroutine doit
END SUBROUTINE check
END MODULE

!> =============================================================================
!> For each loop the code checks for incoming messages and tries to process as
!> many as possible (cf. comments in check procedure).  It then delays for a time
!> that corresponds to a patch update time, and then sends a patch-load of data
!> to a random rank.  Hence each iteration of the loop generates one send, and
!> processes on average one recv.  By tuning the delay time to be larger or 
!> smaller than the typical latency one can study the growth of send and recv
!> queues.
!> =============================================================================
PROGRAM test_probe
  USE mpi_mod
  implicit none
  integer:: ierr
  !-----------------------------------------------------------------------------
  call mp%init
  call MPI_Barrier (mp%comm, ierr)
  wc = omp_get_wtime()
  print *, mp%rank, wc
  !$omp parallel shared(wc)
  do while (omp_get_wtime()-wc < 10.)
    !print *, omp_get_thread_num(), 'check'
    call mp%check
    !print *, omp_get_thread_num(), 'delay'
    call mp%delay (mp%del)
    !print *, omp_get_thread_num(), 'send'
    call mp%send
  end do
  !$omp end parallel
  call mp%end
END PROGRAM test_probe
