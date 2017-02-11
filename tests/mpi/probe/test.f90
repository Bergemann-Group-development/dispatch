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
    integer:: n
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
    real(8), dimension(:), pointer:: buffer => null()
  end type
  integer, save:: nbuf=8*8**3
  real(8), dimension(:), pointer, save:: buffer => null()
  integer, save:: nrecv=0
  real(8):: total=0d0
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init
  integer:: rank
  integer:: ierr
  logical:: wait=.true.
  integer:: n=8
  real(8):: delay=1e-2
  namelist /input/ wait, n, delay
  !-----------------------------------------------------------------------------
  !print*,'init'
  open (1, file='input.nml', form='formatted', status='old')
  read (1,input)
  close (1)
  mp%n = n
  mp%wait = wait
  mp%del = delay
  call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, mp%provided, ierr)
  mp%comm = mpi_comm_world
  call MPI_COMM_SIZE (MPI_COMM_WORLD, mp%size, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, mp%rank, ierr)
  nbuf = 8*mp%n**3
  allocate (buffer(nbuf))
  print*,'init', mp%rank, nbuf/1024.
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE delay (s)
  real:: s
  real(8):: wc
  !-----------------------------------------------------------------------------
  wc = omp_get_wtime()
  do while (omp_get_wtime()-wc < s)
  end do
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE end
  integer:: ierr
  !-----------------------------------------------------------------------------
  call mp%check
  print*,'average latency', total/max(1,nrecv), nrecv
  flush (5)
  !call sleep(3)
  !call MPI_ABORT (mp%comm, 0, ierr)
  call MPI_FINALIZE (ierr)
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE send
  integer:: rank, req, ierr
  integer:: stat(MPI_STATUS_SIZE)
  integer, save:: seed=1
  real(8):: r
  !-----------------------------------------------------------------------------
  rank = mp%rank
  do while (rank==mp%rank)
    r = ran(seed)
    rank = floor((mp%size-1e-7)*r)
    if (rank==mp%rank) cycle
  end do
  !$omp critical (snd_cr)
  buffer(1) = omp_get_wtime()
  !print*,'send to', rank, buffer(1)
  call MPI_ISEND (buffer, nbuf, MPI_REAL8, rank, mp%rank, mp%comm, req, ierr)
  call MPI_REQUEST_FREE (req, ierr)
  !$omp end critical (snd_cr)
END SUBROUTINE

SUBROUTINE test_wait (mesg)
  type(mesg_t):: mesg
  integer:: req
  integer:: stat(MPI_STATUS_SIZE)
  integer:: ierr, from, nr
  logical:: flag
  real:: latency
  !-----------------------------------------------------------------------------
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
  if (flag) then
    latency = omp_get_wtime()-mesg%buffer(1)
    !print*, mp%rank, 'recv from', from, ' latency', latency, mesg%buffer(1)
    flush (6)
    total = total + latency
    nrecv = nrecv + 1
    mesg%msg = -1
  end if
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE check
  logical:: flag
  integer:: msg, ierr, rank, req, from
  integer:: stat(MPI_STATUS_SIZE)
  real, dimension(:), pointer:: buffer => null()
  type (mesg_t), dimension(:), pointer:: recv => null()
  real(8):: wc
  !-----------------------------------------------------------------------------
  allocate (recv(0:mp%size-1))
  do rank=0,mp%size-1
    allocate (recv(rank)%buffer(nbuf))
  end do
  call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mp%comm, flag, msg, stat, ierr)
  !print*,'check:',flag
  do while (flag)
    call MPI_GET_COUNT (stat, MPI_REAL8, nbuf, ierr)
    from = stat(MPI_SOURCE)
    recv(from)%msg = msg
    call MPI_IMRECV (recv(from)%buffer, nbuf, MPI_REAL8, msg, req, ierr)
    recv(from)%req = req
    call test_wait (recv(from))
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mp%comm, flag, msg, stat, ierr)
    !print*,'recv from:',from
  end do
  do from=0,mp%size-1
    if (recv(from)%msg < 0) cycle
    call test_wait (recv(from))
  end do
  do rank=0,mp%size-1
    deallocate (recv(rank)%buffer)
  end do
  deallocate (recv)
END SUBROUTINE
END MODULE

!> =============================================================================
!> =============================================================================
PROGRAM test_probe
  USE mpi_mod
  implicit none
  real(8):: wc
  !-----------------------------------------------------------------------------
  call mp%init
  wc = omp_get_wtime()
  !$omp parallel shared(wc)
  do while (omp_get_wtime()-wc < 10.)
    call mp%check
    call mp%delay (mp%del)
    call mp%send
  end do
  !$omp end parallel
  call mp%end
END PROGRAM test_probe
