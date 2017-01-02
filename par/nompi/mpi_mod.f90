!> =============================================================================
!> Basic MPI calls; initialization, end, barriers, asserts, wallclock for use
!> when MPI is *NOT* activated.
!>
!> =============================================================================
MODULE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_unit_mod
  implicit none
  integer:: mpi_err             ! MPI error status
  logical:: mpi_master          ! true on MPI master (rank 0)
  logical:: mpi_trace=.false.   ! debug trace option
  logical:: master=.true.       ! true on OMP master thread on MPI master process
  integer:: mpi_comm            ! short hand for mpi_comm_world
  integer:: mpi_size            ! number of MPI processes
  integer:: mpi_rank            ! process number (zero-based)
  integer:: output=2
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  ! Object that returns the most important MPI parameters
  !-----------------------------------------------------------------------------
  type mpi_t
    integer:: size=1, rank=0, comm, verbose=0
    integer:: mode
    logical:: master
    logical:: ok = .false.
  contains
    procedure:: init
    procedure:: barrier => barrier_
    procedure:: print => print_
    procedure:: assert
    procedure:: abort => abort_mpi
    procedure:: end => end_mpi
  end type
  type(mpi_t):: mpi
PRIVATE
PUBLIC mpi_t, mpi
CONTAINS

!> =============================================================================
!> Initialize an object with the most important pieces of information about MPI:
!> mpi%size = the number of MPI ranks
!> mpi%rank = the rank of the current MPI process (0-indexed)
!> mpi%coord = the 3-D position of the rank, if a Cartesian arrangement is used
!> =============================================================================
SUBROUTINE init (self, mpi_dims, dims)
  class(mpi_t):: self
  integer, optional, dimension(3):: mpi_dims, dims
  !.............................................................................
  if (self%ok) return
  if (first_time) then
    first_time = .false.
    call init_mpi (self)
    call omp%init
    if (mpi_rank==0) call detect_cores
  end if
  mpi%size   = 1
  mpi%rank   = 0
  mpi%comm   = 0
  mpi%master = (mpi_rank==0)
  mpi%ok     = .true.
  if (mpi%master) then
    print'(a,i4)', ' n_socket =', omp%nsockets
    print'(a,i4)', '   n_core =', omp%ncores
    print'(a,i4)', ' n_thread =', omp_nthreads
  end if
  !-----------------------------------------------------------------------------
  ! For some reason, this must be done element by element.  Turning it around,
  ! and doing mpi = self (which is syntactically correct), values are screwed
  !-----------------------------------------------------------------------------
  self%size   = mpi%size
  self%rank   = mpi%rank
  self%comm   = mpi%comm
  self%master = mpi%master
  self%ok     = mpi%ok
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE print_ (self)
  class(mpi_t):: self
  if (master) print *, 'mpi%size  =', mpi%size
END SUBROUTINE
!> =============================================================================
!> MPI barrier with optional label and delay -- often necessary to arrange serial
!> printout from all ranks, one at a time
!> =============================================================================
SUBROUTINE barrier_ (self, label, delay)
  class(mpi_t):: self
  character(len=*), optional:: label
  real, optional:: delay
  integer:: i
  real(8):: wc
  !.............................................................................
  if (.not.self%ok) return
  if (present(delay)) then
    wc = wallclock()
    do while (wallclock()-wc < delay)
    end do
  end if
  if (present(label)) then
    call barrier_mpi (self, label)
  else
    call barrier_mpi (self, 'mpi%barrier')
  end if
END SUBROUTINE barrier_

!> =============================================================================
SUBROUTINE init_mpi (self)
  implicit none
  class(mpi_t):: self
  character(len=mch):: id='mpi.f90 $Id$'
  character(len=mch):: filename
  integer:: mpi_provided
!-------------------------------------------------------------------------------
!  Start up MPI and get number of nodes and our node rank
!-------------------------------------------------------------------------------
  print*,'MPI is not enabled!'
END SUBROUTINE init_mpi

!> =============================================================================
SUBROUTINE end_mpi (self, label)
  class(mpi_t):: self
  character(len=*), optional:: label
!..............................................................................
  !call finalize_cart
  if (mpi_rank==0) then
    if (present(label)) then
      print*,'NOMPI finalized with message: "'//trim(label)//'"'
    else
      print*,'NOMPI finalized'
    end if
  end if
  call exit        ! avoid "stop" which may produce one output line per process
END SUBROUTINE

!> =============================================================================
SUBROUTINE abort_mpi (self, reason)
  class(mpi_t):: self
  character(len=*), optional:: reason
!...............................................................................
  if (present(reason)) then
    print*,'nompi_abort: '//trim(reason)
    write(io_unit%log,'(3(a,i6,2x))') 'mpi_abort: thread =',omp%thread,trim(reason)
  else
    print*, 'nompi_abort: generic'
    write(io_unit%log,'(3(a,i6,2x))') 'mpi_abort: thread =',omp%thread
  end if
  flush (io_unit%log)
  call exit
END SUBROUTINE abort_mpi

!> =============================================================================
SUBROUTINE clean_stop
!...............................................................................
  call exit
END SUBROUTINE clean_stop

!> =============================================================================
SUBROUTINE barrier_mpi (self, label)
  class(mpi_t):: self
  character(len=*) label
!...............................................................................
  if (mpi_trace) print *,  mpi_rank, 'barrier_nompi: '//label
  if (self%verbose > 2) then
    write (output,*) mpi_rank, 'barrier_nompi: '//label
    flush (output)
  end if
END SUBROUTINE barrier_mpi

!> =============================================================================
SUBROUTINE assert (self, message, code)
  class(mpi_t):: self
  character(len=*):: message
  integer code
  character (len=16):: codemsg
!..............................................................................
  if (code/=0) then
    write (codemsg,'("code =",i10)') code
    write(*,*)      mpi_rank,' error: '//trim(message)
    if (self%verbose > 1) then
      write(output,*) mpi_rank,' error: '//trim(message)
      flush (output)
    end if
  end if
END SUBROUTINE assert

!> =============================================================================
!> Attempt to determine the actual number of cores. If /proc/cpuinfo exists
!> one has to count the number of sockets, but observing hom many different 
!> "physical id" there are.  The number of cores per socket is in the fild
!> "cpu cores".
!> =============================================================================
SUBROUTINE detect_cores
  character(len=32):: file="/proc/cpuinfo"
  character(len=32):: line
  logical:: exist
  integer:: i, iostat, id, oid, n_socket, n_core
  type:: socket_t
    type(socket_t), pointer:: next => null()
    integer:: id
  end type
  type(socket_t), pointer:: p, o, head
  real(8):: wc
  !.............................................................................
  nullify(head)
  n_socket = 0
  oid = -1
  inquire (file=file, exist=exist)
  if (exist) then
    wc = wallclock()
    open (file=file, unit=io_unit%os, form='formatted', status='old')
    do while (.true.)
      read (io_unit%os,'(a32)',iostat=iostat) line
      if (iostat/=0) then
        exit
      else
        i = index(line,':')
        if (line(1:8)=='cpu core') then
          read (line(i+1:),*) n_core
        else if (line(1:8)=='physical') then
          read (line(i+1:),*) id
          if (new(id)) then
            n_socket = n_socket + 1
          end if
        end if
      end if
    end do
    close (io_unit%os)
    wc = wallclock()-wc
    if (wc > 5) print '(a)','############################################################'
    print '(a,f6.1,a)',     'reading /proc/cpuinfo:', wallclock()-wc, ' sec'
    if (wc > 5) print '(a)','############################################################'
  end if
  if (n_socket>0) n_core = n_core*n_socket
  !$omp parallel
  omp%ncores = n_core
  !$omp end parallel
!-----------------------------------------------------------------------------
! For good measure, although totally insignificant, deallocate the list
!-----------------------------------------------------------------------------
  p => head
  do while (associated(p))
    o => p
    p => p%next
    deallocate (o)
  end do
CONTAINS
  logical function new(id)
    integer:: id
    p => head
    o => p
!---------------------------------------------------------------------------
! If the id already exists, return false
!---------------------------------------------------------------------------
    do while (associated(p))
      o => p
      if (p%id == id) then
        new = .false.
        return
      end if
      p => p%next
    end do
!---------------------------------------------------------------------------
! If the id did not exist, add it, and return true
!---------------------------------------------------------------------------
    allocate (p)
    p%id = id
    if (associated(o)) then
      o%next => p
    else
      head => p
    end if
    new = .true.
  end function
END SUBROUTINE detect_cores

END MODULE mpi_mod
