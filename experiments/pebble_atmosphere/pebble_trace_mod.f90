MODULE pebble_trace_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE omp_lock_mod
  USE pebble_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Data types for tracing of selected pebbles
  !-----------------------------------------------------------------------------
  type:: trace_buffer_t
    integer:: n, record=0
    integer:: id, status
    real:: s
    real, dimension(:,:), pointer:: v
    integer, dimension(:), pointer:: ip
  end type
  type pebble_trace_t
    type(lock_t):: lock
    type(trace_buffer_t), dimension(:), pointer:: buffer
    integer:: npebble=10
    integer:: nbuffer=10000                     ! max number of traces
    integer:: every=1000                        ! save frequency
  contains
    procedure init
    procedure trace_pebble                      ! trace pebbles
    procedure trace_output                      ! dump traced pebbles
    procedure check                             ! check if should close
    procedure close                             ! close one trace
    procedure close_all                         ! close all traces
  end type
  type(pebble_trace_t), public:: pebble_trace
  integer, save:: id=0
CONTAINS

!> =============================================================================
SUBROUTINE init (self, trace_every)
  class(pebble_trace_t):: self
  integer:: trace_every, i
  !-----------------------------------------------------------------------------
  pebble_trace%every = trace_every
  if (trace_every > 0) then
    call pebble_trace%lock%init ('trace')
    allocate(pebble_trace%buffer(pebble_trace%nbuffer))
    do i=1,pebble_trace%nbuffer
      allocate (pebble_trace%buffer(i)%v(14,pebble_trace%npebble))
      allocate (pebble_trace%buffer(i)%ip(pebble_trace%npebble))
      pebble_trace%buffer(i)%n = 0
    end do
  end if
END SUBROUTINE init

!> =============================================================================
!> Add a pebble to the trace buffer, and check for output.  A trace pebble has
!> a status word where the first 16 bits is a consecutive trace numbre, starting
!> with 1, and running up to nbuffer.  The 17th and 18th bit are set to mark the
!> pebble as being traced, and as being launched from the boundary, respectively
!> =============================================================================
SUBROUTINE trace_pebble (self, pebble, ip, level)
  class(pebble_trace_t):: self
  class(pebble_t), pointer:: pebble
  integer:: ip, level
  integer:: i, n
  !-----------------------------------------------------------------------------
  call trace_begin ('pebbles_mod::trace_pebble')
  !-----------------------------------------------------------------------------
  ! Check if the pebble has a trace id, and give it one if not
  !-----------------------------------------------------------------------------
  i = iand(pebble%status,not(bits%trace+bits%has_BC))
  if (i == 0) then
    !$omp atomic
    id = id+1
    pebble%status = pebble%status + id
    i = id
  end if
  if (i <= self%nbuffer) then
    self%buffer(i)%id      = pebble%id
    self%buffer(i)%status  = pebble%status
    self%buffer(i)%s       = pebble%s
    n = self%buffer(i)%n+1
    self%buffer(i)%ip(n) = ip
    self%buffer(i)%v(1:3,n) = pebble%p
    self%buffer(i)%v(4:6,n) = pebble%v
    self%buffer(i)%v(7:9,n) = pebble%vg
    self%buffer(i)%v(10,n)  = pebble%c
    self%buffer(i)%v(11,n)  = pebble%d
    self%buffer(i)%v(12,n)  = pebble%t
    self%buffer(i)%v(13,n)  = pebble%s
    self%buffer(i)%v(14,n)  = pebble%w
    self%buffer(i)%n = n
    if (io%verbose > 1) print *,'MK', pebble%id, i, n, pebble%t, pebble%dt
    if (n >= self%npebble) then
      call self%trace_output (i)
    end if
  else
    if (io%verbose > 1) print *,'MK', pebble%id, i
  end if
  call trace_end
END SUBROUTINE trace_pebble

!> =============================================================================
!> Write out a pebble trace buffer
!> =============================================================================
SUBROUTINE trace_output (self, i, close_key)
  class(pebble_trace_t):: self
  integer:: i, close_key
  integer, parameter:: key=-1
  optional:: close_key
  character(len=64):: filename
  integer:: j, n
  !-----------------------------------------------------------------------------
  if (i>self%nbuffer) return
  if (.not.associated(self%buffer)) return
  if (self%buffer(i)%n==0) return
  call trace_begin ('pebbles_mod::trace_output')
  call self%lock%set ('trace_output')
  write (filename,'(a,i4.4,".trc")') trim(io%outputname), i
  if (io%verbose>0) &
    print*,'pebble_trace: ', trim(filename), self%buffer(i)%record, self%buffer(i)%n
  if (self%buffer(i)%record == 0) then
    open (io_unit%trace, file=trim(filename), form='unformatted', &
      status='unknown')
    close (io_unit%trace, status='delete')
    open (io_unit%trace, file=trim(filename), form='unformatted', &
      status='new')
    write (io_unit%trace) key, self%buffer(i)%status
    write (io_unit%trace) self%buffer(i)%s
  else
    open (io_unit%trace, file=trim(filename), form='unformatted', &
      status='unknown', position = 'append')
    backspace (io_unit%trace)
  end if
  n =  self%buffer(i)%n
  self%buffer(i)%record = self%buffer(i)%record+1
  write (io_unit%trace) self%buffer(i)%id, self%buffer(i)%record, self%buffer(i)%n
  write (io_unit%trace) self%buffer(i)%ip(1:n)
  write (io_unit%trace) self%buffer(i)%v(:,1:n)
  self%buffer(i)%n = 0
  write (io_unit%trace) 0,0,0,0
  close (io_unit%trace)
  call self%lock%unset ('trace_output')
  call trace_end
END SUBROUTINE trace_output

!> =============================================================================
!> Check if this pebble is traced, and if so flush and close its file
!> =============================================================================
SUBROUTINE check (self, pebble)
  class(pebble_trace_t):: self
  class(pebble_t), pointer:: pebble
  integer:: i
  !-----------------------------------------------------------------------------
  if (iand(pebble%status,bits%trace)==0) return
  i = iand(pebble%status,not(bits%trace+bits%has_BC))
  call self%close (i)
END SUBROUTINE check

!> =============================================================================
!> Flush and close trace file
!> =============================================================================
SUBROUTINE close (self, i)
  class(pebble_trace_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  call self%trace_output (i, 0)
END SUBROUTINE close

!> =============================================================================
!> Flush and close all pebble trace buffers
!> =============================================================================
SUBROUTINE close_all (self)
  class(pebble_trace_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  do i=1,self%nbuffer
    call self%close (i)
  end do
END SUBROUTINE close_all

END MODULE pebble_trace_mod
