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
  type(trace_buffer_t), pointer :: buffer
  integer:: trace_every, i
  !-----------------------------------------------------------------------------
  pebble_trace%every = trace_every
  if (trace_every > 0) then
    call pebble_trace%lock%init ('trace')
    allocate(pebble_trace%buffer(pebble_trace%nbuffer))
    do i=1,pebble_trace%nbuffer
      buffer => pebble_trace%buffer(i)
      allocate (buffer%v(14,pebble_trace%npebble))
      allocate (buffer%ip(pebble_trace%npebble))
      buffer%n = 0
    end do
  end if
  nullify (buffer)
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
  type(trace_buffer_t), pointer :: buffer
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
    buffer => self%buffer(i)
    buffer%id      = pebble%id
    buffer%status  = pebble%status
    buffer%s       = pebble%s
    n = buffer%n+1
    buffer%ip(n) = ip
    buffer%v(1:3,n) = pebble%p
    buffer%v(4:6,n) = pebble%v
    buffer%v(7:9,n) = pebble%vg
    buffer%v(10,n)  = pebble%c
    buffer%v(11,n)  = pebble%d
    buffer%v(12,n)  = pebble%t
    buffer%v(13,n)  = pebble%s
    buffer%v(14,n)  = pebble%w
    buffer%n = n
    if (io%verbose > 1) print *,'MK', pebble%id, i, n, pebble%t, pebble%dt
    if (n >= self%npebble) then
      call self%trace_output (i)
    end if
  else
    if (io%verbose > 1) print *,'MK', pebble%id, i
  end if
  nullify(buffer)
  call trace_end
END SUBROUTINE trace_pebble

!> =============================================================================
!> Write out a pebble trace buffer
!> =============================================================================
SUBROUTINE trace_output (self, i, close_key)
  class(pebble_trace_t):: self
  integer:: i, close_key
  type(trace_buffer_t), pointer :: buffer
  integer, parameter:: key=-1
  optional:: close_key
  character(len=64):: filename
  integer:: j, n
  !-----------------------------------------------------------------------------
  if (i>self%nbuffer) return
  if (.not.associated(self%buffer)) return
  buffer => self%buffer(i)
  if (buffer%n==0) return
  call trace_begin ('pebbles_mod::trace_output')
  call self%lock%set ('trace_output')
  write (filename,'(a,i4.4,".trc")') trim(io%outputname), i
  if (io%verbose>0) &
    print*,'pebble_trace: ', trim(filename), buffer%record, buffer%n
  if (buffer%record == 0) then
    open (io_unit%trace, file=trim(filename), form='unformatted', &
      status='unknown')
    close (io_unit%trace, status='delete')
    open (io_unit%trace, file=trim(filename), form='unformatted', &
      status='new')
    write (io_unit%trace) key, buffer%status
    write (io_unit%trace) buffer%s
  else
    open (io_unit%trace, file=trim(filename), form='unformatted', &
      status='unknown', position = 'append')
    backspace (io_unit%trace)
  end if
  n =  buffer%n
  buffer%record = buffer%record+1
  write (io_unit%trace) buffer%id, buffer%record, buffer%n
  write (io_unit%trace) buffer%ip(1:n)
  write (io_unit%trace) buffer%v(:,1:n)
  buffer%n = 0
  write (io_unit%trace) 0,0,0,0
  close (io_unit%trace)
  call self%lock%unset ('trace_output')
  nullify(buffer)
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
