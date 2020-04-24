!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Module used for tracing particles, producing one output file for every (save
!> frequency) particles.  Each particle being traced has a buffer here, which is
!> being filled by trace_pebble() calls.  When full, trace_output() is called to
!> write out the buffer and reset it.
!>
!> Calls:
!> pebble_trace%trace_pebble (self, pebble, ip, level)   ! add one trace position
!> pebble_trace%trace_output (self, i)                   ! empty buffer nr. i
!> pebble_trace%check (self, pebble)                     ! close if pebble traced
!> pebble_trace%close (self, i)                          ! close buffer i
!> pebble_trace%close_all (self)                         ! close all buffers
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebble_trace_mod
  USE io_mod
  USE omp_mod
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
    integer:: npebble=20
    integer:: nbuffer=25000                     ! max number of traces
    integer:: every=1000                        ! save frequency
    logical:: do_trace_output=.true.
  contains
    procedure flush_buffer
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

!===============================================================================
SUBROUTINE init (self, trace_every)
  class(pebble_trace_t):: self
  integer:: trace_every, i
  !-----------------------------------------------------------------------------
  if (.not.pebble_trace%do_trace_output) return

  pebble_trace%every = trace_every
  if (trace_every > 0) then
    allocate(pebble_trace%buffer(pebble_trace%nbuffer))
    do i=1,pebble_trace%nbuffer
      allocate (pebble_trace%buffer(i)%v(21,pebble_trace%npebble))
      allocate (pebble_trace%buffer(i)%ip(pebble_trace%npebble))
      pebble_trace%buffer(i)%n = 0
    end do
  end if

END SUBROUTINE init

!===============================================================================
!> Add a pebble to the trace buffer, and check for output.  A trace pebble has
!> a status word where the first 16 bits is a consecutive trace number, starting
!> with 1, and running up to nbuffer.  The 17th and 18th bit are set to mark the
!> pebble as being traced, and as being launched from the boundary, respectively.
!> The 27th bit marks a pebble that has been accreted.
!===============================================================================
SUBROUTINE trace_pebble (self, pebble, ip, level)
  class(pebble_trace_t):: self
  class(pebble_t), pointer:: pebble
  integer:: ip, level
  integer:: i, n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.pebble_trace%do_trace_output) return
  call trace%begin ('pebble_trace_t%trace_pebble', itimer=itimer)

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
    self%buffer(i)%v(15:17,n) = pebble%ad
    self%buffer(i)%v(18:20,n) = pebble%ag
    self%buffer(i)%v(21,n)  = pebble%dt
    self%buffer(i)%n = n
    if (io%verbose > 1) print *,'Pebble ID:', pebble%id, i, n, pebble%t, pebble%dt
    if (n >= self%npebble) then
      call self%trace_output (i)
    end if
  end if

  call trace%end (itimer)
END SUBROUTINE trace_pebble

!===============================================================================
!> Dump trace outputs for all pebbles
!===============================================================================
SUBROUTINE flush_buffer(self, parent_id)
  class(pebble_trace_t):: self
  integer :: i, parent_id
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.pebble_trace%do_trace_output) return
  call trace%begin ('pebble_trace_t%flush_buffer', itimer=itimer)

  if (io%verbose>0) &
    print*,'pebble_trace: parent ID:', parent_id
  do i=1, id
    call self%trace_output (i)
  end do

  call trace%end (itimer)
END SUBROUTINE flush_buffer

!===============================================================================
!> Write out a pebble trace buffer.  Each record consists of:
!>
!> id=trace_id, rec=record_number, n=record_size   (3 ints)
!> ip(1:n)      (n ints)
!> v(:,1:n)     (m*n floats, where m=nn, n=record_size)
!===============================================================================
SUBROUTINE trace_output (self, i, close_key)
  USE os_mod

  class(pebble_trace_t):: self
  integer:: i, close_key
  integer, parameter:: key=-1
  optional:: close_key
  character(len=64):: filename, dirname
  integer:: j, n, unit
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.pebble_trace%do_trace_output) return
  if (i > self%nbuffer) return
  if (.not.associated(self%buffer)) return
  if (self%buffer(i)%n == 0) return
  call trace%begin ('pebbles_mod::trace_output', itimer=itimer)

  !$omp critical (tr_out_cr)
  write (dirname, '(a,"traces/")') trim(io%outputname)
  write (filename,'(a,i5.5,".trc")') trim(dirname), i
  if (io%verbose>0) &
    print*,'pebble_trace: file, record, n: ', trim(filename), self%buffer(i)%record, self%buffer(i)%n
  unit = io_unit%tmp + omp%thread
  if (self%buffer(i)%record == 0) then
    if (io%do_output) call os%mkdir(dirname)
    open (unit, file=trim(filename), form='unformatted', &
      status='unknown')
    close (unit, status='delete')
    open (unit, file=trim(filename), form='unformatted', &
      status='new')
    write (unit) key, self%buffer(i)%status
    write (unit) self%buffer(i)%s
  else
    open (unit, file=trim(filename), form='unformatted', &
      status='unknown', position = 'append')
    backspace (unit)
  end if
  n =  self%buffer(i)%n
  self%buffer(i)%record = self%buffer(i)%record+1
  write (unit) self%buffer(i)%id, self%buffer(i)%record, self%buffer(i)%n
  write (unit) self%buffer(i)%ip(1:n)
  write (unit) self%buffer(i)%v(:,1:n)
  self%buffer(i)%n = 0
  write (unit) 0,0,0,0
  close (unit)
  !$omp end critical (tr_out_cr)

  call trace%end (itimer)
END SUBROUTINE trace_output

!===============================================================================
!> Check if this pebble is traced, and if so flush and close its file
!===============================================================================
SUBROUTINE check (self, pebble)
  class(pebble_trace_t):: self
  class(pebble_t), pointer:: pebble
  integer:: i
  !-----------------------------------------------------------------------------
  if (.not.pebble_trace%do_trace_output) return
  if (iand(pebble%status,bits%trace)==0) return
  i = iand(pebble%status,not(bits%trace+bits%has_BC))
  call self%close (i)
END SUBROUTINE check

!===============================================================================
!> Flush and close trace file
!===============================================================================
SUBROUTINE close (self, i)
  class(pebble_trace_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  call self%trace_output (i, 0)
END SUBROUTINE close

!===============================================================================
!> Flush and close all pebble trace buffers
!===============================================================================
SUBROUTINE close_all (self)
  class(pebble_trace_t):: self
  integer:: i
  !-----------------------------------------------------------------------------
  do i=1,self%nbuffer
    call self%close (i)
  end do
END SUBROUTINE close_all

END MODULE pebble_trace_mod
