!===============================================================================
!> Module to do detailed timing inside a code procedure, adding up contributions
!> and printing the results a given number of times.
!===============================================================================
MODULE used_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  private
  type use_t
    real(8):: time=0d0
    character(len=16):: label=''
  end type
  type, public:: used_t
    integer:: i=0, n=0, set=0, n_max, n_print=10
    real(8):: start, now
    type(use_t), allocatable:: use(:)
  contains
    procedure:: init
    procedure:: update
    procedure:: print
  end type
CONTAINS

!===============================================================================
!> Initialize, setting the max number of entries and number of times to print
!===============================================================================
SUBROUTINE init (self, n_max, n_print)
  class(used_t):: self
  integer, optional:: n_max, n_print
  integer:: n_m, n_p
  !-----------------------------------------------------------------------------
  if (omp%master .and. self%n_print > 0) then
    if (present(n_max)) then
      n_m = n_max
    else
      n_m = 50
    end if
    if (present(n_print)) then
      n_p = n_print
    else
      n_p = 5
    end if
    if (.not.allocated(self%use)) then
      allocate (self%use(n_m))
      self%n_print = n_p
    end if
    self%use(:)%time = 0d0
    self%i = 0
    self%n_max = n_m
    write (stdout,'("used_t%init: n_max, n_print =",2i5)') n_m, self%n_print
    self%start = wallclock()
  end if
END SUBROUTINE init

!===============================================================================
!> Update: Increment slot counter, and add the time since last call there
!===============================================================================
SUBROUTINE update (self, label, set, reset)
  class(used_t):: self
  character(len=*), optional:: label
  logical, optional:: set, reset
  integer:: i
  real(8):: now
  !-----------------------------------------------------------------------------
  if (omp%master .and. self%n_print > 0) then
    if (present(reset)) then
      self%i = self%set
    end if
    now = wallclock()
    i = self%i
    if (i < self%n_max) then
      i = i+1
      if (present(label)) &
        self%use(i)%label = label
      self%use(i)%time = self%use(i)%time + now - self%start
      self%n = max(self%n,i)
      self%i = i
    end if
    self%start = now
    if (present(set)) then
      self%set = self%i
    end if
  end if
END SUBROUTINE update

!===============================================================================
!> Print a table of millisecond use
!===============================================================================
SUBROUTINE print (self)
  class(used_t):: self
  integer:: i, n_max
  real(8):: total
  !-----------------------------------------------------------------------------
  if (omp%master .and. self%n_print > 0) then
    self%n_print = self%n_print-1
    total = 0d0
    do i=1,self%n
      write (stdout,1) i, self%use(i)%time, trim(self%use(i)%label)
      1 format(i3,3p,f12.3,3x,a)
      total = total + self%use(i)%time
    end do
    write (stdout,1) self%n, total, 'TOTAL'
  end if
END SUBROUTINE print

END MODULE used_mod
