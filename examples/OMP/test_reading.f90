MODULE mod
  implicit none
  type exp_t
    real:: a, b
    integer:: thread, i
    contains
    procedure:: init
    procedure:: run
  end type
CONTAINS

SUBROUTINE init (s, i)
  class(exp_t):: s
  integer:: i, omp_get_thread_num, thread
  real:: b
  namelist /b_param/ b, thread
  character(len=32):: filename
  !................................................
  read (*,b_param)
  thread = omp_get_thread_num()
  s%thread = thread
  s%b      = b
  s%i      = i
  write (filename,'(a,i2.2,a)') 'data/', s%i, '.dat'
  open (s%i+10, file=filename, form='formatted', status='unknown')
END SUBROUTINE

SUBROUTINE run (s)
  class(exp_t):: s
  write (s%i+10,*) s%b, s%thread, s%i
  close (s%i+10)
END SUBROUTINE

END MODULE mod

PROGRAM test
  USE mod
  implicit none
  integer:: i, n
  real:: a, b
  namelist /a_param/ a, n
  integer:: omp_get_num_threads, omp_get_thread_num
  type(exp_t), dimension(:), allocatable:: experiment

  read(*,a_param)
  write(*,a_param)
  allocate (experiment(n))

  !$omp parallel do
  do i=1,n
    call experiment(i)%init(i)
    call experiment(i)%run
  end do

END
