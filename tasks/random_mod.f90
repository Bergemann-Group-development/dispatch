!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Simple random number module
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE random_mod
  USE io_mod
  implicit none
  private
  type, public:: random_t
    integer:: seed=11
  contains
    procedure:: init
    procedure:: ran1
    procedure:: ran3
    procedure:: uniform_1
    procedure:: uniform_n
  end type
  interface uniform
    module procedure uniform_1, uniform_n
  end interface
  type(random_t), public:: random
CONTAINS

!> =============================================================================
!> Initialize both random number generators
!> =============================================================================
SUBROUTINE init (self, seed)
  class(random_t):: self
  integer, optional:: seed
  real:: a, b(3)
  !.............................................................................
  if (present(seed)) then
    self%seed = abs(seed)
  else
    self%seed = 11
  end if
  call random_seed
  !call random_number (a)
  !call random_number (b)
  !if (io%master) then
  !  print*,'random_number:',a
  !  print*,'random_number:',b
  !end if
END SUBROUTINE init

!> =============================================================================
!> Ancient Numerical Recipes code.  Will do for now.
!> =============================================================================
FUNCTION ran1(self)
  class(random_t):: self
  integer:: idum,IA,IM,IQ,IR,NTAB,NDIV
  real:: ran1,AM,EPS,RNMX
  parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
    NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  integer:: k,iy
  !.............................................................................
  idum = self%seed
  k =idum/IQ
  idum = IA*(idum-k*IQ)-IR*k
  if (idum<0) idum = idum+IM
  iy = idum
  ran1 = min(AM*iy,RNMX)
  self%seed = idum
END FUNCTION

!> =============================================================================
!> 3-component random vector, components uniform in a spherical distribution,
!> with max radius = 1.0
!> =============================================================================
FUNCTION ran3(self) result (u)
  class(random_t):: self
  real(kind=4):: u(3), u2
  integer:: i
  !-----------------------------------------------------------------------------
  do i=1,3
    u(i) = 2.*self%ran1()-1.
  end do
  u2 = sum(u**2)
  do while (u2 > 1.0 .or. u2 < 1e-8)
    do i=1,3
      u(i) = 2.*self%ran1()-1.
    end do
    u2 = sum(u**2)
  end do
END FUNCTION ran3

!> =============================================================================
!> =============================================================================
FUNCTION uniform_1 (self, a, b) RESULT (out)
  class(random_t):: self
  real, intent(in), optional:: a, b
  real:: out
  !-----------------------------------------------------------------------------
  if (present(a)) then
    if (present(b)) then
      out = a*self%ran1()
    else
      out = a + (b-a)*self%ran1()
    end if
  else
    out = self%ran1()
  end if
END FUNCTION

!> =============================================================================
!> =============================================================================
FUNCTION uniform_n (self, n, a, b) RESULT (out)
  class(random_t):: self
  integer, intent(in):: n
  real, intent(in), optional:: a, b
  real, dimension(n):: out
  integer:: i
  !-----------------------------------------------------------------------------
  if (present(a)) then
    if (present(b)) then
      do i=1,n
        out(i) = a*self%ran1()
      end do
    else
      do i=1,n
        out(i) = a + (b-a)*self%ran1()
      end do
    end if
  else
    do i=1,n
      out(i) = self%ran1()
    end do
  end if
END FUNCTION

END MODULE random_mod
