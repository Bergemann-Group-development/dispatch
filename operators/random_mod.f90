!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Simple random number module
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE random_mod
  USE io_mod
  implicit none
  private
  type, public:: random_t
    integer:: seed=-11
  contains
    procedure:: init
    procedure:: ran1
    procedure:: ran3
  end type
CONTAINS

!> =============================================================================
!> Initialize both random number generators
!> =============================================================================
SUBROUTINE init (self)
  class(random_t):: self
  real:: a, b(3)
  !.............................................................................
  a = self%ran1()
  print*,a
  call random_seed
  call random_number (a)
  call random_number (b)
  if (io%master) then
    print*,'random_number:',a
    print*,'random_number:',b
  end if
END SUBROUTINE init

!> =============================================================================
!> Ancient Numerical Recipes code.  Will do for now.
!> =============================================================================
FUNCTION ran1(self)
  class(random_t):: self
  INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
  REAL ran1,AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER k,iy
  !.............................................................................
  idum = self%seed
  if (idum<=0) then
    idum = max(-idum,1)
  endif
  k = idum/IQ
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

END MODULE random_mod
