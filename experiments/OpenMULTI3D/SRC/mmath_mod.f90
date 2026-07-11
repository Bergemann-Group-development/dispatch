MODULE mmath_mod
  USE refsor_mod

  implicit none

  interface sort
    module procedure sort_i, sort_8, sort_4
  end interface sort

  interface linear1d
    module procedure dlinear1d,rlinear1d ! real(8) and real versions
  end interface linear1d

  interface locate
    module procedure dlocate,rlocate ! real(8) and real versions
  end interface locate

  interface hunt
    module procedure dhunt,rhunt ! real(8) and real versions
  end interface hunt

  interface har_mean_deriv
    module procedure har_mean_deriv_noneq,har_mean_deriv_eq
  end interface har_mean_deriv

  interface hermite_interpolation
    module procedure hermite_interpolation_single, &
      hermite_interpolation_1d, &
      hermite_interpolation_1d4
  end interface hermite_interpolation

  interface ei1
    module procedure dei1,rei1 ! real(8) and real versions
  end interface ei1


contains

  ! *******************************************************************

  ! --------------------------------------------------------------------------
  !
  ! --------------------------------------------------------------------------
  recursive subroutine sort_i(a, first, last)
    ! *************************************************************************
    !
    ! quicksort.f -*-f90-*-
    ! Author: t-nissie
    ! License: GPLv3
    ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !
    !  AJG: - changed name from quicksort to sort
    !       - added double and single precisions and the interface
    !
    ! *************************************************************************

    implicit none

    integer :: a(*), x, t
    integer :: first, last
    integer :: i, j

    ! *************************************************************************

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
        i=i+1
      end do
      do while (x < a(j))
        j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call sort_i(a, first, i-1)
    if (j+1 < last)  call sort_i(a, j+1, last)
  END SUBROUTINE sort_i
  ! --------------------------------------------------------------------------
  !
  ! --------------------------------------------------------------------------
  recursive subroutine sort_8(a, first, last)
    ! *************************************************************************
    !
    ! quicksort.f -*-f90-*-
    ! Author: t-nissie
    ! License: GPLv3
    ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !
    !  AJG: - changed name from quicksort to sort
    !       - added double and single precisions and the interface
    !
    ! *************************************************************************

    implicit none

    real(8) :: a(*), x, t
    integer :: first, last
    integer :: i, j

    ! *************************************************************************

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
        i=i+1
      end do
      do while (x < a(j))
        j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call sort_8(a, first, i-1)
    if (j+1 < last)  call sort_8(a, j+1, last)
  END SUBROUTINE sort_8
  ! --------------------------------------------------------------------------
  !
  ! --------------------------------------------------------------------------
  recursive subroutine sort_4(a, first, last)
    ! *************************************************************************
    !
    ! quicksort.f -*-f90-*-
    ! Author: t-nissie
    ! License: GPLv3
    ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !
    !  AJG: - changed name from quicksort to sort
    !       - added double and single precisions and the interface
    !
    ! *************************************************************************

    implicit none

    real(4) :: a(*), x, t
    integer :: first, last
    integer :: i, j

    ! *************************************************************************

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
        i=i+1
      end do
      do while (x < a(j))
        j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i-1) call sort_4(a, first, i-1)
    if (j+1 < last)  call sort_4(a, j+1, last)
  END SUBROUTINE sort_4
  !----------------------------------------------------------------------------

  ! *******************************************************************

  pure function dei1(x)
    !---------------------------------------------------------------------
    ! --- Gives an approximation to the first exponential
    ! --- integral. Uses polynomial approximations from Abramowitz and
    ! --- Stegun, 1970, handbook of Mathematical Functions, p. 231
    ! --- for x>1 (Absolute accuracy is 5e-5)
    ! --- for 0<x<1 uses polynomial approx comPuted with mpfitfun and expint
    ! --- in IDL. Absolute accuracy is 1E-4
    !---------------------------------------------------------------------
    implicit none

    real(8) :: dei1
    real(8),intent(in) :: x
    real(8) :: c1,c2,x2

    real(8),parameter :: &
      s0 = -0.577117, &
      s1 =  0.997877, &
      s2 = -0.239964, &
      s3 =  0.0386789
    real(8),parameter :: &
      a1=2.334733,&
      a2=0.250621,&
      b1=3.330657,&
      b2=1.681534


    if (x.ge.1.0) then
      x2=x*x
      c1=x2+a1*x+a2
      c2=x2+b1*x+b2
      dei1=c1/(c2*x*exp(x))
    else if ((x.ge.0.0).and.(x.le.1.0)) then
      c1= s0 + x*(s1+x*(s2+x*s3))
      dei1 = c1-log(x)
    else
      dei1=0.0
    end if

  end function dei1

  ! ***********************************************************************

  pure function rei1(x)
    !---------------------------------------------------------------------
    ! --- Gives an approximation to the first exponential
    ! --- integral. Uses polynomial approximations from Abramowitz and
    ! --- Stegun, 1970, Handbook of Mathematical Functions, p. 231
    ! --- Absolute accuracy is 2e-7 for 0<x<1 and 2e-8 for 1<x<infinity.
    ! --- Mar  7 2005 :
    ! ---   Created by Jorrit Leenaarts (j.leenaarts@astro.uu.nl).
    !---------------------------------------------------------------------
    implicit none

    real(kind=4) :: rei1
    real(kind=4),intent(in) :: x
    real(kind=4) :: c1,c2,x2

    real(kind=4),parameter :: &
      s0 = -0.577117, &
      s1 =  0.997877, &
      s2 = -0.239964, &
      s3 =  0.0386789
    real(kind=4),parameter :: &
      a1=2.334733,&
      a2=0.250621,&
      b1=3.330657,&
      b2=1.681534

    if (x.ge.1.0) then
      x2=x*x
      c1=x2+a1*x+a2
      c2=x2+b1*x+b2
      rei1=c1/(c2*x*exp(x))
    else if ((x.ge.0.0).and.(x.le.1.0)) then
      c1= s0 + x*(s1+x*(s2+x*s3))
      rei1 = c1-log(x)
    else
      rei1=0.0
    end if

  end function rei1

  ! *******************************************************************

  PURE SUBROUTINE har_mean_deriv_noneq(wprime,x,y)
    ! *******************************************************************
    !
    ! compute derivative at y(1) using a weighted harmonic mean, see
    ! Ibgui et al. 2013,A&A, "IRIS,a generic 3D radiative transfer code"
    ! x can be nonequidistant
    ! *******************************************************************
    USE const_mod
    implicit none

    real(8),dimension(3),intent(in)  :: x,y
    real(8),             intent(out) :: wprime
    real(8)                          :: fim1,fi,alpha,d32

    fim1 =  y(2)-y(1)
    fi   =  y(3)-y(2)

    if (fim1*fi>zero) then

      d32=x(3)-x(2)
      fim1=fim1/(x(2)-x(1))
      fi=fi/d32
      alpha  = one_third * ( one + d32/(x(3)-x(1)) )
      wprime = fim1*fi / ( (one-alpha)*fim1 + alpha*fi )

    else
      wprime = zero
    endif

  END SUBROUTINE har_mean_deriv_noneq

  ! *******************************************************************

  PURE SUBROUTINE har_mean_deriv_eq(wprime,dx,y)
    ! *******************************************************************
    !
    ! compute derivative at y(1) using a weighted harmonic mean, see
    ! Ibgui et al. 2012,A&A, "IRIS,a generic 3D radiative transfer code"
    ! x axis is assumed equidistant, so only dx is given
    ! *******************************************************************
    USE const_mod
    implicit none

    real(8),dimension(3),intent(in)  :: y
    real(8),             intent(in)  :: dx
    real(8),             intent(out) :: wprime
    real(8)                          :: fim1,fi,p

    fim1 = y(2)-y(1)
    fi   = y(3)-y(2)
    p    = fim1*fi

    if (p>zero) then
      wprime = 2.0*p / ((fim1+fi)*dx)
    else
      wprime = zero
    endif

  END SUBROUTINE har_mean_deriv_eq

  ! *******************************************************************

  PURE SUBROUTINE hermite_derivative(dydx,x,y)
    ! *******************************************************************
    !
    ! compute derivative at y(1) using a weighted harmonic mean, see
    ! Ibgui et al. 2012,A&A, "IRIS,a generic 3D radiative transfer code"
    !
    ! *******************************************************************
    implicit none

    real(8),dimension(:),intent(in)  :: x,y
    real(8),dimension(:),intent(out) :: dydx

    integer :: i,n

    n=size(x)

    ! end points have one sided derivative
    dydx(1)= (y(2)-y(1)   ) /( x(2)-x(1  ) )
    dydx(n)= (y(n)-y(n-1) ) /( x(n)-x(n-1) )

    ! middle points
    if (n>2) then
      do i=2,n-1
        call  har_mean_deriv(dydx(i), x(i-1:i+1), y(i-1:i+1))
      end do
    end if


  END SUBROUTINE hermite_derivative

  ! *******************************************************************

  PURE SUBROUTINE hermite_interpolation_single(yy,xx,x,y,dydx)
    implicit none

    real(kind=8), dimension(:), intent(in)  :: x,y,dydx
    real(8),intent(in) :: xx
    real(8),intent(out):: yy
    integer :: n,i
    real(8) :: q,q2,q3,dx,p

    n=size(x)

    if (xx<=x(1)) then
      yy = y(1)
    else if (xx>=x(n)) then
      yy = y(n)
    else

      call hunt(n, x, xx, i)

      dx = x(i+1)-x(i)
      q = (xx-x(i)) / dx
      q2=q**2
      q3=q2*q
      p = 3.0d0*q2 - 2.0d0*q3

      yy = (1.0d0-p)*y(i) + p*y(i+1) + dx * ((q3-2.0d0*q2 +q)*dydx(i) + (q3-q2)*dydx(i+1))

    end if

  END SUBROUTINE hermite_interpolation_single

  ! *******************************************************************

  PURE SUBROUTINE hermite_interpolation_1d(yy,xx,x,y,dydx)
    implicit none

    real(kind=8), dimension(:), intent(in)  :: x,y,dydx
    real(8), dimension(:),intent(in) :: xx
    real(8), dimension(:),intent(out):: yy
    integer :: n,i,k
    real(8) :: q,q2,q3,dx,p

    n=size(x)

    do k=1,size(xx)

      if (xx(k)<=x(1)) then
        yy(k) = y(1)
      else if (xx(k)>=x(n)) then
        yy(k) = y(n)
      else

        call hunt(n, x, xx(k), i)

        dx = x(i+1)-x(i)
        q = (xx(k)-x(i)) / dx
        q2=q**2
        q3=q2*q
        p = 3.0d0*q2 - 2.0d0*q3

        yy(k) = (1.0d0-p)*y(i) + p*y(i+1) + dx * ((q3-2.0d0*q2 +q)*dydx(i) + (q3-q2)*dydx(i+1))

      end if

    end do

  END SUBROUTINE hermite_interpolation_1d

  ! *******************************************************************

  PURE SUBROUTINE hermite_interpolation_1d4(yy,xx,x,y,dydx)
    implicit none

    real(kind=8), dimension(:), intent(in)  :: x,y,dydx
    real(8), dimension(:),intent(in) :: xx
    real(4), dimension(:),intent(out):: yy
    integer :: n,i,k
    real(8) :: q,q2,q3,dx,p

    n=size(x)

    do k=1,size(xx)

      if (xx(k)<=x(1)) then
        yy(k) = y(1)
      else if (xx(k)>=x(n)) then
        yy(k) = y(n)
      else

        call hunt(n, x, xx(k), i)

        dx = x(i+1)-x(i)
        q = (xx(k)-x(i)) / dx
        q2=q**2
        q3=q2*q
        p = 3.0d0*q2 - 2.0d0*q3

        yy(k) = (1.0d0-p)*y(i) + p*y(i+1) + dx * ((q3-2.0d0*q2 +q)*dydx(i) + (q3-q2)*dydx(i+1))

      end if

    end do

  END SUBROUTINE hermite_interpolation_1d4


  ! ***********************************************************************

  SUBROUTINE dlinear1d(ntable, xtable, ytable, x, y,hnt)
    !---------------------------------------------------------------------
    ! Linear interpolation routine.
    ! Generates interpolate y at value x using the (xtable,ytable) pair
    ! When x is out of the table range, y=ytable(1) or y=ytable(ntable) is
    ! returned.
    ! This version requires double precision input.
    !---------------------------------------------------------------------
    implicit none

    real(8), dimension(:), intent(in)  :: xtable, ytable ! x,y arrays of table
    real(8),               intent(in)  :: x              ! array of points where interpolation is desired
    real(8),               intent(out) :: y              ! interpolated values
    integer,               intent(in)  :: Ntable         ! number of points in table
    logical,intent(in) :: hnt

    logical       :: ascend
    integer,save  :: j
    real(8) :: xmin, xmax, fx, ymin,ymax

    if (xtable(2).ge.xtable(1)) then
      ascend=.true.
      xmin=xtable(1)
      xmax=xtable(Ntable)
      ymin=ytable(1)
      ymax=ytable(Ntable)
    else
      ascend=.false.
      xmin=xtable(Ntable)
      xmax=xtable(1)
      ymin=ytable(Ntable)
      ymax=ytable(1)
    end if

    if (x.le.xmin) then
      y = ymin
    else if (x.ge.xmax) then
      y = ymax
    else

      if (hnt) then
        call dhunt(Ntable, xtable, x, j)
      else
        call dlocate(Ntable, xtable, x, j)
      end if

      fx = (xtable(j+1) - x) / (xtable(j+1) - xtable(j))
      y = fx*ytable(j) + (1.0d0 - fx)*ytable(j+1)

    end if

  END SUBROUTINE dlinear1d

  ! ***********************************************************************

  SUBROUTINE rlinear1d(ntable, xtable, ytable, x, y,hnt)
    !---------------------------------------------------------------------
    ! Linear interpolation routine.
    ! Generates interpolate y at value x using the (xtable,ytable) pair
    ! When x is out of the table range, y=ytable(1) or y=ytable(ntable) is
    ! returned.
    ! This version requires double precision input.
    !---------------------------------------------------------------------
    implicit none

    real(kind=4), dimension(:), intent(in)  :: xtable, ytable ! x,y arrays of table
    real(kind=4),               intent(in)  :: x              ! point where interpolation is desired
    real(kind=4),               intent(out) :: y              ! interpolated value
    integer,                    intent(in)  :: Ntable         ! number of points in table
    logical,intent(in) :: hnt

    logical  :: ascend
    integer  :: j
    real     :: xmin, xmax, fx, ymin,ymax

    if (xtable(2).ge.xtable(1)) then
      ascend=.true.
      xmin=xtable(1)
      xmax=xtable(Ntable)
      ymin=ytable(1)
      ymax=ytable(Ntable)
    else
      ascend=.false.
      xmin=xtable(Ntable)
      xmax=xtable(1)
      ymin=ytable(Ntable)
      ymax=ytable(1)
    end if

    if (x.le.xmin) then
      y = ymin
    else if (x.ge.xmax) then
      y = ymax
    else

      if (hnt) then
        call rhunt(Ntable, xtable, x, j)
      else
        call rlocate(Ntable, xtable, x, j)
      end if
      fx = (xtable(j+1) - x) / (xtable(j+1) - xtable(j))
      y = fx*ytable(j) + (1 - fx)*ytable(j+1)

    end if

  END SUBROUTINE rlinear1d

  ! ***********************************************************************

  PURE SUBROUTINE dlocate(n,xx,x,j)
    !---------------------------------------------------------------------
    ! Given an array xx(1:n) and given a value x, returns a value such
    ! that x is between xx(j) and xx(j+1), xx must be monotonic, j=0
    ! or j=n is returned if x is out of range. See Numerical Recipes.
    ! This version requires double precision input.
    !---------------------------------------------------------------------
    implicit none

    integer             ,intent(in)    :: n
    integer             ,intent(inout) :: j
    real(8)             ,intent(in)    :: x
    real(8),dimension(:),intent(in)    :: xx

    integer :: jl,jm,ju

    if (x.lt.xx(1)) then
      j=0 ! x below table range
    else if (x.gt.xx(n)) then
      j=n ! x above table range
    else
      ! use bisection
      jl=0
      ju=n+1
      do
        if (ju-jl.gt.1) then
          jm=(ju+jl)/2
          if( (xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)) ) then
            jl=jm
          else
            ju=jm
          end if
        else
          j = jl
          exit
        end if
      end do
    endif

  END SUBROUTINE dlocate

  ! ***********************************************************************

  PURE SUBROUTINE rlocate(n,xx,x,j)
    !---------------------------------------------------------------------
    ! Given an array xx(1:n) and given a value x, returns a value such
    ! that x is between xx(j) and xx(j+1), xx must be monotonic, j=0
    ! or j=n is returned if x is out of range. See Numerical Recipes.
    ! This version requires double precision input.
    !---------------------------------------------------------------------
    implicit none

    integer                  ,intent(in)    :: n
    integer                  ,intent(inout) :: j
    real(kind=4)             ,intent(in)    :: x
    real(kind=4),dimension(:),intent(in)    :: xx

    integer :: jl,jm,ju

    if (x.lt.xx(1)) then
      j=0 ! x below table range
    else if (x.gt.xx(n)) then
      j=n ! x above table range
    else
      ! use bisection
      jl=0
      ju=n+1
      do
        if (ju-jl.gt.1) then
          jm=(ju+jl)/2
          if( (xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)) ) then
            jl=jm
          else
            ju=jm
          end if
        else
          j = jl
          exit
        end if
      end do
    endif

  END SUBROUTINE rlocate

  ! ***********************************************************************

  PURE SUBROUTINE dhunt(n,xx,x,jlo)
    !---------------------------------------------------------------------
    ! Given an array xx(1:N), and given a value x, returns a value jlo
    ! such that x is between xx(jlo) and xx(jlo+1). xx must be
    ! monotonicly increasing, jlo=0 or jlo=N is returned if x is out of range.
    ! the value of jlo should be used as initial guess, otherwise HUNT
    ! is not faster than LOCATE.
    ! May 28 2006 : Copied from Numerical Recipes.
    !---------------------------------------------------------------------

    implicit none

    integer              , intent(in)    :: n
    real(8),               intent(in)    :: x
    real(8), dimension(n), intent(in)    :: xx
    integer              , intent(inout) :: jlo

    integer :: inc,jhi,jm
    logical :: ascnd

    ascnd=(xx(n).gt.xx(1))  ! True if ascending order of table, false otherwise.
    if (jlo <= 0 .or. jlo > n) then
      ! Input guess not useful. Go immediately to bisection.
      jlo=0
      jhi=n+1
    else
      inc=1 ! Set the hunting increment.
      if (x.ge.xx(jlo).eqv.ascnd) then ! Hunt up:
        do
          jhi=jlo+inc
          if(jhi.gt.n) then ! Done hunting, since off end of table.
            jhi=n+1
            exit
          else
            if (x.lt.xx(jhi).eqv.ascnd) exit
            jlo=jhi ! Not done hunting,
            inc=inc+inc ! so double the increment
          end if
        end do !and try again.

      else ! Hunt down:

        jhi=jlo
        do
          jlo=jhi-inc
          if(jlo.le.0) then ! Done hunting, since off end of table.
            jlo=0
            exit
          else
            if (x.gt.xx(jlo).eqv.ascnd) exit
            jhi=jlo !Not done hunting,
            inc=inc+inc !so double the increment
          end if
        end do !and try again.
      end if
    end if !Done hunting, value bracketed.

    do !Hunt is done, so begin the final bisection phase:
      if (jhi-jlo == 1) then
        exit
      else
        jm=(jhi+jlo)/2
        if (x>xx(jm) .eqv. ascnd) then
          jlo=jm
        else
          jhi=jm
        end if
      end if
    end do

  END SUBROUTINE dhunt

  ! ***********************************************************************

  PURE SUBROUTINE rhunt(n,xx,x,jlo)
    !---------------------------------------------------------------------
    ! Given an array xx(1:N), and given a value x, returns a value jlo
    ! such that x is between xx(jlo) and xx(jlo+1). xx must be
    ! monotonicly increasing, jlo=0 or jlo=N is returned to indicate
    ! that x is out of range.
    ! the value of jlo should be used as initial guess, otherwise HUNT
    ! is not faster than LOCATE.
    ! May 28 2006 : Copied from Numerical Recipes.
    !---------------------------------------------------------------------

    implicit none

    integer                   , intent(in)    :: n
    real(kind=4)              , intent(in)    :: x
    real(kind=4), dimension(n), intent(in)    :: xx
    integer                   , intent(inout) :: jlo

    integer :: inc,jhi,jm
    logical :: ascnd

    ascnd=(xx(n).ge.xx(1))  ! True if ascending order of table, false otherwise.
    if (jlo <= 0 .or. jlo > n) then
      ! Input guess not useful. Go immediately to bisection.
      jlo=0
      jhi=n+1
    else
      inc=1 ! Set the hunting increment.
      if (x.ge.xx(jlo).eqv.ascnd) then ! Hunt up:
        do
          jhi=jlo+inc
          if(jhi.gt.n) then ! Done hunting, since off end of table.
            jhi=n+1
            exit
          else
            if (x.lt.xx(jhi).eqv.ascnd) exit
            jlo=jhi ! Not done hunting,
            inc=inc+inc ! so double the increment
          end if
        end do !and try again.

      else ! Hunt down:

        jhi=jlo
        do
          jlo=jhi-inc
          if(jlo.lt.1) then ! Done hunting, since off end of table.
            jlo=0
            exit
          else
            if (x.gt.xx(jlo).eqv.ascnd) exit
            jhi=jlo !Not done hunting,
            inc=inc+inc !so double the increment
          end if
        end do !and try again.
      end if
    end if !Done hunting, value bracketed.

    do !Hunt is done, so begin the final bisection phase:
      if (jhi-jlo <= 1) then
        if (x== xx(n)) jlo=n-1
        if (x== xx(1)) jlo=1
        exit
      else
        jm=(jhi+jlo)/2
        if (x>= xx(jm) .eqv. ascnd) then
          jlo=jm
        else
          jhi=jm
        end if
      end if
    end do

  END SUBROUTINE rhunt

  !---------------------------------------------------------------------
  SUBROUTINE cc_kernel(s,k)
    !---------------------------------------------------------------------
    !
    !  Computes the 1D cubic convolution kernel. For an equidistant
    !  grid x_i (i={1,2,3,4}) with spacing h, a coordinate x_2<x<x_3
    !  where the function value is desired, and s= (x-x_2)/h, the
    !  kernel k is computed, so that y(x) approx sigma_1^4 k[i]*y_i
    !  for tabulated values y_i at position x_i
    !
    !---------------------------------------------------------------------
    implicit none

    real(8),intent(in)               :: s
    real(8),intent(out),dimension(4) :: k
    real(8)                          :: s2,s3

    s2=s*s
    s3=s2*s

    k(1) = 0.5*(-s3+2*s2-s)
    k(2) = 0.5*( 3*s3-5*s2+2)
    k(3) = 0.5*(-3*s3+4*s2+s)
    k(4) = 0.5*(s3-s2)

  END SUBROUTINE cc_kernel
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE lineq(a,b,n)
    !***********************************************************************
    !
    !  finds solution of system of linear equations
    !  with gaussian elimination with pivoting
    !
    !: lineq  90-06-05  new routine: (martin j stift)
    !:
    !***********************************************************************
    !
    implicit none
    !
    ! declaration of the supplied variables
    !
    integer :: n
    real(8), dimension(n,n) :: a
    real(8), dimension(n)   :: b

    real(8), dimension(n)  :: c
    integer, dimension(n) :: icol
    integer :: i, ibeg, jbeg, ima, jma, j, it, n1, jmin, i1, imin
    real(8)  :: ama, temp, quot

    ! initialize column count and starting points

    do i = 1,n
      icol(i) = i
    end do

    ibeg = 1
    jbeg = 1

    ! determine pivot

20  ama = abs(a(ibeg,jbeg))
    ima = ibeg
    jma = jbeg
    do i = ibeg,n
      do j = jbeg,n
        if(abs(a(i,j)).le.ama) cycle
        ama = abs(a(i,j))
        ima = i
        jma = j
      end do
    end do

    ! order matrix depending on pivot

    do i = 1,n
      temp = a(i,jma)
      a(i,jma) = a(i,jbeg)
      a(i,jbeg) = temp
    end do

    do j = jbeg,n
      temp = a(ima,j)
      a(ima,j) = a(ibeg,j)
      a(ibeg,j) = temp
    end do

    temp = b(ima)
    b(ima) = b(ibeg)
    b(ibeg) = temp
    it = icol(jbeg)
    icol(jbeg) = icol(jma)
    icol(jma)= it
    ibeg = ibeg + 1
    jbeg = jbeg + 1

    ! eliminate

    imin = ibeg - 1
    jmin = jbeg - 1
    do i = ibeg,n
      quot = a(i,jmin) / a(imin,jmin)
      do j = jbeg,n
        a(i,j) = a(i,j) - quot * a(imin,j)
      end do
      b(i) = b(i) - quot * b(imin)
    end do
    if (ibeg.lt.n) goto 20

    ! determine coefficients

    b(n) = b(n) / a(n,n)
    n1 = n - 1
    do i = n1,1,-1
      i1 = i + 1
      do j = n,i1,-1
        b(i) = b(i) - a(i,j) * b(j)
      end do
      b(i) = b(i) / a(i,i)
    end do

    ! reorder coefficients

    do i = 1,n
      c(icol(i))=b(i)
    end do

    do i = 1,n
      b(i) = c(i)
    end do

  END SUBROUTINE lineq
  !---------------------------------------------------------------------

END MODULE mmath_mod
! *******************************************************************j
