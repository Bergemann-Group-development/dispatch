MODULE mmath
  USE refsor_mod
  USE unista_mod

  implicit none

  interface sort
     module procedure sort_i, sort_8, sort_4
  end interface sort

  interface linint
     module procedure dlinint,rlinint ! real(8) and real versions
  end interface linint

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
  !
  !----------------------------------------------------------------------------
  integer function stat_mode(a, others, olen, ok)
    ! *************************************************************************
    !
    !  Taken from https://rosettacode.org/wiki/Averages/Mode#Fortran 
    !
    !  sort from original code replaced with sort interface from 
    !  https://gist.github.com/t-nissie/479f0f16966925fa29ea, which was 
    !  adapted by A. J. Gallagher for use with Multi3D
    !
    ! *************************************************************************
    
    ! dummy arguments
    integer, dimension(        :),            intent(in ) :: a      ! mode
    logical,                        optional, intent(out) :: ok     ! is a mode
    integer, dimension(size(a, 1)), optional, intent(out) :: others ! more modes
    integer,                        optional, intent(out) :: olen   ! others len

    ! local variables
    integer, dimension(size(a, 1)) :: ta    ! copy of input, a
    integer, dimension(size(a, 1)) :: freq  ! frequencies
    integer, dimension(size(a, 1)) :: idx   ! index of ta
    integer                        :: rs
    integer                        :: i
    integer                        :: tm
    integer                        :: ml
    integer                        :: tf

    ! *************************************************************************

    if ( present(ok) ) ok = .false.

    select case ( size(a, 1) )
    case (0)  ! no mode... ok is false
       return
    case (1)
       if ( present(ok) ) ok = .true.
       stat_mode = a(1)
       return
    case default
       if ( present(ok) ) ok = .true.
       ta   = a         ! copy the array
       call sort(ta,1,size(a, 1))  ! sort it in place (cfr. sort algos on RC)
       freq = 1
       idx  = 0
       rs   = 1         ! rs will be the number of different values

       do i = 2, size(ta, 1)
          if ( ta(i-1) == ta(i) ) then
             freq(rs) = freq(rs) + 1
          else
             idx(rs) = i  - 1
             rs      = rs + 1
          end if
       end do
       idx(rs) = i - 1

       ml = maxloc(freq(1:rs), 1)  ! index of the max value of freq
       tf = freq(ml)               ! the max frequency
       tm = ta(idx(ml))            ! the value with that freq

       ! if we want all the possible modes, we provide others
       if ( present(others) ) then
          i = 1
          others(1) = tm
          do
             freq(ml) = 0
             ml = maxloc(freq(1:rs), 1)
             if ( tf == freq(ml) ) then ! the same freq
                i = i + 1               ! as the max one
                others(i) = ta(idx(ml))
             else
                exit
             end if
          end do

          if ( present(olen) ) then
             olen = i
          end if

       end if
       stat_mode = tm
    end select

  end function stat_mode
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE ipol_cons(x1,y1,x2,y2)
    !+
    !
    !            interpolates (x1,y1) onto new grid x2
    !            conserving the integral y*dx
    !
    ! Nov 21 2014 JL: Fortran implementation of Mats Carlsson's 
    ! ipol_cons.pro IDL routine. 
    !
    !-
    implicit none

    real(8),dimension(:),intent(in) :: x1,y1,x2
    real(8),dimension(:),intent(out) :: y2

    real(8),dimension(:),allocatable :: ynew,xmerge,yint,xnew,xint
    integer,dimension(:),allocatable :: imerge,iold,inew

    real(8) :: int
    integer :: n1,n2,i,ii,nint

    n1=size(x1)
    n2=size(x2)

    allocate( ynew(n1+n2+1) )
    allocate( inew(n1+n2) )
    allocate(xmerge(n1+n2))
    xmerge(1:n1)=x1
    xmerge(n1+1:n1+n2)=x2

    allocate(imerge(n1+n2))
    do i=1,n1+n2 
       imerge(i)=i
    end do

    allocate(iold(n1+n2))
    call qsortd(xmerge,iold,n1+n2)
    inew(iold)=imerge

    allocate(xnew(n1+n2))
    xnew=xmerge(iold)


    ynew(n1+n2+1)=0.0d0
    ynew( inew(1:n1) ) = y1

    do i=n1+n2,n1+1,-1
       ynew(inew(i)) = ynew(inew(i)+1)
    end do

    do i=2,n2

       ii=i+n1
       nint=inew(ii)-inew(ii-1)+1

       if (allocated(yint)) deallocate(yint)
       allocate(yint(nint))
       yint = ynew( inew(ii-1):inew(ii) )

       if (allocated(xint)) deallocate(xint)
       allocate(xint(nint))
       xint = xnew( inew(ii-1):inew(ii) )

       ! integral of y*dx over new cell
       int=sum( ( xint(2:nint) - xint(1:nint-1) ) * yint(2:nint))
       y2(i)=int/(x2(i)-x2(i-1))  ! new value
    end do

  END SUBROUTINE ipol_cons

  ! ***********************************************************************

  pure function trapez(x,y)
    ! trapezoidal integration of tabulated values y as function of x
    implicit none
    real(8),dimension(:),intent(in) :: x,y

    real(8) :: trapez
    integer :: i,n

    n=size(x)

    trapez = 0.0d0
    do i=1,n-1
       trapez = trapez  + (y(i)+y(i+1)) * (x(i+1)-x(i))
    end do
    trapez=0.5*trapez

  end function trapez

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
    USE mconst    
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
    USE mconst
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

  SUBROUTINE dlinint(ntable,xtable,ytable,x,y,dxtable)
    !---------------------------------------------------------------------
    ! Gets interpolates value from a table with an equidistant axis, with
    ! xtable(ntable) > xtable(1).
    ! When x is out of the table range, y=ytable(1) or y=ytable(ntable) is
    ! returned.
    ! Use linear1d if table axis is not equidistant.
    ! This version needs double precision input.
    !---------------------------------------------------------------------
    implicit none
    real(kind=8), dimension(:), intent(in)  :: xtable, ytable
    real(kind=8),               intent(in)  :: x,dxtable
    real(kind=8),               intent(out) :: y
    integer,                    intent(in)  :: Ntable
    real(kind=8)                            :: ifr,tt,m1t
    integer                                 :: il

    ! get fractional x index
    ifr=(x-xtable(1))/dxtable+1.0d0
    il =aint(ifr)
    tt =ifr-il
    m1t=1.0d0-tt
    ! check ranges
    if (il.lt.1) then
       il =1
       tt =0.0d0
       m1t=1.0d0
    end if
    if (il.gt.ntable) then
       il =ntable-1
       tt =1.0d0
       m1t=0.0d0
    end if
    y=tt*ytable(il) + m1t*ytable(il+1)

  END SUBROUTINE dlinint

  ! ***********************************************************************

  SUBROUTINE rlinint(ntable,xtable,ytable,x,y,dxtable)
    !---------------------------------------------------------------------
    ! Gets interpolates value from a table with an equidistant axis, with
    ! xtable(ntable) > xtable(1).
    ! When x is out of the table range, y=ytable(1) or y=ytable(ntable) is 
    ! returned.
    ! Use linear1d if table axis is not equidistant
    ! This version has whatever is defined as real as input 
    !---------------------------------------------------------------------
    implicit none
    real(kind=4),dimension(:), intent(in)  :: xtable, ytable
    real(kind=4),              intent(in)  :: x,dxtable
    real(kind=4),              intent(out) :: y
    integer,                   intent(in)  :: Ntable
    real(kind=4)                           :: ifr,tt,m1t
    integer                                :: il

    ! get fractional x index
    ifr=(x-xtable(1))/dxtable+1.0d0
    il =aint(ifr)
    tt =ifr-il
    m1t=1.0d0-tt
    ! check ranges
    if (il.lt.1) then
       il =1
       tt =0.0d0
       m1t=1.0d0
    end if
    if (il.gt.ntable) then
       il =ntable-1
       tt =1.0d0
       m1t=0.0d0
    end if
    y=tt*ytable(il) + m1t*ytable(il+1)

  END SUBROUTINE rlinint

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

  ! ***********************************************************************

  SUBROUTINE cc_interp(x,kern,xx)
    implicit none
    real(8),dimension(4,4),intent(in):: kern,xx
    real(8),intent(out) :: x
    x=sum( kern * xx )
    x=max(x,minval(xx(2:3,2:3)))
  END SUBROUTINE cc_interp

  ! ***********************************************************************

  SUBROUTINE qsortd(x,ind,n)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2002-12-18  Time: 11:55:47

    implicit none
    real(8), intent(IN)  :: x(:)
    integer, intent(OUT)   :: ind(:)
    integer, intent(IN)    :: n

    !***************************************************************************

    !                                                         ROBERT RENKA
    !                                                 OAK RIDGE NATL. LAB.

    !   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
    ! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
    ! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL 
    ! INTERCHANGES
    ! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
    ! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
    ! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
    ! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
    ! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
    ! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
    ! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
    ! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
    ! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
    ! UNSORTED PORTION.

    ! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

    !                      X - VECTOR OF LENGTH N TO BE SORTED.

    !                    IND - VECTOR OF LENGTH >= N.

    ! N AND X ARE NOT ALTERED BY THIS ROUTINE.

    ! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
    !                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
    !                          X IS DEFINED BY Y(I) = X(IND(I)).

    !*********************************************************************

    ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

    !*********************************************************************

    integer   :: iu(21), il(21)
    integer   :: m, i, j, k, l, ij, it, itt, indx
    real      :: r
    real(8) :: t

    ! LOCAL PARAMETERS -

    ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
    !            INDICES OF PORTIONS OF THE ARRAY X
    ! M =      INDEX FOR IU AND IL
    ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
    ! K,L =    INDICES IN THE RANGE I,...,J
    ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
    ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
    ! INDX =   TEMPORARY INDEX FOR X
    ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
    ! T =      CENTRAL ELEMENT OF X

    if (n <= 0) return

    ! INITIALIZE IND, M, I, J, AND R

    do  i = 1, n
       ind(i) = i
    end do
    m = 1
    i = 1
    j = n
    r = .375

    ! TOP OF LOOP

20  if (i >= j) GO TO 70
    if (r <= .5898437) then
       r = r + .0390625
    else
       r = r - .21875
    end if

    ! INITIALIZE K

30  k = i

    ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

    ij = i + r*(j-i)
    it = ind(ij)
    t = x(it)

    ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
    !   INTERCHANGE IT WITH T

    indx = ind(i)
    if (x(indx) > t) then
       ind(ij) = indx
       ind(i) = it
       it = indx
       t = x(it)
    end if

    ! INITIALIZE L

    l = j

    ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
    !   INTERCHANGE IT WITH T

    indx = ind(j)
    if (x(indx) >= t) GO TO 50
    ind(ij) = indx
    ind(j) = it
    it = indx
    t = x(it)

    ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
    !   INTERCHANGE IT WITH T

    indx = ind(i)
    if (x(indx) <= t) GO TO 50
    ind(ij) = indx
    ind(i) = it
    it = indx
    t = x(it)
    GO TO 50

    ! INTERCHANGE ELEMENTS K AND L

40  itt = ind(l)
    ind(l) = ind(k)
    ind(k) = itt

    ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
    !   NOT LARGER THAN T

50  l = l - 1
    indx = ind(l)
    if (x(indx) > t) GO TO 50

    ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60  k = k + 1
    indx = ind(k)
    if (x(indx) < t) GO TO 60

    ! IF K <= L, INTERCHANGE ELEMENTS K AND L

    if (k <= l) GO TO 40

    ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
    !   ARRAY YET TO BE SORTED

    if (l-i > j-k) then
       il(m) = i
       iu(m) = l
       i = k
       m = m + 1
       GO TO 80
    end if

    il(m) = k
    iu(m) = j
    j = l
    m = m + 1
    GO TO 80

    ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70  m = m - 1
    if (m == 0) return
    i = il(m)
    j = iu(m)

80  if (j-i >= 11) GO TO 30
    if (i == 1) GO TO 20
    i = i - 1

    ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90  i = i + 1
    if (i == j) GO TO 70
    indx = ind(i+1)
    t = x(indx)
    it = indx
    indx = ind(i)
    if (x(indx) <= t) GO TO 90
    k = i

100 ind(k+1) = ind(k)
    k = k - 1
    indx = ind(k)
    if (t < x(indx)) GO TO 100

    ind(k+1) = it
    GO TO 90
  END SUBROUTINE qsortd
  !---------------------------------------------------------------------


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
  SUBROUTINE bicubic_convol(fi,xl,yl,f,xker,yker)
    !---------------------------------------------------------------------
    !
    ! performs bicubic convolution interpolation
    !
    ! idea is as follows: for 2D tabulated data f, with constant grid spacing dx
    ! in x and dy in y:
    !      1. find xl,yl such that   x(xl)<x<x(xl+1) and likewise for y
    !      2. compute normalized offset sx= (x-x(xl))/dx and likewise for y
    !      3. use cc_kernel to get the interpolation kernel for x and likewise for y
    !      4. use bicubic_convol (this routine) to get fi = f(x,y)
    !
    ! note that this method gives spurious extrema in case of large
    ! gradients (e.g. step function). Code has been added to remove
    ! this. The first derivative of f is in this case no longer
    ! continuous.
    !
    ! f: tabulated data
    ! xl,yl: (x,y) coordinates of stencil
    ! fi: contains interpolated value on output
    ! xker,yker : 1D cubic convolution kernels (should be computed with cc_kernel)
    !
    !---------------------------------------------------------------------
    implicit none

    real(8),               intent(out) :: fi
    real(8),dimension(:,:),intent(in)  :: f
    real(8),dimension(4),  intent(in)  :: xker,yker    
    integer(4),            intent(in)  :: xl,yl

    real(8),dimension(4) :: yt
    integer(4)           :: i

    do i=1,4
       yt(i)=sum(xker*f(xl-1:xl+2,i))
    end do
    fi=sum(yker*yt)

    ! check for and remove spurious extrema
    fi=min(fi,maxval(f(xl:xl+1,yl:yl+1)))
    fi=max(fi,minval(f(xl:xl+1,yl:yl+1)))

  END SUBROUTINE bicubic_convol
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



  function expint2(n,x,ex)
    implicit none
    !  computes the n'th exponential integral of x
    !  input - x,  independent variable (-100. .le. x .le. +100.)
    !         n,  order of desired exponential integral (1 .le. n .le. 8)
    !  output - expint,  the desired result
    !           ex,  expf(-x)
    !  note   returns with e1(0)=0, (not infinity).
    !  based on the share routine nuexpi, written by j. w. cooley,
    !  courant institute of mathematical sciences, new york university
    !  obtained from rudolf loeser
    !-----general compilation of 1 august 1967.

    integer,intent(in)  :: n
    real(8),intent(in)  :: x
    real(8),intent(out) :: ex
    real(8)             :: expint2

    real(8),dimension( 7),parameter :: xint=(/1.,2.,3.,4.,5.,6.,7./)
    real(8),dimension(20),parameter :: tab = (/ & 
         0.2707662555,.2131473101,.1746297218,.1477309984, &
         1.1280843565,.1131470205,.1014028126,.0919145454,.0840790292, &
         1.0774922515,.0718735405,.0670215610,.0627878642,.0590604044, &
         1.0557529077,.0527977953,.0501413386,.0477402600,.0455592945, &
         1.0435694088/)
    real(8) :: xsave,u,xm,emx,eisave,arg,exsave,xzero,delta,y,deltax,power

    integer :: i

    xsave=0.0

    u=x
    if(u)603,602,603
602 ex=1.
    if(n-1)800,800,801
800 expint2=0.
    goto 777
801 expint2=1./xint(n-1)
    goto 777
603 if(u-xsave)604,503,604
604 xsave=u
    xm=-u
    emx=exp(xm)

    !  select method for computing ei(xm)

    if(xm-24.5)501,400,400
501 if(xm-5.)502,300,300
502 if(xm+1.)100,200,200
503 eisave=-arg
    exsave=emx

    !  now recurse to higher orders

    if(n-1)507,507,505
505 do 506 i=2,n
       eisave=(u*eisave-exsave)/(-xint(i-1))
506    continue
507    expint2=eisave
       ex=exsave
777    return

       !  ei(xm) for xm .lt. -1.0
       !  hastings polynomial approximation

100    arg=((((((u+8.573328740 )*u+18.05901697  )*u+8.634760893 )*u &
            +.2677737343)/xm)*emx)/((((u+9.573322345 )*u+25.63295615  )*u &
            +21.09965308  )*u+3.958496923 )
       goto 503
       !     ei(xm) for -1. .le. xm .lt. 5.0
       !     power series expansion about zero
200    arg=log(abs(xm))
       arg=((((((((((((((((.41159050e-14*xm+.71745406e-13)*xm+.76404637e-12) &
            *xm+.11395905e-10)*xm+.17540077e-9)*xm+.23002666e-8)*xm+.27536018e-7) &
            *xm+.30588626e-6)*xm+.31003842e-5)*xm+.28346991e-4)*xm+.23148057e-3) &
            *xm+.0016666574)*xm+.010416668)*xm+.055555572)*xm+.25)*xm+ &
            .99999999)*xm+.57721566)+arg
       goto 503

       !  ei(xm) for 5.0 .le. xm .lt. 24.5
       !  table look-up and interpolation

300    i=xm+.5
       xzero=i
       delta=xzero-xm
       arg=tab(i-4)
       if(delta)303,305,303
303    y=arg
       deltax=delta/xzero
       power=1./deltax
       do 304 i=1,7
          power=power*deltax
          y=((y-power/xzero)*delta)/xint(i)
          arg=arg+y
          if(abs(y/arg)-1.e-8)305,304,304
304       continue
305       arg=emx*arg
          goto 503
          !     ei(xm) for 24.5 .le. xm
          !     truncated continued fraction
400       arg=((((xm-15.)*xm+58.)*xm-50.)*emx)/((((xm-16.)*xm+72.)*xm-96.) &
               *xm+24.)
          goto 503
        end function expint2

        !**********************************************************************

        SUBROUTINE gauleg(x1, x2, x, w, n)
          !**********************************************************************
          ! Gives Gauss-Legendre integration ordinates and weigths, see Numerical
          ! Recipes.
          !
          !  n     : integration order       (in)
          !  x1     : lower quadrature limit (in)
          !  x2     : upper quadrature limit (in)
          !  w    : quadrature weights       (out)
          !  x    : quadrature points        (out)
          !
          !**********************************************************************
          USE mconst, only: pi
          implicit none
          integer, intent(in) :: n
          real(8), intent(in) :: x1, x2
          real(8), dimension(n), intent(out) :: x, w
          integer :: i, j, m
          real(8) :: p1, p2, p3, pp, xl, xm, z, z1
          real(8), parameter :: eps=3.d-14

          m = (n+1)/2
          xm = 0.5d0*(x2+x1)
          xl = 0.5d0*(x2-x1)

          do i=1,m

             z = cos(pi*(i-0.25d0)/(n+0.5d0))
             z1 = 0.0

             do while(abs(z-z1) .gt. eps)
                p1 = 1.0d0
                p2 = 0.0d0
                do j=1,n
                   p3 = p2
                   p2 = p1
                   p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
                end do
                pp = n*(z*p1-p2)/(z*z-1.0d0)
                z1 = z
                z = z1 - p1/pp
             end do

             x(i) = xm - xl*z
             x(n+1-i) = xm + xl*z
             w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
             w(n+1-i) = w(i)

          end do

        END SUBROUTINE gauleg

        !**********************************************************************

        !---------------------------------------------------------------------
        SUBROUTINE gausi(k,a,b,ai,xi)
          !---------------------------------------------------------------------
          !
          !  supplies quadrature weights and points for gaussian quadrature
          !  between a and b.   source of data:    lowan,  davids,
          !  levenson,  bull. amer. math. soc. 48, p 739 (1942).
          !
          !  k     : integration order, k between 1 and 10 (in)
          !  a     : lower quadrature limit                (in)
          !  b     : upper quadrature limit                (in)
          !  ai    : quadrature weights                   (out)
          !  xi    : quadrature points                    (out)
          !
          !  coded by : b. gustafsson (1974).
          !
          !---------------------------------------------------------------------
          implicit none
          !
          ! declaration of the supplied variables
          !
          integer,intent(in) :: k
          real(8), intent(in) :: a, b
          real(8), intent(out), dimension(k) :: ai, xi
          !
          ! declaration of the local variables
          !
          real(8), dimension(29) :: ap, xp
          integer, dimension(9) :: indov
          integer :: kud, k2, ioev, ined, i, ip  
          real(8) :: flk, fk
          !
          data ap/1.0d0 &
               ,0.55555555555555d0,.88888888888888d0,.347854845137d0 &
               ,0.65214515486254d0,0.23692688505618d0,0.47862867049936d0,&
               0.56888888888888d0,0.17132449237917d0,0.36076157304813d0,&
               0.46791393457269d0,0.12948496616887d0,0.27970539148927d0,&
               0.38183005050511d0,0.41795918367346d0,0.10122853629037d0,&
               0.22238103445337d0,0.31370664587788d0,0.36268378337836d0,&
               0.08127438836157d0,0.18064816069485d0,0.26061069640293d0,&
               0.31234707704000d0,0.33023935500126d0,0.06667134430868d0,&
               0.14945134915058d0,0.21908636251598d0,0.26926671930999d0,&
               0.29552422471475d0/
          data xp/ &
               0.57735026918962d0,.77459666924148d0,.0d0,0.86113631159405d0,&
               0.33998104358485d0,.90617984593866d0,.53846931010568d0,.0d0,&
               0.93246951420315d0,.66120938646626d0,.23861918608319d0,&
               0.94910791234275d0,.74153118559939d0,.40584515137739d0,.0d0,&
               0.96028985649753d0,.79666647741362d0,.52553240991632d0,&
               0.18343464249565d0,.96816023950762d0,.83603110732663d0,&
               0.61337143270059d0,.32425342340380d0,.0d0,0.97390652851717d0,&
               0.86506336668898d0,.67940956829902d0,.43339539412924d0,&
               0.14887433898163d0/
          data indov/1,3,5,8,11,15,19,24,29/
          !
          if (k.eq.1) go to 7
          !
          kud=0
          flk=float(k)/2.
          k2=k/2
          fk=float(k2)
          if (abs(flk-fk)-1.e-7) 2,1,1
1         k2=k2+1
          kud=1
2         ioev=indov(k-1)
          ined=ioev-k2
          !
          do 3 i=1,k2
             ip=ined+i
             xi(i)=-xp(ip)*(b-a)*0.5+(b+a)*0.5
3            ai(i)=(b-a)*0.5*ap(ip)
             !
             k2=k2+1
             do 4 i=k2,k
                ip=ioev+k2-i
                if (kud) 6,6,5
5               ip=ip-1
6               continue
                xi(i)=xp(ip)*(b-a)*0.5+(b+a)*0.5
4               ai(i)=(b-a)*0.5*ap(ip)
                return
                !
7               xi(1)=(b+a)*0.5
                ai(1)=b-a
                return
                !
              END SUBROUTINE gausi
              !---------------------------------------------------------------------

              !---------------------------------------------------------------------
              SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y)
                !---------------------------------------------------------------------
                !     Interpolates to obtain the value at position (x1,x2)
                !    
                !     REF : Numerical Recipes in C p.124-125
                ! 
                !---------------------------------------------------------------------

                implicit none
                !
                ! Declaration of the supplied variables
                ! 
                integer :: m, n
                real(8) :: dy, x1, x2, y, x1a(m), x2a(n), ya(m,n)
                !
                ! Declaration of the local variables
                !
                integer :: j, k
                real(8) :: ymtmp(20), yntmp(20)
                !
                do j=1,m
                   do k=1,n
                      yntmp(k)=ya(j,k)
                   end do
                   call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
                end do
                call polint(x1a,ymtmp,m,x1,y,dy)
                !
              END SUBROUTINE polin2
              !---------------------------------------------------------------------

              !---------------------------------------------------------------------
              SUBROUTINE polint(xa,ya,n,x,y,dy)
                !---------------------------------------------------------------------
                !     Interpolates to obtain the value at position (x)
                !    
                !     REF : Numerical Recipes in C p.124-125
                !
                !---------------------------------------------------------------------

                implicit none
                !
                ! Declaration of the supplied variables
                ! 
                integer :: n
                real(8) :: dy, x, y, xa(n), ya(n)
                !
                ! Declaration of the local variables
                !
                integer :: i, m, ns
                real :: den, dif, dift, ho, hp, w, c(10), d(10)
                !
                ns=1
                dif=abs(x-xa(1))
                do i=1,n
                   dift=abs(x-xa(i))
                   if (dift.lt.dif) then
                      ns=i
                      dif=dift
                   endif
                   c(i)=ya(i)
                   d(i)=ya(i)
                end do
                y=ya(ns)
                ns=ns-1
                do m=1,n-1
                   do i=1,n-m
                      ho=xa(i)-x
                      hp=xa(i+m)-x
                      w=c(i+1)-d(i)
                      den=ho-hp
                      if(den.eq.0.) stop 'failure in polint'
                      den=w/den
                      d(i)=hp*den
                      c(i)=ho*den
                   end do
                   if (2*ns.lt.n-m)then
                      dy=c(ns+1)
                   else
                      dy=d(ns)
                      ns=ns-1
                   endif
                   y=y+dy
                end do
                !
              END SUBROUTINE polint
              !---------------------------------------------------------------------

              !---------------------------------------------------------------------
              function ppvalu (break, coef, l, k, x, jderiv ,mflag)
                !---------------------------------------------------------------------
                !
                !: ppvalu 94-07-20 mats carlsson
                !:        sets value to 0 if extrapolation
                !: 
                !  from  * a practical guide to splines *  by c. de boor
                ! calls  interv
                ! calculates value at  x  of  jderiv-th derivative of pp fct from pp-repr
                !
                !******  i n p u t  ******
                !  break, coef, l, k.....forms the pp-representation of the function  f
                !        to be evaluated. specifically, the j-th derivative of  f  is
                !        given by
                !
                !     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
                !                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
                !
                !        with  h = x - break(i),  and
                !
                !       i  =  max( 1 , max( j ,  break(j) .le. x , 1 .le. j .le. l ) ).
                !
                !  x.....the point at which to evaluate.
                !  jderiv.....integer giving the order of the derivative to be evaluat-
                !        ed.  a s s u m e d  to be zero or positive.
                !
                !******  o u t p u t  ******
                !  ppvalu.....the value of the (jderiv)-th derivative of  f  at  x.
                !
                !******  m e t h o d  ******
                !     the interval index  i , appropriate for  x , is found through a
                !  call to  interv . the formula above for the  jderiv-th derivative
                !  of  f  is then evaluated (by nested multiplication).
                ! 
                !: ppvalu 95-11-06 Mats Carlsson
                !:        returns value of flag from interv
                !:
                !---------------------------------------------------------------------
                !
                implicit none
                !
                ! declaration of the supplied variables
                !
                integer :: jderiv, k, l, mflag 
                real(8)  :: x 
                real(8), dimension(l)   :: break
                real(8), dimension(k,l) :: coef
                !
                ! declaration of the local variables
                !
                real(8) :: ppvalu, fmmjdr, h
                integer :: m, i
                !
                ppvalu = 0.
                fmmjdr = k - jderiv
                !              derivatives of order  k  or higher are identically zero.
                if (fmmjdr .le. 0.)               go to 99
                !
                !              find index  i  of largest breakpoint to the left of  x .
                call interv ( break, l, x, i, mflag )
                !
                if(mflag.ne.0) then
                   ppvalu=0.0
                   return
                endif
                !
                !      evaluate  jderiv-th derivative of  i-th polynomial piece at  x .
                h = x - break(i)
                m = k
9               ppvalu = (ppvalu/fmmjdr)*h + coef(m,i)
                m = m - 1
                fmmjdr = fmmjdr - 1.
                if (fmmjdr .gt. 0.)            go to 9
99              return
              end function ppvalu
              !---------------------------------------------------------------------


              !---------------------------------------------------------------------
              SUBROUTINE interv ( xt, lxt, x, left, mflag )
                !---------------------------------------------------------------------
                !  from  * a practical guide to splines *  by c. de boor
                ! computes  left = max( i , 1 .le. i .le. lxt  .and.  xt(i) .le. x )  .
                !
                !******  i n p u t  ******
                !  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
                !  lxt.....number of terms in the sequence  xt .
                !  x.....the point whose location with respect to the sequence  xt  is
                !        to be determined.
                !
                !******  o u t p u t  ******
                !  left, mflag.....both integers, whose value is
                !
                !   1     -1      if               x .lt.  xt(1)
                !   i      0      if   xt(i)  .le. x .lt. xt(i+1)
                !  lxt     1      if  xt(lxt) .le. x
                !
                !        in particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
                !        indicates that  x  lies outside the halfopen interval
                !        xt(1) .le. y .lt. xt(lxt) . the asymmetri! treatment of the
                !        interval is due to the decision to make all pp functions cont-
                !        inuous from the right.
                !
                !******  m e t h o d  ******
                !  the program is designed to be efficient in the common situation that
                !  it is called repeatedly, with  x  taken from an increasing or decrea-
                !  sing sequence. this will happen, e.g., when a pp function is to be
                !  graphed. the first guess for  left  is therefore taken to be the val-
                !  ue returned at the previous call and stored in the  l o ! a l  varia-
                !  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
                !  essary since the present call may have nothing to do with the previ-
                !  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
                !  ilo  and are done after just three comparisons.
                !     otherwise, we repeatedly double the difference  istep = ihi - ilo
                !  while also moving  ilo  and  ihi  in the direction of  x , until
                !                      xt(ilo) .le. x .lt. xt(ihi) ,
                !  after which we use bisection to get, in addition, ilo+1 = ihi .
                !  left = ilo  is then returned.
                !
                !:
                !: interv 88-08-30 Mats Carlsson
                !:        stop if more than 20 extrapolations
                !:
                !:        88-10-26 Mats Carlsson
                !:        save statement inserted, extrapolation check changed
                !:        to work for the case x=xt
                !:
                !:        94-01-06 mats carlsson
                !:        iextr in common to enable zeroing outside routine
                !:        initialization with data statement in violation of F77 standard
                !:
                !:        94-07-20 mats carlsson
                !:        extrapolation results in zero value and warning message
                !:        only 10 warnings printed but program continues
                !:
                !:        95-08-16 mats carlsson
                !:        initialization of iextr moved to block data subprogram
                !:
                !:        95-11-27 mats carlsson
                !:        mflag=1 only set if x.gt.xt
                !:
                !---------------------------------------------------------------------
                !
                USE omp_lib
                implicit none
                !
                ! declaration of the supplied variables
                !
                integer :: lxt, left, mflag
                real(8), dimension(lxt) :: xt
                real(8)  :: x 
                !
                ! declaration of the local variables
                !
                integer :: ihi, ilo, istep, middle, iextr
                !
                common /cintv/ iextr
                !$omp threadprivate(/cintv/)
                data ilo /1/

                save ilo  

                !
                ihi = ilo + 1
                if (ihi .lt. lxt)                 go to 20
                if (x .ge. xt(lxt))            go to 110
                if (lxt .le. 1)                go to 90
                ilo = lxt - 1
                ihi = lxt
                !
20              if (x .ge. xt(ihi))               go to 40
                if (x .ge. xt(ilo))               go to 100
                !
                !              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
                istep = 1
31              ihi = ilo
                ilo = ihi - istep
                if (ilo .le. 1)                go to 35
                if (x .ge. xt(ilo))            go to 50
                istep = istep*2
                go to 31
35              ilo = 1
                if (x .lt. xt(1))                 go to 90
                go to 50
                !              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
40              istep = 1
41              ilo = ihi
                ihi = ilo + istep
                if (ihi .ge. lxt)              go to 45
                if (x .lt. xt(ihi))            go to 50
                istep = istep*2
                go to 41
45              if (x .ge. xt(lxt))               go to 110
                ihi = lxt

                !           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
50              middle = (ilo + ihi)/2
                if (middle .eq. ilo)              go to 100
                !     note. it is assumed that middle = ilo in case ihi = ilo+1 .
                if (x .lt. xt(middle))            go to 53
                ilo = middle
                go to 50
53              ihi = middle
                go to 50
                !**** set output and return.
90              mflag = -1
                left = 1
                if(iextr.le.10) then
                   write(*,92) x,xt(1)
92                 format(' interv: x outside xt, value set to zero, x,xt=',&
                        1p,2e14.6)
                endif
                iextr=iextr+1
                return
100             mflag = 0
                left = ilo
                return
110             continue
                left = lxt-1
                mflag = 0
                if(x.gt.xt(lxt)) then
                   mflag = 1
                   if(iextr.le.10) then
                      write(*,92) x,xt(lxt)
                   endif
                   iextr=iextr+1
                endif
                return
              END SUBROUTINE interv
              !---------------------------------------------------------------------

              !---------------------------------------------------------------------
              SUBROUTINE ratio(dx,xmax,numx,r)
                !---------------------------------------------------------------------
                !     Taken from sc2d.f by L.Auer and F.Paletou                 
                !     Original may be in the MAM-code
                !
                !     Uses Newton-Rapson iteration to find 'r' in the following :
                !          Sigma=dx+dx*r+dx*r^2+......+dx*r^nx=xmax
                !          
                !          Sigma=(r^nx-1)/(r-1)=xmax
                !
                ! NOTE : The convergence kriteria of this routine requires that it is
                !        compiled with the option -r8.Changing (abs(dr/r).le.1.e-9)
                !        will alter this.
                !---------------------------------------------------------------------
                implicit none
                !
                ! Declaration of the supplied variables
                !
                real(8) :: dx, xmax, r
                integer :: numx
                !
                ! Declaration the local variables
                !
                real(8) :: fr, fpr, c, dr, alph, alph1
                integer :: itr
                !
                alph=numx-2.
                alph1=alph+1.
                c=xmax/dx
                !
                ! Make a first guess
                !
                r=c**(1./alph)
                r=(c*(r-1.)+1.)**(1./alph1)
                itr=0
                !
                ! Start iterating
                !
1               fr=(r**alph1-1.)/(r-1.)-c
                fpr=(alph1*(r**alph)*(r-1.)-r**alph1+1.)/(r-1.)**2
                dr=-fr/fpr
                r=r+dr
                if(abs(dr/r).le.1.e-9) return
                itr=itr+1
                if(itr.le.10) goto 1
                !
                return
              END SUBROUTINE ratio
              !---------------------------------------------------------------------

              !---------------------------------------------------------------------
              SUBROUTINE voigt(a,vs,h)
                !---------------------------------------------------------------------
                !
                !  computes a voigt function  h = h(a,v)
                !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
                !  traving (landolt-b\rnstein, p. 449).
                !
                !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
                !
                !---------------------------------------------------------------------
                !
                implicit none
                real(8),intent(in) :: a,vs
                real(8),intent(out) :: h

                integer :: k,m,i

                real(8),dimension(5) :: a1 
                real(8) :: v,u,v2,ex,quo,h1,pqs,h1p,h2p,h3p,h4p,psi
                real(8) :: a2,u2

                real(8),dimension(19),parameter :: ak = (/ &
                     -1.12470432d0, -0.15516677d0,  3.28867591d0, -2.34357915d0,             &
                     0.42139162d0, -4.48480194d0,  9.39456063d0, -6.61487486d0,  1.98919585d0,&
                     -0.22041650d0, 0.554153432d0, 0.278711796d0,-0.188325687d0, 0.042991293d0,&
                     -0.003278278d0, 0.979895023d0,-0.962846325d0, 0.532770573d0,-0.122727278d0 /)

                real(8),parameter :: sqp = 1.772453851d0
                real(8),parameter :: sq2 = 1.414213562d0

                v = abs(vs)
                u = a + v
                v2 = v*v
                if (a.eq.0.0) go to 140
                if (a.gt.0.2) go to 120
                if (v.ge.5.0) go to 121
                !
                ex=0.
                if(v2.lt.100.)ex = exp(-v2)
                k = 1
                !
100             quo = 1.
                if (v.lt.2.4) go to 101
                quo = 1.0d0/(v2 - 1.5d0)
                m = 11
                go to 102
                !
101             m = 6
                if (v.lt.1.3) m = 1
102             do 103 i=1,5
                   a1(i) = ak(m)
                   m = m + 1
103                continue
                   h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))
                   if (k.gt.1) go to 110
                   !
                   ! a le 0.2  and v lt 5.
                   !
                   h = h1*a + ex*(1. + a*a*(1. - 2.*v2))
                   return
                   !
110                pqs = 2.0d0/sqp
                   h1p = h1 + pqs*ex
                   h2p = pqs*h1p - 2.0d0*v2*ex
                   h3p = (pqs*(1.0d0 - ex*(1.0d0 - 2.0d0*v2)) - 2.0d0*v2*h1p)/3.0d0 + pqs*h2p
                   h4p = (2.0d0*v2*v2*ex - pqs*h1p)/3.0d0 + pqs*h3p
                   psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
                   !
                   ! 0.2 lt a le 1.4  and  a + v le 3.2
                   !
                   h = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
                   return
                   !
120                if (a.gt.1.4.or.u.gt.3.2) go to 130
                   ex=0.
                   if(v2.lt.100.)ex = exp(-v2)
                   k = 2
                   go to 100
                   !
                   ! a le 0.2  and  v ge 5.
                   !
121                h = a*(15.0d0 + 6.0d0*v2 + 4.0d0*v2*v2)/(4.0d0*v2*v2*v2*sqp)
                   return
                   !
130                a2 = a*a
                   u = sq2*(a2 + v2)
                   u2 = 1.0d0/(u*u)
                   !
                   ! a gt 1.4  or  a + v gt 3.2
                   !
                   h = sq2/sqp*a/u*(1.0d0 + u2*(3.0d0*v2 - a2) + &
                        u2*u2*(15.0d0*v2*v2 - 30.0d0*v2*a2 + 3.0d0*a2*a2))
                   return
                   !
                   ! a eq 0.
                   !
140                h=0.0d0
                   if(v2.lt.100.)h=exp(-v2)

                 END SUBROUTINE voigt
                 !---------------------------------------------------------------------

                 !---------------------------------------------------------------------
                 SUBROUTINE splin0(n,dx,f,d,d2,d3,ww)
                   !---------------------------------------------------------------------
                   !
                   !  calculates standard cubic spline with continous second derivative
                   !
                   !  coded by g.scharmer, 1982
                   !
                   !  dx(k)=x(k)-x(k-1)
                   !  ww is a working array
                   !***********************************************************************
                   !

                   implicit none
                   !
                   ! Declaration of the supplied parameters
                   !
                   integer :: n 
                   real(8) :: dx(n), f(n)
                   real(8) :: d(n), d2(n), d3(n), ww(n)
                   !
                   ! Declaration of the local variables
                   !
                   real(8) :: fac, an, c2
                   integer :: k 
                   !
                   fac=-dx(2)/dx(3)
                   d(2)=(2.-fac)*(dx(2)+dx(3))
                   c2=dx(3)+fac*dx(2)
                   ww(2)=(f(2)-f(1))/dx(2)
                   ww(3)=(f(3)-f(2))/dx(3)
                   d3(2)=6.*(ww(3)-ww(2))
                   do k=3,n-1
                      d(k)=2.*(dx(k)+dx(k+1))
                      ww(k+1)=(f(k+1)-f(k))/dx(k+1)
                      fac=-dx(k)/d(k-1)
                      d(k)=d(k)+fac*dx(k)
                      if(k.eq.3) d(k)=d(k)+fac*(c2-dx(k))
                      d3(k)=6.*(ww(k+1)-ww(k))+fac*d3(k-1)
                   end do
                   fac=-dx(n)/d(n-2)
                   an=-dx(n-1)-dx(n)+fac*dx(n-1)
                   d3(n)=fac*d3(n-2)
                   fac=-an/d(n-1)
                   d2(n)=(d3(n)+fac*d3(n-1))/(dx(n-1)+fac*dx(n))
                   do k=n-1,3,-1
                      d2(k)=(d3(k)-dx(k+1)*d2(k+1))/d(k)
                   end do
                   d2(2)=(d3(2)-c2*d2(3))/d(2)
                   d2(1)=((dx(2)+dx(3))*d2(2)-dx(2)*d2(3))/dx(3)
                   do k=1,n-1
                      d(k)=ww(k+1)-(d2(k+1)+d2(k)+d2(k))*dx(k+1)/6.
                      d3(k)=(d2(k+1)-d2(k))/dx(k+1)
                   end do
                   d3(n)=d3(n-1)
                   d(n)=d(n-1)+dx(n)*(d2(n-1)+0.5*dx(n)*d3(n-1))
                   !
                 END SUBROUTINE splin0
                 !---------------------------------------------------------------------

                 !---------------------------------------------------------------------
                 SUBROUTINE splin1(n, dx, f, d, d2, d3)
                   !---------------------------------------------------------------------
                   !
                   ! Calculates qubic splines with locally determined first derivative.
                   !
                   ! Coded by AA.Nordlund, mar-83
                   !
                   !---------------------------------------------------------------------
                   !
                   implicit none
                   !
                   ! Declaration of the supplied parameters
                   !      
                   integer :: n
                   real(8) :: dx(n), f(n), d(n), d2(n), d3(n)
                   !
                   ! Local variables
                   !
                   integer :: k
                   real(8) :: dfdx, cx, cxn, dfdxn
                   !
                   ! First derivative by centered difference
                   ! 
                   do k=2, n-1
                      d(k)=(f(k+1)-f(k-1))/(dx(k+1)+dx(k))
                   end do
                   !
                   d(1)=(f(2)-f(1))/dx(2)
                   d(n)=(f(n)-f(n-1))/dx(n)
                   !
                   do k=1, n-1
                      cx=1./dx(k+1)
                      dfdx=(f(k+1)-f(k))*cx
                      d2(k)=(6.*dfdx-4.*d(k)-2.*d(k+1))*cx
                      d3(k)=6.*(d(k)+d(k+1)-2.*dfdx)*cx*cx
                   end do
                   !
                   cxn=1./dx(n)
                   dfdxn=(f(n)-f(n-1))*cxn
                   d2(n)=(4.*d(n)+2.*d(n-1)-6.*dfdxn)*cxn
                   d3(n)=d3(n-1)
                   !
                 END SUBROUTINE splin1
                 !---------------------------------------------------------------------
                 !
                 function splin(t, x, y, n, nnn)
                   USE mconst
                   !
                   ! Splin interpolation of t in a given table of points (x[i],y[i]) 
                   !
                   ! INPUT:         t        ARGUMENT   (real)
                   !                x        ARRAY OF ARGUMENTS        (real)
                   !                y        ARRAY OF function VALUES (real)
                   !                n        LENGTH OF THE ARRAYS X,Y
                   !                        if N = 2 ---> LINEAR INTERPOLATION
                   !                nnn        NUMBER OF POINTS OF THE CUBIC SPLINE
                   !                        if NNN > N   ---> NNN = N
                   !                        if NNN > 100 ---> NNN = 100
                   !                        if NNN = 2 LINEAR INTERPOLATION ADOPTED.
                   !
                   !                        if NNN IS ZERO NNN=7  (ALTERED BY P JUDGE)
                   !                       (USED TO BE if NNN OMITTED )
                   ! OUTPUT:        SPLIN INTERPOLATED value AT T
                   !

                   real(8),intent(in) :: t
                   integer,intent(in) :: n,nnn
                   real(8),dimension(n) :: x,y

                   integer :: nn,k,i,i1,i2,ishift,n1
                   real(8) :: s2(100), s3(100), dely(100), h(100)
                   real(8) :: h2(100), delsqy(100), c(100), b(100), xx(100), yy(100)
                   real(8) :: iover,int,rmin,rmax,xnorm,ynorm,eta,w,ht1,ht2,error,epsln,omega,ss2,prod,delsqs,splin
                   logical search, found
                   character*4 init
                   !        data INIT/'NOTO'/
                   !        save
                   !
                   !  nn = number of points of the spline interval
                   !
                   nn = min(7,n)
                   if (nnn .eq. 0) go to 10
                   nn = nnn
                   nn = min(n,nnn,100)
10                 error = 1.0e-12
                   if (t.ne.0) error = abs(1.0e-12*t)
                   found = search(x,t,1,n,k,error)
                   if (.not.found) goto 20
                   splin = y(k)
                   return
20                 if (n.ne.2.and.nn.ne.2) goto 30
                   !
                   !           if n=2 OR NNN=2 LINEAR INTERPOLATION IS ADOPTED
                   !
                   splin = (y(k+1)-y(k))*(t-x(k))/(x(k+1)-x(k)) + y(k)
                   return
30                 if (n.ne.1) goto 40
                   splin = y(1)
                   return
                   !
                   !           interval of the spline is computed
                   !
40                 iover = log10(float(nn))*3
                   int = (k-iover)/(nn-2*iover) + 1
                   i1 = (int-1)*(nn-2*iover)
                   i2 = i1 + nn -1
                   if (i1.ge.1) goto 50
                   i2 = nn
                   i1 = 1
50                 if (i2.le.n) goto 60
                   i2 = n
                   i1 = n - nn + 1
60                 ishift = i1 - 1
                   n1 = nn - 1
                   !
                   !           initialization check
                   !
                   if (init.ne.'okay') goto 80
                   !
                   !           interval check
                   !
                   do 70 i=1,nn
                      if (xx(i).ne.x(i+ishift)) goto 90
                      if (yy(i).ne.y(i+ishift)) goto 90
70                    continue
                      goto 200
80                    init = 'okay'
                      !
                      !           start of the computation of spline coefficients
                      !
90                    rmin = 1.0e37
                      rmax = -1.0e37
                      do 100 i=1,nn
                         if (y(i+ishift).lt.rmin)  rmin = y(i+ishift)
                         if (y(i+ishift).gt.rmax)  rmax = y(i+ishift)
100                      continue
                         !
                         !           computation of the norm factors xnorm and ynorm
                         !
                         xnorm = 1/(x(nn+ishift)-x(1+ishift))
                         ynorm = rmax - rmin
                         if (ynorm.eq.0) ynorm = rmax
                         if (ynorm.eq.0) ynorm = 1
                         ynorm = 1 / ynorm
                         do 110 i=1,n1
                            h(i) = (x(i+1+ishift)-x(i+ishift))*xnorm
                            dely(i) = ((y(i+1+ishift)-y(i+ishift))/h(i))*ynorm
110                         continue
                            do 120 i=2,n1
                               h2(i) = h(i-1) + h(i)
                               b(i) = 0.5*h(i-1)/h2(i)
                               delsqy(i) = (dely(i)-dely(i-1))/h2(i)
                               s2(i) = 2.*delsqy(i)
                               c(i) = 3.*delsqy(i)
120                            continue
                               s2(1) = 0
                               s2(n1+1) = 0
                               !
                               !        solution of the linear system of the spline coefficients
                               !        by successive overrelaxation.
                               !        constants: eps = error criterion in the iterative solution.
                               !                   omega = relaxation coefficient.
                               !
                               epsln = 1.0e-8
                               omega = 1.0717968
130                            eta = 0.
                               do 160 i=2,n1
                                  w = (c(i)-b(i)*s2(i-1)-(0.5-b(i))*s2(i+1)-s2(i))*omega
                                  if (abs(w)-eta) 150, 150, 140
140                               eta = abs(w)
150                               s2(i) = s2(i) + w
160                               continue
                                  if (eta-epsln) 170, 130, 130
170                               do 180 i=1,n1
                                     s3(i) = (s2(i+1)-s2(i))/h(i)
180                                  continue
                                     !        x and y stored in xx and yy for interval check
                                     !
                                     do 190 i=1,nn
                                        xx(i) = x(i+ishift)
                                        yy(i) = y(i+ishift)
190                                     continue
200                                     i = k - i1 + 1
                                        ht1 = (t-x(i+ishift))*xnorm
                                        ht2 = (t-x(i+ishift+1) )*xnorm
                                        if ((t-x(1))*xnorm.gt.0) goto 210
                                        !
                                        !        extrapolation   t < X(1)
                                        !
                                        splin = y(1)
                                        return
210                                     if ((t-x(n))*xnorm.gt.0) goto 220
                                        !
                                        !        interpolation by means of the spline coefficients
                                        !
                                        ss2 = s2(i) + ht1*s3(i)
                                        prod = ht1*ht2
                                        delsqs = (s2(i)+s2(i+1)+ss2)/6.
                                        splin = y(i+ishift) + (ht1*dely(i)+prod*delsqs)/ynorm
                                        return
                                        !
                                        !        extrapolation   t > x(n)
                                        !
220                                     splin =  y(n)
                                        return
                                      end function splin
                                      !
                                      !---------------------------------------------------------------------


                                      SUBROUTINE tautsp ( tau, gtau, ntau, gamma, s, break, coef, l, k, iflag )
                                        !---------------------------------------------------------------------
                                        ! from  * a practical guide to splines *  by c. de boor
                                        ! constructs cubic spline interpolant to given data
                                        !        tau(i), gtau(i), i=1,...,ntau.
                                        ! if  gamma .gt. 0., additional knots are introduced where needed to
                                        ! make the interpolant more flexible locally. this avoids extraneous
                                        ! inflection points typical of cubic spline interpolation at knots to
                                        ! rapidly changing data.
                                        !
                                        ! parameters
                                        !           input
                                        ! tau      sequence of data points. must be strictly increasing.
                                        ! gtau     corresponding sequence of function values.
                                        ! ntau     number of data points. must be at least  4 .
                                        ! gamma    indicates whether additional flexibility is desired.
                                        !         = 0., no additional knots
                                        !         in (0.,3.), under certain conditions on the given data at
                                        !               points i-1, i, i+1, and i+2, a knot is added in the
                                        !               i-th interval, i=2,...,ntau-2. see description of meth-
                                        !               od below. the interpolant gets rounded with increasing
                                        !               gamma. a value of  2.5  for gamma is typical.
                                        !         in (3.,6.), same , except that knots might also be added in
                                        !               intervals in which an inflection point would be permit-
                                        !               ted.  a value of  5.5  for gamma is typical.
                                        !           output
                                        ! break, coef, l, k  give the pp-representation of the interpolant.
                                        !         specifically, for break(i) .le. x .le. break(i+1), the
                                        !       interpolant has the form
                                        ! f(x) = coef(1,i) +dx(coef(2,i) +(dx/2)(coef(3,i) +(dx/3)coef(4,i)))
                                        !       with  dx = x - break(i) and i=1,...,l .
                                        ! iflag   = 1, ok
                                        !         = 2, input was incorrect. a printout specifying the mistake
                                        !           was made.
                                        !           workspace
                                        ! s     is required, of size (ntau,6). the individual columns of this
                                        !       array contain the following quantities mentioned in the write-
                                        !       up and below.
                                        !    s(.,1) = dtau = tau(.+1) - tau
                                        !    s(.,2) = diag = diagonal in linear system
                                        !    s(.,3) = u = upper diagonal in linear system
                                        !    s(.,4) = r = right side for linear system (initially)
                                        !           = fsecnd = solution of linear system , namely the second
                                        !                      derivatives of interpolant at  tau
                                        !    s(.,5) = z = indicator of additional knots
                                        !    s(.,6) = 1/hsecnd(1,x) with x = z or = 1-z. see below.
                                        !
                                        ! ------  m e t h o d  ------
                                        ! on the i-th interval, (tau(i), tau(i+1)), the interpolant is of the
                                        ! form
                                        ! (*)  f(u(x)) = a + b*u + c*h(u,z) + d*h(1-u,1-z) ,
                                        ! with  u = u(x) = (x - tau(i))/dtau(i). here,
                                        !      z = z(i) = addg(i+1)/(addg(i) + addg(i+1))
                                        ! (= .5, in case the denominator vanishes). with
                                        !      addg(j) = abs(ddg(j)), ddg(j) = dg(j+1) - dg(j),
                                        !      dg(j) = divdif(j) = (gtau(j+1) - gtau(j))/dtau(j)
                                        ! and
                                        !      h(u,z) = alpha*u**3 + (1 - alpha)*(max(((u-zeta)/(1-zeta)),0)**3
                                        ! with
                                        !      alpha(z) = (1-gamma/3)/zeta
                                        !      zeta(z) = 1 - gamma*min((1 - z), 1/3)
                                        ! thus, for 1/3 .le. z .le. 2/3,  f  is just a cubic polynomial on
                                        ! the interval i. otherwise, it has one additional knot, at
                                        !        tau(i) + zeta*dtau(i) .
                                        ! as  z  approaches  1, h(.,z) has an increasingly sharp bend  near 1,
                                        ! thus allowing  f  to turn rapidly near the additional knot.
                                        !    in terms of f(j) = gtau(j) and
                                        !      fsecnd(j) = 2.derivative of  f  at  tau(j),
                                        ! the coefficients for (*) are given as
                                        !      a = f(i) - d
                                        !      b = (f(i+1) - f(i)) - (c - d)
                                        !      c = fsecnd(i+1)*dtau(i)**2/hsecnd(1,z)
                                        !      d = fsecnd(i)*dtau(i)**2/hsecnd(1,1-z)
                                        ! hence can be computed once fsecnd(i),i=1,...,ntau, is fixed.
                                        !  f  is automatically continuous and has a continuous second derivat-
                                        ! ive (except when z = 0 or 1 for some i). we determine fscnd(.) from
                                        ! the requirement that also the first derivative of  f  be contin-
                                        ! uous. in addition, we require that the third derivative be continuous
                                        ! across  tau(2) and across  tau(ntau-1) . this leads to a strictly
                                        ! diagonally dominant tridiagonal linear system for the fsecnd(i)
                                        ! which we solve by gauss elimination without pivoting.
                                        !
                                        !:
                                        !: tautsp 89-08-30 mats carlsson
                                        !:        test of onemzt changed from 0. to 1.e-11 to avoid division
                                        !:        by underflowed expression (onemzt**3)
                                        !:
                                        !:        89-09-05 mats carlsson
                                        !:        dimension of break and coef changed from 1 to * to make
                                        !:        possible check of index out of bounds
                                        !:
                                        !:        98-01-20 mats carlsson
                                        !:        statement function alph eliminated
                                        !:
                                        !---------------------------------------------------------------------
                                        implicit none
                                        integer ::  iflag,k,l,ntau,i,method,ntaum1
                                        real(8) :: gamma
                                        real(8),dimension(ntau) :: break,tau,gtau
                                        real(8),dimension(4,ntau) :: coef
                                        real(8),dimension(ntau,6) :: s
                                        real(8),parameter :: one=1.0
                                        real(8) :: sixth, gam, z, temp, denom,zeta, onemzt,zt2,alpha,onemg3
                                        real(8) :: factor, entry3,factr2, ratio, entry, divdif,c,d,del

                                        !     alph(x,onemg3) = min(one,onemg3/x)
                                        !
                                        ! there must be at least  4  interpolation points.
                                        if (ntau .ge. 4)                  go to 5
                                        print 600,ntau
600                                     format('ntau = ',i4,'. should be .ge. 4')
                                        go to 999
                                        !
                                        ! construct delta tau and first and second (divided) differences of data
                                        !
5                                       ntaum1 = ntau - 1
                                        do 6 i=1,ntaum1
                                           s(i,1) = tau(i+1) - tau(i)
                                           if (s(i,1) .gt. 0.0)            go to 6
                                           print 610,i,tau(i),tau(i+1)
610                                        format(7h point ,i3,13h and the next,2e15.6,15h are disordered)
                                           go to 999
6                                          s(i+1,4) = (gtau(i+1)-gtau(i))/s(i,1)
                                           do 7 i=2,ntaum1
7                                             s(i,4) = s(i+1,4) - s(i,4)
                                              !
                                              ! construct system of equations for second derivatives at  tau. at each
                                              ! interior data point, there is one continuity equation, at the first
                                              ! and the last interior data point there is an additional one for a
                                              ! total of  ntau  equations in  ntau  unknowns.
                                              !
                                              i = 2
                                              s(2,2) = s(1,1)/3.0
                                              sixth = 1.0/6.0
                                              method = 2
                                              gam = gamma
                                              if (gam .le. 0.0)   method = 1
                                              if ( gam .le. 3.0)                 go to 9
                                              method = 3
                                              gam = gam - 3.0
9                                             onemg3 = 1.0 - gam/3.0
                                              !                ------ loop over i ------
10                                            continue
                                              !         construct z(i) and zeta(i)
                                              z = .50
                                              go to (19,11,12),method
11                                            if (s(i,4)*s(i+1,4) .lt. 0.0)      go to 19
12                                            temp = abs(s(i+1,4))
                                              denom = abs(s(i,4)) + temp
                                              if (denom .eq. 0.0)                go to 19
                                              z = temp/denom
                                              if (abs(z - .50) .le. sixth)  z = .50
19                                            s(i,5) = z
                                              !  ******set up part of the i-th equation which depends on
                                              !        the i-th interval
                                              if (z - .50)                       21,22,23
21                                            zeta = gam*z
                                              onemzt = 1.0 - zeta
                                              zt2 = zeta**2
                                              alpha = min(one,onemg3/onemzt)
                                              factor = zeta/(alpha*(zt2-1.0) + 1.0)
                                              s(i,6) = zeta*factor/6.0
                                              s(i,2) = s(i,2) + s(i,1)*((1.0-alpha*onemzt)*factor/ &
                                                   2.0 - s(i,6))
                                              !    if z = 0 and the previous z = 1, then d(i) = 0. since then
                                              !    also u(i-1) = l(i+1) = 0, its value does not matter. reset
                                              !    d(i) = 1 to insure nonzero pivot in elimination.
                                              if (s(i,2) .le. 0.0) s(i,2) = 1.0
                                              s(i,3) = s(i,1)/6.0
                                              go to 25
22                                            s(i,2) = s(i,2) + s(i,1)/3.0
                                              s(i,3) = s(i,1)/6.0
                                              go to 25
23                                            onemzt = gam*(1.0 - z)
                                              zeta = 1.0 - onemzt
                                              alpha = min(one,onemg3/zeta)
                                              factor = onemzt/(1.0 - alpha*zeta*(1.0+onemzt))
                                              s(i,6) = onemzt*factor/6.0
                                              s(i,2) = s(i,2) + s(i,1)/3.0
                                              s(i,3) = s(i,6)*s(i,1)
25                                            if (i .gt. 2)                     go to 30
                                              s(1,5) = .50
                                              ! ******the first two equations enforce continuity of the first and of
                                              !       the third derivative across tau(2).
                                              s(1,2) = s(1,1)/6.0
                                              s(1,3) = s(2,2)
                                              entry3 = s(2,3)
                                              if (z - .50)                       26,27,28
26                                            factr2 = zeta*(alpha*(zt2-1.0) + 1.0)/(alpha*(zeta* &
                                                   zt2-1.0)+1.0)
                                              ratio = factr2*s(2,1)/s(1,2)
                                              s(2,2) = factr2*s(2,1) + s(1,1)
                                              s(2,3) = -factr2*s(1,1)
                                              go to 29
27                                            ratio = s(2,1)/s(1,2)
                                              s(2,2) = s(2,1) + s(1,1)
                                              s(2,3) = -s(1,1)
                                              go to 29
28                                            ratio = s(2,1)/s(1,2)
                                              s(2,2) = s(2,1) + s(1,1)
                                              s(2,3) = -s(1,1)*6.0*alpha*s(2,6)
                                              !      at this point, the first two equations read
                                              !             diag(1)*x1 + u(1)*x2 + entry3*x3 = r(2)
                                              !      -ratio*diag(1)*x1 + diag(2)*x2 + u(2)*x3 = 0.
                                              !      eliminate first unknown from second equation
29                                            s(2,2) = ratio*s(1,3) + s(2,2)
                                              s(2,3) = ratio*entry3 + s(2,3)
                                              s(1,4) = s(2,4)
                                              s(2,4) = ratio*s(1,4)
                                              go to 35
30                                            continue
                                              ! ******the i-th equation enforces continuity of the first derivative
                                              !       across tau(i). it has been set up in statements 35 up to 40
                                              !       and 21 up to 25 and reads now
                                              !        -ratio*diag(i-1)*xi-1 + diag(i)*xi + u(i)*xi+1 = r(i) .
                                              !       eliminate (i-1)st unknown from this equation
                                              s(i,2) = ratio*s(i-1,3) + s(i,2)
                                              s(i,4) = ratio*s(i-1,4) + s(i,4)
                                              !
                                              ! ******set up the part of the next equation which depends on the
                                              !       i-th interval.
35                                            if (z - .50)                       36,37,38
36                                            ratio = -s(i,6)*s(i,1)/s(i,2)
                                              s(i+1,2) = s(i,1)/3.0
                                              go to 40
37                                            ratio = -(s(i,1)/6.0)/s(i,2)
                                              s(i+1,2) = s(i,1)/3.0
                                              go to 40
38                                            ratio = -(s(i,1)/6.0)/s(i,2)
                                              s(i+1,2) = s(i,1)*((1.0-zeta*alpha)*factor/2.0 - s(i,6))
                                              !        ------  end of i loop ------
40                                            i = i+1
                                              if (i .lt. ntaum1)                go to 10
                                              s(i,5) = .50
                                              !
                                              !       ------  last two equations  ------
                                              ! the last two equations enforce continuity of third derivative and
                                              ! of first derivative across  tau(ntau-1).
                                              entry = ratio*s(i-1,3) + s(i,2) + s(i,1)/3.0
                                              s(i+1,2) = s(i,1)/6.0
                                              s(i+1,4) = ratio*s(i-1,4) + s(i,4)
                                              if (z - .50)                       41,42,43
41                                            ratio = s(i,1)*6.0*s(i-1,6)*alpha/s(i-1,2)
                                              s(i,2) = ratio*s(i-1,3) + s(i,1) + s(i-1,1)
                                              s(i,3) = -s(i-1,1)
                                              go to 45
42                                            ratio = s(i,1)/s(i-1,2)
                                              s(i,2) = ratio*s(i-1,3) + s(i,1) + s(i-1,1)
                                              s(i,3) = -s(i-1,1)
                                              go to 45
43                                            factr2 = onemzt*(alpha*(onemzt**2-1.0)+1.0)/ &
                                                   (alpha*(onemzt**3-1.0)+1.0)
                                              ratio = factr2*s(i,1)/s(i-1,2)
                                              s(i,2) = ratio*s(i-1,3) + factr2*s(i-1,1) + s(i,1)
                                              s(i,3) = -factr2*s(i-1,1)
                                              !    at this point, the last two equations read
                                              !            diag(i)*xi + u(i)*xi+1 = r(i)
                                              !     -ratio*diag(i)*xi + diag(i+1)*xi+1 = r(i+1)
                                              !    eliminate xi from last equation
45                                            s(i,4) = ratio*s(i-1,4)
                                              ratio = -entry/s(i,2)
                                              s(i+1,2) = ratio*s(i,3) + s(i+1,2)
                                              s(i+1,4) = ratio*s(i,4) + s(i+1,4)
                                              !
                                              !       ------ back substitution ------
                                              !
                                              s(ntau,4) = s(ntau,4)/s(ntau,2)
50                                            s(i,4) = (s(i,4) - s(i,3)*s(i+1,4))/s(i,2)
                                              i = i - 1
                                              if (i .gt. 1)                  go to 50
                                              s(1,4) = (s(1,4)-s(1,3)*s(2,4)-entry3*s(3,4))/s(1,2)
                                              !
                                              !      ------ construct polynomial pieces ------
                                              !
                                              break(1) = tau(1)
                                              l = 1
                                              do 70 i=1,ntaum1
                                                 coef(1,l) = gtau(i)
                                                 coef(3,l) = s(i,4)
                                                 divdif = (gtau(i+1)-gtau(i))/s(i,1)
                                                 z = s(i,5)
                                                 if (z - .50)                    61,62,63
61                                               if (z .eq. 0.0)                 go to 65
                                                 zeta = gam*z
                                                 onemzt = 1.0 - zeta
                                                 c = s(i+1,4)/6.0
                                                 d = s(i,4)*s(i,6)
                                                 l = l + 1
                                                 del = zeta*s(i,1)
                                                 break(l) = tau(i) + del
                                                 zt2 = zeta**2
                                                 alpha = min(one,onemg3/onemzt)
                                                 factor = onemzt**2*alpha
                                                 coef(1,l) = gtau(i) + divdif*del &
                                                      + s(i,1)**2*(d*onemzt*(factor-1.0)+c*zeta* &
                                                      (zt2-1.0))
                                                 coef(2,l) = divdif + s(i,1)*(d*(1.0-3.0*factor)+c*(3.0* &
                                                      zt2-1.0))
                                                 coef(3,l) = 6.0*(d*alpha*onemzt + c*zeta)
                                                 coef(4,l) = 6.0*(c - d*alpha)/s(i,1)
                                                 coef(4,l-1) = coef(4,l) - 6.0*d*(1.0-alpha)/(del*zt2)
                                                 coef(2,l-1) = coef(2,l) - del*(coef(3,l) -(del/2.0)* &
                                                      coef(4,l-1))
                                                 go to 68
62                                               coef(2,l) = divdif - s(i,1)*(2.0*s(i,4) + s(i+1,4))/6.0
                                                 coef(4,l) = (s(i+1,4)-s(i,4))/s(i,1)
                                                 go to 68
63                                               onemzt = gam*(1.0 - z)
                                                 if (abs(onemzt) .lt. 1.d-11)            go to 65
                                                 zeta = 1.0 - onemzt
                                                 alpha = min(one,onemg3/zeta)
                                                 c = s(i+1,4)*s(i,6)
                                                 d = s(i,4)/6.0
                                                 del = zeta*s(i,1)
                                                 break(l+1) = tau(i) + del
                                                 coef(2,l) = divdif - s(i,1)*(2.0*d + c)
                                                 coef(4,l) = 6.0*(c*alpha - d)/s(i,1)
                                                 l = l + 1
                                                 coef(4,l) = coef(4,l-1) + 6.0*(1.0-alpha)*c/(s(i,1)* &
                                                      onemzt**3)
                                                 coef(3,l) = coef(3,l-1) + del*coef(4,l-1)
                                                 coef(2,l) = coef(2,l-1)+del*(coef(3,l-1)+(del/2.0)* &
                                                      coef(4,l-1))
                                                 coef(1,l) = coef(1,l-1)+del*(coef(2,l-1)+(del/2.0)* &
                                                      (coef(3,l-1) &
                                                      +(del/3.0)*coef(4,l-1)))
                                                 go to 68
65                                               coef(2,l) = divdif
                                                 coef(3,l) = 0.0
                                                 coef(4,l) = 0.0
68                                               l = l + 1
70                                               break(l) = tau(i+1)
                                                 !     l = l - 1   deletion of this statement makes extrapol warning possible
                                                 k = 4
                                                 iflag = 1
                                                 return
999                                              iflag = 2
                                                 return

                                               END SUBROUTINE tautsp

                                               !---------------------------------------------------------------------

                                               SUBROUTINE interpol(y1,x1,y2,x2)
                                                 implicit none
                                                 real(8),dimension(:),intent(in)  :: y1,x1,x2
                                                 real(8),dimension(:),intent(out) :: y2

                                                 ! work space for tautsp
                                                 integer :: iflag,ntmp,i,k,n
                                                 real(8),dimension(2*size(y1),6         ) :: work
                                                 real(8),dimension(4         ,2*size(y1)) ::coef
                                                 real(8),dimension(2*size(y1)           ) :: break

                                                 real(8),dimension(size(y1)) :: yy1, xx1
                                                 real(8) :: gamma=5.5d0


                                                 n=size(y1)
                                                 ! x axis must be ascending
                                                 if (x1(2) < x1(1)) then
                                                    do i=1,n
                                                       yy1(i)=y1(n+1-i)
                                                    end do
                                                    do i=1,n
                                                       xx1(i)=x1(n+1-i)
                                                    end do
                                                 else
                                                    yy1=y1
                                                    xx1=x1
                                                 endif

                                                 call tautsp(xx1,yy1,n,gamma,work,break,coef,ntmp,k,iflag)
                                                 do i=1,n
                                                    y2(i)= ppvalu(break(1:ntmp),coef(1:4,1:ntmp),ntmp,4,x2(i),0,iflag)
                                                 end do

                                               END SUBROUTINE interpol
                                               !---------------------------------------------------------------------

! ES
! ES 18.08.2020
! transfered from mul23lus
REAL(8) FUNCTION PROFLL(ALAMD,ALAMB,GA,GW,GQ,ION,PION,CHIEV,CHJEV, AWGT,TOTHI,QNORM, T, N_E, toth2in)
    REAL(8) :: ALAMD,ALAMB,GA,GW,GQ,PION,CHIEV,CHJEV,AWGT,TOTHI,QNORM, T, N_E, VEL
    REAL(8) :: GV, GR, GS, SIGMA, XALPHA, GX, GAMMAF, VBAR, GAM_GW, C625, GSLG
    INTEGER(4) :: ION, ZZ
    real(8),optional :: toth2in
    real(8) :: toth2
! some constants
    REAL(8) :: EE=1.602189E-12,HH=6.626176E-27, CC=2.99792458E10, EM=9.109534E-28, UU=1.6605655E-24, BK=1.380662E-16, PI=3.14159265359
    REAL(8) :: DNYD, DOP, ADAMP, Q, V, PHI, H, HNY4P, BIJ
    toth2 = 0.0
    if (present(toth2in)) toth2 = toth2in
! N_E      electron density at a point in the atmosphere
! T       temperature at a point in the atmosphere
!  CALCULATE DOPPLER WIDTH
    DNYD=SQRT(2.*BK*T/AWGT)*1.E-5/QNORM
! NO MICROTURBULENCE INCLUDED
! GW as given in the model atom
      IF(GW.GE.20.0) THEN
        GV=GW
        SIGMA=INT(GV)*2.80028E-21
        XALPHA=GV-INT(GV)
        GX=2.-XALPHA*.5
        GX=GX-1.0
        GAMMAF=1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX)*GX)*GX)*GX)*GX
        GV=(4./PI)**(XALPHA*0.5)*GAMMAF*1.E4*SIGMA
        VBAR=SQRT(21173.*T*(1./1.008+1./(AWGT/UU)))
        GV=GV*((VBAR/1.E4)**(1.-XALPHA))
!       TOTHI    NUMBER DENSITY OF NEUTRAL HYDROGEN
        GV=GV*TOTHI*1.E6*2.
      ! MOD PE: only for positive values below 20
      ELSEIF(GW.GT.0.0) THEN
        ZZ=ION
!       PION     IONIZATION POTENTIAL FROM GROUND STAGE (ELECTRON VOLTS)
!       CHJEV    EXCITATION POTENTIAL OF UPPER LEVEL (ELECTRON VOLTS)
        C625=1.283984E-12*ZZ**0.8*(ABS(1./(PION-CHJEV)**2-1./(PION-CHIEV)**2))**.4
        GV=GW*8.411*(8.*BK*T/PI*(1./(1.008*UU)+1./AWGT))**0.3*TOTHI*C625
      ! PE 15-09-21:
      ! Add van der Waals broadening according to Turbospectrum 
      ! if the damping constant is negative
      elseif (GW.lt.0.0) then
         ! this is copied from turobospectrum2019 Bpz 05/04-2013
         GV = 10.0**GW * (T/10000.)**0.3 * tothi
         GV = GV * (1.0 + toth2/tothi*0.85) ! Missing He part, set it to 0 for now
      ELSE
        GV=0.0
      ENDIF
      GR=GA
      IF(GQ.GE.0.0) THEN
        GS=GQ*N_E
      ELSE
!       FORMULA FROM GRAY, P.237
        GSLG=19.4+2./3.*GQ+LOG10(N_E*BK)+LOG10(T)/6.
        GS=10.**GSLG
      ENDIF
!
      GAM_GW=GR+GV+GS

!  DNYD*QNORM IS IN UNITS OF KM/S
!  DOP        IS IN UNITS HZ

      DOP=DNYD*QNORM/ALAMB*1.E13
      ADAMP=GAM_GW/(4.*PI*DOP)
!
!  Q IS DELTA(NY) IN UNITS OF A TYPICAL DOPPLER WIDTH
!  V IS DELTA(NY)/DNYD, UNITLESS
!
      Q = (CC/ALAMD - CC/ALAMB)*ALAMB*1.E-5/QNORM
! ES VEL =0
! macroturbulence set to zero for now
! the same way as in mulv2.3, 1D version
      V=(Q-0.)/DNYD
!     
      CALL VOIGT(ADAMP,V,H)
!      print*, 'H=',  H
 !---------------------------------------------------------------------
 !
 !  computes a voigt function  h = h(a,v)
 !  a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
 !  traving (landolt-b\rnstein, p. 449).
 !
 !  coded by: kjell eriksson  (apr-1974). revised: (july-1975).
 !
 !---------
      PHI=H/(DNYD*1.772453851)

      HNY4P=HH*CC/QNORM/4./PI*1.E-5
      BIJ=1.00*ALAMB/(2.*HH*CC*1.E24)*6.671E15
!
      PROFLL=HNY4P*BIJ*PHI
      RETURN
end


END MODULE mmath
! *******************************************************************j
