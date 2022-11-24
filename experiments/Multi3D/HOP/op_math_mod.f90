MODULE op_math_mod
   USE prec_mod
   implicit none

   interface poly1d
       module procedure poly3, poly4, poly5, poly6
   end interface poly1d
    
contains

   FUNCTION interp1d(xx,yy,newx, order, sorted) result(newy)
         real(dp), dimension(:), intent(in) :: xx, yy, newx
         real(dp), dimension(:), allocatable :: newy
         integer, value :: order
         logical, value :: sorted
         integer :: i, j, nx, ny
         procedure(interp1d_linear), pointer :: interp

         nx = size(xx)
         ny = size(yy)

         if (nx /= ny) then
            print *, 'X array is of length', nx
            print *, 'but Y array is of lenght', ny
            stop
         endif

         interp => interp1d_linear

         !if (.not. present(sorted)) sorted = .false.
         !if (.not. present(order))  order  = 1

         if (order == 2) then
            if (nx < 3) then
               print *, 'too few points for quadratic interpolation'
               stop
            else
               interp => interp1d_quadratic
            endif
         elseif (order > 2) then
            print *, 'only linear and quadratic are supported'
            stop
         endif

         allocate(newy, mold=newx)

         if (sorted) then
            j = 1
            i = 1

            do while (i <= size(newy))
               if (j >= nx) then
                  newy(i) = yy(ny)

               elseif (xx(j+1) < newx(i)) then
                  j = j+1
                  cycle

               elseif (xx(j) < newx(i)) then
                  newy(i) = interp(j,xx,yy,newx(i))
               else
                  newy(i) = yy(j)
               endif

               i = i+1
            end do
         else
            do i=1, size(newy)
               ! xx and yy are still sorted
               j = searchsorted(xx, newx(i))
               newy(i) = interp(j,xx,yy,newx(i))
            end do
         endif

   end FUNCTION interp1d

   pure FUNCTION interp1d_linear(j,xx,yy,newx) result(newy)
      real(dp), dimension(:), intent(in) :: xx, yy
      real(dp), intent(in) :: newx
      integer, intent(in) :: j
      real(dp) :: m, newy

      m = (yy(j+1) - yy(j)) / (xx(j+1) - xx(j))
      newy = m * (newx - xx(j)) + yy(j)
   end FUNCTION interp1d_linear

   pure FUNCTION interp1d_quadratic(j,xx,yy,newx) result(newy)
      real(dp), dimension(:), intent(in) :: xx, yy
      real(dp), intent(in) :: newx
      integer, intent(in) :: j
      real(dp) :: newy, w(3)
      integer :: j1, j2, j3

      if (size(xx) - j > 1) then
         j1 = j
         j2 = j+1
         j3 = j+2
      else
         j1 = j-1
         j2 = j
         j3 = j+1
      endif

      w(1) = (newx-xx(j2)) * (newx-xx(j3)) / ((xx(j1)-xx(j2)) * (xx(j1)-xx(j3)))
      w(2) = (newx-xx(j1)) * (newx-xx(j3)) / ((xx(j2)-xx(j1)) * (xx(j2)-xx(j3)))
      w(3) = (newx-xx(j1)) * (newx-xx(j2)) / ((xx(j3)-xx(j1)) * (xx(j3)-xx(j2)))

      newy = sum(yy(j1:j3)*w(1:3))

   end FUNCTION interp1d_quadratic

   ! Finds index of element after which value should be inserted to maintain order.
   pure FUNCTION searchsorted(arr, val, outbounds_ok) result(i)
      real(8), intent(in) :: arr(:), val
      logical, optional, intent(in) :: outbounds_ok
      logical :: ob_ok, above, below
      integer :: i, n

      if (present(outbounds_ok)) then
         ob_ok = outbounds_ok
      else
         ob_ok = .false.
      endif

      n = size(arr)
      below = val < arr(1)
      above = val > arr(n)

      if (below .or. above) then
         if (ob_ok) then
            if (below) i = 0
            if (above) i = n
         else
            i = -1
         endif

      else
         i = BinarySearch(arr, 1, n, val)
      
      endif
   end FUNCTION searchsorted

   pure recursive FUNCTION BinarySearch(arr, l, r, x) result(idx)
      real(8), intent(in) :: arr(:)
      real(8), intent(in) :: x
      integer, intent(in) :: l, r
      integer :: idx, mid

      if (r - l > 1) then
          mid = l + (r - l) / 2

          if(arr(mid) >= x) then
              idx = BinarySearch(arr, l, mid, x)

          elseif(arr(mid) < x) then
              idx = BinarySearch(arr, mid, r, x)

          endif
      else
          idx = l
      endif

   end FUNCTION BinarySearch

   elemental real(dp) FUNCTION poly3(x, a, b, c) result(poly)
      real(dp), intent(in) :: x, a, b, c
      poly = a*x*x + b*x + c
   end FUNCTION poly3

   elemental real(dp) FUNCTION poly4(x, a, b, c, d) result(poly)
      real(dp), intent(in) :: x, a, b, c, d
      poly = x*poly3(x, a, b, c) + d
   end FUNCTION poly4

   elemental real(dp) FUNCTION poly5(x, a, b, c, d, e) result(poly)
      real(dp), intent(in) :: x, a, b, c, d, e
      poly = x*poly4(x, a, b, c, d) + e
   end FUNCTION poly5

   elemental real(dp) FUNCTION poly6(x, a, b, c, d, e, f) result(poly)
      real(dp), intent(in) :: x, a, b, c, d, e, f
      poly = x*poly5(x, a, b, c, d, e) + f
   end FUNCTION poly6


!   elemental real(dp) FUNCTION solve_quadratic(p, q, avoid_negative) result(x)
!      real(dp), intent(in) :: p, q
!      logical, intent(in), optional :: avoid_negative
!      real(dp) :: cam, root
!
!      cam = - 0.5*p
!      root = sqrt(abs(0.25*p*p - q))
!
!      x = cam + root
!
!      if (present(avoid_negative)) then
!         if (avoid_negative .and. x<0) x = cam - root
!      endif
!   end FUNCTION solve_quadratic

   !===============================================================================
   elemental real(dp) FUNCTION solve_quadratic(p, q, avoid_negative) result(x)
      real(dp), intent(in) :: p, q
      logical, intent(in), optional :: avoid_negative
      real(dp) :: cam, root
      !----------------------------------------------------------------------------
      if (p > 0) then
         root = sqrt(abs(p*p - 4.d0*q))
         x = - 2.d0*q / (p + root)
      else
         cam = - 0.5d0*p
         root = sqrt(abs(0.25d0*p*p - q))
         x = cam + root
      endif
      if (x < 0) then
         if (present(avoid_negative)) then
            if (avoid_negative) x = cam - root
         endif
      endif
   end FUNCTION solve_quadratic
    
   !===============================================================================
   !  finds solution of system of linear equations
   !  with gaussian elimination with pivoting
   !: lineq  90-06-05  new routine: (martin j stift)
   !===============================================================================
   SUBROUTINE lineq(a,b,n)
      integer :: n
      real(dp), dimension(n,n) :: a
      real(dp), dimension(n)   :: b

      real(dp), dimension(n)  :: c
      integer, dimension(n) :: icol
      integer :: i, ibeg, jbeg, ima, jma, j, it, n1, jmin, i1, imin
      real(dp)  :: ama, temp, quot
      !-----------------------------------------------------------------------------
      ! initialize column count and starting points

      do i = 1,n
         icol(i) = i
      end do

      ibeg = 1
      jbeg = 1

      ! determine pivot
      do while (ibeg.lt.n)
         ama = abs(a(ibeg,jbeg))
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
      end do
      !if (ibeg.lt.n) goto 20

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

   end SUBROUTINE lineq
   !---------------------------------------------------------------------

   !===============================================================================
   SUBROUTINE lineq_no_pivot (a, b, n)
      real(8):: a(n,n), b(n)
      integer :: i, j, k, n
      !-----------------------------------------------------------------------------
      do k=1,n-1
        do i=k+1,n
          a(i,k) = a(i,k)/a(k,k)
        end do
        do j=k+1,n
          if (a(k,j) /= 0) then
            do i=k+1,n
              a(i,j) = a(i,j) - a(i,k)*a(k,j)
            end do
          endif
        end do
      end do
      do i=1,n
        do j=1,i-1
          b(i) = b(i) - a(i,j)*b(j)
        end do
      end do
      do i=n,1,-1
        do j=i+1,n
          b(i) = b(i) - a(i,j)*b(j)
        end do
        b(i) = b(i)/a(i,i)
      end do
   END SUBROUTINE lineq_no_pivot

   ! **************************************************************************
   SUBROUTINE AINV2(A,M)
      implicit none
      integer :: M
      real(8), dimension(M,M) :: A
      integer, dimension(M) :: IND
      real(8), dimension(M) :: C
      real(8) :: AMAX,STO,W1,W
      integer :: I, I4, MM, J, ISTO, K, J1, I1, J2, I2
      !***** THIS SUBROUTINE EVALUATES THE INVERSE OF A
      !***** SQUARE M*M MATRIX A
      !***** THE DIMENSIONS OF THE ARRAYS MAY BE INTEGER VARIABLES WHEN USED*********
      !***** ON THE IBM 7094,BUT THEY MUST BE INTEGER CONSTANTS ON THE IBM 1130******
      
      AMAX=0.0

      do I=1,M
         IND(I)=I
         if(abs(A(I,1))-AMAX > 0) then
            AMAX=abs(A(I,1))
            I4=I
         endif
      end do

      MM=M-1

      do J=1,MM
         if(I4-J > 0) then
            ISTO=IND(J)
            IND(J)=IND(I4)
            IND(I4)=ISTO
            do K=1,M
               STO=A(I4,K)
               A(I4,K)=A(J,K)
               A(J,K)=STO
            end do
         endif
         AMAX=0.0
         J1=J+1
         do I=J1,M
            A(I,J)=A(I,J)/A(J,J)
            do K=J1,M
               A(I,K)=A(I,K)-A(I,J)*A(J,K)
               IF (K-J1 <= 0) then
                  IF(abs(A(I,K))-AMAX > 0) then
                     AMAX=abs(A(I,K))
                     I4=I
                  endif
               endif
            end do
         end do
      end do

      do I1=1,MM
         I=M+1-I1
         I2=I-1
         do J1=1,I2
            J=I2+1-J1
            J2=J+1
            W1=-A(I,J)
            IF(I2-J2 >= 0) then
               do K=J2,I2
                  W1=W1-A(K,J)*C(K)
               end do
            endif
            C(J)=W1
         end do
         do K=1,I2
            A(I,K)=C(K)
         end do
      end do

      do I1=1,M
         I=M+1-I1
         I2=I+1
         W=A(I,I)
         do J=1,M
            if (I-J < 0) then
               W1=0.0
            elseif(I-J == 0) then
               W1=1.0
            else
               W1=A(I,J)
            endif
            if (I1-1 > 0) then
               do K=I2,M
                  W1=W1-A(I,K)*A(K,J)
               end do
            endif
            C(J)=W1
         end do
         do J=1,M
            A(I,J)=C(J)/W
         end do
      end do

      do I=1,M
         do while (IND(I)-I /=0)
            J=IND(I)
            do K=1,M
               STO=A(K,I)
               A(K,I)=A(K,J)
               A(K,J)=STO
            end do
            ISTO=IND(J)
            IND(J)=J
            IND(I)=ISTO
         end do
      end do
      return
   END SUBROUTINE
   ! **************************************************************************

    !===============================================================================
    SUBROUTINE argSort(n, a, order)
      integer :: n, i, j
      integer, intent(in) :: a(n)
      integer, intent(out) :: order(n)
      real :: scr(n)
      !-----------------------------------------------------------------------------
      scr = a
      order = 1
   
      do i = 2, n
         j = i - 1
         do while (j >= 1)
            if (scr(j) <= a(i)) exit
            scr(j + 1) = scr(j)
            order(j+1) = order(j)
            j = j - 1
         end do
         scr(j + 1) = a(i)
         order(j+1) = i
      end do
    END SUBROUTINE argSort

    !===============================================================================
    elemental FUNCTION air2vac (air) result(vac)
      real(8), intent(in) :: air
      real(8) :: vac, n, s, s2
      !-----------------------------------------------------------------------------
      s = 1d4 / air
      s2 = s*s
      n = 1d0 + 0.00008336624212083d0 + 0.02408926869968d0 / (130.1065924522d0 - s2) &
            + 0.0001599740894897d0 / (38.92568793293d0 - s2)
      vac = air*n
    END FUNCTION air2vac
END MODULE op_math_mod