!
!  Test program for radiative transfer.  The solutions Q_p for source
!  functions S_p = S(tau)=tau**p are
!
!  Q_1 = 0.5*exp(-tau)
!  Q_2 = 2 - exp(-tau)
!  Q_3 = 6*tau + 3*exp(-tau)
!
!  However, for a logarithmicallyu spaced tau-scale this assumes that
!  one inputs the correct forward intensity at tau=tau(1).  These are:
!
!  I_1 = -1 + tau + exp(-tau)
!  I_2 = 2. - 2.*tau + tau**2 - 2.*exp(-tau)
!  I_3 = -6. + 6.*tau - 3.*tau**2 + tau**3 + 6.*exp(-tau)
!
 implicit none
 real, allocatable, dimension(:):: s, tau, q, q0, dtau
 integer :: i, n, m, niter, iter, nmod
 real :: logtaumin, logtaumax, dlogtau, logtau, cput, cpu(2), dtime
 real :: Ifw1, IfwN, Irv1, IrvN, P, error, sigma
n=91
 logtaumin = -6.
 logtaumax = +3.
 dlogtau = (logtaumax-logtaumin)/float(n-1)!0.1
 sigma = 1.4

 !n = (logtaumax-logtaumin)/dlogtau + 1.5
 allocate (s(n), tau(n), q(n), q0(n), dtau(n))

 do m=1,4
 do i=1,n
   logtau  = logtaumin+dlogtau*(i-1)
   tau(i) = 10.**logtau
   s(i)   = tau(i)**m
   select case (m)
   case (1)
     q0(i)  = 0.5*exp(-tau(i))
   case (2) 
     q0(i)  = 2.-exp(-tau(i))
   case (3)
     q0(i)  = 6.*tau(i)+3.*exp(-tau(i))
   case (4)
     s(i)   = (erf(logtau/sigma)+1.)*0.5
     q0(i)  = 0.
   end select
   if (i > 1) dtau(i) = tau(i)-tau(i-1)
 end do

 select case (m)
 case (1)
   Ifw1 = -1. + tau(1) + exp(-tau(1))
   IrvN = +1. + tau(n)
 case (2) 
   Ifw1 = 2. - 2.*tau(1) + tau(1)**2 - 2.*exp(-tau(1))
   IrvN = 2. + 2.*tau(n) + tau(n)**2
 case (3)
   Ifw1 = -6. + 6.*tau(1) - 3.*tau(1)**2 + tau(1)**3 + 6.*exp(-tau(1))
   IrvN = +6. + 6.*tau(n) + 3.*tau(n)**2 + tau(n)**3
 case (4)
   Ifw1 = 0.
   IrvN = s(n)
    tau(1) = 1e-10
   s(1)   = (erf(alog10(tau(1))/sigma)+1.)*0.5
   dtau(2) = tau(2)-tau(1)
 end select
 !IrvN = S(n)

 niter=50000/n*50
 cput=dtime(cpu)
 do iter=1,niter
   call transfer (n,dtau,s,q,ifwN,.false.)
 end do
 cput=dtime(cpu)
 print *,'cpu =',cput
 niter = niter*2./cput                                                  ! number of iterations for 2 seconds CPU

 cput=dtime(cpu)
 do iter=1,niter 
   call transfer (n,dtau,s,q,ifwN,.false.)
 end do
 cput=dtime(cpu)
 print *, 'Ifw1=',Ifw1
 print *,'  i  log(tau)         S           P           Q       Q0        err'
 nmod=0.1/dlogtau*5.+0.1
 print*,'nmod',nmod
 do i=1,n
   P = s(i)+q(i)
   error = (q(i)-q0(i))/P
  ! if (mod(i-1,nmod) == 0) then
     print '(i4,5e12.4,f12.6)',i, alog10(tau(i)), s(i), P, q(i), q0(i), error
  ! endif
 end do
 print *,'cpu =',cput
 print '(1x,a,f6.1)','ns/pt =',cput*1e9/(n*niter)
 write(11) m,n,tau,s,q,q0

 end do

 END
