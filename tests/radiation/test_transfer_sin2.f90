!
!  Test program for radiative transfer.  The solutions Q_p for source
!  functions S(tau)=tau**p are
!
!  Q_1 = 0.5*exp(-tau)
!  Q_2 = 2 - exp(-tau)
!  Q_3 = 6*tau - 3*exp(-tau)
!
!  However, for a logarithmicallyu spaced tau-scale this assumes that
!  one inputs the correct forward intensity at tau=tau(1).  These are:
!
!  I_1 = -1 + tau + exp(-tau)
!  I_2 = 2. - 2.*tau + tau**2 - 2.*exp(-tau)
!  I_3 = -6. + 6.*tau - 3.*tau**2 + tau**3 + 6.*exp(-tau)
!
 implicit none
 real, allocatable, dimension(:,:) :: s, q, dtau
 real, allocatable, dimension(:)   :: q0, tau, Ifw1, IfwN, Irv1, IrvN
 integer :: i, j, k, n, m, niter, iter, nn(6), nmod
 real :: taumin, taumax, cput, cpu(2), dtime, dlogtau
 real :: JfwN, Jrv1, P, error, pi, maxerr
 real(kind=8) dtau1

 pi = 4.*atan(1.)
 taumin = -pi
 taumax =  pi
 
 k = 1
 m = 50

 nn = (/10,20,50,100,200,500/)
 do j=1,6
 n = nn(j)
 print*,''
 print*,'N =',n

 allocate (s(m,n), tau(n), q(m,n), q0(n), dtau(m,n))
 allocate (Ifw1(m), IfwN(m), Irv1(m), IrvN(m))

 dtau1 = (taumax-taumin)/(n-1)
 dtau  = dtau1
 do i=1,n
   tau(i)   = taumin+(i-1)*dtau1
   s(:,i)   = 1.+sin(k*tau(i))
   q0(i)  = -k**2*sin(k*tau(i))/(1.+k**2)
 end do

 Ifw1 = 1.+(-k*cos(k*tau(1))+sin(k*tau(1)))/(1.+k**2)
 JfwN = 1.+(-k*cos(k*tau(n))+sin(k*tau(n)))/(1.+k**2)

 Jrv1 = 1.+(+k*cos(k*tau(1))+sin(k*tau(1)))/(1.+k**2)
 IrvN = 1.+(+k*cos(k*tau(n))+sin(k*tau(n)))/(1.+k**2)

 niter=50000/n*50/m
 cput=dtime(cpu)
 do iter=1,niter
   call transfer2 (m,n,dtau,s,Ifw1,IrvN,q,Irv1,IfwN)
 end do
 cput=dtime(cpu)
 print *,'cpu =',cput
 niter = niter*2./cput                                                  ! number of iterations for 2 seconds CPU

 cput=dtime(cpu)
 do iter=1,niter
   call transfer2 (m,n,dtau,s,Ifw1,IrvN,q,Irv1,IfwN)
 end do
 cput=dtime(cpu)

! print *,'  i      tau         S           P           Q          err'
 maxerr = 0.
nmod=0.1/n*5.+0.1
 print*,'nmod',nmod
 do i=1,n
   P = s(i)+q(i)
   error = (q(i)-q0(i))/P
   maxerr = max(maxerr,error)
   if (mod(i-1,nmod) == 0) then
     print '(i4,4e12.4,f12.6)',i, alog10(tau(i)), s(i), P, q(i), error
   endif
 end do
 print *,(1.-Irv1(1)/Jrv1),(1.-IfwN(1)/JfwN),maxerr
 print *,'cpu =',cput
 print '(1x,a,f6.1)','ns/pt =',cput*1e9/(n*m*niter)

 deallocate (s, tau, q, q0, dtau)
 deallocate (Ifw1, IfwN, Irv1, IrvN)
 end do

 END
