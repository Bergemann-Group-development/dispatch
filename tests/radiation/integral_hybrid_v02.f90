! $Id$
!***********************************************************************
SUBROUTINE transfer (n, dtau, S, I_in, I_in2, q, I_out, I_out2)
  implicit none
  real, dimension(n) :: dtau, S, q
  real  :: I_in, I_out, I_in2, I_out2
  real, dimension(n) :: qfw, qre  
  
  !-----------------------------------------------------------------------
  real(kind=8), dimension(n) :: ex1,ex2, ds1, d1s, d2s, ex0
  real(kind=8) :: ctau
  integer i, n
  
  do i=2,n
    ex0(i) = exp(-dtau(i))
    ex1(i) = 1.-ex0(i)
    ex2(i) = ex1(i) - dtau(i)*ex0(i)
    ds1(i) = (s(i) - s(i-1)) / dtau(i)                  ! half-point centered dS/dtau
  end do 
  
  do i=2,n-1
    ctau = 1. / (dtau(i+1) + dtau(i))
    d1s(i) =(ds1(i+1) * dtau(i) &
           + ds1(i  ) * dtau(i+1)) * ctau                   ! dS/dtau, 2nd. order accurate
    d2s(i)= (ds1(i+1) - ds1 (i  )) * 2. * ctau              ! d2S/dtau2, 2nd. order accurate
  end do
  d2s(1) = d2s(2)                                               ! assume constant d2S/dtau2
  d1s(1) = d1s(2) - dtau(2)*d2s(2)                      ! extrapolate to 1
  d2s(n) = d2s(n-1)
  d1s(n) = d1s(n-1) - dtau(n)*d2s(n)

  qfw(1) = I_in - s(1)                                            ! initial forward Q
  do i=2,n
    qfw(i) = qfw(i-1) *ex0(i) &
             - d1s(i) *ex1(i) &
             + d2s(i) *ex2(i)                                   ! formal solution forward
  end do

    I_out = qfw(n) + s(n)                                         ! forward boundary intensity
    
    
  qre(1) = I_in2 - s(n)                                            ! initial reverse Q
  
  do i=n-1,1,-1
    qre(i) = qre(i+1) *ex0(i+1) &
             + d1s(i) *ex1(i+1) &
             + d2s(i) *ex2(i+1)                                  ! formal solution reverse
  end do
    

    I_out2 = q(1) + s(1)                                         ! reverse boundary intensity
q = 0.5*(qre+qfw)

END
