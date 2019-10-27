! cooling_Rosen_etal.f90 $Id$
!***********************************************************************
FUNCTION cooling_Rosen_etal (self, temp)
  implicit none
  class(cooling_t):: self
  real(kind=dp) cooling_Rosen_etal, temp, qq
  real(kind=dp), parameter:: &
    TE1=10.0, &                ! extended down to 10 K from 300 K.
    TE2=2000.0, &
    TE3=8000.0, &
    TE4=1.0e5, &
    TE5=4.0e7, &
    C1=2.0, &
    C2=1.5, &
    C3=2.867, &
    C4=-0.65, &
    C5=0.5, &
! The coefficients below are the ones from Rosen, Bregman & Normman 1993,
! with a factor of 1e20 multiplied on the coefficients, which is cancelled
! by dividing not with 2e-24**2, but with 4e-28.  The factor 2e-24 is the
! approximate molecular weight of ISM particles (without electrons), and
! converts between hydrogen number density and mass density.  With the 
! original coefficients the cooling is obtained by multiplying the cooling 
! function with CGS number density squared, while with these much larger
! coefficients one multiplies instead with CHS mass density squared. The
! expression in the paper has a multiplication with "rho^2", but this should
! be number density squared (N^2). The expression in Rosen & Bregman 1995,
! which has a factor pressure squared, is also incorrect.  Rosen et al state
! that the expressions take account of the varying degree of ionization (the
! cooling is really proportional to N_H * N_e).
    H1=2.2380e-32, &
    H2=1.0012e-30, &
    H3=4.6240e-16, &
    H4=1.7800e-18, &
    H5=3.2217e-27
!.......................................................................
    if      (temp .le. T_MC) then
      qq = 0.
    else if (temp .le. TE2) then
      qq = H1*temp**C1
    else if (temp .le. TE3) then
      qq = H2*temp**C2
    else if (temp .le. TE4) then
      qq = H3*temp**C3*1e-20
    else if (temp .le. TE5) then
      qq = H4*temp**C4
    else
      qq = H5*temp**C5
    end if
    cooling_Rosen_etal=qq
END FUNCTION
