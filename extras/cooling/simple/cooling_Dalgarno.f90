! cooling_Dalgarno.f90 $Id$
!***********************************************************************
FUNCTION cooling_Dalgarno (self, temp)
!
! Fit to cooling profile, merged between Dalgarno & McCray(1972)
! for low temp, and Raymond, Cox & Smith(1976) for high temp.
! C4/H4 corrected relative to sncool.pro to make it smooth at T=1e5.

  implicit none
  class(cooling_t):: self
  real(kind=dp) cooling_Dalgarno, temp, qq  
  real(kind=dp), parameter:: &
    TE1=10.0, &
    TE2=2000.0, &
    TE3=8000.0, &
    TE4=2.5e4, &
    TE5=1.0e5, &
    TE6=2.5e5, &
    TE7=4.0e5, &
    TE8=4.5e7, &
    C1=2.0, &
    C2=1.5, &
    C3=5.0, &
    C4=0.83, &
    C5=0.0, &
    C6=-2.54, &
    C7=-0.49, &
    C8=0.5, &
    H1=8.00e15, &
    H2=3.58e17, &
    H3=7.79e03, &
    H4=1.70e22, &
    H5=2.38e26, &
    H6=1.22e20, &
    H7=4.00e28, &
    H8=1.07e21

    if      (temp .le. T_MC) then
      qq=0.
    else if (temp .le. TE2) then
      qq=H1*temp**C1
    else if (temp .le. TE3) then
      qq=H2*temp**C2
    else if (temp .le. TE4) then
      qq=H3*temp**C3
    else if (temp .le. TE5) then
      qq=H4*temp**C4
    else if (temp .le. TE6) then
      qq=H5*temp**C5
    else if (temp .le. TE7) then
      qq=H6*temp**C6*1e20
    else if (temp .le. TE8) then
      qq=H7*temp**C7
    else
      qq=H8*temp**C8
    end if
    cooling_Dalgarno=qq*4e-28
END FUNCTION
