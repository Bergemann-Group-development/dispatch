! cooling_Casiana.f90 $Id$
!***********************************************************************
FUNCTION cooling_Casiana (self, temp)
  implicit none
  class(cooling_t):: self
  real(kind=dp) cooling_Casiana, temp, qq  
  real(kind=dp), parameter:: &
    TEm=1e1, &
    TE0=1e3, &
    TE1=10964.781e0, &
    TE2=19952.632e0, & 
    TE3=29512.078e0, & 
    TE4=39810.708e0, & 
    TE5=69183.121e0, &
    TE6=125892.51e0, &
    TE7=251188.7e0, &
    TE8=501187.01e0, &
    TE9=891250.55e0, &
    TE10=1.412537e6, & 
    TE11=2.511887e6, &
    TE12=6.309576e6, &
    TE13=1.e7, &
    TE14=1.5848925e7, &
    TE15=1.e9, &

    Cm=2.,   &
    C0=2.596,   &
    C1=3.846, &
    C2=-0.1764, &
    C3=1.00001, &
    C4=1.6666, &
    C5=1.15384, &
    C6=0.1, &
    C7= -3.933, &
    C8= 0.72, &
    C9= -0.65,  &
    C10= 0.4, &
    C11= -1.25, &
    C12= 0., &
    C13= -0.4,  &
    C14= 0.43333, &
    
    Hm=4.98282e-05*4e-28, &
    H0=3.2473989e-34, &
    H2=5.739e-22, &
    H3=3.1620277e-27, &
    H4=2.7123742e-30, &
    H5=8.2298862e-28, &
    H6=1.9497e-22, &
    H7=1.1749e+0*4e-28, &
    H8=3.5155e-27, &
    H9=4.982757e-19, &
    H10=1.7379e-25, &
    H11=6.309e-15, &
    H12=1.99525e-22, &
    H13=1.2589e-20, &
    H14=1.2589e-26

  real(kind=8), parameter:: &
    H1=2.894d-39

    if      (temp .le. T_MC) then
      qq = 0.
    else if (temp .le. TE0) then
      qq = Hm*temp**Cm
    else if (temp .le. TE1) then
      qq = H0*temp**C0
    else if (temp .le. TE2) then
      qq = H1*temp**C1
    else if (temp .le. TE3) then
      qq = H2*temp**C2
    else if (temp .le. TE4) then
      qq = H3*temp**C3
    else if (temp .le. TE5) then
      qq = H4*temp**C4
    else if (temp .le. TE6) then
      qq = H5*temp**C5
    else if (temp .le. TE7) then
      qq = H6*temp**C6
    else if (temp .le. TE8) then
      qq = (H7*temp**C7)*2.5e27
    else if (temp .le. TE9) then
      qq = H8*temp**C8
    else if (temp .le. TE10) then
      qq = H9*temp**C9
    else if (temp .le. TE11) then
      qq = H10*temp**C10
    else if (temp .le. TE12) then
      qq = H11*temp**C11
    else if (temp .le. TE13) then
      qq = H12*temp**C12
    else if (temp .le. TE14) then
      qq = H13*temp**C13
    else
      qq = H14*temp**C14
    end if
    cooling_Casiana=qq
END FUNCTION
