module prec_mod
    use, intrinsic :: iso_fortran_env
    use mconst
    implicit none

    !integer, parameter :: sp = REAL32 ! in conflict with spectrum object
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128
     
    ! physics
    real(dp), parameter :: Ryd = 13.598408    ! Rydberg constant [eV]
    real(dp), parameter :: EK  = EE/BK
    real(dp), parameter :: XMH = 1.67339e-24
    !real(dp), parameter :: op_saha_const = (2*PI*EM*BK/HH**2)**1.5
    real(dp), parameter :: op_saha_const = 1.202e9/BK/5040**2.5/2
    !real(dp), parameter :: temp_to_theta = log10(exp(1.d0))*EK
    real(dp), parameter :: temp_to_theta = 5040
    real(dp), parameter :: prec_err = temp_to_theta / (log10(exp(1.d0))*EK)
    real(dp), parameter :: uop_BK = 1.38053e-16
    !real(dp), parameter :: uop_BK = BK
    
contains
    
end module prec_mod