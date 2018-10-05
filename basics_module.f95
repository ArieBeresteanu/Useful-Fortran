MODULE basics

! This module contains important definitions and constants to be
! used in other modules or programs.

IMPLICIT NONE
INTEGER, PARAMETER 		:: dp = kind(1.0d0)
INTEGER, PARAMETER 		:: regInt = selected_int_kind (8)
INTEGER, PARAMETER 		:: longInt = selected_int_kind (16)

REAL(dp), PARAMETER		:: ZERO=0.00_dp
REAL(dp), PARAMETER		:: ONE=1.00_dp
REAL(dp), PARAMETER		:: negONE=-1.0_dp
REAL(dp), PARAMETER		:: TWO=2.0_dp
REAL(dp), PARAMETER     :: THREE=3.0_dp
REAL(dp), PARAMETER     :: HALF=0.5_dp

REAL(dp), PARAMETER     :: PI=3.1415926535897932_dp
REAL(dp), PARAMETER	    :: Two_PI=6.2831853071795864_dp
REAL(dp), PARAMETER     :: Half_PI=1.5707963267948966_dp
REAL(dp), PARAMETER     :: sqrtOfTwoPI=2.50662827463100051_dp
REAL(dp), PARAMETER     :: invSqrtOfTwoPi=0.398942280401432678_dp

END MODULE basics
