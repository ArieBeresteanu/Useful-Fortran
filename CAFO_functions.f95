MODULE CAFO_functions

IMPLICIT NONE

!!! global parameters
REAL(8), PARAMETER  :: Earth_Radius=6371000.0D0 
REAL(8), PARAMETER  :: Pi=3.14159265359D0 
REAL(8), PARAMETER	:: Two_PI=6.28318530718D0
REAL(8), PARAMETER	:: HALF=0.50D0, ONE=1.00D0, ZERO=0.00D0

REAL(8), PARAMETER 	:: bandwidth=3000.0D0

CONTAINS

REAL(8) FUNCTION sphere_dist(lon1,lat1,lon2,lat2)

!	=======================================
!	distance between two points on a sphere
!	=======================================

REAL(8), INTENT(IN) :: lon1,lat1,lon2,lat2
REAL(8)				:: phi1,phi2,delta_phi,delta_lambda,a,c

phi1=lat1/(Two_PI) 
phi2=lat2/(Two_PI) 
delta_phi=(lat1-lat2)/(Two_PI) 
delta_lambda=(lon1-lon2)/(Two_PI) 
a=SIN(HALF*delta_phi)**2 + COS(phi1)*COS(phi2)*(SIN(HALF*delta_lambda))**2 
c=ATAN2(SQRT(a),SQRT(1.0-a)) 
sphere_dist=Earth_Radius*c

END FUNCTION sphere_dist



REAL(8) FUNCTION epanechnikov(x0,bandwidth)

!	============================
! 	Epanechnikov kernel function
!	============================

REAL(8), INTENT(IN)	:: x0, bandwidth

    IF (ABS(x0).LE.bandwidth) THEN
        epanechnikov=0.750D0*(ONE-(x0/bandwidth)**2)
    ELSE
        epanechnikov=ZERO
    END IF
END FUNCTION epanechnikov



FUNCTION to_lower (str) RESULT (string)

!   ==============================
!   Changes a string to lower case
!   ==============================

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: str
    CHARACTER(LEN(str))      :: string

    INTEGER :: ic, i

    CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'

!   lowercase each letter if it is capitalized
    string = str
    DO i = 1, LEN_TRIM(str)
        ic = INDEX(cap, str(i:i))
        IF (ic .GT. 0) string(i:i) = low(ic:ic)
    END DO

END FUNCTION to_lower



! REAL(8) FUNCTION haversine(lat1,lon1,lat2,lon2)

!	================================
!	Haversine distance over a sphere
!	================================

! REAL(8), INTENT(IN) 	:: lat1,lon1,lat2,lon2

! haversine=2.0D0*Earth_Radius*ASIN(SQRT(SIN(HALF*(lat2-lat1))**2 + COS(lat1)*COS(lat2)*SIN(HALF*(lon2- lon1))**2))

! END FUNCTION haversine

FUNCTION to_radian(degree) RESULT(rad)
! 	===================
! 	degrees to radians
! 	===================

    REAL(8),INTENT(IN) :: degree
    REAL(8), PARAMETER :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
    REAL(8) :: rad
 
	rad = degree*deg_to_rad
END FUNCTION to_radian
 
FUNCTION haversine(deglat1,deglon1,deglat2,deglon2) RESULT (dist)
!	=============================================
! 	great circle distance -- adapted from Matlab 
!	source: https://rosettacode.org/wiki/Haversine_formula#Fortran
!	=============================================

    REAL(8),INTENT(IN) 	:: deglat1,deglon1,deglat2,deglon2
    REAL(8) 			:: a,c,dist,dlat,dlon,lat1,lat2
    REAL(8),PARAMETER 	:: radius = 6372800.0 !in meters
 
    dlat = to_radian(deglat2-deglat1)
    dlon = to_radian(deglon2-deglon1)
    lat1 = to_radian(deglat1)
    lat2 = to_radian(deglat2)
    a = (SIN(dlat/2))**2 + COS(lat1)*COS(lat2)*(SIN(dlon/2))**2
    c = 2*ASIN(SQRT(a))
    dist = radius*c
END FUNCTION haversine
 
FUNCTION compassDirection(deglat1,deglon1,deglat2,deglon2) RESULT(deg)
!	=============================================
! 	receives the latitude and longitude of two  
!	points on heart and returns the degree a 
!	compass will show if you look from point 1 to
!	point 2.
!	=============================================
	REAL(8), INTENT(IN)	:: deglat1,deglat2,deglon1,deglon2
	REAL(8)				:: deg
	REAL(8)				: a,c
    radlat1=to_radian(deglat1)
    radlat2=to_radian(deglat2)
    radlon1=to_radian(deglon1)
    radlon2=to_radian(deglon2)
    a=acos(SIN(radlat1)*SIN(radlat2)+COS(radlat1)*COS(radlat2)*COS(radlon2-radlon1))
    c=asin((SIN(radlon2-radlon1)*SIN(pi/2-radlat2))/SIN(a))
    IF (deglat1>deglat2) THEN
        deg=180.0D0-(180.0D0*c/pi)
    ELSEIF (c<0)
            deg=360.0D0+(180.0D0*c/pi)
    ELSE
            deg= 180.0D0*c/pi
    END IF
END FUNCTION compassDirection


REAL(8) FUNCTION heat(N,y,x1,x2,lat,lon,bandwidth)

!	==========================================
!	heat at point (lat,lon) from y,x1,x2
!	computing the influence of the y's
!	located in (x1,x2)'s on location (lat,lon)
!	==========================================

	INTEGER, INTENT(IN)					:: N
	REAL(8), INTENT(IN), DIMENSION(N)	:: y,x1,x2
	REAL(8), INTENT(IN)					:: lon,lat,bandwidth

	REAL(8)								:: numerator 
	REAL(8)								:: dist,k
    INTEGER								:: i
	!denominator=ZERO
    numerator=ZERO
    DO i=1,N
        dist=haversine(lat,lon,x1(i),x2(i))  !using haversine distance
        k=epanechnikov(dist,bandwidth) 
        !denominator = denominator+k 
        numerator = numerator+k*y(i) 
	END DO
	heat=numerator
END FUNCTION heat

FUNCTION kernreg1(N,y,x,bandwidth)

!	======================================
!	kernel regression E(y|x) evaluated for
!	each x using a global bandwidth.
!	======================================

	INTEGER, INTENT(IN)					:: N
	REAL(8), INTENT(IN), DIMENSION(N)	:: y,x
	REAL(8), INTENT(IN)					:: bandwidth

	REAL(8), DIMENSION(N)	:: kernreg1

	! local variables
	REAL(8), PARAMETER	:: tol=1.0D-06
	REAL(8)				:: denominator, numerator, temp
	INTEGER				:: i,j
	
	DO i=1,N
		denominator=ZERO
		numerator=ZERO
		DO j=1,N 
			temp=epanechnikov(x(i)-x(j),bandwidth)
			denominator=denominator+temp
			numerator=numerator+temp*y(j)
		END DO
		IF (denominator.GE.tol) THEN
			kernreg1(i)=numerator/denominator
		ELSE
			kernreg1(i)=ZERO
		END IF
	END DO
END FUNCTION kernreg1


SUBROUTINE studentize(rows,columns,X_in,X_out,mu,sigma)

!
!   ===================================================================
!   Studentizing the variables in the matrix X_in buy substructing the 
!   mean and then dividing by the standard deviation of each column.
!   ===================================================================
!

INTEGER, INTENT(IN)							:: rows,columns
REAL(8), DIMENSION(rows,columns),INTENT(IN)	:: X_in
REAL(8), DIMENSION(rows,columns),INTENT(OUT):: X_out
REAL(8), DIMENSION(columns), INTENT(OUT)	:: mu,sigma

INTEGER				:: i
REAL(8)				:: rN

rN=REAL(rows,kind=8)

mu=SUM(X_in)/rN
sigma=SQRT(SUM(X_in**2)/rN - mu**2)

DO i=1,columns
	X_out(:,i)=(X_in(:,i)-mu(i))/sigma(i)
END DO

END SUBROUTINE studentize


				



!	End of Module

END MODULE CAFO_functions