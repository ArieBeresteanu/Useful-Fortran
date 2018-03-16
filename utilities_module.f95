MODULE utilities
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
IMPLICIT NONE
INTEGER, PARAMETER 		:: dp = kind(1.0d0)
REAL(dp), PARAMETER		:: ZERO=0.00_dp
REAL(dp), PARAMETER		:: ONE=1.00_dp

CONTAINS
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION LP2(vec) RESULT(dist)

REAL(dp), INTENT(IN)		:: vec(:)
REAL(dp)					:: dist

dist=SQRT(SUM(vec**2))

END FUNCTION LP2

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION LP1(vec) RESULT(dist)
REAL(dp), INTENT(IN)		:: vec(:)
REAL(dp)					:: dist

dist=SUM(ABS(vec))

END FUNCTION LP1

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION LPinf(vec) Result(dist)
REAL(dp), INTENT(IN)		:: vec(:)
REAL(dp)					:: dist

dist=MAXVAL(ABS(vec))

END FUNCTION LPinf

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
SUBROUTINE descriptive(Nobs,x,mu,sigma)
INTEGER, INTENT(IN)						:: Nobs
REAL(dp), INTENT(IN), DIMENSION(Nobs)	:: X
REAL(dp), INTENT(OUT)					:: mu,sigma

REAL(dp)		:: rNobs


rNobs=Nobs+0.0D0
mu=SUM(X)/rNobs
sigma=SQRT(SUM((X-mu)**2)/(rNobs-ONE))
END SUBROUTINE descriptive

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
SUBROUTINE printRealMatrix(A,f)
REAL(dp), INTENT(IN), DIMENSION(:,:)	:: A
CHARACTER(LEN=*), INTENT(IN), OPTIONAL	:: f

INTEGER	:: i,m

m=SIZE(A,1) !number of rows
DO i=1,m
	IF (PRESENT(f)) THEN 
		PRINT f, A(i,:)
	ELSE
		PRINT *,A(i,:)
	END IF
END DO

END SUBROUTINE printRealMatrix

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
SUBROUTINE printRealVector(vec,f)
REAL(dp), INTENT(IN), DIMENSION(:)		:: vec
CHARACTER(LEN=*), INTENT(IN), OPTIONAL	:: f

INTEGER	:: i,m

m=SIZE(vec) !number of elements
DO i=1,m
	IF (PRESENT(f)) THEN 
		PRINT f, vec(i)
	ELSE
		PRINT *,vec(i)
	END IF
END DO

END SUBROUTINE printRealVector




END MODULE utilities