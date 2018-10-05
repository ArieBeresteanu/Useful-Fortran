MODULE regression

USE basics
USE utilities, ONLY: LP2

IMPLICIT NONE

CONTAINS

SUBROUTINE studentize(rows,columns,X_in,X_out,mu,sigma)

	!
	!   ===================================================================
	!   Studentizing the variables in the matrix X_in buy substructing the
	!   mean and then dividing by the standard deviation of each column.
	!   ===================================================================
	!

	INTEGER(regInt), INTENT(IN)					:: rows,columns
	REAL(dp), DIMENSION(rows,columns),INTENT(IN):: X_in
	REAL(dp), DIMENSION(rows,columns),INTENT(OUT):: X_out
	REAL(dp), DIMENSION(columns), INTENT(OUT)	:: mu,sigma

	INTEGER(regInt)			:: i
	REAL(dp)				:: rN

	rN=REAL(rows,kind=8)

	mu=SUM(X_in)/rN

	DO i=1,columns
		sigma(i)=SQRT(SUM((x_in(:,i)-mu(i))**2)/(rN-ONE))
		X_out(:,i)=(X_in(:,i)-mu(i))/sigma(i)
	END DO

END SUBROUTINE studentize
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

SUBROUTINE makeGrid(X,mesh,padding,Xgrid)
	!
	! An equidistance grid is built from the lower to the upper
	! value of X. The user can ask to add some "padding" which
	! means that the grid will not start and end at the minimum
	! and maximum but will leave 1% padding from each side.
	!
		REAL(dp), INTENT(IN)			:: X(:)
		INTEGER(regInt), INTENT(IN)		:: mesh
		LOGICAL, INTENT(IN), OPTIONAL	:: padding
		REAL(dp), INTENT(OUT)			:: Xgrid(mesh)

		INTEGER(regInt)	:: i
		REAL(dp)	:: low,high,length,skip,pad
		REAL(dp), PARAMETER	:: tol=1.0D-06
		REAL(dp), PARAMETER	:: padding_percent=0.005D0

		high=MAXVAL(X)
		low =MINVAL(X)
		length=high-low
		IF (length.LT.tol .OR. mesh.EQ.1) THEN
			Xgrid=high
		ELSE
			pad=ZERO
			IF (PRESENT(padding)) pad=padding_percent*length

			skip=(high-low-TWO*pad)/(mesh-1)
			DO i=1,mesh
				Xgrid(i)=low+pad+skip*(i-1)
			END DO
		END IF
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
END SUBROUTINE makeGrid

SUBROUTINE makeGrid2(low,high,mesh,Xgrid)
    REAL(dp), INTENT(IN)        :: low,high
    INTEGER(regInt), INTENT(IN) :: mesh
    REAL(dp), INTENT(OUT)       :: Xgrid(mesh)

    INTEGER(regINT) :: i

    DO i=1,mesh
        Xgrid(i)=low+(i-1)*(high-low)/(mesh-1)
    END DO

END SUBROUTINE makeGrid2
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


REAL(dp) FUNCTION epanechnikov(x0,bandwidth)

!	============================
! 	Epanechnikov kernel function
!	============================

REAL(dp), INTENT(IN)	:: x0, bandwidth

   IF (ABS(x0).LE.bandwidth) THEN
       epanechnikov=0.750D0*(ONE-(x0/bandwidth)**2)
   ELSE
       epanechnikov=ZERO
   END IF
END FUNCTION epanechnikov
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


REAL(dp) FUNCTION silverman(N,vec)
	!
	! Computing the Silverman rule of thumb bandwidth
	!
	INTEGER(regInt), INTENT(IN)		:: N
	REAL(dp), INTENT(IN)	        :: vec(N)

	REAL(dp)		:: mu,sigma,rN

	rN=N+0.0_dp
	mu=SUM(vec)/rN
	sigma=SQRT(SUM((vec-mu)**2)/(rN-ONE))

	silverman = 1.060_dp*sigma*(N**(-0.20_dp))
END FUNCTION silverman
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


END MODULE regression
