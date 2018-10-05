MODULE kernelRegression

USE regression


IMPLICIT NONE


INTERFACE kernreg
    MODULE PROCEDURE kernreg_kD_global
    MODULE PROCEDURE kernreg_1D_global
END INTERFACE

INTERFACE kernreg_sorted
    MODULE PROCEDURE kernreg_kD_global_sorted
    MODULE PROCEDURE kernreg_1D_global_sorted
END INTERFACE

CONTAINS


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! Kernel Regression !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE kernreg_kD_global(y,X,Xgrid,bandwidth,y_fit)

	!
	! Estimating a Kernel regression with:
	!   a. k dimensional covariatesa
	!   b. global bandwidth
	! E(y|X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to y_fit
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!   X and Xgrid (if exists) are not necessarily sorted

	REAL(dp), INTENT(IN),  DIMENSION(:)					:: y
	REAL(dp), INTENT(IN),  DIMENSION(:,:)				:: X
	REAL(dp), INTENT(IN),  DIMENSION(:,:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN),  DIMENSION(:),   OPTIONAL		:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:),   ALLOCATABLE	:: y_fit

	! additional variables
	INTEGER(regInt)				        :: Nobs !number of obsevations
	INTEGER(regInt)				        :: Ncov !Number of covariates (Note: do not include a constant)
	REAL(dp), DIMENSION(:), ALLOCATABLE	:: bwidth !the bandwidth vector we will use
	REAL(dp), DIMENSION(:), ALLOCATABLE	:: point
	INTEGER(regInt)				        :: i,j,N
	REAL(dp)					        :: dist,denom,numer,epa
	REAL(dp), PARAMETER			        :: tol=1.0D-06
	LOGICAL                             :: xgridFLAG

	Nobs=SIZE(y) !number of observations
	IF (Nobs .NE. SIZE(X,1)) THEN !dimension mismatch
		y_fit=ZERO
		RETURN
	END IF

	Ncov=SIZE(X,2) !number of covariates
	ALLOCATE(bwidth(Ncov),point(Ncov))

	IF (PRESENT(bandwidth) .AND. (SIZE(bandwidth) .EQ. Ncov)) THEN
		bwidth=bandwidth
	ELSE !use silverman rule of thumb
		DO i=1,Ncov
			bwidth(i)=silverman(Nobs,X(:,i))
		END DO
	END IF

	!if Xgrid is an input
	IF (PRESENT(Xgrid)) THEN
		N=SIZE(Xgrid,1)
		xgridFLAG=.TRUE.
	ELSE
		N=Nobs
		xgridFLAG=.FALSE.
	END IF

	ALLOCATE(y_fit(N))
	!if Xgrid is present and its 2nd dimention matches that of the covariates,
	DO i=1,N
		denom=ZERO
		numer=ZERO
		IF (xgridFLAG) THEN
			point=Xgrid(i,:)
		ELSE
			point=X(i,:)
		END IF
		DO j=1,Nobs
			dist=LP2((point-X(j,:))/bwidth)
			IF (ABS(dist).LT.ONE) THEN
				epa=0.75D0*(ONE-dist**2)
				denom=denom+epa
                numer=numer+epa*y(j)
			END IF
		END DO
		IF (denom.LT.tol) THEN
			y_fit(i)=ZERO
		ELSE
			y_fit(i)=numer/denom
		END IF
	END DO

END SUBROUTINE kernreg_kD_global
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


SUBROUTINE kernreg_1D_global(y,X,Xgrid,bandwidth,y_fit)

	!
	! Estimating a Kernel regression with:
	!   a. 1 dimensional covariatesa
	!   b. global bandwidth
	! E(y|X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to y_fit
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!
	! (this is a specialized version of kernreg_kD_global)

	REAL(dp), INTENT(IN), DIMENSION(:)					:: y
	REAL(dp), INTENT(IN), DIMENSION(:)					:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL						:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: y_fit

	! additional variables
	INTEGER(regInt)				:: Nobs !number of obsevations
	REAL(dp)					:: bwidth !the bandwidth vector we will use
	INTEGER(regInt)				:: i,j,N
	REAL(dp)					:: dist,denom,numer,epa,point
	REAL(dp), PARAMETER			:: tol=1.0D-06
	LOGICAL                     :: xgridFLAG

	Nobs=SIZE(y)
	IF (Nobs .NE. SIZE(X,1)) THEN !dimension mismatch
		y_fit=ZERO
		RETURN
	END IF

	IF (PRESENT(bandwidth)) THEN
		bwidth=bandwidth
	ELSE
		bwidth=silverman(Nobs,X)
	END IF

	!if Xgrid is an input
	IF (PRESENT(Xgrid)) THEN
		N=SIZE(Xgrid)
		xgridFLAG=.TRUE.
	ELSE
		N=Nobs
		xgridFLAG=.FALSE.
	END IF

	ALLOCATE(y_fit(N))

	DO i=1,N
		denom=ZERO
		numer=ZERO
		IF (xgridFLAG) THEN
			point=Xgrid(i)
		ELSE
			point=X(i)
		END IF

		DO j=1,Nobs
			dist=((point-X(j))/bwidth)
			IF (ABS(dist).LT.ONE) THEN
				epa=0.75D0*(ONE-dist**2)
				denom=denom+epa
                numer=numer+epa*y(j)
			END IF
		END DO
		IF (denom.LT.tol) THEN
			y_fit(i)=ZERO
		ELSE
			y_fit(i)=numer/denom
		END IF
	END DO
END SUBROUTINE kernreg_1D_global
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

SUBROUTINE kernreg_1D_global_sorted(y,X,Xgrid,bandwidth,y_fit)

	!
	! Estimating a Kernel regression with:
	!   a. 1 dimensional covariate
	!   b. global bandwidth
	! E(y|X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to y_fit
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!
	! WARNING: The procedure assumes that (y,X) is sorted by X in an ascending way!!
	!  If the data is not sorted, you'll get a wrong estimate!
	!	If Xgrid is present, Xgrid must be sorted in an ascending way as well!!
	!	If not, you'll get a wrong estimate!
	! (this is a specialized version of kernreg_1D_global)

	REAL(dp), INTENT(IN), DIMENSION(:)					:: y
	REAL(dp), INTENT(IN), DIMENSION(:)					:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL						:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: y_fit

	! additional variables
	INTEGER(regInt)				:: Nobs !number of observations
	REAL(dp)					:: bwidth !the bandwidth vector we will use
	INTEGER(regInt)				:: i,j,N
	REAL(dp)					:: dist,denom,numer,epa,point
	REAL(dp), PARAMETER			:: tol=1.0D-06
	INTEGER(regInt)				:: startIndex

	Nobs=SIZE(y)
	IF (Nobs .NE. SIZE(X,1)) THEN !dimension mismatch
		y_fit=ZERO
		RETURN
	END IF

	IF (PRESENT(bandwidth)) THEN
		bwidth=bandwidth
	ELSE
		bwidth=silverman(Nobs,X)
	END IF

	!if Xgrid is an input
	IF (PRESENT(Xgrid)) THEN
		N=SIZE(Xgrid)
	ELSE
		N=Nobs
	END IF

	ALLOCATE(y_fit(N))
	startIndex=1


outer:	DO i=1,N
		denom=ZERO
		numer=ZERO
		IF (PRESENT(Xgrid)) THEN
			point=Xgrid(i)
		ELSE
			point=X(i)
		END IF

inner:		DO j=startIndex,N
			dist=((X(j)-point)/bwidth)
			IF (dist.LT.-ONE) THEN
				startIndex=j  ! didn't reach the relevant points yet
			ELSE IF (dist.GT.ONE) THEN
				EXIT inner
			ELSE
				epa=0.75D0*(ONE-dist**2)
				denom=denom+epa
				numer=numer+epa*y(j)
			END IF
!end inner
		END DO inner
		IF (denom.LT.tol) THEN
			y_fit(i)=ZERO
		ELSE
			y_fit(i)=numer/denom
		END IF
!end outer
	END DO outer

END SUBROUTINE kernreg_1D_global_sorted
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

SUBROUTINE kernreg_kD_global_sorted(y,X,Xgrid,bandwidth,y_fit)

	!
	! Estimating a Kernel regression with:
	!   a. k dimensional covariates
	!   b. global bandwidth
	! E(y|X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to y_fit
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!
	! WARNING: The procedure assumes that (y,X) is sorted by X's first coordinate (column)
	!   in an ascending way!!
	!   If the data is not sorted, you'll get a wrong estimate!
	!	If Xgrid is present, Xgrid must be sorted by its first coordinate
	!   in an ascending way as well!!
	!	If not, you'll get a wrong estimate!
	! (this is a specialized version of kernreg_kD_global)

	REAL(dp), INTENT(IN), DIMENSION(:)					:: y
	REAL(dp), INTENT(IN), DIMENSION(:,:)				:: X
	REAL(dp), INTENT(IN), DIMENSION(:,:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: y_fit

	! additional variables
	INTEGER(regInt)				            :: Nobs !number of obsevations
	INTEGER(regInt)				            :: Ncov
	REAL(dp), ALLOCATABLE, DIMENSION(:)		:: bwidth !the bandwidth vector we will use
	REAL(dp), ALLOCATABLE, DIMENSION(:)		:: point
	INTEGER(regInt)				            :: i,j,N
	REAL(dp)					            :: dist1, dist,denom,numer,epa
	REAL(dp), PARAMETER			            :: tol=1.0D-06
	INTEGER(regInt)				            :: startIndex

	Nobs=SIZE(y) !number of observations
	IF (Nobs .NE. SIZE(X,1)) THEN !dimension mismatch
		y_fit=ZERO
		RETURN
	END IF

	Ncov=SIZE(X,2) !number of covariates
    ALLOCATE(bwidth(Ncov),point(Ncov))

	IF (PRESENT(bandwidth)) THEN
		bwidth=bandwidth
	ELSE
		bwidth=silverman(Nobs,X)
	END IF

	!if Xgrid is an input
	IF (PRESENT(Xgrid)) THEN
		N=SIZE(Xgrid,1)
	ELSE
		N=Nobs
	END IF

	ALLOCATE(y_fit(N))
	startIndex=1


outer:	DO i=1,N
		denom=ZERO
		numer=ZERO
		IF (PRESENT(Xgrid)) THEN
			point=Xgrid(i,:)
		ELSE
			point=X(i,:)
		END IF

inner:		DO j=startIndex,N
			dist1=((X(j,1)-point(1))/bwidth(1))
			IF (dist1.LT.-ONE) THEN
				startIndex=j
			ELSE IF (dist1.GT.ONE) THEN
				EXIT inner
			ELSE
				dist=LP2((point-X(j,:))/bwidth)
				epa=0.75D0*(ONE-dist**2)
				denom=denom+epa
				numer=numer+epa*y(j)
			END IF
!end inner
		END DO inner
		IF (denom.LT.tol) THEN
			y_fit(i)=ZERO
		ELSE
			y_fit(i)=numer/denom
		END IF
!end outer
	END DO outer

END SUBROUTINE kernreg_kD_global_sorted
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8




END MODULE kernelRegression