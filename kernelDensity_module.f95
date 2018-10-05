MODULE kernelDensity

USE regression



IMPLICIT NONE


INTERFACE kerndensity
    MODULE PROCEDURE kerndensity_kD_global
    MODULE PROCEDURE kerndensity_1D_global
END INTERFACE

INTERFACE kerndensity_sorted
    MODULE PROCEDURE kerndensity_kD_global_sorted
    MODULE PROCEDURE kerndensity_1D_global_sorted
END INTERFACE

CONTAINS


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! Kernel Regression !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE kerndensity_kD_global(X,Xgrid,bandwidth,X_density)
	!
	! Estimating a Kernel density with:
	!   a. k dimensional X
	!   b. global bandwidth
	! f(X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to X_density
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!   X and Xgrid (if exists) are not necessarily sorted

	REAL(dp), INTENT(IN),  DIMENSION(:,:)				:: X
	REAL(dp), INTENT(IN),  DIMENSION(:,:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN),  DIMENSION(:),   OPTIONAL		:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: X_density

	! additional variables
	REAL(dp), DIMENSION(:), ALLOCATABLE	:: bwidth !the bandwidth vector we will use
	REAL(dp), DIMENSION(:), ALLOCATABLE	:: point
	INTEGER(regInt)				        :: i,j,N,Nvar,Nobs
	REAL(dp)					        :: dist
    REAL(dp)                            :: denom
	LOGICAL                             :: xgridFLAG

    Nobs=SIZE(X,1)
	Nvar=SIZE(X,2) 
    
	ALLOCATE(bwidth(Nvar),point(Nvar))

	IF (PRESENT(bandwidth) .AND. (SIZE(bandwidth) .EQ. Nvar)) THEN
		bwidth=bandwidth
	ELSE !use silverman rule of thumb
		DO i=1,Nvar
			bwidth(i)=silverman(Nobs,X(:,i))
		END DO
	END IF

	!if Xgrid is an input
	IF (PRESENT(Xgrid)) THEN
		N=SIZE(Xgrid)
		xgridFLAG=.TRUE.
	ELSE
		N=Nobs
		xgridFLAG=.FALSE.
	END IF

	ALLOCATE(X_density(N))
    X_density=ZERO
    
    denom=Nobs*LP2(bwidth)
	!if Xgrid is present and its 2nd dimention matches that of the covariates,
	DO i=1,N
		IF (xgridFLAG) THEN
			point=Xgrid(i,:)
		ELSE
			point=X(i,:)
		END IF
		DO j=1,Nobs
			dist=LP2((point-X(j,:))/bwidth)
			IF (ABS(dist).LT.ONE) X_density(i)=X_density(i)+0.75D0*(ONE-dist**2)
		END DO
		X_density(i)=X_density(i)/denom
	END DO

END SUBROUTINE kerndensity_kD_global
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


SUBROUTINE kerndensity_1D_global(X,Xgrid,bandwidth,X_density)
	!
	! Estimating a Kernel density with:
	!   a. 1 dimensional X
	!   b. global bandwidth
	! f(X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to X_density
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!   X and Xgrid (if exists) are not necessarily sorted
	!
	! (this is a specialized version of kerndensity_kD_global)

	REAL(dp), INTENT(IN), DIMENSION(:)					:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL						:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: X_density

	! additional variables
	INTEGER(regInt)				:: Nobs !number of obsevations
	REAL(dp)					:: bwidth !the bandwidth vector we will use
	INTEGER(regInt)				:: i,j,N
	REAL(dp)					:: dist,point,denom
	LOGICAL                     :: xgridFLAG

	Nobs=SIZE(X)
    

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

	ALLOCATE(X_density(N))
    X_density=ZERO
    
    denom=Nobs*bwidth
	DO i=1,N
		IF (xgridFLAG) THEN
			point=Xgrid(i)
		ELSE
			point=X(i)
		END IF

		DO j=1,Nobs
			dist=((point-X(j))/bwidth)
			IF (ABS(dist).LT.ONE) X_density(i)=X_density(i)+0.75D0*(ONE-dist**2)
		END DO
		X_density(i)=X_density(i)/denom
	END DO
END SUBROUTINE kerndensity_1D_global
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

SUBROUTINE kerndensity_1D_global_sorted(X,Xgrid,bandwidth,X_density)

	!
	! Estimating a Kernel density with:
	!   a. 1 dimensional covariate
	!   b. global bandwidth
	! f(X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to X_density
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!
	! WARNING: The procedure assumes that (X) is sorted in an ascending way!!
	!  If the data is not sorted, you'll get a wrong estimate!
	!	If Xgrid is present, Xgrid must be sorted in an ascending way as well!!
	!	If not, you'll get a wrong estimate!
	! (this is a specialized version of kerndensity_1D_global)

	REAL(dp), INTENT(IN), DIMENSION(:)					:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL						:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: X_density

	! additional variables
	INTEGER(regInt)				:: Nobs !number of observations
	REAL(dp)					:: bwidth !the bandwidth vector we will use
	INTEGER(regInt)				:: i,j,N
	REAL(dp)					:: dist,point,denom
	INTEGER(regInt)				:: startIndex
    LOGICAL                     :: xgridFLAG

	Nobs=SIZE(X)

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

	ALLOCATE(X_density(N))
    X_density=ZERO
	startIndex=1

denom=Nobs*bwidth
outer:	DO i=1,N
		IF (xgridFLAG) THEN
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
				X_density(i)=X_density(i)+0.75D0*(ONE-dist**2)
			END IF
!end inner
		END DO inner
		X_density(i)=X_density(i)/denom
!end outer
	END DO outer

END SUBROUTINE kerndensity_1D_global_sorted
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

SUBROUTINE kerndensity_kD_global_sorted(X,Xgrid,bandwidth,X_density)

	!
	! Estimating a Kernel density with:
	!   a. k dimensional X
	!   b. global bandwidth
	! f(X) is computed for each value of X or Xgrid (if present)
	! The estimator is outputed to X_density
	! The epanechnikov kernel is hard coded into the subroutine for speed
	!
	! WARNING: The procedure assumes that X is sorted by X's first coordinate (column)
	!   in an ascending way!!
	!   If the data is not sorted, you'll get a wrong estimate!
	!	If Xgrid is present, Xgrid must be sorted by its first coordinate
	!   in an ascending way as well!!
	!	If not, you'll get a wrong estimate!
	! (this is a specialized version of kerndensity_kD_global)

	REAL(dp), INTENT(IN), DIMENSION(:,:)				:: X
	REAL(dp), INTENT(IN), DIMENSION(:,:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: X_density

	! additional variables
	INTEGER(regInt)				            :: Nobs !number of obsevations
	INTEGER(regInt)				            :: Nvar
	REAL(dp), ALLOCATABLE, DIMENSION(:)		:: bwidth !the bandwidth vector we will use
	REAL(dp), ALLOCATABLE, DIMENSION(:)		:: point
	INTEGER(regInt)				            :: i,j,N
	REAL(dp)					            :: dist1,dist,denom
	INTEGER(regInt)				            :: startIndex
    LOGICAL                                 :: xgridFLAG

	Nobs=SIZE(X,1) !number of observations
	Nvar=SIZE(X,2) !number of covariates
    ALLOCATE(bwidth(Nvar),point(Nvar))

	IF (PRESENT(bandwidth)) THEN
		bwidth=bandwidth
	ELSE
		bwidth=silverman(Nobs,X)
	END IF

	!if Xgrid is an input
	IF (PRESENT(Xgrid)) THEN
		N=SIZE(Xgrid,1)
        xgridFLAG=.TRUE.
	ELSE
		N=Nobs
        xgridFLAG=.FALSE.
	END IF

	ALLOCATE(X_density(N))
    X_density=ZERO
	startIndex=1

denom=Nobs*LP2(bwidth)
outer:	DO i=1,N
		IF (xgridFLAG) THEN
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
				X_density(i)=X_density(i)+0.75D0*(ONE-dist**2)
			END IF
!end inner
		END DO inner
		X_density(i)=X_density(i)/denom
!end outer
	END DO outer

END SUBROUTINE kerndensity_kD_global_sorted
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


END MODULE kernelDensity