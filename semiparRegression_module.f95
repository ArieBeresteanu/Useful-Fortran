MODULE semiparametricRegression


USE linearRegression
USE kernelRegression
USE kernelDensity

IMPLICIT NONE




INTERFACE semiLinear
    MODULE PROCEDURE semiLinear_kD_kD_global
END INTERFACE

CONTAINS



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! 3. Local Linear Regression !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE localLinear_1D_global(y,X,Xgrid,bandwidth,y_fit,y_derivative)
	!
	! Estimating a local linear regression with a global bandwidth
	! E(y|X) is computed for each value of X and outputed to y_fit
	! the procedure also outputs the derivative of the regression function
	! The epanechnikov kernel is hard coded into the subroutine
	! This subroutine assumes that X is a one dimensional vector!!!
	!

	REAL(dp), INTENT(IN), DIMENSION(:)  			:: y
	REAL(dp), INTENT(IN), DIMENSION(:)  			:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL	:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL					:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE:: y_fit, y_derivative


	! additional variables
	INTEGER(regInt)			:: Nobs !number of obsevations
	REAL(dp)				:: bwidth !the bandwidth vector we will use
	INTEGER(regInt)			:: i,j,N
	REAL(dp)				:: dist,epa,s0, s1,s2, m0,m1,delta,point
	REAL(dp), PARAMETER		:: tol=1.0D-06
	REAL(dp), DIMENSION(2,2):: mat, mat_inv
	REAL(dp), DIMENSION(2)	:: v

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

	ALLOCATE(y_fit(N),y_derivative(N))

	DO i=1,N
		S0=ZERO
		S1=ZERO
		S2=ZERO
		m0=ZERO
		m1=ZERO
		IF (PRESENT(Xgrid)) THEN
			point=Xgrid(i)
		ELSE
			point=X(i)
		END IF

		DO j=1,Nobs
			delta=(X(j)-point)
			dist=(delta/bwidth)
			IF ( ABS(dist) .LT. ONE ) THEN
				epa=0.75D0*(ONE-dist**2)
			ELSE
				epa=ZERO
			END IF
			S0=S0+epa
			S1=S1+epa*delta
			S2=S2+epa*delta*delta
			m0=m0+y(j)*epa
			m1=m1+y(j)*epa*delta
		END DO
		mat=RESHAPE((/s0,s1,s1,s2/),(/2,2/))
		mat_inv=matinv2(mat)
		v=MATMUL(mat_inv,(/m0,m1/))
		y_fit(i)=v(1)
		y_derivative=v(2)
	END DO

END SUBROUTINE localLinear_1D_global
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! 4. Semi-linear Regression !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE semiLinear_kD_kD_global(y,Xlin,Xnonpar,bandwidth,trimming,verbose, &
                                    statsOut,y_fit,y_fit_lin,y_fit_nonpar,beta)
    !
    ! Semi-linear regression E(y|X)=Xlin'beta+f(Xnonpar)
    !   a. Xlin has Ncov_lin>1 covariates
    !   b. Xnonpar has Ncov_nonpar>1 covariates
    !   c. a global bandwidth is used for the nonparametric part
    !
    REAL(dp), INTENT(IN), DIMENSION(:)                  :: y
    REAL(dp), INTENT(IN), DIMENSION(:,:)                :: Xlin
    REAL(dp), INTENT(IN), DIMENSION(:,:)                :: Xnonpar
    REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL        :: bandwidth
    LOGICAL, INTENT(IN), OPTIONAL                       :: trimming
    LOGICAL, INTENT(IN), OPTIONAL                       :: verbose
    TYPE(linearStats), INTENT(OUT), OPTIONAL            :: statsOut
    REAL(dp), INTENT(OUT), DIMENSION(:)                 :: y_fit
    REAL(dp), INTENT(OUT), DIMENSION(:)                 :: y_fit_lin
    REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE    :: y_fit_nonpar
    REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE    :: beta

   	! additional variables
	INTEGER(regInt)				        :: Nobs !number of obsevations
	INTEGER(regInt)				        :: Ncov_lin !Number of linear covariates (No constant!)
	INTEGER(regInt)				        :: Ncov_nonpar !Number of nonparametric covariates

	REAL(dp), DIMENSION(:), ALLOCATABLE	    :: bwidth !the bandwidth vector we will use
	!REAL(dp), DIMENSION(:), ALLOCATABLE     :: beta_temp, trimming_indicator
	REAL(dp), DIMENSION(:), ALLOCATABLE     :: ytemp
	REAL(dp), DIMENSION(:,:), ALLOCATABLE   :: Xtemp
	INTEGER(regInt)				            :: i,j,N
    LOGICAL                                 :: verb

	Nobs=SIZE(y) !number of observations
	IF ( (Nobs .NE. SIZE(Xlin,1)) .OR. (Nobs .NE. SIZE(Xnonpar,1)) ) THEN !dimension mismatch
		y_fit=ZERO
		RETURN
	END IF

	Ncov_lin=SIZE(Xlin,2) !number of linear covariates
	Ncov_nonpar=SIZE(Xnonpar,2) !number of nonparametric covariates

	! Step 1: Nonparametric estimation of E(y|Xnonpar) and E(Xlin|Xnonpar)
	! ytemp=Y-E(Y|Xnonpar) and Xtemp=[ 1s Xlin-E(Xlin|Xnonpar) ]
	ALLOCATE(beta(Ncov_lin)) !,beta_temp(Ncov_lin+1),Xtemp(Nobs,Ncov_lin+1))
	ALLOCATE(ytemp(Nobs),Xtemp(Nobs,Ncov_lin),bwidth(Ncov_nonpar))
    verb=.FALSE.
    IF (PRESENT(verbose) .AND. verbose) verb=.TRUE.
	IF (PRESENT(bandwidth) .AND. (SIZE(bandwidth,1) .EQ. Ncov_nonpar)) THEN
		bwidth=bandwidth
	ELSE !use silverman rule of thumb
		DO i=1,Ncov_nonpar
			bwidth(i)=silverman(Nobs,Xnonpar(:,i))
		END DO
	END IF

    !preparing the ytemp and Xtemp arrays
    !Xtemp(:,1)=ONE
    DO i=1,Ncov_lin
       !CALL kernreg(y,X,Xgrid,bandwidth,y_fit)
       CALL kernreg(y=Xlin(:,i),X=Xnonpar,y_fit=ytemp,bandwidth=bwidth)
       Xtemp(:,i)=Xlin(:,i)-ytemp !residual
       IF (verb) PRINT *,"Done for covariate ",i, "out of ", Ncov_lin
    END DO
    CALL kernreg(y=y,X=Xnonpar,y_fit=ytemp,bandwidth=bwidth)
    ytemp=y-ytemp !residual
    IF (verb) PRINT *,"Done for y"

    ! Step 2: Linear regression of ytemp on Xtemp:
    !OLS(Nobs,Ncov,y,X,y_fit,beta,R2,robustSE,verbose)

    IF (PRESENT(statsOut)) THEN
            CALL OLS(Nobs=Nobs,Ncov=Ncov_lin,y=ytemp,X=Xtemp,y_fit=y_fit_lin,beta=beta,statsOut=statsOut,verbose=verb)
        ELSE
            CALL OLS(Nobs=Nobs,Ncov=Ncov_lin,y=ytemp,X=Xtemp,y_fit=y_fit_lin,beta=beta,verbose=verb)
    END IF



    !Step 3: Calculating the nonparametric component
    ! Regress y-Xlin'beta on Xnonpar
    ytemp=y
    DO i=1,Ncov_lin
        ytemp=ytemp-y_fit_lin
    END DO

    !step 4: computing the nonparametric part
    CALL kernreg(y=ytemp,X=Xnonpar,y_fit=y_fit_nonpar,bandwidth=bwidth)

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
END SUBROUTINE semiLinear_kD_kD_global


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! 5. Additive separable Regression !!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE addSepReg(Nobs,Klin,Knonpar,Niter,y,Xlin,constIncluded,&
						Xnonpar,bandwidth,y_fit_lin,y_fit_nonpar,y_fit,beta)

	!
	! Additive separable regression. y is the outcome variable, Xlin is the covariate matrix
	! of the linear part (include a constant!) and Xnonlin is the covariate matrix of the nonparametric
	!
	INTEGER(regInt), INTENT(IN)						:: Nobs 	! number of observations
	INTEGER(regInt), INTENT(IN)						:: Klin 	! number of covariates in the linear part
	INTEGER(regInt), INTENT(IN)						:: Knonpar 	! number of covariates in the nonparametric part
	INTEGER(regInt), INTENT(IN)						:: Niter 	! number of iterations
	REAL(dp), INTENT(IN), DIMENSION(:)				:: y 		! response variable
	REAL(dp), INTENT(IN), DIMENSION(:,:)			:: Xlin		! the matrix of linear covariates
	LOGICAL, INTENT(IN),OPTIONAL					:: constIncluded ! .TRUE. if Xlin has a column of 1's
	REAL(dp), INTENT(IN), DIMENSION(:,:)			:: Xnonpar 	! the matrix of nonparametric covariates
	REAL(dp), INTENT(IN), DIMENSION(:),OPTIONAL		:: bandwidth 	! optional bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:)				:: y_fit_lin 	! output linear fit
	REAL(dp), INTENT(OUT), DIMENSION(:,:)			:: y_fit_nonpar	! output nonparametric fits (one column per each covariate)
	REAL(dp), INTENT(OUT), DIMENSION(:)				:: y_fit		! output total fit
	REAL(dp), INTENT(OUT), DIMENSION(:)				:: beta			! output beta of linear form

	!additional variables
	INTEGER(regInt)						:: i,j,n
	REAL(dp), DIMENSION(Nobs)			:: y_resid
	REAL(dp), DIMENSION(Nobs)			:: y_nonpar
	REAL(dp), DIMENSION(:),ALLOCATABLE	:: y_fit_temp
	REAL(dp), DIMENSION(Klin)			:: mu_lin,sigma_lin
	REAL(dp), DIMENSION(Knonpar)		:: mu_nonpar,sigma_nonpar
	REAL(dp), DIMENSION(Nobs,Klin)		:: Xlin_st
	REAL(dp), DIMENSION(Nobs,Knonpar)	:: Xnonpar_st
	REAL(dp), DIMENSION(Nobs,1)			:: y_st
	REAL(dp), DIMENSION(1)				:: mu_y,sigma_y
	REAL(dp)							:: avg_y, meanfit
	TYPE(linearStats)                   :: statsOut

	!studentizing the variables
	!CALL studentize(Nobs,Klin,Xlin,Xlin_st,mu_lin,sigma_lin)
	Xlin_st=Xlin
	CALL studentize(Nobs,Knonpar,Xnonpar,Xnonpar_st,mu_nonpar,sigma_nonpar)

	IF (PRESENT(constIncluded) .AND. constIncluded) THEN
		y_st(:,1)=y
		avg_y=SUM(y)/Nobs
	ELSE
		CALL studentize(Nobs,1,y,y_st,mu_y,sigma_y)
		avg_y=0.0_dp
	END IF

	!initializing
	y_fit_nonpar=ZERO
	y_nonpar=ZERO
	ALLOCATE(y_fit_temp(Nobs))

	DO i=1,Niter
	PRINT *,"iteration ", i
		y_resid=y_st(:,1)-y_nonpar !substract the nonparametric fit so far
		CALL OLS(Nobs=Nobs,Ncov=Klin,y=y_resid,X=Xlin_st,y_fit=y_fit_lin,beta=beta,statsOut=statsOut) !the linear fit

		PRINT '(9F16.8)',statsOut%R2
		DO j=1,Knonpar
			y_resid=y_st(:,1)-y_fit_lin-y_nonpar+y_fit_nonpar(:,j) !leaving only the residual from the j`th nonparametric covariate
			!univariate nonparametric regression
			CALL kernreg(y=y_resid,X=Xnonpar_st(:,j),y_fit=y_fit_temp)
			meanfit=SUM(y_fit_temp)/Nobs
			y_fit_nonpar(:,j)=y_fit_temp-meanfit

		END DO
		y_nonpar=SUM(y_fit_nonpar,DIM=2) !total fit coming from the nonparametric part
	END DO
	y_fit=y_nonpar+y_fit_lin
	statsOut%R2=SUM((y_fit-avg_y)**2,DIM=1)/SUM((y-avg_y)**2,DIM=1)
	PRINT '(F16.8)',statsOut%R2
END SUBROUTINE addSepReg
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

END MODULE semiparametricRegression