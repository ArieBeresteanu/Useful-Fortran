MODULE regression

USE utilities, ONLY: LP2,printRealMatrix
USE linAlgebra
IMPLICIT NONE

REAL(dp), PARAMETER		:: ZERO=0.00_dp
REAL(dp), PARAMETER		:: ONE=1.00_dp

CONTAINS

	SUBROUTINE OLS(Nobs,Ncov,y,X,y_fit,beta,R2)
	
	!
	! a simple OLS regression assuming the covariates matrix includes a constant
	!
	
	INTEGER, INTENT(IN)							:: Nobs !number of obsevations
	INTEGER, INTENT(IN)							:: Ncov !number of covariates +1 for the constant
	REAL(dp), DIMENSION(Nobs), INTENT(IN)		:: y !vector of outcome variables
	REAL(dp), DIMENSION(Nobs,Ncov), INTENT(IN)	:: X !covariates matrix(assumes the constant is already included)
	REAL(dp), DIMENSION(Nobs), INTENT(OUT)		:: y_fit !fited values	
	REAL(dp), DIMENSION(Ncov), INTENT(OUT)		:: beta !beta coefficient
	REAL(dp), INTENT(OUT), OPTIONAL				:: R2

	
	!additional variables
	REAL(dp), DIMENSION(Ncov,Ncov)	:: XX, XXinv
	REAL(dp), DIMENSION(Ncov,Nobs)	:: Xtrans
	REAL(dp), DIMENSION(Ncov)		:: Xy
	INTEGER							:: i,j
	REAL(dp)						:: avg_y
	
	
	Xtrans=TRANSPOSE(X)
	XX=MATMUL(Xtrans,X)
	Xy=MATMUL(Xtrans,y)
	CALL matinv(Ncov,XX,XXinv)
	beta=MATMUL(XXinv,Xy)
	y_fit=PACK(MATMUL(X,beta),.TRUE.)
	
	IF (PRESENT(R2)) THEN
		avg_y=SUM(y)/Nobs
		R2=SUM((y_fit-avg_y)**2,DIM=1)/SUM((y-avg_y)**2,DIM=1)
	END IF
	
	END SUBROUTINE OLS
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8	
	
	SUBROUTINE kernreg(Nobs,Ncov,y,X,Xgrid,bandwidth,y_fit)
	
	!
	! Estimating a simple Kernel regression with a global bandwidth
	! E(y|X) is computed for each value of X and outputed to y_fit
	! The epanechnikov kernel is hard coded into the subroutine
	!
	
	INTEGER, INTENT(IN)									:: Nobs !number of obsevations
	INTEGER, INTENT(IN)									:: Ncov !Number of covariates (Note: do not include a constant)
	REAL(dp), INTENT(IN), DIMENSION(Nobs)				:: y
	REAL(dp), INTENT(IN), DIMENSION(Nobs,Ncov)			:: X
	REAL(dp), INTENT(IN), DIMENSION(:,:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), DIMENSION(:),OPTIONAL			:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: y_fit
	
	! additional variables
	REAL(dp), DIMENSION(Ncov)	:: bwidth !the bandwidth vector we will use
	INTEGER						:: i,j,N
	REAL(dp)					:: dist,denom,numer,epa
	REAL(dp), PARAMETER			:: tol=1.0D-06
	
	IF (PRESENT(bandwidth)) THEN
		IF (SIZE(bandwidth) .EQ. Ncov) THEN
			bwidth=bandwidth
		ELSE
			bwidth=bandwidth(1)
		END IF
	ELSE !use silverman rule of thumb
		DO i=1,Ncov
			bwidth(i)=silverman(Nobs,X(:,i))
		END DO
	END IF	
		
	!if Xgrid is an input and its 2nd dimention matches that of the covariates,
	IF (PRESENT(Xgrid) .AND. SIZE(Xgrid,2).EQ.Ncov) THEN
		N=SIZE(Xgrid,1)
		ALLOCATE(y_fit(N))
		DO i=1,N
			denom=ZERO
			numer=ZERO
			DO j=1,Nobs
				dist=LP2(Xgrid(i,:)-X(j,:))/bwidth(j)
				IF (ABS(dist).LT.ONE) THEN
					epa=0.75D0*(ONE-dist**2)
				ELSE
					epa=ZERO
				END IF
				denom=denom+epa
				numer=numer+epa*y(j)
			END DO
			IF (denom.LT.tol) THEN
				y_fit(i)=ZERO
			ELSE
				y_fit(i)=numer/denom
			END IF
		END DO
	!if Xgrid is not present or not of the right dimension,
	ELSE
		ALLOCATE(y_fit(Nobs))
		DO i=1,Nobs
			denom=ZERO
			numer=ZERO
			DO j=1,Nobs
				dist=LP2(X(i,:)-X(j,:))/bwidth(j)
				IF (ABS(dist).LT.ONE) THEN
					epa=0.75D0*(ONE-dist**2)
				ELSE
					epa=ZERO
				END IF
				denom=denom+epa
				numer=numer+epa*y(j)
			END DO
			IF (denom.LT.tol) THEN
				y_fit(i)=ZERO
			ELSE
				y_fit(i)=numer/denom
			END IF
		END DO
	END IF
	
	END SUBROUTINE kernreg
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


SUBROUTINE kernreg1D(Nobs,y,X,Xgrid,bandwidth,y_fit)
	
	!
	! Estimating a simple Kernel regression with a global bandwidth
	! E(y|X) is computed for each value of X and outputed to y_fit
	! The epanechnikov kernel is hard coded into the subroutine
	! This subroutine assumes that X is a one dimensional vector!!!
	!
	
	INTEGER, INTENT(IN)									:: Nobs !number of obsevations
	REAL(dp), INTENT(IN), DIMENSION(Nobs)				:: y
	REAL(dp), INTENT(IN), DIMENSION(Nobs)				:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL		:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL						:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE	:: y_fit
	
	! additional variables
	REAL(dp)				:: bwidth !the bandwidth vector we will use
	INTEGER					:: i,j,N
	REAL(dp)				:: dist,denom,numer,epa,point
	REAL(dp), PARAMETER		:: tol=1.0D-06
	
	
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
	
	DO i=1,N
		denom=ZERO
		numer=ZERO
		IF (PRESENT(Xgrid)) THEN
			point=Xgrid(i)
		ELSE
			point=X(i)
		END IF
			
		DO j=1,Nobs
			dist=((point-X(j))/bwidth)
			IF (ABS(dist).LT.ONE) THEN
				epa=0.75D0*(ONE-dist**2)
			ELSE
				epa=ZERO
			END IF
			denom=denom+epa
			numer=numer+epa*y(j)
		END DO
		IF (denom.LT.tol) THEN
			y_fit(i)=ZERO
		ELSE
			y_fit(i)=numer/denom
		END IF
	END DO
	END SUBROUTINE kernreg1D

	SUBROUTINE localLinear1D(Nobs,y,X,Xgrid,bandwidth,y_fit,y_derivative)
	!
	! Estimating a local linear regression with a global bandwidth
	! E(y|X) is computed for each value of X and outputed to y_fit
	! the procedure also outputs the derivative of the regression function
	! The epanechnikov kernel is hard coded into the subroutine
	! This subroutine assumes that X is a one dimensional vector!!!
	!
	
	INTEGER, INTENT(IN)								:: Nobs !number of obsevations
	REAL(dp), INTENT(IN), DIMENSION(Nobs)			:: y
	REAL(dp), INTENT(IN), DIMENSION(Nobs)			:: X
	REAL(dp), INTENT(IN), DIMENSION(:), OPTIONAL	:: Xgrid
	REAL(dp), INTENT(IN), OPTIONAL					:: bandwidth
	REAL(dp), INTENT(OUT), DIMENSION(:), ALLOCATABLE:: y_fit, y_derivative
	
	
	! additional variables
	REAL(dp)				:: bwidth !the bandwidth vector we will use
	INTEGER					:: i,j,N
	REAL(dp)				:: dist,epa,s0, s1,s2, m0,m1,delta
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
	
	END SUBROUTINE localLinear1D
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


	SUBROUTINE studentize(rows,columns,X_in,X_out,mu,sigma)

	!
	!   ===================================================================
	!   Studentizing the variables in the matrix X_in buy substructing the 
	!   mean and then dividing by the standard deviation of each column.
	!   ===================================================================
	!

	INTEGER, INTENT(IN)							:: rows,columns
	REAL(dp), DIMENSION(rows,columns),INTENT(IN)	:: X_in
	REAL(dp), DIMENSION(rows,columns),INTENT(OUT):: X_out
	REAL(dp), DIMENSION(columns), INTENT(OUT)	:: mu,sigma

	INTEGER				:: i	
	REAL(dp)				:: rN

	rN=REAL(rows,kind=8)
	
	mu=SUM(X_in)/rN
	
	DO i=1,columns
		sigma(i)=SQRT(SUM((x_in(:,i)-mu(i))**2)/(rN-ONE))
		X_out(:,i)=(X_in(:,i)-mu(i))/sigma(i)
	END DO

	END SUBROUTINE studentize
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

	
	SUBROUTINE addSepReg(Nobs,Klin,Knonpar,Niter,y,Xlin,constIncluded,&
						Xnonpar,bandwidth,y_fit_lin,y_fit_nonpar,y_fit,beta)
		
	!
	! Additive separable regression. y is the outcome variable, Xlin is the covariate matrix
	! of the linear part (include a constant!) and Xnonlin is the covariate matrix of the nonparametric
	!
	INTEGER, INTENT(IN)								:: Nobs 	! number of observations
	INTEGER, INTENT(IN)								:: Klin 	! number of covariates in the linear part
	INTEGER, INTENT(IN)								:: Knonpar 	! number of covariates in the nonparametric part
	INTEGER, INTENT(IN)								:: Niter 	! number of iterations
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
	INTEGER								:: i,j,n
	REAL(dp), DIMENSION(Nobs)			:: y_resid
	REAL(dp), DIMENSION(Nobs)			:: y_nonpar
	REAL(dp), DIMENSION(:),ALLOCATABLE	:: y_fit_temp
	REAL(dp), DIMENSION(Klin)			:: mu_lin,sigma_lin
	REAL(dp), DIMENSION(Knonpar)		:: mu_nonpar,sigma_nonpar
	REAL(dp), DIMENSION(Nobs,Klin)		:: Xlin_st
	REAL(dp), DIMENSION(Nobs,Knonpar)	:: Xnonpar_st
	REAL(dp), DIMENSION(Nobs,1)			:: y_st
	REAL(dp), DIMENSION(1)				:: mu_y,sigma_y
	REAL(dp)							:: R2, avg_y, meanfit
	
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
		CALL OLS(Nobs,Klin,y_resid,Xlin_st,y_fit_lin,beta,R2=R2) !the linear fit
		
		PRINT '(9F16.8)',R2
		!$OMP PARALLEL DO
		DO j=1,Knonpar
			y_resid=y_st(:,1)-y_fit_lin-y_nonpar+y_fit_nonpar(:,j) !leaving only the residual from the j`th nonparametric covariate
			!univariate nonparametric regression
			CALL kernreg1D(Nobs=Nobs,y=y_resid,X=Xnonpar_st(:,j),y_fit=y_fit_temp)
			meanfit=SUM(y_fit_temp)/Nobs
			y_fit_nonpar(:,j)=y_fit_temp-meanfit
			
		END DO
		!$OMP END PARALLEL DO
		y_nonpar=SUM(y_fit_nonpar,DIM=2) !total fit coming from the nonparametric part
	END DO	
	y_fit=y_nonpar+y_fit_lin
	R2=SUM((y_fit-avg_y)**2,DIM=1)/SUM((y-avg_y)**2,DIM=1)
	PRINT '(F16.8)',R2	
	END SUBROUTINE addSepReg
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
	
	
	REAL(dp) FUNCTION silverman(N,vec)
	!
	! Computing the Silverman rule of thumb bandwidth
	!
	INTEGER, INTENT(IN)		:: N
	REAL(dp), INTENT(IN)	:: vec(N)
	
	REAL(dp)		:: mu,sigma,rN
	
	rN=N+0.0_dp
	mu=SUM(vec)/rN
	sigma=SQRT(SUM((vec-mu)**2)/(rN-ONE))
	
	silverman = 1.060_dp*sigma*(N**(-0.20_dp))
	END FUNCTION silverman
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

	SUBROUTINE makeGrid(X,mesh,padding,Xgrid)
		REAL(dp), INTENT(IN)			:: X(:)
		INTEGER, INTENT(IN)				:: mesh
		LOGICAL, INTENT(IN), OPTIONAL	:: padding
		REAL(dp), INTENT(OUT)			:: Xgrid(mesh)
		
		INTEGER	:: i 
		REAL(dp)	:: low,high,length,skip,pad
		REAL(dp), PARAMETER	:: tol=1.0D-06
		REAL(dp), PARAMETER	:: padding_percent=0.01D0
		
		high=MAXVAL(X)
		low =MINVAL(X)
		length=high-low
		IF (length.LT.tol .OR. mesh.EQ.1) THEN
			Xgrid=high
		ELSE
			pad=0.0D0
			IF (PRESENT(padding)) pad=padding_percent*length
				
			skip=(high-low-2.0D0*pad)/(mesh-1)
			DO i=1,mesh
				Xgrid(i)=low+pad+skip*(i-1)
			END DO
		END IF
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
	END SUBROUTINE makeGrid	
	
	
END MODULE regression