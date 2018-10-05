MODULE linearRegression

USE utilities, ONLY: printRealMatrix, printRealVector,remove_dups
USE statistics, ONLY: freq
USE linAlgebra
IMPLICIT NONE

TYPE linearStats
    INTEGER(regInt)                     :: Npars
    REAL(dp)                            :: R2
    REAL(dp), DIMENSION(:,:), POINTER   :: robustSE
END TYPE



CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Linear Least Squares     !!!!!
    !!!!!   1.1 OLS                !!!!!
    !!!!!   1,2 WLS                !!!!!
    !!!!!   1.3 createFixedEffects !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE OLS(Nobs,Ncov,y,X,y_fit,beta,statsOut,verbose)
	!
	! a simple OLS regression assuming the covariates matrix includes a constant
	!

	INTEGER(regInt), INTENT(IN)							    :: Nobs !number of obsevations
	INTEGER(regInt), INTENT(IN)							    :: Ncov !number of covariates +1 for the constant
	REAL(dp), DIMENSION(Nobs), INTENT(IN)		            :: y !vector of outcome variables
	REAL(dp), DIMENSION(Nobs,Ncov), INTENT(IN)	            :: X !covariates matrix(assumes the constant is already included)
	REAL(dp), DIMENSION(Nobs), INTENT(OUT)		            :: y_fit !fited values
	REAL(dp), DIMENSION(Ncov), INTENT(OUT)		            :: beta !beta coefficient
	TYPE(linearStats), INTENT(OUT), OPTIONAL                :: statsOut
	LOGICAL, INTENT(IN), OPTIONAL                           :: verbose


	!additional variables
	REAL(dp), DIMENSION(Ncov,Ncov)	:: XX, XXinv,temp
	REAL(dp), DIMENSION(Ncov,Nobs)	:: Xtrans
	REAL(dp), DIMENSION(Ncov)		:: Xy
	INTEGER(regInt)					:: i,j
	REAL(dp)						:: avg_y, s
	REAL(dp), DIMENSION(Nobs)       :: eps2


	Xtrans=TRANSPOSE(X)
	XX=MATMUL(Xtrans,X)
	Xy=MATMUL(Xtrans,y)
	CALL matinv(Ncov,XX,XXinv)
	beta=MATMUL(XXinv,Xy)
	y_fit=PACK(MATMUL(X,beta),.TRUE.)

	IF (PRESENT(statsOut)) THEN
        ALLOCATE (statsOut%robustSE(Ncov,Ncov))
		avg_y=SUM(y)/Nobs
		statsOut%R2=SUM((y_fit-avg_y)**2,DIM=1)/SUM((y-avg_y)**2,DIM=1)

        eps2=(y-y_fit)**2
        DO i=1,Ncov
            DO j=1,Ncov
                temp(i,j)=DOT_PRODUCT(X(:,i)*eps2,X(:,j))
            END DO
        END DO
        statsOut%robustSE=MATMUL(XXinv,MATMUL(temp,XXinv))
        statsOUT%Npars=Ncov

        IF (PRESENT(verbose) .AND. verbose ) THEN
            PRINT *," "
            PRINT *,"----------------------------------------"
            PRINT *,"       beta            (SE)         t "
            PRINT *,"----------------------------------------"
            DO i=1,Ncov
                s=SQRT(statsOut%robustSE(i,i))
                PRINT '(2F16.8,F8.3)', beta(i), s, beta(i)/s
            END DO
            PRINT *,"----------------------------------------"
            PRINT '(A10,F8.6)',"R^2 = ",statsOut%R2
            PRINT *,"----------------------------------------"
            PRINT *," "
        END IF

    END IF


END SUBROUTINE OLS
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8


SUBROUTINE WLS(Nobs,Ncov,y,X,weights,y_fit,beta,statsOut,verbose)
	!
	! weighted LS regression assuming the covariates matrix includes a constant
	!

	INTEGER(regInt), INTENT(IN)							    :: Nobs !number of obsevations
	INTEGER(regInt), INTENT(IN)							    :: Ncov !number of covariates +1 for the constant
	REAL(dp), DIMENSION(Nobs), INTENT(IN)		            :: y !vector of outcome variables
	REAL(dp), DIMENSION(Nobs,Ncov), INTENT(IN)	            :: X !covariates matrix(assumes the constant is already included)
	REAL(dp), DIMENSION(Nobs), INTENT(IN)                   :: weights ! the weights vector
	REAL(dp), DIMENSION(Nobs), INTENT(OUT)		            :: y_fit !fited values
	REAL(dp), DIMENSION(Ncov), INTENT(OUT)		            :: beta !beta coefficient
	TYPE(linearStats), INTENT(OUT), OPTIONAL                :: statsOut
	LOGICAL, INTENT(IN), OPTIONAL                           :: verbose


	!additional variables
	REAL(dp), DIMENSION(Nobs,Ncov)  :: Xw
	REAL(dp), DIMENSION(Nobs)       :: Yw, sqrtWeights
	INTEGER(regInt)					:: i
	REAL(dp)						:: s

    ! rewight the X matrix and y vector
    sqrtWeights=SQRT(weights)
    DO i=1,Nobs
        Yw(i)=y(i)*sqrtWeights(i)
        Xw(i,:)=X(i,:)*sqrtWeights(i)
    END DO

    !now it's an OLS of Yw on Xw:

	IF (PRESENT(statsOut)) THEN
        ALLOCATE (statsOut%robustSE(Ncov,Ncov))
        CALL OLS(Nobs=Nobs,Ncov=Ncov,y=Yw,X=Xw,y_fit=y_fit,beta=beta,statsOut=statsOut)

        IF (PRESENT(verbose) .AND. verbose ) THEN
            PRINT *," "
            PRINT *,"----------------------------------------"
            PRINT *,"       beta            (SE)         t "
            PRINT *,"----------------------------------------"
            DO i=1,Ncov
                s=SQRT(statsOut%robustSE(i,i))
                PRINT '(2F16.8,F8.3)', beta(i), s, beta(i)/s
            END DO
            PRINT *,"----------------------------------------"
            PRINT '(A10,F8.6)',"R^2 = ",statsOut%R2
            PRINT *,"----------------------------------------"
            PRINT *," "
        END IF

    ELSE
        CALL OLS(Nobs=Nobs,Ncov=Ncov,y=Yw,X=Xw,y_fit=y_fit,beta=beta)
    END IF


END SUBROUTINE WLS
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8




FUNCTION createFixedEffects(N,x,threshold,verbose) RESULT(FE)
	!
	! The vector x of integers is used to create fixed effects.
	! First, the unique values in x are detected and their frequency.
	! Then a binary matrix of the appropriate dimension is created.
	! By default the first unique value is used as the benchmark category.
	! Therefore, one less fixed effect is created than the number of
	! unique values in x.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! threshold is an optional input. FEs will be created only for values
	! that appear at least threshold times in the data. It is to avoid
	! creating a FE for a value that appears only once or a very low number
	! of times in the data.
    !
		INTEGER(regInt), INTENT(IN)						:: N !number of observations
		INTEGER(regInt), INTENT(IN), DIMENSION(N)		:: x !a discrete covariate
        INTEGER(regInt), INTENT(IN), OPTIONAL           :: threshold
		LOGICAL, OPTIONAL, INTENT(IN)                   :: verbose
		REAL(dp), ALLOCATABLE, DIMENSION(:,:)	        :: FE !the output matrix of fixed effects

		INTEGER(regInt), DIMENSION(:), ALLOCATABLE		:: unique,f
		LOGICAL, DIMENSION(:), ALLOCATABLE              :: above_threshold

		INTEGER(regInt)									:: i,k,j

        !frequency table for x
        CALL freq(N=N,list=x,unique=unique,frequency=f,verbose=.FALSE.)
!

        ALLOCATE(above_threshold(size(unique)))
        IF (PRESENT(threshold)) THEN
            above_threshold=.FALSE.
            k=0
            DO i=1,size(unique)
                IF (f(i) .GE. threshold) THEN
                    above_threshold(i)=.TRUE.
                    k=k+1
                END IF
            END DO
            ALLOCATE(FE(N,k-1))
            FE=ZERO
            j=0
            DO i=1,size(unique)
                IF (f(i).GE. threshold) THEN
                    j=j+1
                    IF (j .LT. k) WHERE (x.EQ.unique(i)) FE(:,j)=ONE
                END IF
            END DO
        ELSE
            k=size(unique)
            ALLOCATE(FE(N,k-1))
            FE=ZERO
            DO i=1,k-1
                WHERE (x.EQ.unique(i)) FE(:,i)=ONE
            END DO
            j=k
            above_threshold=.TRUE.
        END IF

		IF (PRESENT(verbose) .AND. verbose) THEN
            PRINT *,' '
            PRINT *, k-1, ' fixed effects created.'
            PRINT *,' '
            PRINT *,'Here are the FE categories and the frequency at which they appear:'
            DO i=1,size(unique)
                IF (above_threshold(i)) THEN
                    PRINT *,unique(i),f(i)
                    j=i ! to keep track of the last category
                END IF
            END DO
            PRINT *,' '
            PRINT *, 'benchmark category (no FE created):'
            PRINT *, unique(j)
            PRINT *,' '
        END IF
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
END FUNCTION createFixedEffects

END MODULE linearRegression
