MODULE statistics

USE quickSort

IMPLICIT NONE

TYPE descriptiveStatistics
	REAL(dp)			:: mean,var,skewness,kurtosis,std
	REAL(dp)			:: minimum,maximum
	REAL(dp)			:: q01,q10,q25,median,q75,q90,q99
	REAL(dp)            :: silverman1,silverman2
	INTEGER(regInt)		:: N
END TYPE
!some useful constants
REAL(dp), PARAMETER, PRIVATE :: epsilon=0.05D0, tol=1.0D-5

CONTAINS


!********** wighted mean *******************
REAL(dp) FUNCTION mean(x,weight,mask)
! this procedure computes the eighted mean of the vector x
! if weights are not given, we assume equal weights of 1

REAL(dp), INTENT(IN), DIMENSION(:)           	:: x
REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:) 	:: weight
LOGICAL, OPTIONAL, INTENT(IN), DIMENSION(:) 	:: mask

REAL(dp)			:: denom,numer
INTEGER(regInt)		:: N,i

N=size(x)
denom=ZERO
numer=ZERO


DO i=1,N
	IF (PRESENT(mask) .AND. (.NOT. mask(i))) CYCLE
	IF (PRESENT(weight)) THEN
		numer=numer+weight(i)
		denom=denom+x(i)*weight(i)
	ELSE
		numer=numer+ONE
		denom=denom+x(i)
	END IF
END DO

mean=numer/denom

END FUNCTION mean




!*********** weighted variance *********************
REAL(dp) FUNCTION var(x,weight,mask)
! this function computes the weighted variance of a vector x
! if weights are not given, we assume equal weights.
! The function can handle a mask variable.
! A typical call can look like:
! v=var(x,WEIGHTS=w,MASK=m)

REAL(dp), INTENT(IN),DIMENSION(:) 				:: x
REAL(dp), INTENT(IN), OPTIONAL, DIMENSION(:) 	:: weight
LOGICAL, OPTIONAL, INTENT(IN), DIMENSION(:) 	:: mask

REAL(dp)			:: Swx2, Swx, Sw
INTEGER(regInt)		:: N,i

N=size(x)
Swx2=ZERO
Swx=ZERO
Sw=ZERO

IF (PRESENT(mask)) THEN
	DO i=1,N
		IF (.NOT. mask(i)) CYCLE
		IF (PRESENT(weight)) THEN
			Swx2=Swx2+weight(i)*x(i)*x(i)
			Swx=Swx+weight(i)*x(i)
			Sw=Sw+weight(i)
		ELSE
			Swx2=Swx2+x(i)*x(i)
			Swx=Swx+x(i)
			Sw=Sw+ONE
		END IF
	END DO
ELSE
	DO i=1,N
		IF (PRESENT(weight)) THEN
			Swx2=Swx2+weight(i)*x(i)*x(i)
			Swx=Swx+weight(i)*x(i)
			Sw=Sw+weight(i)
		ELSE
			Swx2=Swx2+x(i)*x(i)
			Swx=Swx+x(i)
			Sw=Sw+ONE
		END IF
	END DO
END IF

var=(Swx2-Swx)/Sw

END FUNCTION var



!*************** covariance ***********************
REAL(dp) FUNCTION cov (x,y,weight)
! this function computes the weighted covariance between x and y
! if weights are not given, we assume equal weights
! if x and y are not of the same length the covariance is set to -HUGE

REAL(dp), INTENT(IN),DIMENSION(:) 				:: x,y
REAL(dp), INTENT(IN), OPTIONAL, DIMENSION(:) 	:: weight

IF (SIZE(x) .EQ. SIZE(y)) THEN
	IF (PRESENT(weight)) THEN
		cov= mean(x*y,weight)-mean(x,weight)*mean(y,weight)
	ELSE
		cov=mean(x**y) - mean(x)*mean(y)
	ENDIF
ELSE
	cov=-HUGE(ONE)
ENDIF
END FUNCTION cov


REAL(dp) FUNCTION norm(x,p)
! this function computes the L_p norm of vector x
! if p is not given, we assume p=2

REAL(dp), INTENT(IN), DIMENSION(:) 		:: x
INTEGER(regInt), INTENT(IN), OPTIONAL 	:: p

REAL(dp)	:: pr

IF (PRESENT(p)) THEN
	pr=REAL(p,dp)
ELSE
	pr=TWO
ENDIF

norm=(SUM(x**pr))**(ONE/pr)

END FUNCTION norm


!**************** descriptive statistics ********************
SUBROUTINE descriptiveStats(x_in,mask,verbose,descriptive)
! this subroutine returns a set of descriptive statistics for the univariate
! variable x. A mask variable is optional.

REAL(dp), INTENT(IN), DIMENSION(:)           	:: x_in
LOGICAL, OPTIONAL, INTENT(IN), DIMENSION(:) 	:: mask
LOGICAL, OPTIONAL								:: verbose
TYPE(descriptiveStatistics), INTENT(OUT)		:: descriptive

INTEGER(regInt)		:: i,N
REAL(dp)            :: x(SIZE(x_in))
REAL(dp)			:: s1,s2,s3,s4,rN
REAL(dp)            :: s1sq,s1cu,rN2,rN3

N=size(x)
IF (N .LT. 2) RETURN
rN=REAL(N,dp)
rN2=rN*rN
rN3=rN2*rN
s1=ZERO
s2=ZERO
s3=ZERO
s4=ZERO

CALL quick_sort(x_in, x)
descriptive%q01=x(INT(N/100))
descriptive%q10=x(INT(N/10))
descriptive%q25=x(INT(N/4))
descriptive%median=x(INT(N/2))
descriptive%q75=x(INT(3*N/4))
descriptive%q90=x(INT(9*N/10))
descriptive%q99=x(INT(99*N/100))

IF (PRESENT(mask)) THEN
	descriptive%minimum=MINVAL(x,mask)
	descriptive%maximum=MAXVAL(x,mask)
	DO i=1,N
		IF (.NOT. mask(i)) CYCLE
		s1=s1+x(i)
		s2=s2+x(i)*x(i)
		s3=s3+x(i)*x(i)*x(i)
		s4=s4+x(i)*x(i)*x(i)*x(i)
	END DO
ELSE
	descriptive%minimum=x(1)
	descriptive%maximum=x(N)
	DO i=1,N
		s1=s1+x(i)
		s2=s2+x(i)*x(i)
		s3=s3+x(i)*x(i)*x(i)
		s4=s4+x(i)*x(i)*x(i)*x(i)
	END DO
END IF
s1sq=s1*s1
s1cu=s1sq*s1
rN2=rN*rN
rN3=rN2*rN
descriptive%mean=s1/rN
descriptive%var=(s2-s1sq/rn)/(rN-ONE)
descriptive%std=SQRT(descriptive%var)
descriptive%skewness=((rN*SQRT(rN-ONE))/(rN-TWO)) * &
                    (s3-3.0*s2*descriptive%mean+3.0*s1*(descriptive%mean)**2 - rN*descriptive%mean**3) / &
                    ((SQRT(s2-rN*descriptive%mean**2))**3)
descriptive%kurtosis=( (rN*(rN+ONE))*(s4-4.0*s3*s1/rN+6.0*s2*s1sq/rN2-3.0*s1*s1cu/rN3)/ &
                            ((rN-ONE)*descriptive%var**2) -3.0*(rN-ONE)**2 ) / &
                            ((rN-TWO)*(rN-3.0))
descriptive%N=N
descriptive%silverman1=1.06_dp*descriptive%std/(rN**0.2)
descriptive%silverman2=0.90_dp*MIN(descriptive%q75-descriptive%q25,descriptive%std)/(rN**0.2)

500 FORMAT (A10,F15.6)
IF (verbose) THEN
	PRINT *,"--------------------------------------"
	PRINT *, "Number of observations: ", descriptive%N
	PRINT *, " "
	PRINT *, "Moments:"
	PRINT *, "--------"
	PRINT 500, "Mean:     ", descriptive%mean
	PRINT 500, "Variance: ", descriptive%var
	PRINT 500, "Stdev:    " , descriptive%std
	PRINT 500, "Skewness: ",descriptive%skewness
	PRINT 500, "Kurtosis: ", descriptive%kurtosis
	PRINT *, " "
	PRINT *, "Location:"
    PRINT *, "---------"
	PRINT 500, "Minimum: ", descriptive%minimum
	PRINT 500, "Q01:     ", descriptive%q01
	PRINT 500, "Q10:     ", descriptive%q10
	PRINT 500, "Q25:     ", descriptive%q25
	PRINT 500, "Median:  ", descriptive%median
	PRINT 500, "Q75:     ", descriptive%q75
	PRINT 500, "Q90:     ", descriptive%q90
	PRINT 500, "Q99:     ", descriptive%q99
	PRINT 500, "Maximum: ", descriptive%maximum
	PRINT *, " "
	PRINT *, "Silverman rules of thumb:"
    PRINT *, "-------------------------"
	PRINT 500, "version 1:", descriptive%silverman1
	PRINT 500, "version 2:", descriptive%silverman2
	PRINT *,"--------------------------------------"
END IF

END SUBROUTINE descriptiveStats

!****************** sub-sample *******************************
SUBROUTINE subsample(original_list,m,with,sub_list)
! this subroutine subsamples m elements from the vector original_list and put them in sub_list
! if with=.TRUE. we sample with replacement and if with=.FALSE. we sample without replacement
! if we sample without replacement m needs to be smaller than size(original_list) or we get an error

REAL(dp), INTENT(IN), DIMENSION(:)  	:: original_list
INTEGER(regInt), INTENT(IN)            	:: m
LOGICAL, INTENT(IN)                		:: with
REAL(dp), INTENT(OUT), DIMENSION(m) 	:: sub_list

REAL(dp), DIMENSION(SIZE(original_list))   	:: duplicate_list
INTEGER(regInt)                        		:: n, j, rand_indx
REAL(dp)                            		:: h

duplicate_list=original_list

IF ((.NOT. with) .AND. m.GT.n) THEN
	sub_list=-HUGE(ONE)
ELSE
	DO j=1,m
		CALL RANDOM_NUMBER(h)
		rand_indx=INT(h*REAL(n-j+1,dp))+1
		sub_list(j)=duplicate_list(rand_indx)
		!swapping:
		duplicate_list(rand_indx)=duplicate_list(n-j+1)-duplicate_list(rand_indx)
		duplicate_list(n-j+1)=duplicate_list(n-j+1)-duplicate_list(rand_indx)
		duplicate_list(rand_indx)=duplicate_list(n-j+1)+duplicate_list(rand_indx)
	END DO
ENDIF

END SUBROUTINE subsample


!******************** Gini Coefficient ***************************
REAL(dp)  FUNCTION realGini(N,list,mask)
!
! computes the Gini index for a list of reals
!
INTEGER(regInt), INTENT(IN)					:: N
REAL(dp), INTENT(IN), DIMENSION(N) 			:: list
LOGICAL, INTENT(IN), DIMENSION(N), OPTIONAL	:: mask

REAL(dp)	:: denom,numer
INTEGER(regInt)		:: i,j

denom=ZERO
numer=ZERO
IF (PRESENT(mask)) THEN
	DO i=1,N
		IF (.NOT. mask(i)) CYCLE
		denom=denom+list(i)
		DO j=1,N
			IF (.NOT. mask(j)) CYCLE
			numer=numer+ABS(list(i)-list(j))
		END DO
	END DO
ELSE
	DO i=1,N
		denom=denom+list(i)
		DO j=1,N
			numer=numer+ABS(list(i)-list(j))
		END DO
	END DO
END IF
realGini=numer/(2*N*denom)

END FUNCTION realGini


!******************** Frequency table ******************
SUBROUTINE freq(N,list,unique,frequency,verbose)
!   Finds the unique elements of a list	of integers (unsorted)
!   and returns them in another list (allocatable) as well as
!	a count of occurrence of each unique value (allocatable)
    INTEGER(regInt), INTENT(IN) 				:: N
	INTEGER(regInt), INTENT(IN)					:: list(N)
	INTEGER(regInt), INTENT(OUT), ALLOCATABLE	:: unique(:),frequency(:)
	LOGICAL, INTENT(IN), OPTIONAL               :: verbose

	INTEGER(regInt) 	:: temp_list(N), temp_count(N)
    INTEGER(regInt) 	:: k      ! The number of unique elements
    INTEGER(regInt) 	:: i,j,item
	LOGICAL				:: found

    k = 1
    temp_list(1) = list(1)
	temp_count(1)=1
    DO i=2,N
		item=list(i)
		found=.FALSE.
		DO j=1,k
		    ! if the number already exist in res check next
			IF (item .EQ. temp_list(j)) THEN
				temp_count(j)=temp_count(j)+1
				found=.TRUE.
				EXIT
			END IF
		END DO
		IF (.NOT.found) THEN
			k=k+1
			temp_list(k)=item
			temp_count(k)=1
		END IF
    END DO

	ALLOCATE(unique(k),frequency(k))
	unique=temp_list(1:k)
	frequency=temp_count(1:k)
	IF (PRESENT(verbose) .AND. verbose) THEN
        DO i=1,k
            PRINT *,unique(i),frequency(i)
        END DO
    END IF

END SUBROUTINE freq


REAL(dp) FUNCTION entropy(N,freq)
! compute the Shannon information entropy
! the input is a frequency vector (or probabilities)
INTEGER(regInt), INTENT(IN)		:: N
REAL(dp), INTENT(IN)	:: freq(N)

REAL(dp)	:: temp(N)


temp=freq/SUM(freq)
temp=temp*LOG(temp)

entropy=-SUM(temp)/LOG(TWO)
END FUNCTION entropy

END MODULE statistics
