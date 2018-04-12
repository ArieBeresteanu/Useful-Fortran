MODULE statistics

IMPLICIT NONE
SAVE

!some useful constants
INTEGER, PARAMETER :: dp = kind(1.0d0)
REAL(dp), PARAMETER :: HALF=0.50D0, epsilon=0.05D0, tol=1.0D-5, negONE=-1.0D0, ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0

INCLUDE


!********** wighted mean *******************
REAL(dp) FUNCTION mean( x,weights,mask) RETURN(mean)
! this procedure computes the eighted mean of the vector x
! if weights are not given, we assume equal weights of 1

REAL(dp), INTENT(IN), DIMENSION(:)           	:: x
REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:) 	:: weights
LOGICAL, OPTIONAL, INTENT(IN), DIMENSION(:) 	:: mask

REAL(dp)	:: denom,numer
INTEGER		:: N,i

N=size(x)
denom=ZERO
numer=ZERO
	
	
DO i=1,N
	IF (PRESENT(mask) .AND. (.NOT, mask)) CYCLE
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
REAL(dp) FUNCTION var(x,weights,mask) RESULT(var)
! this function computes the weighted variance of a vector x
! if weights are not given, we assume equal weights.
! The function can handle a mask variable.
! A typical call can look like:
! v=var(x,WEIGHTS=w,MASK=m)

REAL(dp), INTENT(IN),DIMENSION(:) :: x
REAL(dp), INTENT(IN), OPTIONAL, DIMENSION(:) :: weights
LOGICAL, OPTIONAL, INTENT(IN), DIMENSION(:) 	:: mask

REAL(dp)	:: Swx2, Swx, Sw
INTEGER		:: N,i

N=size(x)
Swx2=ZERO
Swx=ZERO
Sw=ZERO
	
IF (PRESENT(mask))) THEN
	DO i=1,N
		IF (.NOT, mask) CYCLE
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
REAL(dp) FUNCTION covariance (x,y,weights)
! this function computes the weighted covariance between x and y
! if weights are not given, we assume equal weights
! if x and y are not of the same length the covariance is set to -HUGE

REAL(dp), INTENT(IN),DIMENSION(:) :: x,y
REAL(dp), INTENT(IN), OPTIONAL, DIMENSION(:) :: weights

IF (SIZE(x)=SIZE(y)) THEN
	IF (PRESENT(weighted)) THEN
		covariance= mean(x*y,weights)-mean(x,weights)*mean(y,weights)
	ELSE
		covariance=mean(x**y) - mean(x)*mean(y)
	ENDIF
ELSE
	covariance=-HUGE(ONE)
ENDIF

REAL(dp) FUNCTION norm(x,p) RETURN(norm)
! this function computes the L_p norm of vector x
! if p is not given, we assume p=2

REAL(dp), INTENT(IN), DIMENSION(:) :: x
INTEGER, INTENT(IN), OPTIONAL :: p

IF (PRESENT(p)) THEN
	pr=REAL(p,dp)
ELSE
	pr=TWO
ENDIF

norm=(SUM(x**pr))**(ONE/pr)

END FUNCTION norm

!****************** sub-sample *******************************
SUBROUTINE subsample(original_list,m,with,sub_list)
! this subroutine subsamples m elements from the vector original_list and put them in sub_list
! if with=.TRUE. we sample with replacement and if with=.FALSE. we sample without replacement
! if we sample without replacement m needs to be smaller than size(original_list) or we get an error

REAL(dp), INTENT(IN), DIMENSION(:)  	:: original_list
INTEGER, INTENT(IN)                		:: m
LOGICAL, INTENT(IN)                		:: with
REAL(dp), INTENT(OUT), DIMENSION(m) 	:: sub_list

REAL(dp), DIMENSION(:)              	:: duplicate_list
INTEGER                            		:: n, j, rand_indx 
REAL(dp)                            	:: h

duplicate_list=original_list

IF (.NOT.with .AND. m.GT.n) THEN
	sub_list=-HUGE(ONE)
ELSE
	DO j=1,m 
		CALL RANDOM_NUMBER(h)
		rand_indx=INT(h*REAL(n-j+1,8)+1
		sub_list(j)=duplicate_list(rand_indx)
		!swapping:
		duplicate_list(rand_indx)=duplicate_list(n-j+1)-duplicate_list(rand_indx)
		duplicate_list(n-j+1)=duplicate_list(n-j+1)-duplicate_list(rand_indx) 
		duplicate_list(rand_indx)=duplicate_list(n-j+1)+duplicate_list(rand_indx) 
	END DO
ENDIF		

END SUBROUTINE subsample			


!******************** Gini Coefficient ***************************
 FUNCTION realGini(N,list,mask) RETURN(G)
!
! computes the Gini index for a list of reals
!
INTEGER, INTENT(IN)							:: N
REAL(dp), INTENT(IN), DIMENSION(N) 			:: list
LOGICAL, INTENT(IN), DIMENSION(N), OPTIONAL	:: mask
REAL(dp), INTNET(OUT)						:: G

REAL(dp)	:: denom,numer
INTEGER		:: i,j

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
G=numer/(2*N*denom)

END FUNCTION realGini


!******************** Frequency table ******************
SUBROUTINE freq(N,list,mask,unique,frequency)
!   Finds the unique elements of a list	of integers (unsorted) 
!   and returns them in another list (allocatable) as well as
!	a count of occurrence of each unique value (allocatable)
    INTEGER, INTENT(IN) 				:: N
	INTEGER, INTENT(IN)					:: list(N)
	INTEGER, INTENT(OUT), ALLOCATABLE	:: unique,frequency
	
    INTEGER 	:: temp_list(N), temp_count(N)
    INTEGER 	:: k      ! The number of unique elements
    INTEGER 	:: i

    k = 1
    temp_list(1) = list(1)
    DO i=2,N
        ! if the number already exist in res check next
        IF (ANY( temp_list(1:k) .EQ. list(i) )) 
			! a match was found
			WHERE (temp_list(1:k) .EQ. list(i)) temp_count=temp_count+1
		ELSE
			! No match found so add it to the output
			k = k + 1
			temp_list(k) = list(i)
		END IF
    END DO
	
	ALLOCATE(unique(k),frequency(k))
	unique=temp_list(1:k)
	frequency=temp_count(1:k)
	
END SUBROUTINE freq

FUNCTION entropy(N,freq) RETURN(shannon)
! compute the Shannon information entropy
! the input is a frequency vector (or probabilities)
INTEGER, INTENT(IN)		:: N
REAL(dp), INTENT(IN)	:: freq(N)
REAL(dp), INTENT(OUT)	:: shannon

REAL(dp)	:: temp(N)


temp=freq/SUM(feq)
temp=temp*LOG(temp)

shannon=-SUM(temp)/LOG(TWO)
END FUNCTION shannon

END MODULE statistics