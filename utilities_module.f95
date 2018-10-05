MODULE utilities
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

USE basics
IMPLICIT NONE

CONTAINS

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION LP2(vec) RESULT(dist)

REAL(dp), INTENT(IN)	:: vec(:)
REAL(dp)							:: dist

dist=SQRT(SUM(vec**2))

END FUNCTION LP2

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION LP1(vec) RESULT(dist)
REAL(dp), INTENT(IN)		:: vec(:)
REAL(dp)					:: dist

dist=SUM(ABS(vec))

END FUNCTION LP1

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION LPinf(vec) RESULT(dist)
REAL(dp), INTENT(IN)		:: vec(:)
REAL(dp)					:: dist

dist=MAXVAL(ABS(vec))

END FUNCTION LPinf

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
SUBROUTINE descriptive(Nobs,x,mu,sigma)
INTEGER(regInt), INTENT(IN)				:: Nobs
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

INTEGER(regInt)	:: i,m

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

INTEGER(regInt)	:: i,m

m=SIZE(vec) !number of elements
DO i=1,m
	IF (PRESENT(f)) THEN
		PRINT f, vec(i)
	ELSE
		PRINT *,vec(i)
	END IF
END DO

END SUBROUTINE printRealVector

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION intToString(i) result(res)
  CHARACTER(LEN=256)            :: res
  INTEGER(regInt),INTENT(IN)    :: i
  CHARACTER(RANGE(i)+2)         :: tmp

  WRITE(tmp,'(i0)') i
  res = TRIM(tmp)
END FUNCTION intToString

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION to_lower (str) RESULT (string)
!   ==============================
!   Changes a string to lower case
!   ==============================

    CHARACTER(*), INTENT(IN) :: str
    CHARACTER(LEN(str))      :: string

    INTEGER(regInt) :: ic, i

    CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'

!   lowercase each letter IF it is capitalized
    string = str
    DO i = 1, LEN_TRIM(str)
        ic = INDEX(cap, str(i:i))
        IF (ic .GT. 0) string(i:i) = low(ic:ic)
    END DO

END FUNCTION to_lower


!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8
FUNCTION remove_dups(N,list) RESULT(unique)
!   ====================================
!   Finds the unique elements of a list
!	of INTEGER(regInt)s (unsorted) and returns
!	them in another list (allocatable)
!	Frequencies are not computed!
!   ====================================

    INTEGER(regInt), INTENT(IN) 				:: N
	INTEGER(regInt), INTENT(IN)					:: list(N)
	INTEGER(regInt), DIMENSION(:), ALLOCATABLE	:: unique

    INTEGER(regInt) 	:: temp(N)
    INTEGER(regInt) 	:: k      ! The number of unique elements
    INTEGER(regInt) 	:: i

    k = 1
    temp(1) = list(1)
    DO i=2,N
        ! IF the number already exist in res check next
        IF (ANY( temp(1:k) .EQ. list(i) )) CYCLE
        ! No match found so add it to the output
        k = k + 1
        temp(k) = list(i)
    END DO

	ALLOCATE(unique(k))
	unique=temp(1:k)

END FUNCTION remove_dups

SUBROUTINE hms_current()

!*****************************************************************************80
! based on: https://people.sc.fsu.edu/~jburkardt/f_src/timestamp/timestamp.f90
!

  CHARACTER ( LEN = 2 )     :: ampm
  INTEGER(regInt)      :: h
  INTEGER(regInt)      :: mm
  INTEGER(regInt)      :: n
  INTEGER(regInt)      :: s
  INTEGER(regInt)      :: values(8)

  CALL date_and_time ( values = values )

  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  IF ( h < 12 ) THEN
    ampm = 'AM'
  ELSE IF ( h == 12 ) THEN
    IF ( n == 0 .and. s == 0 ) THEN
      ampm = 'Nn'
    ELSE
      ampm = 'PM'
    END IF
  ELSE
    h = h - 12
    IF ( h < 12 ) THEN
      ampm = 'PM'
    ELSE IF ( h == 12 ) THEN
      IF ( n == 0 .and. s == 0 ) THEN
        ampm = 'Md'
      ELSE
        ampm = 'AM'
      END IF
    END IF
  END IF

  PRINT '(i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)', &
    h, ':', n, ':', s, '.', mm, trim ( ampm )

  RETURN
END SUBROUTINE hms_current

SUBROUTINE array2json(arr,varNames,jsonOutFile)
    REAL(dp), DIMENSION(:,:),INTENT(IN)         :: arr
    CHARACTER(len=*), DIMENSION(:), INTENT(IN)  :: varNames
    CHARACTER(len=*), INTENT(IN)                :: jsonOutFile

    INTEGER(regINT)     :: N,k,i,j,ios

    N=SIZE(arr,1)
    k=SIZE(arr,2)
    IF (SIZE(varNames) .NE. k) THEN
        PRINT *,"dimension mismatch"
        RETURN
    END IF

    OPEN (Unit=88,file=jsonOutFile,IOSTAT=ios)
    IF (ios.ne.0) THEN
        PRINT *,'error opening file'
        RETURN
    END IF

    WRITE(unit=88,FMT=*) "["
    DO i=1,N
        WRITE(unit=88,FMT=*) " {"
        DO j=1,k-1
            WRITE(unit=88,FMT=*) " """,varNames(j),""":",arr(i,j),","
        END DO
        WRITE(unit=88,FMT=*) " """,varNames(k),""":",arr(i,k)

        IF (i.NE. N) THEN
            WRITE(unit=88,FMT=*) "},"
        ELSE
            WRITE(unit=88,FMT=*) "}"
        END IF
    END DO
    WRITE(unit=88,FMT=*) "]"
    CLOSE(unit=88)
    PRINT *,"data saved at ",jsonOutFile

END SUBROUTINE array2json

END MODULE utilities
