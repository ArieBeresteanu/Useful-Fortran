MODULE linAlgebra
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8

IMPLICIT NONE
! This module contains the REAL(dp) version of the code posted at:
! http://fortranwiki.org/fortran/show/Matrix+inversion
INTEGER, PARAMETER :: dp = kind(1.0d0)

CONTAINS

  PURE FUNCTION matinv2(A) RESULT(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    REAL(dp), INTENT(IN) :: A(2,2)   !! Matrix
    REAL(dp)             :: B(2,2)   !! Inverse matrix
    REAL(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0D0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  END FUNCTION

  PURE FUNCTION matinv3(A) RESULT(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    REAL(dp), INTENT(IN) :: A(3,3)   !! Matrix
    REAL(dp)             :: B(3,3)   !! Inverse matrix
    REAL(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0D0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  END FUNCTION

  PURE FUNCTION matinv4(A) RESULT(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    REAL(dp), INTENT(IN) :: A(4,4)   !! Matrix
    REAL(dp)             :: B(4,4)   !! Inverse matrix
    REAL(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
	1.0D0/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
	- A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
	+ A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
	- A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  END FUNCTION
  
  
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8 
  SUBROUTINE matinv(n,matIn,matOut)
  
	INTEGER, INTENT(IN)		:: n !matrix dimension
	REAL(dp), INTENT(IN)	:: matIN(n,n)  !input matrix to be inverted
	REAL(dp), INTENT(OUT)	:: matOut(n,n) !the inverse of the above
	
	SELECT CASE(n)
		CASE (1)
			matOut=1.0D0/matIn
		CASE (2)
			matOut=matinv2(matIn)
		CASE (3)
			matOut=matinv3(matIn)
		CASE (4)
			matOut=matinv4(matIn)
		CASE DEFAULT
			CALL matinvL(n,matIn,matOut)
	END SELECT
	
  END SUBROUTINE matinv

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8   

SUBROUTINE matinvL(n,a,minverse)
!matrix inverse, used if n>=5
INTEGER, INTENT(IN)					:: n
REAL(dp),INTENT(IN),DIMENSION(n,n)	:: a
REAL(dp),INTENT(OUT),DIMENSION(n,n)	:: minverse
!more variables
REAL(dp),DIMENSION(n,n)				:: y,b
INTEGER,DIMENSION(n)				:: indx
REAL(dp)							:: d
INTEGER								:: i

y=0.0_dp
DO i=1,n
	y(i,i)=1.0_dp
END DO

b=a

CALL ludcmp(b,indx,d)

do i=1,n
  CALL lubksb(b,indx,y(:,i))
END DO

minverse=y

END SUBROUTINE matinvL

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8   


SUBROUTINE ludcmp(a,indx,d)
! LU decomposition of matrix a 

REAL(dp),DIMENSION(:,:),INTENT(INOUT) 	:: a
INTEGER,DIMENSION(:),INTENT(OUT) 		:: indx
REAL(dp),INTENT(OUT) 					:: d
!additional variables
REAL(dp),DIMENSION(SIZE(a,1)) :: vv
REAL(dp),PARAMETER :: TINY=1.0e-20_dp
INTEGER:: j,n,imax

n=SIZE(indx)
d=1.0_dp
vv=MAXVAL(ABS(a),DIM=2)
IF (ANY(vv.LT.TINY)) THEN
   d=-99.0_dp !error
   RETURN
 END IF
   
vv=1.0_dp/vv
DO j=1,n
        imax=(j-1)+MAXLOC(vv(j:n)*ABS(a(j:n,j)), dim=1)
        IF (j /=imax) THEN
           CALL swap_rv(a(imax,:),a(j,:))
           d=-d
           vv(imax)=vv(j)
        END IF
        indx(j)=imax
        if (a(j,j).EQ.0.0_dp) a(j,j)=TINY
        a(j+1:n,j)=a(j+1:n,j)/a(j,j)
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
END DO
END SUBROUTINE ludcmp
 
!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8   
 SUBROUTINE lubksb(a,indx,b)
!solving linear system A*x=b when A is a matrix after a LU decomposition
!the solution is returned into b
REAL(dp),DIMENSION(:,:),INTENT(IN) 	:: a
INTEGER,DIMENSION(:),INTENT(IN) 	:: indx
REAL(dp),DIMENSION(:),INTENT(INOUT) :: b
!additional variables
INTEGER:: i,n,ii,ll
REAL(dp) :: summ

n=SIZE(indx)
ii=0
do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii /= 0)then
                 summ=summ-DOT_PRODUCT(a(i,ii:i-1),b(ii:i-1))
        else if (summ /= 0.0_dp) THEN
                 ii=i
        endif
        b(i)=summ
END DO
DO i=n,1,-1
        b(i)=(b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/a(i,i)
END DO

END SUBROUTINE lubksb

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8   
FUNCTION outerprod(a,b)
REAL(dp),DIMENSION(:),INTENT(IN)	::a,b
REAL(dp),DIMENSION(SIZE(a),SIZE(b))	::outerprod
outerprod=SPREAD(a,DIM=2,ncopies=SIZE(b))*&
          SPREAD(b,DIM=1,ncopies=SIZE(a))
END function outerprod

!--------1!--------2!--------3!--------4!--------5!--------6!--------7!--------8   
SUBROUTINE swap_rv(a,b)
REAL(dp),DIMENSION(:),INTENT(INOUT)	::a,b
REAL(dp),DIMENSION(SIZE(a))			::dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rv


END MODULE linAlgebra