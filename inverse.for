      subroutine inverse(a,c,n)
c============================================================
c Inverse matrix
c Method: Based on Doolittle LU factorization for Ax=b
c Alex G. December 2009
c-----------------------------------------------------------
c input ...
c a(n,n) - array of coefficients for matrix A
c n      - dimension
c output ...
c c(n,n) - inverse matrix of A
c comments ...
c the original matrix a(n,n) will be destroyed 
c during the calculation
c===========================================================
      implicit none 
      integer n
c      double precision a(n,n), c(n,n)
c      double precision L(n,n), U(n,n), b(n), d(n), x(n)
c      double precision coeff
      real*8 a(n,n), c(n,n)
      real*8 L(n,n), U(n,n), b(n), d(n), x(n)
      real*8 coeff
      integer i, j, k

c step 0: initialization for matrices L and U and b
c Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

c step 1: forward elimination
      do k=1, n-1
        do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
          end do
        end do
      end do

c Step 2: prepare L and U matrices 
c L matrix is a matrix of the elimination coefficient
c + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
c U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do

c Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
c Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
c Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
c Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse
