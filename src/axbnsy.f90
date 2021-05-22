subroutine axbnsy(n,ntot,nterm,iat,ja,coef_A,xvec,bvec)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: axbnsy
!
!  Coded by Carlo Janna, Massimiliano Ferronato and Nicola Castelletto
!  December 2010
!
!  Purpose: compute bvec = [A]*xvec
!           
!  Variables:
!
!  n      : # of matrix A rows 
!  ntot   : # of xvec components
!  nterm  : # of matrix A non-zeroes 
!  iat    : integer array of the pointers to the beginning of each row of matrix A
!  ja     : integer array of the column indices for matrix A
!  coef_A : real(kind=dp) array of the matrix coefficients
!  xvec   : real(kind=dp) array
!  bvec   : matrix-vector product [A]*xvec
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer, intent(in)           :: n,ntot,nterm
real(kind=dp), intent(in)     :: coef_A(nterm),xvec(ntot)
integer, intent(in)           :: ja(nterm),iat(n+1)

! Ouput variables
real(kind=dp), intent(out) :: bvec(n)

! Local variables
integer i,k,m,mm

!-----------------------------------------------------------------------------------------

do k = 1,n
   m = iat(k)
   mm = iat(k+1)-1
   bvec(k) = coef_A(m)*xvec(ja(m))
   do i = m+1,mm
      bvec(k) = bvec(k) + coef_A(i)*xvec(ja(i))
   end do
end do

!-----------------------------------------------------------------------------------------

end subroutine axbnsy
