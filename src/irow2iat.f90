subroutine irow2iat(n,nterm,irow,iat)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: irow2iat
!
!  Coded by Carlo Janna
!  April 2012
!
!  Purpose: Build the topology vector iat from a row indices list irow
!
!  Variables:
!
!  n     : # of rows of a given matrix stored in CSR format
!  nterm : # of non-zeroes of a given matrix stored in CSR format
!  irow  : row indices of non-zeroes of a given matrix stored in CSR format
!  iat   : integer array of the pointers to the beginning of each row
!
!-----------------------------------------------------------------------------------------
 
implicit none

!Input variables
integer, intent(in)  :: n,nterm
integer, intent(in)  :: irow(nterm)

!Output variables
integer, intent(out) :: iat(n+1)
 
!Local variables
integer              :: j,k,irow_old,irow_new

!-----------------------------------------------------------------------------------------

irow_old = 0
do k=1,nterm
   irow_new = irow(k)
   if ( irow_new .gt. irow_old ) then
      do j = irow_old+1,irow_new
         iat(j)=k
      enddo
      irow_old = irow_new
   end if
end do
k = nterm+1
do j = irow_old+1,n+1
   iat(j) = k
enddo
 
!-----------------------------------------------------------------------------------------

end subroutine irow2iat
