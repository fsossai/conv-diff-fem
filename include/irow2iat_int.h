interface irow2iat_int
   subroutine irow2iat(n,nterm,irow,iat)
 
   implicit none

   !Input variables
   integer, intent(in) :: n,nterm
   integer, intent(in) :: irow(nterm)

   !Output variables
   integer, intent(out) :: iat(n+1)

   end subroutine irow2iat
end interface irow2iat_int
