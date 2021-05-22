interface axbnsy_int
   subroutine axbnsy(n,ntot,nterm,iat,ja,coef_A,xvec,bvec)

   use class_precision

   implicit none

   ! Input variables
   integer, intent(in)           :: n,ntot,nterm
   real(kind=dp), intent(in)     :: coef_A(nterm),xvec(ntot)
   integer, intent(in)           :: ja(nterm),iat(n+1)
   ! Ouput variables
   real(kind=dp), intent(out) :: bvec(n)

   end subroutine axbnsy
end interface axbnsy_int
