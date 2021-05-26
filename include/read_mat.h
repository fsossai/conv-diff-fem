interface readmat_int
   subroutine readmat(iunit,BINREAD,mat,info)

   use class_CSRMAT

   implicit none

   ! Input variables
   integer, intent(in)         :: iunit
   logical, intent(in)         :: BINREAD

   ! Input/Output variables
   type(CSRMAT), intent(inout) :: mat

   ! Output variables
   integer, intent(out)        :: info

   end subroutine readmat
end interface readmat_int
