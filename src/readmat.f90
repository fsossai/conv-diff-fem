subroutine readmat(iunit,BINREAD,mat,info)
!*****************************************************************************************
!
!  Subroutine: readmat
!
!  Coded by Carlo Janna, Massimiliano Ferronato and Nicola Castelletto
!  December 2010
!
!  Purpouse: read a matrix in coordinate format and store it in CSR type format
!
!  Variables
!
!  iunit    : Input unit
!  mat      : matrix
!  info     : error code 0 ---> successful run
!                        1 ---> error in allocation/deallocation of local scratch
!                        2 ---> end of file encountered
!
!*****************************************************************************************

use class_CSRMAT

implicit none

! Input variables
logical, intent(in)         :: BINREAD
integer, intent(in)         :: iunit

! Input/Output variables
type(CSRMAT), intent(inout) :: mat

! Output variables
integer, intent(out)        :: info

! Local variables
integer                     :: i,alinfo,readinfo
integer, allocatable        :: irow(:)

!-----------------------------------------------------------------------------------------

info = 0

! Allocate the local scratch
allocate(irow(mat%patt%nterm),stat=alinfo)
if (alinfo .ne. 0) then
   info = 1
   return
endif

if (BINREAD) then

   ! Read matrix in binary format
   do i = 1,mat%patt%nterm
      read(iunit,iostat=readinfo) irow(i),mat%patt%ja(i),mat%coef(i)
      if (readinfo .ne. 0) then
         info = 2
         deallocate(irow)
         return
      endif
   enddo

else

   ! Read matrix in coordinate format
   do i = 1,mat%patt%nterm
      read(iunit,*,iostat=readinfo) irow(i),mat%patt%ja(i),mat%coef(i)
      if (readinfo .ne. 0) then
         info = 2
         deallocate(irow)
         return
      endif
   enddo

endif

! Build topology vector mat%iat from the row indices list
call irow2iat(mat%patt%nrows,mat%patt%nterm,irow,mat%patt%iat)

! Deallocate the local scratch
deallocate(irow,stat=alinfo)
if (alinfo .ne. 0) then
   info = 1
endif

!-----------------------------------------------------------------------------------------

end subroutine readmat
