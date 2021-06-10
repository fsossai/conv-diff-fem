subroutine read_mat(filename, binary, mat, info)
!*****************************************************************************************
!
!  Subroutine: readmat
!
!  Coded by Carlo Janna, Massimiliano Ferronato and Nicola Castelletto
!  December 2010
!
!  Read a matrix in coordinate format and store it in CSR type format.
!
!  Variables
!
!  mat      : matrix
!  info     : error code 0 ---> successful run
!                        1 ---> error in allocation/deallocation of local scratch
!                        2 ---> end of file encountered
!
!*****************************************************************************************

use class_CSRMAT

implicit none

character(len=*), intent(in)   :: filename
logical, intent(in)              :: binary
type(CSRMAT), intent(inout) :: mat
integer, intent(out)        :: info

! Local variables
integer                     :: i, iunit, ierr, alinfo, readinfo, nn, nt
integer, allocatable        :: irow(:)

info = 0
iunit = 11

! Open and read the matrix file
if (binary) then
   ! Open in binary format and read matrix header
   open(iunit, file=filename, status='old', form='unformatted', access='stream')
   read(iunit) nn, nt
else
   ! Open as ascii and read matrix header
   open(iunit, file=filename, status='old')
   read(iunit, *) nn, nt
endif

! Allocate the matrix
ierr = new_CSRMAT(nn, nt, mat)
if (ierr /= 0) stop 'Error in allocating mat'

! Allocate the local scratch
allocate(irow(mat%patt%nterm), stat=alinfo)
if (alinfo .ne. 0) then
   info = 1
   return
endif

! Read the matrix 

if (binary) then
   ! Read matrix in binary format
   do i = 1,mat%patt%nterm
      read(iunit, iostat=readinfo) irow(i), mat%patt%ja(i), mat%coef(i)
      if (readinfo .ne. 0) then
         info = 2
         deallocate(irow)
         return
      endif
   enddo
else
   ! Read matrix in coordinate format
   do i = 1,mat%patt%nterm
      read(iunit, *, iostat=readinfo) irow(i), mat%patt%ja(i), mat%coef(i)
      if (readinfo .ne. 0) then
         info = 2
         deallocate(irow)
         return
      endif
   enddo
endif

! Build topology vector mat%iat from the row indices list
call irow2iat(mat%patt%nrows, mat%patt%nterm, irow, mat%patt%iat)

! Deallocate the local scratch
deallocate(irow, stat=alinfo)
if (alinfo .ne. 0) then
   info = 1
endif

end subroutine read_mat
