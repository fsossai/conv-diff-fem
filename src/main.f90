program main

use class_CSRMAT

implicit none

! Interfaces
include 'axbnsy_int.h'
include 'readmat_int.h'

! Input variables
logical                        :: BINREAD
type(CSRMAT)                   :: mat_A
character(len=100)             :: mat_name, temp

! Local variables
integer                        :: nn,nt
integer                        :: i,ierr

! Local allocatable variables
! Local allocatable variables
real(kind=double), allocatable :: vec_x(:),vec_y(:)

! Handles
integer, pointer               :: iat(:),ja(:)
real(kind=double), pointer     :: coef(:)

! Open the input file
open(1,file='inputs/test1.param',status='old')
read(1,*) BINREAD

read(1,*) mat_name
close(1)

! Open and read the matrix file
if (BINREAD) then
   ! Open in binary format and read matrix header
   open(11,file=mat_name,status='old',form='unformatted',access='stream')
   read(11) nn,nt
else
   ! Open as ascii and read matrix header
   open(11,file=mat_name,status='old')
   read(11,*) nn,nt
endif

! Allocate the matrix
ierr = new_CSRMAT(nn,nt,mat_A)
if (ierr /= 0) stop 'Error in allocating mat_A'

! Read the matrix 
call readmat(11,BINREAD,mat_A,ierr)
if (ierr /= 0) stop 'Error in reading mat_A'

! Allocate vec_x and vec_y
allocate(vec_x(nn),vec_y(nn),stat=ierr)
if (ierr /= 0) stop 'Error in allocating vec_x and vec_y'

! Set all ones in vec_x
vec_x = 1._double

! Set handles to mat_A
iat  => mat_A%patt%iat
ja   => mat_A%patt%ja
coef => mat_A%coef

! Compute matrix by vector product
call axbnsy(nn,nn,nt,iat,ja,coef,vec_x,vec_y)

! Print the result
open(10,file='vec_y.txt')
write(10,'(i10,e15.6)') (i,vec_y(i),i=1,nn)
close(10)

! Deallocate the matrix
ierr = dlt_CSRMAT(mat_A)
if (ierr /= 0) stop 'Error in deallocating mat_A'

! Deallocate vec_x and vec_y
deallocate(vec_x,vec_y,stat=ierr)
if (ierr /= 0) stop 'Error in deallocating vec_x and vec_y'

write(*,'(a)') 'Done'
end program main
