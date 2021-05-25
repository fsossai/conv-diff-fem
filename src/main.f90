program main

use class_CSRMAT
use class_utils
use class_bicgstab
use omp_lib

implicit none

! Input variables
logical                        :: BINREAD
type(CSRMAT)                   :: mat_A
character(len=100)             :: mat_name

! Local variables
integer                        :: nn,nt
integer                        :: i,ierr

! Local allocatable variables
real(dp), allocatable :: b(:), vec_y(:), x(:)

! Handles
integer, pointer               :: iat(:), ja(:)
real(dp), pointer     :: coef(:)
real(dp) :: timer

! Open the input file
open(1, file='inputs/test1.param', status='old')
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

! Allocate b and vec_y
allocate(b(nn),vec_y(nn),x(nn),stat=ierr)
if (ierr /= 0) stop 'Error in allocating b and vec_y'

! Set handles to mat_A
iat  => mat_A%patt%iat
ja   => mat_A%patt%ja
coef => mat_A%coef


! Core computations
print *, 'Matrix size:', nn
b = 1.0_dp
timer = omp_get_wtime()
call bicgstab(mat_A, b, x, tol=1e-6_dp)
timer = omp_get_wtime() - timer
!call print_vec_compact(x, 5)
call write_vec('solution.txt', x)
print '(a,1f10.6)', 'Elapsed time (s):', timer

! Deallocate the matrix
ierr = dlt_CSRMAT(mat_A)
if (ierr /= 0) stop 'Error in deallocating mat_A'

! Deallocate b and vec_y
deallocate(b,vec_y,stat=ierr)
if (ierr /= 0) stop 'Error in deallocating b and vec_y'

write(*,'(a)') 'Done'
end program main
