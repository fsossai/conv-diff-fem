program main
use class_precision

use class_CSRMAT
use omp_lib

include 'bicgstab.h'

integer                             :: nn, nt, nthreads
type(CSRMAT)                        :: A
real(dp), dimension(:), allocatable :: x, b
real(dp)                            :: timer
character(len=100)                  :: filename

read(*,'(a)') filename

! reading matrix
open(1, file=filename, status='old')
read(1, *) nn, nt
ierr = new_CSRMAT(nn, nt, A)
call readmat(1, .false., A, ierr)

! allocation and initialization
allocate(x(nn), b(nn))
b = 1.0_dp
x = 0.0_dp

print '(a25,i15)', 'Matrix size:',      nn
print '(a25,i15)', 'Non-zero terms:',   nt

timer = omp_get_wtime()
call bicgstab(A, b, x, tol=1e-6_dp)
timer = omp_get_wtime() - timer

print '(a25,i15)',    'Number of threads:', nthreads
print '(a25,en15.3)', 'Elapsed time [s]:',  timer

ierr = dlt_CSRMAT(A)
deallocate(x, b)

end program