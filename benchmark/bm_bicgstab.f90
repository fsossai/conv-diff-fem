program main

use class_precision
use class_CSRMAT
use omp_lib

include 'bicgstab.h'
include 'readers.h'

integer                             :: nn, nt, max_it
type(CSRMAT)                        :: A
real(dp), dimension(:), allocatable :: x, b
real(dp)                            :: timer
character(len=100)                  :: filename, temp

! Checking input arguments
if (iargc().lt.2) then
    stop 'Missing arguments. Please specify the matrix filename and max iterations.'
end if

! Getting command line arguments
call getarg(1, filename)
call getarg(2, temp)
read(temp, '(i100)') max_it
print '(a25,i15)', 'Max iterations:',   max_it

! Reading sparse matrix from file
print '(a)', 'Reading matrix from file...'
call read_mat(filename, .false., A, ierr)

nn = A%patt%nrows
nt = A%patt%nterm

! Allocation and initialization
allocate(x(nn), b(nn))
b = 1.0_dp
x = 0.0_dp

print '(a25,i15)', 'Matrix size:',      nn
print '(a25,i15)', 'Non-zero terms:',   nt

timer = omp_get_wtime()
call bicgstab(A, b, x, tol=1e-6_dp, max_it=max_it)
timer = omp_get_wtime() - timer

print '(a25,en15.3)', 'Elapsed time [s]:',  timer

ierr = dlt_CSRMAT(A)
deallocate(x, b)

end program