program main
use BLAS
use class_precision
use omp_lib

real(dp), parameter                 :: flops_per_term = 2
real(dp), parameter                 :: flops_per_row = 2
integer                             :: nn, nt, iterations, nthreads
type(CSRMAT)                        :: A
real(dp), dimension(:), allocatable :: x, y, z
real(dp)                            :: timer, total_flops

read(*,*) iterations

! reading matrix
open(1, file='../inputs/mat1500.txt', status='old')
read(1, *) nn, nt
ierr = new_CSRMAT(nn, nt, A)
call readmat(1, .false., A, ierr)

! allocation and initialization
allocate(x(nn), y(nn), z(nn))
do i = 1,nn
    x = dble(i)
    y = nn - dble(i) + 1
end do

print '(a25,i15)', 'Iterations:',       iterations
print '(a25,i15)', 'Matrix size:',      nn
print '(a25,i15)', 'Non-zero terms:',   nt

timer = omp_get_wtime()
!$omp parallel private(i)
!$omp single
nthreads = omp_get_num_threads()
!$omp end single
do i = 1, iterations
    call amxpby(z, dble(i), A, x)
end do
!$omp end parallel
timer = omp_get_wtime() - timer

total_flops = (flops_per_term * nt + flops_per_row * nn) * iterations
print '(a25,i15)',    'Number of threads:', nthreads
print '(a25,en15.3)', 'Elapsed time [s]:',  timer
print '(a25,en15.3)', 'Avg time [s]:',      timer / iterations
print '(a25,en15.3)', 'Performance [GFlops]:', total_flops / timer * 1e-9

ierr = dlt_CSRMAT(A)
deallocate(x, y, z)

end program