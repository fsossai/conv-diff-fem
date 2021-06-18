program main
use BLAS
use class_precision
use omp_lib

include 'readers.h'

real(dp), parameter                 :: flops_per_term = 2
real(dp), parameter                 :: flops_per_row = 2
integer                             :: nn, nt, iterations
type(CSRMAT)                        :: A
real(dp), dimension(:), allocatable :: x, y, z
real(dp)                            :: timer, total_flops, delta
character(len=100)                  :: filename

read(*,'(a)') filename
iterations = 10000

! reading matrix
call read_mat(filename, .false., A, ierr)

nn = A%patt%nrows
nt = A%patt%nterm

! allocation and initialization
allocate(x(nn), y(nn), z(nn))
do i = 1,nn
    x = dble(i)
    y = nn - dble(i) + 1
end do

print '(a25,i15)', 'Iterations:',       iterations
print '(a25,i15)', 'Matrix size:',      nn
print '(a25,i15)', 'Non-zero terms:',   nt

timer = 0.0_dp
!$omp parallel private(i)
do i = 1, iterations
    !$omp master
    delta = omp_get_wtime()
    !$omp end master
    
    call amxpby_set(z, dble(i), A, x)
    
    !$omp master
    delta = omp_get_wtime() - delta
    timer = timer + delta
    !$omp end master
    
    !$omp barrier
end do
!$omp end parallel

total_flops = (flops_per_term * nt + flops_per_row * nn) * iterations
print '(a25,en15.3)', 'Elapsed time [s]:',  timer
print '(a25,en15.3)', 'Avg time [s]:',      timer / iterations
print '(a25,en15.3)', 'Performance [GFlops]:', total_flops / timer * 1e-9

ierr = dlt_CSRMAT(A)
deallocate(x, y, z)

end program