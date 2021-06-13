program main

use class_CSRMAT
use fem
use utils
use omp_lib
   
implicit none

include 'readers.h'

! Local variables
real(dp), allocatable   :: x(:), x0(:)
real(dp)                :: timer, clock_t_start, clock_t_end

! Topology and Nodal coordinates
integer                 :: nnodes, nelem
integer, allocatable    :: topo(:,:)
real(dp), allocatable   :: coord(:,:)

! Core computations
call cpu_time(clock_t_start)
timer = omp_get_wtime()

call read_coord('inputs/grid1.coord.txt', coord)
call read_topo('inputs/grid1.topo.txt', topo)
nnodes = size(coord, 1)
nelem = size(topo, 1)
print '(a20,i10)', 'Number of nodes:',      nnodes
print '(a20,i10)', 'Number of elements:',   nelem
allocate(x(nnodes), x0(nnodes))
x0 = 0.0_dp
call solve(coord, topo, x0, x)

timer = omp_get_wtime() - timer
call cpu_time(clock_t_end)

!call print_vec_compact(x, 5)
!call write_vec('solution.txt', x)

print *
print '(a20,1en20.3)', 'Elapsed time (s):', timer
print '(a20,1en20.3)', 'Clock time (s):', clock_t_end - clock_t_start

deallocate(coord, topo, x0, x)

write(*,'(a)') 'Done'
end program main
