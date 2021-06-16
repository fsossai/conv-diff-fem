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
character(100)          :: arg_coord, arg_topo

! Checking input arguments
if (iargc().lt.2) then
    stop 'Missing arguments. Please specify coordinates and topology files.'
end if

call getarg(1, arg_coord)
call getarg(2, arg_topo)

! Core computations
call cpu_time(clock_t_start)
timer = omp_get_wtime()

print '(a)', 'Reading coordinates file...'
call read_coord(arg_coord, coord)
print '(a)', 'Reading topology file...'
call read_topo(arg_topo, topo)

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
call write_vec('solution.txt', x)

print *
print '(a20,1en20.3)', 'Elapsed time (s):', timer
print '(a20,1en20.3)', 'Clock time (s):', clock_t_end - clock_t_start

deallocate(coord, topo, x0, x)

write(*,'(a)') 'Done'
end program main
