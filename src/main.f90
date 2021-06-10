program main

use class_CSRMAT
use utils
use omp_lib
   
implicit none

include 'bicgstab.h'
include 'readers.h'

! Input variables
logical                        :: binary
type(CSRMAT)                   :: mat_A
character(len=100)             :: mat_name

! Local variables
integer                        :: nn
integer                        :: i, ierr

! Local allocatable variables
real(dp), allocatable :: b(:), vec_y(:), x(:)

! Handles
integer, pointer        :: iat(:), ja(:)
real(dp), pointer       :: coef(:)
real(dp)                :: timer, clock_t_start, clock_t_end

! Topology and Nodal coordinates
integer, allocatable        :: topo(:,:)
integer, allocatable        :: bnodes(:)
real(dp), allocatable       :: coord(:,:)

! Open the input file
open(1, file='inputs/test1.param', status='old')
read(1,*) binary
read(1,*) mat_name
close(1)

call read_mat(mat_name, binary, mat_A, ierr)
nn = mat_A%patt%nrows

! Allocate b and vec_y
allocate(b(nn),vec_y(nn),x(nn),stat=ierr)
if (ierr .ne. 0) stop 'Error in allocating b and vec_y'

! Set handles to mat_A
iat  => mat_A%patt%iat
ja   => mat_A%patt%ja
coef => mat_A%coef


! Core computations
print *, 'Matrix size:', nn
b = 1.0_dp
call cpu_time(clock_t_start)
timer = omp_get_wtime()
!call bicgstab(mat_A, b, x, tol=1e-6_dp)
timer = omp_get_wtime() - timer
call cpu_time(clock_t_end)


call read_coord('inputs/grid1.coord.txt', coord)
call read_topo('inputs/grid1.topo.txt', topo)
call get_boundaries(coord, 1e-5_dp, bnodes)
print *, 'bnodes:'
print *, (bnodes(i), i=1,10)

!call print_vec_compact(x, 5)
call write_vec('solution.txt', x)
print '(a20,1en20.3)', 'Elapsed time (s):', timer
print '(a20,1en20.3)', 'Clock time (s):', clock_t_end - clock_t_start

! Deallocate the matrix
ierr = dlt_CSRMAT(mat_A)
if (ierr /= 0) stop 'Error in deallocating mat_A'

! Deallocate b and vec_y
deallocate(b, vec_y, coord, topo, bnodes, stat=ierr)
if (ierr .ne. 0) stop 'ERROR: deallocation error in main.f90'

write(*,'(a)') 'Done'
end program main
