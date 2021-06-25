program main
use class_precision
use class_CSRMAT
use omp_lib
use fem, only: assemble, create_pattern, compute_workloads
use utils, only: write_CSRMAT

include 'readers.h'

real(dp)                :: timer, delta
real(dp), allocatable   :: coord(:,:)
integer, allocatable    :: topo(:,:)
real(dp), parameter     :: boundary_cond = 5.0_dp, dt = 0.001_dp
type(CSRMAT)            :: H, P
real(dp), allocatable   :: q(:)
integer                 :: nnodes, nelem, ierr, it, i
integer, parameter      :: max_it = 1
character(100)          :: arg_coord, arg_topo
integer, pointer        :: i_start_of(:), i_end_of(:), offset(:), duty_of(:), el_idx(:)


! Checking input arguments
if (iargc().lt.2) then
    stop 'Missing arguments. Please specify coordinates and topology files.'
end if

! Getting command line arguments
call getarg(1, arg_coord)
call getarg(2, arg_topo)

! Reading input files
print '(a)', 'Reading coordinates file...'
call read_coord(arg_coord, coord)
print '(a)', 'Reading topology file...'
call read_topo(arg_topo, topo)

nnodes = size(coord, 1)
nelem = size(topo, 1)

allocate(q(nnodes))

print '(a25,i15)', 'Number of nodes:',      nnodes
print '(a25,i15)', 'Number of elements:',   nelem
print '(a25,i15)', 'Number of iterations:', max_it

call create_pattern(nnodes, topo, H)
call copy_Pattern(P, H)
call compute_workloads(topo, nnodes, i_start_of, i_end_of, offset, duty_of, el_idx)

! Assembling the system matrix
timer = 0.0_dp
do it = 1, max_it
    ! Reset
    H%coef = 0.0_dp
    P%coef = 0.0_dp
    q = 0.0_dp
    
    ! Assembly
    delta = omp_get_wtime()
    call assemble(coord, topo, dt, H, P, q, i_start_of, i_end_of, offset, duty_of, el_idx)
    timer = timer + omp_get_wtime() - delta
end do


print '(a25,en15.3)', 'Elapsed time [s]:',  timer

! Writing to file
print '(a)', 'Exporting matrix to file...'
call write_CSRMAT('mat.txt', H)

! Showing the first row
print *, 'First row:'
print '(*(en20.8,1x))', (H%coef(i), i=H%patt%iat(1),H%patt%iat(2)-1)

! Deleting/deallocating all matrices
ierr = 0
ierr = ierr + dlt_CSRMAT(H)
ierr = ierr + dlt_CSRMAT(P)
if (ierr.ne.0) stop 'ERROR: failed to delete one or more CSRMAT.'
deallocate(q, i_start_of, i_end_of, offset, duty_of, el_idx)

end program
