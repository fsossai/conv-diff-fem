module fem
    use class_precision
    use class_CSRMAT
    use blas
    use utils
    implicit none

contains


subroutine assemble(coord, topo, dt, H, P, q, i_start_of, i_end_of, offset, duty_of, el_idx)
    real(dp), intent(in)                :: coord(:,:)
    integer, intent(in)                 :: topo(:,:)
    real(dp), intent(in)                :: dt
    type(CSRMAT), intent(inout)         :: H, P
    real(dp), intent(out)               :: q(:)
    integer, intent(in), dimension(:)   :: i_start_of, i_end_of, offset, duty_of, el_idx

    integer                             :: i, ne, nn, nthreads, tid, i_start, i_end
    integer                             :: nodes(3), idx(9), regions(3)
    real(dp), dimension(3,3)            :: T, C
    real(dp), target, dimension(3,3)    :: He, Pe
    real(dp), pointer, dimension(:)     :: He_flat, Pe_flat
    real(dp)                            :: area, grad(2,3), diff(2,2), qe, &
                                           ones3x1(3,1), vel(1,2), time
    logical                             :: frontier
    real(dp), allocatable               :: local_H(:), local_P(:)
    integer, pointer                    :: iat(:) => null()

    diff = 0.0_dp               ! Diffusivity coefficients
    diff(1,1) = 1.0_dp
    diff(2,2) = 1.0_dp
    ones3x1 = 1.0_dp
    vel = 1e-2_dp               ! (constant) velocity field

    ne = size(topo, 1)          ! Number of elements
    nn = size(coord, 1)         ! Number of nodes

    if (size(topo, 2).ne.3) stop 'ERROR: elements are not triangles.'

    !$omp parallel shared(H,P,q,ne,coord,topo,duty_of,el_idx,i_start_of,i_end_of,offset) &
    !$omp& firstprivate(diff,vel,ones3x1,dt,nn) & 
    !$omp& private(i,nodes,T,C,area,grad,He,Pe,qe,He_flat,Pe_flat,idx,regions,nthreads,tid) &
    !$omp& private(frontier,local_H,local_P,iat,i_start,i_end)
    
    !$omp workshare
    q = 0.0_dp
    !$omp end workshare

    nthreads = omp_get_num_threads()
    tid = omp_get_thread_num()

    if (size(i_start_of).ne.nthreads) then
        stop 'ERROR: unexpected number of threads'
    end if

    ! Flattening local matrices, so that we will be able to update
    ! the array 'coef' in a single assignment.
    He_flat(1:9) => He(:,:)
    Pe_flat(1:9) => Pe(:,:)

    ! local_H and local_P will serve as private partitions of the
    ! global matrices, kept locally for performance reasons.
    ! The total memory footprint will be twice the memory needed
    ! to store the sparse matrices only.
    ! These subdomains will contribute to the global matrices only
    ! at the end.
    iat => H%patt%iat
    allocate(local_H( iat(offset(tid + 1)) - iat(offset(tid)) ))
    allocate(local_P( iat(offset(tid + 1)) - iat(offset(tid)) ))
    !$omp workshare
    local_H = 0.0_dp
    local_P = 0.0_dp
    !$omp end workshare
    
    !$omp single
    time = omp_get_wtime()
    !$omp end single
    
    ! Computing bounds 
    i_start = i_start_of(tid)
    i_end = i_end_of(tid)

    ! Every element is processed by exactly one thread
    do i = i_start, i_end
        ! el_idx(i) is the index of the current element    
        nodes = topo(el_idx(i), :)
        
        T = 1.0_dp
        T(:, 2:3) = coord(nodes, :)
        area = abs(det3x3(T)) / 2.0_dp
        call inv3x3(T, C)
        grad = C(2:3, :)
        
        Pe = area / 12.0_dp
        Pe(1, 1) = area / 6.0_dp
        Pe(2, 2) = area / 6.0_dp
        Pe(3, 3) = area / 6.0_dp
        
        He = area * matmul(transpose(grad), matmul(diff, grad)) &
           + matmul(ones3x1, matmul(vel, grad) / 6.0_dp)        &
           + Pe / dt
        
        qe = area / 3.0_dp
        
        ! Getting the index of the entry to be updated
        ! in the sparse matrix
        idx = get_idx(H, nodes)
        
        ! update the global matrices
        if (duty_of(i).gt.0) then   ! concurrent update
            idx = idx - iat(offset(tid)) + 1
            local_H(idx) = local_H(idx) + He_flat
            local_P(idx) = local_P(idx) + Pe_flat
            q(nodes) = q(nodes) + qe
        else                        ! atomic update
            !$omp critical
            H%coef(idx) = H%coef(idx) + He_flat        
            P%coef(idx) = P%coef(idx) + Pe_flat
            q(nodes) = q(nodes) + qe
            !$omp end critical
        end if
       
    end do

    ! Final update of the local subdomain of the matrices into
    ! the global ones.
    !$omp critical
    i_start = iat(offset(tid))
    i_end = iat(offset(tid + 1)) - 1
    H%coef(i_start:i_end) = H%coef(i_start:i_end) + local_H
    P%coef(i_start:i_end) = P%coef(i_start:i_end) + local_P
    !$omp end critical

    deallocate(local_H, local_P)

    !$omp end parallel

    time = omp_get_wtime() - time

    print '(a25,en20.3)', 'Elements proc time [s]:', time

end subroutine


subroutine compute_workloads(topo, nn, i_start_of, i_end_of, offset, duty_of, el_idx)
    integer, intent(in)                 :: topo(:,:)
    integer, intent(in)                 :: nn
    integer, intent(out), allocatable   :: i_start_of(:), i_end_of(:), offset(:)
    integer, intent(out), dimension(:)  :: duty_of, el_idx

    integer                             :: i, ne, nthreads, tid, tmp
    integer                             :: nodes(3), regions(3)
    integer, allocatable                :: workload(:)

    ne = size(topo, 1)          ! Number of elements

    if (size(topo, 2).ne.3) stop 'ERROR: elements are not triangles.'

    !$omp parallel shared(topo,duty_of,workload,el_idx,nn,i_start_of,i_end_of,offset) &
    !$omp& private(i,nodes,regions,nthreads,tid,tmp)

    nthreads = omp_get_num_threads()
    tid = omp_get_thread_num()

    !$omp single
    allocate(workload(0:nthreads), i_start_of(0:nthreads-1), i_end_of(0:nthreads-1))
    allocate(offset(0:nthreads))
    workload = 0
    !$omp end single

    ! Calculating, given an element, which is the thread that will
    ! take care of the computation.
    !$omp do
    do i = 1, ne
        nodes = topo(i, :)

        ! Getting which thread every node of 'nodes' belongs to
        call get_regions(nodes, nn, nthreads, regions)
        
        if (regions(1).eq.regions(2).and.regions(2).eq.regions(3)) then
            ! Element i will be computed by thread 'regions(1)'
            duty_of(i) = regions(1) + 1
        else
            ! Here, the minus sign indicated the need for a critical section
            duty_of(i) = - (minval(regions) + 1)
        end if

        tmp = abs(duty_of(i))

        !$omp atomic
        workload(tmp) = workload(tmp) + 1
        !$omp end atomic
    end do
    
    !$omp workshare
    el_idx = [(i, i = 1, ne)]
    !$omp end workshare

    !$omp single
    call paired_quicksort_abs(duty_of, el_idx, 1, ne)
    
    ! Calculating offsets.
    ! offset(i) is the first internal node of the subdomain
    ! associated to the i-th thread.
    ! Additionally, offset(i+1) - offset(i) counts how many
    ! nodes belongs to the subdomain of the i-th thread.
    offset(0) = 1
    do i = 1, nthreads
        offset(i) = offset(i - 1) + (nn / nthreads)
        if (i.le.modulo(nn, nthreads)) then
            offset(i) = offset(i) + 1
        end if
    end do
    !$omp end single

    ! Computing bounds
    i_start_of(tid) = 1 + sum(workload(0:tid))
    i_end_of(tid) = i_start_of(tid) + workload(tid + 1) - 1

    !$omp end parallel

    deallocate(workload)
    
end subroutine


subroutine spmat_update(A, indices, Ae)
    type(CSRMAT), intent(inout) :: A
    integer, intent(in)         :: indices(:)
    real(dp), intent(in)        :: Ae(:,:)

    integer                     :: i, j, k, n, offset, ii, jj
    real(dp), pointer           :: coef(:) => null()
    integer, pointer            :: iat(:) => null(), ja(:) => null()

    n = size(indices)

    ! Sanity check
    if (size(Ae, 1) .ne. n .or. size(Ae, 2) .ne. n) then
        stop 'ERROR: sanity check failed in spmat_update'
    end if

    ! handles to the sparse matrix
    coef => A%coef
    iat => A%patt%iat
    ja => A%patt%ja
    
    do i = 1, n
        do j = 1, n
            ! Finding the position of (indices(i), indices(j)) in
            ! the sparse matrix coefficient array

            ii = indices(i)
            jj = indices(j)
            offset = iat(ii)

            ! We assume that the sparse matrix has a suitable
            ! pattern. If this condition is not met, the following
            ! loop will cause a segmentation fault.
            k = 0
            do while (ja(offset + k) .ne. jj)
                k = k + 1
            end do

            ! Finally, offset + k is the position of (ii, jj)
            coef(offset + k) = coef(offset + k) + Ae(i, j)
        end do        
    end do
    
end subroutine


subroutine create_pattern(nnodes, topo, A)
    integer, intent(in)             :: nnodes
    integer, intent(in)             :: topo(:,:)
    type(CSRMAT), intent(inout)     :: A

    integer                         :: i, nelem, elem(3), ncols, nz
    integer, parameter              :: max_degree = 8
    integer, allocatable            :: conn(:,:)
    integer, pointer, dimension(:)  :: iat, ja
        
    nelem = size(topo, 1)
    allocate(conn(max_degree, nnodes))
    conn = -1
    
    do i = 1, nelem
        elem = topo(i, :)

        ! storing that elem(1) is connected with elem(2) and elem(3)
        call set_insert(conn(:, elem(1)), elem(2))
        call set_insert(conn(:, elem(1)), elem(3))

        ! storing that elem(2) is connected with elem(1) and elem(3)
        call set_insert(conn(:, elem(2)), elem(1))
        call set_insert(conn(:, elem(2)), elem(3))

        ! storing that elem(3) is connected with elem(1) and elem(2)
        call set_insert(conn(:, elem(3)), elem(1))
        call set_insert(conn(:, elem(3)), elem(2))
    end do

    ! At this point, the i-th column of 'conn' is the column pointer
    ! of the i-th row in the sparse matrix.
    ! Now the actual pattern is ready to be built.
  
    iat => A%patt%iat
    ja => A%patt%ja
    
    ! Calculating the number of nonzeros
    nz = max_degree * nnodes - count(conn .eq. -1) + nnodes

    if (associated(iat)) deallocate(iat)
    allocate(iat(nnodes + 1))

    if (associated(ja)) deallocate(ja)
    allocate(ja(nz))

    A%patt%nrows = nnodes
    A%patt%nterm = nz

    iat(1) = 1
    do i = 1, nnodes
        ncols = index_of_first(conn(:, i), -1)
        if (ncols.eq.-1) ncols = size(conn, 1) + 1
        iat(i + 1) = iat(i) + ncols
        ja(iat(i) : iat(i) + ncols - 2) = conn(1 : ncols - 1, i)
        ja(iat(i + 1) - 1) = i
        call quicksort(ja, iat(i), iat(i+1)-1)
    end do

    if (associated(A%coef)) deallocate(A%coef)
    allocate(A%coef(nz))

    A%coef = 0.0_dp
    A%patt%iat => iat
    A%patt%ja => ja

end subroutine


subroutine solve(coord, topo, x0, x)
    real(dp), intent(in)    :: coord(:,:)
    integer, intent(in)     :: topo(:,:)
    real(dp), intent(in)    :: x0(:)
    real(dp), intent(out)   :: x(:)

    integer, parameter      :: max_it = 100
    integer, parameter      :: bicgstab_max_it = 200
    real(dp), parameter     :: boundary_cond = 5.0_dp
    type(CSRMAT)            :: H, P
    real(dp), allocatable   :: q(:), rhs(:), diag(:)
    integer, allocatable    :: bnodes(:), i_start_of(:), i_end_of(:), &
                               offset(:), duty_of(:), el_idx(:)
    integer                 :: nnodes, nelem, i, ierr
    real(dp)                :: dt = 0.01

    include 'bicgstab.h'
    
    nnodes = size(coord, 1)
    nelem = size(topo, 1)

    allocate(q(nnodes), rhs(nnodes), diag(nnodes))

    call compute_workloads(topo, nnodes, i_start_of, i_end_of, offset, duty_of, el_idx)
    
    call get_boundaries(coord, 1e-5_dp, bnodes)
    call create_pattern(nnodes, topo, H)
    call copy_Pattern(P, H)

    ! Setting up the the matrix and the RHS of the linear system
    call get_diag(H, diag)
    call jacobi_precond_mat(H)

    x = x0

    do i = 1, max_it
        print *, 'Iteration:', i

        call assemble(coord, topo, dt, H, P, q, i_start_of, i_end_of, offset, duty_of, el_idx)

        ! Setting up the the matrix and the RHS of the linear system
        call get_diag(H, diag)
        call jacobi_precond_mat(H)

        ! imposing the Dirichlet boundary conditions
        x(bnodes) = boundary_cond

        ! solving the linear system
        call amxpby_set(rhs, 1.0_dp / dt, P, x, 1.0_dp, q)

        ! preconditioning 
        rhs = rhs / diag
        call bicgstab(H, rhs, x, 1e-5_dp, bicgstab_max_it)
    end do

    ierr = 0
    ierr = ierr + dlt_CSRMAT(H)
    ierr = ierr + dlt_CSRMAT(P)
    if (ierr .ne. 0) stop 'ERROR: failed to delete one or more CSRMAT in "solve".'
    
    deallocate(q, rhs, diag, i_start_of, i_end_of, offset, duty_of, el_idx)

end subroutine


function has_diagonal(A, r) result(flag)
    ! Checks whether the diagonal of the sparse matrix A 
    ! belongs to the matrix.
    ! Every integer 'i' in the output array 'r' represents
    ! a row index such that A(i,i) is NOT part of the pattern.
    
    type(CSRMAT), intent(in)    :: A
    integer, intent(out), optional, allocatable  :: r(:)
    logical                     :: flag
    integer                     :: i, j, k, n
    integer, allocatable        :: temp(:)

    n = A%patt%nrows
    allocate(temp(n))
    temp = -1
    k = 0
    flag = .false.

    do i = 1, n
        flag = .false.
        do j = A%patt%iat(i), A%patt%iat(i+1) - 1
            if (A%patt%ja(j) .eq. i) then
                flag = .true.
                exit
            end if
        end do
        if (.not.flag) then
            k = k + 1
            temp(k) = i
        end if
    end do

    if (present(r)) then
        allocate(r(k))
        r = temp(1:k)
    end if
    
    if (temp(1) .eq. -1) then
        flag = .true.
    else
        flag = .false.
    end if
    deallocate(temp)

end function


function get_idx(A, nodes) result(idx)
    type(CSRMAT), intent(inout) :: A
    integer, intent(in)         :: nodes(3)
    integer                     :: idx(9)

    integer                     :: i, j, k, offset, ii, jj, h
    real(dp), pointer           :: coef(:) => null()
    integer, pointer            :: iat(:) => null(), ja(:) => null()

    ! handles to the sparse matrix
    coef => A%coef
    iat => A%patt%iat
    ja => A%patt%ja

    h = 0
    do j = 1, 3
        do i = 1, 3
            ! Finding the position of (indices(i), indices(j)) in
            ! the sparse matrix coefficient array

            ii = nodes(i)
            jj = nodes(j)
            offset = iat(ii)

            ! We assume that the sparse matrix has a suitable
            ! pattern. If this condition is not met, the following
            ! loop will cause a segmentation fault.
            k = 0
            do while (ja(offset + k).ne.jj)
                k = k + 1
            end do

            h = h + 1
            idx(h) = offset + k
        end do
    end do

end function


subroutine get_regions(nodes, nnodes, nthreads, regions)
    integer, intent(in)     :: nodes(3), nnodes
    integer, intent(in)     :: nthreads
    integer, intent(out)    :: regions(3)

    integer                 :: i, d, wload, rem, s

    wload = nnodes / nthreads
    rem = modulo(nnodes, nthreads)
    s = 0
    if (rem.ne.0) s = 1

    do i = 1, 3
        d = (nodes(i) - 1) / (wload + s)
        if (d.ge.rem) then
            regions(i) = (nodes(i) - 1 - (wload + s) * rem) / wload + rem
        else
            regions(i) = d
        end if
    end do

end subroutine


end module fem