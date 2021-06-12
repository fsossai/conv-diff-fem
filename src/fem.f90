subroutine solve(coord, topo)
    use class_CSRMAT
    use blas, only: det3x3, inv3x3
    implicit none

    real(dp), intent(in)        :: coord(:,:)
    integer, intent(in)         :: topo(:,:)

    type(CSRMAT)                :: H, B, P
    real(dp), allocatable       :: q(:)
    integer                     :: i, ne, nn
    integer                     :: nodes(3)
    real(dp), dimension(3,3)    :: T, C, He, Be, Pe
    real(dp)                    :: area, grad(2,3), diff(2,2), qe, &
                                   ones3x1(3,1), vel(1,2)

    diff = 0.0_dp               ! Diffusivity coefficients
    diff(1,1) = 1.0_dp
    diff(2,2) = 1.0_dp
    ones3x1 = 1.0_dp
    vel = 1e-2_dp               ! (constant) velocity field

    ne = size(topo, 1)          ! Number of elements
    nn = size(coord, 1)         ! Number of nodes

    if (size(topo, 2) .ne. 3) stop 'ERROR: elements are not triangles.'

    allocate(q(nn))             ! RHS vector
    q = 0.0_dp

    do i = 1, ne
        nodes = topo(i, :)
        T = 1.0_dp
        T(:, 2:3) = coord(nodes, :)
        area = abs(det3x3(T)) / 2.0_dp
        call inv3x3(T, C)
        grad = C(2:3, :)
        He = area * matmul(transpose(grad), matmul(diff, grad))
        qe = area / 3.0_dp
        Be = matmul(ones3x1, matmul(vel, grad) / 6.0_dp)

        ! update global matrices
        ! H(nodes, nodes) = H(nodes, nodes) + He
        ! B(nodes, nodes) = B(nodes, nodes) + Be
        ! q(nodes) = q(nodes) + qe
    end do

    
    deallocate(q)

end subroutine

subroutine spmat_update(A, indices, Ae)
    use class_CSRMAT

    type(CSRMAT), intent(inout) :: A
    integer, intent(in)         :: indices(:)
    real(dp), intent(in)        :: Ae(:,:)

    integer :: i, j, k, n, offset, ii, jj
    real(dp), pointer :: coef(:) => null()
    integer, dimension(:), pointer :: iat => null(), ja => null()

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
            ! loop is likely to cycle forever.
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
    use class_CSRMAT
    use utils

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
        iat(i + 1) = iat(i) + ncols
        ja(iat(i) : iat(i) + ncols - 2) = conn(1 : ncols - 1, i)
        ja(iat(i + 1) - 1) = i
        call quicksort(ja, iat(i), iat(i+1)-1)
    end do

    if (associated(A%coef)) deallocate(A%coef)
    allocate(A%coef(nz))

    A%coef = 1.0_dp
    A%patt%iat => iat
    A%patt%ja => ja

end subroutine