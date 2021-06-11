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