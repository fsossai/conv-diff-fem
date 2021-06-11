subroutine solve(coord, topo)
    use class_CSRMAT
    use blas, only: det3x3, inv3x3
    implicit none

    real(dp), intent(in)        :: coord(:,:)
    integer, intent(in)         :: topo(:,:)

    !type(CSRMAT)                :: H, B, P
    real(dp), allocatable       :: q(:)
    integer                     :: i, ne, nn
    integer                     :: nodes(3)
    real(dp), dimension(3,3)    :: T, C
    real(dp)                    :: area

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
    end do

    deallocate(q)

end subroutine