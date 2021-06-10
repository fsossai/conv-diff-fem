subroutine read_coord(filename, C)
    ! Read a 2D matrix of nodal coordinates.
    ! Each node is a point in two dimensions.
    ! The first row must contain the total number of nodes,
    ! which must corresponds to the number of remaining rows.
    ! Each row has the following syantax:
    ! NODE_ID X_COORD Y_COORD

    use class_precision
    
    implicit none
    
    character(len=*), intent(in)      :: filename
    real(dp), allocatable, intent(out)  :: C(:,:)
    
    integer :: i, iunit, ierr, nnodes, row
    
    iunit = 1

    ! Open and read the matrix file
    open(iunit, file=filename, status='old')
    read(iunit, *) nnodes

    allocate(C(nnodes,2), stat=ierr)
    if (ierr .ne. 0) then
        stop 'ERROR: nodal coordinates matrix allocation failed.'
    end if

    do i=1,nnodes
        read(iunit, *) row, C(i,1), C(i,2)
        if (row .ne. i) stop 'ERROR: unexpected syntax in topology data.'
    end do

end subroutine read_coord
    