subroutine read_topo(filename, T)
    ! Read a 2D topology matrix.
    ! Each element is a triangle described by three node IDs.
    ! The first row must contain the total number of elements,
    ! which must corresponds to the number of remaining rows.
    ! Each row has the following syantax:
    ! ELEMENT_ID NODE1_ID NODE2_ID NODE3_ID

    use class_precision
    
    implicit none
    
    character(len=100), intent(in)      :: filename
    integer, allocatable, intent(out)   :: T(:,:)
    
    integer :: i, iunit, ierr, ne, el
    
    iunit = 1

    ! Open and read the matrix file
    open(iunit, file=filename, status='old')
    read(iunit, *) ne

    allocate(T(ne,3), stat=ierr)
    if (ierr .ne. 0) then
        stop 'ERROR: topology matrix allocation failed.'
    end if

    do i=1,ne
        read(iunit, *) el, T(i,1), T(i,2), T(i,3)
        if (el .ne. i) stop 'ERROR: unexpected syntax in topology data.'
    end do

end subroutine read_topo
    