interface readers_int

    subroutine read_mat(filename, binary, mat, info)
        use class_CSRMAT
        implicit none

        character(len=*), intent(in)   :: filename
        logical, intent(in)              :: binary
        type(CSRMAT), intent(inout) :: mat
        integer, intent(out)        :: info

    end subroutine

    subroutine read_topo(filename, T)
        implicit none
        
        character(len=*), intent(in)      :: filename
        integer, allocatable, intent(out)   :: T(:,:)
    end subroutine

    subroutine read_coord(filename, C)
        use class_precision
        implicit none
        
        character(len=*), intent(in)      :: filename
        real(dp), allocatable, intent(out)  :: C(:,:)
    end subroutine

end interface readers_int
