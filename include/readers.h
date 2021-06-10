interface read_mat_int
    subroutine read_mat(filename, binary, mat, info)
        use class_CSRMAT
        implicit none

        character(len=100), intent(in)   :: filename
        logical, intent(in)              :: binary
        type(CSRMAT), intent(inout) :: mat
        integer, intent(out)        :: info

    end subroutine
end interface read_mat_int
