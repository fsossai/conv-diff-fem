interface fem_int

    subroutine solve(coord, topo)
        use class_CSRMAT
        implicit none
        
        real(dp), intent(in)    :: coord(:,:)
        integer, intent(in)    :: topo(:,:)
    end subroutine

end interface fem_int