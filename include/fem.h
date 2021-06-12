interface solve_int

    subroutine solve(coord, topo)
        use class_CSRMAT
        real(dp), intent(in)    :: coord(:,:)
        integer, intent(in)    :: topo(:,:)
    end subroutine

end interface solve_int

include 'create_pattern.h'
include 'spmat_update.h'