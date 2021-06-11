interface fem_int

    subroutine solve(coord, topo)
        use class_CSRMAT
        real(dp), intent(in)    :: coord(:,:)
        integer, intent(in)    :: topo(:,:)
    end subroutine


    subroutine spmat_update(A, indices, Ae)
        use class_CSRMAT
        type(CSRMAT), intent(inout) :: A
        integer, intent(in)         :: indices(:)
        real(dp), intent(in)        :: Ae(:,:)
    end subroutine

    
    subroutine create_pattern(nnodes, topo, A)
        use class_CSRMAT
        integer, intent(in)             :: nnodes
        integer, intent(in)             :: topo(:,:)
        type(CSRMAT), intent(inout)     :: A
    end subroutine
    
end interface fem_int