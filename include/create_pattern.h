
interface create_pattern_int

    subroutine create_pattern(nnodes, topo, A)
        use class_CSRMAT
        integer, intent(in)             :: nnodes
        integer, intent(in)             :: topo(:,:)
        type(CSRMAT), intent(inout)     :: A
    end subroutine

end interface create_pattern_int