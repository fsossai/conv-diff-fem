interface spmat_update_int

    subroutine spmat_update(A, indices, Ae)
        use class_CSRMAT
        use class_precision
        type(CSRMAT), intent(inout) :: A
        integer, intent(in)         :: indices(:)
        real(dp), intent(in)        :: Ae(:,:)
    end subroutine
    
end interface spmat_update_int

