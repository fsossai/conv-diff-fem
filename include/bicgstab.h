interface bicgstab_int
    subroutine bicgstab(A, b, x, tol, max_it)
        use class_CSRMAT

        type(CSRMAT), intent(in)            :: A
        real(dp), intent(in)                :: b(:)
        real(dp), intent(out)               :: x(:)
        real(dp), optional, intent(in)      :: tol
        integer, optional, intent(in)       :: max_it
    end subroutine
end interface