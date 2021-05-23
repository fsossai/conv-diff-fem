! Basic Linear Algebra Subroutines

module class_BLAS
    use class_precision
    use class_CSRMAT
    implicit none

contains
    
    subroutine axpy(alpha, x, y, z)
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), intent(in) :: alpha
        real(dp), allocatable, intent(out) :: z(:)
        integer :: n

        n = size(x)

        ! dimension check
        if (n .ne. size(y)) then
            print *, 'ERROR: dimension mismatch (in axpy).'
            stop
        end if

        ! automatic memory allocation
        if (.not.allocated(z)) allocate(z(n))

        ! computation

    end subroutine

    subroutine mxv(A, x, z)
        type(CSRMAT), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: z(:)
        integer :: n
        ! handles for the matrix
        real(dp), pointer :: coef(:) => null()
        integer, dimension(:), pointer :: iat => null(), ja => null()

        n = A%patt%nrows

        ! dimension check
        if (n .ne. size(x)) then
            print *, 'ERROR: dimension mismatch (in mxv).'
            stop
        end if
        
        ! automatic memory allocation
        if (.not.allocated(z)) allocate(z(n))

        ! setting handles 
        coef => A%coef
        iat => A%patt%iat
        ja => A%patt%ja

        ! computation


    end subroutine

end module class_BLAS