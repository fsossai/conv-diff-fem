! Basic Linear Algebra Subroutines

module class_BLAS
    use class_precision
    use class_CSRMAT
    implicit none

contains

    subroutine reset(x)
        real(dp), intent(out) :: x(:)
        integer :: i,n
        n = size(x)
        do i = 1,n
            x(i) = 0.0_dp
        end do
    end subroutine

    subroutine axpby(alpha, x, beta, y, z)
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), intent(in) :: alpha, beta
        real(dp), allocatable, intent(out) :: z(:)
        integer :: i,n

        n = size(x)

        ! dimension check
        if (n .ne. size(y)) then
            print *, 'ERROR: dimension mismatch (in axpby).'
            stop
        end if

        ! automatic memory allocation
        if (.not.allocated(z)) then
            allocate(z(n))
            call reset(z)
        end if

        ! computation
        do i = 1,n
            z(i) = alpha * x(i) + beta * y(i)
        end do
    end subroutine

    subroutine amxpby(alpha, A, x, beta, y, z)
        type(CSRMAT), intent(in) :: A
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), allocatable, intent(out) :: z(:)
        real(dp), intent(in) :: alpha, beta
        integer :: i,j,n,c_start,c_end

        ! handles for the matrix
        real(dp), pointer :: coef(:) => null()
        integer, dimension(:), pointer :: iat => null(), ja => null()

        n = A%patt%nrows

        ! dimension check
        if (n .ne. size(x)) then
            print *, 'ERROR: dimension mismatch (in amxpby).'
            stop
        end if
        
        ! automatic memory allocation
        if (.not.allocated(z)) then
            allocate(z(n))
            call reset(z)
        end if

        ! setting handles 
        coef => A%coef
        iat => A%patt%iat
        ja => A%patt%ja

        ! computation
        do i=1,n
            c_start = iat(i)
            c_end = iat(i+1) - 1
            do j = c_start,c_end
                z(i) = z(i) + coef(j) * x(ja(j))
            end do
            z(i) = alpha * z(i) + beta * y(i)
        end do
    end subroutine

end module class_BLAS