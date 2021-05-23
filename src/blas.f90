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

    subroutine axpby(z, alpha, x, beta, y)
        real(dp), allocatable, intent(inout) :: z(:)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: x(:)
        real(dp), optional, intent(in) :: beta
        real(dp), optional, intent(in) :: y(:)
        real(dp) :: beta_local
        integer :: i,n

        n = size(x)

        ! dimension check
        if (present(y) .and. n .ne. size(y)) then
            print *, 'ERROR: dimension mismatch (in axpby).'
            stop
        end if

        ! automatic memory allocation
        if (.not.allocated(z)) then
            allocate(z(n))
            call reset(z)
        end if

        ! computation
        if (present(y)) then
            if (present(beta)) then
                beta_local = beta
            else
                beta_local = 1.0_dp
            end if
            do i = 1,n
                z(i) = z(i) + alpha * x(i) + beta_local * y(i)
            end do
        else
            do i = 1,n
                z(i) = z(i) + alpha * x(i)
            end do
        end if

        
    end subroutine

    subroutine amxpby(z, alpha, A, x, beta, y)
        real(dp), allocatable, intent(out) :: z(:)
        real(dp), intent(in) :: alpha
        type(CSRMAT), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), optional, intent(in) :: beta
        real(dp), optional, intent(in) :: y(:)
        real(dp) :: beta_local
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
        if (present(y)) then
            if (present(beta)) then
                beta_local = beta
            else
                beta_local = 1.0_dp
            end if
            do i=1,n
                c_start = iat(i)
                c_end = iat(i+1) - 1
                do j = c_start,c_end
                    z(i) = z(i) + coef(j) * x(ja(j))
                end do
                z(i) = alpha * z(i) + beta * y(i)
            end do
        else
            do i=1,n
                c_start = iat(i)
                c_end = iat(i+1) - 1
                do j = c_start,c_end
                    z(i) = z(i) + coef(j) * x(ja(j))
                end do
                z(i) = alpha * z(i)
            end do
        end if
    end subroutine

    function inner_prod(x, y) result(z)
        real(dp) :: z
        real(dp), dimension(:), intent(in) :: x, y
        integer :: i,n

        n = size(x)

        ! dimension check
        if (n .ne. size(y)) then
            print *, 'ERROR: dimension mismatch (in axpby).'
            stop
        end if

        ! computation
        z = 0.0_dp
        do i = 1,n
            z = z + x(i) * y(i)
        end do
    end function

end module class_BLAS