! Basic Linear Algebra Subroutines

module BLAS
    use class_precision
    use class_CSRMAT
    use omp_lib
    implicit none

contains

    subroutine reset(x)
        real(dp), intent(out) :: x(:)
        integer :: i,n
        n = size(x)
        !$omp parallel do
        do i = 1,n
            x(i) = 0.0_dp
        end do
    end subroutine

    subroutine axpby(z, alpha, x, beta, y)
        real(dp), intent(inout) :: z(:)
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

        ! computation
        if (present(y)) then
            if (present(beta)) then
                beta_local = beta
            else
                beta_local = 1.0_dp
            end if
            !$omp parallel do
            do i = 1,n
                z(i) = z(i) + alpha * x(i) + beta_local * y(i)
            end do
        else
            !$omp parallel do
            do i = 1,n
                z(i) = z(i) + alpha * x(i)
            end do
        end if
        
    end subroutine

    subroutine axpby_set(z, alpha, x, beta, y)
        real(dp), intent(inout) :: z(:)
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

        ! computation
        if (present(y)) then
            if (present(beta)) then
                beta_local = beta
            else
                beta_local = 1.0_dp
            end if
            !$omp parallel do
            do i = 1,n
                z(i) = alpha * x(i) + beta_local * y(i)
            end do
        else
            !$omp parallel do
            do i = 1,n
                z(i) = alpha * x(i)
            end do
        end if
        
    end subroutine

    subroutine amxpby(z, alpha, A, x, beta, y)
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: alpha
        type(CSRMAT), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), optional, intent(in) :: beta
        real(dp), optional, intent(in) :: y(:)
        real(dp) :: beta_local, partial
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
            !$omp parallel do private(j,partial)
            do i = 1,n
                c_start = iat(i)
                c_end = iat(i+1) - 1
                partial = z(i)
                do j = c_start,c_end
                    partial = partial + coef(j) * x(ja(j))
                end do
                z(i) = alpha * partial + beta * y(i)
            end do
        else
            !$omp parallel do private(j,partial)
            do i = 1,n
                c_start = iat(i)
                c_end = iat(i+1) - 1
                partial = z(i)
                do j = c_start,c_end
                    partial = partial + coef(j) * x(ja(j))
                end do
                z(i) = alpha * partial
            end do
        end if
    end subroutine

    subroutine amxpby_set(z, alpha, A, x, beta, y)
        real(dp), intent(inout) :: z(:)
        real(dp), intent(in) :: alpha
        type(CSRMAT), intent(in) :: A
        real(dp), intent(in) :: x(:)
        real(dp), optional, intent(in) :: beta
        real(dp), optional, intent(in) :: y(:)
        real(dp) :: beta_local, partial
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
            !$omp parallel do private(j,partial,c_start,c_end)
            do i = 1,n
                c_start = iat(i)
                c_end = iat(i+1) - 1
                partial = 0.0_dp
                do j = c_start,c_end
                    partial = partial + coef(j) * x(ja(j))
                end do
                z(i) = alpha * partial + beta * y(i)
            end do
        else
            ! this loop is the bottleneck of the whole BiCGSTAB
            ! unrolled version
            !!$omp do private(j,partial,c_start,c_end)
            !do i = 1,n
            !    c_start = iat(i)
            !    c_end = iat(i+1) - 1
            !    partial = 0.0_dp
            !    do j = c_start, c_start + mod(c_end - c_start + 1, 4) - 1
            !        partial = partial + coef(j) * x(ja(j))
            !    end do
            !    do j = j,c_end,4  ! 4-unrolled
            !        partial = partial &
            !                + coef(j) * x(ja(j)) &
            !                + coef(j+1) * x(ja(j+1)) &
            !                + coef(j+2) * x(ja(j+2)) &
            !                + coef(j+3) * x(ja(j+3))
            !    end do
            !    z(i) = alpha * partial
            !end do
            !$omp parallel do private(j,partial,c_start,c_end)
            do i = 1,n
                c_start = iat(i)
                c_end = iat(i+1) - 1
                partial = 0.0_dp
                do j = c_start,c_end
                    partial = partial + coef(j) * x(ja(j))
                end do
                z(i) = alpha * partial
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
        !$omp parallel do reduction(+:z)
        do i = 1,n
            z = z + x(i) * y(i)
        end do
    end function

    function norm(x) result(z)
        real(dp), intent(in) :: x(:)
        real(dp) :: z
        z = sqrt(inner_prod(x, x))
    end function

end module BLAS