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
    !!$omp parallel do
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
        !!$omp parallel do
        do i = 1,n
            z(i) = z(i) + alpha * x(i) + beta_local * y(i)
        end do
    else
        !!$omp parallel do
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
        !!$omp parallel do
        do i = 1,n
            z(i) = alpha * x(i) + beta_local * y(i)
        end do
    else
        !!$omp parallel do
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
        !!$omp parallel do private(j,partial)
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
        !!$omp parallel do private(j,partial)
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
        !!$omp parallel do private(j,partial,c_start,c_end)
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
        !!!$omp do private(j,partial,c_start,c_end)
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
    !!$omp parallel do reduction(+:z)
    do i = 1,n
        z = z + x(i) * y(i)
    end do
end function


function norm(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp) :: z
    z = sqrt(inner_prod(x, x))
end function


function det3x3(A) result(det)
    real(dp), intent(in)    :: A(3,3)
    real(dp)                :: det
    det = A(1,1) * A(2,2) * A(3,3) &
        - A(1,1) * A(2,3) * A(3,2) &
        - A(1,2) * A(2,1) * A(3,3) &
        + A(1,2) * A(2,3) * A(3,1) &
        + A(1,3) * A(2,1) * A(3,2) &
        - A(1,3) * A(2,2) * A(3,1)
end function


subroutine inv3x3(A, A_inv)
    real(dp), intent(in)    :: A(3,3)
    real(dp), intent(out)   :: A_inv(3,3)
    
    real(dp), parameter     :: eps = 1.0e-10_dp
    real(dp)                :: cofactor(3,3), det
    
    det = A(1,1) * A(2,2) * A(3,3) &
        - A(1,1) * A(2,3) * A(3,2) &
        - A(1,2) * A(2,1) * A(3,3) &
        + A(1,2) * A(2,3) * A(3,1) &
        + A(1,3) * A(2,1) * A(3,2) &
        - A(1,3) * A(2,2) * A(3,1)
    
    if (abs(det) .le. eps) then
        A_inv = 0.0_dp
        return
    end if
    
    cofactor(1,1) = + ( A(2,2)*A(3,3) - A(2,3)*A(3,2) )
    cofactor(1,2) = - ( A(2,1)*A(3,3) - A(2,3)*A(3,1) )
    cofactor(1,3) = + ( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
    cofactor(2,1) = - ( A(1,2)*A(3,3) - A(1,3)*A(3,2) )
    cofactor(2,2) = + ( A(1,1)*A(3,3) - A(1,3)*A(3,1) )
    cofactor(2,3) = - ( A(1,1)*A(3,2) - A(1,2)*A(3,1) )
    cofactor(3,1) = + ( A(1,2)*A(2,3) - A(1,3)*A(2,2) )
    cofactor(3,2) = - ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )
    cofactor(3,3) = + ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )
    
    A_inv = transpose(cofactor) / det        
end subroutine inv3x3

subroutine jacobi_precond_mat(A)
    ! Divide rows by the diagonal element

    type(CSRMAT), intent(inout)         :: A

    real(dp), pointer               :: coef(:)
    integer, pointer, dimension(:)  :: iat, ja
    integer                         :: i, j, n
    real(dp)                        :: pivot

    ! setting handles
    coef => A%coef
    iat => A%patt%iat
    ja => A%patt%ja

    pivot = 1.0_dp
    n = A%patt%nrows

    do i = 1, n
        ! finding pivot, a.k.a. the diagonal element
        do j = iat(i), iat(i+1) - 1
            if (ja(j).eq.i) then
                if (coef(j).eq.0.0_dp) then
                    stop 'ERROR: diagonal element is zero (jacobi_precond)'
                end if
                pivot = 1.0_dp / coef(j)
                exit
            end if
        end do

        ! changing the i-th row
        do j = iat(i), iat(i+1) - 1
            coef(j) = pivot * coef(j)
        end do
    end do
end subroutine


subroutine jacobi_precond_rhs(A, b)
    ! Apply the preconditioning only to the right hand side

    type(CSRMAT), intent(in)    :: A
    real(dp), intent(inout)     :: b(:)

    real(dp), pointer           :: coef(:)
    integer, pointer            :: iat(:), ja(:)
    integer                     :: i, j, n
    real(dp)                    :: pivot

    ! setting handles
    coef => A%coef
    iat => A%patt%iat
    ja => A%patt%ja

    pivot = 1.0_dp
    n = A%patt%nrows
    if (n.ne.size(b)) stop 'ERROR: wrong system shape (jacobi_precond)'

    do i = 1, n
        ! finding pivot, a.k.a. the diagonal element
        do j = iat(i), iat(i+1) - 1
            if (ja(j).eq.i) then
                if (coef(j).eq.0.0_dp) then
                    stop 'ERROR: diagonal element is zero (jacobi_precond)'
                end if
                pivot = 1.0_dp / coef(j)
                exit
            end if
        end do

        b(i) = pivot * b(i)
    end do
end subroutine

end module BLAS