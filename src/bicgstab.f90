! Bi-conjugate Gradient Stabilized

module class_bicgstab
    implicit none
contains

subroutine bicgstab(A, b, x)
    use class_precision
    use class_CSRMAT
    use class_BLAS

    type(CSRMAT), intent(in)            :: A
    real(dp), intent(in)                :: b(:)
    real(dp), intent(out)               :: x(:)
    integer :: j,n
    real(dp) :: alpha, beta, omega, temp1
    real(dp), allocatable, dimension(:) :: r, r0, r_new, p, s, temp2, temp3
    integer, parameter :: max_it = 1

    ! handles for the matrix
    real(dp), pointer :: coef(:) => null()
    integer, dimension(:), pointer :: iat => null(), ja => null()

    n = size(b)
    ! dimension check
    if (A%patt%nrows .ne. n) then
        print *, 'ERROR: dimension mismatch (in BiCGSTAB).'
        stop
    end if

    ! setting handles 
    coef => A%coef
    iat => A%patt%iat
    ja => A%patt%ja

    allocate(r(n), r0(n), r_new(n), p(n), s(n), temp2(n), temp3(n))

    ! Start of BiCGSTAB algorithm
    ! r_0^ = b - A x_0
    call amxpby(r, -1.0_dp, A, x, 1.0_dp, b)
    
    ! r_0^ arbitrary
    r0 = r
    
    ! p0 = r0
    p = r
    
    do j = 1,max_it
        ! alpha_j = (r_j,r_0^) / (A p_j, r_0^)
        temp1 = inner_prod(r, r0)
        call amxpby(temp2, 1.0_dp, A, p)
        alpha = temp1 / inner_prod(temp2, r0)
        
        ! s_j = r_j - alpha_j A p_j
        s = r - alpha * temp2
        
        ! omega_j = (A s_j, s_j) / (A s_j, A s_j)
        call amxpby(temp3, 1.0_dp, A, s)
        omega = inner_prod(temp3, s) / inner_prod(temp3, temp3)

        ! x_j+1 = x_j + alpha_j p_j + omega_j s_j
        call axpby(x, alpha, p, 1.0_dp, x)
        call axpby(x, omega, s, 1.0_dp, x)
        
        ! r_j+1 = s_j - omega_j A s_j
        call axpby(r_new, -omega, temp3, 1.0_dp, s)
        
        ! beta_j = (r_j+1, r_0^) / (r_j, r_0^) * (alpha_j / omega_j)
        beta = inner_prod(r_new, r0) / temp1 * (alpha / omega)

        ! p_j+1 = r_j+1 + beta_j (p_j - omega(j) A p_j)
        p = r_new + beta * (p - omega * temp2)

        r = r_new
    end do
    
end subroutine


end module