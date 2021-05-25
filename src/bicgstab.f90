! Bi-conjugate Gradient Stabilized

module class_bicgstab
    implicit none
contains

subroutine bicgstab(A, b, x, tol, max_it)
    use class_precision
    use class_CSRMAT
    use class_BLAS
    use class_utils
    use omp_lib

    type(CSRMAT), intent(in)            :: A
    real(dp), intent(in)                :: b(:)
    real(dp), intent(out)               :: x(:)
    real(dp), optional, intent(in)      :: tol
    integer, optional, intent(in)       :: max_it
    integer :: j, n, max_iterations
    real(dp) :: alpha, beta, omega, r_r0, tolerance
    real(dp), allocatable, dimension(:) :: r, r0, r_new, p, s, A_p, A_s

    ! handles for the matrix
    real(dp), pointer :: coef(:) => null()
    integer, dimension(:), pointer :: iat => null(), ja => null()

    ! setting default values
    max_iterations = size(b) * 10
    tolerance = 1e-5_dp
    if (present(tol)) tolerance = tol
    if (present(max_it)) max_iterations = max_it
    
    ! dimension check
    n = size(b)
    if (A%patt%nrows .ne. n) then
        print *, 'ERROR: dimension mismatch (in BiCGSTAB).'
        stop
    end if

    ! setting handles 
    coef => A%coef
    iat => A%patt%iat
    ja => A%patt%ja

    allocate(r(n), r0(n), r_new(n), p(n), s(n), A_p(n), A_s(n))

    ! Start of BiCGSTAB algorithm
    
    ! the initial guess is the zero-vector
    x = 0.0_dp
    
    ! r_0 = b - A x_0
    call amxpby(r, -1.0_dp, A, x, 1.0_dp, b)
    
    ! r_0^ arbitrary
    r0 = r
    
    ! p0 = r0
    p = r
    
    do j = 1,max_iterations
        if (norm(r) <= tolerance) exit
        
        ! alpha_j = (r_j,r_0^) / (A p_j, r_0^)
        r_r0 = inner_prod(r, r0)
        call reset(A_p)
        call amxpby(A_p, 1.0_dp, A, p)
        alpha = r_r0 / inner_prod(A_p, r0)

        ! s_j = r_j - alpha_j A p_j
        s = r - alpha * A_p
        
        ! omega_j = (A s_j, s_j) / (A s_j, A s_j)
        call reset(A_s)
        call amxpby(A_s, 1.0_dp, A, s)
        omega = inner_prod(A_s, s) / inner_prod(A_s, A_s)

        ! x_j+1 = x_j + alpha_j p_j + omega_j s_j
        call axpby(x, alpha, p, omega, s)
        
        ! r_j+1 = s_j - omega_j A s_j
        call reset(r_new)
        call axpby(r_new, -omega, A_s, 1.0_dp, s)
        
        ! beta_j = (r_j+1, r_0^) / (r_j, r_0^) * (alpha_j / omega_j)
        beta = inner_prod(r_new, r0) / r_r0 * (alpha / omega)

        ! p_j+1 = r_j+1 + beta_j (p_j - omega_j A p_j)
        p = r_new + beta * (p - omega * A_p)

        r = r_new
    end do
    
    print *, 'BiCGSTAB iterations:', j-1

    deallocate(r, r0, r_new, p, s, A_p, A_s)
end subroutine


end module