! Bi-conjugate Gradient Stabilized

! subroutine axbnsy(n,ntot,nterm,iat,ja,coef_A,xvec,bvec)
!subroutine bicgstab(n, ntot, nterm, iat, ja, coef, xvec, 
subroutine bicgstab(A, b, x)
    use class_precision
    use class_CSRMAT

    type(CSRMAT), intent(in)            :: A
    real(dp), intent(in)                :: b(:)
    real(dp), allocatable, intent(out)  :: x(:)

!    allocate(x(1))
    
end subroutine