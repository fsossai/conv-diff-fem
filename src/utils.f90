module utils
    use class_precision
    use class_CSRMAT
    
contains

    subroutine write_vec(name, x)
        integer :: n
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: x(:)
        
        n = size(x)
        open(unit=1, file=name)
        write(unit=1, fmt='(i10,1x,1d15.8)') (i, x(i), i=1,n)
        close(unit=1)
    end subroutine

    subroutine write_csrmat(name, A)
        character(len=*), intent(in) :: name
        type(CSRMAT), intent(in) :: A
        integer, pointer :: iat(:), ja(:)
        real(dp), pointer :: coef(:)
        character(len=19) :: frmt = '(i7,1x,i7,1x,1d15.8)'
        integer :: n,i,j,c_start,c_end

        ! creating handles
        iat => A%patt%iat
        ja => A%patt%ja
        coef => A%coef
        n = A%patt%nrows

        open(unit=1, file=name)
        write(1,*) n,A%patt%nterm
        do i=1,n
            c_start = iat(i)
            c_end = iat(i+1) - 1
            do j = c_start,c_end
                write(1,frmt) i,ja(j),coef(j)
            end do
        end do
        
        close(unit=1)
    end subroutine

    subroutine print_vec_compact(x, cols)
        real(dp), intent(in) :: x(:)
        integer, optional, intent(in) :: cols
        integer :: i,j,n,c
        n = size(x)
        
        c = 5 ! default value
        if (present(cols)) c = cols

        do i = 1,n,c
            do j = 0,c-1
                if (i+j <= n) write(*, '(1f15.8x)', advance='no') x(i+j)
            end do
            print *
        end do
    end subroutine
end module utils