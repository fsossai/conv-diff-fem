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
        write(unit=1, fmt='(i10,1x,e15.8)') (i, x(i), i=1,n)
        close(unit=1)
    end subroutine

    subroutine write_csrmat(name, A)
        character(len=*), intent(in) :: name
        type(CSRMAT), intent(in) :: A
        integer, pointer :: iat(:), ja(:)
        real(dp), pointer :: coef(:)
        character(len=19) :: frmt = '(i7,1x,i7,1x,e15.8)'
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

end module utils