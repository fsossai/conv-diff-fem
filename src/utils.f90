module utils
    use class_precision
    
contains

    subroutine write_vec(name, x)
        integer :: n
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: x(:)
        
        n = size(x)
        open(unit=1, file=name)
        write(unit=1, fmt='(i10,e15.8)') (i, x(i), i=1,n)
        close(unit=1)
    end subroutine

end module utils