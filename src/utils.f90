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
        write(unit=1, fmt='(i10,1x,1e15.8)') (i, x(i), i=1,n)
        close(unit=1)
    end subroutine

    subroutine write_csrmat(name, A)
        character(len=*), intent(in) :: name
        type(CSRMAT), intent(in) :: A
        integer, pointer :: iat(:), ja(:)
        real(dp), pointer :: coef(:)
        character(len=20) :: frmt = '(i7,1x,i7,1x,1e15.8)'
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
                if (i+j <= n) write(*, '(e15.8,1x)', advance='no') x(i+j)
            end do
            print *
        end do
    end subroutine

    subroutine get_boundaries(coord, tolerance, bnodes)
        ! Given a triangulation, finds all the nodes lying on the borders 
        ! of a square defined by the following two corners
        ! [min(x), min(y)] [max(x), max(y)]
        ! called point_a and point_b respectively

        real(dp), intent(in)                :: coord(:,:)
        real(dp), intent(in)                :: tolerance
        integer, allocatable, intent(out)   :: bnodes(:)
        
        logical, allocatable    :: otb(:) ! stands for 'On The Border'
        integer                 :: i, j, nn, notb
        real(dp), dimension(2)  :: point_a, point_b

        nn = size(coord, 1)
        allocate(otb(nn))
        otb = .false.

        ! Finding point_a and point_b (x,y) coordinates
        point_a = minval(coord, 1)
        point_b = maxval(coord, 1)

        ! Scanning x and y separately for cache-friendliness
        
        ! scanning x 
        do i = 1, nn
            if (abs(coord(i,1) - point_a(1)) .le. tolerance .or. &
                abs(coord(i,1) - point_b(1)) .le. tolerance) then
                otb(i) = .true.
            end if
        end do

        ! scanning y
        do i = 1, nn
            if (abs(coord(i,2) - point_a(2)) .le. tolerance .or. &
                abs(coord(i,2) - point_b(2)) .le. tolerance) then
                otb(i) = .true.
            end if
        end do

        notb = count(otb)
        allocate(bnodes(notb))
        
        j = 0
        do i = 1, nn
            if (otb(i)) then
                j = j + 1
                bnodes(j) = i
            end if
        end do

        print *, 'j',j,'notb',notb
        if (j .ne. notb) stop 'ERROR: corrupted'

        end subroutine
end module utils