module utils
    use class_precision
    use class_CSRMAT
    
contains

subroutine write_vec(name, x)
    character(len=*), intent(in)    :: name
    real(dp), intent(in)            :: x(:)
    integer                         :: n
    
    n = size(x)
    open(unit=1, file=name)
    write(unit=1, fmt='(i10,1x,f30.10)') (i, x(i), i=1,n)
    close(unit=1)
end subroutine


subroutine write_CSRMAT(name, A)
    character(len=*), intent(in)    :: name
    type(CSRMAT), intent(in)        :: A

    integer, pointer                :: iat(:), ja(:)
    real(dp), pointer               :: coef(:)
    character(len=20)               :: frmt = '(i7,1x,i7,1x,1e15.8)'
    integer                         :: n, i, j, c_start, c_end

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
    real(dp), intent(in)            :: x(:)
    integer, optional, intent(in)   :: cols
    integer                         :: i,j,n,c
    n = size(x)
    
    c = 5 ! default value
    if (present(cols)) c = cols

    do i = 1,n,c
        do j = 0,c-1
            if (i+j <= n) write(*, '(f15.8,1x)', advance='no') x(i+j)
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

    deallocate(otb)

end subroutine


subroutine sort2(x)
    integer, intent(inout) :: x(2)
    integer :: temp

    if (x(1) .gt. x(2)) then
        temp = x(1)
        x(1) = x(2)
        x(2) = temp
    end if
end subroutine


subroutine sort3(x)
    integer, intent(inout) :: x(3)

    call sort2(x(1:2))
    call sort2(x(2:3))
    call sort2(x(1:2))
end subroutine


function index_of_first(x, val) result(index)
    integer, intent(in) :: x(:), val
    integer :: i, n, index
    n = size(x, 1)

    do i = 1, n
        if (x(i) .eq. val) then
            index = i
            return
        end if
    end do

    index = -1
end function


subroutine print_int_array_2d(x, frmt)
    integer, intent(in)             :: x(:,:)
    character(len=*), intent(in)    :: frmt
    integer                         :: i, j

    do i = 1, size(x, 1)
        do j = 1, size(x, 2)
            write(*, fmt=frmt, advance='no') x(i, j)
        end do
        print *
    end do
end subroutine


subroutine print_real_array_2d(x, frmt)
    real(dp), intent(in)            :: x(:,:)
    character(len=*), intent(in)    :: frmt
    integer                         :: i, j

    do i = 1, size(x, 1)
        do j = 1, size(x, 2)
            write(*, fmt=frmt, advance='no') x(i, j)
        end do
        print *
    end do
end subroutine


subroutine set_insert(x, val)
    ! Replaces the first occurence of -1 in x only if val
    ! is not in x

    integer, intent(inout)  :: x(:)
    integer, intent(in)     :: val

    integer                 :: i, n

    n = size(x, 1)

    do i = 1, n
        if (x(i) .eq. val) return
        if (x(i) .eq. -1) then
            x(i) = val
            return
        end if
    end do
end subroutine  


recursive subroutine quicksort(a, first, last)
    implicit none
    integer, intent(inout) :: a(:)
    integer :: first, last
    integer :: i, j, x, t

    x = a( (first + last) / 2 )
    i = first
    j = last
    do
        do while (a(i) < x)
            i = i + 1
        end do
        do while (x < a(j))
            j = j - 1
        end do
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        i = i + 1
        j = j - 1
    end do
    if (first < i - 1) call quicksort(a, first, i - 1)
    if (j + 1 < last)  call quicksort(a, j + 1, last)
end subroutine quicksort


recursive subroutine paired_quicksort_abs(a, b, first, last)
    implicit none
    integer, intent(inout)  :: a(:), b(:)
    integer                 :: first, last
    integer                 :: i, j, x, t
    
    x = abs(a( (first + last) / 2 ))
    i = first
    j = last
    do
        do while (abs(a(i)) < x)
            i = i + 1
        end do
        do while (x < abs(a(j)))
            j = j - 1
        end do
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        t = b(i);  b(i) = b(j);  b(j) = t
        i = i + 1
        j = j - 1
    end do
    if (first < i - 1) call paired_quicksort_abs(a, b, first, i - 1)
    if (j + 1 < last)  call paired_quicksort_abs(a, b, j + 1, last)
end subroutine paired_quicksort_abs

end module utils