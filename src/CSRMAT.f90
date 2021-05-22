!*****************************************************************************************
!
!  MODULE: Pattern
!
!> @brief Data type to store the pattern of a CSR matrix
!
!> @details
!! This module is used to create, store and destroy the non-zero pattern of a matrix
!! stored in the Compact Sparse Row format.
!
!> @author Carlo Janna
!
!> @version 1.0
!
!> @date January 2013
!
!> @par License:
!! This program is intended for internal research only and can not be distributed
!! elsewhere without authors' consent.
!*****************************************************************************************
module class_Pattern

use class_precision

implicit none

! Public member function/subroutines
public  :: new_Pattern,dlt_Pattern,wrPATT_unit

!> @brief Pattern of a matrix in CSR format
type Pattern
   !> @brief number of rows
   integer :: nrows = 0
   !> @brief number of non-zeroes
   integer :: nterm = 0
   !> @brief pointers to the beginning of each row
   integer, pointer :: iat(:) => null()
   !> @brief column indices of each row
   integer, pointer :: ja(:) => null()
end type

contains

   !**************************************************************************************
   ! FUNCTION: new_Pattern
   !
   !> @brief Allocates the Pattern data structure
   !
   !> @author Carlo Janna
   !
   !> @param[in]    nrows # of rows
   !! @param[in]    nterm # of non zeroes
   !! @param[inout] Pattern_inout Pattern variable
   !! @param[out]   ierr error code
   !! @param        == 0, successful allocation
   !! @param        /= 0, allocation error
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   function new_Pattern(nrows,nterm,Pattern_inout) result(ierr)

   implicit none

   ! Input variables
   integer, intent(in)          :: nrows,nterm
   ! Input/Output variables
   type(Pattern), intent(inout) :: Pattern_inout
   ! Output variables
   integer                      :: ierr
   ! Local variables
   integer                      :: info

   ! Initialize error code
   ierr = 0

   Pattern_inout%nrows = nrows
   Pattern_inout%nterm = nterm
   if ( associated(Pattern_inout%iat) ) then
      deallocate(Pattern_inout%iat,stat=info)
      if (info .ne. 0) ierr = 1
   endif
   if ( associated(Pattern_inout%ja) ) then
      deallocate(Pattern_inout%ja,stat=info)
      if (info .ne. 0) ierr = 1
   endif
   if (nrows.gt.0) allocate(Pattern_inout%iat(nrows+1),stat=info)
   if (info .ne. 0) ierr = 1
   if (nterm.gt.0) allocate(Pattern_inout%ja(nterm),stat=info)
   if (info .ne. 0) ierr = 1

   end function new_Pattern

   !**************************************************************************************
   ! FUNCTION: dlt_Pattern
   !
   !> @brief deallocates the Pattern data structure
   !
   !> @author Carlo Janna
   !
   !! @param[inout] Pattern_inout Pattern variable
   !! @param[out]   ierr error code
   !! @param        == 0, successful allocation
   !! @param        /= 0, allocation error
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   function dlt_Pattern(Pattern_inout) result(ierr)

   implicit none

   ! Input/Output variables
   type(Pattern), intent(inout) :: Pattern_inout
   ! Output variables
   integer                      :: ierr
   ! Local variables
   integer                      :: info

   ! Initialize error code
   ierr = 0

   Pattern_inout%nrows = 0
   Pattern_inout%nterm = 0
   if ( associated(Pattern_inout%ja) ) then
      deallocate(Pattern_inout%ja,stat=info)
      if (info .ne. 0) ierr = 1
      Pattern_inout%ja => null()
   endif
   if ( associated(Pattern_inout%iat) ) then
      deallocate(Pattern_inout%iat,stat=info)
      if (info .ne. 0) ierr = 1
      Pattern_inout%iat => null()
   endif

   end function dlt_Pattern

   !-----------------------------------------------------------------------------------------
   !  SUBROUTINE: wrPATT_unit
   !
   !> @brief Write a matrix PATTERN to a given unit.
   !
   !> @details
   !! Write a matrix PATTERN to a given unit.
   !
   !  Parameters
   !
   !! @param[in]    ounit Output unit
   !! @param[in]    PATT Input pattern
   !
   !> @author Carlo Janna
   !
   !> @version 1.0
   !
   !> @date January 2013
   !
   !> @par License:
   !! This program is intended for internal research only and can not be distributed
   !! elsewhere without authors' consent.
   !-----------------------------------------------------------------------------------------

   subroutine wrPATT_unit(ounit,PATT)

   implicit none

   ! Input variables
   integer, intent(in)           :: ounit
   type(Pattern), intent(in)     :: PATT

   ! Local variables
   integer                       :: i,j,mm,nn
   integer                       :: nrows,nterm
   integer, pointer              :: iat(:),ja(:)
   character*16, parameter       :: frmt='(i7,1x,i7,1x,i1)'

   !-----------------------------------------------------------------------------------------

   ! Set handles
   nrows = PATT%nrows
   nterm = PATT%nterm
   iat => PATT%iat
   ja => PATT%ja
   ! Print matrix
   do i=1,nrows
      mm = iat(i)
      nn = iat(i+1) - 1
      do j=mm,nn
         write(ounit,frmt) i,ja(j),1
      end do
   end do
   flush(ounit)

   !-----------------------------------------------------------------------------------------

   end subroutine wrPATT_unit

end module class_Pattern

!*****************************************************************************************
!
!  MODULE: CSRMAT
!
!> @brief Data type to store a CSR matrix
!
!> @details
!! This module is used to define the Compact Sparse Row (CSR) data type for matrices
!! storage (see Y. Saad -- Iterative Methods for Sparse Linear Systems)
!
!> @author Carlo Janna
!
!> @version 1.0
!
!> @date January 2013
!
!> @par License:
!! This program is intended for internal research only and can not be distributed
!! elsewhere without authors' consent.
!*****************************************************************************************
module class_CSRMAT

use class_Pattern

implicit none

! Public member function/subroutines
public  :: new_CSRMAT,dlt_CSRMAT,new_CSRCOEF,dlt_CSRCOEF,copy_CSRMAT,wrCSR_unit
public  :: errchk_CSRMAT

!> @brief matrix stored in CSR format
type CSRMAT
   !> @brief matrix pattern
   type(Pattern)              :: patt
   !> @brief matrix coefficients
   real(kind=dp), pointer :: coef(:) => null()
end type

interface assignment (=)
   module procedure copy_CSRMAT
end interface

contains

   !**************************************************************************************
   ! FUNCTION: new_CSRCOEF
   !
   !> @brief Allocates the CSR matrix coefficients
   !
   !> @author Carlo Janna
   !
   !! @param[inout] mat_inout CSRMAT variable
   !! @param[out]   ierr error code
   !! @param        == 0, successful allocation
   !! @param        /= 0, allocation error
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   function new_CSRCOEF(CSRMAT_inout) result(ierr)

   implicit none

   ! Input/Output variables
   type(CSRMAT), intent(inout) :: CSRMAT_inout
   ! Output variables
   integer                     :: ierr
   ! Local variables
   integer                     :: nterm,info

   ierr = 0
   info = 0
   nterm = CSRMAT_inout%patt%nterm
   if ( associated( CSRMAT_inout%coef ) ) then
      deallocate(CSRMAT_inout%coef,stat=info)
      if (info .ne. 0) ierr = 1
   endif
   if (nterm .gt. 0) allocate(CSRMAT_inout%coef(nterm),stat=info)
   if (info .ne. 0) ierr = 1

   end function new_CSRCOEF

   !**************************************************************************************
   ! FUNCTION: dlt_CSRCOEF
   !
   !> @brief Deallocates the CSR matrix coefficients
   !
   !> @author Carlo Janna
   !
   !! @param[inout] mat_inout CSRMAT variable
   !! @param[out]   ierr error code
   !! @param        == 0, successful allocation
   !! @param        /= 0, allocation error
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   function dlt_CSRCOEF(CSRMAT_inout) result(ierr)

   implicit none

   ! Input/Output variables
   type(CSRMAT), intent(inout) :: CSRMAT_inout
   ! Output variables
   integer                     :: ierr

   ! Initialize error code
   ierr = 0

   if ( associated( CSRMAT_inout%coef ) ) then
      deallocate(CSRMAT_inout%coef,stat=ierr)
      CSRMAT_inout%coef => null()
   endif

   end function dlt_CSRCOEF

   !**************************************************************************************
   ! FUNCTION: new_CSRMAT
   !
   !> @brief Allocates the CSR matrix data structure
   !
   !> @author Carlo Janna
   !
   !> @param[in]    nrows # of rows
   !> @param[in]    nterm # of non zeroes
   !! @param[inout] mat_inout CSRMAT variable
   !! @param[out]   ierr error code
   !! @param        == 0, successful allocation
   !! @param        /= 0, allocation error
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   function new_CSRMAT(nrows,nterm,CSRMAT_inout) result(ierr)

   implicit none

   ! Input variables
   integer, intent(in)         :: nrows,nterm
   ! Input/Output variables
   type(CSRMAT), intent(inout) :: CSRMAT_inout
   ! Output variables
   integer                     :: ierr
   ! Local variables 
   integer                     :: info

   ! Initialize the error code
   ierr = 0

   info = new_Pattern(nrows,nterm,CSRMAT_inout%Patt)
   if (info .ne. 0) ierr = 1
   info = new_CSRCOEF(CSRMAT_inout)
   if (info .ne. 0) ierr = 1

   end function new_CSRMAT

   !**************************************************************************************
   ! FUNCTION: dlt_CSRMAT
   !
   !> @brief Deallocates the CSR matrix data structure
   !
   !> @author Carlo Janna
   !
   !! @param[inout] mat_inout CSRMAT variable
   !! @param[out]   ierr error code
   !! @param        == 0, successful allocation
   !! @param        /= 0, allocation error
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   function dlt_CSRMAT(CSRMAT_inout) result(ierr)

   implicit none

   ! Input/Output variables
   type(CSRMAT), intent(inout) :: CSRMAT_inout
   ! Output variables
   integer                     :: ierr
   ! Local variables 
   integer                     :: info

   ! Initialize the error code
   ierr = 0

   info = dlt_CSRCOEF(CSRMAT_inout)
   if (info .ne. 0) ierr = 1
   info = dlt_Pattern(CSRMAT_inout%Patt)
   if (info .ne. 0) ierr = 1

   end function dlt_CSRMAT

   !**************************************************************************************
   ! SUBROUTINE: copy_CSRMAT
   !
   !> @brief Copy a CSR matrix data structure
   !
   !> @author Carlo Janna
   !
   !> @param[in]    mat_in CSRMAT variable
   !> @param[inout] mat_out CSRMAT variable
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   subroutine copy_CSRMAT(mat_out,mat_in)

   implicit none

   ! Input variables
   type(CSRMAT), intent(in)    :: mat_in
   ! Output variables
   type(CSRMAT), intent(inout) :: mat_out
   ! Local variables
   integer                     :: nn,nt,ierr

   ! Get current size
   nn = mat_in%patt%nrows
   nt = mat_in%patt%nterm

   ! Allocate space for mat_out
   ierr = new_CSRMAT(nn,nt,mat_out)
   if (ierr .ne. 0) stop 'ERROR during copy_CSRMAT!'

   ! Copy the content of mat_in in mat_out
   if (nn .gt. 0) call SCOPY(nn+1,mat_in%patt%iat,1,mat_out%patt%iat,1)
   if (nt .gt. 0) call SCOPY(nt,mat_in%patt%ja,1,mat_out%patt%ja,1)
   if (nt .gt. 0) call DCOPY(nt,mat_in%coef,1,mat_out%coef,1)

   end subroutine copy_CSRMAT

   !-----------------------------------------------------------------------------------------
   !  SUBROUTINE: wrCSR_unit
   !
   !> @brief Write a CSR matrix to a given unit.
   !
   !> @details
   !! Write a CSR matrix to a given unit.
   !
   !  Parameters
   !
   !! @param[in]    ounit Output unit
   !! @param[in]    CSR_mat Input matrix
   !
   !> @author Carlo Janna
   !
   !> @version 1.0
   !
   !> @date January 2013
   !
   !> @par License:
   !! This program is intended for internal research only and can not be distributed
   !! elsewhere without authors' consent.
   !-----------------------------------------------------------------------------------------

   subroutine wrCSR_unit(ounit,CSR_mat)

   implicit none

   ! Input variables
   integer, intent(in)           :: ounit
   type(CSRMAT), intent(in)      :: CSR_mat

   ! Local variables
   integer                       :: i,j,mm,nn
   integer                       :: nrows,nterm
   integer, pointer              :: iat(:),ja(:)
   real(kind=dp), pointer    :: coef(:)
   character*20, parameter       :: frmt='(i7,1x,i7,1x,e20.10)'

   !-----------------------------------------------------------------------------------------

   ! Set handles
   nrows = CSR_mat%patt%nrows
   nterm = CSR_mat%patt%nterm
   iat => CSR_mat%patt%iat
   ja => CSR_mat%patt%ja
   coef => CSR_mat%coef
   ! Print matrix
   do i=1,nrows
      mm = iat(i)
      nn = iat(i+1) - 1
      do j=mm,nn
         write(ounit,frmt) i,ja(j),coef(j)
      end do
   end do
   flush(ounit)

   !-----------------------------------------------------------------------------------------

   end subroutine wrCSR_unit

   !**************************************************************************************
   ! FUNCTION: errchk_CSRMAT
   !
   !> @brief Error handling routine
   !
   !> @author Carlo Janna
   !
   !! @param[in] ounit output unit
   !! @param[in] sub string for case selection
   !! @param[in] ierr error to interpret
   !
   !> @version 1.0
   !
   !> @date January 2013
   !**************************************************************************************
   subroutine errchk_CSRMAT(ounit,sub,ierr)

   implicit none 
! Input variables
   character*(*), intent(in)    :: sub
   integer, intent(in)          :: ierr,ounit
! Local variables
   integer, parameter           :: shift=(ichar('a')-ichar('A'))
   integer, parameter           :: indA=ichar('A'),indZ=ichar('Z')
   integer                      :: i,ind,lenstr
   character*100                :: string

   lenstr = len(sub)
   do i = 1,lenstr
      ind = ichar(sub(i:i))
      if (ind.le.indZ .and. ind.ge.indA) then
         string(i:i) = char(ind+shift) 
      else
         string(i:i) = sub(i:i)
      endif
   enddo
   select case(string(1:lenstr))
      case('new_pattern')
         if (ierr.ne.0) then 
            write(ounit,'(a)') 'Pattern data type not allocated in new_Pattern'
         endif
      case('dlt_pattern')
         if (ierr.ne.0) then
            write(ounit,'(a)') 'Pattern data type not deallocated in dlt_Pattern'
         endif
      case('new_csrcoef')
         if (ierr.ne.0) then 
            write(ounit,'(a)') 'CSRCOEF not allocated in new_CSRCOEF'
         endif
      case('dlt_csrcoef')
         if (ierr.ne.0) then
            write(ounit,'(a)') 'CSRCOEF not deallocated in dlt_CSRCOEF'
         endif

      case('new_csrmat')
         if (ierr.ne.0) then 
            write(ounit,'(a)') 'CSRMAT data type not allocated in new_CSRMAT'
         endif
      case('dlt_csrmat')
         if (ierr.ne.0) then
            write(ounit,'(a)') 'CSRMAT data type not deallocated in dlt_CSRMAT'
         endif
      case default
         write(ounit,'(3a)') 'check string ''',sub,''' does not exist'
   end select

   end subroutine errchk_CSRMAT

end module class_CSRMAT
