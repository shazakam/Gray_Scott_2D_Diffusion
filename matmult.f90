!------------------------------------------------------------------------
!
!  Subroutine to add a matrix by vector product in compressed row storage 
!  format, i.e to calculate b = b + alpha * A*u (Parallel version).
!
!------------------------------------------------------------------------
! OPTIMISED IMPLEMENTATION Q10
subroutine Mat_Mult(A,alpha,u,b)

  use header
  implicit None
  include "mpif.h"

  type(Matrix), intent(in) :: A     ! Matrix to multiply with
  real(kind=rk), intent(in):: alpha ! Multiple of matrix
  type(Vector), intent(in) :: u     ! Input vector u
  type(Vector), intent(inout) :: b  ! Output vector b = b + alpha*A*u
  real(kind=rk) :: Aurow            ! temp variable for (Au)(i)
  integer :: i, i_j, j

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate each component of b by taking the scalar product of the
!     i-th row of A and the vector u.
!
!       Note that in the compressed row storage format the nonzero 
!       entries of row i are stored in 
!
!         A%aa(A%ii(i)), A%aa(A%ii(i)+1), ..., A%aa(A%ii(i+1)-1)
!
!       the according (global) column numbers are stored in
!
!         A%jj(A%ii(i)), A%jj(A%ii(i)+1), ..., A%jj(A%ii(i+1)-1)   
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!-------------------------------------------------------------------------------!
! Begin Matrix Multiplication with gathered vector values and stored CRS Matrix !
!-------------------------------------------------------------------------------!

  do i = b%ibeg,b%iend
     Aurow = 0.0_rk
     b%xx(i) = 0.0_8

     do i_j = A%ii(i-b%ibeg+1), A%ii(i-b%ibeg+2)-1
        j = A%jj(i_j)
        Aurow = Aurow + A%aa(i_j) * u%xx(j)
     end do

     b%xx(i) = b%xx(i) + alpha*Aurow
  end do

end subroutine Mat_Mult
