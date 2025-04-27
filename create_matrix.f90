subroutine create_matrix(m, Delta)
  !===========================================
  ! A subroutine for assembling Laplace matrix
  ! in the compressed row storage format
  ! (Parallel version)
  ! Delta is a portion of the final Delta Matrix
  !===========================================
  use header
  implicit none
  include "mpif.h"

  ! Variables as per assignment
  integer, intent(in) :: m
  type(Matrix), intent(inout) :: Delta

  ! Loop counters: total, horizontal, vertical, nonzeros, neighbour
  integer ::       irow,  i,          j,        inz,      next
  ! Total problem size
  integer :: n
  
  integer :: myid     ! processor rank                                                                                                                                         
  integer :: numprocs ! number of processors       
  integer :: ierr
  integer :: rows_per_process, start_row, end_row

  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
! Assume Delta is properly declared prior to input  
  n = m*m

  rows_per_process = n / numprocs

  start_row = myid*rows_per_process + 1
  end_row = (myid + 1)*rows_per_process
  
  inz = 1  ! index of current nonzero element in processor

  ! Iterate over rows local to our process
  do irow = start_row, end_row

     ! Calculate Cartesian index splitting
     j = (irow-1)/m + 1
     i = irow - (j-1)*m

     ! Init current row position
     ! Delta%ii(irow) = inz
     Delta%ii(irow - start_row + 1) = inz
     
     ! We want Delta%jj to contain GLOBAL COLUMN POSITIONS
     ! Diagonal element
     Delta%aa(inz) = -4.0_rk * m**2 ! -4/h^2
     Delta%jj(inz) = irow           ! column index = row index on diagonal
     ! Done with diagonal element, shift the counter of elements
     inz = inz + 1

     ! Off-diagonal elements, including periodically looped ones
     ! left
     Delta%aa(inz) = m**2           ! 1/h^2
     next = i-1
     if (next<1) next = m           ! loop over if we are on a boundary
     Delta%jj(inz) = next + (j-1)*m ! column index back in global range
     ! shift the nnz counter
     inz = inz + 1
     
     ! right
     Delta%aa(inz) = m**2
     next = i+1
     if (next>m) next = 1
     Delta%jj(inz) = next + (j-1)*m
     inz = inz + 1
     
     ! bottom
     Delta%aa(inz) = m**2
     next = j-1
     if (next<1) next = m
     Delta%jj(inz) = i + (next-1)*m
     inz = inz + 1
     
     ! top
     Delta%aa(inz) = m**2
     next = j+1
     if (next>m) next = 1
     Delta%jj(inz) = i + (next-1)*m
     inz = inz + 1
  end do

  ! Finalise row positions
  Delta%ii(rows_per_process + 1) = inz
  Delta%nnz = inz-1

  ! Sanity check to see sizes of Delta
  if (myid == 1) then 
     print*, "myID: ", myid
     print*, "Delta%aa",size(Delta%aa)
     print*, "Delta%jj",size(Delta%jj)
     print*, "Delta%ii",size(Delta%ii)
  end if 

end subroutine create_matrix
