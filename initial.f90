subroutine initial(m, u, v)
  !========================================
  ! A subroutine for assembling 2-Gaussians
  ! initial condition
  !========================================
  use header
  implicit none
  include "mpif.h"

  ! Variables as per assignment
  integer, intent(in) :: m
  type(Vector), intent(inout) :: u,v

  ! Loop counters: total, horizontal, vertical
  integer ::       irow,  i,          j
  ! mesh step, horizontal, vertical coordinates
  real(kind=rk) :: h, x, y

  ! 2*variance of Gaussians
  real(kind=rk), parameter :: v2 = 0.02_rk
  integer :: myid     ! processor rank
  integer :: numprocs ! number of processors                                                                         
  integer :: ierr
  integer :: rows_per_process, start_row, end_row

  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)

  ! Calculate starting and end indices for the current processor
  rows_per_process = (m*m) / numprocs
  start_row = myid*rows_per_process + 1
  end_row = (myid + 1)*rows_per_process
  h = 1.0_rk/m
      
  do irow=start_row,end_row
     ! Calculate Cartesian index splitting
     j = (irow-1)/m + 1
     i = irow - (j-1)*m
     y = j*h
     x = i*h

     u%xx(irow) = 1.0_rk - 0.5_rk*exp(-(x-0.5_rk)**2/v2 -(y-0.5_rk)**2/v2) &
                         - 0.5_rk*exp(-(x-0.4_rk)**2/v2 -(y-0.6_rk)**2/v2)
     v%xx(irow) = 0.25_rk*exp(-(x-0.4_rk)**2/v2 -(y-0.4_rk)**2/v2)  &
                + 0.25_rk*exp(-(x-0.5_rk)**2/v2 -(y-0.6_rk)**2/v2)     
  end do

end subroutine initial
     
