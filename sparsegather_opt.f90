subroutine sparsegather_opt(u, v)
   use header
   implicit none
   include "mpif.h"
   type(Vector), intent(inout) :: u, v
   integer :: myid, num_proc, ibeg, iend, ierr, m ! Processor id, number of processors, associated processor beginning and ending indices                                      
   real(kind = 8), allocatable, dimension(:) :: prev_proc_vals, next_proc_vals, send_forward_vals, send_backward_vals ! Arrays to be sent and received                              
   
   !---------------------------------------------------------------------------------------------------------------------------!
   ! Difference with initial sparsegather subroutine is that this one communicates u and v in a single vector across processors!
   !---------------------------------------------------------------------------------------------------------------------------!

   call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
   call mpi_comm_size(MPI_COMM_WORLD, num_proc, ierr)

   !-----------------------------------!                                                                                                                                            
   ! Allocate memory and assign values !                                                                                                                                            
   !-----------------------------------!                                                                                                                                           
   ibeg = u%ibeg
   iend = u%iend
   m = INT(SQRT(REAL(u%n)))
   allocate(prev_proc_vals(2*m))
   allocate(next_proc_vals(2*m))
   allocate(send_forward_vals(2*m))
   allocate(send_backward_vals(2*m))

   send_forward_vals(1:m) = u%xx(iend-m+1:iend)
   send_forward_vals(m+1:2*m) = v%xx(iend-m+1:iend)

   send_backward_vals(1:m) = u%xx(ibeg:ibeg+m-1)
   send_backward_vals(m+1:2*m) = v%xx(ibeg:ibeg+m-1)

   !----------------------------------------------------------------------------------------------------------------------------------!                                             
   ! An assumption with the following gathering implementation is that the number of rows each processor attends to is greater than m !                                             
   ! Hence we can check whether we are in the first or last processor to deal with the edge cases in the Laplace Matrix.               !                                            
   !----------------------------------------------------------------------------------------------------------------------------------!                                             

   ! Check whether on first processor                                                                                                                                               
   if (myid == 0) then

      ! Send iend to next processor (i.e. to ibeg-1 for the processor) and receive from last processor as we are at the boundary case                                              \
                                                                                                                                                                                    
      call MPI_Sendrecv(send_forward_vals, 2*m, MPI_DOUBLE_PRECISION, myid+1, 0, &
                  prev_proc_vals,2*m, MPI_DOUBLE_PRECISION, num_proc-1, 0, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      ! send ibeg to previous processor (i.e. to iend+1 for the prev processor) and receive from next proeccor                                                                     \
                                                                                                                                                                                    
      call MPI_Sendrecv(send_backward_vals,2*m, MPI_DOUBLE_PRECISION, num_proc-1, 0, &
                  next_proc_vals, 2*m, MPI_DOUBLE_PRECISION, myid+1, 0, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      ! Fill the vector with the additional values received                                                                                                                         
      u%xx(u%n - m+1:u%n) = prev_proc_vals(1:m)
      v%xx(v%n - m+1:v%n) = prev_proc_vals(m+1:2*m)

      u%xx(iend+1:iend+m) = next_proc_vals(1:m)
      v%xx(iend+1:iend+m) = next_proc_vals(m+1:2*m)

   ! Check whether on last processor                                                                                                                                                
   else if (myid == num_proc-1) then

      ! Send the forward values to Processor zero and receive from myid-1 processor            
      call MPI_Sendrecv(send_forward_vals, 2*m, MPI_DOUBLE_PRECISION, 0, 0, &
                  prev_proc_vals, 2*m, MPI_DOUBLE_PRECISION, myid-1, 0, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      ! Send backward values to procssor myid-1 and receive from processor zero                                                                                                     
      call MPI_Sendrecv(send_backward_vals, 2*m, MPI_DOUBLE_PRECISION, myid-1, 0, &
                  next_proc_vals, 2*m, MPI_DOUBLE_PRECISION, 0, 0, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      ! Fill new values received                                                                                                                                                    
      u%xx(ibeg - m: ibeg-1) = prev_proc_vals(1:m)
      v%xx(ibeg - m: ibeg-1) = prev_proc_vals(m+1:2*m)
      
      u%xx(1:m) = next_proc_vals(1:m)
      v%xx(1:m) = next_proc_vals(m+1:2*m)

   else

      ! Send forward values to processor myid+1 and receive from myid-1                                                                                                             
      call MPI_Sendrecv(send_forward_vals, 2*m, MPI_DOUBLE_PRECISION, myid+1, 0, &
                  prev_proc_vals, 2*m, MPI_DOUBLE_PRECISION, myid-1, 0, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      ! Send backward values to processor myid-1 and receive from myid+1                                                                                                            
      call MPI_Sendrecv(send_backward_vals, 2*m, MPI_DOUBLE_PRECISION, myid-1, 0, &
                  next_proc_vals, 2*m, MPI_DOUBLE_PRECISION, myid+1, 0, &
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      ! Fill new values received                                                                                                                                                    
      u%xx(ibeg-m:ibeg-1) = prev_proc_vals(1:m)
      v%xx(ibeg-m:ibeg-1) = prev_proc_vals(m+1: 2*m)

      u%xx(iend+1:iend+m) = next_proc_vals(1:m)
      v%xx(iend+1:iend+m) = next_proc_vals(m+1: 2*m)

   end if

   deallocate(prev_proc_vals)
   deallocate(next_proc_vals)
   deallocate(send_forward_vals)
   deallocate(send_backward_vals)

end subroutine
