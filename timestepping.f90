subroutine timestepping(u,v,tau,T,F,k,Du,Dv,Delta)
  !=========================================
  ! Implement Heun method (Algorithm 1) here
  !=========================================
  use header
  implicit None
  include "mpif.h"

  type(Vector), intent(inout) :: u, v
  type(Matrix), intent(in) :: Delta
  real(kind = 8), intent(in) :: tau, T, F, k, Du, Dv
  integer :: num_steps, i, myid, ierr 
  type(Vector) :: ru, rv, ru_hat, rv_hat

  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  num_steps = INT(T/tau)

  !---------------------------------------!
  ! Allocate and assigning required values!
  !---------------------------------------!

  ru%n = u%n
  rv%n = v%n

  ru%ibeg = u%ibeg
  ru%iend = u%iend

  rv%ibeg = v%ibeg
  rv%iend = v%iend

  allocate(ru%xx(ru%n))
  allocate(rv%xx(rv%n))

  ru_hat%n = u%n
  rv_hat%n = v%n

  ru_hat%ibeg = u%ibeg
  ru_hat%iend = u%iend

  rv_hat%ibeg = v%ibeg
  rv_hat%iend = v%iend

  allocate(ru_hat%xx(ru_hat%n))
  allocate(rv_hat%xx(rv_hat%n))
  
  ! ======= Sanity check on allocated vectors =======!
  if (.not. allocated(ru%xx)) print*, "ID:", myid, "ru not allocated!"
  if (.not. allocated(rv%xx)) print*, "ID:", myid, "rv not allocated!"
  if (.not. allocated(ru_hat%xx)) print*, "ID:", myid, "ru_hat not allocated!"
  if (.not. allocated(rv_hat%xx)) print*, "ID:", myid, "rv_hat not allocated!"
  
  print*, "ID: ", myid, " About to start loop"
  print*,"ID: ", myid, "NUMBER OF STEPS: ",  num_steps


  if (all(u%xx(u%ibeg:u%iend) == 0.0)) then
    print*,"ID ", myid, " u is zero prior to LOOP"
  end if

  ! ======================================== !

  !----- START OF TIME STEPPING ------!
  do i = 1, num_steps
     
     ! Calculate ru and rv
     call RHS(u,v,F,k,Du,Dv,Delta,ru,rv)
     
     ! Update to get u_hat and v_hat in algorithm
     u%xx(u%ibeg:u%iend) = u%xx(u%ibeg:u%iend) + tau*ru%xx(ru%ibeg:ru%iend)
     v%xx(v%ibeg:v%iend) = v%xx(v%ibeg:v%iend) + tau*rv%xx(rv%ibeg:rv%iend)

     ! Calculate ru_hat and rv_hat
     call RHS(u,v,F,k,Du,Dv,Delta,ru_hat,rv_hat)
     
     ! Final Correction Update 
     u%xx(u%ibeg:u%iend) = u%xx(u%ibeg:u%iend) + tau*(ru_hat%xx(ru_hat%ibeg:ru_hat%iend) - ru%xx(ru%ibeg:ru%iend))*0.5_8
     v%xx(v%ibeg:v%iend) = v%xx(v%ibeg:v%iend) + tau*(rv_hat%xx(rv_hat%ibeg:rv_hat%iend) - rv%xx(rv%ibeg:rv%iend))*0.5_8
  end do

end subroutine timestepping

