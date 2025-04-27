subroutine RHS(u,v,F,k,Du,Dv,Delta,ru,rv)
  ! OPTIMISED IMPLEMENTATION Q10
  !==========================================
  ! Implement Gray-Scott right hand side here
  !==========================================
  use header
  implicit None
  include "mpif.h"

  type(Vector), intent(inout) :: u, v
  real(kind=8), intent(in) :: F, k, Du, Dv
  type(Matrix), intent(in) :: Delta
  type(Vector), intent(inout) :: ru, rv
  real(kind=8), allocatable, dimension(:) :: u_v_v
  integer :: nprocs, ierr

  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! Allocate memory for temporary variable
  allocate(u_v_v(u%n))

  if (nprocs > 1) then                                    
      call sparsegather_opt(u, v)
  end if       

 ! Call Parallel Matrix Multiplcation and overwrite to ru and rv
  call Mat_Mult(Delta, Du, u, ru)
  call Mat_Mult(Delta, Dv, v, rv)

  ! Calculate u*v*v once since it can be reused
  u_v_v(u%ibeg:u%iend) = u%xx(u%ibeg:u%iend) * v%xx(v%ibeg:v%iend) * v%xx(v%ibeg:v%iend)

  ! Up to here ru%xx is populated with Du*Delta*u and likewise for v
  ru%xx(ru%ibeg:ru%iend) = -u_v_v(u%ibeg:u%iend) + F - F*u%xx(u%ibeg:u%iend) + ru%xx(ru%ibeg:ru%iend) 
  rv%xx(v%ibeg:v%iend) = u_v_v(v%ibeg:v%iend) - (F + k)*v%xx(v%ibeg:v%iend) + rv%xx(rv%ibeg:rv%iend)
 
  ! Deallocate
  deallocate(u_v_v)

end subroutine RHS
