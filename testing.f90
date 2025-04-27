program testing
  ! ===================================                                                                                                                                            
  ! Testing Program
                                                                                                                                                           
  ! ===================================                                                                                                                                            

  use header
  implicit none
  include "mpif.h"

  ! Variables as per assignment                                                                                                                                                    
  type(Matrix) :: Delta
  type(Vector) :: u, v, ru 
  real(kind=rk) :: Du
  integer       :: myid, nprocs, nrows, ibeg, iend, rows_per_process, ierr, m

  ! 2D problem size                                                                                                                                                               
  integer :: n

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  
  m = 4

  if (m<=0) then
    call MPI_Finalize(ierr)
    stop
  end if
 
  n = m*m
  rows_per_process = n / nprocs
  ibeg = myid*rows_per_process + 1
  iend = (myid +1)*rows_per_process
  nrows = iend - ibeg + 1
  
  ! Allocate Memory to Delta portion on processor                                                                                                                                  
  allocate(Delta%ii(rows_per_process+1))
  allocate(Delta%aa(5*rows_per_process))
  allocate(Delta%jj(5*rows_per_process))

  Delta%n = n
  Delta%ibeg = ibeg
  Delta%iend = iend

  ! Assemble Laplace matrix for this processor                                                                                                                                     
  print*, "PROCESSOR ID: ", myid

  ! LAPLACE TESTING                                                                                                                                                                
  call create_matrix(4, Delta)
  print*, "Delta Values: ", Delta%aa(1:8)
  print*, "Delta Column indices: ", Delta%jj(1:8)
  print*, " Delta Row start index: ", Delta%ii(1:2)

  ! Initial guess                                                                                                                                                                  
  allocate(u%xx(n))
  allocate(v%xx(n))

  u%n = n
  u%ibeg = ibeg
  u%iend = iend

  v%n = n
  v%ibeg = ibeg
  v%iend = iend

  allocate(ru%xx(n))

  ru%n = n
  ru%ibeg = ibeg
  ru%iend = iend

  print*, "INIITAL TEST"
  call initial(m,u,v)

  print*,"ID ", myid,  " First 4 values of U", u%xx(1:4)
  print*,"ID ", myid,  " First 4 values of V", v%xx(1:4)

  print*, "MAT MULT TEST"
 
  ! Setting to easier values to calculate and compare
  Du = 1.0_8
  u%xx = 1.0_8
  call Mat_Mult(Delta, Du, u, ru)
  if (all(ru%xx ==0.0_8)) then
     print*, "Matrix Multiplcation test passed!"
  else
     print*, "Matrix Multiplcation test not passed!"
  end if
  
  ! Shut down MPI                                                                                                                                                                 
  call mpi_finalize(ierr)

end program testing
