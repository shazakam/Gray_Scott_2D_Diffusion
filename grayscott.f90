program grayscott
  ! ===================================
  ! Main program for solving
  ! 2-dimensional Gray-Scott equations
  ! (Parallel version)
  ! ===================================
  
  use header
  implicit none
  include "mpif.h"

  ! Variables as per assignment
  type(Matrix) :: Delta
  type(Vector) :: u, v 
  real(kind=rk) :: tau, T, F, k, Du, Dv
  integer       :: myid, nprocs, nrows, ibeg, iend, rows_per_process, ierr, m
  real(kind = 8), allocatable, dimension(:) ::  total_u, total_v
  
  ! 2D problem size
  integer :: n

  ! Ticks for CPU timing
  real(kind=rk) :: t_start, t_finish

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  
  ! Read in the model parameters                                                                                                                                         
  if (myid == 0) then
     open(unit=10, file="input.dat", status="old", action="read")
     read(10,*) m
     read(10,*) tau
     read(10,*) T
     read(10,*) F
     read(10,*) k
     read(10,*) Du
     read(10,*) Dv
     close(10)
     print*, "m", m
     print*, "tau", tau
     print*, "T", T
     print*, "F", F
     print*, "k", k
     print*, "Du", Du
     print*, "Dv", Dv
     print*, ""
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                        
!      Broadcast parameters to the other processes                                                                                        
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                        
  
   call mpi_bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(F,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(k,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(Du,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(Dv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  
   if (m<=0) then
      call MPI_Finalize(ierr)
      stop
   end if

  ! Local range of the current process
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
  call create_matrix(m, Delta)

  ! Initial guess
  allocate(u%xx(n))
  allocate(v%xx(n))

  u%n = n
  u%ibeg = ibeg
  u%iend = iend

  v%n = n 
  v%ibeg = ibeg
  v%iend = iend
  
  ! Allocate memory for final solution vectors
  allocate(total_u(u%n))  
  allocate(total_v(v%n))  
 
  call initial(m,u,v)

  ! Run the Heun method - NOTE THAT THE FINAL SOLUTION IS ONLY WRITTEN TO PROCESSOR ZERO

  call cpu_time(t_start)
    call timestepping(u, v, tau, T, F, k, Du, Dv, Delta)
  call cpu_time(t_finish)

  ! Gather solutions onto root processor                                                                                                              
  call MPI_Gather(u%xx(u%ibeg:u%iend), iend-ibeg+1, MPI_DOUBLE_PRECISION, &
                total_u, iend-ibeg+1, MPI_DOUBLE_PRECISION, &
                0, MPI_COMM_WORLD, ierr)

  call MPI_Gather(v%xx(v%ibeg:v%iend), iend-ibeg+1, MPI_DOUBLE_PRECISION, &
                total_v, iend-ibeg+1, MPI_DOUBLE_PRECISION, &
                0, MPI_COMM_WORLD, ierr)


  if (myid == 0) then
     u%xx = total_u
     v%xx = total_v
  end if

  ! Sanity check
  if (all(u%xx == 0.0) .and. (myid == 1)) then
    print *, "FINAL U Array is all zero!"
  end if
  
  ! Print out solution and time from processor zero
  if (myid == 0) then
     write(*, '(A,F17.12)') 'u_{n/2,n}     = ', u%xx(m/2+(m-1)*m)
     write(*, '(A,F17.12)') 'v_{n/2,n}     = ', v%xx(m/2+(m-1)*m)
     write(*, '(A,F8.4)')  'total cpu time = ', t_finish - t_start

     print*, "PROCESSOR 0 is SAVING THE SOLUTION"
     call save_fields(u,v,'solution.dat')
  end if

  ! Free up memory
  deallocate(Delta%ii)
  deallocate(Delta%aa)
  deallocate(Delta%jj)
  deallocate(u%xx, v%xx)
  deallocate(total_u, total_v)

  ! Shut down MPI                                                                                                                  
  call mpi_finalize(ierr)

end program grayscott
