program test_solver
  use class_mesh,   only : mesh, generate, write_to_file
  use class_solver, only : solver, initialize, solve, write_results
  implicit none
  double precision              :: x_max, x_min
  integer                       :: npts
  character (len=30)            :: soln_file
  double precision, allocatable :: ic(:,:)
  integer                       :: ntsteps,i
  double precision              :: dt,gam

  ! Creating mesh object
  type(mesh)   :: grid
  type(solver) :: euler_solver

  ! Defining inputs
  x_max     = 1.0d0
  x_min     = 0.0d0
  npts      = 201
  ntsteps   = 800
  dt        = 0.0005d0
  gam       = 1.4d0
  soln_file = "solution.dat"
  
  ! Generating mesh
  call generate(grid,npts,x_min,x_max)

  ! Creating initial condition -- including ghost cells
  allocate(ic(3,0:grid%num_elements+1))
!  do i=0,grid%num_elements+1
!    if (grid%xc(i).le.2.0d0) then
!      ic(1,i) = 10.0d0   ! Density
!      ic(2,i) = 0.0d0    ! Velocity
!      ic(3,i) = 100.0d0  ! Pressure
!    else
!      ic(1,i) = 1.0d0    ! Density
!      ic(2,i) = 0.0d0    ! Velocity
!      ic(3,i) = 1.0d0    ! Pressure
!    end if
!  end do
  do i=0,grid%num_elements+1
    if (grid%xc(i).le.0.5d0) then
      ic(1,i) = 1.0d0    ! Density
      ic(2,i) = 0.0d0    ! Velocity
      ic(3,i) = 1.0d0    ! Pressure
    else
      ic(1,i) = 0.125d0  ! Density
      ic(2,i) = 0.0d0    ! Velocity
      ic(3,i) = 0.1d0    ! Pressure
    end if
  end do

  ! Initializing solver
  call initialize(euler_solver,grid,ntsteps,dt,gam,ic)
  
  ! Solving the problem
  call solve(euler_solver)

  ! Writing results to file
  call write_results(euler_solver,soln_file)

end program test_solver
