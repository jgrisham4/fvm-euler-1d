program test_mesh
  use class_mesh, only : mesh, generate, write_to_file
  implicit none
  double precision  :: x_max, x_min
  integer           :: npts
  character (len=30) :: mesh_file

  ! Creating mesh object
  type(mesh) :: grid

  ! Defining inputs
  x_max = 1.0d0
  x_min = 0.0d0
  npts  = 100
  mesh_file = "grid.dat"
  
  ! Generating mesh
  call generate(grid,npts,x_min,x_max)

  ! Writing mesh to file
  call write_to_file(grid,mesh_file)

end program 
