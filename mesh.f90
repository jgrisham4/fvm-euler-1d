!===========================================================
! This module contains the definition for the mesh class.
!===========================================================
module class_mesh
  implicit none
  private
  public :: mesh,generate,write_to_file

  !---------------------------------------------------------
  ! Mesh class definition
  !---------------------------------------------------------
  type mesh
    integer                       :: num_points,num_elements
    double precision              :: xmin,xmax
    double precision              :: dx
    double precision, allocatable :: x(:)
    double precision, allocatable :: xc(:)
  end type mesh

  contains

    ! Subroutine for generating the mesh
    subroutine generate(this,npts,x_min,x_max)
      implicit none
      type(mesh), intent(inout) :: this
      integer, intent(in) :: npts
      double precision, intent(in) :: x_min, x_max
      integer :: i, aerr

      ! Assigning members
      this%num_points = npts
      this%num_elements = npts-1
      this%xmin = x_min
      this%xmax = x_max

      ! Allocating memory for node coordinates
      allocate(this%x(this%num_points),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: Couldn't allocate memory for mesh."
        stop
      end if

      ! Allocating memory for cell centers
      allocate(this%xc(this%num_elements),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: Couldn't allocate memory for cell center coordinates."
        stop
      end if

      ! Determining step size
      this%dx = (this%xmax - this%xmin)/(dble(this%num_points)-1.0d0)

      ! Generating the grid
      do i=1,this%num_points
        this%x(i) = dble(i-1)*this%dx
      end do

      ! Computing coordinates of the cell centers
      this%xc = 0.5d0*(this%x(1:this%num_points-1)+this%x(2:this%num_points))

      print *, "Done generating mesh." 

    end subroutine generate
        
    ! Subroutine for writing the mesh to a Tecplot file
    subroutine write_to_file(this,file_name)
      implicit none
      type(mesh), intent(in) :: this
      character (len=*) :: file_name

      ! Declaring local variables
      integer :: j

      ! Opening file
      open(2,file=file_name)
      write(2,*) 'variables=x'
      write(2,*) 'zone i=',this%num_points
      do j=1,this%num_points
        write(2,*) this%x(j)
      end do 
      close(2)

    end subroutine write_to_file
end module class_mesh
