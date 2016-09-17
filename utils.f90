!===========================================================
! This module contains some utility functions
!
! NOTE: The assumption of a calorically perfect gas is 
!       made in the below conversions to and from 
!       conservative and primitive variable vectors.
!===========================================================

module utils
  implicit none
  private
  public :: w_to_u, u_to_w

  contains

    !------------------------------------------------------
    ! Function for converting a vector of primitive 
    ! variables to a vector of conserved variables.  
    ! Assuming a calorically perfect gas.
    !------------------------------------------------------
    function w_to_u(w,g) result(u)
      implicit none
      double precision, intent(in) :: w(3),g
      double precision             :: u(3)
      u(1) = w(1)
      u(2) = w(1)*w(2)
      u(3) = w(3)/(g-1.0d0)+0.5d0*w(1)*w(2)**2
    end function w_to_u
    
    !------------------------------------------------------
    ! Function for converting from vector of conserved 
    ! variables to a vector of primitive variables.
    ! Assuming a calorically perfect gas.
    !------------------------------------------------------
    function u_to_w(u,g) result(w)
      implicit none
      double precision, intent(in) :: u(3),g
      double precision             :: w(3)
      w(1) = u(1)
      w(2) = u(2)/u(1)
      w(3) = (g-1.0d0)*(u(3)-0.5d0*u(2)**2/u(1))
    end function u_to_w

end module utils
