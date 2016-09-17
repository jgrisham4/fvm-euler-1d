!===========================================================
! This module contains different limiters for use in the 
! FVM solver
!===========================================================

module limiters
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: minmod, vanleer

  contains

    !------------------------------------------------------
    ! minmod limiter
    !------------------------------------------------------
    double precision function minmod(r) result(phi)
      implicit none
      double precision, intent(in) :: r
      if (r.le.0.0d0) then
        phi = 0.0d0
      else if (r.lt.1.0d0) then
        phi = r
      else 
        phi = 1.0d0
      end if 
    end function minmod

    !------------------------------------------------------
    ! van Leer limiter
    !------------------------------------------------------
    double precision function vanleer(r) result(phi)
      implicit none
      double precision, intent(in) :: r
      if (r.le.0.0d0) then
        phi = 0.0d0
      else if (isnan(r).or..not.ieee_is_finite(r)) then
        phi = 2.0d0
      else 
        phi = 2.0d0*r/(r+1.0d0)
      end if 
    end function vanleer

end module limiters
