!===========================================================
! This module contains functions for solving the Riemann 
! problem.  It also contains a function which is used to 
! compute the actual value of the flux when provided with 
! the vector of primitive variables.
!===========================================================
module flux
  use utils, only : w_to_u
  implicit none
  private
  public :: roe_flux, euler_flux, roe_flux2

  contains

    !------------------------------------------------------
    ! Function for computing the flux term for the Euler
    ! equations.
    !
    ! parameters:
    ! - w: vector of primitive variables.
    ! - g: ratio of specific heats
    !------------------------------------------------------
    function euler_flux(w,g) result(f)
      implicit none
      double precision, intent(in) :: w(3),g
      double precision             :: f(3),u(3)
      u = w_to_u(w,g)
      f(1) = w(1)*w(2)
      f(2) = w(1)*w(2)**2 + w(3)
      f(3) = w(2)*(u(3)+w(3))
    end function euler_flux

    !------------------------------------------------------
    ! Function for approximately solving the Riemann 
    ! problem using Roe's solver.  This works by linearizing
    ! the flux Jacobian.
    !
    ! parameters:
    ! - wL: vector of primitive variables on LHS of interface.
    ! - wR: vector of primitive variables on RHS of interface.
    ! -  g: ratio of specific heats, gamma.
    !------------------------------------------------------
    function roe_flux(wL,wR,g) result(f)
      implicit none
      double precision, intent(in) :: wL(3), wR(3) ! Left and right states 
      double precision, intent(in) :: g            ! cp/cv
      double precision             :: uL(3), uR(3) ! Left and right states
      double precision             :: a(3)         ! Wave speeds
      double precision             :: l(3)         ! Eigenvalues of averaged J
      double precision             :: R(3,3)       ! Eigenvectors
      double precision             :: ut,Ht,at     ! Roe averaged variables
      double precision             :: rhoL,vL,pL   ! Left primitives
      double precision             :: rhoR,vR,pR   ! Right primitives
      double precision             :: HL,HR,aL,aR  ! Prims to be computed
      double precision             :: du1,du2,du3  ! Jumps in conserved vars
      double precision             :: rhohat       ! average density
      integer                      :: i            ! dummy indices
      double precision             :: fL(3),fR(3)  ! Fluxes from Euler eqs
      double precision             :: sum_term(3)  ! Term used in final calc
      double precision             :: f(3)         ! Approximate flux

      ! Converting states to conservative form for convenience
      uL = w_to_u(wL,g)
      uR = w_to_u(wR,g)

      ! Computing primitive variables
      rhoL = wL(1)
      rhoR = wR(1)
      vL   = wL(2)
      vR   = wR(2)
      pL   = wL(3)
      pR   = wR(3)
      aL   = sqrt(g*pL/rhoL)
      aR   = sqrt(g*pR/rhoR)
      HL   = (uL(3)+pL)/rhoL
      HR   = (uR(3)+pR)/rhoR

      ! Computing the Roe average values for u,H and a
      ut = (sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
      Ht = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
      at = sqrt((g-1.0d0)*(Ht - 0.5d0*ut**2))

      ! Computing the averaged eigenvalues 
      l(1) = ut - at
      l(2) = ut
      l(3) = ut + at

      ! Computing the averaged right eigenvectors
      R(1,1) = 1.0d0
      R(2,1) = ut - at
      R(3,1) = Ht - ut*at
      R(1,2) = 1.0d0
      R(2,2) = ut
      R(3,2) = 0.5d0*ut**2
      R(1,3) = 1.0d0
      R(2,3) = ut + at
      R(3,3) = Ht + ut*at

      ! Computing the jumps in the conserved variables
      du1 = uR(1) - uL(1)
      du2 = uR(2) - uL(2)
      du3 = uR(3) - uL(3)

      ! Computing the wave strengths
      a(2) = (g-1.0d0)/at**2*(du1*(Ht-ut**2)+ut*du2-du3)
      a(1) = 1.0d0/(2.0d0*at)*(du1*(ut+at)-du2-at*a(2))
      a(3) = du1 - (a(1) + a(2))

      ! Computing the approximate flux at i+1/2 (i.e., the RHS)
      fL = euler_flux(wL,g)
      fR = euler_flux(wR,g)
      do i=1,3
        sum_term = sum_term + a(i)*abs(l(i))*R(:,i)
      end do
      f = 0.5d0*(fL+fR) - 0.5d0*sum_term

    end function roe_flux

 function roe_flux2(wL,wR,gamma) result(flux)

 implicit none
 double precision, intent(in) :: wL(3), wR(3), gamma !  Input (conservative variables rho*[1, v, E])
 double precision             :: flux(3)             ! Output (numerical flux across L and R states)

!Local parameters
 double precision, parameter ::    zero = 0.0d0
 double precision, parameter ::     one = 1.0d0
 double precision, parameter ::    four = 4.0d0
 double precision, parameter ::    half = 0.5d0
 double precision, parameter :: quarter = 0.25d0
!Local variables
 double precision :: uL(3), uR(3)
 double precision :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 double precision :: aL, aR, HL, HR               ! Speeds of sound.
 double precision :: RT,rho,v,H,a                 ! Roe-averages
 double precision :: drho,du,dP,dV(3)
 double precision :: ws(3),Da, R(3,3)
 integer :: j, k

    uL = w_to_u(wL,gamma)
    uR = w_to_u(wR,gamma)

!Primitive and other variables.
!  Left state
    rhoL = wL(1)
      vL = wL(2)
      pL = wL(3)
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = wR(1)
      vR = wR(2)
      pR = wR(3)
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

!First compute the Roe Averages **************************
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*v*v) )

!Differences in primitive variables.
   drho = rhoR - rhoL
     du =   vR - vL
     dP =   pR - pL

!Wave strength (Characteristic Variables).
   dV(1) =  half*(dP-rho*a*du)/(a*a)
   dV(2) = -( dP/(a*a) - drho )
   dV(3) =  half*(dP+rho*a*du)/(a*a)

!Absolute values of the wave speeds (Eigenvalues)
   ws(1) = abs(v-a)
   ws(2) = abs(v  )
   ws(3) = abs(v+a)

!Modified wave speeds for nonlinear fields (the so-called entropy fix, which
!is often implemented to remove non-physical expansion shocks).
!There are various ways to implement the entropy fix. This is just one
!example. Try turn this off. The solution may be more accurate.
   Da = max(zero, four*((vR-aR)-(vL-aL)) )
   if (ws(1) < half*Da) ws(1) = ws(1)*ws(1)/Da + quarter*Da
   Da = max(zero, four*((vR+aR)-(vL+aL)) )
   if (ws(3) < half*Da) ws(3) = ws(3)*ws(3)/Da + quarter*Da

!Right eigenvectors
   R(1,1) = one
   R(2,1) = v - a
   R(3,1) = H - v*a

   R(1,2) = one
   R(2,2) = v
   R(3,2) = half*v*v

   R(1,3) = one
   R(2,3) = v + a
   R(3,3) = H + v*a

!Compute the average flux.
   flux = half*( euler_flux(wL,gamma) + euler_flux(wR,gamma) )

!Add the matrix dissipation term to complete the Roe flux.
  do j = 1, 3
   do k = 1, 3
    flux(j) = flux(j) - half*ws(k)*dV(k)*R(j,k) 
   end do
  end do

 end function roe_flux2


end module flux
