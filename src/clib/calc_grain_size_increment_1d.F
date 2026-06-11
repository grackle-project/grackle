#include "phys_const.def"

!=======================================================================
!////////////////  SUBROUTINE SOLVE_CUBIC_EQUATION \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_cubic_equation(a, b, c, root)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      real*8 a, b, c
      real*8 root
      real*8 root1, root2, root3

!  Locals

      real*8 q, r, m
      real*8 th
      real*8 s,t
      real*8, parameter :: pi = pi_val

      q = (a*a - 3.d0*b)/9.d0
      r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c)/54.d0
      m = r*r - q*q*q
    
      if (m .lt. 0.d0) then !! three real roots 
        th = acos( r / sqrt(q*q*q) )
        root1 = -(2.d0*sqrt(q)*cos( th         /3.d0))-a/3.d0
        root2 = -(2.d0*sqrt(q)*cos((th+2.d0*pi)/3.d0))-a/3.d0
        root3 = -(2.d0*sqrt(q)*cos((th-2.d0*pi)/3.d0))-a/3.d0
        write(*,*) "three real roots", root1, root2, root3
        stop
      else !! one real root
        if(r .gt. 0.d0) then
          s = -( r+sqrt(m))**(1.d0/3.d0)
        else
          s =  (-r+sqrt(m))**(1.d0/3.d0)
        endif
        t = q/s
        root = s+t-a/3.d0
      endif

      return
      end
