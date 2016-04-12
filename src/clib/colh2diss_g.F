#include "phys_const.def"
c=======================================================================
c/////////////////////  SUBROUTINE COLH2DISS_G  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine colh2diss_g(tgas, idt, f1, f2, f3, f4, f5, f6, f7)
c
c  COMPUTE DENSITY-DEPENDENT COLLISIONAL H2 DISSOCIATION RATE
c
c  written by: Tom Abel
c  date:       
c  modified1: Feb, 2000 by Greg Bryan; adapted to AMR
c  modified2: November, 2015 by Simon Glover; added dissociative tunnelling
c
c  PURPOSE:
c    Computes the 7 temperature-dependent functions required to
c    generate the density-dependent rate coefficients for the
c    collisional dissociation of H2 by HI. If the flag idt=0, 
c    the routine generates the functions corresponding
c    to direct collisional dissociation:  
c
c      H2 + H  -->  H + H + H
c
c    If the flag idt=1, the routine generates instead the 
c    functions corresponding to the destruction of H2 by
c    dissociative tunnelling:
c
c      H2 + H  --> (H2)* + H  --> H + H + H
c
c    where (H2)* represents a quasi-bound state of H2.
c
c    In both cases, the data for the functions is taken
c    from Martin, Schwartz & Mandy (1996, ApJ 461, 265).
c
c  UNITS:
c
c  PARAMETERS:
c
c  INPUTS:
c     T is the gas temperature in Kelvin
c     idt -- flag determining whether to return direct collisional 
c       dissociation rate (k_CID in Martin et al, 1996, selected 
c       if idt=0) or the dissociative tunnelling rate (k_dt in 
c       Martin et al, 1996, selected if idt = 1)
c         
c  OUTPUTS:
c     f1-7: rates as given below
c
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "grackle_fortran_types.def"
c
c  Arguments
c
      real*8 f1,f2,f3,f4,f5,f6,f7
      real*8 tgas
      integer idt
c
c  Parameters
c
c
c  Locals
c
      real*8 t, tl, a1, a, b, b1, c, c1, d
      real*8 y(21)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c
      t  = tgas
      f1 = tiny
      f2 = tiny
      f3 = tiny
      f4 = tiny
      f5 = 1._DKIND
      f6 = 1._DKIND
      f7 = 0._DKIND
c
c     do not use tiny values for temperatures below 500 K also do not
c     return values for temperatures above 1 million Kelvin. Note that     
c     data and fits are only accurate for t< 5 10^5. However,
c     collisional dissociation by electrons will be dominant above this
c     temperature, anyhow. 
c
      if (tgas .le. 500._DKIND) then
         t=500._DKIND
      endif
      if (t .ge. 1.e6_DKIND) then
         t = 1.e6_DKIND
      endif
c
c     fitting parameters
c
      if (idt .eq. 0) then
        y(1)   =   -1.784239e2_DKIND
        y(2)   =   -6.842243e1_DKIND
        y(3)   =    4.320243e1_DKIND
        y(4)   =   -4.633167e0_DKIND
        y(5)   =    6.970086e1_DKIND
        y(6)   =    4.087038e4_DKIND
        y(7)   =   -2.370570e4_DKIND
        y(8)   =    1.288953e2_DKIND
        y(9)   =   -5.391334e1_DKIND
        y(10)  =    5.315517e0_DKIND
        y(11)  =   -1.973427e1_DKIND
        y(12)  =    1.678095e4_DKIND
        y(13)  =   -2.578611e4_DKIND
        y(14)  =    1.482123e1_DKIND
        y(15)  =   -4.890915e0_DKIND
        y(16)  =    4.749030e-1_DKIND
        y(17)  =   -1.338283e2_DKIND
        y(18)  =   -1.164408e0_DKIND
        y(19)  =    8.227443e-1_DKIND
        y(20)  =    5.864073e-1_DKIND
        y(21)  =   -2.056313e0_DKIND
      elseif (idt .eq. 1) then
        y(1)   =   -1.427664e+02_DKIND
        y(2)   =    4.270741e+01_DKIND
        y(3)   =   -2.027365e+00_DKIND
        y(4)   =   -2.582097e-01_DKIND
        y(5)   =    2.136094e+01_DKIND
        y(6)   =    2.753531e+04_DKIND
        y(7)   =   -2.146779e+04_DKIND
        y(8)   =    6.034928e+01_DKIND
        y(9)   =   -2.743096e+01_DKIND
        y(10)  =    2.676150e+00_DKIND
        y(11)  =   -1.128215e+01_DKIND
        y(12)  =    1.425455e+04_DKIND
        y(13)  =   -2.312520e+04_DKIND
        y(14)  =    9.305564e+00_DKIND
        y(15)  =   -2.464009e+00_DKIND
        y(16)  =    1.985955e-01_DKIND
        y(17)  =    7.430600e+02_DKIND
        y(18)  =   -1.174242e+00_DKIND
        y(19)  =    7.502286e-01_DKIND
        y(20)  =    2.358848e-01_DKIND
        y(21)  =    2.937507e+00_DKIND
      else
        write(0,*) "Unknown value for idt!"
        write(0,*) "Expecting 0 or 1; got ", idt
        stop
      endif
c
      tl=log10(t)
c high density limit
      a =   y(1)+y(2)*tl+y(3)*tl*tl+y(4)*tl*tl*tl
     $     +y(5)*dlog10(1._DKIND+y(6)/t)
      a1=   y(7)/t
c low density limit
      b =   y(8)+y(9)*tl+y(10)*tl*tl+y(11)*dlog10(1._DKIND+y(12)/t)       
      b1=   y(13)/t
c critical density
      c =   y(14)+y(15)*tl+y(16)*tl*tl+y(17)/t
      c1 =  y(18)+c
      d =   y(19)+y(20)*exp(-t/1850._DKIND)+y(21)*exp(-t/440._DKIND)    
c     tabulate the following temperature dependent coefficients:
      f1 = a
      f2 = (a-b)
      f3 = a1   
      f4 = (a1-b1)
      f5 = 10._DKIND**c
      f6 = 10._DKIND**c1
      f7 = d
c
      return
      end
