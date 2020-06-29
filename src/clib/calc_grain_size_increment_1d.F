#ifdef GRACKLE_MD
#include "phys_const.def"

!     Compute grain size increment

      subroutine calc_grain_size_increment_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , metal, metal_SN1, metal_SN2
     &              , SN1_fSiM, SN1_fFeM, SN1_fMg2SiO4, SN1_fMgSiO3
     &              , SN1_fFe3O4, SN1_fAC, SN1_fSiO2D, SN1_fMgO
     &              , SN1_fFeS, SN1_fAl2O3
     &              , SN1_r0SiM, SN1_r0FeM, SN1_r0Mg2SiO4, SN1_r0MgSiO3
     &              , SN1_r0Fe3O4, SN1_r0AC, SN1_r0SiO2D, SN1_r0MgO
     &              , SN1_r0FeS, SN1_r0Al2O3
     &              , SN1_a0SiM, SN1_a0FeM, SN1_a0Mg2SiO4, SN1_a0MgSiO3
     &              , SN1_a0Fe3O4, SN1_a0AC, SN1_a0SiO2D, SN1_a0MgO
     &              , SN1_a0FeS, SN1_a0Al2O3
     &              , SN1_v0SiM, SN1_v0FeM, SN1_v0Mg2SiO4, SN1_v0MgSiO3
     &              , SN1_v0Fe3O4, SN1_v0AC, SN1_v0SiO2D, SN1_v0MgO
     &              , SN1_v0FeS, SN1_v0Al2O3
     &              , SN1_N
     &              , SN1_RSiM, SN1_RFeM, SN1_RMg2SiO4, SN1_RMgSiO3
     &              , SN1_RFe3O4, SN1_RAC, SN1_RSiO2D, SN1_RMgO
     &              , SN1_RFeS, SN1_RAl2O3
     &              , SN2_fSiM, SN2_fFeM, SN2_fMg2SiO4, SN2_fMgSiO3
     &              , SN2_fFe3O4, SN2_fAC, SN2_fSiO2D, SN2_fMgO
     &              , SN2_fFeS, SN2_fAl2O3
     &              , SN2_r0SiM, SN2_r0FeM, SN2_r0Mg2SiO4, SN2_r0MgSiO3
     &              , SN2_r0Fe3O4, SN2_r0AC, SN2_r0SiO2D, SN2_r0MgO
     &              , SN2_r0FeS, SN2_r0Al2O3
     &              , SN2_a0SiM, SN2_a0FeM, SN2_a0Mg2SiO4, SN2_a0MgSiO3
     &              , SN2_a0Fe3O4, SN2_a0AC, SN2_a0SiO2D, SN2_a0MgO
     &              , SN2_a0FeS, SN2_a0Al2O3
     &              , SN2_v0SiM, SN2_v0FeM, SN2_v0Mg2SiO4, SN2_v0MgSiO3
     &              , SN2_v0Fe3O4, SN2_v0AC, SN2_v0SiO2D, SN2_v0MgO
     &              , SN2_v0FeS, SN2_v0Al2O3
     &              , SN2_N
     &              , SN2_RSiM, SN2_RFeM, SN2_RMg2SiO4, SN2_RMgSiO3
     &              , SN2_RFe3O4, SN2_RAC, SN2_RSiO2D, SN2_RMgO
     &              , SN2_RFeS, SN2_RAl2O3
     &              , FSN1_SiM, FSN1_FeM, FSN1_Mg2SiO4, FSN1_MgSiO3
     &              , FSN1_Fe3O4, FSN1_AC, FSN1_SiO2D, FSN1_MgO
     &              , FSN1_FeS, FSN1_Al2O3
     &              , FSN2_SiM, FSN2_FeM, FSN2_Mg2SiO4, FSN2_MgSiO3
     &              , FSN2_Fe3O4, FSN2_AC, FSN2_SiO2D, FSN2_MgO
     &              , FSN2_FeS, FSN2_Al2O3
     &            )

      implicit NONE
#include "grackle_fortran_types.def"

!     in
      logical itmask(in)
      integer igrgr
      integer in, jn, kn, is, ie, j, k
      R_PREC  SiM(in,jn,kn)
     &      , FeM(in,jn,kn)  , Mg2SiO4(in,jn,kn), MgSiO3(in,jn,kn)
     &      , Fe3O4(in,jn,kn), AC(in,jn,kn)     , SiO2D(in,jn,kn)
     &      , MgO(in,jn,kn)  , FeS(in,jn,kn)    , Al2O3(in,jn,kn)
     &      , metal(in,jn,kn)
     &      , metal_SN1(in,jn,kn), metal_SN2(in,jn,kn)
!     table
      real*8  SN1_fSiM, SN1_fFeM, SN1_fMg2SiO4, SN1_fMgSiO3
     &      , SN1_fFe3O4, SN1_fAC, SN1_fSiO2D, SN1_fMgO
     &      , SN1_fFeS, SN1_fAl2O3
      real*8  SN1_r0SiM, SN1_r0FeM, SN1_r0Mg2SiO4, SN1_r0MgSiO3
     &      , SN1_r0Fe3O4, SN1_r0AC, SN1_r0SiO2D, SN1_r0MgO
     &      , SN1_r0FeS, SN1_r0Al2O3
      real*8  SN1_a0SiM, SN1_a0FeM, SN1_a0Mg2SiO4, SN1_a0MgSiO3
     &      , SN1_a0Fe3O4, SN1_a0AC, SN1_a0SiO2D, SN1_a0MgO
     &      , SN1_a0FeS, SN1_a0Al2O3
      real*8  SN1_v0SiM, SN1_v0FeM, SN1_v0Mg2SiO4, SN1_v0MgSiO3
     &      , SN1_v0Fe3O4, SN1_v0AC, SN1_v0SiO2D, SN1_v0MgO
     &      , SN1_v0FeS, SN1_v0Al2O3
      integer SN1_N(3)
      real*8  SN1_RSiM(SN1_N(3)), SN1_RFeM(SN1_N(3))
     &      , SN1_RMg2SiO4(SN1_N(3)), SN1_RMgSiO3(SN1_N(3))
     &      , SN1_RFe3O4(SN1_N(3)), SN1_RAC(SN1_N(3))
     &      , SN1_RSiO2D(SN1_N(3)), SN1_RMgO(SN1_N(3))
     &      , SN1_RFeS(SN1_N(3)), SN1_RAl2O3(SN1_N(3))
      real*8  SN2_fSiM, SN2_fFeM, SN2_fMg2SiO4, SN2_fMgSiO3
     &      , SN2_fFe3O4, SN2_fAC, SN2_fSiO2D, SN2_fMgO
     &      , SN2_fFeS, SN2_fAl2O3
      real*8  SN2_r0SiM, SN2_r0FeM, SN2_r0Mg2SiO4, SN2_r0MgSiO3
     &      , SN2_r0Fe3O4, SN2_r0AC, SN2_r0SiO2D, SN2_r0MgO
     &      , SN2_r0FeS, SN2_r0Al2O3
      real*8  SN2_a0SiM, SN2_a0FeM, SN2_a0Mg2SiO4, SN2_a0MgSiO3
     &      , SN2_a0Fe3O4, SN2_a0AC, SN2_a0SiO2D, SN2_a0MgO
     &      , SN2_a0FeS, SN2_a0Al2O3
      real*8  SN2_v0SiM, SN2_v0FeM, SN2_v0Mg2SiO4, SN2_v0MgSiO3
     &      , SN2_v0Fe3O4, SN2_v0AC, SN2_v0SiO2D, SN2_v0MgO
     &      , SN2_v0FeS, SN2_v0Al2O3
      integer SN2_N(3)
      real*8  SN2_RSiM(SN2_N(3)), SN2_RFeM(SN2_N(3))
     &      , SN2_RMg2SiO4(SN2_N(3)), SN2_RMgSiO3(SN2_N(3))
     &      , SN2_RFe3O4(SN2_N(3)), SN2_RAC(SN2_N(3))
     &      , SN2_RSiO2D(SN2_N(3)), SN2_RMgO(SN2_N(3))
     &      , SN2_RFeS(SN2_N(3)), SN2_RAl2O3(SN2_N(3))
!     out
      real*8  FSN1_SiM(in), FSN1_FeM(in), FSN1_Mg2SiO4(in)
     &      , FSN1_MgSiO3(in), FSN1_Fe3O4(in), FSN1_AC(in)
     &      , FSN1_SiO2D(in), FSN1_MgO(in), FSN1_FeS(in)
     &      , FSN1_Al2O3(in)
      real*8  FSN2_SiM(in), FSN2_FeM(in), FSN2_Mg2SiO4(in)
     &      , FSN2_MgSiO3(in), FSN2_Fe3O4(in), FSN2_AC(in)
     &      , FSN2_SiO2D(in), FSN2_MgO(in), FSN2_FeS(in)
     &      , FSN2_Al2O3(in)
      integer i;


      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , SiM, metal_SN1, metal_SN2
     &              , SN1_fSiM
     &              , SN1_r0SiM, SN1_a0SiM, SN1_v0SiM
     &              , SN1_N, SN1_RSiM
     &              , SN2_fSiM
     &              , SN2_r0SiM, SN2_a0SiM, SN2_v0SiM
     &              , SN2_N, SN2_RSiM
     &              , FSN1_SiM, FSN2_SiM
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , FeM, metal_SN1, metal_SN2
     &              , SN1_fFeM
     &              , SN1_r0FeM, SN1_a0FeM, SN1_v0FeM
     &              , SN1_N, SN1_RFeM
     &              , SN2_fFeM
     &              , SN2_r0FeM, SN2_a0FeM, SN2_v0FeM
     &              , SN2_N, SN2_RFeM
     &              , FSN1_FeM, FSN2_FeM
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , Mg2SiO4, metal_SN1, metal_SN2
     &              , SN1_fMg2SiO4
     &              , SN1_r0Mg2SiO4, SN1_a0Mg2SiO4, SN1_v0Mg2SiO4
     &              , SN1_N, SN1_RMg2SiO4
     &              , SN2_fMg2SiO4
     &              , SN2_r0Mg2SiO4, SN2_a0Mg2SiO4, SN2_v0Mg2SiO4
     &              , SN2_N, SN2_RMg2SiO4
     &              , FSN1_Mg2SiO4, FSN2_Mg2SiO4
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , MgSiO3, metal_SN1, metal_SN2
     &              , SN1_fMgSiO3
     &              , SN1_r0MgSiO3, SN1_a0MgSiO3, SN1_v0MgSiO3
     &              , SN1_N, SN1_RMgSiO3
     &              , SN2_fMgSiO3
     &              , SN2_r0MgSiO3, SN2_a0MgSiO3, SN2_v0MgSiO3
     &              , SN2_N, SN2_RMgSiO3
     &              , FSN1_MgSiO3, FSN2_MgSiO3
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , Fe3O4, metal_SN1, metal_SN2
     &              , SN1_fFe3O4
     &              , SN1_r0Fe3O4, SN1_a0Fe3O4, SN1_v0Fe3O4
     &              , SN1_N, SN1_RFe3O4
     &              , SN2_fFe3O4
     &              , SN2_r0Fe3O4, SN2_a0Fe3O4, SN2_v0Fe3O4
     &              , SN2_N, SN2_RFe3O4
     &              , FSN1_Fe3O4, FSN2_Fe3O4
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , AC, metal_SN1, metal_SN2
     &              , SN1_fAC
     &              , SN1_r0AC, SN1_a0AC, SN1_v0AC
     &              , SN1_N, SN1_RAC
     &              , SN2_fAC
     &              , SN2_r0AC, SN2_a0AC, SN2_v0AC
     &              , SN2_N, SN2_RAC
     &              , FSN1_AC, FSN2_AC
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , SiO2D, metal_SN1, metal_SN2
     &              , SN1_fSiO2D
     &              , SN1_r0SiO2D, SN1_a0SiO2D, SN1_v0SiO2D
     &              , SN1_N, SN1_RSiO2D
     &              , SN2_fSiO2D
     &              , SN2_r0SiO2D, SN2_a0SiO2D, SN2_v0SiO2D
     &              , SN2_N, SN2_RSiO2D
     &              , FSN1_SiO2D, FSN2_SiO2D
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , MgO, metal_SN1, metal_SN2
     &              , SN1_fMgO
     &              , SN1_r0MgO, SN1_a0MgO, SN1_v0MgO
     &              , SN1_N, SN1_RMgO
     &              , SN2_fMgO
     &              , SN2_r0MgO, SN2_a0MgO, SN2_v0MgO
     &              , SN2_N, SN2_RMgO
     &              , FSN1_MgO, FSN2_MgO
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , FeS, metal_SN1, metal_SN2
     &              , SN1_fFeS
     &              , SN1_r0FeS, SN1_a0FeS, SN1_v0FeS
     &              , SN1_N, SN1_RFeS
     &              , SN2_fFeS
     &              , SN2_r0FeS, SN2_a0FeS, SN2_v0FeS
     &              , SN2_N, SN2_RFeS
     &              , FSN1_FeS, FSN2_FeS
     &            )

      call calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , Al2O3, metal_SN1, metal_SN2
     &              , SN1_fAl2O3
     &              , SN1_r0Al2O3, SN1_a0Al2O3, SN1_v0Al2O3
     &              , SN1_N, SN1_RAl2O3
     &              , SN2_fAl2O3
     &              , SN2_r0Al2O3, SN2_a0Al2O3, SN2_v0Al2O3
     &              , SN2_N, SN2_RAl2O3
     &              , FSN1_Al2O3, FSN2_Al2O3
     &            )

      return
      end

!=======================================================================
!//////////////  CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D \\\\\\\\\\\\\\\\\\

      subroutine calc_grain_size_increment_species_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , gr, metal_SN1, metal_SN2
     &              , SN1_fgr
     &              , SN1_r0gr, SN1_a0gr, SN1_v0gr
     &              , SN1_N, SN1_Rgr
     &              , SN2_fgr
     &              , SN2_r0gr, SN2_a0gr, SN2_v0gr
     &              , SN2_N, SN2_Rgr
     &              , FSN1_gr, FSN2_gr
     &            )


      implicit NONE
! input
      logical itmask(in)
      integer igrgr
      integer in, jn, kn, is, ie, j, k
      R_PREC gr(in,jn,kn)
      R_PREC metal_SN1(in,jn,kn), metal_SN2(in,jn,kn)
      integer SN1_N(3)
      real*8  SN1_fgr, SN1_r0gr, SN1_a0gr, SN1_v0gr, SN1_Rgr(SN1_N(3))
      integer SN2_N(3)
      real*8  SN2_fgr, SN2_r0gr, SN2_a0gr, SN2_v0gr, SN2_Rgr(SN2_N(3))
! output
      real*8 FSN1_gr(in), FSN2_gr(in)
! local
      integer i
      real*8 coef0, coef1, coef2, coef3
      real*8 dr
      real*8 d_SN1, n_SN1
      real*8 d_SN2, n_SN2
      real*8 dRgr, dRgr2, dRgr3
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)
      real*8 a2r, m2r
      parameter (a2r = 4.d0 * pi)
      parameter (m2r = 4.d0 * pi / 3.d0)

      do i = is+1, ie+1
         if ( itmask(i) ) then

         if(igrgr .eq. 0) then

            FSN1_gr(i) = 0.d0
            FSN2_gr(i) = 0.d0

         else

            if( SN1_fgr .gt. 0.d0) then
               d_SN1 = SN1_fgr * metal_SN1(i,j,k)
               n_SN1 = d_SN1 / SN1_v0gr
            else
               d_SN1 = 0.d0
               n_SN1 = 0.d0
            endif
            if( SN2_fgr .gt. 0.d0) then
               d_SN2 = SN2_fgr * metal_SN2(i,j,k)
               n_SN2 = d_SN2 / SN2_v0gr
            else
               d_SN2 = 0.d0
               n_SN2 = 0.d0
            endif
            
            if ( (d_SN1 + d_SN2) .gt. 1.e-20 ) then
               coef0 = d_SN1 + d_SN2 - gr(i,j,k)
               coef1 = n_SN1 * SN1_a0gr
     &               + n_SN2 * SN2_a0gr
               coef2 =(n_SN1 * SN1_r0gr
     &               + n_SN2 * SN2_r0gr
     &                ) * a2r
               coef3 =(n_SN1 + n_SN2) * m2r
               call solve_cubic_equation(
     &           coef2/coef3, coef1/coef3, coef0/coef3, dr)
               dRgr = min(max(dr, 0.d0), SN1_Rgr(SN1_N(3)))
            else
               dRgr = 0.d0
            endif

            dRgr2 = dRgr * dRgr
            dRgr3 = dRgr2 * dRgr

            if (SN1_v0gr .gt. 0.d0) then
               FSN1_gr(i) = log10(
     &                      (   SN1_v0gr
     &                    +     SN1_a0gr*dRgr
     &                    + a2r*SN1_r0gr*dRgr2
     &                    + m2r         *dRgr3
     &                      ) / SN1_v0gr )
            else
               FSN1_gr(i) = 0.d0
            endif
            if (SN2_v0gr .gt. 0.d0) then
               FSN2_gr(i) = log10(
     &                      (   SN2_v0gr
     &                    +     SN2_a0gr*dRgr
     &                    + a2r*SN2_r0gr*dRgr2
     &                    + m2r         *dRgr3
     &                      ) / SN2_v0gr )
            else
               FSN2_gr(i) = 0.d0
            endif

         endif

         endif
      enddo

      return
      end

!=======================================================================
!////////////////  SUBROUTINE SOLVE_CUBIC_EQUATION \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_cubic_equation(a, b, c, root)

      implicit NONE

!  General Arguments

      real*8 a, b, c
      real*8 root
      real*8 root1, root2, root3

!  Locals

      real*8 q, r, m
      real*8 th
      real*8 s,t
      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      q = (a*a - 3.d0*b)/9.d0
      r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c)/54.d0
      m = r*r - q*q*q
    
      if (m .lt. 0.d0) then !! three real roots 
        write(*,*) "something's wrong"
        stop
        th = acos( r / sqrt(q*q*q) )
        root1 = -(2.d0*sqrt(q)*cos( th         /3.d0))-a/3.d0
        root2 = -(2.d0*sqrt(q)*cos((th+2.d0*pi)/3.d0))-a/3.d0
        root3 = -(2.d0*sqrt(q)*cos((th-2.d0*pi)/3.d0))-a/3.d0
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
#endif
