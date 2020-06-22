#ifdef GRACKLE_MD
#include "phys_const.def"

!     Compute grain size increment

      subroutine lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FSN_SiM, FSN_FeM, FSN_Mg2SiO4, FSN_MgSiO3
     &              , FSN_Fe3O4, FSN_AC, FSN_SiO2D, FSN_MgO
     &              , FSN_FeS, FSN_Al2O3
     &              , SN_fSiM, SN_fFeM, SN_fMg2SiO4, SN_fMgSiO3
     &              , SN_fFe3O4, SN_fAC, SN_fSiO2D, SN_fMgO
     &              , SN_fFeS, SN_fAl2O3
     &              , SN_N, SN_Size
     &              , SN_D, SN_T, SN_dD, SN_dT
     &              , SN_RSiM, SN_RFeM, SN_RMg2SiO4, SN_RMgSiO3
     &              , SN_RFe3O4, SN_RAC, SN_RSiO2D, SN_RMgO
     &              , SN_RFeS, SN_RAl2O3
     &              , SN_dRSiM, SN_dRFeM, SN_dRMg2SiO4, SN_dRMgSiO3
     &              , SN_dRFe3O4, SN_dRAC, SN_dRSiO2D, SN_dRMgO
     &              , SN_dRFeS, SN_dRAl2O3
     &              , RSN_SiM, RSN_FeM, RSN_Mg2SiO4, RSN_MgSiO3
     &              , RSN_Fe3O4, RSN_AC, RSN_SiO2D, RSN_MgO
     &              , RSN_FeS, RSN_Al2O3
     &              , RSN
     &            )

      implicit NONE
#include "grackle_fortran_types.def"
!     in
      logical itmask(in)
      integer in, jn, kn, is, ie, j, k
      real*8 logrho(in), logT(in)
      real*8  FSN_SiM(in), FSN_FeM(in), FSN_Mg2SiO4(in), FSN_MgSiO3(in)
     &      , FSN_Fe3O4(in), FSN_AC(in), FSN_SiO2D(in), FSN_MgO(in)
     &      , FSN_FeS(in), FSN_Al2O3(in)
      real*8  SN_fSiM, SN_fFeM, SN_fMg2SiO4, SN_fMgSiO3
     &      , SN_fFe3O4, SN_fAC, SN_fSiO2D, SN_fMgO
     &      , SN_fFeS, SN_fAl2O3
!     table
      integer SN_N(3), SN_Size
      real*8  SN_D(SN_N(1)), SN_T(SN_N(2)), SN_dD, SN_dT
      real*8  SN_RSiM(SN_N(3)), SN_RFeM(SN_N(3))
     &      , SN_RMg2SiO4(SN_N(3)), SN_RMgSiO3(SN_N(3))
     &      , SN_RFe3O4(SN_N(3)), SN_RAC(SN_N(3))
     &      , SN_RSiO2D(SN_N(3)), SN_RMgO(SN_N(3))
     &      , SN_RFeS(SN_N(3)), SN_RAl2O3(SN_N(3))
      real*8  SN_dRSiM, SN_dRFeM, SN_dRMg2SiO4, SN_dRMgSiO3
     &      , SN_dRFe3O4, SN_dRAC, SN_dRSiO2D, SN_dRMgO
     &      , SN_dRFeS, SN_dRAl2O3
      real*8  RSN_SiM(SN_Size), RSN_FeM(SN_Size)
     &      , RSN_Mg2SiO4(SN_Size), RSN_MgSiO3(SN_Size)
     &      , RSN_Fe3O4(SN_Size), RSN_AC(SN_Size)
     &      , RSN_SiO2D(SN_Size), RSN_MgO(SN_Size)
     &      , RSN_FeS(SN_Size), RSN_Al2O3(SN_Size)
!     out
      real*8  RSN(in)
!     local
      integer i
      real*8  log_R

      do i = is+1, ie+1
         if ( itmask(i) ) then

            RSN(i) = 0.d0

            if (SN_fSiM .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_SiM(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RSiM, SN_dRSiM,
     &        int8(SN_Size), RSN_SiM, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fFeM .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_FeM(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RFeM, SN_dRFeM,
     &        int8(SN_Size), RSN_FeM, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fMg2SiO4 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_Mg2SiO4(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RMg2SiO4, SN_dRMg2SiO4,
     &        int8(SN_Size), RSN_Mg2SiO4, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fMgSiO3 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_MgSiO3(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RMgSiO3, SN_dRMgSiO3,
     &        int8(SN_Size), RSN_MgSiO3, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fFe3O4 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_Fe3O4(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RFe3O4, SN_dRFe3O4,
     &        int8(SN_Size), RSN_Fe3O4, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fAC .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_AC(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RAC, SN_dRAC,
     &        int8(SN_Size), RSN_AC, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fSiO2D .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_SiO2D(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RSiO2D, SN_dRSiO2D,
     &        int8(SN_Size), RSN_SiO2D, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fMgO .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_MgO(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RMgO, SN_dRMgO,
     &        int8(SN_Size), RSN_MgO, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fFeS .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_FeS(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RFeS, SN_dRFeS,
     &        int8(SN_Size), RSN_FeS, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif
            if (SN_fAl2O3 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_Al2O3(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RAl2O3, SN_dRAl2O3,
     &        int8(SN_Size), RSN_Al2O3, log_R)
              RSN(i) = RSN(i) + 1.d1**(log_R)
            endif

         endif
      enddo

      return
      end


! ------------------------------------------------------------------
!     Compute grain growth rate

      subroutine lookup_grain_growth_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FSN_SiM, FSN_FeM, FSN_Mg2SiO4, FSN_MgSiO3
     &              , FSN_Fe3O4, FSN_AC, FSN_SiO2D, FSN_MgO
     &              , FSN_FeS, FSN_Al2O3
     &              , metal_SN
     &              , SN_fSiM, SN_fFeM, SN_fMg2SiO4, SN_fMgSiO3
     &              , SN_fFe3O4, SN_fAC, SN_fSiO2D, SN_fMgO
     &              , SN_fFeS, SN_fAl2O3
     &              , SN_N, SN_Size
     &              , SN_D, SN_T, SN_dD, SN_dT
     &              , SN_RSiM, SN_RFeM, SN_RMg2SiO4, SN_RMgSiO3
     &              , SN_RFe3O4, SN_RAC, SN_RSiO2D, SN_RMgO
     &              , SN_RFeS, SN_RAl2O3
     &              , SN_dRSiM, SN_dRFeM, SN_dRMg2SiO4, SN_dRMgSiO3
     &              , SN_dRFe3O4, SN_dRAC, SN_dRSiO2D, SN_dRMgO
     &              , SN_dRFeS, SN_dRAl2O3
     &              , KSN_SiM, KSN_FeM, KSN_Mg2SiO4, KSN_MgSiO3
     &              , KSN_Fe3O4, KSN_AC, KSN_SiO2D, KSN_MgO
     &              , KSN_FeS, KSN_Al2O3
     &              , kd0, kd1, kd2, kd3, kd4
     &              , kd5, kd6, kd7, kd8, kd9
     &            )

      implicit NONE
#include "grackle_fortran_types.def"
!     in
      logical itmask(in)
      integer in, jn, kn, is, ie, j, k
      real*8 logrho(in), logT(in)
      real*8  FSN_SiM(in), FSN_FeM(in), FSN_Mg2SiO4(in), FSN_MgSiO3(in)
     &      , FSN_Fe3O4(in), FSN_AC(in), FSN_SiO2D(in), FSN_MgO(in)
     &      , FSN_FeS(in), FSN_Al2O3(in)
      R_PREC  metal_SN(in,jn,kn)
      real*8  SN_fSiM, SN_fFeM, SN_fMg2SiO4, SN_fMgSiO3
     &      , SN_fFe3O4, SN_fAC, SN_fSiO2D, SN_fMgO
     &      , SN_fFeS, SN_fAl2O3
!     table
      integer SN_N(3), SN_Size
      real*8  SN_D(SN_N(1)), SN_T(SN_N(2)), SN_dD, SN_dT
      real*8  SN_RSiM(SN_N(3)), SN_RFeM(SN_N(3))
     &      , SN_RMg2SiO4(SN_N(3)), SN_RMgSiO3(SN_N(3))
     &      , SN_RFe3O4(SN_N(3)), SN_RAC(SN_N(3))
     &      , SN_RSiO2D(SN_N(3)), SN_RMgO(SN_N(3))
     &      , SN_RFeS(SN_N(3)), SN_RAl2O3(SN_N(3))
      real*8  SN_dRSiM, SN_dRFeM, SN_dRMg2SiO4, SN_dRMgSiO3
     &      , SN_dRFe3O4, SN_dRAC, SN_dRSiO2D, SN_dRMgO
     &      , SN_dRFeS, SN_dRAl2O3
      real*8  KSN_SiM(SN_Size), KSN_FeM(SN_Size)
     &      , KSN_Mg2SiO4(SN_Size), KSN_MgSiO3(SN_Size)
     &      , KSN_Fe3O4(SN_Size), KSN_AC(SN_Size)
     &      , KSN_SiO2D(SN_Size), KSN_MgO(SN_Size)
     &      , KSN_FeS(SN_Size), KSN_Al2O3(SN_Size)
!     out
      real*8  kd0(in), kd1(in), kd2(in), kd3(in), kd4(in)
     &      , kd5(in), kd6(in), kd7(in), kd8(in), kd9(in)
!     local
      integer i
      real*8  log_kgg

      do i = is+1, ie+1
         if ( itmask(i) ) then
         if ( metal_SN(i,j,k) .gt. 0._DKIND ) then

            if (SN_fSiM .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_SiM(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RSiM, SN_dRSiM,
     &        int8(SN_Size), KSN_SiM, log_kgg)
              kd0(i) = 1.d1**(log_kgg)
            else
              kd0(i) = 0.d0
            endif

            if (SN_fFeM .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_FeM(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RFeM, SN_dRFeM,
     &        int8(SN_Size), KSN_FeM, log_kgg)
              kd1(i) = 1.d1**(log_kgg)
            else
              kd1(i) = 0.d0
            endif

            if (SN_fMg2SiO4 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_Mg2SiO4(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RMg2SiO4, SN_dRMg2SiO4,
     &        int8(SN_Size), KSN_Mg2SiO4, log_kgg)
              kd2(i) = 1.d1**(log_kgg)
            else
              kd2(i) = 0.d0
            endif

            if (SN_fMgSiO3 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_MgSiO3(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RMgSiO3, SN_dRMgSiO3,
     &        int8(SN_Size), KSN_MgSiO3, log_kgg)
              kd3(i) = 1.d1**(log_kgg)
            else
              kd3(i) = 0.d0
            endif

            if (SN_fFe3O4 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_Fe3O4(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RFe3O4, SN_dRFe3O4,
     &        int8(SN_Size), KSN_Fe3O4, log_kgg)
              kd4(i) = 1.d1**(log_kgg)
            else
              kd4(i) = 0.d0
            endif

            if (SN_fAC .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_AC(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RAC, SN_dRAC,
     &        int8(SN_Size), KSN_AC, log_kgg)
              kd5(i) = 1.d1**(log_kgg)
            else
              kd5(i) = 0.d0
            endif

            if (SN_fSiO2D .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_SiO2D(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RSiO2D, SN_dRSiO2D,
     &        int8(SN_Size), KSN_SiO2D, log_kgg)
              kd6(i) = 1.d1**(log_kgg)
            else
              kd6(i) = 0.d0
            endif

            if (SN_fMgO .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_MgO(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RMgO, SN_dRMgO,
     &        int8(SN_Size), KSN_MgO, log_kgg)
              kd7(i) = 1.d1**(log_kgg)
            else
              kd7(i) = 0.d0
            endif

            if (SN_fFeS .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_FeS(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RFeS, SN_dRFeS,
     &        int8(SN_Size), KSN_FeS, log_kgg)
              kd8(i) = 1.d1**(log_kgg)
            else
              kd8(i) = 0.d0
            endif

            if (SN_fAl2O3 .gt. 0.d0) then
            call interpolate_3D_g(
     &        logrho(i), logT(i), FSN_Al2O3(i), int8(SN_N), 
     &        SN_D, SN_dD, SN_T, SN_dT, SN_RAl2O3, SN_dRAl2O3,
     &        int8(SN_Size), KSN_Al2O3, log_kgg)
              kd9(i) = 1.d1**(log_kgg)
            else
              kd9(i) = 0.d0
            endif

         endif
         endif
      enddo

      return
      end

#endif
