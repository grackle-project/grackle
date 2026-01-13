#include "phys_const.def"

! This file originally defined the solve_rate_cool_g subroutine and a
! number of helper subroutines
! - solve_rate_cool_g has now been transcribed to C++
! - the remaining helper subroutines still need to be transcribed

c -----------------------------------------------------------
!   This routine ensures that the species aren't below tiny.

      subroutine ceiling_species_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, dust,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, idustfield
     &                   , DM, HDII, HeHII
     &                   , imabund, imchem, idspecies, immulti
     &                   , igrgr, idsub
     &                   , CI, CII, CO, CO2
     &                   , OI, OH, H2O, O2
     &                   , SiI, SiOI, SiO2I
     &                   , CH, CH2, COII, OII
     &                   , OHII, H2OII, H3OII, O2II
     &                   , Mg, Al, S, Fe
     &                   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                   , AC, SiO2D, MgO, FeS, Al2O3
     &                   , reforg, volorg, H2Oice
     &                   , metal_loc
     &                   , metal_C13, metal_C20, metal_C25, metal_C30
     &                   , metal_F13, metal_F15, metal_F50, metal_F80
     &                   , metal_P170, metal_P200, metal_Y19)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
     &      , idustfield
      R_PREC  d(in,jn,kn),
     &        de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn), dust(in,jn,kn)
      integer imabund, imchem, idspecies, immulti, igrgr, idsub
       R_PREC DM(in,jn,kn), HDII(in,jn,kn), HeHII(in,jn,kn)
     &      , CI(in,jn,kn)   , CII(in,jn,kn)    , CO(in,jn,kn)
     &      , CO2(in,jn,kn)  , OI(in,jn,kn)     , OH(in,jn,kn)
     &      , H2O(in,jn,kn)  , O2(in,jn,kn)     , SiI(in,jn,kn)
     &      , SiOI(in,jn,kn) , SiO2I(in,jn,kn)  , CH(in,jn,kn)
     &      , CH2(in,jn,kn)  , COII(in,jn,kn)   , OII(in,jn,kn)
     &      , OHII(in,jn,kn) , H2OII(in,jn,kn)  , H3OII(in,jn,kn)
     &      , O2II(in,jn,kn) , Mg(in,jn,kn)     , Al(in,jn,kn)
     &      , S(in,jn,kn)    , Fe(in,jn,kn)     
      R_PREC  SiM(in,jn,kn), FeM(in,jn,kn), Mg2SiO4(in,jn,kn)
     &      , MgSiO3(in,jn,kn), Fe3O4(in,jn,kn), AC(in,jn,kn)
     &      , SiO2D(in,jn,kn), MgO(in,jn,kn), FeS(in,jn,kn)
     &      , Al2O3(in,jn,kn)
     &      , reforg(in,jn,kn), volorg(in,jn,kn), H2Oice(in,jn,kn)
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C13(in,jn,kn), metal_C20(in,jn,kn)
     &      , metal_C25(in,jn,kn), metal_C30(in,jn,kn)
     &      , metal_F13(in,jn,kn), metal_F15(in,jn,kn)
     &      , metal_F50(in,jn,kn), metal_F80(in,jn,kn)
     &      , metal_P170(in,jn,kn), metal_P200(in,jn,kn)
     &      , metal_Y19(in,jn,kn)

!     locals

      integer i, j, k

      if (ispecies .gt. 0) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  de(i,j,k)    = max(de(i,j,k), tiny)
                  HI(i,j,k)    = max(HI(i,j,k), tiny)
                  HII(i,j,k)   = max(HII(i,j,k), tiny)
                  HeI(i,j,k)   = max(HeI(i,j,k), tiny)
                  HeII(i,j,k)  = max(HeII(i,j,k), tiny)
                  HeIII(i,j,k) = max(HeIII(i,j,k), 1e-5_RKIND*tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = max(HM(i,j,k), tiny)
                  H2I(i,j,k)  = max(H2I(i,j,k), tiny)
                  H2II(i,j,k) = max(H2II(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = max(DI(i,j,k), tiny)
                  DII(i,j,k) = max(DII(i,j,k), tiny)
                  HDI(i,j,k) = max(HDI(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 3) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DM(i,j,k)      = max(DM(i,j,k), tiny)
                  HDII(i,j,k)    = max(HDII(i,j,k), tiny)
                  HeHII(i,j,k)   = max(HeHII(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = max(metal(i,j,k), tiny)
                  if (metal(i,j,k) .gt. d(i,j,k)) then
                     write(6, *) 'WARNING: metal density exceeds ',
     &                    'total density!'
                     write(6, *) 'i, j, k, metal, density = ',
     &                    i, j, k, metal(i,j,k), d(i,j,k)
                  endif
C                 if( immulti .gt. 0 ) then
C                    metal_loc(i,j,k) = max(metal_loc(i,j,k), tiny)
C                    metal_C13(i,j,k) = max(metal_C13(i,j,k), tiny)
C                    metal_C20(i,j,k) = max(metal_C20(i,j,k), tiny)
C                    metal_C25(i,j,k) = max(metal_C25(i,j,k), tiny)
C                    metal_C30(i,j,k) = max(metal_C30(i,j,k), tiny)
C                    metal_F13(i,j,k) = max(metal_F13(i,j,k), tiny)
C                    metal_F15(i,j,k) = max(metal_F15(i,j,k), tiny)
C                    metal_F50(i,j,k) = max(metal_F50(i,j,k), tiny)
C                    metal_F80(i,j,k) = max(metal_F80(i,j,k), tiny)
C                    metal_P170(i,j,k)= max(metal_P170(i,j,k),tiny)
C                    metal_P200(i,j,k)= max(metal_P200(i,j,k),tiny)
C                    metal_Y19(i,j,k) = max(metal_Y19(i,j,k), tiny)
C                 endif
!!                if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
                  if (imchem .eq. 1) then
                     CI(i,j,k)      = max(CI(i,j,k), tiny)
                     CII(i,j,k)     = max(CII(i,j,k), tiny)
                     CO(i,j,k)      = max(CO(i,j,k), tiny)
                     CO2(i,j,k)     = max(CO2(i,j,k), tiny)
                     OI(i,j,k)      = max(OI(i,j,k), tiny)
                     OH(i,j,k)      = max(OH(i,j,k), tiny)
                     H2O(i,j,k)     = max(H2O(i,j,k), tiny)
                     O2(i,j,k)      = max(O2(i,j,k), tiny)
                     SiI(i,j,k)     = max(SiI(i,j,k), tiny)
                     SiOI(i,j,k)    = max(SiOI(i,j,k), tiny)
                     SiO2I(i,j,k)   = max(SiO2I(i,j,k), tiny)
                     CH(i,j,k)      = max(CH(i,j,k), tiny)
                     CH2(i,j,k)     = max(CH2(i,j,k), tiny)
                     COII(i,j,k)    = max(COII(i,j,k), tiny)
                     OII(i,j,k)     = max(OII(i,j,k), tiny)
                     OHII(i,j,k)    = max(OHII(i,j,k), tiny)
                     H2OII(i,j,k)   = max(H2OII(i,j,k), tiny)
                     H3OII(i,j,k)   = max(H3OII(i,j,k), tiny)
                     O2II(i,j,k)    = max(O2II(i,j,k), tiny)
                     if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
                     if (idspecies .gt. 0) then
                     Mg(i,j,k)      = max(Mg(i,j,k), tiny)
                     endif
                     if (idspecies .gt. 1) then
                     Al(i,j,k)      = max(Al(i,j,k), tiny)
                     S(i,j,k)       = max(S(i,j,k), tiny)
                     Fe(i,j,k)      = max(Fe(i,j,k), tiny)
                     endif
                     endif
                  endif
!!                endif
               enddo
            enddo
         enddo
      endif
      if (idustfield .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  dust(i,j,k)       = max(dust(i,j,k), tiny)
                  if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
!!                if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
                  if (idspecies .gt. 0) then
                     MgSiO3(i,j,k)  = max(MgSiO3(i,j,k), tiny)
                     AC(i,j,k)      = max(AC(i,j,k), tiny)
                  endif
                  if (idspecies .gt. 1) then
                     SiM(i,j,k)     = max(SiM(i,j,k), tiny)
                     FeM(i,j,k)     = max(FeM(i,j,k), tiny)
                     Mg2SiO4(i,j,k) = max(Mg2SiO4(i,j,k), tiny)
                     Fe3O4(i,j,k)   = max(Fe3O4(i,j,k), tiny)
                     SiO2D(i,j,k)   = max(SiO2D(i,j,k), tiny)
                     MgO(i,j,k)     = max(MgO(i,j,k), tiny)
                     FeS(i,j,k)     = max(FeS(i,j,k), tiny)
                     Al2O3(i,j,k)   = max(Al2O3(i,j,k), tiny)
                  endif
                  if (idspecies .gt. 2) then
                     reforg(i,j,k)  = max(reforg(i,j,k), tiny)
                     volorg(i,j,k)  = max(volorg(i,j,k), tiny)
                     H2Oice(i,j,k)  = max(H2Oice(i,j,k), tiny)
                  endif
!!                endif
                  endif
               enddo
            enddo
         enddo
      endif

      return
      end




! -------------------------------------------------------------------
!  This routine calculates the electron and HI rates of change in
!    order to determine the maximum permitted timestep

      subroutine rate_timestep_g(
     &                     dedot, HIdot, ispecies, anydust,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58, 
     &                     h2dust, ncrn, ncrd1, ncrd2, rhoH, 
     &                     k24shield, k25shield, k26shield, 
     &                     k28shield, k29shield, k30shield, k31shield,
     &                     iradtrans, irt_honly, 
     &                     kphHI, kphHeI, kphHeII,
     &                     itmask, edot, chunit, dom, metal
     &                   , HDI, imchem, CI, OI, OH, CO, H2O
     &                   , idissHDI, kdissHDI, iionZ, kphCI, kphOI
     &                   , idissZ, kdissCO, kdissOH, kdissH2O
     &                          )

! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     arguments

      integer ispecies, is, ie, j, k, in, jn, kn,
     &        iradtrans, irt_honly
      real*8 dedot(in), HIdot(in), dom
      real*8 edot(in)
      MASK_TYPE itmask(in), anydust

!     Density fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC  metal(in,jn,kn)

!      Radiative Transfer Fields
      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn)

      R_PREC  HDI(in,jn,kn)
      integer imchem
      R_PREC  CI(in,jn,kn), OI(in,jn,kn),
     &        CO(in,jn,kn), OH(in,jn,kn), H2O(in,jn,kn)
      integer idissHDI, iionZ, idissZ
      R_PREC  kdissHDI(in,jn,kn), kphCI(in,jn,kn), kphOI(in,jn,kn),
     &        kdissCO(in,jn,kn), kdissOH(in,jn,kn), kdissH2O(in,jn,kn)

      real*8 chunit

!     Rate values

      real*8 k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in), h2dust(in), 
     &       ncrn(in), ncrd1(in), ncrd2(in), rhoH(in), 
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k28shield(in), k29shield(in), k30shield(in),
     &       k31shield(in),
     &       k24, k25, k26, k27, k28, k29, k30

!     locals

      integer i
      real*8 h2heatfac(in), H2delta(in), H2dmag, atten, tau

!     Debug
      integer i_max
      R_PREC  d_max

      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i) .ne. MASK_FALSE) then
!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k3(i)*HeI(i,j,k)*de(i,j,k)/4._DKIND
     &               + k5(i)*HeII(i,j,k)*de(i,j,k)/4._DKIND
     &               - k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k4(i)*HeII(i,j,k)*de(i,j,k)/4._DKIND
     &               - k6(i)*HeIII(i,j,k)*de(i,j,k)/4._DKIND
     &               + k57(i)*HI(i,j,k)*HI(i,j,k)
     &               + k58(i)*HI(i,j,k)*HeI(i,j,k)/4._DKIND
     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4._DKIND
     &               + k26shield(i)*HeI(i,j,k)/4._DKIND)

!     Compute the HI density rate-of-change

            HIdot(i) =
     &               - k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k57(i)*HI(i,j,k)*HI(i,j,k)
     &               - k58(i)*HI(i,j,k)*HeI(i,j,k)/4._DKIND
     &               -      k24shield(i)*HI(i,j,k)

         endif                  ! itmask
         enddo
      else

!         Include molecular hydrogen rates for HIdot

         do i = is+1, ie+1
            if (itmask(i) .ne. MASK_FALSE) then
               HIdot(i) = 
     &               -      k1(i) *de(i,j,k)    *HI(i,j,k)  
     &               -      k7(i) *de(i,j,k)    *HI(i,j,k)
     &               -      k8(i) *HM(i,j,k)    *HI(i,j,k)
     &               -      k9(i) *HII(i,j,k)   *HI(i,j,k)
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2._DKIND
     &               - 2._DKIND*k22(i)*HI(i,j,k)**2 *HI(i,j,k)
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) 
     &               + 2._DKIND*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2._DKIND
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._DKIND
     &               + 2._DKIND*k12(i)*de(i,j,k)    *H2I(i,j,k)/2._DKIND
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k)
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k)
     &               + 2._DKIND*k16(i)*HM(i,j,k)    *HII(i,j,k)
     &               + 2._DKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._DKIND
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._DKIND
     &               -      k57(i)*HI(i,j,k)    *HI(i,j,k)
     &               -      k58(i)*HI(i,j,k)    *HeI(i,j,k)/4._DKIND
     &               -      k24shield(i)*HI(i,j,k)
     &               +   2.0_DKIND*k31shield(i) * H2I(i,j,k)/2.0_DKIND

!     Add H2 formation on dust grains

            if (anydust .ne. MASK_FALSE) then
            if (metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k)) then
               HIdot(i) = HIdot(i) 
     &              - 2._DKIND * h2dust(i) * rhoH(i) * HI(i,j,k)
            endif            ! correct GC20200701
            endif

!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i) * HI(i,j,k)   * de(i,j,k)
     &               + k3(i) * HeI(i,j,k)  * de(i,j,k)/4._DKIND
     &               + k5(i) * HeII(i,j,k) * de(i,j,k)/4._DKIND
     &               + k8(i) * HM(i,j,k)   * HI(i,j,k)
     &               + k15(i)* HM(i,j,k)   * HI(i,j,k)
     &               + k17(i)* HM(i,j,k)   * HII(i,j,k)
     &               + k14(i)* HM(i,j,k)   * de(i,j,k)
     &               - k2(i) * HII(i,j,k)  * de(i,j,k)
     &               - k4(i) * HeII(i,j,k) * de(i,j,k)/4._DKIND
     &               - k6(i) * HeIII(i,j,k)* de(i,j,k)/4._DKIND
     &               - k7(i) * HI(i,j,k)   * de(i,j,k)
     &               - k18(i)* H2II(i,j,k) * de(i,j,k)/2._DKIND
     &               + k57(i)* HI(i,j,k)   * HI(i,j,k)
     &               + k58(i)* HI(i,j,k)   * HeI(i,j,k)/4._DKIND
     &               + (k24shield(i)*HI(i,j,k)
     &               +  k25shield(i)*HeII(i,j,k)/4._DKIND
     &               +  k26shield(i)*HeI(i,j,k)/4._DKIND)

!     HII, HeII, HeIII recombination heating

            edot(i) = edot(i) - chunit * (
     &    13.6_DKIND*( k1(i) * HI(i,j,k)   * de(i,j,k)
     &               - k2(i) * HII(i,j,k)  * de(i,j,k)          )
     &  + 24.6_DKIND*( k3(i) * HeI(i,j,k)  * de(i,j,k)/4._DKIND
     &               - k4(i) * HeII(i,j,k) * de(i,j,k)/4._DKIND )
     &  + 79.0_DKIND*( k5(i) * HeII(i,j,k) * de(i,j,k)/4._DKIND
     &               - k6(i) * HeIII(i,j,k)* de(i,j,k)/4._DKIND )
     &            )

!     H2 formation heating

!     Equation 23 from Omukai (2000)
            h2heatfac(i) = (1._DKIND + (ncrn(i) / (dom *
     &           (HI(i,j,k) * ncrd1(i) +
     &           H2I(i,j,k) * 0.5_DKIND * ncrd2(i)))))**(-1._DKIND)

            ! We only want to apply this if the formation dominates, but we
            ! need to apply it outside the delta calculation.

            H2delta(i) = 
     &          HI(i,j,k) *
     &           ( (3.53_DKIND * k8 (i) * HM(i,j,k)
     &            + 4.48_DKIND * k22(i) * HI(i,j,k)**2._DKIND)
     &            * h2heatfac(i)
     &            - 4.48_DKIND * k13(i) * H2I(i,j,k)/2._DKIND)
            !! corrected by GC 202002

!!          if(H2delta(i).gt.0._DKIND) then
!!            H2delta(i) = H2delta(i) * h2heatfac(i)
!!          endif


            if (anydust .ne. MASK_FALSE) then
            if (metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k)) then
               H2delta(i) = H2delta(i) + 
     &              h2dust(i) * HI(i,j,k) * rhoH(i) * 
     &              (0.2_DKIND + 4.2_DKIND * h2heatfac(i))
            endif
            endif

!            H2dmag = abs(H2delta)/(
!     &          HI(i,j,k)*( k22(i) * HI(i,j,k)**2._DKIND
!     &                    + k13(i) * H2I(i,j,k)/2._DKIND))
!            tau = (H2dmag/1e-5_DKIND)**-1.0_DKIND
!            tau = max(tau, 1.e-5_DKIND)
!            atten = min((1.-exp(-tau))/tau,1._DKIND)
            atten = 1._DKIND
            edot(i) = edot(i) + chunit * H2delta(i) * atten
!     &       + H2I(i,j,k)*( k21(i) * HI(i,j,k)**2.0_DKIND
!     &                    - k23(i) * H2I(i,j,k))
!H * (k22 * H^2 - k13 * H_2) + H_2 * (k21 * H^2 - k23 * H_2) */
         endif                  ! itmask
         enddo

!!!!! output cooling rate for debug
!        d_max = 0._DKIND
!        do i = is+1, ie+1
!        if (itmask(i)) then
!           if (d(i,j,k) .gt. d_max) then
!              i_max = i
!              d_max = d(i,j,k)
!           endif
!        endif
!        enddo

!        if (d_max * dom .gt. 1.e9_DKIND) then
C           open(11,file='nH.dat',status='unknown',form='formatted'
C    &         , access='append'
C    &          )
C           do i = is+1, ie+1
C              if (itmask(i)) then
C                 write(11,'(6E13.5)')
C    &            d(i,j,k)*dom
C    &          , chunit * H2delta(i) * atten / d(i,j,k)
C    &          , 3.53_DKIND * k8 (i) * HM(i,j,k) * HI(i,j,k)
C    &             * h2heatfac(i) * chunit / d(i,j,k)
C    &          , 4.48_DKIND * k22(i) * HI(i,j,k)**2._DKIND * HI(i,j,k)
C    &             * h2heatfac(i) * chunit / d(i,j,k)
C    &          ,-4.48_DKIND * k13(i) * H2I(i,j,k)/2._DKIND * HI(i,j,k)
C    &                            *chunit / d(i,j,k)
C    &          ,(0.2_DKIND + 4.2_DKIND * h2heatfac(i)) * h2dust(i)
C    &             * HI(i,j,k) * rhoH(i) * chunit / d(i,j,k)
C              endif
C           enddo
C           close(11)
!        endif

!        if (d_max * dom .gt. 1.e11_DKIND) then
!           stop
!        endif
        

      endif

!     Add photo-ionization rates if needed

      if (iradtrans .eq. 1) then
         if (irt_honly .eq. 0) then
            do i = is+1, ie+1
               if (itmask(i) .ne. MASK_FALSE) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     &                 + kphHeI(i,j,k) * HeI(i,j,k) / 4._DKIND
     &                 + kphHeII(i,j,k) * HeII(i,j,k) / 4._DKIND
               endif
            enddo
         else
            do i = is+1, ie+1
               if (itmask(i) .ne. MASK_FALSE) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
               endif
            enddo
         endif
         if ((ispecies .gt. 2).and.(idissHDI .gt. 0)) then
            do i = is+1, ie+1
               if (itmask(i) .ne. MASK_FALSE) then
                  HIdot(i) = HIdot(i)
     &              + kdissHDI(i,j,k) * HDI(i,j,k)/3.0_DKIND
               endif
            enddo
         endif
         if ((imchem .gt. 0).and.(iionZ .gt. 0)) then
            do i = is+1, ie+1
               if (itmask(i) .ne. MASK_FALSE) then
                  dedot(i) = dedot(i)
     &              + kphCI(i,j,k) * CI(i,j,k)/12.0_DKIND
     &              + kphOI(i,j,k) * OI(i,j,k)/16.0_DKIND
               endif
            enddo
         endif
         if ((imchem .gt. 0).and.(idissZ .gt. 0)) then
            do i = is+1, ie+1
               if (itmask(i) .ne. MASK_FALSE) then
                  HIdot(i) = HIdot(i)
     &              + kdissOH (i,j,k) * OH(i,j,k) /17.0_DKIND
     &              + kdissH2O(i,j,k) * H2O(i,j,k)/18.0_DKIND
               endif
            enddo
         endif
      endif

      


      return
      end

