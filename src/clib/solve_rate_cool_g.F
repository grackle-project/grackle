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


! ------------------------------------------------------------------
!   This routine correct the highest abundence species to
!     insure conservation of particle number and charge.

      subroutine make_consistent_g(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, dust,
     &                     d, is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, fh, dtoh
     &                   , idustfield, imchem, igrgr, dom
     &                   , DM, HDII, HeHII
     &                   , CI, CII, CO, CO2
     &                   , OI, OH, H2O, O2
     &                   , SiI, SiOI, SiO2I
     &                   , CH, CH2, COII, OII
     &                   , OHII, H2OII, H3OII, O2II
     &                   , Mg, Al, S, Fe
     &                   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                   , AC, SiO2D, MgO, FeS, Al2O3
     &                   , reforg, volorg, H2Oice
     &                   , immulti, imabund, idspecies, itdmulti, idsub
     &                   , metal_loc
     &                   , metal_C13, metal_C20, metal_C25, metal_C30
     &                   , metal_F13, metal_F15, metal_F50, metal_F80
     &                   , metal_P170, metal_P200, metal_Y19
     &                   , SN0_N
     &                   , SN0_XC, SN0_XO, SN0_XMg, SN0_XAl, SN0_XSi
     &                   , SN0_XS, SN0_XFe
     &                   , SN0_fC, SN0_fO, SN0_fMg, SN0_fAl, SN0_fSi
     &                   , SN0_fS, SN0_fFe
     &                      )
! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      integer idustfield, imchem, igrgr
      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn), metal(in,jn,kn), dust(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real*8 fh, dtoh
      real*8 dom
       R_PREC DM(in,jn,kn)   , HDII(in,jn,kn)   , HeHII(in,jn,kn)
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
      integer immulti, imabund, idspecies, itdmulti, idsub
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C13(in,jn,kn), metal_C20(in,jn,kn)
     &      , metal_C25(in,jn,kn), metal_C30(in,jn,kn)
     &      , metal_F13(in,jn,kn), metal_F15(in,jn,kn)
     &      , metal_F50(in,jn,kn), metal_F80(in,jn,kn)
     &      , metal_P170(in,jn,kn), metal_P200(in,jn,kn)
     &      , metal_Y19(in,jn,kn)
      integer SN0_N
      real*8  SN0_XC (SN0_N), SN0_XO(SN0_N), SN0_XMg(SN0_N)
     &      , SN0_XAl(SN0_N), SN0_XSi(SN0_N), SN0_XS(SN0_N)
     &      , SN0_XFe(SN0_N)
      real*8  SN0_fC (SN0_N), SN0_fO(SN0_N), SN0_fMg(SN0_N)
     &      , SN0_fAl(SN0_N), SN0_fSi(SN0_N), SN0_fS(SN0_N)
     &      , SN0_fFe(SN0_N)

!     locals

      integer i, j, k
      real*8 totalH(in), totalHe(in),
     &       totalD, metalfree(in)
      R_PREC correctH, correctHe, correctD
      real*8 totalZ
      real*8 totalC, totalO, totalMg, totalAl
     &     , totalSi, totalS, totalFe
      real*8 totalCg, totalOg, totalMgg, totalAlg
     &     , totalSig, totalSg, totalFeg
      real*8 totalCd, totalOd, totalMgd, totalAld
     &     , totalSid, totalSd, totalFed
      R_PREC correctC, correctO, correctMg, correctAl
     &     , correctSi, correctS, correctFe
      R_PREC correctCg, correctOg, correctMgg, correctAlg
     &     , correctSig, correctSg, correctFeg
      R_PREC correctCd, correctOd, correctMgd, correctAld
     &     , correctSid, correctSd, correctFed
      R_PREC correctZ
      integer iSN, nSN, iSN0
      integer SN_i(SN0_N)
      R_PREC  SN_metal(in, SN0_N)
      real*8 Ct(in), Ot(in), Mgt(in), Alt(in)
     &     , Sit(in), St(in), Fet(in)
      real*8 Cg(in), Og(in), Mgg(in), Alg(in)
     &     , Sig(in), Sg(in), Feg(in)
      real*8 Cd(in), Od(in), Mgd(in), Ald(in)
     &     , Sid(in), Sd(in), Fed(in)

!     Loop over all zones

      do k = ks+1, ke+1
      do j = js+1, je+1

!     Compute total densities of H and He
!         (ensure non-negativity)

      if (imetal .eq. 1) then
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k) - metal(i,j,k)
         enddo
      else
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k)
         enddo
      endif

      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
      enddo

!     include molecular hydrogen

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
         enddo
      endif

      if(ispecies .gt. 2) then
         do i = is+1, ie+1
            HDI(i,j,k) = abs(HDI(i,j,k))
            totalH (i) = totalH (i)
     &           + 1._DKIND/3._DKIND*HDI(i,j,k)
         enddo
      endif
      !! GC202005

      if(ispecies .gt. 3) then
         do i = is+1, ie+1
            HDII (i,j,k) = abs(HDII (i,j,k))
            HeHII(i,j,k) = abs(HeHII(i,j,k))
            totalH (i) = totalH (i)
     &           + 1._DKIND/3._DKIND*HDII (i,j,k)
     &           + 1._DKIND/5._DKIND*HeHII(i,j,k)
            totalHe(i) = totalHe(i)
     &           + 4._DKIND/5._DKIND*HeHII(i,j,k)
         enddo
      endif

!     Iteration mask for metal-rich cells

!     do i = is+1, ie + 1
!        itmask_metal(i) = .false.
!     enddo
!     if (imetal .eq. 1) then
!         do i = is+1, ie + 1
!            if (metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k)) then
!               itmask_metal(i) = .true.
!            endif
!         enddo
!     endif

      if(imchem .gt. 0) then
         if(immulti .eq. 0) then
            iSN0 = imabund + 1
            do i = is+1, ie+1
                Ct(i) = SN0_XC (iSN0) * metal(i,j,k)
                Ot(i) = SN0_XO (iSN0) * metal(i,j,k)
               Mgt(i) = SN0_XMg(iSN0) * metal(i,j,k)
               Alt(i) = SN0_XAl(iSN0) * metal(i,j,k)
               Sit(i) = SN0_XSi(iSN0) * metal(i,j,k)
                St(i) = SN0_XS (iSN0) * metal(i,j,k)
               Fet(i) = SN0_XFe(iSN0) * metal(i,j,k)
         
                Cg(i) = SN0_fC (iSN0) * metal(i,j,k)
                Og(i) = SN0_fO (iSN0) * metal(i,j,k)
               Mgg(i) = SN0_fMg(iSN0) * metal(i,j,k)
               Alg(i) = SN0_fAl(iSN0) * metal(i,j,k)
               Sig(i) = SN0_fSi(iSN0) * metal(i,j,k)
                Sg(i) = SN0_fS (iSN0) * metal(i,j,k)
               Feg(i) = SN0_fFe(iSN0) * metal(i,j,k)
            enddo
         
         else

C           do i = is+1, ie+1
C              totalZ = metal_loc(i,j,k)
C    &           + metal_C13(i,j,k) + metal_C20(i,j,k)
C    &           + metal_C25(i,j,k) + metal_C30(i,j,k)
C    &           + metal_F13(i,j,k) + metal_F15(i,j,k)
C    &           + metal_F50(i,j,k) + metal_F80(i,j,k)
C    &           + metal_P170(i,j,k)+ metal_P200(i,j,k)
C    &           + metal_Y19(i,j,k)
C              correctZ = metal(i,j,k) / totalZ
C              metal_loc(i,j,k) = metal_loc(i,j,k) * correctZ
C              metal_C13(i,j,k) = metal_C13(i,j,k) * correctZ
C              metal_C20(i,j,k) = metal_C20(i,j,k) * correctZ
C              metal_C25(i,j,k) = metal_C25(i,j,k) * correctZ
C              metal_C30(i,j,k) = metal_C30(i,j,k) * correctZ
C              metal_F13(i,j,k) = metal_F13(i,j,k) * correctZ
C              metal_F15(i,j,k) = metal_F15(i,j,k) * correctZ
C              metal_F50(i,j,k) = metal_F50(i,j,k) * correctZ
C              metal_F80(i,j,k) = metal_F80(i,j,k) * correctZ
C              metal_P170(i,j,k)= metal_P170(i,j,k)* correctZ
C              metal_P200(i,j,k)= metal_P200(i,j,k)* correctZ
C              metal_Y19(i,j,k) = metal_Y19(i,j,k) * correctZ
C           enddo

            nSN = 12
            SN_i( 1) = 1
            SN_metal(:, 1) = metal_loc(:,j,k)
            SN_i( 2) = 2
            SN_metal(:, 2) = metal_C13(:,j,k)
            SN_i( 3) = 3
            SN_metal(:, 3) = metal_C20(:,j,k)
            SN_i( 4) = 4
            SN_metal(:, 4) = metal_C25(:,j,k)
            SN_i( 5) = 5
            SN_metal(:, 5) = metal_C30(:,j,k)
            SN_i( 6) = 6
            SN_metal(:, 6) = metal_F13(:,j,k)
            SN_i( 7) = 7
            SN_metal(:, 7) = metal_F15(:,j,k)
            SN_i( 8) = 8
            SN_metal(:, 8) = metal_F50(:,j,k)
            SN_i( 9) = 9
            SN_metal(:, 9) = metal_F80(:,j,k)
            SN_i(10) =10
            SN_metal(:,10) = metal_P170(:,j,k)
            SN_i(11) =11
            SN_metal(:,11) = metal_P200(:,j,k)
            SN_i(12) =12
            SN_metal(:,12) = metal_Y19(:,j,k)
         
            do i = is+1, ie+1
                Ct(i) = 0._DKIND
                Cg(i) = 0._DKIND
                Ot(i) = 0._DKIND
                Og(i) = 0._DKIND
               Mgt(i) = 0._DKIND
               Mgg(i) = 0._DKIND
               Alt(i) = 0._DKIND
               Alg(i) = 0._DKIND
               Sit(i) = 0._DKIND
               Sig(i) = 0._DKIND
                St(i) = 0._DKIND
                Sg(i) = 0._DKIND
               Fet(i) = 0._DKIND
               Feg(i) = 0._DKIND
               do iSN = 1, nSN
                  iSN0 = SN_i(iSN)

                   Ct(i) =  Ct(i) + SN0_XC (iSN0) * SN_metal(i,iSN)
                   Ot(i) =  Ot(i) + SN0_XO (iSN0) * SN_metal(i,iSN)
                  Mgt(i) = Mgt(i) + SN0_XMg(iSN0) * SN_metal(i,iSN)
                  Alt(i) = Alt(i) + SN0_XAl(iSN0) * SN_metal(i,iSN)
                  Sit(i) = Sit(i) + SN0_XSi(iSN0) * SN_metal(i,iSN)
                   St(i) =  St(i) + SN0_XS (iSN0) * SN_metal(i,iSN)
                  Fet(i) = Fet(i) + SN0_XFe(iSN0) * SN_metal(i,iSN)

                   Cg(i) =  Cg(i) + SN0_fC (iSN0) * SN_metal(i,iSN)
                   Og(i) =  Og(i) + SN0_fO (iSN0) * SN_metal(i,iSN)
                  Mgg(i) = Mgg(i) + SN0_fMg(iSN0) * SN_metal(i,iSN)
                  Alg(i) = Alg(i) + SN0_fAl(iSN0) * SN_metal(i,iSN)
                  Sig(i) = Sig(i) + SN0_fSi(iSN0) * SN_metal(i,iSN)
                   Sg(i) =  Sg(i) + SN0_fS (iSN0) * SN_metal(i,iSN)
                  Feg(i) = Feg(i) + SN0_fFe(iSN0) * SN_metal(i,iSN)
               enddo
            enddo
            
         endif
            
         do i = is+1, ie+1
             Cd(i) =  Ct(i) -  Cg(i)
             Od(i) =  Ot(i) -  Og(i)
            Mgd(i) = Mgt(i) - Mgg(i)
            Ald(i) = Alt(i) - Alg(i)
            Sid(i) = Sit(i) - Sig(i)
             Sd(i) =  St(i) -  Sg(i)
            Fed(i) = Fet(i) - Feg(i)
         enddo

         do i = is+1, ie+1
!        if (itmask_metal(i)) then
            OH   (i,j,k) = abs(OH   (i,j,k))
            H2O  (i,j,k) = abs(H2O  (i,j,k))
            CH   (i,j,k) = abs(CH   (i,j,k))
            CH2  (i,j,k) = abs(CH2  (i,j,k))
            OHII (i,j,k) = abs(OHII (i,j,k))
            H2OII(i,j,k) = abs(H2OII(i,j,k))
            H3OII(i,j,k) = abs(H3OII(i,j,k))
            totalH(i) = totalH(i)
     &        + OH   (i,j,k)/17._DKIND
     &        + H2O  (i,j,k)/18._DKIND*2._DKIND
     &        + CH   (i,j,k)/13._DKIND
     &        + CH2  (i,j,k)/14._DKIND*2._DKIND
     &        + OHII (i,j,k)/17._DKIND
     &        + H2OII(i,j,k)/18._DKIND*2._DKIND
     &        + H3OII(i,j,k)/19._DKIND*3._DKIND
!        endif
         enddo
      endif

!     Correct densities by keeping fractions the same

      do i = is+1, ie+1
         correctH = real(fh*metalfree(i)/totalH(i), RKIND)
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH

         correctHe = real((1._DKIND - fh)*
     &        metalfree(i)/totalHe(i), RKIND)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe

!     Correct molecular hydrogen-related fractions

         if (ispecies .gt. 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
         if(ispecies .gt. 3) then
!!          HDII (i,j,k) = HDII (i,j,k)*correctH
            HeHII(i,j,k) = HeHII(i,j,k)*correctHe
         endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) +
     &           2._DKIND/3._DKIND*HDI(i,j,k)
            if(ispecies .gt. 3) then
               DM   (i,j,k) = abs(DM   (i,j,k))
               HDII (i,j,k) = abs(HDII (i,j,k))
               totalD = totalD + DM(i,j,k) + 
     &           2._DKIND/3._DKIND*HDII(i,j,k)
            endif
            correctD = real(fh*dtoh*metalfree(i)/totalD, RKIND)
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
            if(ispecies .gt. 3) then
               DM   (i,j,k) = DM   (i,j,k)*correctD
               HDII (i,j,k) = HDII (i,j,k)*correctD
            endif
         enddo
      endif

!     Do the same thing for metal species

      if (imchem .eq. 1) then
      do i = is+1, ie+1
!        if (itmask_metal(i)) then
            CI(i,j,k)      = abs(CI(i,j,k)     )
            CII(i,j,k)     = abs(CII(i,j,k)    )
            CO(i,j,k)      = abs(CO(i,j,k)     )
            CO2(i,j,k)     = abs(CO2(i,j,k)    )
            OI(i,j,k)      = abs(OI(i,j,k)     )
            OH(i,j,k)      = abs(OH(i,j,k)     )
            H2O(i,j,k)     = abs(H2O(i,j,k)    )
            O2(i,j,k)      = abs(O2(i,j,k)     )
            SiI(i,j,k)     = abs(SiI(i,j,k)    )
            SiOI(i,j,k)    = abs(SiOI(i,j,k)   )
            SiO2I(i,j,k)   = abs(SiO2I(i,j,k)  )
            CH(i,j,k)      = abs(CH(i,j,k)     )
            CH2(i,j,k)     = abs(CH2(i,j,k)    )
            COII(i,j,k)    = abs(COII(i,j,k)   )
            OII(i,j,k)     = abs(OII(i,j,k)    )
            OHII(i,j,k)    = abs(OHII(i,j,k)   )
            H2OII(i,j,k)   = abs(H2OII(i,j,k)  )
            H3OII(i,j,k)   = abs(H3OII(i,j,k)  )
            O2II(i,j,k)    = abs(O2II(i,j,k)   )
            if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 0) then
                  Mg(i,j,k)      = abs(Mg(i,j,k)     )
               endif
               if (idspecies .gt. 1) then
                  Al(i,j,k)      = abs(Al(i,j,k)     )
                  S(i,j,k)       = abs(S(i,j,k)      )
                  Fe(i,j,k)      = abs(Fe(i,j,k)     )
               endif
            endif
!        endif
      enddo
      endif

      if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
      do i = is+1, ie+1
!        if (itmask_metal(i)) then
         if (idspecies .gt. 0) then
            MgSiO3(i,j,k)  = abs(MgSiO3(i,j,k) )
            AC(i,j,k)      = abs(AC(i,j,k)     )
         endif
         if (idspecies .gt. 1) then
            SiM(i,j,k)     = abs(SiM(i,j,k)    )
            FeM(i,j,k)     = abs(FeM(i,j,k)    )
            Mg2SiO4(i,j,k) = abs(Mg2SiO4(i,j,k))
            Fe3O4(i,j,k)   = abs(Fe3O4(i,j,k)  )
            SiO2D(i,j,k)   = abs(SiO2D(i,j,k)  )
            MgO(i,j,k)     = abs(MgO(i,j,k)    )
            FeS(i,j,k)     = abs(FeS(i,j,k)    )
            Al2O3(i,j,k)   = abs(Al2O3(i,j,k)  )
         endif
         if (idspecies .gt. 2) then
            reforg(i,j,k)  = abs(reforg(i,j,k)  )
            volorg(i,j,k)  = abs(volorg(i,j,k)  )
            H2Oice(i,j,k)  = abs(H2Oice(i,j,k)  )
         endif
!        endif
      enddo
      endif

      if (imchem .eq. 1) then
      do i = is+1, ie+1
!     if (itmask_metal(i)) then
c        if (d(i,j,k)*dom .lt. 1.e-2_DKIND) then
!!       if (d(i,j,k)*dom .lt. 
!!   &    min(1.e6_DKIND/(metal(i,j,k)/d(i,j,k)/0.02d-4)**2
!!   &       ,1.e6_DKIND)) then
         if ( ( (imetal .eq. 0)
     &    .and. (d(i,j,k)*dom .lt. 1.e8_DKIND) )
     &   .or. ( (imetal .eq. 1)
     &    .and. ( ( (metal(i,j,k) .le. 1.e-9_DKIND * d(i,j,k))
     &        .and. (d(i,j,k)*dom .lt. 1.e8_DKIND) )
     &       .or. ( (metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k))
     &        .and. (d(i,j,k)*dom .lt. 1.e6_DKIND) ) ) ) ) then

            totalOg = 16._DKIND/28._DKIND*   CO(i,j,k)
     &              + 32._DKIND/44._DKIND*  CO2(i,j,k)
     &              +                        OI(i,j,k)
     &              + 16._DKIND/17._DKIND*   OH(i,j,k)
     &              + 16._DKIND/18._DKIND*  H2O(i,j,k)
     &              +                        O2(i,j,k)
     &              + 16._DKIND/44._DKIND* SiOI(i,j,k)
     &              + 32._DKIND/60._DKIND*SiO2I(i,j,k)
     &              + 16._DKIND/28._DKIND* COII(i,j,k)
     &              +                       OII(i,j,k)
     &              + 16._DKIND/17._DKIND* OHII(i,j,k)
     &              + 16._DKIND/18._DKIND*H2OII(i,j,k)
     &              + 16._DKIND/19._DKIND*H3OII(i,j,k)
     &              +                      O2II(i,j,k)
            correctOg = real(Og(i)/totalOg, RKIND)
               CO(i,j,k) =    CO(i,j,k)*correctOg
              CO2(i,j,k) =   CO2(i,j,k)*correctOg
               OI(i,j,k) =    OI(i,j,k)*correctOg
               OH(i,j,k) =    OH(i,j,k)*correctOg
              H2O(i,j,k) =   H2O(i,j,k)*correctOg
               O2(i,j,k) =    O2(i,j,k)*correctOg
             SiOI(i,j,k) =  SiOI(i,j,k)*correctOg
            SiO2I(i,j,k) = SiO2I(i,j,k)*correctOg
             COII(i,j,k) =  COII(i,j,k)*correctOg
              OII(i,j,k) =   OII(i,j,k)*correctOg
             OHII(i,j,k) =  OHII(i,j,k)*correctOg
            H2OII(i,j,k) = H2OII(i,j,k)*correctOg
            H3OII(i,j,k) = H3OII(i,j,k)*correctOg
             O2II(i,j,k) =  O2II(i,j,k)*correctOg
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
         if (idspecies .gt. 0) then
            totalOd = 48._DKIND/100._DKIND* MgSiO3(i,j,k)
         endif
         if (idspecies .gt. 1) then
            totalOd = totalOd
     &              +  64._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 64._DKIND/232._DKIND*  Fe3O4(i,j,k)  
     &              + 32._DKIND/ 60._DKIND*  SiO2D(i,j,k)  
     &              + 16._DKIND/ 40._DKIND*    MgO(i,j,k)    
     &              + 48._DKIND/102._DKIND*  Al2O3(i,j,k)  
         endif
         if (idspecies .gt. 2) then
            totalOd = totalOd
     &              +  8._DKIND/22.68_DKIND*reforg(i,j,k)
     &              + 16._DKIND/32._DKIND  *volorg(i,j,k)
     &              + 16._DKIND/18._DKIND  *H2Oice(i,j,k)
         endif
            correctOd = real(Od(i)/totalOd, RKIND)
         if (idspecies .gt. 0) then
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctOd
         endif
         if (idspecies .gt. 1) then
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctOd
              Fe3O4(i,j,k) =   Fe3O4(i,j,k)*correctOd
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctOd
                MgO(i,j,k) =     MgO(i,j,k)*correctOd
              Al2O3(i,j,k) =   Al2O3(i,j,k)*correctOd
         endif
         if (idspecies .gt. 2) then
             reforg(i,j,k) =  reforg(i,j,k)*correctOd
             volorg(i,j,k) =  volorg(i,j,k)*correctOd
             H2Oice(i,j,k) =  H2Oice(i,j,k)*correctOd
         endif
         endif
  
            totalCg =                       CI(i,j,k)
     &              +                      CII(i,j,k)
     &              + 12._DKIND/28._DKIND*  CO(i,j,k)
     &              + 12._DKIND/44._DKIND* CO2(i,j,k)
     &              + 12._DKIND/13._DKIND*  CH(i,j,k)
     &              + 12._DKIND/14._DKIND* CH2(i,j,k)
     &              + 12._DKIND/28._DKIND*COII(i,j,k)
            correctCg = real(Cg(i)/totalCg, RKIND)
               CI(i,j,k) =   CI(i,j,k)*correctCg
              CII(i,j,k) =  CII(i,j,k)*correctCg
               CO(i,j,k) =   CO(i,j,k)*correctCg
              CO2(i,j,k) =  CO2(i,j,k)*correctCg
               CH(i,j,k) =   CH(i,j,k)*correctCg
              CH2(i,j,k) =  CH2(i,j,k)*correctCg
             COII(i,j,k) = COII(i,j,k)*correctCg
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
         if (idspecies .gt. 0) then
            totalCd =                           AC(i,j,k)
         endif
         if (idspecies .gt. 2) then
            totalCd = totalCd
     &              + 12._DKIND/22.68_DKIND*reforg(i,j,k)
     &              + 12._DKIND/32._DKIND  *volorg(i,j,k)
         endif
            correctCd = real(Cd(i)/totalCd, RKIND)
         if (idspecies .gt. 0) then
                AC(i,j,k) =     AC(i,j,k)*correctCd
         endif
         if (idspecies .gt. 2) then
            reforg(i,j,k) = reforg(i,j,k)*correctCd
            volorg(i,j,k) = volorg(i,j,k)*correctCd
         endif
         endif
            
            totalSig =                        SiI(i,j,k)
     &               + 28._DKIND/ 44._DKIND* SiOI(i,j,k)
     &               + 28._DKIND/ 60._DKIND*SiO2I(i,j,k)
            correctSig = real(Sig(i)/totalSig, RKIND)
              SiI(i,j,k) =   SiI(i,j,k)*correctSig
             SiOI(i,j,k) =  SiOI(i,j,k)*correctSig
            SiO2I(i,j,k) = SiO2I(i,j,k)*correctSig
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
         if (idspecies .gt. 0) then
            totalSid = 28._DKIND/100._DKIND* MgSiO3(i,j,k)
         endif
         if (idspecies .gt. 1) then
            totalSid = totalSid
     &               +                          SiM(i,j,k)
     &               + 28._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &               + 28._DKIND/ 60._DKIND*  SiO2D(i,j,k)
         endif
            correctSid = real(Sid(i)/totalSid, RKIND)
         if (idspecies .gt. 0) then
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctSid
         endif
         if (idspecies .gt. 1) then
                SiM(i,j,k) =     SiM(i,j,k)*correctSid
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctSid
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctSid
         endif
         endif

         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
         if (idspecies .gt. 1) then
            totalFeg = Fe(i,j,k)
            correctFeg = real(Feg(i)/totalFeg, RKIND)
            Fe(i,j,k) =    Fe(i,j,k)*correctFeg

            totalFed =                        FeM(i,j,k)
     &               +168._DKIND/232._DKIND*Fe3O4(i,j,k)
     &               + 56._DKIND/ 88._DKIND*  FeS(i,j,k)
            correctFed = real(Fed(i)/totalFed, RKIND)
              FeM(i,j,k) =   FeM(i,j,k)*correctFed
            Fe3O4(i,j,k) = Fe3O4(i,j,k)*correctFed
              FeS(i,j,k) =   FeS(i,j,k)*correctFed
         endif

         if (idspecies .gt. 0) then
            totalMgg =                   Mg(i,j,k)
            correctMgg = real( Mgg(i)/totalMgg, RKIND)
            Mg(i,j,k)      = Mg(i,j,k)     *correctMgg
            totalMgd = 24._DKIND/100._DKIND* MgSiO3(i,j,k) 
            if (idspecies .gt. 1) then
            totalMgd = totalMgd
     &               + 48._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &               + 24._DKIND/ 40._DKIND*    MgO(i,j,k)    
            endif
            correctMgd = real( Mgd(i)/totalMgd, RKIND)
            MgSiO3(i,j,k)  = MgSiO3(i,j,k) *correctMgd
            if (idspecies .gt. 1) then
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctMgd
            MgO(i,j,k)     = MgO(i,j,k)    *correctMgd
            endif
         endif

         if (idspecies .gt. 1) then
              S(i,j,k) =                         Sg(i)
            FeS(i,j,k) = 88._DKIND / 32._DKIND * Sd(i)

               Al(i,j,k) =                         Alg(i)
            Al2O3(i,j,k) =  102._DKIND/54._DKIND * Ald(i)
         endif
         endif

         else ! d(i,j,k)

            totalO  = 16._DKIND/28._DKIND*   CO(i,j,k)
     &              + 32._DKIND/44._DKIND*  CO2(i,j,k)
     &              +                        OI(i,j,k)
     &              + 16._DKIND/17._DKIND*   OH(i,j,k)
     &              + 16._DKIND/18._DKIND*  H2O(i,j,k)
     &              +                        O2(i,j,k)
     &              + 16._DKIND/44._DKIND* SiOI(i,j,k)
     &              + 32._DKIND/60._DKIND*SiO2I(i,j,k)
     &              + 16._DKIND/28._DKIND* COII(i,j,k)
     &              +                       OII(i,j,k)
     &              + 16._DKIND/17._DKIND* OHII(i,j,k)
     &              + 16._DKIND/18._DKIND*H2OII(i,j,k)
     &              + 16._DKIND/19._DKIND*H3OII(i,j,k)
     &              +                      O2II(i,j,k)
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
         if (idspecies .gt. 0) then
            totalO  = totalO
     &              + 48._DKIND/100._DKIND* MgSiO3(i,j,k)
         endif
         if (idspecies .gt. 1) then
            totalO  = totalO
     &              + 64._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 64._DKIND/232._DKIND*  Fe3O4(i,j,k)  
     &              + 32._DKIND/ 60._DKIND*  SiO2D(i,j,k)  
     &              + 16._DKIND/ 40._DKIND*    MgO(i,j,k)    
     &              + 48._DKIND/102._DKIND*  Al2O3(i,j,k)  
         endif
         if (idspecies .gt. 2) then
            totalO  = totalO
     &              +  8._DKIND/22.68_DKIND*reforg(i,j,k)
     &              + 16._DKIND/32._DKIND  *volorg(i,j,k)
     &              + 16._DKIND/18._DKIND  *H2Oice(i,j,k)
         endif
         endif
         if ( ( igrgr .eq. 0 ) .and. ( idsub .eq. 0) ) then
            correctO  = real(Og(i)/ totalO, RKIND)
               CO(i,j,k)     =    CO(i,j,k)*correctO
              CO2(i,j,k)     =   CO2(i,j,k)*correctO
               OI(i,j,k)     =    OI(i,j,k)*correctO
               OH(i,j,k)     =    OH(i,j,k)*correctO
              H2O(i,j,k)     =   H2O(i,j,k)*correctO
               O2(i,j,k)     =    O2(i,j,k)*correctO
             SiOI(i,j,k)     =  SiOI(i,j,k)*correctO
            SiO2I(i,j,k)     = SiO2I(i,j,k)*correctO
             COII(i,j,k)     =  COII(i,j,k)*correctO
              OII(i,j,k)     =   OII(i,j,k)*correctO
             OHII(i,j,k)     =  OHII(i,j,k)*correctO
            H2OII(i,j,k)     = H2OII(i,j,k)*correctO
            H3OII(i,j,k)     = H3OII(i,j,k)*correctO
             O2II(i,j,k)     =  O2II(i,j,k)*correctO
         else
            correctO  = real(Ot(i)/ totalO, RKIND)
               CO(i,j,k)     =    CO(i,j,k)*correctO
              CO2(i,j,k)     =   CO2(i,j,k)*correctO
               OI(i,j,k)     =    OI(i,j,k)*correctO
               OH(i,j,k)     =    OH(i,j,k)*correctO
              H2O(i,j,k)     =   H2O(i,j,k)*correctO
               O2(i,j,k)     =    O2(i,j,k)*correctO
             SiOI(i,j,k)     =  SiOI(i,j,k)*correctO
            SiO2I(i,j,k)     = SiO2I(i,j,k)*correctO
             COII(i,j,k)     =  COII(i,j,k)*correctO
              OII(i,j,k)     =   OII(i,j,k)*correctO
             OHII(i,j,k)     =  OHII(i,j,k)*correctO
            H2OII(i,j,k)     = H2OII(i,j,k)*correctO
            H3OII(i,j,k)     = H3OII(i,j,k)*correctO
             O2II(i,j,k)     =  O2II(i,j,k)*correctO
         if (idspecies .gt. 0) then
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctO 
         endif
         if (idspecies .gt. 1) then
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctO 
              Fe3O4(i,j,k) =   Fe3O4(i,j,k)*correctO 
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctO 
                MgO(i,j,k) =     MgO(i,j,k)*correctO 
              Al2O3(i,j,k) =   Al2O3(i,j,k)*correctO 
         endif
         if (idspecies .gt. 2) then
             reforg(i,j,k) =  reforg(i,j,k)*correctO 
             volorg(i,j,k) =  volorg(i,j,k)*correctO 
             H2Oice(i,j,k) =  H2Oice(i,j,k)*correctO 
         endif
         endif
  
            totalC  =                           CI(i,j,k)
     &              +                          CII(i,j,k)
     &              + 12._DKIND/28._DKIND    *  CO(i,j,k)
     &              + 12._DKIND/44._DKIND    * CO2(i,j,k)
     &              + 12._DKIND/13._DKIND    *  CH(i,j,k)
     &              + 12._DKIND/14._DKIND    * CH2(i,j,k)
     &              + 12._DKIND/28._DKIND    *COII(i,j,k)
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
         if (idspecies .gt. 0) then
            totalC  = totalC
     &              +                           AC(i,j,k)
         endif
         if (idspecies .gt. 2) then
            totalC  = totalC
     &              + 12._DKIND/22.68_DKIND*reforg(i,j,k)
     &              + 12._DKIND/32._DKIND  *volorg(i,j,k)
         endif
         endif
         if ( ( igrgr .eq. 0 ) .and. ( idsub .eq. 0 ) ) then
            correctC = real(Cg(i)/ totalC, RKIND)
               CI(i,j,k) =   CI(i,j,k)*correctC
              CII(i,j,k) =  CII(i,j,k)*correctC
               CO(i,j,k) =   CO(i,j,k)*correctC
              CO2(i,j,k) =  CO2(i,j,k)*correctC
               CH(i,j,k) =   CH(i,j,k)*correctC
              CH2(i,j,k) =  CH2(i,j,k)*correctC
             COII(i,j,k) = COII(i,j,k)*correctC
         else
            correctC  = real(  Ct(i)/ totalC, RKIND)
               CI(i,j,k) =     CI(i,j,k)*correctC
              CII(i,j,k) =    CII(i,j,k)*correctC
               CO(i,j,k) =     CO(i,j,k)*correctC
              CO2(i,j,k) =    CO2(i,j,k)*correctC
               CH(i,j,k) =     CH(i,j,k)*correctC
              CH2(i,j,k) =    CH2(i,j,k)*correctC
             COII(i,j,k) =   COII(i,j,k)*correctC
         if (idspecies .gt. 0) then
               AC(i,j,k) =     AC(i,j,k)*correctC
         endif
         if (idspecies .gt. 2) then
           reforg(i,j,k) = reforg(i,j,k)*correctC
           volorg(i,j,k) = volorg(i,j,k)*correctC
         endif
         endif
            
            totalSi =                       SiI(i,j,k)
     &              + 28._DKIND/ 44._DKIND*SiOI(i,j,k)
     &              + 28._DKIND/ 60._DKIND*SiO2I(i,j,k)
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
         if (idspecies .gt. 0) then
            totalSi = totalSi
     &              + 28._DKIND/100._DKIND* MgSiO3(i,j,k)
         endif
         if (idspecies .gt. 1) then
            totalSi = totalSi
     &              +                          SiM(i,j,k)
     &              + 28._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 28._DKIND/ 60._DKIND*  SiO2D(i,j,k)
         endif
         endif
         if ( ( igrgr .eq. 0 ) .and. ( idsub .eq. 0 ) ) then
            correctSi = real(Sig(i)/totalSi, RKIND)
                SiI(i,j,k) =   SiI(i,j,k)*correctSi
               SiOI(i,j,k) =  SiOI(i,j,k)*correctSi
              SiO2I(i,j,k) = SiO2I(i,j,k)*correctSi
         else
            correctSi = real(Sit(i)/totalSi, RKIND)
                SiI(i,j,k) =   SiI(i,j,k)*correctSi
               SiOI(i,j,k) =  SiOI(i,j,k)*correctSi
              SiO2I(i,j,k) = SiO2I(i,j,k)*correctSi
         if (idspecies .gt. 0) then
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctSi
         endif
         if (idspecies .gt. 1) then
                SiM(i,j,k) =     SiM(i,j,k)*correctSi
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctSi
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctSi
         endif
         endif

         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
         if (idspecies .gt. 1) then
            totalFe =                         Fe(i,j,k)     
     &              +                        FeM(i,j,k)  
     &              +168._DKIND/232._DKIND*Fe3O4(i,j,k)
     &              + 56._DKIND/ 88._DKIND*  FeS(i,j,k)  
            correctFe = real( Fet(i)/totalFe, RKIND)
            Fe(i,j,k)    =    Fe(i,j,k)*correctFe
            FeM(i,j,k)   =   FeM(i,j,k)*correctFe
            Fe3O4(i,j,k) = Fe3O4(i,j,k)*correctFe
            FeS(i,j,k)   =   FeS(i,j,k)*correctFe
         endif

         if (idspecies .gt. 0) then
            totalMg =                           Mg(i,j,k)     
     &              + 24._DKIND/100._DKIND* MgSiO3(i,j,k) 
            if (idspecies .gt. 1) then
            totalMg = totalMg
     &              + 48._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 24._DKIND/ 40._DKIND*    MgO(i,j,k)    
            endif
            correctMg = real( Mgt(i)/totalMg, RKIND)
                 Mg(i,j,k) =      Mg(i,j,k)*correctMg
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctMg
            if (idspecies .gt. 1) then
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctMg
                MgO(i,j,k) =     MgO(i,j,k)*correctMg
            endif
         endif

         if (idspecies .gt. 1) then
            totalS  =                        S(i,j,k)  
     &              + 32._DKIND/ 88._DKIND*FeS(i,j,k)
            correctS  = real(  St(i)/totalS , RKIND)
              S(i,j,k) =   S(i,j,k)*correctS
            FeS(i,j,k) = FeS(i,j,k)*correctS

            totalAl =                         Al(i,j,k)  
     &              + 54._DKIND/102._DKIND*Al2O3(i,j,k)
            correctAl = real( Alt(i)/totalAl, RKIND)
               Al(i,j,k) =    Al(i,j,k)*correctAl
            Al2O3(i,j,k) = Al2O3(i,j,k)*correctAl
         endif
         endif

         endif ! d(i,j,k)

!           CI(i,j,k)      = max(CI(i,j,k), tiny)
!           CII(i,j,k)     = max(CII(i,j,k), tiny)
!           CO(i,j,k)      = max(CO(i,j,k), tiny)
!           CO2(i,j,k)     = max(CO2(i,j,k), tiny)
!           OI(i,j,k)      = max(OI(i,j,k), tiny)
!           OH(i,j,k)      = max(OH(i,j,k), tiny)
!           H2O(i,j,k)     = max(H2O(i,j,k), tiny)
!           O2(i,j,k)      = max(O2(i,j,k), tiny)
!           SiI(i,j,k)     = max(SiI(i,j,k), tiny)
!           SiOI(i,j,k)    = max(SiOI(i,j,k), tiny)
!           SiO2I(i,j,k)   = max(SiO2I(i,j,k), tiny)
!           CH(i,j,k)      = max(CH(i,j,k), tiny)
!           CH2(i,j,k)     = max(CH2(i,j,k), tiny)
!           COII(i,j,k)    = max(COII(i,j,k), tiny)
!           OII(i,j,k)     = max(OII(i,j,k), tiny)
!           OHII(i,j,k)    = max(OHII(i,j,k), tiny)
!           H2OII(i,j,k)   = max(H2OII(i,j,k), tiny)
!           H3OII(i,j,k)   = max(H3OII(i,j,k), tiny)
!           O2II(i,j,k)    = max(O2II(i,j,k), tiny)
!        if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
!        if (idspecies .gt. 0) then
!           Mg(i,j,k)      = max(Mg(i,j,k), tiny)
!           MgSiO3(i,j,k)  = max(MgSiO3(i,j,k), tiny)
!           AC(i,j,k)      = max(AC(i,j,k), tiny)
!        endif
!        if (idspecies .gt. 1) then
!           Al(i,j,k)      = max(Al(i,j,k), tiny)
!           S(i,j,k)       = max(S(i,j,k), tiny)
!           Fe(i,j,k)      = max(Fe(i,j,k), tiny)
!           SiM(i,j,k)     = max(SiM(i,j,k), tiny)
!           FeM(i,j,k)     = max(FeM(i,j,k), tiny)
!           Mg2SiO4(i,j,k) = max(Mg2SiO4(i,j,k), tiny)
!           Fe3O4(i,j,k)   = max(Fe3O4(i,j,k), tiny)
!           SiO2D(i,j,k)   = max(SiO2D(i,j,k), tiny)
!           MgO(i,j,k)     = max(MgO(i,j,k), tiny)
!           FeS(i,j,k)     = max(FeS(i,j,k), tiny)
!           Al2O3(i,j,k)   = max(Al2O3(i,j,k), tiny)
!        endif
!        if (idspecies .gt. 2) then
!           reforg(i,j,k)  = max(reforg(i,j,k), tiny)
!           volorg(i,j,k)  = max(volorg(i,j,k), tiny)
!           H2Oice(i,j,k)  = max(H2Oice(i,j,k), tiny)
!        endif
!        endif

!     endif
      enddo
      endif ! imchem

C     if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then
C        do i = is+1, ie+1
C!          if ( itmask_metal(i) ) then
C           if (idspecies .gt. 0) then
C              dust(i,j,k) = MgSiO3  (i,j,k)
C    &                     + AC      (i,j,k)
C           endif
C           if (idspecies .gt. 1) then
C              dust(i,j,k) = dust(i,j,k)
C    &                     + SiM     (i,j,k)
C    &                     + FeM     (i,j,k)
C    &                     + Mg2SiO4 (i,j,k)
C    &                     + Fe3O4   (i,j,k)
C    &                     + SiO2D   (i,j,k)
C    &                     + MgO     (i,j,k)
C    &                     + FeS     (i,j,k)
C    &                     + Al2O3   (i,j,k)
C           endif
C           if (idspecies .gt. 2) then
C              dust(i,j,k) = dust(i,j,k)
C    &                     + reforg  (i,j,k)
C    &                     + volorg  (i,j,k)
C    &                     + H2Oice  (i,j,k)
C           endif
C!          endif
C        enddo
C     endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND +
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &        - HM(i,j,k) + H2II(i,j,k)/2._RKIND
         if (ispecies .gt. 3) de(i,j,k) = de(i,j,k)
     &        - DM   (i,j,k)/2._RKIND
     &        + HDII (i,j,k)/3._RKIND
     &        + HeHII(i,j,k)/5._RKIND
         if (imchem .eq. 1) then
!        if (itmask_metal(i)) then
            de(i,j,k) = de(i,j,k)
     &        + CII  (i,j,k)/12._RKIND
     &        + COII (i,j,k)/28._RKIND
     &        + OII  (i,j,k)/16._RKIND
     &        + OHII (i,j,k)/17._RKIND
     &        + H2OII(i,j,k)/18._RKIND
     &        + H3OII(i,j,k)/19._RKIND
     &        + O2II (i,j,k)/32._RKIND
!        endif
         endif
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k

      return
      end
