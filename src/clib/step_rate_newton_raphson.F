      subroutine step_rate_newton_raphson(icool, d, e, u, v, w, de, HI,
     &                HII, HeI, HeII, HeIII, in, jn, kn, nratec,
     &                iexpand, ispecies, imetal, imcool, idust,
     &                idustall, idustfield, is, ie, ih2co, ipiht,
     &                idustrec, igammah, aye, temstart, temend, utem,
     &                uxyz, uaye, urho, utim, gamma, fh, z_solar, fgr,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a,
     &                k10a, k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a, k24, k25, k26, k27, k28,
     &                k29, k30, k31, k50a, k51a, k52a, k53a, k54a,
     &                k55a, k56a, k57a, k58a, ndratec, dtemstart,
     &                dtemend, h2dusta, ncrna, ncrd1a, ncrd2a, ceHIa,
     &                ceHeIa, ceHeIIa, ciHIa, ciHeIa, ciHeISa, ciHeIIa,
     &                reHIIa, reHeII1a, reHeII2a, reHeIIIa, brema,
     &                compa, gammaha, isrf, regra, gamma_isrfa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII, HM,
     &                H2I, H2II, DI, DII, HDI, metal, dust, hyd01ka,
     &                h2k01a, vibha, rotha, rotla, gpldla, gphdla,
     &                hdltea, hdlowa, gaHIa, gaH2a, gaHea, gaHpa,
     &                gaela, h2ltea, gasgra, iH2shield, iradshield,
     &                avgsighi, avgsighei, avgsigheii, iradtrans,
     &                irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                photogamma, xH2shield, ih2optical, iciecool,
     &                ithreebody, ih2cr, ihdcr, ciecoa, icmbTfloor,
     &                iClHeat, clEleFra, priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim, metPar1, metPar2,
     &                metPar3, metPar4, metPar5, metDataSize,
     &                metCooling, metHeating, clnew, iVheat, iMheat,
     &                Vheat, Mheat, iTfloor, Tfloor_scalar, Tfloor,
     &                imchem, igrgr, ipcont, tmcool, DM, HDII, HeHII,
     &                CI, CII, CO, CO2, OI, OH, H2O, O2, SiI, SiOI,
     &                SiO2I, CH, CH2, COII, OII, OHII, H2OII, H3OII,
     &                O2II, Mg, Al, S, Fe, SiM, FeM, Mg2SiO4, MgSiO3,
     &                Fe3O4, AC, SiO2D, MgO, FeS, Al2O3, reforg,
     &                volorg, H2Oice, k125a, k129a, k130a, k131a,
     &                k132a, k133a, k134a, k135a, k136a, k137a, k148a,
     &                k149a, k150a, k151a, k152a, k153a, kz15a, kz16a,
     &                kz17a, kz18a, kz19a, kz20a, kz21a, kz22a, kz23a,
     &                kz24a, kz25a, kz26a, kz27a, kz28a, kz29a, kz30a,
     &                kz31a, kz32a, kz33a, kz34a, kz35a, kz36a, kz37a,
     &                kz38a, kz39a, kz40a, kz41a, kz42a, kz43a, kz44a,
     &                kz45a, kz46a, kz47a, kz48a, kz49a, kz50a, kz51a,
     &                kz52a, kz53a, kz54a, cieY06a, LH2_N, LH2_Size,
     &                LH2_D, LH2_T, LH2_H, LH2_dD, LH2_dT, LH2_dH,
     &                LH2_L, LHD_N, LHD_Size, LHD_D, LHD_T, LHD_H,
     &                LHD_dD, LHD_dT, LHD_dH, LHD_L, LCI_N, LCI_Size,
     &                LCI_D, LCI_T, LCI_H, LCI_dD, LCI_dT, LCI_dH,
     &                LCI_L, LCII_N, LCII_Size, LCII_D, LCII_T, LCII_H,
     &                LCII_dD, LCII_dT, LCII_dH, LCII_L, LOI_N,
     &                LOI_Size, LOI_D, LOI_T, LOI_H, LOI_dD, LOI_dT,
     &                LOI_dH, LOI_L, LCO_N, LCO_Size, LCO_D, LCO_T,
     &                LCO_H, LCO_dD, LCO_dT, LCO_dH, LCO_L, LOH_N,
     &                LOH_Size, LOH_D, LOH_T, LOH_H, LOH_dD, LOH_dT,
     &                LOH_dH, LOH_L, LH2O_N, LH2O_Size, LH2O_D, LH2O_T,
     &                LH2O_H, LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L,
     &                alphap_N, alphap_Size, alphap_D, alphap_T,
     &                alphap_dD, alphap_dT, alphap_Data, immulti,
     &                imabund, idspecies, itdmulti, idsub, metal_loc,
     &                metal_C13, metal_C20, metal_C25, metal_C30,
     &                metal_F13, metal_F15, metal_F50, metal_F80,
     &                metal_P170, metal_P200, metal_Y19, SN0_N,
     &                SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3,
     &                SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO,
     &                SN0_fFeS, SN0_fAl2O3, SN0_freforg, SN0_fvolorg,
     &                SN0_fH2Oice, SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4,
     &                SN0_r0MgSiO3, SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D,
     &                SN0_r0MgO, SN0_r0FeS, SN0_r0Al2O3, SN0_r0reforg,
     &                SN0_r0volorg, SN0_r0H2Oice, gr_N, gr_Size, gr_dT,
     &                gr_Td, SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4,
     &                SN0_kpMgSiO3, SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D,
     &                SN0_kpMgO, SN0_kpFeS, SN0_kpAl2O3, SN0_kpreforg,
     &                SN0_kpvolorg, SN0_kpH2Oice, h2dustSa, h2dustCa,
     &                gasgr2a, gamma_isrf2a, grogra, idissHDI,
     &                kdissHDI, iionZ, kphCI, kphOI, idissZ, kdissCO,
     &                kdissOH, kdissH2O, iuseH2shield, iisrffield,
     &                isrf_habing, iH2shieldcustom, f_shield_custom,
     &                ierror, j, k, iter, dom, comp1,
     &                comp2, coolunit, tbase1, xbase1, chunit, dx_cgs,
     &                c_ljeans, indixe, t1, t2, logtem, tdef, dtit,
     &                p2d, tgas, tgasold, tdust, metallicity, dust2gas,
     &                rhoH, mmw, mynh, myde, gammaha_eff, gasgr_tdust,
     &                regr, h2dust, ncrn, ncrd1, ncrd2, tSiM, tFeM,
     &                tMg2SiO4, tMgSiO3, tFe3O4, tAC, tSiO2D, tMgO,
     &                tFeS, tAl2O3, treforg, tvolorg, tH2Oice, ceHI,
     &                ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem, edot,
     &                hyd01k, h2k01, vibh, roth, rotl, gpldl, gphdl,
     &                hdlte, hdlow, cieco, anydust, itmask_nr,
     &                itmask_metal, itr, imp_eng
     &              )

!  PURPOSE:
!    An alternative to step_rate_g that employs Newton-Raphson
!    rather than using Gauss-Seidel

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer icool, in, jn, kn, is, ie, nratec, iexpand, ih2co, ipiht,
     &        ispecies, imetal, ierror, imcool, idust, idustall,
     &        idustfield, idustrec, igammah, ih2optical, iciecool,
     &        ithreebody, ih2cr, ihdcr, ndratec, clnew, iVheat, iMheat,
     &        iTfloor, iH2shield, iradshield, iradtrans, irt_honly,
     &        imchem, igrgr, ipcont, iisrffield, iH2shieldcustom

      real*8 aye, temstart, temend, gamma, utim, uxyz, uaye, urho,
     &       utem, fh, z_solar, fgr, dtemstart, dtemend, clEleFra,
     &       Tfloor_scalar, tmcool

!  Density, energy and velocity fields fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        d(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn), dust(in,jn,kn),
     &        Vheat(in,jn,kn), Mheat(in,jn,kn), Tfloor(in,jn,kn)
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

!  Radiative transfer fields

      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn), photogamma(in,jn,kn)

      integer idissHDI, iionZ, idissZ
      R_PREC  kdissHDI(in,jn,kn), kphCI(in,jn,kn), kphOI(in,jn,kn),
     &        kdissCO(in,jn,kn), kdissOH(in,jn,kn), kdissH2O(in,jn,kn)
      integer iuseH2shield

!  H2 self-shielding length-scale field

      R_PREC  xH2shield(in,jn,kn)

!  Interstellar radiation field for dust heating

      R_PREC  isrf_habing(in,jn,kn)

!  Custom H2 shielding factor

      R_PREC f_shield_custom(in, jn, kn)

!  Cooling tables (coolings rates as a function of temperature)

      real*8  hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     &        gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), h2ltea(nratec),
     &        gasgra(nratec), ciecoa(nratec),
     &        ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec),
     &        compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        gammaha, isrf, regra(nratec), gamma_isrfa
       real*8  cieY06a(nratec)
       integer LH2_N(3), LH2_Size
       real*8  LH2_D(LH2_N(1)), LH2_T(LH2_N(2)), LH2_H(LH2_N(3))
     &       , LH2_dD, LH2_dT, LH2_dH, LH2_L(LH2_Size)
       integer LHD_N(3), LHD_Size
       real*8  LHD_D(LHD_N(1)), LHD_T(LHD_N(2)), LHD_H(LHD_N(3))
     &       , LHD_dD, LHD_dT, LHD_dH, LHD_L(LHD_Size)
       integer LCI_N(3), LCI_Size
       real*8  LCI_D(LCI_N(1)), LCI_T(LCI_N(2)), LCI_H(LCI_N(3))
     &       , LCI_dD, LCI_dT, LCI_dH, LCI_L(LCI_Size)
       integer LCII_N(3), LCII_Size
       real*8  LCII_D(LCII_N(1)), LCII_T(LCII_N(2)), LCII_H(LCII_N(3))
     &       , LCII_dD, LCII_dT, LCII_dH, LCII_L(LCII_Size)
       integer LOI_N(3), LOI_Size
       real*8  LOI_D(LOI_N(1)), LOI_T(LOI_N(2)), LOI_H(LOI_N(3))
     &       , LOI_dD, LOI_dT, LOI_dH, LOI_L(LOI_Size)
       integer LCO_N(3), LCO_Size
       real*8  LCO_D(LCO_N(1)), LCO_T(LCO_N(2)), LCO_H(LCO_N(3))
     &       , LCO_dD, LCO_dT, LCO_dH, LCO_L(LCO_Size)
       integer LOH_N(3), LOH_Size
       real*8  LOH_D(LOH_N(1)), LOH_T(LOH_N(2)), LOH_H(LOH_N(3))
     &       , LOH_dD, LOH_dT, LOH_dH, LOH_L(LOH_Size)
       integer LH2O_N(3), LH2O_Size
       real*8  LH2O_D(LH2O_N(1)), LH2O_T(LH2O_N(2)), LH2O_H(LH2O_N(3))
     &       , LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L(LH2O_Size)
       integer alphap_N(2), alphap_Size
       real*8  alphap_D(alphap_N(1)), alphap_T(alphap_N(2))
     &       , alphap_dD, alphap_dT
     &       , alphap_Data(alphap_Size)
      integer immulti, imabund, idspecies, itdmulti, idsub
      integer SN0_N
      real*8  SN0_fSiM(SN0_N), SN0_fFeM(SN0_N), SN0_fMg2SiO4(SN0_N)
     &      , SN0_fMgSiO3(SN0_N), SN0_fFe3O4(SN0_N), SN0_fAC(SN0_N)
     &      , SN0_fSiO2D(SN0_N), SN0_fMgO(SN0_N), SN0_fFeS(SN0_N)
     &      , SN0_fAl2O3(SN0_N)
     &      , SN0_freforg(SN0_N), SN0_fvolorg(SN0_N), SN0_fH2Oice(SN0_N)
      real*8  SN0_r0SiM(3,SN0_N), SN0_r0FeM(3,SN0_N)
     &      , SN0_r0Mg2SiO4(3,SN0_N), SN0_r0MgSiO3(3,SN0_N)
     &      , SN0_r0Fe3O4(3,SN0_N), SN0_r0AC(3,SN0_N)
     &      , SN0_r0SiO2D(3,SN0_N), SN0_r0MgO(3,SN0_N)
     &      , SN0_r0FeS(3,SN0_N), SN0_r0Al2O3(3,SN0_N)
     &      , SN0_r0reforg(3,SN0_N)
     &      , SN0_r0volorg(3,SN0_N), SN0_r0H2Oice(3,SN0_N)
!     opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  SN0_kpSiM(gr_Size,SN0_N), SN0_kpFeM(gr_Size,SN0_N)
     &      , SN0_kpMg2SiO4(gr_Size,SN0_N), SN0_kpMgSiO3(gr_Size,SN0_N)
     &      , SN0_kpFe3O4(gr_Size,SN0_N), SN0_kpAC(gr_Size,SN0_N)
     &      , SN0_kpSiO2D(gr_Size,SN0_N), SN0_kpMgO(gr_Size,SN0_N)
     &      , SN0_kpFeS(gr_Size,SN0_N), SN0_kpAl2O3(gr_Size,SN0_N)
     &      , SN0_kpreforg(gr_Size,SN0_N)
     &      , SN0_kpvolorg(gr_Size,SN0_N), SN0_kpH2Oice(gr_Size,SN0_N)
      real*8  gasgr2a(nratec), gamma_isrf2a

      real*8  avgsighi, avgsighei, avgsigheii

!  Chemistry tables (rates as a function of temperature)

      real*8 k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &       k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &       k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &       k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &       k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &       k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &       k54a(nratec), k55a(nratec), k56a(nratec),
     &       k57a(nratec), k58a(nratec),
     &       k13dda(nratec, 14), h2dusta(nratec, ndratec),
     &       ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &       k24, k25, k26, k27, k28, k29, k30, k31
       real*8  k125a(nratec), k129a(nratec), k130a(nratec)
     &       , k131a(nratec), k132a(nratec), k133a(nratec)
     &       , k134a(nratec), k135a(nratec), k136a(nratec)
     &       , k137a(nratec), k148a(nratec), k149a(nratec)
     &       , k150a(nratec), k151a(nratec), k152a(nratec)
     &       , k153a(nratec)
     &       , kz15a(nratec), kz16a(nratec), kz17a(nratec)
     &       , kz18a(nratec), kz19a(nratec), kz20a(nratec)
     &       , kz21a(nratec), kz22a(nratec), kz23a(nratec)
     &       , kz24a(nratec), kz25a(nratec), kz26a(nratec)
     &       , kz27a(nratec), kz28a(nratec), kz29a(nratec)
     &       , kz30a(nratec), kz31a(nratec), kz32a(nratec)
     &       , kz33a(nratec), kz34a(nratec), kz35a(nratec)
     &       , kz36a(nratec), kz37a(nratec), kz38a(nratec)
     &       , kz39a(nratec), kz40a(nratec), kz41a(nratec)
     &       , kz42a(nratec), kz43a(nratec), kz44a(nratec)
     &       , kz45a(nratec), kz46a(nratec), kz47a(nratec)
     &       , kz48a(nratec), kz49a(nratec), kz50a(nratec)
     &       , kz51a(nratec), kz52a(nratec), kz53a(nratec)
     &       , kz54a(nratec)

!  Cloudy cooling data

      integer icmbTfloor, iClHeat
      integer*8 priGridRank, priDataSize,
     &     metGridRank, metDataSize,
     &     priGridDim(priGridRank), metGridDim(metGridRank)
      real*8 priPar1(priGridDim(1)), priPar2(priGridDim(2)), 
     &     priPar3(priGridDim(3)), priPar4(priGridDim(4)),
     &     priPar5(priGridDim(5)),
     &     metPar1(metGridDim(1)), metPar2(metGridDim(2)), 
     &     metPar3(metGridDim(3)), metPar4(metGridDim(4)),
     &     metPar5(metGridDim(5)),
     &     priCooling(priDataSize), priHeating(priDataSize),
     &     priMMW(priDataSize),
     &     metCooling(metDataSize), metHeating(metDataSize)

!  Local variable
      integer i

!  Assorted arguments
      integer j, k, iter
      real*8 dom, comp1, comp2
      real*8 coolunit, tbase1, xbase1, chunit
      real*8 dx_cgs, c_ljeans

!  Assorted row buffers

      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in), dtit(in), p2d(in),
     &       tgas(in), tgasold(in), tdust(in), metallicity(in),
     &       dust2gas(in), rhoH(in), mmw(in), mynh(in), myde(in),
     &       gammaha_eff(in), gasgr_tdust(in), regr(in)

!  Rate equation row buffers

      real*8 h2dust(in), ncrn(in), ncrd1(in), ncrd2(in)

      real*8  h2dustSa(nratec*ndratec), h2dustCa(nratec*ndratec)
      real*8  grogra(nratec)
!     grain temperature
      real*8  tSiM(in), tFeM(in), tMg2SiO4(in)
     &      , tMgSiO3(in), tFe3O4(in), tAC(in)
     &      , tSiO2D(in), tMgO(in), tFeS(in)
     &      , tAl2O3(in)
     &      , treforg(in), tvolorg(in), tH2Oice(in)

!  Cooling/heating row buffers

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), edot(in)
      real*8 hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in), cieco(in)

!  Iteration mask

      MASK_TYPE anydust
      MASK_TYPE itmask_nr(in)
      MASK_TYPE itmask_metal(in)
      integer itr, imp_eng(in)
!  Local variable
      integer itr_time
      integer nsp, isp, jsp, id
      real*8  dspj, err, err_max
      integer,parameter   :: i_eng = 52
! There may be an argument for allocating the following at a higher
! level function, but we will leave that for after transcription
      real*8  dsp(i_eng), dsp0(i_eng), dsp1(i_eng)
     &      , dspdot(i_eng), dspdot1(i_eng)
     &      , ddsp(i_eng), der(i_eng,i_eng)

! (In the future, we may want to reconsider when/how we allocate
! the following 3 variables)
      integer,allocatable :: idsp(:)
      real*8, allocatable :: mtrx(:,:), vec(:)
! Another parameter
      real*8, parameter   :: eps = 1.d-4

! The following was extracted from another subroutine
            do i = is+1, ie+1
            if (itmask_nr(i) .ne. MASK_FALSE) then

!           If density and temperature are low, update gas energy explicitly

               if (icool .eq. 1) then
                  if (imp_eng(i) .eq. 0) then
                     e(i,j,k)  = e(i,j,k) +
     &                    real(edot(i)/d(i,j,k)*dtit(i), RKIND)
                  endif
               endif

!              initialize arrays
               if (ispecies .gt. 0) nsp = 6
               if (ispecies .gt. 1) nsp = nsp + 3
               if (ispecies .gt. 2) nsp = nsp + 3
               if (ispecies .gt. 3) nsp = nsp + 3
               if (itmask_metal(i) .ne. MASK_FALSE) then
                  if (imchem .eq. 1) then
                     nsp = nsp + 19
                     if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
                        if (idspecies .gt. 0) nsp = nsp + 1
                        if (idspecies .gt. 1) nsp = nsp + 3
                     endif
                  endif
                  if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
                     if (idspecies .gt. 0) nsp = nsp + 2
                     if (idspecies .gt. 1) nsp = nsp + 8
                     if (idspecies .gt. 2) nsp = nsp + 3
                  endif
               endif
               nsp = nsp + imp_eng(i)
               allocate(idsp(nsp))
               allocate(mtrx(nsp,nsp))
               allocate(vec(nsp))

               if ( ispecies .gt. 0 ) then
                  dsp( 1) = de(i,j,k)
                  dsp( 2) = HI(i,j,k)
                  dsp( 3) = HII(i,j,k)
                  dsp( 4) = HeI(i,j,k)
                  dsp( 5) = HeII(i,j,k)
                  dsp( 6) = HeIII(i,j,k)
               endif
               if ( ispecies .gt. 1 ) then
                  dsp( 7) = HM(i,j,k)
                  dsp( 8) = H2I(i,j,k)
                  dsp( 9) = H2II(i,j,k)
               endif
               if ( ispecies .gt. 2 ) then
                  dsp(10) = DI(i,j,k)
                  dsp(11) = DII(i,j,k)
                  dsp(12) = HDI(i,j,k)
               endif
               if ( ispecies .gt. 3 ) then
                  dsp(13) = DM(i,j,k)
                  dsp(14) = HDII(i,j,k)
                  dsp(15) = HeHII(i,j,k)
               endif
               if ( itmask_metal(i) .ne. MASK_FALSE ) then
                  if ( imchem .eq. 1 ) then
                     dsp(16) = CI(i,j,k)
                     dsp(17) = CII(i,j,k)
                     dsp(18) = CO(i,j,k)
                     dsp(19) = CO2(i,j,k)
                     dsp(20) = OI(i,j,k)
                     dsp(21) = OH(i,j,k)
                     dsp(22) = H2O(i,j,k)
                     dsp(23) = O2(i,j,k)
                     dsp(24) = SiI(i,j,k)
                     dsp(25) = SiOI(i,j,k)
                     dsp(26) = SiO2I(i,j,k)
                     dsp(27) = CH(i,j,k)
                     dsp(28) = CH2(i,j,k)
                     dsp(29) = COII(i,j,k)
                     dsp(30) = OII(i,j,k)
                     dsp(31) = OHII(i,j,k)
                     dsp(32) = H2OII(i,j,k)
                     dsp(33) = H3OII(i,j,k)
                     dsp(34) = O2II(i,j,k)
                     if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
                        if (idspecies .gt. 0) then
                           dsp(35) = Mg(i,j,k)
                        endif
                        if (idspecies .gt. 1) then
                           dsp(36) = Al(i,j,k)
                           dsp(37) = S(i,j,k)
                           dsp(38) = Fe(i,j,k)
                        endif
                     endif
                  endif
                  if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
                     if (idspecies .gt. 0) then
                        dsp(39) = MgSiO3(i,j,k)
                        dsp(40) = AC(i,j,k)
                     endif
                     if (idspecies .gt. 1) then
                        dsp(41) = SiM(i,j,k)
                        dsp(42) = FeM(i,j,k)
                        dsp(43) = Mg2SiO4(i,j,k)
                        dsp(44) = Fe3O4(i,j,k)
                        dsp(45) = SiO2D(i,j,k)
                        dsp(46) = MgO(i,j,k)
                        dsp(47) = FeS(i,j,k)
                        dsp(48) = Al2O3(i,j,k)
                     endif
                     if (idspecies .gt. 2) then
                        dsp(49) = reforg(i,j,k)
                        dsp(50) = volorg(i,j,k)
                        dsp(51) = H2Oice(i,j,k)
                     endif
                  endif
               endif ! if ( itmask_metal(i) )
               dsp(i_eng) = e(i,j,k)

               id = 0
               if (ispecies .gt. 0) then
                  do isp = 1, 6
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               if (ispecies .gt. 1) then
                  do isp = 7, 9
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               if (ispecies .gt. 2) then
                  do isp = 10, 12
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               if (ispecies .gt. 3) then
                  do isp = 13, 15
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               if (itmask_metal(i) .ne. MASK_FALSE) then
                  if (imchem .eq. 1) then
                     do isp = 16, 34
                        id = id + 1
                        idsp(id) = isp
                     enddo
                     if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
                        if (idspecies .gt. 0) then
                           do isp = 35, 35
                              id = id + 1
                              idsp(id) = isp
                           enddo
                        endif
                        if (idspecies .gt. 1) then
                           do isp = 36, 38
                              id = id + 1
                              idsp(id) = isp
                           enddo
                        endif
                     endif
                  endif
                  if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
                     if (idspecies .gt. 0) then
                        do isp = 39, 40
                           id = id + 1
                           idsp(id) = isp
                        enddo
                     endif
                     if (idspecies .gt. 1) then
                        do isp = 41, 48
                           id = id + 1
                           idsp(id) = isp
                        enddo
                     endif
                     if (idspecies .gt. 2) then
                        do isp = 49, 51
                           id = id + 1
                           idsp(id) = isp
                        enddo
                     endif
                  endif
               endif
               if ( imp_eng(i) .eq.1 ) then
                  id = id + 1
                  idsp(id) = i_eng
               endif

!              Save arrays at ttot(i)

               dsp0 = dsp
               ddsp = 0.d0

!              Search for the timestep for which chemistry converges

               ierror=1
               itr_time=0
               do while ((ierror.eq.1))

!                 If not converge, restore arrays at ttot(i)

                  dsp = dsp0
                  ddsp = 0.d0

!                 Iteration to solve ODEs

                  err_max=1.d2
                  itr=0
                  do while (err_max.gt.1.d-8)
                    if(itr.ge.20) then
                       ierror = 1
                       go to 9996
                    endif

            call lookup_cool_rates0d(dtit(i),
     &   d(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), 
     &   nsp, dsp, dspdot, nratec, 
     &   iexpand, ispecies, imetal, imcool,
     &   idust, idustall, idustfield, idustrec,
     &   ih2co, ipiht, iter, igammah,
     &   aye, temstart, temend, z_solar, fgr,
     &   utem, uxyz, uaye, urho, utim,
     &   gamma, fh,
     &   ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &   ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &   reHeII2a, reHeIIIa, brema, compa, gammaha,
     &   isrf, regra, gamma_isrfa, comp_xraya, comp_temp,
     &   piHI, piHeI, piHeII, comp1, comp2,
     &   metal(i,j,k), dust(i,j,k),
     &   hyd01ka, h2k01a, vibha, rotha, rotla,
     &   hyd01k(i), h2k01(i), vibh(i), roth(i), rotl(i),
     &   gpldla, gphdla, gpldl(i), gphdl(i),
     &   hdltea, hdlowa, hdlte(i), hdlow(i),
     &   gaHIa, gaH2a, gaHea, 
     &   gaHpa, gaela,
     &   h2ltea, gasgra,
     &   ceHI(i), ceHeI(i), ceHeII(i), ciHI(i), ciHeI(i), 
     &   ciHeIS(i), ciHeII(i),
     &   reHII(i), reHeII1(i), reHeII2(i), reHeIII(i), brem(i),
     &   indixe(i), t1(i), t2(i), logtem(i), tdef(i), edot(i),
     &   tgas(i), tgasold(i), mmw(i), p2d(i), tdust(i), metallicity(i),
     &   dust2gas(i), rhoH(i), mynh(i), myde(i),
     &   gammaha_eff(i), gasgr_tdust(i), regr(i),
     &   iradshield, avgsighi, avgsighei, avgsigheii,
     &   iradtrans, irt_honly,
     &   ih2optical, iciecool, ih2cr, ihdcr, ciecoa, cieco(i),
     &   icmbTfloor, iClHeat, clEleFra,
     &   priGridRank, priGridDim,
     &   priPar1, priPar2, priPar3, priPar4, priPar5,
     &   priDataSize, priCooling, priHeating, priMMW,
     &   metGridRank, metGridDim,
     &   metPar1, metPar2, metPar3, metPar4, metPar5,
     &   metDataSize, metCooling, metHeating, clnew,
     &   iVheat, iMheat, Vheat(i,j,k), Mheat(i,j,k),
     &   iTfloor, Tfloor_scalar, Tfloor(i,j,k),
     &   iisrffield, isrf_habing(i,j,k),
     &   ithreebody, anydust, iH2shield, 
     &   k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &   k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &   k17a, k18a, k19a, k22a,
     &   k24, k25, k26, k27, k28, k29, k30, k31,
     &   k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &   k57a, k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &   ncrna, ncrd1a, ncrd2a, 
     &   h2dust(i), ncrn(i), ncrd1(i), ncrd2(i), 
     &   dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
     &   kphHI(i,j,k), kphHeI(i,j,k), kphHeII(i,j,k), kdissH2I(i,j,k), 
     &   photogamma(i,j,k), xH2shield(i,j,k), chunit, itmask_nr(i),
     &   itmask_metal(i)
     & , imchem, igrgr, ipcont, tmcool
     & , k125a, k129a, k130a, k131a, k132a
     & , k133a, k134a, k135a, k136a, k137a
     & , k148a, k149a, k150a, k151a, k152a
     & , k153a
     & , kz15a, kz16a, kz17a, kz18a, kz19a
     & , kz20a, kz21a, kz22a, kz23a, kz24a
     & , kz25a, kz26a, kz27a, kz28a, kz29a
     & , kz30a, kz31a, kz32a, kz33a, kz34a
     & , kz35a, kz36a, kz37a, kz38a, kz39a
     & , kz40a, kz41a, kz42a, kz43a, kz44a
     & , kz45a, kz46a, kz47a, kz48a, kz49a
     & , kz50a, kz51a, kz52a, kz53a, kz54a
     & , cieY06a
     & , LH2_N, LH2_Size
     & , LH2_D, LH2_T, LH2_H
     & , LH2_dD, LH2_dT, LH2_dH, LH2_L
     & , LHD_N, LHD_Size
     & , LHD_D, LHD_T, LHD_H
     & , LHD_dD, LHD_dT, LHD_dH, LHD_L
     & , LCI_N, LCI_Size
     & , LCI_D, LCI_T, LCI_H
     & , LCI_dD, LCI_dT, LCI_dH, LCI_L
     & , LCII_N, LCII_Size
     & , LCII_D, LCII_T, LCII_H
     & , LCII_dD, LCII_dT, LCII_dH, LCII_L
     & , LOI_N, LOI_Size
     & , LOI_D, LOI_T, LOI_H
     & , LOI_dD, LOI_dT, LOI_dH, LOI_L
     & , LCO_N, LCO_Size
     & , LCO_D, LCO_T, LCO_H
     & , LCO_dD, LCO_dT, LCO_dH, LCO_L
     & , LOH_N, LOH_Size
     & , LOH_D, LOH_T, LOH_H
     & , LOH_dD, LOH_dT, LOH_dH, LOH_L
     & , LH2O_N, LH2O_Size
     & , LH2O_D, LH2O_T, LH2O_H
     & , LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L
     & , alphap_N, alphap_Size
     & , alphap_D, alphap_T, alphap_dD, alphap_dT
     & , alphap_Data
     & , immulti, imabund, idspecies, itdmulti, idsub
     & , metal_loc(i,j,k), metal_C13(i,j,k), metal_C20(i,j,k)
     & , metal_C25(i,j,k), metal_C30(i,j,k), metal_F13(i,j,k)
     & , metal_F15(i,j,k), metal_F50(i,j,k), metal_F80(i,j,k)
     & , metal_P170(i,j,k), metal_P200(i,j,k), metal_Y19(i,j,k)
     & , SN0_N
     & , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     & , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     & , SN0_fFeS, SN0_fAl2O3
     & , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     & , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     & , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     & , SN0_r0FeS, SN0_r0Al2O3
     & , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     & , gr_N, gr_Size, gr_dT, gr_Td
     & , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     & , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     & , SN0_kpFeS, SN0_kpAl2O3
     & , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     & , h2dustSa, h2dustCa, grogra
     & , tSiM(i), tFeM(i), tMg2SiO4(i), tMgSiO3(i), tFe3O4(i)
     & , tAC(i), tSiO2D(i), tMgO(i), tFeS(i), tAl2O3(i)
     & , treforg(i), tvolorg(i), tH2Oice(i)
     & , gasgr2a, gamma_isrf2a
     & , imp_eng(i)
     & , idissHDI, kdissHDI(i,j,k), iionZ, kphCI(i,j,k), kphOI(i,j,k)
     & , idissZ, kdissCO(i,j,k), kdissOH(i,j,k), kdissH2O(i,j,k)
     & , iuseH2shield, iH2shieldcustom, f_shield_custom(i,j,k)
     &   )

                    do jsp = 1, nsp
                      dspj = eps * dsp(idsp(jsp))
                      do isp = 1, nsp
                        if(isp .eq. jsp) then
                          dsp1(idsp(isp)) = dsp(idsp(isp)) + dspj
                        else
                          dsp1(idsp(isp)) = dsp(idsp(isp))
                        endif
                      enddo

            call lookup_cool_rates0d(dtit(i),
     &   d(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), 
     &   nsp, dsp1, dspdot1, nratec, 
     &   iexpand, ispecies, imetal, imcool,
     &   idust, idustall, idustfield, idustrec,
     &   ih2co, ipiht, iter, igammah,
     &   aye, temstart, temend, z_solar, fgr,
     &   utem, uxyz, uaye, urho, utim,
     &   gamma, fh,
     &   ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &   ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &   reHeII2a, reHeIIIa, brema, compa, gammaha,
     &   isrf, regra, gamma_isrfa, comp_xraya, comp_temp,
     &   piHI, piHeI, piHeII, comp1, comp2,
     &   metal(i,j,k), dust(i,j,k),
     &   hyd01ka, h2k01a, vibha, rotha, rotla,
     &   hyd01k(i), h2k01(i), vibh(i), roth(i), rotl(i),
     &   gpldla, gphdla, gpldl(i), gphdl(i),
     &   hdltea, hdlowa, hdlte(i), hdlow(i),
     &   gaHIa, gaH2a, gaHea, 
     &   gaHpa, gaela,
     &   h2ltea, gasgra,
     &   ceHI(i), ceHeI(i), ceHeII(i), ciHI(i), ciHeI(i), 
     &   ciHeIS(i), ciHeII(i),
     &   reHII(i), reHeII1(i), reHeII2(i), reHeIII(i), brem(i),
     &   indixe(i), t1(i), t2(i), logtem(i), tdef(i), edot(i),
     &   tgas(i), tgasold(i), mmw(i), p2d(i), tdust(i), metallicity(i),
     &   dust2gas(i), rhoH(i), mynh(i), myde(i),
     &   gammaha_eff(i), gasgr_tdust(i), regr(i),
     &   iradshield, avgsighi, avgsighei, avgsigheii,
     &   iradtrans, irt_honly,
     &   ih2optical, iciecool, ih2cr, ihdcr, ciecoa, cieco(i),
     &   icmbTfloor, iClHeat, clEleFra,
     &   priGridRank, priGridDim,
     &   priPar1, priPar2, priPar3, priPar4, priPar5,
     &   priDataSize, priCooling, priHeating, priMMW,
     &   metGridRank, metGridDim,
     &   metPar1, metPar2, metPar3, metPar4, metPar5,
     &   metDataSize, metCooling, metHeating, clnew,
     &   iVheat, iMheat, Vheat(i,j,k), Mheat(i,j,k),
     &   iTfloor, Tfloor_scalar, Tfloor(i,j,k),
     &   iisrffield, isrf_habing(i,j,k),
     &   ithreebody, anydust, iH2shield, 
     &   k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &   k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &   k17a, k18a, k19a, k22a,
     &   k24, k25, k26, k27, k28, k29, k30, k31,
     &   k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &   k57a, k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &   ncrna, ncrd1a, ncrd2a, 
     &   h2dust(i), ncrn(i), ncrd1(i), ncrd2(i), 
     &   dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
     &   kphHI(i,j,k), kphHeI(i,j,k), kphHeII(i,j,k), kdissH2I(i,j,k), 
     &   photogamma(i,j,k), xH2shield(i,j,k), chunit, itmask_nr(i),
     &   itmask_metal(i)
     & , imchem, igrgr, ipcont, tmcool
     & , k125a, k129a, k130a, k131a, k132a
     & , k133a, k134a, k135a, k136a, k137a
     & , k148a, k149a, k150a, k151a, k152a
     & , k153a
     & , kz15a, kz16a, kz17a, kz18a, kz19a
     & , kz20a, kz21a, kz22a, kz23a, kz24a
     & , kz25a, kz26a, kz27a, kz28a, kz29a
     & , kz30a, kz31a, kz32a, kz33a, kz34a
     & , kz35a, kz36a, kz37a, kz38a, kz39a
     & , kz40a, kz41a, kz42a, kz43a, kz44a
     & , kz45a, kz46a, kz47a, kz48a, kz49a
     & , kz50a, kz51a, kz52a, kz53a, kz54a
     & , cieY06a
     & , LH2_N, LH2_Size
     & , LH2_D, LH2_T, LH2_H
     & , LH2_dD, LH2_dT, LH2_dH, LH2_L
     & , LHD_N, LHD_Size
     & , LHD_D, LHD_T, LHD_H
     & , LHD_dD, LHD_dT, LHD_dH, LHD_L
     & , LCI_N, LCI_Size
     & , LCI_D, LCI_T, LCI_H
     & , LCI_dD, LCI_dT, LCI_dH, LCI_L
     & , LCII_N, LCII_Size
     & , LCII_D, LCII_T, LCII_H
     & , LCII_dD, LCII_dT, LCII_dH, LCII_L
     & , LOI_N, LOI_Size
     & , LOI_D, LOI_T, LOI_H
     & , LOI_dD, LOI_dT, LOI_dH, LOI_L
     & , LCO_N, LCO_Size
     & , LCO_D, LCO_T, LCO_H
     & , LCO_dD, LCO_dT, LCO_dH, LCO_L
     & , LOH_N, LOH_Size
     & , LOH_D, LOH_T, LOH_H
     & , LOH_dD, LOH_dT, LOH_dH, LOH_L
     & , LH2O_N, LH2O_Size
     & , LH2O_D, LH2O_T, LH2O_H
     & , LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L
     & , alphap_N, alphap_Size
     & , alphap_D, alphap_T, alphap_dD, alphap_dT
     & , alphap_Data
     & , immulti, imabund, idspecies, itdmulti, idsub
     & , metal_loc(i,j,k), metal_C13(i,j,k), metal_C20(i,j,k)
     & , metal_C25(i,j,k), metal_C30(i,j,k), metal_F13(i,j,k)
     & , metal_F15(i,j,k), metal_F50(i,j,k), metal_F80(i,j,k)
     & , metal_P170(i,j,k), metal_P200(i,j,k), metal_Y19(i,j,k)
     & , SN0_N
     & , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     & , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     & , SN0_fFeS, SN0_fAl2O3
     & , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     & , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     & , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     & , SN0_r0FeS, SN0_r0Al2O3
     & , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     & , gr_N, gr_Size, gr_dT, gr_Td
     & , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     & , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     & , SN0_kpFeS, SN0_kpAl2O3
     & , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     & , h2dustSa, h2dustCa, grogra
     & , tSiM(i), tFeM(i), tMg2SiO4(i), tMgSiO3(i), tFe3O4(i)
     & , tAC(i), tSiO2D(i), tMgO(i), tFeS(i), tAl2O3(i)
     & , treforg(i), tvolorg(i), tH2Oice(i)
     & , gasgr2a, gamma_isrf2a
     & , imp_eng(i)
     & , idissHDI, kdissHDI(i,j,k), iionZ, kphCI(i,j,k), kphOI(i,j,k)
     & , idissZ, kdissCO(i,j,k), kdissOH(i,j,k), kdissH2O(i,j,k)
     & , iuseH2shield, iH2shieldcustom, f_shield_custom(i,j,k)
     &   )

                      do isp = 1, nsp
                        if ( (dsp(idsp(isp)).eq.0.d0)
     &                 .and. (dspdot1(idsp(isp))
     &                   .eq.  dspdot(idsp(isp))) ) then
                          der(idsp(isp),idsp(jsp)) = 0.d0
                        else
                          der(idsp(isp),idsp(jsp)) = 
     &                       (dspdot1(idsp(isp))
     &                       - dspdot(idsp(isp))) / dspj
                        endif
                      enddo

                    enddo ! jsp

                    do isp = 1, nsp
                      do jsp = 1, nsp
                        if(isp .eq. jsp) then
                          mtrx(isp,jsp) = 1.d0 - dtit(i) 
     &                       * der(idsp(isp),idsp(jsp))
                        else
                          mtrx(isp,jsp) =      - dtit(i)
     &                       * der(idsp(isp),idsp(jsp))
                        endif
                      enddo
                    enddo

                    do isp = 1, nsp
                      vec(isp) = dspdot(idsp(isp)) * dtit(i)
     &                       - ddsp(idsp(isp))
                    enddo

!                   to get more accuracy
                    do isp = 1, nsp
                      vec(isp) = vec(isp)/d(i,j,k)
                    enddo

                    call gaussj_g(nsp, mtrx, vec, ierror)
                    if(ierror .eq. 1) then
                      go to 9998
                    endif

!                   multiply with density again
                    do isp = 1, nsp
                      vec(isp) = vec(isp)*d(i,j,k)
                    enddo

                    do isp = 1, nsp
                      ddsp(idsp(isp)) = ddsp(idsp(isp)) + vec(isp)
                      dsp(idsp(isp))  = dsp(idsp(isp))  + vec(isp)
                    enddo

                    if (imp_eng(i) .eq. 1) then
                    if( (ispecies .gt. 0) .and. (icool .eq. 1) ) then
                       do isp = 1, nsp
                          if ( (dsp(idsp(isp)) .ne. dsp(idsp(isp)))
     &                    .or. (dsp(idsp(isp)) .le. 0._DKIND) ) then
                             ierror = 1
                             go to 9997
                          endif
                       enddo
                    endif
                    endif

                    err_max = 0.d0
                    do isp = 1, nsp
                      if(dsp(idsp(isp)) .gt. tiny8) then
                        err = dabs(vec(isp) / dsp(idsp(isp)))
                      else
                        err = 0.d0
                      endif
                      if(err .gt. err_max) then
                        err_max = err
                      endif
                    enddo

                    itr=itr+1
                  enddo

9998              continue
9997              continue
9996              continue


!             Check if the fractions are valid after an iteration

                  if( (ispecies .gt. 0) .and. (icool .eq. 1) ) then
                     do isp = 1, nsp
                        if ( (dsp(idsp(isp)) .ne. dsp(idsp(isp)))
     &                  .or. (dsp(idsp(isp)) .le. 0._DKIND) ) then
                           ierror = 1
                        endif
                     enddo
                  endif
                  if(ierror .eq. 1) then
                     dtit(i) = 0.5d0*dtit(i)
                  endif

                  itr_time=itr_time+1
               enddo

               if ( ispecies .gt. 0 ) then
                  de(i,j,k)      = dsp( 1)
                  HI(i,j,k)      = dsp( 2)
                  HII(i,j,k)     = dsp( 3)
                  HeI(i,j,k)     = dsp( 4)
                  HeII(i,j,k)    = dsp( 5)
                  HeIII(i,j,k)   = dsp( 6)
               endif
               if ( ispecies .gt. 1 ) then
                  HM(i,j,k)      = dsp( 7)
                  H2I(i,j,k)     = dsp( 8)
                  H2II(i,j,k)    = dsp( 9)
               endif
               if ( ispecies .gt. 2 ) then
                  DI(i,j,k)      = dsp(10)
                  DII(i,j,k)     = dsp(11)
                  HDI(i,j,k)     = dsp(12)
               endif
               if ( ispecies .gt. 3 ) then
                  DM(i,j,k)      = dsp(13)
                  HDII(i,j,k)    = dsp(14)
                  HeHII(i,j,k)   = dsp(15)
               endif
               if ( itmask_metal(i) .ne. MASK_FALSE ) then
                  if ( imchem .eq. 1 ) then
                     CI(i,j,k)      = dsp(16)
                     CII(i,j,k)     = dsp(17)
                     CO(i,j,k)      = dsp(18)
                     CO2(i,j,k)     = dsp(19)
                     OI(i,j,k)      = dsp(20)
                     OH(i,j,k)      = dsp(21)
                     H2O(i,j,k)     = dsp(22)
                     O2(i,j,k)      = dsp(23)
                     SiI(i,j,k)     = dsp(24)
                     SiOI(i,j,k)    = dsp(25)
                     SiO2I(i,j,k)   = dsp(26)
                     CH(i,j,k)      = dsp(27)
                     CH2(i,j,k)     = dsp(28)
                     COII(i,j,k)    = dsp(29)
                     OII(i,j,k)     = dsp(30)
                     OHII(i,j,k)    = dsp(31)
                     H2OII(i,j,k)   = dsp(32)
                     H3OII(i,j,k)   = dsp(33)
                     O2II(i,j,k)    = dsp(34)
                     if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
                        if (idspecies .gt. 0) then
                           Mg(i,j,k)      = dsp(35)
                        endif
                        if (idspecies .gt. 1) then
                           Al(i,j,k)      = dsp(36)
                           S(i,j,k)       = dsp(37)
                           Fe(i,j,k)      = dsp(38)
                        endif
                     endif
                  endif
                  if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
                     if (idspecies .gt. 0) then
                        MgSiO3(i,j,k)  = dsp(39)
                        AC(i,j,k)      = dsp(40)
                     endif
                     if (idspecies .gt. 1) then
                        SiM(i,j,k)     = dsp(41)
                        FeM(i,j,k)     = dsp(42)
                        Mg2SiO4(i,j,k) = dsp(43)
                        Fe3O4(i,j,k)   = dsp(44)
                        SiO2D(i,j,k)   = dsp(45)
                        MgO(i,j,k)     = dsp(46)
                        FeS(i,j,k)     = dsp(47)
                        Al2O3(i,j,k)   = dsp(48)
                     endif
                     if (idspecies .gt. 2) then
                        reforg(i,j,k)  = dsp(49)
                        volorg(i,j,k)  = dsp(50)
                        H2Oice(i,j,k)  = dsp(51)
                     endif
                  endif
               endif

               e(i,j,k)     = dsp(i_eng)

               deallocate(idsp)
               deallocate(vec)
               deallocate(mtrx)

            endif
            enddo ! end Newton-Raphson scheme
      end
