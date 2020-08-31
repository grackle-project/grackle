#include "phys_const.def"
#define CORRECTED_BY_CHIAKI2020MAY

! calculate rates

      subroutine lookup_cool_rates0d(output, dtit,
     &               d, u, v, w, 
     &               nsp, dsp, dspdot, nratec, 
     &               iexpand, ispecies, imetal, imcool,
     &               idust, idustall, idustfield,
     &               ih2co, ipiht, iter, igammah,
     &               aye, temstart, temend, z_solar, fgr,
     &               utem, uxyz, uaye, urho, utim,
     &               gamma, fh,
     &               ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &               ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &               reHeII2a, reHeIIIa, brema, compa, gammaha,
     &               isrf, regra, gamma_isrfa, comp_xraya, comp_temp,
     &               piHI, piHeI, piHeII, comp1, comp2,
     &               metal, dust,
     &               hyd01ka, h2k01a, vibha, rotha, rotla,
     &               hyd01k, h2k01, vibh, roth, rotl,
     &               gpldla, gphdla, gpldl, gphdl,
     &               hdltea, hdlowa, hdlte, hdlow,
     &               gaHIa, gaH2a, gaHea, 
     &               gaHpa, gaela,
     &               h2ltea, gasgra,
     &               ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &               reHII, reHeII1, reHeII2, reHeIII, brem,
     &               indixe, t1, t2, logtem, tdef, edot,
     &               tgas, tgasold, mmw, p2d, tdust, metallicity,
     &               dust2gas, rhoH, mynh, myde,
     &               gammaha_eff, gasgr_tdust, regr,
     &               iradshield, avgsighi, avgsighei, avgsigheii,
     &               iradtrans, irt_honly,
     &               ih2optical, iciecool, ciecoa, cieco,
     &               icmbTfloor, iClHeat, clEleFra,
     &               priGridRank, priGridDim,
     &               priPar1, priPar2, priPar3, priPar4, priPar5,
     &               priDataSize, priCooling, priHeating, priMMW,
     &               metGridRank, metGridDim,
     &               metPar1, metPar2, metPar3, metPar4, metPar5,
     &               metDataSize, metCooling, metHeating, clnew,
     &               iVheat, iMheat, Vheat, Mheat,
     &               iisrffield, isrf_habing,
     &               ithreebody, anydust, iH2shield, 
     &               k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &               k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &               k17a, k18a, k19a, k22a,
     &               k24, k25, k26, k27, k28, k29, k30, k31,
     &               k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &               k57a, k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &               ncrna, ncrd1a, ncrd2a, 
     &               h2dust, ncrn, ncrd1, ncrd2, 
     &               dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
     &               kphHI, kphHeI, kphHeII, kdissH2I, 
     &               photogamma, xH2shield, chunit, itmask
     &             , imchem, igrgr
     &             , k125a, k129a, k130a, k131a, k132a
     &             , k133a, k134a, k135a, k136a, k137a
     &             , k148a, k149a, k150a, k151a, k152a
     &             , k153a
     &             , kz15a, kz16a, kz17a, kz18a, kz19a
     &             , kz20a, kz21a, kz22a, kz23a, kz24a
     &             , kz25a, kz26a, kz27a, kz28a, kz29a
     &             , kz30a, kz31a, kz32a, kz33a, kz34a
     &             , kz35a, kz36a, kz37a, kz38a, kz39a
     &             , kz40a, kz41a, kz42a, kz43a, kz44a
     &             , kz45a, kz46a, kz47a, kz48a, kz49a
     &             , kz50a, kz51a, kz52a, kz53a, kz54a
     &             , cieY06a
     &             , LH2_N, LH2_Size
     &             , LH2_D, LH2_T, LH2_H
     &             , LH2_dD, LH2_dT, LH2_dH, LH2_L
     &             , LHD_N, LHD_Size
     &             , LHD_D, LHD_T, LHD_H
     &             , LHD_dD, LHD_dT, LHD_dH, LHD_L
     &             , LCI_N, LCI_Size
     &             , LCI_D, LCI_T, LCI_H
     &             , LCI_dD, LCI_dT, LCI_dH, LCI_L
     &             , LCII_N, LCII_Size
     &             , LCII_D, LCII_T, LCII_H
     &             , LCII_dD, LCII_dT, LCII_dH, LCII_L
     &             , LOI_N, LOI_Size
     &             , LOI_D, LOI_T, LOI_H
     &             , LOI_dD, LOI_dT, LOI_dH, LOI_L
     &             , LCO_N, LCO_Size
     &             , LCO_D, LCO_T, LCO_H
     &             , LCO_dD, LCO_dT, LCO_dH, LCO_L
     &             , LOH_N, LOH_Size
     &             , LOH_D, LOH_T, LOH_H
     &             , LOH_dD, LOH_dT, LOH_dH, LOH_L
     &             , LH2O_N, LH2O_Size
     &             , LH2O_D, LH2O_T, LH2O_H
     &             , LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L
     &             , alphap_N, alphap_Size
     &             , alphap_D, alphap_T, alphap_dD, alphap_dT
     &             , alphap_Data
     &             , grain_N, grain_Size
     &             , grain_D, grain_T, grain_dD, grain_dT
     &             , Hgrain, Tgrain, Ograin, Lgrain
     &             , immulti, imabund, idspecies, itdmulti, idsub
     &             , metal_loc
     &             , metal_C13, metal_C20, metal_C25, metal_C30
     &             , metal_F13, metal_F15, metal_F50, metal_F80
     &             , metal_P170, metal_P200, metal_Y19
     &             , SN0_N
     &             , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &             , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &             , SN0_fFeS, SN0_fAl2O3
     &             , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &             , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &             , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &             , SN0_r0FeS, SN0_r0Al2O3
     &             , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &             , gr_N, gr_Size, gr_dT, gr_Td
     &             , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &             , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &             , SN0_kpFeS, SN0_kpAl2O3
     &             , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &             , h2dustSa, h2dustCa, grogra
     &             , tSiM, tFeM, tMg2SiO4, tMgSiO3, tFe3O4
     &             , tAC, tSiO2D, tMgO, tFeS, tAl2O3
     &             , treforg, tvolorg, tH2Oice
     &             , gasgr2a, gamma_isrf2a
     &             , imp_eng)

c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      integer output, nratec, nsp,
     &        iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        imcool, idust, idustall, idustfield,
     &        igammah, ih2optical, iciecool, ithreebody,
     &        ndratec, clnew, iVheat, iMheat, iH2shield, iradshield,
     &        iradtrans, irt_honly, iisrffield
     &       ,imchem, igrgr
      logical itmask, anydust

      real*8  dx_cgs, aye, temstart, temend, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, z_solar, 
     &        fgr, dtemstart, dtemend, clEleFra

!  Density, energy and velocity fields fields
      integer, parameter :: i_eng = 52
      real*8  dsp(i_eng), dspdot(i_eng)
      R_PREC  d, e, u, v, w,
     &        metal, dust,
     &        Vheat, Mheat
       R_PREC metal_loc
     &      , metal_C13          , metal_C20          
     &      , metal_C25          , metal_C30          
     &      , metal_F13          , metal_F15          
     &      , metal_F50          , metal_F80          
     &      , metal_P170          , metal_P200          
     &      , metal_Y19          

!  Radiative transfer fields

      R_PREC  kphHI, kphHeI, kphHeII,
     &        kdissH2I, photogamma

!  H2 self-shielding length-scale field

      R_PREC  xH2shield

!  Interstellar radiation field for dust heating

      R_PREC  isrf_habing

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

!  Parameters

      integer itmax
      parameter (itmax = 10000)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.0e-05_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.0e-10_RKIND)
#endif

      real*8 mh, pi
      parameter (mh = mass_h, pi = pi_val)

!  Locals

      integer i, j, k, iter
      integer t, dj, dk
      real*8 ttmin, dom, energy, comp1, comp2
      real*8 coolunit, dbase1, tbase1, xbase1, chunit, uvel
      real*8 heq1, heq2, eqk221, eqk222, eqk131, eqk132,
     &       eqt1, eqt2, eqtdef, dheq, heq, dlogtem,
     &       c_ljeans

!  row temporaries

      integer*8 indixe
      real*8 t1, t2, logtem, tdef,
     &       dtit, ttot, p2d, tgas, tgasold,
     &       tdust, metallicity, dust2gas,
     &       rhoH, mmw, mynh, myde, gammaha_eff,
     &       gasgr_tdust, regr, olddtit

!  Rate equation row temporaries

      real*8 HI , HII , HeI , HeII , HeIII ,
     &       HM , H2I , H2II ,
     &       de , dedot,HIdot, dedot_prev,
     &       DI , DII , HDI , HIdot_prev,
     &       k24shield, k25shield, k26shield,
     &       k28shield, k29shield, k30shield,
     &       k31shield,
     &       k1 , k2 , k3 , k4 , k5 ,
     &       k6 , k7 , k8 , k9 , k10,
     &       k11, k12, k13, k14, k15,
     &       k16, k17, k18, k19, k22,
     &       k50, k51, k52, k53, k54,
     &       k55, k56, k57, k58,
     &       k13dd(14), h2dust,
     &       ncrn, ncrd1, ncrd2
      R_PREC DM , HDII , HeHII 
     &     , CI   , CII    , CO
     &     , CO2  , OI     , OH
     &     , H2O  , O2     , SiI
     &     , SiOI , SiO2I  , CH
     &     , CH2  , COII   , OII
     &     , OHII , H2OII  , H3OII
     &     , O2II , Mg     , Al
     &     , S    , Fe     
      R_PREC  SiM, FeM, Mg2SiO4
     &      , MgSiO3, Fe3O4, AC
     &      , SiO2D, MgO, FeS
     &      , Al2O3
     &      , reforg, volorg, H2Oice
      real*8 k125, k129, k130, k131, k132
     &     , k133, k134, k135, k136, k137
     &     , k148, k149, k150, k151, k152
     &     , k153
     &     , kz15, kz16, kz17, kz18, kz19
     &     , kz20, kz21, kz22, kz23, kz24
     &     , kz25, kz26, kz27, kz28, kz29
     &     , kz30, kz31, kz32, kz33, kz34
     &     , kz35, kz36, kz37, kz38, kz39
     &     , kz40, kz41, kz42, kz43, kz44
     &     , kz45, kz46, kz47, kz48, kz49
     &     , kz50, kz51, kz52, kz53, kz54

!  Cooling/heating row locals

      real*8 ceHI, ceHeI, ceHeII,
     &       ciHI, ciHeI, ciHeIS, ciHeII,
     &       reHII, reHeII1, reHeII2, reHeIII,
     &       brem, edot
      real*8 hyd01k, h2k01, vibh, roth, rotl,
     &       gpldl, gphdl, hdlte, hdlow, cieco
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
       integer grain_N(2), grain_Size
       real*8  grain_D(grain_N(1)), grain_T(grain_N(2))
     &       , grain_dD, grain_dT
     &       , Hgrain(grain_Size), Tgrain(grain_Size)
     &       , Ograin(grain_Size), Lgrain(grain_Size)
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
      real*8  h2dustSa(nratec, ndratec), h2dustCa(nratec, ndratec)
      real*8  grogra(nratec)
      real*8  kdSiM    , kdFeM    , kdMg2SiO4    
     &      , kdMgSiO3    , kdFe3O4    , kdAC    
     &      , kdSiO2D    , kdMgO    , kdFeS    
     &      , kdAl2O3
     &      , kdreforg    , kdvolorg    , kdH2Oice    
!     grain temperature
      real*8  tSiM    , tFeM    , tMg2SiO4    
     &      , tMgSiO3    , tFe3O4    , tAC    
     &      , tSiO2D    , tMgO    , tFeS    
     &      , tAl2O3    
     &      , treforg    , tvolorg    , tH2Oice    
      real*8  gasgr2a(nratec), gamma_isrf2a

!     locals

      integer imp_eng, isp
      real*8 scoef, acoef
      real*8 atten, H2delta, h2heatfac
      logical itmask_metal

!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Iteration mask for metal-rich cells

      itmask_metal = .false.
      if (imetal .eq. 1) then
         itmask_metal = itmask
         if (metal .lt. 1.e-9_DKIND * d) then
            itmask_metal = .false.
         endif
      endif

      if ( ispecies .gt. 0 ) then
         de      = dsp( 1)
         HI      = dsp( 2)
         HII     = dsp( 3)
         HeI     = dsp( 4)
         HeII    = dsp( 5)
         HeIII   = dsp( 6)
      endif
      if ( ispecies .gt. 1 ) then
         HM      = dsp( 7)
         H2I     = dsp( 8)
         H2II    = dsp( 9)
      endif
      if ( ispecies .gt. 2 ) then
         DI      = dsp(10)
         DII     = dsp(11)
         HDI     = dsp(12)
      endif
      if ( ispecies .gt. 3 ) then
         DM      = dsp(13)
         HDII    = dsp(14)
         HeHII   = dsp(15)
      endif
      if ( itmask_metal ) then
      if (imchem .eq. 1) then
         CI      = dsp(16)
         CII     = dsp(17)
         CO      = dsp(18)
         CO2     = dsp(19)
         OI      = dsp(20)
         OH      = dsp(21)
         H2O     = dsp(22)
         O2      = dsp(23)
         SiI     = dsp(24)
         SiOI    = dsp(25)
         SiO2I   = dsp(26)
         CH      = dsp(27)
         CH2     = dsp(28)
         COII    = dsp(29)
         OII     = dsp(30)
         OHII    = dsp(31)
         H2OII   = dsp(32)
         H3OII   = dsp(33)
         O2II    = dsp(34)
         if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
         if (idspecies .gt. 0) then
         Mg      = dsp(35)
         endif
         if (idspecies .gt. 1) then
         Al      = dsp(36)
         S       = dsp(37)
         Fe      = dsp(38)
         endif
         endif
      endif
      if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
      if (idspecies .gt. 0) then
         MgSiO3  = dsp(39)
         AC      = dsp(40)
      endif
      if (idspecies .gt. 1) then
         SiM     = dsp(41)
         FeM     = dsp(42)
         Mg2SiO4 = dsp(43)
         Fe3O4   = dsp(44)
         SiO2D   = dsp(45)
         MgO     = dsp(46)
         FeS     = dsp(47)
         Al2O3   = dsp(48)
      endif
      if (idspecies .gt. 2) then
         reforg  = dsp(49)
         volorg  = dsp(50)
         H2Oice  = dsp(51)
      endif
      endif
      endif
      e     = dsp(i_eng)

!     Compute the cooling rate, tgas, tdust, and metallicity for this row

      if (imp_eng .eq. 1) then
            call cool1d_multi_g(
     &                 d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                 1,  1,  1, nratec, 
     &                iexpand, ispecies, imetal, imcool,
     &                idust, idustall, idustfield,
     &                   1,  0,  0, 1, 1, ih2co, ipiht, iter, igammah,
     &                aye, temstart, temend, z_solar, fgr,
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                isrf, regra, gamma_isrfa, comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal, dust,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                h2ltea, gasgra,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, mmw, p2d, tdust, metallicity,
     &                dust2gas, rhoH, mynh, myde,
     &                gammaha_eff, gasgr_tdust, regr,
     &                iradshield, avgsighi, avgsighei, avgsigheii,
     &                k24, k26,
     &                iradtrans, photogamma,
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat,
     &                iisrffield, isrf_habing, itmask
     &              , imchem, igrgr
     &              , DM, HDII, HeHII
     &              , CI, CII, CO, CO2
     &              , OI, OH, H2O, O2
     &              , SiI, SiOI, SiO2I
     &              , CH, CH2, COII, OII
     &              , OHII, H2OII, H3OII, O2II
     &              , Mg, Al, S, Fe
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , reforg, volorg, H2Oice
     &              , cieY06a
     &              , LH2_N, LH2_Size
     &              , LH2_D, LH2_T, LH2_H
     &              , LH2_dD, LH2_dT, LH2_dH, LH2_L
     &              , LHD_N, LHD_Size
     &              , LHD_D, LHD_T, LHD_H
     &              , LHD_dD, LHD_dT, LHD_dH, LHD_L
     &              , LCI_N, LCI_Size
     &              , LCI_D, LCI_T, LCI_H
     &              , LCI_dD, LCI_dT, LCI_dH, LCI_L
     &              , LCII_N, LCII_Size
     &              , LCII_D, LCII_T, LCII_H
     &              , LCII_dD, LCII_dT, LCII_dH, LCII_L
     &              , LOI_N, LOI_Size
     &              , LOI_D, LOI_T, LOI_H
     &              , LOI_dD, LOI_dT, LOI_dH, LOI_L
     &              , LCO_N, LCO_Size
     &              , LCO_D, LCO_T, LCO_H
     &              , LCO_dD, LCO_dT, LCO_dH, LCO_L
     &              , LOH_N, LOH_Size
     &              , LOH_D, LOH_T, LOH_H
     &              , LOH_dD, LOH_dT, LOH_dH, LOH_L
     &              , LH2O_N, LH2O_Size
     &              , LH2O_D, LH2O_T, LH2O_H
     &              , LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L
     &              , alphap_N, alphap_Size
     &              , alphap_D, alphap_T, alphap_dD, alphap_dT
     &              , alphap_Data
     &              , grain_N, grain_Size
     &              , grain_D, grain_T, grain_dD, grain_dT
     &              , Tgrain, Ograin, Lgrain
     &              , immulti, imabund, idspecies, itdmulti, idsub
     &              , metal_loc
     &              , metal_C13, metal_C20, metal_C25, metal_C30
     &              , metal_F13, metal_F15, metal_F50, metal_F80
     &              , metal_P170, metal_P200, metal_Y19
     &              , SN0_N
     &              , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &              , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &              , SN0_fFeS, SN0_fAl2O3
     &              , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &              , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &              , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &              , SN0_r0FeS, SN0_r0Al2O3
     &              , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &              , gr_N, gr_Size, gr_dT, gr_Td
     &              , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &              , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &              , SN0_kpFeS, SN0_kpAl2O3
     &              , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &              , tSiM, tFeM, tMg2SiO4, tMgSiO3, tFe3O4
     &              , tAC, tSiO2D, tMgO, tFeS, tAl2O3
     &              , treforg, tvolorg, tH2Oice
     &              , gasgr2a, gamma_isrf2a
     &  )
      endif

! -----------------------------------------------------------
! This routine uses the temperature to look up the chemical
!   rates which are tabulated in a log table as a function
!   of temperature.

            call lookup_cool_rates1d_g(temstart, temend, nratec, 1, 1,
     &                 0,  0, ithreebody,
     &                 1,  1,  1, ispecies, anydust,
     &                iH2shield, iradshield,
     &                tgas, mmw, d, HI, HII, HeI, HeII, HeIII, 
     &                HM, H2I, H2II, DI, DII, HDI,
     &                tdust, dust2gas,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &                k57a, k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &                ncrna, ncrd1a, ncrd2a, 
     &                avgsighi, avgsighei, avgsigheii, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26, k28, k29, k30, k31,
     &                k50, k51, k52, k53, k54, k55, k56, k57,
     &                k58, k13dd, k24shield, k25shield, k26shield,
     &                k28shield, k29shield, k30shield,
     &                k31shield, h2dust, ncrn, ncrd1, ncrd2, 
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
     &                iradtrans, kdissH2I, xH2shield, itmask
     &              , fh, metal
     &              , DM, HDII, HeHII, imetal, imchem, igrgr
     &              , CI, CII, CO, CO2
     &              , OI, OH, H2O, O2
     &              , SiI, SiOI, SiO2I
     &              , CH, CH2, COII, OII
     &              , OHII, H2OII, H3OII, O2II
     &              , Mg, Al, S, Fe
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , reforg, volorg, H2Oice
     &              , k125a, k129a, k130a, k131a, k132a
     &              , k133a, k134a, k135a, k136a, k137a
     &              , k148a, k149a, k150a, k151a, k152a
     &              , k153a
     &              , kz15a, kz16a, kz17a, kz18a, kz19a
     &              , kz20a, kz21a, kz22a, kz23a, kz24a
     &              , kz25a, kz26a, kz27a, kz28a, kz29a
     &              , kz30a, kz31a, kz32a, kz33a, kz34a
     &              , kz35a, kz36a, kz37a, kz38a, kz39a
     &              , kz40a, kz41a, kz42a, kz43a, kz44a
     &              , kz45a, kz46a, kz47a, kz48a, kz49a
     &              , kz50a, kz51a, kz52a, kz53a, kz54a
     &              , k125 , k129 , k130 , k131 , k132 
     &              , k133 , k134 , k135 , k136 , k137 
     &              , k148 , k149 , k150 , k151 , k152 
     &              , k153 
     &              , kz15 , kz16 , kz17 , kz18 , kz19 
     &              , kz20 , kz21 , kz22 , kz23 , kz24 
     &              , kz25 , kz26 , kz27 , kz28 , kz29 
     &              , kz30 , kz31 , kz32 , kz33 , kz34 
     &              , kz35 , kz36 , kz37 , kz38 , kz39 
     &              , kz40 , kz41 , kz42 , kz43 , kz44 
     &              , kz45 , kz46 , kz47 , kz48 , kz49 
     &              , kz50 , kz51 , kz52 , kz53 , kz54 
     &              , grain_N, grain_Size
     &              , grain_D, grain_T, grain_dD, grain_dT
     &              , Hgrain
     &              , immulti, imabund, idspecies, itdmulti, idsub
     &              , metal_loc
     &              , metal_C13, metal_C20, metal_C25, metal_C30
     &              , metal_F13, metal_F15, metal_F50, metal_F80
     &              , metal_P170, metal_P200, metal_Y19
     &              , SN0_N
     &              , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &              , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &              , SN0_fFeS, SN0_fAl2O3
     &              , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &              , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &              , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &              , SN0_r0FeS, SN0_r0Al2O3
     &              , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &              , gr_N, gr_Size, gr_dT, gr_Td
     &              , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &              , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &              , SN0_kpFeS, SN0_kpAl2O3
     &              , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &              , h2dustSa, h2dustCa, rhoH, grogra, dtit
     &              , kdSiM, kdFeM, kdMg2SiO4
     &              , kdMgSiO3, kdFe3O4, kdAC, kdSiO2D, kdMgO, kdFeS
     &              , kdAl2O3, kdreforg, kdvolorg, kdH2Oice
     &              , tSiM, tFeM, tMg2SiO4, tMgSiO3, tFe3O4
     &              , tAC, tSiO2D, tMgO, tFeS, tAl2O3
     &              , treforg, tvolorg, tH2Oice
     &      )

!           Compute dedot and HIdot, the rates of change of de and HI
!             (should add itmask to this call)

      if (imp_eng .eq. 1) then
            call rate_timestep_g(
     &                     dedot, HIdot, ispecies, anydust,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                      1,  1,  1,  0,  0, 1, 1,
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
     &                         )
      endif

!     Initialize

      dspdot = 0.d0


!   Heating/cooling rate (per unit volume -> gas mass)

      dspdot(i_eng) = edot / d

!   A) the 6-species integrator
!      
      if (ispecies .eq. 1) then




!        1) HI

            scoef  = k2   *HII       *de       
            acoef  = k1   *de       
     &             + k57   *HI       
     &             + k58   *HeI       /4._DKIND
     &             + k24shield   
            if (iradtrans .eq. 1) acoef = acoef + kphHI       
          dspdot(2) = dspdot(2) + (scoef - acoef * HI)














!        2) HII
c 
            scoef  = k1   *HI    *de       
     &             + k57   *HI    *HI    
     &             + k58   *HI    *HeI       /4._DKIND
     &             + k24shield   *HI    
            if (iradtrans .eq. 1)
     &          scoef = scoef + kphHI       *HI    
            acoef  = k2   *de        
          dspdot(3) = dspdot(3) + (scoef - acoef * HII)
















!        3) Electron density

            scoef = 0._DKIND
     &                 + k57   *HI    *HI    
     &                 + k58   *HI    *HeI       /4._DKIND
     &                 + k24shield   *HI       
     &                 + k25shield   *HeII       /4._DKIND
     &                 + k26shield   *HeI       /4._DKIND

            if ( (iradtrans .eq. 1) .and. ( irt_honly .eq. 0) )
     &          scoef = scoef + kphHI        * HI       
     &                + kphHeI         * HeI         / 4._DKIND
     &                + kphHeII        * HeII        / 4._DKIND
            if ( (iradtrans .eq. 1) .and. ( irt_honly .eq. 1) )
     &          scoef = scoef + kphHI        * HI       



            acoef = -(k1   *HI             - k2   *HII       
     &              + k3   *HeI       /4._DKIND -
     &           k6   *HeIII       /4._DKIND
     &              + k5   *HeII       /4._DKIND -
     &           k4   *HeII       /4._DKIND)
          dspdot(1) = dspdot(1) + (scoef - acoef * de)





      endif ! if (ispecies .eq. 1) then

!  --- (B) Do helium chemistry in any case: (for all ispecies values) ---




!        4) HeI

         scoef  = k4   *HeII       *de       
         acoef  = k3   *de       
     &                + k26shield   

         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       acoef = acoef + kphHeI       
         if (ispecies .gt. 3) then
            scoef = scoef +  4._DKIND * ( 0._DKIND
     &          + k152    * HeHII        *    HI        /  5._DKIND
     &          + k153    * HeHII        *    de        /  5._DKIND
     &         )
            acoef = acoef
     &          + k148    *   HII       
     &          + k149    *   HII       
     &          + k150    *  H2II        /  2._DKIND
         endif
          dspdot(4) = dspdot(4) + (scoef - acoef * HeI)


!        5) HeII

         scoef  = k3   *HeI    *de       
     &          + k6   *HeIII       *de       
     &          + k26shield   *HeI    

         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       scoef = scoef + kphHeI       *HeI    

         acoef  = k4   *de        + k5   *de
     &          + k25shield   

         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       acoef = acoef + kphHeII       
         if (ispecies .gt. 3) then
            acoef = acoef
     &          + k151    *    HI       
         endif
          dspdot(5) = dspdot(5) + (scoef - acoef * HeII)


!       6) HeIII

         scoef   = k5   *HeII    *de       
     &           + k25shield   *HeII    
         if ((iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       scoef = scoef + kphHeII        * HeII    
         acoef   = k6   *de       
          dspdot(6) = dspdot(6) + (scoef - acoef * HeIII)





c --- (C) Now do extra 3-species for molecular hydrogen ---

      if (ispecies .gt. 1) then

!        First, do HI/HII with molecular hydrogen terms




!        1) HI
!     
            scoef  =      k2    * HII        * de       
     &             + 2._DKIND*k13   * HI         * H2I       /2._DKIND
     &             +      k11   * HII        * H2I       /2._DKIND
     &             + 2._DKIND*k12   * de         * H2I       /2._DKIND
     &             +      k14   * HM         * de       
     &             +      k15   * HM         * HI       
     &             + 2._DKIND*k16   * HM         * HII       
     &             + 2._DKIND*k18   * H2II       * de       /2._DKIND
     &             +      k19   * H2II       * HM       /2._DKIND
     &             + 2._DKIND*k31shield      * H2I       /2._DKIND

            acoef  =      k1    * de       
     &             +      k7    * de       
     &             +      k8    * HM       
     &             +      k9    * HII       
     &             +      k10   * H2II       /2._DKIND
     &             + 2._DKIND*k22   * HI       **2
     &             +      k57   * HI       
     &             +      k58   * HeI       /4._DKIND
     &             + k24shield   

            if (iradtrans .eq. 1) acoef = acoef + kphHI       

            if (anydust) then
            if(itmask_metal   ) then
               acoef = acoef + 2._DKIND * h2dust    * rhoH   
            endif
            endif
!           contribution of minor species
            if (ispecies .gt. 2) then
               scoef = scoef 
     &               + k50    * HII        * DI         / 2._DKIND
     &               + k54    * H2I        * DI         / 4._DKIND
               acoef = acoef
     &               + k51    * DII        / 2._DKIND
     &               + k55    * HDI        / 3._DKIND
            endif

            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k131    *  HDII        *    de        /  3._DKIND
     &             + k134    *   HII        *    DM        /  2._DKIND
     &             + k135    *    HM        *    DI        /  2._DKIND
     &             + k150    *   HeI        *  H2II        /  8._DKIND
     &             + k153    * HeHII        *    de        /  5._DKIND
               acoef = acoef
     &             + k125    *  HDII        /  3._DKIND
     &             + k130    *   DII        /  2._DKIND
     &             + k136    *    DM        /  2._DKIND
     &             + k137    *    DM        /  2._DKIND
     &             + k151    *  HeII        /  4._DKIND
     &             + k152    * HeHII        /  5._DKIND
            endif

            if ( (imchem .eq. 1) .and. (itmask_metal   ) ) then
               scoef = scoef
     &             + kz20    *    CI        *   H2I        / 24._DKIND
     &             + kz21    *    OI        *   H2I        / 32._DKIND
     &             + kz22    *   HII        *    OI        / 16._DKIND
     &             + kz23    *   H2I        *    CH        / 26._DKIND
     &             + kz24    *   H2I        *    OH        / 34._DKIND
     &             + kz26    *    OH        *    CO        / 476._DKIND
     &             + kz28    *    CI        *    OH        / 204._DKIND
     &             + kz32    *    OI        *    CH        / 208._DKIND
     &             + kz33    *    OI        *    OH        / 272._DKIND
     &             + kz34    *   HII        *    OH        / 17._DKIND
     &             + kz35    *   HII        *   H2O        / 18._DKIND
     &             + kz36    *   HII        *    O2        / 32._DKIND
     &             + kz37    *   CII        *    OH        / 204._DKIND
     &             + kz40    *   OII        *   H2I        / 32._DKIND
     &             + kz41    *  OHII        *   H2I        / 34._DKIND
     &             + kz42    * H2OII        *   H2I        / 36._DKIND
     &             + kz46    * H2OII        *    de        / 18._DKIND
     &             + kz48    * H3OII        *    de        / 19._DKIND
     &             + kz49    * H3OII        *    de        / 9.5_DKIND
     &             + kz52    *   SiI        *    OH        / 476._DKIND
     &             + kz54    *  SiOI        *    OH        / 748._DKIND
               acoef = acoef
     &             + kz15    *    CH        / 13._DKIND
     &             + kz16    *   CH2        / 14._DKIND
     &             + kz17    *    OH        / 17._DKIND
     &             + kz18    *   H2O        / 18._DKIND
     &             + kz19    *    O2        / 32._DKIND
     &             + kz27    *    CI        / 12._DKIND
     &             + kz30    *    OI        / 16._DKIND
     &             + kz39    *   OII        / 16._DKIND
     &             + kz43    *  COII        / 28._DKIND
            endif
          dspdot(2) = dspdot(2) + (scoef - acoef * HI)













!          2) HII

            scoef  =    k1     * HI        * de       
     &             +    k10    * H2II       *HI       /2._DKIND
     &             +    k57    * HI        * HI       
     &             +    k58    * HI        * HeI       /4._DKIND
     &             + k24shield   *HI       

            if (iradtrans .eq. 1)
     &          scoef = scoef + kphHI        * HI       

            acoef  =    k2     * de       
     &             +    k9     * HI       
     &             +    k11    * H2I       /2._DKIND
     &             +    k16    * HM       
     &             +    k17    * HM       
!           contribution of minor species
            if (ispecies .gt. 2) then
               scoef = scoef
     &               + k51    * HI         * DII        / 2._DKIND
     &               + k52    * H2I        * DII        / 4._DKIND
               acoef = acoef
     &               + k50    * DI         / 2._DKIND
     &               + k53    * HDI        / 3._DKIND
            endif

            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k125    *  HDII        *    HI        /  3._DKIND
               acoef = acoef
     &             + k129    *    DI        /  2._DKIND
     &             + k134    *    DM        /  2._DKIND
     &             + k148    *   HeI        /  4._DKIND
     &             + k149    *   HeI        /  4._DKIND
            endif

            if ( (imchem .eq. 1) .and. (itmask_metal   ) ) then
               scoef = scoef
     &             + kz39    *   OII        *    HI        / 16._DKIND
     &             + kz43    *  COII        *    HI        / 28._DKIND
               acoef = acoef
     &             + kz22    *    OI        / 16._DKIND
     &             + kz34    *    OH        / 17._DKIND
     &             + kz35    *   H2O        / 18._DKIND
     &             + kz36    *    O2        / 32._DKIND
            endif
          dspdot(3) = dspdot(3) + (scoef - acoef * HII)

!     
!          3) electrons:

            scoef =   k8    * HM        * HI       
     &             +  k15   * HM        * HI       
     &             +  k17   * HM        * HII       
     &             +  k57   * HI        * HI       
     &             +  k58   * HI        * HeI       /4._DKIND
!
     &             + k24shield   *HI    
     &             + k25shield   *HeII    /4._DKIND
     &             + k26shield   *HeI    /4._DKIND

            if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0) )
     &          scoef = scoef + kphHI        * HI    
     &                + kphHeI         * HeI      / 4._DKIND
     &                + kphHeII        * HeII     / 4._DKIND
            if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 1) )
     &          scoef = scoef + kphHI        * HI    

            acoef = - (k1    *HI           - k2   *HII       
     &              +  k3    *HeI       /4._DKIND -
     &           k6   *HeIII       /4._DKIND
     &              +  k5    *HeII       /4._DKIND -
     &           k4   *HeII       /4._DKIND
     &              +  k14   *HM       
     &              -  k7    *HI       
     &              -  k18   *H2II       /2._DKIND)
!           contribution of minor species
            if (ispecies .gt. 2) then
               scoef = scoef
     &               + k56    * DI         * HM        / 2._DKIND
               acoef = acoef
     &               - k1     * DI         / 2._DKIND
     &               + k2     * DII        / 2._DKIND
            endif

            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k137    *    DM        *    HI        /  2._DKIND
               acoef = acoef
     &             + k131    *  HDII        /  3._DKIND
     &             + k132    *    DI        /  2._DKIND
     &             + k153    * HeHII        /  5._DKIND
            endif

            if ( (imchem .eq. 1) .and. (itmask_metal   ) ) then
               scoef = scoef
               acoef = acoef
     &             + kz44    *   CII        / 12._DKIND
     &             + kz45    *   OII        / 16._DKIND
     &             + kz46    * H2OII        / 18._DKIND
     &             + kz47    * H2OII        / 18._DKIND
     &             + kz48    * H3OII        / 19._DKIND
     &             + kz49    * H3OII        / 19._DKIND
     &             + kz50    *  O2II        / 32._DKIND
            endif
          dspdot(1) = dspdot(1) + (scoef - acoef * de)


!           7) H2

            scoef = 2._DKIND*(k8     * HM          * HI       
     &            +       k10    * H2II        * HI       /2._DKIND
     &            +       k19    * H2II        * HM       /2._DKIND
     &            +       k22    * HI        * (HI       )**2._DKIND)
            acoef = ( k13   *HI        + k11   *HII       
     &              + k12   *de        )
     &              + k29shield    + k31shield   

            if (anydust) then
            if(itmask_metal   ) then
               scoef = scoef + 2._DKIND * h2dust    *
     &              HI        * rhoH   
            endif
            endif
!           contribution of minor species
            if (ispecies .gt. 2) then
               scoef = scoef + 2._DKIND * (
     &                 k53    * HDI        * HII        / 3._DKIND
     &               + k55    * HDI        * HI         / 3._DKIND
     &                  )
               acoef = acoef
     &               + k52    * DII        / 2._DKIND
     &               + k54    * DI         / 2._DKIND
            endif

            if ( (imchem .eq. 1) .and. (itmask_metal   ) ) then
               scoef = scoef +  2._DKIND * ( 0._DKIND
     &             + kz15    *    HI        *    CH        / 13._DKIND
     &             + kz16    *    HI        *   CH2        / 14._DKIND
     &             + kz17    *    HI        *    OH        / 17._DKIND
     &             + kz18    *    HI        *   H2O        / 18._DKIND
     &             + kz47    * H2OII        *    de        / 18._DKIND
     &            )
               acoef = acoef
     &             + kz20    *    CI        / 12._DKIND
     &             + kz21    *    OI        / 16._DKIND
     &             + kz23    *    CH        / 13._DKIND
     &             + kz24    *    OH        / 17._DKIND
     &             + kz40    *   OII        / 16._DKIND
     &             + kz41    *  OHII        / 17._DKIND
     &             + kz42    * H2OII        / 18._DKIND
     &             + kz51    *    CI        / 12._DKIND
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 0) then
                  scoef = scoef + 2._DKIND *
     &                  kdMgSiO3      * 2._DKIND

               endif
               if (idspecies .gt. 1) then
                  scoef = scoef + 2._DKIND * (
     &                  kdMg2SiO4     * 3._DKIND
     &                + kdFe3O4       * 4._DKIND
     &                + kdMgO        
     &                + kdAl2O3       * 3._DKIND
     &              )
               endif
               if (idspecies .gt. 2) then
                    acoef = acoef
     &            + kdvolorg      / H2I        * 2._DKIND * 2._DKIND
               endif
               endif
            endif
          dspdot(8) = dspdot(8) + (scoef - acoef * H2I)


!           8) H-

            scoef = k7    * HI        * de       
            acoef = (k8     + k15   )  * HI        +
     &              (k16    + k17   )  * HII        +
     &              k14    * de        + k19    * H2II       /2.0 +
     &              k27
!           contribution of minor species
            if (ispecies .gt. 2) then
               acoef = acoef
     &               + k56    * DI         / 2._DKIND
            endif

            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k136    *    DM        *    HI        /  2._DKIND
               acoef = acoef
     &             + k135    *    DI        /  2._DKIND
            endif
          dspdot(7) = dspdot(7) + (scoef - acoef * HM)



!           9) H2+

            scoef =    2._DKIND*( k9    *HI    *HII    
     &                    +   k11   *H2I    /2._DKIND*HII    
     &                    +   k17   *HM    *HII    
     &                    + k29shield   *H2I    /2._DKIND
     &                    )
            acoef =         k10   *HI     + k18   *de    
     &                    + k19   *HM    
     &                    + (k28shield   +k30shield   )
            if (ispecies .gt. 3) then
               scoef = scoef +  2._DKIND * ( 0._DKIND
     &             + k152    * HeHII        *    HI        /  5._DKIND
     &            )
               acoef = acoef
     &             + k150    *   HeI        /  4._DKIND
            endif
          dspdot(9) = dspdot(9) + (scoef - acoef * H2II)









      endif ! (ispecies .gt. 1) then

!  --- (D) Now do extra 3-species for molecular HD ---
!     
      if (ispecies .gt. 2) then


!     
!         1) DI
!     
            scoef =   (       k2    * DII        * de       
     &                 +      k51   * DII        * HI       
     &                 + 2._DKIND*k55   * HDI        *
     &              HI       /3._DKIND
     &                 )
            acoef  =    k1    * de       
     &             +    k50    * HII       
     &             +    k54    * H2I       /2._DKIND
     &             +    k56    * HM       
     &             + k24shield   
            if (iradtrans .eq. 1) acoef = acoef + kphHI       
            if (ispecies .gt. 3) then
               scoef = scoef +  2._DKIND * ( 0._DKIND
     &             + k131    *  HDII        *    de        /  3._DKIND
     &             + k133    *   DII        *    DM        /  2._DKIND
     &             + k134    *   HII        *    DM        /  2._DKIND
     &             + k136    *    DM        *    HI        /  2._DKIND
     &             )
               acoef = acoef
     &             + k129    *   HII       
     &             + k132    *    de       
     &             + k135    *    HM       
            endif
            dspdot(10) = dspdot(10) + (scoef - acoef * DI)
                                                    

!         2) DII
c 
            scoef =   (   k1     * DI        * de       
     &            +       k50    * HII       * DI       
     &            +  2._DKIND*k53    * HII       * HDI       /3._DKIND
     &            )
     &            + k24shield   *DI       
            acoef = 0._DKIND
            !! initialize GC202002
            if (iradtrans .eq. 1) scoef = scoef + kphHI       *DI       
            acoef =    k2     * de       
     &            +    k51    * HI       
     &            +    k52    * H2I       /2._DKIND
            if (ispecies .gt. 3) then
               acoef = acoef
     &             + k130    *    HI       
     &             + k133    *    DM        /  2._DKIND
            endif
            dspdot(11) = dspdot(11) + (scoef - acoef * DII)


!          3) HDI
c 
            scoef = 3._DKIND*(k52    * DII       *
     &           H2I       /2._DKIND/2._DKIND
     &           + k54    * DI        * H2I       /2._DKIND/2._DKIND
!!   &           + 2._DKIND*k56    * DI        * HM       /2._DKIND
                       !! corrected by GC202005
     &           +          k56    * DI        * HM       /2._DKIND
     &                 )
            acoef  =    k53    * HII       
     &             +    k55    * HI       
            if (ispecies .gt. 3) then
               scoef = scoef +  3._DKIND * ( 0._DKIND
     &             + k125    *  HDII        *    HI        /  3._DKIND
     &             + k137    *    DM        *    HI        /  2._DKIND
     &             )
            endif
            dspdot(12) = dspdot(12) + (scoef - acoef * HDI)




      endif

!  --- (D2) Now do extra 3-species for minor primordial species ---
!     
      if (ispecies .gt. 3) then



!         1) DM

            scoef =
     &            k132    *    DI        *    de
     &          + k135    *    HM        *    DI
            acoef =
     &            k133    *   DII        /  2._DKIND
     &          + k134    *   HII       
     &          + k136    *    HI       
     &          + k137    *    HI       

            dspdot(13) = dspdot(13) + (scoef - acoef * DM)


!         2) HDII

            scoef = 3._DKIND * (
     &            k129    *    DI        *   HII        /  2._DKIND
     &          + k130    *   DII        *    HI        /  2._DKIND
     &         )
            acoef = 
     &            k125    *    HI       
     &          + k131    *    de       

            dspdot(14) = dspdot(14) + (scoef - acoef * HDII)


!         3) HeHII

            scoef = 5._DKIND * (
     &            k148    *   HeI        *   HII        /  4._DKIND
     &          + k149    *   HeI        *   HII        /  4._DKIND
     &          + k150    *   HeI        *  H2II        /  8._DKIND
     &          + k151    *  HeII        *    HI        /  4._DKIND
     &         )
            acoef =
     &            k152    *    HI       
     &          + k153    *    de       

            dspdot(15) = dspdot(15) + (scoef - acoef * HeHII)




      endif

!  --- (D3) Now do metal species ---
!     
      if (imchem .eq. 1) then

            if (itmask_metal   ) then

C***** CI **********
               scoef = 0._DKIND + 12._DKIND * ( 0._DKIND
     &             + kz15    *    HI        *    CH        / 13._DKIND
     &             + kz44    *   CII        *    de        / 12._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz20    *   H2I        /  2._DKIND
     &             + kz27    *    HI       
     &             + kz28    *    OH        / 17._DKIND
     &             + kz29    *    O2        / 32._DKIND
     &             + kz51    *   H2I        /  2._DKIND
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 0) then
                    acoef = acoef 
     &            + kdAC          / CI        * 12._DKIND      
               endif
               endif

            dspdot(16) = dspdot(16) + (scoef - acoef * CI)



C***** CII **********
               scoef = 0._DKIND + 12._DKIND * ( 0._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz37    *    OH        / 17._DKIND
     &             + kz38    *    O2        / 32._DKIND
     &             + kz44    *    de       

            dspdot(17) = dspdot(17) + (scoef - acoef * CII)



C***** CO **********
               scoef = 0._DKIND + 28._DKIND * ( 0._DKIND
     &             + kz28    *    CI        *    OH        / 204._DKIND
     &             + kz29    *    CI        *    O2        / 384._DKIND
     &             + kz32    *    OI        *    CH        / 208._DKIND
     &             + kz38    *   CII        *    O2        / 384._DKIND
     &             + kz43    *  COII        *    HI        / 28._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz26    *    OH        / 17._DKIND
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 2) then
                    acoef = acoef
     &            + kdreforg      / CO        * 17._DKIND * 0.5_DKIND
     &            + kdvolorg      / CO        * 17._DKIND
               endif
               endif

            dspdot(18) = dspdot(18) + (scoef - acoef * CO)



C***** CO2 **********
               scoef = 0._DKIND + 44._DKIND * ( 0._DKIND
     &             + kz26    *    OH        *    CO        / 476._DKIND
     &            )
               acoef = 0._DKIND

            dspdot(19) = dspdot(19) + (scoef - acoef * CO2)



C***** OI **********
               scoef = 0._DKIND + 16._DKIND * ( 0._DKIND
     &             + kz17    *    HI        *    OH        / 17._DKIND
     &             + kz19    *    HI        *    O2        / 32._DKIND
     &             + kz25    *    OH        *    OH        / 289._DKIND
     &             + kz29    *    CI        *    O2        / 384._DKIND
     &             + kz39    *   OII        *    HI        / 16._DKIND
     &             + kz45    *   OII        *    de        / 16._DKIND
     &             + kz47    * H2OII        *    de        / 18._DKIND
     &             + kz50    *  O2II        *    de        / 16._DKIND
     &             + kz53    *   SiI        *    O2        / 896._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz21    *   H2I        /  2._DKIND
     &             + kz22    *   HII       
     &             + kz30    *    HI       
     &             + kz31    *    OI        / 8._DKIND
     &             + kz32    *    CH        / 13._DKIND
     &             + kz33    *    OH        / 17._DKIND

            dspdot(20) = dspdot(20) + (scoef - acoef * OI)



C***** OH **********
               scoef = 0._DKIND + 17._DKIND * ( 0._DKIND
     &             + kz18    *    HI        *   H2O        / 18._DKIND
     &             + kz19    *    HI        *    O2        / 32._DKIND
     &             + kz21    *    OI        *   H2I        / 32._DKIND
     &             + kz30    *    OI        *    HI        / 16._DKIND
     &             + kz46    * H2OII        *    de        / 18._DKIND
     &             + kz49    * H3OII        *    de        / 19._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz17    *    HI       
     &             + kz24    *   H2I        /  2._DKIND
     &             + kz25    *    OH        / 8.5_DKIND
     &             + kz26    *    CO        / 28._DKIND
     &             + kz28    *    CI        / 12._DKIND
     &             + kz33    *    OI        / 16._DKIND
     &             + kz34    *   HII       
     &             + kz37    *   CII        / 12._DKIND
     &             + kz52    *   SiI        / 28._DKIND
     &             + kz54    *  SiOI        / 44._DKIND

            dspdot(21) = dspdot(21) + (scoef - acoef * OH)



C***** H2O **********
               scoef = 0._DKIND + 18._DKIND * ( 0._DKIND
     &             + kz24    *   H2I        *    OH        / 34._DKIND
     &             + kz25    *    OH        *    OH        / 289._DKIND
     &             + kz48    * H3OII        *    de        / 19._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz18    *    HI       
     &             + kz35    *   HII       
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 0) then
                    acoef = acoef
     &            + kdMgSiO3      / H2O        * 18._DKIND * 2._DKIND
               endif
               if (idspecies .gt. 1) then
                    acoef = acoef
     &            + kdMg2SiO4     / H2O        * 18._DKIND * 3._DKIND
     &            + kdFe3O4       / H2O        * 18._DKIND * 4._DKIND
     &            + kdMgO         / H2O        * 18._DKIND
     &            + kdAl2O3       / H2O        * 18._DKIND * 3._DKIND
               endif
               if (idspecies .gt. 2) then
                    acoef = acoef
     &            + kdH2Oice      / H2O        * 18._DKIND
               endif
               endif

            dspdot(22) = dspdot(22) + (scoef - acoef * H2O)



C***** O2 **********
               scoef = 0._DKIND + 32._DKIND * ( 0._DKIND
     &             + kz31    *    OI        *    OI        / 256._DKIND
     &             + kz33    *    OI        *    OH        / 272._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz19    *    HI       
     &             + kz29    *    CI        / 12._DKIND
     &             + kz36    *   HII       
     &             + kz38    *   CII        / 12._DKIND
     &             + kz53    *   SiI        / 28._DKIND

            dspdot(23) = dspdot(23) + (scoef - acoef * O2)



C***** SiI **********
               scoef = 0._DKIND + 28._DKIND * ( 0._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz52    *    OH        / 17._DKIND
     &             + kz53    *    O2        / 32._DKIND
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 1) then
                     acoef = acoef 
     &               + kdSiM         / SiI        * 28._DKIND
               endif
               endif

            dspdot(24) = dspdot(24) + (scoef - acoef * SiI)



C***** SiOI **********
               scoef = 0._DKIND + 44._DKIND * ( 0._DKIND
     &             + kz52    *   SiI        *    OH        / 476._DKIND
     &             + kz53    *   SiI        *    O2        / 896._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz54    *    OH        / 17._DKIND
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 0) then
                       acoef = acoef 
     &               + kdMgSiO3      / SiOI        * 44._DKIND
               endif
               if (idspecies .gt. 1) then
                       acoef = acoef 
     &               + kdMg2SiO4     / SiOI        * 44._DKIND
               endif
               endif

            dspdot(25) = dspdot(25) + (scoef - acoef * SiOI)



C***** SiO2I **********
               scoef = 0._DKIND + 60._DKIND * ( 0._DKIND
     &             + kz54    *  SiOI        *    OH        / 748._DKIND
     &            )
               acoef = 0._DKIND
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 1) then
                       acoef = acoef 
     &               + kdSiO2D       / SiO2I        * 60._DKIND
               endif
               endif

            dspdot(26) = dspdot(26) + (scoef - acoef * SiO2I)



C***** CH **********
               scoef = 0._DKIND + 13._DKIND * ( 0._DKIND
     &             + kz16    *    HI        *   CH2        / 14._DKIND
     &             + kz20    *    CI        *   H2I        / 24._DKIND
     &             + kz27    *    CI        *    HI        / 12._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz15    *    HI       
     &             + kz23    *   H2I        /  2._DKIND
     &             + kz32    *    OI        / 16._DKIND

            dspdot(27) = dspdot(27) + (scoef - acoef * CH)



C***** CH2 **********
               scoef = 0._DKIND + 14._DKIND * ( 0._DKIND
     &             + kz23    *   H2I        *    CH        / 26._DKIND
     &             + kz51    *   H2I        *    CI        / 24._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz16    *    HI       
               if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
               if (idspecies .gt. 2) then
                    acoef = acoef 
     &            + kdreforg      / CH2        * 14._DKIND * 0.5_DKIND
               endif
               endif

            dspdot(28) = dspdot(28) + (scoef - acoef * CH2)



C***** COII **********
               scoef = 0._DKIND + 28._DKIND * ( 0._DKIND
     &             + kz37    *   CII        *    OH        / 204._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz43    *    HI       

            dspdot(29) = dspdot(29) + (scoef - acoef * COII)



C***** OII **********
               scoef = 0._DKIND + 16._DKIND * ( 0._DKIND
     &             + kz22    *   HII        *    OI        / 16._DKIND
     &             + kz38    *   CII        *    O2        / 384._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz39    *    HI       
     &             + kz40    *   H2I        /  2._DKIND
     &             + kz45    *    de       

            dspdot(30) = dspdot(30) + (scoef - acoef * OII)



C***** OHII **********
               scoef = 0._DKIND + 17._DKIND * ( 0._DKIND
     &             + kz34    *   HII        *    OH        / 17._DKIND
     &             + kz40    *   OII        *   H2I        / 32._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz41    *   H2I        /  2._DKIND

            dspdot(31) = dspdot(31) + (scoef - acoef * OHII)



C***** H2OII **********
               scoef = 0._DKIND + 18._DKIND * ( 0._DKIND
     &             + kz35    *   HII        *   H2O        / 18._DKIND
     &             + kz41    *  OHII        *   H2I        / 34._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz42    *   H2I        /  2._DKIND
     &             + kz46    *    de       
     &             + kz47    *    de       

            dspdot(32) = dspdot(32) + (scoef - acoef * H2OII)



C***** H3OII **********
               scoef = 0._DKIND + 19._DKIND * ( 0._DKIND
     &             + kz42    * H2OII        *   H2I        / 36._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz48    *    de       
     &             + kz49    *    de       

            dspdot(33) = dspdot(33) + (scoef - acoef * H3OII)



C***** O2II **********
               scoef = 0._DKIND + 32._DKIND * ( 0._DKIND
     &             + kz36    *   HII        *    O2        / 32._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz50    *    de       

            dspdot(34) = dspdot(34) + (scoef - acoef * O2II)



            if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then
            if (idspecies .gt. 0) then
C***** Mg **********
               scoef = 0._DKIND
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kdMgSiO3      / Mg        * 24._DKIND
               if (idspecies .gt. 1) then
                    acoef = acoef 
     &            + kdMg2SiO4     / Mg        * 24._DKIND * 2._DKIND
     &            + kdMgO         / Mg        * 24._DKIND
               endif

            dspdot(35) = dspdot(35) + (scoef - acoef * Mg)


            endif

            if (idspecies .gt. 1) then
C***** Al **********
               scoef = 0._DKIND 
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kdAl2O3       / Al        * 27._DKIND * 2._DKIND

            dspdot(36) = dspdot(36) + (scoef - acoef * Al)



C***** S  **********
               scoef = 0._DKIND 
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kdFeS         / S        * 32._DKIND

            dspdot(37) = dspdot(37) + (scoef - acoef * S)



C***** Fe **********
               scoef = 0._DKIND 
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kdFeM         / Fe        * 56._DKIND
     &         + kdFe3O4       / Fe        * 56._DKIND * 3._DKIND
     &         + kdFeS         / Fe        * 56._DKIND

            dspdot(38) = dspdot(38) + (scoef - acoef * Fe)


            endif
            endif

            endif ! itmask_metal

      endif ! imchem

!  --- (D4) Now do dust species ---
!     
      if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1) ) then

            if (itmask_metal   ) then

            if (idspecies .gt. 0) then
C***** MgSiO3 **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdMgSiO3      * 100._DKIND
               acoef = 0._DKIND

            dspdot(39) = dspdot(39) + (scoef - acoef * MgSiO3)



C***** AC **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdAC          * 12._DKIND
               acoef = 0._DKIND

            dspdot(40) = dspdot(40) + (scoef - acoef * AC)


            endif

            if (idspecies .gt. 1) then
C***** SiM **********
               scoef = 0._DKIND
                 scoef = scoef 
     &         + kdSiM         * 28._DKIND
               acoef = 0._DKIND

            dspdot(41) = dspdot(41) + (scoef - acoef * SiM)



C***** FeM **********
               scoef = 0._DKIND
                 scoef = scoef 
     &         + kdFeM         * 56._DKIND
               acoef = 0._DKIND

            dspdot(42) = dspdot(42) + (scoef - acoef * FeM)



C***** Mg2SiO4 **********
               scoef = 0._DKIND
                 scoef = scoef 
     &         + kdMg2SiO4     * 140._DKIND
               acoef = 0._DKIND

            dspdot(43) = dspdot(43) + (scoef - acoef * Mg2SiO4)



C***** Fe3O4 **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdFe3O4       * 232._DKIND
               acoef = 0._DKIND

            dspdot(44) = dspdot(44) + (scoef - acoef * Fe3O4)



C***** SiO2D **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdSiO2D       * 60._DKIND
               acoef = 0._DKIND

            dspdot(45) = dspdot(45) + (scoef - acoef * SiO2D)



C***** MgO **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdMgO         * 40._DKIND
               acoef = 0._DKIND

            dspdot(46) = dspdot(46) + (scoef - acoef * MgO)



C***** FeS **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdFeS         * 88._DKIND
               acoef = 0._DKIND

            dspdot(47) = dspdot(47) + (scoef - acoef * FeS)



C***** Al2O3 **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdAl2O3       * 102._DKIND
               acoef = 0._DKIND

            dspdot(48) = dspdot(48) + (scoef - acoef * Al2O3)


            endif

            if (idspecies .gt. 2) then
C***** reforg **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdreforg      * 22.68_DKIND
               acoef = 0._DKIND

            dspdot(49) = dspdot(49) + (scoef - acoef * reforg)



C***** volorg **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdvolorg      * 32._DKIND
               acoef = 0._DKIND

            dspdot(50) = dspdot(50) + (scoef - acoef * volorg)



C***** H2Oice **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kdH2Oice      * 18._DKIND
               acoef = 0._DKIND

            dspdot(51) = dspdot(51) + (scoef - acoef * H2Oice)


            endif

            endif ! itmask_metal

      endif ! igrgr or idsub

      return
      end
