#include "phys_const.def"

!=======================================================================
!/////////////////////  SUBROUTINE COOL1D_MULTI_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_multi_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool,
     &                idust, idustall, idustfield, idustrec,
     &                idim, is, ie, j, k, ih2co, ipiht, iter, igammah,
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
     &                iradshield, avgsighi, avgsighei,
     &                avgsigheii,
     &                k24, k26, iradtrans, photogamma,
     &                ih2optical, iciecool, ih2cr, ihdcr, ciecoa, cieco,
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat,
     &                iTfloor, Tfloor_scalar, Tfloor,
     &                iisrffield, isrf_habing,
     &                itmask, itmask_metal
     &              , imchem, igrgr, ipcont, tmcool
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
     &            )

!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
!  modified6: September, 2009 by BDS to include cloudy cooling
!
!  PURPOSE:
!    Solve the energy cooling equations.
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!  PARAMETERS:
!

!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!  Arguments

      integer in, jn, kn, is, ie, j, k, nratec, idim,
     &        iexpand, ih2co, ipiht, ispecies, imcool,
     &        idust, idustall, idustfield, idustrec,
     &        imetal, igammah, ih2optical, iciecool, ih2cr, ihdcr,
     &        clnew,
     &        iVheat, iMheat, iTfloor, iradtrans, iradshield,
     &        iisrffield
     &       ,imchem, igrgr, ipcont

      real*8  aye, temstart, temend, z_solar, fgr,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma, fh, clEleFra, Tfloor_scalar, tmcool
      R_PREC  d(in,jn,kn),    e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn), dust(in,jn,kn),
     &        Vheat(in,jn,kn), Mheat(in,jn,kn), Tfloor(in,jn,kn),
     &        photogamma(in,jn,kn), isrf_habing(in,jn,kn)
     &      , DM(in,jn,kn)   , HDII(in,jn,kn)   , HeHII(in,jn,kn)
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
     &        gammaha, isrf, regra(nratec), gamma_isrfa,
     &        avgsighi, avgsighei, avgsigheii, k24, k26

!  Cloudy cooling data

      integer icmbTfloor, iClHeat, iZscale, mycmbTfloor
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

      real*8, parameter :: mh = mass_h      !DPC
      ! approx. mean molecular weight of metals
      real*8, parameter :: mu_metal = 16._DKIND
      integer, parameter :: ti_max = 20

!  Locals

      integer i, iter, ti, iradfield
      real*8 dom, qq, vibl, logtem0, logtem9, dlogtem, zr,
     &       comp1, comp2, 
     &       hdlte1, hdlow1, gamma2, x, fudge, fH2,
     &       gphdl1, dom_inv, tau, ciefudge,
     &       coolunit, dbase1, tbase1, xbase1,
     &       nH2, nother, nSSh, nratio, nssh_he, nratio_he,
     &       fSShHI, fSShHeI, pe_eps, pe_X, grbeta, ih2cox,
     &       min_metallicity
      integer itd

!  Slice locals
 
      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &       tgas(in), tgasold(in), mmw(in), tdust(in), rhoH(in),
     &       mynh(in), metallicity(in), dust2gas(in), edot(in),
     &       myde(in), gammaha_eff(in)

!  Cooling/heating slice locals

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), cieco(in),
     &       hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in),
     &       gaHI(in), gaH2(in), gaHe(in), gaHp(in), gael(in),
     &       h2lte(in), galdl(in), gasgr(in), gasgr_tdust(in),
     &       regr(in), myisrf(in)
       integer iden, item, itab
       real*8  cieY06a(nratec), cieY06(in)
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
       real*8 logdom, logT(in), logTcmb(in), logrho(in)
     &      , logH(in), logH2I(in), logHDI(in)
     &      , logH2(in), logCI(in), logCII(in), logOI(in)
     &      , logCO(in), logOH(in), logH2O(in)
       real*8 logdvdr(in), lognhat
       real*8 log_Linv, log_Ginv, L, G
       real*8 Lpri(in)
     &      , LH2(in), LCIE(in), LHD(in)
     &      , LCI(in), LCII(in), LOI(in)
     &      , LCO(in), LOH(in), LH2O(in)
       real*8 Ldst(in), alpha(in), alphad(in)
       real*8 lshield_con(in), tau_con(in)
       real*8 log_a, log_L, log_G
!     grain growth
      real*8  sgSiM(in), sgFeM(in), sgMg2SiO4(in)
     &      , sgMgSiO3(in), sgFe3O4(in), sgAC(in)
     &      , sgSiO2D(in), sgMgO(in), sgFeS(in)
     &      , sgAl2O3(in)
     &      , sgreforg(in), sgvolorg(in), sgH2Oice(in)
     &      , sgtot(in)
      real*8  alSiM(gr_N(2),in), alFeM(gr_N(2),in)
     &      , alMg2SiO4(gr_N(2),in), alMgSiO3(gr_N(2),in)
     &      , alFe3O4(gr_N(2),in), alAC(gr_N(2),in)
     &      , alSiO2D(gr_N(2),in), alMgO(gr_N(2),in)
     &      , alFeS(gr_N(2),in), alAl2O3(gr_N(2),in)
     &      , alreforg(gr_N(2),in)
     &      , alvolorg(gr_N(2),in), alH2Oice(gr_N(2),in)
     &      , altot(gr_N(2),in)
      real*8  kpSiM(in), kpFeM(in), kpMg2SiO4(in)
     &      , kpMgSiO3(in), kpFe3O4(in), kpAC(in)
     &      , kpSiO2D(in), kpMgO(in), kpFeS(in)
     &      , kpAl2O3(in)
     &      , kpreforg(in), kpvolorg(in), kpH2Oice(in)
     &      , kptot(in)
!     grain temperature
      real*8  tSiM(in), tFeM(in), tMg2SiO4(in)
     &      , tMgSiO3(in), tFe3O4(in), tAC(in)
     &      , tSiO2D(in), tMgO(in), tFeS(in)
     &      , tAl2O3(in)
     &      , treforg(in), tvolorg(in), tH2Oice(in)
      real*8  gasgr2a(nratec), gamma_isrf2a
      real*8  gasSiM(in), gasFeM(in), gasMg2SiO4(in)
     &      , gasMgSiO3(in), gasFe3O4(in), gasAC(in)
     &      , gasSiO2D(in), gasMgO(in), gasFeS(in)
     &      , gasAl2O3(in)
     &      , gasreforg(in), gasvolorg(in), gasH2Oice(in)
      real*8  gasgr_tSiM(in), gasgr_tFeM(in), gasgr_tMg2SiO4(in)
     &      , gasgr_tMgSiO3(in), gasgr_tFe3O4(in), gasgr_tAC(in)
     &      , gasgr_tSiO2D(in), gasgr_tMgO(in), gasgr_tFeS(in)
     &      , gasgr_tAl2O3(in)
     &      , gasgr_treforg(in), gasgr_tvolorg(in), gasgr_tH2Oice(in)
      real*8  mygisrf(in), fv2k, fac
      real*8  gisrfSiM(in), gisrfFeM(in), gisrfMg2SiO4(in)
     &      , gisrfMgSiO3(in), gisrfFe3O4(in), gisrfAC(in)
     &      , gisrfSiO2D(in), gisrfMgO(in), gisrfFeS(in)
     &      , gisrfAl2O3(in)
     &      , gisrfreforg(in), gisrfvolorg(in), gisrfH2Oice(in)
!  Iteration mask

      MASK_TYPE itmask(in), anydust, interp
      MASK_TYPE itmask_metal(in), itmask_tab(in)
!!#define CALCULATE_TGAS_SELF_CONSISTENTLY
!#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
      integer iter_tgas
      real*8  tgas_err, tgas0
!#endif /* NOT important */
!     debug
      real*8  edotunit
      integer i_max
      R_PREC  d_max

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set flag for dust-related options

      if ((idust .gt. 0) .or. (idustall .gt. 0) .or.
     &    (idustrec .gt. 0)) then
         anydust = MASK_TRUE
      else
         anydust = MASK_FALSE
      endif

!     Set flag for needing interpolation variables

      if ((ispecies .gt. 0) .or. (idustall .gt. 0)) then
         interp = MASK_TRUE
      else
         interp = MASK_FALSE
      endif
!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1, DKIND)

!     Set units

      dom      = urho*(aye**3)/mh
      dom_inv  = 1._DKIND/dom
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1._DKIND/(aye*uaye) - 1._DKIND
      fudge    = 1._DKIND
      iradfield = -1

!     Set compton cooling coefficients (and temperature)

      comp1 = compa  * (1._DKIND + zr)**4
      comp2 = 2.73_DKIND * (1._DKIND + zr)

!     multiplicative factor for including/excluding H2 cooling
      ih2cox = real(ih2co, DKIND)

!     ignore metal chemistry/cooling below this metallicity
      min_metallicity = 1.d-9 / z_solar

!     Initialize edot

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            edot(i) = 0._DKIND
         end if
      enddo

!     Compute Pressure

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            p2d(i) = (gamma - 1._DKIND)*d(i,j,k)*e(i,j,k)
         end if
      enddo

!     Compute Temperature

!     If no chemistry, use a tabulated mean molecular weight
!     and iterate to convergence.

      if (ispecies .eq. 0) then

!     fh is H mass fraction in metal-free gas.

         if (imetal .eq. 1) then
            do i = is+1, ie+1
               if ( itmask(i) .ne. MASK_FALSE ) then
                  rhoH(i) = fh * (d(i,j,k) - metal(i,j,k))
               endif
            enddo
         else
            do i = is+1, ie+1
               if ( itmask(i) .ne. MASK_FALSE ) then
                  rhoH(i) = fh * d(i,j,k)
               endif
            enddo
         endif

         call calc_temp1d_cloudy_g(d, metal, e, rhoH,
     &        in, jn, kn, is, ie, j, k,
     &        tgas, mmw, dom, zr, 
     &        temstart, temend,
     &        gamma, utem, imetal,
     &        priGridRank, priGridDim,
     &        priPar1, priPar2, priPar3,
     &        priDataSize, priMMW,
     &        itmask)

      else

!     Compute mean molecular weight (and temperature) directly

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               mmw(i) =
     &              (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4._DKIND +
     &              HI(i,j,k) + HII(i,j,k) + de(i,j,k)
               rhoH(i) = HI(i,j,k) + HII(i,j,k)
               myde(i) = de(i,j,k)
            end if
         enddo

!     (include molecular hydrogen, but ignore deuterium)

         if (ispecies .gt. 1) then
            do i = is+1, ie+1
               if ( itmask(i) .ne. MASK_FALSE ) then
                  mmw(i) = mmw(i) +
     &                 HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2._DKIND
                  rhoH(i) = rhoH(i) + H2I(i,j,k) + H2II(i,j,k)
               end if
            enddo
         endif

!     Include metal species

         if (imetal .eq. 1) then
            do i = is+1, ie+1
               if ( itmask(i) .ne. MASK_FALSE ) then
                  mmw(i) = mmw(i) + metal(i,j,k)/mu_metal
               end if
            enddo
         endif

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               tgas(i) = max(p2d(i)*utem/mmw(i), temstart)
               mmw(i) = d(i,j,k) / mmw(i)
            end if
         enddo

!     Correct temperature for gamma from H2

         if (ispecies .gt. 1) then
            do i = is+1, ie+1
               if ( itmask(i) .ne. MASK_FALSE ) then
                  nH2 = 0.5_DKIND*(H2I(i,j,k) + H2II(i,j,k))
                  nother = (HeI(i,j,k) + HeII(i,j,k) +
     &                 HeIII(i,j,k))/4._DKIND +
     &                 HI(i,j,k) + HII(i,j,k) + de(i,j,k)

                  iter_tgas = 0
                  tgas_err = huge8
                  do while ((iter_tgas .lt. 100)
     &                 .and.(tgas_err .gt. 1.d-3))
                     tgas0 = tgas(i)
                     if (nH2/nother .gt. 1.0e-3_DKIND) then
                        x = 6100._DKIND/tgas(i) ! not quite self-consistent
                        if (x .gt. 10._DKIND) then
                           gamma2 = 0.5_DKIND*5._DKIND
                        else
                           gamma2 = 0.5_DKIND*(5._DKIND + 2._DKIND*x**2
     &                       * exp(x)/(exp(x)-1)**2)
                        endif
                     else
                        gamma2 = 2.5_DKIND
                     endif
                     gamma2 = 1._DKIND + (nH2 + nother)/
     &                    (nH2*gamma2 + nother/(gamma-1._DKIND))
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
                     tgas(i) = max((gamma2 - 1._DKIND)*mmw(i)*e(i,j,k)*
     &                    utem, temstart)
                     tgas_err = dabs(tgas0 - tgas(i)) / tgas0
                     iter_tgas = iter_tgas + 1
#else
                     tgas(i) = tgas(i) * (gamma2 - 1._DKIND)/
     &                    (gamma - 1._DKIND)
                     iter_tgas = 101
#endif
                  end do
               end if
            enddo
         endif

      endif

!     Skip if below the temperature floor

      if (iTfloor .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               if (tgas(i) .le. Tfloor_scalar) then
                  edot(i) = tiny
                  itmask(i) = MASK_FALSE
               endif
            endif
         enddo
      else if (iTfloor .eq. 2) then
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               if (tgas(i) .le. Tfloor(i,j,k)) then
                  edot(i) = tiny
                  itmask(i) = MASK_FALSE
               endif
            endif
         enddo
      endif

!     Calculate metallicity and H number density

      if (imetal .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar
            endif
         enddo
      else
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               metallicity(i) = tiny
            endif
         enddo
      endif

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            mynh(i) = rhoH(i) * dom
         end if
      enddo

!     If this is the first time through, just set tgasold to tgas

      if (iter .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
            tgasold(i) = tgas(i)
            end if
         enddo
      endif

!     Compute log densities

      logdom = log10(dom)
      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            logT(i)   = log10(tgas(i))
            if(icmbTfloor .eq. 1)
     &         logTcmb(i) = log10(comp2)
            logrho(i) = log10(d(i,j,k) * dom*mh)
            if( ispecies .gt. 0) then
               logH(i)   = log10(HI(i,j,k) * dom)
               logH2(i)  = log10(HI(i,j,k) * dom)
            endif
            if( ispecies .gt. 1 ) then
              logH2(i)  = log10((HI(i,j,k) + H2I(i,j,k) / 2.d0) * dom)
              logH2I(i) = log10(H2I(i,j,k) * dom /  2.d0)
            endif
            if( ispecies .gt. 2) then
              logHDI(i) = log10(HDI(i,j,k) * dom /  3.d0)
            endif
            if( imcool .eq. 1 ) then
            if( imchem .eq. 1 ) then
              logCI(i)  = log10(CI (i,j,k) * dom / 12.d0)
              logCII(i) = log10(CII(i,j,k) * dom / 12.d0)
              logOI(i)  = log10(OI (i,j,k) * dom / 16.d0)
              logCO(i)  = log10(CO (i,j,k) * dom / 28.d0)
              logOH(i)  = log10(OH (i,j,k) * dom / 17.d0)
              logH2O(i) = log10(H2O(i,j,k) * dom / 18.d0)
            endif
            endif

!     From Chiaki & Wise (2019), approximate dv/dr as 1/(3 * t_ff)
            logdvdr(i) = -8.79947961814d0 + 0.5d0 * logrho(i) ! km/s / cm
            lshield_con(i) = 
     &         sqrt((gamma * pi_val * kboltz * tgas(i)) /
     &         (GravConst *  mmw(i)*mh * d(i,j,k)*dom*mh))

          endif
      enddo
     
      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then

!        Compute log temperature and truncate if above/below table max/min

         logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

         endif
      enddo

!     Compute interpolation indices

      if (interp .ne. MASK_FALSE) then
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then

!     Compute index into the table and precompute parts of linear interp

         indixe(i) = min(nratec-1,
     &           max(1,int((logtem(i)-logtem0)/dlogtem, DIKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

            endif
         enddo
      endif

!     --- 6 species cooling ---

      if (ispecies .gt. 0) then

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then

!        Lookup cooling values and do a linear temperature in log(T)

         ceHI(i) = ceHIa(indixe(i)) + tdef(i)
     &         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))
         ceHeI(i) = ceHeIa(indixe(i)) + tdef(i)
     &         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))
         ceHeII(i) = ceHeIIa(indixe(i)) + tdef(i)
     &         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))
         ciHI(i) = ciHIa(indixe(i)) + tdef(i)
     &         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))
         ciHeI(i) = ciHeIa(indixe(i)) + tdef(i)
     &         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))
         ciHeIS(i) = ciHeISa(indixe(i)) + tdef(i)
     &         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))
         ciHeII(i) = ciHeIIa(indixe(i)) + tdef(i)
     &         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))
         reHII(i) = reHIIa(indixe(i)) + tdef(i)
     &         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))
         reHeII1(i)=reHeII1a(indixe(i)) + tdef(i)
     &        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))
         reHeII2(i)=reHeII2a(indixe(i)) + tdef(i)
     &        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))
         reHeIII(i)=reHeIIIa(indixe(i)) + tdef(i)
     &        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))
         brem(i) = brema(indixe(i)) + tdef(i)
     &         *(brema(indixe(i)+1) -brema(indixe(i)))

         end if
      enddo

!     Compute the cooling function

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
         edot(i) = (

!                    Collisional excitations

     &             - ceHI  (i)*HI  (i,j,k)*de(i,j,k)              ! ce of HI
     &             - ceHeI (i)*HeII(i,j,k)*de(i,j,k)**2*dom/4._DKIND  ! ce of HeI
     &             - ceHeII(i)*HeII(i,j,k)*de(i,j,k)/4._DKIND         ! ce of HeII

!                    Collisional ionizations

     &             - ciHI  (i)*HI  (i,j,k)*de(i,j,k)              ! ci of HI
     &             - ciHeI (i)*HeI (i,j,k)*de(i,j,k)/4._DKIND         ! ci of HeI
     &             - ciHeII(i)*HeII(i,j,k)*de(i,j,k)/4._DKIND         ! ci of HeII
     &             - ciHeIS(i)*HeII(i,j,k)*de(i,j,k)**2*dom/4._DKIND  ! ci of HeIS

!                    Recombinations

     &             - reHII  (i)*HII  (i,j,k)*de(i,j,k)           ! re of HII
     &             - reHeII1(i)*HeII (i,j,k)*de(i,j,k)/4._DKIND      ! re of HeII
     &             - reHeII2(i)*HeII (i,j,k)*de(i,j,k)/4._DKIND      ! re of HeII
     &             - reHeIII(i)*HeIII(i,j,k)*de(i,j,k)/4._DKIND      ! re of HeIII

!                    Bremsstrahlung

     &             - brem(i)*(HII(i,j,k)+HeII(i,j,k)/4._DKIND +
     &           HeIII(i,j,k)) * de(i,j,k)

     &             )
         Lpri(i) = edot(i)

         if (edot(i) .ne. edot(i)) then
#ifdef _OPENMP
!$omp critical
#endif
            write(6,*) 'NaN in edot[1]: ', i, j, k, edot(i), 
     &           HI(i,j,k), HII(i,j,k), HeI(i,j,k), HeII(i,j,k), 
     &           HeIII(i,j,k), de(i,j,k), d(i,j,k), 
     &           tgas(i), p2d(i)
#ifdef _OPENMP
!$omp end critical
#endif
         endif
         
         end if
      enddo

      endif

!     --- H2 cooling ---

      if (ispecies .gt. 1) then
!        Chiaki & Wise (2019) H2 cooling rate
         if (ih2cr .eq. 3) then

      do i = is+1, ie+1
      if ( itmask(i) .ne. MASK_FALSE ) then

         lognhat = logH2I(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int(LH2_N,8),
     &     LH2_D, LH2_dD, LH2_T, LH2_dT, LH2_H, LH2_dH,
     &     int(LH2_Size,8), LH2_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int(LH2_N,8),
     &     LH2_D, LH2_dD, LH2_T, LH2_dT, LH2_H, LH2_dH,
     &     int(LH2_Size,8), LH2_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LH2(i) = ih2cox * (G - L) / dom * H2I(i,j,k)/2.d0
         if (LH2(i) .ne. LH2(i)) LH2(i) = 0.d0
         edot(i) = edot(i) + LH2(i)

      endif
      enddo

!     Glover & Abel (2008) H2 cooling rate
      else if (ih2cr .eq. 2) then

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
            gaHI(i) = gaHIa(indixe(i)) + tdef(i)
     &         *(gaHIa(indixe(i)+1) - gaHIa(indixe(i)))
            gaH2(i) = gaH2a(indixe(i)) + tdef(i)
     &         *(gaH2a(indixe(i)+1) - gaH2a(indixe(i)))
            gaHe(i) = gaHea(indixe(i)) + tdef(i)
     &         *(gaHea(indixe(i)+1) - gaHea(indixe(i)))
            gaHp(i) = gaHpa(indixe(i)) + tdef(i)
     &         *(gaHpa(indixe(i)+1) - gaHpa(indixe(i)))
            gael(i) = gaela(indixe(i)) + tdef(i)
     &         *(gaela(indixe(i)+1) - gaela(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            h2lte(i) = h2ltea(indixe(i)) + tdef(i)
     &         *(h2ltea(indixe(i)+1) - h2ltea(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5_DKIND*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + 
     &           HeIII(i,j,k))/4._DKIND +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40._DKIND * 10._DKIND**(4.8_DKIND * 
     &           sqrt(max(log10(tgas(i)),2._DKIND)-2._DKIND)) / fH2**2)/
     &           ((nH2 + nother)*dom) )
            fudge = min(fudge, 1._DKIND)
#endif /* OPTICAL_DEPTH_FUDGE */
            ! Note that this optical depth approximation comes from
            ! RA04.
            if (ih2optical.eq.1) then
                fudge = (0.76_DKIND*d(i,j,k)*dom/
     &              8.e9_DKIND)**(-0.45_DKIND)
                fudge = min(fudge, 1._DKIND)
            else
                fudge = 1._DKIND
            endif
            galdl(i) = gaHI(i) * HI(i,j,k)  
     &               + gaH2(i) * H2I(i,j,k) / 2._DKIND
     &               + gaHe(i) * HeI(i,j,k) / 4._DKIND
     &               + gaHp(i) * HII(i,j,k)
     &               + gael(i) * de(i,j,k)
c            gphdl1 = gphdl(i)/dom
            gphdl1 = h2lte(i)/dom
            edot(i) = edot(i) - ih2cox*fudge*H2I(i,j,k)*
     &           h2lte(i)/(1._DKIND + gphdl1/galdl(i)) / (2._DKIND*dom)

            end if
         enddo

!     Galli & Palla (1998) H2 cooling rate
      else if (ih2cr .eq. 1) then

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
            gpldl(i) = gpldla(indixe(i)) + tdef(i)
     &         *(gpldla(indixe(i)+1) - gpldla(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then

#define NO_OPTICAL_DEPTH_FUDGE
#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5_DKIND*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) +
     &           HeIII(i,j,k))/4._DKIND +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40._DKIND * 10._DKIND**(4.8_DKIND * 
     &           sqrt(max(log10(tgas(i)),2._DKIND)-2._DKIND)) / fH2**2)/
     &           ((nH2 + nother)*dom) )
            fudge = min(fudge, 1._DKIND)
#endif /* OPTICAL_DEPTH_FUDGE */
            ! Note that this optical depth approximation comes from
            ! RA04.
            if (ih2optical.eq.1) then
                fudge = (0.76_DKIND*d(i,j,k)*dom/
     &              8.e9_DKIND)**(-0.45_DKIND)
                fudge = min(fudge, 1._DKIND)
            else
                fudge = 1._DKIND
            endif
            gphdl1 = gphdl(i)/(HI(i,j,k)*dom)
            edot(i) = edot(i) - ih2cox*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1._DKIND + gphdl1/gpldl(i)) / (2._DKIND*dom)

            end if
         enddo

!      Lepp & Shull (1983) H2 cooling rate
       else if (ih2cr .eq. 0) then

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
            hyd01k(i) = hyd01ka(indixe(i)) + tdef(i)
     &         *(hyd01ka(indixe(i)+1)-hyd01ka(indixe(i)))
            h2k01(i) = h2k01a(indixe(i)) + tdef(i)
     &         *(h2k01a(indixe(i)+1) - h2k01a(indixe(i)))
            vibh(i) = vibha(indixe(i)) + tdef(i)
     &         *(vibha(indixe(i)+1) - vibha(indixe(i)))
            roth(i) = rotha(indixe(i)) + tdef(i)
     &         *(rotha(indixe(i)+1) - rotha(indixe(i)))
            rotl(i) = rotla(indixe(i)) + tdef(i)
     &         *(rotla(indixe(i)+1) - rotla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
            qq   = 1.2_DKIND*(HI(i,j,k)*dom)**0.77_DKIND + 
     &                (H2I(i,j,k)*dom/2._DKIND)**0.77_DKIND
            vibl = (HI(i,j,k)*hyd01k(i) + 
     &             H2I(i,j,k)/2._DKIND*h2k01(i))
     &             *dom*8.18e-13_DKIND

#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5_DKIND*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) +
     &           HeIII(i,j,k))/4._DKIND +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40._DKIND * 10._DKIND**(4.8_DKIND * 
     &           sqrt(max(log10(tgas(i)),2._DKIND)-2._DKIND)) / fH2**2)/
     &           ((nH2 + nother)*dom) )
            fudge = min(fudge, 1._DKIND)
#endif /* OPTICAL_DEPTH_FUDGE */

            edot(i) = edot(i) - ih2cox*fudge*H2I(i,j,k)*(
     &           vibh(i)/(1._DKIND+vibh(i)/max(   vibl     ,tiny)) +
     &           roth(i)/(1._DKIND+roth(i)/max(qq*rotl(i),tiny))     
     &           )/2._DKIND/dom
            end if
         enddo

      endif ! ih2cr

c     CIE
c     cooling from H2-H2 and He-H2 collisional induced emission comes
C     with its own radiative transfer correction as discussed in
C     Ripamonti & Abel 2003
         if (iciecool.eq.1) then
            do i = is+1, ie+1
            if (itmask(i) .ne. MASK_FALSE) then
c     Only calculate if H2I(i) is a substantial fraction
              if (d(i,j,k)*dom.gt.1e10_DKIND) then
                ciefudge = 1._DKIND
                tau = ((d(i,j,k)/2e16_DKIND)*dom)**2.8_DKIND  ! 2e16 is in units of cm^-3
                tau = max(tau, 1.e-5_DKIND)
                ciefudge = min((1._DKIND-exp(-tau))/tau,1._DKIND)
c               Matt's attempt at a second exponentialier cutoff
                tau = ((d(i,j,k)/2.e18_DKIND)*dom)**8._DKIND  ! 2e18 is in units of cm^-3
                tau = max(tau, 1.e-5_DKIND)
                ciefudge = ciefudge*min((1.-exp(-tau))/tau,1._DKIND)
c               ciefudge, which is applied to the continuum, is applied to edot
                edot(i) = ciefudge*(edot(i) - 
     &                  H2I(i,j,k)*(d(i,j,k)*cieco(i)))
              endif
            endif
            enddo
c     CIE H2 cooling using Yoshida et al. (2006)
         else if (iciecool .eq. 2) then
            do i = is+1, ie+1
               if (itmask(i) .ne. MASK_FALSE) then
                cieY06(i) = cieY06a(indixe(i)) + tdef(i)
     &             *(cieY06a(indixe(i)+1) - cieY06a(indixe(i)))
                LCIE(i) = - cieY06(i) * (H2I(i,j,k)/2.d0)**2
                edot(i) = edot(i) + LCIE(i)
               endif
            enddo
         endif

      endif

!     --- Cooling from HD ---

      if (ispecies .gt. 2) then

!     Chiaki & Wise (2019) HD cooling rate
      if (ihdcr .eq. 1 ) then

      do i = is+1, ie+1
      if ( itmask(i) .ne. MASK_FALSE ) then

         lognhat = logHDI(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int(LHD_N,8),
     &     LHD_D, LHD_dD, LHD_T, LHD_dT, LHD_H, LHD_dH,
     &     int(LHD_Size,8), LHD_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int(LHD_N,8),
     &     LHD_D, LHD_dD, LHD_T, LHD_dT, LHD_H, LHD_dH,
     &     int(LHD_Size,8), LHD_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LHD(i) = (G - L) / dom * HDI(i,j,k)/3.d0
         if (LHD(i) .ne. LHD(i)) LHD(i) = 0.d0
         edot(i) = edot(i) + LHD(i)

      endif
      enddo

!     Coppola et al (2011) and Wrathmall, Gusdorf, & Flower (2007) HD cooling rate
      else if (ihdcr .eq. 0) then

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
c CMB cooling floor
               if (tgas(i) .gt. comp2) then
                  hdlte(i) = hdltea(indixe(i)) + tdef(i)
     &            *(hdltea(indixe(i)+1) - hdltea(indixe(i)))
                  hdlow(i) = hdlowa(indixe(i)) + tdef(i)
     &            *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))
               else
                  hdlte(i) = tiny
                  hdlow(i) = tiny
               endif
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
c  old (incorrect) way:
c               hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2._DKIND)
c               hdlow1 = max(hdlow(i), tiny)
c               edot(i) = edot(i) - HDI(i,j,k)*
c     .                     (hdlte1/(1._DKIND + hdlte1/hdlow1)/(2._DKIND*dom))
c  new (correct) way: (april 4, 2007)
               hdlte1 = hdlte(i)/(HI(i,j,k)*dom)
               hdlow1 = max(hdlow(i), tiny)
               edot(i) = edot(i) - HDI(i,j,k)*
     &              (hdlte(i)/(1._DKIND + hdlte1/hdlow1)) /
     &              (3._DKIND*dom)
            end if
         enddo

      endif ! ihdcr
      endif

!     Iteration mask for metal-rich cells
      if (imetal .eq. 1) then
         do i = is+1, ie+1
            if (metallicity(i) .ge. min_metallicity) then
               itmask_metal(i) = itmask(i)
            else
               itmask_metal(i) = MASK_FALSE
            endif
         enddo
      else
         do i = is+1, ie+1
            itmask_metal(i) = MASK_FALSE
         enddo
      endif

!     Compute grain size increment

      if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then

         call calc_grain_size_increment_1d(
     &             immulti, imabund, idspecies, igrgr, itmask_metal
     &           , in, jn, kn, is, ie, j, k, dom, d
     &           , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &           , AC, SiO2D, MgO, FeS, Al2O3
     &           , reforg, volorg, H2Oice
     &           , metal, metal_loc
     &           , metal_C13, metal_C20, metal_C25, metal_C30
     &           , metal_F13, metal_F15, metal_F50, metal_F80
     &           , metal_P170, metal_P200, metal_Y19
     &           , SN0_N
     &           , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &           , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &           , SN0_fFeS, SN0_fAl2O3
     &           , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &           , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &           , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &           , SN0_r0FeS, SN0_r0Al2O3
     &           , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &           , gr_N, gr_Size, gr_dT, gr_Td
     &           , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &           , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &           , SN0_kpFeS, SN0_kpAl2O3
     &           , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &           , sgSiM, sgFeM, sgMg2SiO4, sgMgSiO3, sgFe3O4, sgAC
     &           , sgSiO2D, sgMgO, sgFeS, sgAl2O3
     &           , sgreforg, sgvolorg, sgH2Oice, sgtot
     &           , alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC
     &           , alSiO2D, alMgO, alFeS, alAl2O3
     &           , alreforg, alvolorg, alH2Oice, altot
     &         )

      endif

!     Calculate dust to gas ratio AND interstellar radiation field
!     -> an earlier version of this logic would store values @ indices
!        where `itmask_metal(i) .ne. MASK_FALSE`
!     -> this was undesirable, b/c these quantities are required for
!        photo-electric heating, which can occur when
!        `itmask_metal(i) .eq. MASK_FALSE` (we can revisit this choice
!        later). Moreover, in most cases, these calculations will be
!        faster when there is no branching

      if ((anydust .ne. MASK_FALSE) .or. (igammah .gt. 0)) then
         if (idustfield .gt. 0) then
            do i = is+1, ie+1
               ! REMINDER: use of `itmask` over `itmask_metal` is
               !   currently required by Photo-electric heating
               if (itmask(i) .ne. MASK_FALSE) then
                  ! it may be faster to remove this branching
                  dust2gas(i) = dust(i,j,k) / d(i,j,k)
               endif
            enddo
         else
            do i = is+1, ie+1
               dust2gas(i) = fgr * metallicity(i)
            enddo
         endif
      endif

      if ((anydust .ne. MASK_FALSE) .or. (igammah .gt. 1)) then
         if (iisrffield .gt. 0) then
            do i = is+1, ie+1
               myisrf(i) = isrf_habing(i,j,k)
            enddo
         else
            do i = is+1, ie+1
               myisrf(i) = isrf
            enddo
         endif
      endif

      ! compute dust temperature and cooling due to dust
      if (anydust .ne. MASK_FALSE) then

         call calc_all_tdust_gasgr_1d_g(in, jn, kn, nratec,
     &                idustfield, is, ie, j, k, fgr, gamma_isrfa,
     &                comp2, gasgra, indixe, tdef, tgas, tdust,
     &                metallicity, dust2gas, mynh, gasgr_tdust,
     &                priGridRank, priGridDim, priDataSize,
     &                metGridRank, metGridDim, metDataSize,
     &                itmask_metal, LH2_N, LH2_Size, LHD_N, LHD_Size,
     &                LCI_N, LCI_Size, LCII_N, LCII_Size, LOI_N,
     &                LOI_Size, LCO_N, LCO_Size, LOH_N, LOH_Size,
     &                LH2O_N, LH2O_Size, alphap_N, alphap_Size,
     &                idspecies, itdmulti, SN0_N, gr_N, gr_Size, gr_dT,
     &                gr_Td, tSiM, tFeM, tMg2SiO4, tMgSiO3, tFe3O4,
     &                tAC, tSiO2D, tMgO, tFeS, tAl2O3, treforg,
     &                tvolorg, tH2Oice, gasgr2a, gamma_isrf2a,
     &                coolunit, gasgr, myisrf, sgSiM, sgFeM, sgMg2SiO4,
     &                sgMgSiO3, sgFe3O4, sgAC, sgSiO2D, sgMgO, sgFeS,
     &                sgAl2O3, sgreforg, sgvolorg, sgH2Oice, sgtot,
     &                alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC,
     &                alSiO2D, alMgO, alFeS, alAl2O3, alreforg,
     &                alvolorg, alH2Oice, altot, kpSiM, kpFeM,
     &                kpMg2SiO4, kpMgSiO3, kpFe3O4, kpAC, kpSiO2D,
     &                kpMgO, kpFeS, kpAl2O3, kpreforg, kpvolorg,
     &                kpH2Oice, kptot, gasSiM, gasFeM, gasMg2SiO4,
     &                gasMgSiO3, gasFe3O4, gasAC, gasSiO2D, gasMgO,
     &                gasFeS, gasAl2O3, gasreforg, gasvolorg,
     &                gasH2Oice, gasgr_tSiM, gasgr_tFeM,
     &                gasgr_tMg2SiO4, gasgr_tMgSiO3, gasgr_tFe3O4,
     &                gasgr_tAC, gasgr_tSiO2D, gasgr_tMgO, gasgr_tFeS,
     &                gasgr_tAl2O3, gasgr_treforg, gasgr_tvolorg,
     &                gasgr_tH2Oice, mygisrf, fv2k, fac, gisrfSiM,
     &                gisrfFeM, gisrfMg2SiO4, gisrfMgSiO3, gisrfFe3O4,
     &                gisrfAC, gisrfSiO2D, gisrfMgO, gisrfFeS,
     &                gisrfAl2O3, gisrfreforg, gisrfvolorg,
     &                gisrfH2Oice, anydust
     &              )

         endif

!     Calculate dust cooling rate
         if (anydust .ne. MASK_FALSE) then


         do i = is+1, ie+1
            if ( itmask_metal(i) .ne. MASK_FALSE ) then

            if (idspecies .eq. 0) then

               Ldst(i) =         - 
     &              gasgr(i) * (tgas(i) - tdust(i)) * 
     &              dust2gas(i) * rhoH(i) * rhoH(i)

            else ! idspecies

               if (itdmulti .eq. 0) then
                   Ldst(i) = - gasgr(i) * (tgas(i) - tdust(i))
     &                * d(i,j,k) * rhoH(i)
               else

               if (idspecies .gt. 0) then
                  Ldst(i) = - (
     &                gasMgSiO3  (i) * (tgas(i) - tMgSiO3  (i))
     &              + gasAC      (i) * (tgas(i) - tAC      (i))
     &              ) * d(i,j,k) * rhoH(i)
               endif

               if (idspecies .gt. 1) then
                  Ldst(i) = Ldst(i) - (
     &                gasSiM     (i) * (tgas(i) - tSiM     (i))
     &              + gasFeM     (i) * (tgas(i) - tFeM     (i))
     &              + gasMg2SiO4 (i) * (tgas(i) - tMg2SiO4 (i))
     &              + gasFe3O4   (i) * (tgas(i) - tFe3O4   (i))
     &              + gasSiO2D   (i) * (tgas(i) - tSiO2D   (i))
     &              + gasMgO     (i) * (tgas(i) - tMgO     (i))
     &              + gasFeS     (i) * (tgas(i) - tFeS     (i))
     &              + gasAl2O3   (i) * (tgas(i) - tAl2O3   (i))
     &              ) * d(i,j,k) * rhoH(i)
               endif

               if (idspecies .gt. 2) then
                  Ldst(i) = Ldst(i) - (
     &                gasreforg  (i) * (tgas(i) - treforg  (i))
     &              + gasvolorg  (i) * (tgas(i) - tvolorg  (i))
     &              + gasH2Oice  (i) * (tgas(i) - tH2Oice  (i))
     &              ) * d(i,j,k) * rhoH(i)
               endif
               endif


            endif ! idspecies

            edot(i) = edot(i) + Ldst(i)

            endif
         enddo

      endif

!     Compute continuum opacity

      if ( ipcont .eq. 1 ) then

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then

            !! primordial continuum opacity !!
            call interpolate_2D_g(
     &        logrho(i), logT(i), int(alphap_N,8), alphap_D, alphap_dD,
     &        alphap_T, alphap_dT, int(alphap_Size,8),
     &        alphap_Data, log_a)
              alpha(i) = 1.d1**log_a
         endif
      enddo

      else
         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then
               alpha(i) = 0.
         endif
      enddo

      endif

      ! Add dust opacity.
      ! if (idspecies .eq. 0), dust opacity is overestimated at Td > 50 K
      ! We better not include dust opacity.
      if ((anydust .ne. MASK_FALSE).and.(idspecies .gt. 0)) then
         do i = is+1, ie+1
            if ( itmask_metal(i) .ne. MASK_FALSE ) then

            if (itdmulti .eq. 0) then

               alphad(i) = kptot(i)

            else

               if (idspecies .gt. 0) then
               alphad(i) = kpMgSiO3  (i)
     &                   + kpAC      (i)
               endif
               if (idspecies .gt. 1) then
               alphad(i) = alphad(i)
     &                   + kpSiM     (i)
     &                   + kpFeM     (i)
     &                   + kpMg2SiO4 (i)
     &                   + kpFe3O4   (i)
     &                   + kpSiO2D   (i)
     &                   + kpMgO     (i)
     &                   + kpFeS     (i)
     &                   + kpAl2O3   (i)
               endif
               if (idspecies .gt. 2) then
               alphad(i) = alphad(i)
     &                   + kpreforg  (i)
     &                   + kpvolorg  (i)
     &                   + kpH2Oice  (i)
               endif
            endif

            alpha(i) = alpha(i) + alphad(i) * d(i,j,k)*dom*mh

            endif
         enddo

      endif !! anydust

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
         tau_con(i) = alpha(i) * lshield_con(i)
         endif
      enddo

!     --- Compute (external) radiative heating terms ---
!     Photoionization heating

      if (ispecies .gt. 0) then

      if (iradshield == 0) then ! no shielding
        do i = is+1, ie+1
           if ( itmask(i) .ne. MASK_FALSE ) then
              edot(i) = edot(i) + real(ipiht, DKIND)*(
     &               piHI  *HI  (i,j,k) ! pi of HI
     &             + piHeI *HeI (i,j,k)*0.25_DKIND ! pi of HeI
     &             + piHeII*HeII(i,j,k)*0.25_DKIND ! pi of HeII
     &             )/dom
           end if
        enddo

      else if (iradshield == 1) then
!
!     approximate self shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     to shield HI, while leaving HeI and HeII optically thin
!

        do i = is+1, ie+1
           if (itmask(i) .ne. MASK_FALSE) then
              if (k24 .lt. tiny8) then
                 fSShHI = 1._DKIND
              else
                 nSSh = 6.73e-3_DKIND *
     &            (avgsighi /2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &            (tgas(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &            (k24/tbase1/1.0e-12_DKIND)**(2._DKIND/3._DKIND)
                 nratio = (HI(i,j,k) + HII(i,j,k))*dom/nSSh
                 fSShHI =
     &            0.98_DKIND*(1._DKIND+
     &            nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &            0.02_DKIND*(1._DKIND+
     &            nratio)**(-0.84_DKIND)
              endif

             edot(i) = edot(i) + real(ipiht,DKIND)*(
     &              piHI  *HI  (i,j,k)* fSShHI
     &            + piHeI * HeI(i,j,k)*0.25_DKIND
     &            + piHeII*HeII(i,j,k)*0.25_DKIND
     &             )/dom
           endif
        enddo

      else if (iradshield == 2)  then
!
!     Better self-shielding in HI using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     approximate self shielding in HeI and HeII
!

        do i = is+1, ie+1
           if ( itmask(i) .ne. MASK_FALSE ) then
!
!            HI self shielding ratio
!
              if (k24 .lt. tiny8) then
                 fSShHI = 1._DKIND
              else
                 nSSh = 6.73e-3_DKIND *
     &            (avgsighi/2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &            (tgas(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &            (k24/tbase1/1.0e-12_DKIND)**(2._DKIND/3._DKIND)
                 nratio = (HI(i,j,k) + HII(i,j,k))*dom/nSSh
                 fSShHI =
     &            0.98_DKIND*(1._DKIND+
     &             nratio**(1.64_DKIND))**(-2.28_DKIND)+
     &            0.02_DKIND*(1._DKIND+
     &             nratio)**(-0.84_DKIND)
              endif
!
!            HeI self shielding ratio
!
              if (k26 .lt. tiny8) then
                 fSShHeI = 1._DKIND
              else
                 nSSh_he = 6.73e-3_DKIND *
     &            (avgsighei/ 2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &            (tgas(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &            (k26/tbase1/1.0e-12_DKIND)**(2._DKIND/3._DKIND)
                 nratio_he = 0.25_DKIND*
     &            (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))*dom/nSSh_he
                 fSShHeI =
     &            0.98_DKIND*(1._DKIND+
     &             nratio_he**(1.64_DKIND))**(-2.28_DKIND)+
     &            0.02_DKIND*(1._DKIND+
     &             nratio_he)**(-0.84_DKIND)
              endif

             edot(i) = edot(i) + real(ipiht, DKIND)*(
     &              piHi * HI(i,j,k)* fSShHI
     &            + piHeI * HeI(i,j,k)*0.25_DKIND* fSShHeI
     &            + piHeII*HeII(i,j,k)*0.25_DKIND
     &              )/dom
           endif
        enddo

      else if (iradshield == 3) then
!
!     shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     in HI and HeI, but ignoring HeII heating entirely
!

        do i = is+1, ie+1
           if ( itmask(i) .ne. MASK_FALSE ) then
!
!            HI self shielding ratio
!
              if (k24 .lt. tiny8) then
                 fSShHI = 1._DKIND
              else
                 nSSh = 6.73e-3_DKIND *
     &            (avgsighi /2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &            (tgas(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &            (k24/tbase1/1.0e-12_DKIND)**(2._DKIND/3._DKIND)
                 nratio = (HI(i,j,k) + HII(i,j,k))*dom/nSSh
                 fSShHI =
     &            0.98_DKIND*(1._DKIND+
     &             nratio**(1.64_DKIND))**(-2.28_DKIND)+
     &            0.02_DKIND*(1._DKIND+
     &             nratio)**(-0.84_DKIND)
              endif
!
!            HeI self shielding ratio
!
              if (k26 .lt. tiny8) then
                 fSShHeI = 1._DKIND
              else
                 nSSh_he = 6.73e-3_DKIND *
     &            (avgsighei /2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &            (tgas(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &            (k26/tbase1/1.0e-12_DKIND)**(2._DKIND/3._DKIND)
                 nratio_he = 0.25_DKIND*
     &            (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))*dom/nSSh_he
                 fSShHeI =
     &            0.98_DKIND*(1._DKIND+
     &             nratio_he**(1.64_DKIND))**(-2.28_DKIND)+
     &            0.02_DKIND*(1._DKIND+
     &             nratio_he)**(-0.84_DKIND)
              endif

             edot(i) = edot(i) + real(ipiht, DKIND)*(
     &              piHi * HI (i,j,k)* fSShHI
     &           + piHeI * HeI(i,j,k)* fSShHeI
     &           )/dom
!
!          Ignoring HeII heating (HeII heating rate -> 0)
!
           endif
        enddo

      endif

      endif

!     --- Cloudy primordial cooling and heating ---

      if (ispecies .eq. 0) then

         iZscale = 0
         mycmbTfloor = 0
         call cool1d_cloudy_g(d, rhoH, metallicity,
     &        in, jn, kn, is, ie, j, k,
     &        logtem, edot, comp2, dom, zr,
     &        mycmbTfloor, iClHeat, iZscale,
     &        priGridRank, priGridDim,
     &        priPar1, priPar2, priPar3,
     &        priDataSize, priCooling, priHeating,
     &        itmask)

!     Calculate electron density from mean molecular weight

         do i = is+1, ie+1
            if ( itmask(i) .ne. MASK_FALSE ) then

               myde(i) = 1 - mmw(i) * (3.0_DKIND * fh + 1.0_DKIND) /
     &              4.0_DKIND
               if (imetal .eq. 1) then
                  myde(i) = myde(i) - mmw(i) * metal(i,j,k) /
     &                 (d(i,j,k) * mu_metal)
               endif
               myde(i) = d(i,j,k) * myde(i) / mmw(i)
               myde(i) = max(myde(i), 0._DKIND)

            end if
         enddo

      endif

!     Photo-electric heating by UV-irradiated dust

      if (igammah .eq. 1) then

          do i = is + 1, ie + 1
             if (itmask(i) .ne. MASK_FALSE) then
                if ( tgas(i) > 2.d4 ) then
                   gammaha_eff(i) = 0._DKIND
                else
                   gammaha_eff(i) = gammaha
                endif
            endif
         enddo

!     Use eqn. 1 of Wolfire et al. (1995)
      else if (igammah .eq. 2) then

          do i = is + 1, ie + 1
             if (itmask(i) .ne. MASK_FALSE) then
                if ( tgas(i) > 2.d4 ) then
                   gammaha_eff(i) = 0._DKIND
                else
!                  Assume constant epsilon = 0.05.
                   gammaha_eff(i) = gammaha * 0.05_DKIND * myisrf(i)
                endif
            endif
         enddo

!     Full calculation of epsilon (eqn. 2 of Wolfire 1995)
      else if (igammah .eq. 3) then

          do i = is + 1, ie + 1
             if (itmask(i) .ne. MASK_FALSE) then
                pe_X = myisrf(i) * dom_inv * sqrt(tgas(i)) / myde(i)
                pe_eps =
     &               (4.9d-2 /
     &                (1._DKIND + (pe_X / 1925._DKIND)**0.73_DKIND)) +
     &               ((3.7d-2 * (tgas(i) / 1.d4)**0.7_DKIND) /
     &                (1._DKIND + (pe_X / 5000._DKIND)))
                gammaha_eff(i) = gammaha * pe_eps * myisrf(i)
             endif
          enddo

      endif

      if (igammah .gt. 0) then
          do i = is + 1, ie + 1
             if (itmask(i) .ne. MASK_FALSE) then
                edot(i) = edot(i) + gammaha_eff(i) * rhoH(i) *
     &               dom_inv * dust2gas(i) / fgr
            endif
         enddo
      endif

!     Electron recombination onto dust grains (eqn. 9 of Wolfire 1995)

      if ((idustall .gt. 0) .or. (idustrec .gt. 0)) then

          do i = is + 1, ie + 1
             if (itmask(i) .ne. MASK_FALSE) then
                regr(i) = regra(indixe(i)) + tdef(i)
     &               *(regra(indixe(i)+1) -regra(indixe(i)))
             endif
          enddo

          do i = is + 1, ie + 1
             if (itmask(i) .ne. MASK_FALSE) then
                grbeta = 0.74_DKIND / tgas(i)**0.068_DKIND
                edot(i) = edot(i) -
     &               regr(i) * (myisrf(i)*dom_inv / myde(i))**grbeta *
     &               myde(i) * rhoH(i) * dust2gas(i) / fgr
             endif
          enddo

      endif

!     Compton cooling or heating and X-ray compton heating

      do i = is + 1, ie + 1
         if (itmask(i) .ne. MASK_FALSE) then

            edot(i) = edot(i)

!                  Compton cooling or heating

     &           - comp1      * (tgas(i) - comp2)     * myde(i)*dom_inv

!                  X-ray compton heating

     &           - comp_xraya * (tgas(i) - comp_temp) * myde(i)*dom_inv

         endif
      enddo
 
!     Photoheating from radiative transfer

      if (iradtrans .eq. 1) then
          do i = is + 1, ie + 1
            if (itmask(i) .ne. MASK_FALSE) then
              edot(i) = edot(i) + real(ipiht, DKIND) * photogamma(i,j,k)
     &                          / coolunit * HI(i,j,k) / dom

              if (edot(i) .ne. edot(i)) then
#ifdef _OPENMP
!$omp critical
#endif
                  write(6,*) 'NaN in edot[2]: ', i,j,k, edot(i),
     &                photogamma(i,j,k), HI(i,j,k), de(i,j,k), d(i,j,k),
     &                e(i,j,k), p2d(i), tgas(i), dom, urho, aye, mh
#ifdef _OPENMP
!$omp end critical
#endif
              endif

            endif
          enddo
        endif

!     --- Cloudy metal cooling and heating ---

      if (imcool .eq. 1) then

!     Determine if the temperature is above the threshold to do tabulated cooling.
         do i = is+1, ie+1
            itmask_tab(i) = itmask_metal(i)
            if ( itmask_tab(i) .ne. MASK_FALSE ) then
               if (( tmcool .gt. 0.0d0 ) .and.
     &             ( tgas(i) .lt. tmcool )) then
                  itmask_tab(i) = MASK_FALSE
               endif
            endif
         enddo

         if (clnew .eq. 1) then

            iZscale = 1
            call cool1d_cloudy_g(d, rhoH, metallicity,
     &           in, jn, kn, is, ie, j, k,
     &           logtem, edot, comp2, dom, zr,
     &           icmbTfloor, iClHeat, iZscale,
     &           metGridRank, metGridDim,
     &           metPar1, metPar2, metPar3,
     &           metDataSize, metCooling, metHeating,
     &           itmask_tab)

         else

            call cool1D_cloudy_old_tables_g(
     &           d, de, rhoH, metallicity,
     &           in, jn, kn, is, ie, j, k,
     &           logtem, edot, comp2, ispecies, dom, zr,
     &           icmbTfloor, iClHeat, 
     &           clEleFra, metGridRank, metGridDim,
     &           metPar1, metPar2, metPar3, metPar4, metPar5,
     &           metDataSize, metCooling, metHeating, 
     &           itmask_tab)

         endif

      if (imchem .eq. 1) then

!        --- C/O fine-structure, metal molecular rotational cooling for low temperatures ---
      
!        C/O fine-structure cooling
         do i = is+1, ie+1
         if ( itmask_metal(i) .ne. MASK_FALSE ) then
      
!           CI
            lognhat = logCI(i) - logdvdr(i)

            call interpolate_3D_g(
     &        lognhat, logT(i), logH(i), int(LCI_N,8),
     &        LCI_D, LCI_dD, LCI_T, LCI_dT, LCI_H, LCI_dH,
     &        int(LCI_Size,8), LCI_L, log_Linv)
            L = 1.d1**(-log_Linv)
      
            if (icmbTfloor .eq. 1) then
            call interpolate_3D_g(
     &        lognhat, logTcmb(i), logH(i), int(LCI_N,8),
     &        LCI_D, LCI_dD, LCI_T, LCI_dT, LCI_H, LCI_dH,
     &        int(LCI_Size,8), LCI_L, log_Ginv)
            G = 1.d1**(-log_Ginv)
            else
            G = tiny8
            endif
      
            LCI(i) = (G - L) / dom * CI(i,j,k)/12.d0
            if (LCI(i) .ne. LCI(i)) LCI(i) = 0.d0
            edot(i) = edot(i) + LCI(i)
      
      
!           CII
            lognhat = logCII(i) - logdvdr(i)
      
            call interpolate_3D_g(
     &        lognhat, logT(i), logH(i), int(LCII_N,8),
     &        LCII_D, LCII_dD, LCII_T, LCII_dT, LCII_H, LCII_dH,
     &        int(LCII_Size,8), LCII_L, log_Linv)
            L = 1.d1**(-log_Linv)
      
            if (icmbTfloor .eq. 1) then
            call interpolate_3D_g(
     &        lognhat, logTcmb(i), logH(i), int(LCII_N,8),
     &        LCII_D, LCII_dD, LCII_T, LCII_dT, LCII_H, LCII_dH,
     &        int(LCII_Size,8), LCII_L, log_Ginv)
            G = 1.d1**(-log_Ginv)
            else
            G = tiny8
            endif
      
            LCII(i) = (G - L) / dom * CII(i,j,k)/12.d0
            if (LCII(i) .ne. LCII(i)) LCII(i) = 0.d0
            edot(i) = edot(i) + LCII(i)
      
      
!           OI
            lognhat = logOI(i) - logdvdr(i)
      
            call interpolate_3D_g(
     &        lognhat, logT(i), logH(i), int(LOI_N,8),
     &        LOI_D, LOI_dD, LOI_T, LOI_dT, LOI_H, LOI_dH,
     &        int(LOI_Size,8), LOI_L, log_Linv)
            L = 1.d1**(-log_Linv)
      
            if (icmbTfloor .eq. 1) then
            call interpolate_3D_g(
     &        lognhat, logTcmb(i), logH(i), int(LOI_N,8),
     &        LOI_D, LOI_dD, LOI_T, LOI_dT, LOI_H, LOI_dH,
     &        int(LOI_Size,8), LOI_L, log_Ginv)
            G = 1.d1**(-log_Ginv)
            else
            G = tiny8
            endif
      
            LOI(i) = (G - L) / dom * OI(i,j,k)/16.d0
            if (LOI(i) .ne. LOI(i)) LOI(i) = 0.d0
            edot(i) = edot(i) + LOI(i)
      
      
!           Metal molecules rotational cooling

!           CO
            lognhat = logCO(i) - logdvdr(i)
      
            call interpolate_3D_g(
     &        lognhat, logT(i), logH2(i), int(LCO_N,8),
     &        LCO_D, LCO_dD, LCO_T, LCO_dT, LCO_H, LCO_dH,
     &        int(LCO_Size,8), LCO_L, log_Linv)
            L = 1.d1**(-log_Linv)
      
            if (icmbTfloor .eq. 1) then
            call interpolate_3D_g(
     &        lognhat, logTcmb(i), logH2(i), int(LCO_N,8),
     &        LCO_D, LCO_dD, LCO_T, LCO_dT, LCO_H, LCO_dH,
     &        int(LCO_Size,8), LCO_L, log_Ginv)
            G = 1.d1**(-log_Ginv)
            else
            G = tiny8
            endif
      
            LCO(i) = (G - L) / dom * CO(i,j,k)/28.d0
            if (LCO(i) .ne. LCO(i)) LCO(i) = 0.d0
            edot(i) = edot(i) + LCO(i)
      

!           OH
            lognhat = logOH(i) - logdvdr(i)
      
            call interpolate_3D_g(
     &        lognhat, logT(i), logH2(i), int(LOH_N,8),
     &        LOH_D, LOH_dD, LOH_T, LOH_dT, LOH_H, LOH_dH,
     &        int(LOH_Size,8), LOH_L, log_Linv)
            L = 1.d1**(-log_Linv)
      
            if (icmbTfloor .eq. 1) then
            call interpolate_3D_g(
     &        lognhat, logTcmb(i), logH2(i), int(LOH_N,8),
     &        LOH_D, LOH_dD, LOH_T, LOH_dT, LOH_H, LOH_dH,
     &        int(LOH_Size,8), LOH_L, log_Ginv)
            G = 1.d1**(-log_Ginv)
            else
            G = tiny8
            endif
      
            LOH(i) = (G - L) / dom * OH(i,j,k)/17.d0
            if (LOH(i) .ne. LOH(i)) LOH(i) = 0.d0
            edot(i) = edot(i) + LOH(i)


!           H2O
            lognhat = logH2O(i) - logdvdr(i)
      
            call interpolate_3D_g(
     &        lognhat, logT(i), logH2(i), int(LH2O_N,8),
     &        LH2O_D, LH2O_dD, LH2O_T, LH2O_dT, LH2O_H, LH2O_dH,
     &        int(LH2O_Size,8), LH2O_L, log_Linv)
            L = 1.d1**(-log_Linv)
      
            if (icmbTfloor .eq. 1) then
            call interpolate_3D_g(
     &        lognhat, logTcmb(i), logH2(i), int(LH2O_N,8),
     &        LH2O_D, LH2O_dD, LH2O_T, LH2O_dT, LH2O_H, LH2O_dH,
     &        int(LH2O_Size,8), LH2O_L, log_Ginv)
            G = 1.d1**(-log_Ginv)
            else
            G = tiny8
            endif
      
            LH2O(i) = (G - L) / dom * H2O(i,j,k)/18.d0
            if (LH2O(i) .ne. LH2O(i)) LH2O(i) = 0.d0
            edot(i) = edot(i) + LH2O(i)
      
         endif
         enddo

      endif ! imchem

      endif ! imcool

!     Add user-provided volumetric and/or specific heating terms

      if (iVheat .eq. 1) then

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            edot(i) = edot(i) + Vheat(i,j,k) / coolunit / dom**2
         end if
      enddo

      endif

      if (iMheat .eq. 1) then

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            edot(i) = edot(i) + Mheat(i,j,k) * d(i,j,k) * mh
     &          / coolunit / dom
         end if
      enddo

      endif

!     Continuum opacity

      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
            if ( tau_con(i) .gt. 1.d0 ) then
               if ( tau_con(i) .lt. 1.d2 ) then
                  edot(i) = edot(i) * tau_con(i)**(-2.d0)
               else
                  edot(i) = 0.d0
               endif
            endif
         end if
      enddo

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
         tgasold(i) = tgas(i)
         end if
      enddo

      return
      end
