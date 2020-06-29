#include "phys_const.def"

!!#define CHEM_DEBUG
!=======================================================================
!/////////////////////  SUBROUTINE COOL1D_MULTI_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_multi_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool,
     &                idust, idustall, idustfield,
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
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat,
     &                iisrffield, isrf_habing,
     &                itmask
#ifdef GRACKLE_MD
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
     &              , impop3
     &              , metal_loc, metal_C30, metal_F13
     &              , C30_XC , C30_XO, C30_XMg, C30_XAl
     &              , C30_XSi, C30_XS, C30_XFe
     &              , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &              , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     &              , C30_fFeS, C30_fAl2O3
     &              , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4, C30_r0MgSiO3
     &              , C30_r0Fe3O4, C30_r0AC, C30_r0SiO2D, C30_r0MgO
     &              , C30_r0FeS, C30_r0Al2O3
     &              , C30_a0SiM, C30_a0FeM, C30_a0Mg2SiO4, C30_a0MgSiO3
     &              , C30_a0Fe3O4, C30_a0AC, C30_a0SiO2D, C30_a0MgO
     &              , C30_a0FeS, C30_a0Al2O3
     &              , C30_v0SiM, C30_v0FeM, C30_v0Mg2SiO4, C30_v0MgSiO3
     &              , C30_v0Fe3O4, C30_v0AC, C30_v0SiO2D, C30_v0MgO
     &              , C30_v0FeS, C30_v0Al2O3
     &              , C30_N, C30_Size
     &              , C30_D, C30_T, C30_dD, C30_dT
     &              , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     &              , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     &              , C30_RFeS, C30_RAl2O3
     &              , C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     &              , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     &              , C30_dRFeS, C30_dRAl2O3
     &              , OC30_SiM, OC30_FeM, OC30_Mg2SiO4, OC30_MgSiO3
     &              , OC30_Fe3O4, OC30_AC, OC30_SiO2D, OC30_MgO
     &              , OC30_FeS, OC30_Al2O3
     &              , LC30_SiM, LC30_FeM, LC30_Mg2SiO4, LC30_MgSiO3
     &              , LC30_Fe3O4, LC30_AC, LC30_SiO2D, LC30_MgO
     &              , LC30_FeS, LC30_Al2O3
     &              , F13_XC , F13_XO, F13_XMg, F13_XAl
     &              , F13_XSi, F13_XS, F13_XFe
     &              , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &              , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     &              , F13_fFeS, F13_fAl2O3
     &              , F13_r0SiM, F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     &              , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0MgO
     &              , F13_r0FeS, F13_r0Al2O3
     &              , F13_a0SiM, F13_a0FeM, F13_a0Mg2SiO4, F13_a0MgSiO3
     &              , F13_a0Fe3O4, F13_a0AC, F13_a0SiO2D, F13_a0MgO
     &              , F13_a0FeS, F13_a0Al2O3
     &              , F13_v0SiM, F13_v0FeM, F13_v0Mg2SiO4, F13_v0MgSiO3
     &              , F13_v0Fe3O4, F13_v0AC, F13_v0SiO2D, F13_v0MgO
     &              , F13_v0FeS, F13_v0Al2O3
     &              , F13_N, F13_Size
     &              , F13_D, F13_T, F13_dD, F13_dT
     &              , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     &              , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     &              , F13_RFeS, F13_RAl2O3
     &              , F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     &              , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     &              , F13_dRFeS, F13_dRAl2O3
     &              , OF13_SiM, OF13_FeM, OF13_Mg2SiO4, OF13_MgSiO3
     &              , OF13_Fe3O4, OF13_AC, OF13_SiO2D, OF13_MgO
     &              , OF13_FeS, OF13_Al2O3
     &              , LF13_SiM, LF13_FeM, LF13_Mg2SiO4, LF13_MgSiO3
     &              , LF13_Fe3O4, LF13_AC, LF13_SiO2D, LF13_MgO
     &              , LF13_FeS, LF13_Al2O3
#endif
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
     &        idust, idustall, idustfield,
     &        imetal, igammah, ih2optical, iciecool, clnew,
     &        iVheat, iMheat, iradtrans, iradshield,
     &        iisrffield
#ifdef GRACKLE_MD
     &       ,imchem, igrgr
#endif

      real*8  aye, temstart, temend, z_solar, fgr,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma, fh, clEleFra
      R_PREC  d(in,jn,kn),    e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn), dust(in,jn,kn),
     &        Vheat(in,jn,kn), Mheat(in,jn,kn),
     &        photogamma(in,jn,kn), isrf_habing(in,jn,kn)
#ifdef GRACKLE_MD
     &      , DM(in,jn,kn)   , HDII(in,jn,kn)   , HeHII(in,jn,kn)
     &      , CI(in,jn,kn)   , CII(in,jn,kn)    , CO(in,jn,kn)
     &      , CO2(in,jn,kn)  , OI(in,jn,kn)     , OH(in,jn,kn)
     &      , H2O(in,jn,kn)  , O2(in,jn,kn)     , SiI(in,jn,kn)
     &      , SiOI(in,jn,kn) , SiO2I(in,jn,kn)  , CH(in,jn,kn)
     &      , CH2(in,jn,kn)  , COII(in,jn,kn)   , OII(in,jn,kn)
     &      , OHII(in,jn,kn) , H2OII(in,jn,kn)  , H3OII(in,jn,kn)
     &      , O2II(in,jn,kn) , Mg(in,jn,kn)     , Al(in,jn,kn)
     &      , S(in,jn,kn)    , Fe(in,jn,kn)     , SiM(in,jn,kn)
     &      , FeM(in,jn,kn)  , Mg2SiO4(in,jn,kn), MgSiO3(in,jn,kn)
     &      , Fe3O4(in,jn,kn), AC(in,jn,kn)     , SiO2D(in,jn,kn)
     &      , MgO(in,jn,kn)  , FeS(in,jn,kn)    , Al2O3(in,jn,kn)
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
#endif
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

      integer ti_max
      real*8 mh, mu_metal
      parameter (mh = mass_h)      !DPC
      parameter (mu_metal = 16._DKIND)    ! approx. mean molecular weight of metals
      parameter (ti_max = 20)

!  Locals

      integer i, iter, ti, iradfield
      real*8 dom, qq, vibl, logtem0, logtem9, dlogtem, zr,
     &       comp1, comp2, 
     &       hdlte1, hdlow1, gamma2, x, fudge, fH2,
     &       gphdl1, dom_inv, tau, ciefudge,
     &       coolunit, dbase1, tbase1, xbase1,
     &       nH2, nother, nSSh, nratio, nssh_he, nratio_he,
     &       fSShHI, fSShHeI, pe_eps, pe_X, grbeta

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
#ifdef GRACKLE_MD
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
       integer grain_N(2), grain_Size
       real*8  grain_D(grain_N(1)), grain_T(grain_N(2))
     &       , grain_dD, grain_dT, Tgrain(grain_Size)
     &       , Ograin(grain_Size), Lgrain(grain_Size)
      integer impop3
      real*8  C30_XC , C30_XO, C30_XMg, C30_XAl
     &      , C30_XSi, C30_XS, C30_XFe
      real*8  C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &      , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     &      , C30_fFeS, C30_fAl2O3
      real*8  C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4, C30_r0MgSiO3
     &      , C30_r0Fe3O4, C30_r0AC, C30_r0SiO2D, C30_r0MgO
     &      , C30_r0FeS, C30_r0Al2O3
      real*8  C30_a0SiM, C30_a0FeM, C30_a0Mg2SiO4, C30_a0MgSiO3
     &      , C30_a0Fe3O4, C30_a0AC, C30_a0SiO2D, C30_a0MgO
     &      , C30_a0FeS, C30_a0Al2O3
      real*8  C30_v0SiM, C30_v0FeM, C30_v0Mg2SiO4, C30_v0MgSiO3
     &      , C30_v0Fe3O4, C30_v0AC, C30_v0SiO2D, C30_v0MgO
     &      , C30_v0FeS, C30_v0Al2O3
      integer C30_N(3), C30_Size
      real*8  C30_D(C30_N(1)), C30_T(C30_N(2)), C30_dD, C30_dT
      real*8  C30_RSiM(C30_N(3)), C30_RFeM(C30_N(3))
     &      , C30_RMg2SiO4(C30_N(3)), C30_RMgSiO3(C30_N(3))
     &      , C30_RFe3O4(C30_N(3)), C30_RAC(C30_N(3))
     &      , C30_RSiO2D(C30_N(3)), C30_RMgO(C30_N(3))
     &      , C30_RFeS(C30_N(3)), C30_RAl2O3(C30_N(3))
      real*8  C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     &      , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     &      , C30_dRFeS, C30_dRAl2O3
      real*8  OC30_SiM(C30_Size), OC30_FeM(C30_Size)
     &      , OC30_Mg2SiO4(C30_Size), OC30_MgSiO3(C30_Size)
     &      , OC30_Fe3O4(C30_Size), OC30_AC(C30_Size)
     &      , OC30_SiO2D(C30_Size), OC30_MgO(C30_Size)
     &      , OC30_FeS(C30_Size), OC30_Al2O3(C30_Size)
      real*8  LC30_SiM(C30_Size), LC30_FeM(C30_Size)
     &      , LC30_Mg2SiO4(C30_Size), LC30_MgSiO3(C30_Size)
     &      , LC30_Fe3O4(C30_Size), LC30_AC(C30_Size)
     &      , LC30_SiO2D(C30_Size), LC30_MgO(C30_Size)
     &      , LC30_FeS(C30_Size), LC30_Al2O3(C30_Size)
      real*8  F13_XC , F13_XO, F13_XMg, F13_XAl
     &      , F13_XSi, F13_XS, F13_XFe
      real*8  F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &      , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     &      , F13_fFeS, F13_fAl2O3
      real*8  F13_r0SiM, F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     &      , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0MgO
     &      , F13_r0FeS, F13_r0Al2O3
      real*8  F13_a0SiM, F13_a0FeM, F13_a0Mg2SiO4, F13_a0MgSiO3
     &      , F13_a0Fe3O4, F13_a0AC, F13_a0SiO2D, F13_a0MgO
     &      , F13_a0FeS, F13_a0Al2O3
      real*8  F13_v0SiM, F13_v0FeM, F13_v0Mg2SiO4, F13_v0MgSiO3
     &      , F13_v0Fe3O4, F13_v0AC, F13_v0SiO2D, F13_v0MgO
     &      , F13_v0FeS, F13_v0Al2O3
      integer F13_N(3), F13_Size
      real*8  F13_D(F13_N(1)), F13_T(F13_N(2)), F13_dD, F13_dT
      real*8  F13_RSiM(F13_N(3)), F13_RFeM(F13_N(3))
     &      , F13_RMg2SiO4(F13_N(3)), F13_RMgSiO3(F13_N(3))
     &      , F13_RFe3O4(F13_N(3)), F13_RAC(F13_N(3))
     &      , F13_RSiO2D(F13_N(3)), F13_RMgO(F13_N(3))
     &      , F13_RFeS(F13_N(3)), F13_RAl2O3(F13_N(3))
      real*8  F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     &      , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     &      , F13_dRFeS, F13_dRAl2O3
      real*8  OF13_SiM(F13_Size), OF13_FeM(F13_Size)
     &      , OF13_Mg2SiO4(F13_Size), OF13_MgSiO3(F13_Size)
     &      , OF13_Fe3O4(F13_Size), OF13_AC(F13_Size)
     &      , OF13_SiO2D(F13_Size), OF13_MgO(F13_Size)
     &      , OF13_FeS(F13_Size), OF13_Al2O3(F13_Size)
      real*8  LF13_SiM(F13_Size), LF13_FeM(F13_Size)
     &      , LF13_Mg2SiO4(F13_Size), LF13_MgSiO3(F13_Size)
     &      , LF13_Fe3O4(F13_Size), LF13_AC(F13_Size)
     &      , LF13_SiO2D(F13_Size), LF13_MgO(F13_Size)
     &      , LF13_FeS(F13_Size), LF13_Al2O3(F13_Size)
       real*8 logdom, logT(in), logTcmb(in), logrho(in)
     &      , logH(in), logH2I(in), logHDI(in)
     &      , logH2(in), logCI(in), logCII(in), logOI(in)
     &      , logCO(in), logOH(in), logH2O(in)
       real*8 logdvdr(in), lognhat
       real*8 log_Linv, log_Ginv, L, G
       real*8 LH2(in), LCIE(in), LHD(in)
     &      , LCI(in), LCII(in), LOI(in)
     &      , LCO(in), LOH(in), LH2O(in)
       real*8 Ldst(in), alpha(in), alphad(in)
       real*8 lshield_con(in), tau_con(in)
       real*8 log_a, log_L, log_G
       real*8 LC30(in), GC30(in), OC30(in)
       real*8 LF13(in), GF13(in), OF13(in)
!     grain growth
      real*8  FC30_SiM(in), FC30_FeM(in), FC30_Mg2SiO4(in)
     &      , FC30_MgSiO3(in), FC30_Fe3O4(in), FC30_AC(in)
     &      , FC30_SiO2D(in), FC30_MgO(in), FC30_FeS(in)
     &      , FC30_Al2O3(in)
      real*8  FF13_SiM(in), FF13_FeM(in), FF13_Mg2SiO4(in)
     &      , FF13_MgSiO3(in), FF13_Fe3O4(in), FF13_AC(in)
     &      , FF13_SiO2D(in), FF13_MgO(in), FF13_FeS(in)
     &      , FF13_Al2O3(in)
#endif
!  Iteration mask

      logical itmask(in), anydust, interp
#ifdef GRACKLE_MD
      logical itmask_line(in), itmask_tmp(in)
#endif
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
      integer iter_tgas
      real*8  tgas_err, tgas0
#endif /* NOT important */
#ifdef CHEM_DEBUG
      integer, parameter :: iref = 1
      integer, parameter :: jref = 1
      integer, parameter :: kref = 1
#endif

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set flag for dust-related options

      anydust = (idust .gt. 0) .or. (idustall .gt. 0)

!     Set flag for needing interpolation variables

      interp = (ispecies .gt. 0) .or. (idustall .gt. 0)

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

!     Initialize edot

      do i = is+1, ie+1
         if ( itmask(i) ) then
            edot(i) = 0._DKIND
         end if
      enddo

!     Compute Pressure

      do i = is+1, ie+1
         if ( itmask(i) ) then
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
               if ( itmask(i) ) then
                  rhoH(i) = fh * (d(i,j,k) - metal(i,j,k))
               endif
            enddo
         else
            do i = is+1, ie+1
               if ( itmask(i) ) then
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
            if ( itmask(i) ) then
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
               if ( itmask(i) ) then
                  mmw(i) = mmw(i) +
     &                 HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2._DKIND
                  rhoH(i) = rhoH(i) + H2I(i,j,k) + H2II(i,j,k)
               end if
            enddo
         endif

!     Include metal species

         if (imetal .eq. 1) then
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  mmw(i) = mmw(i) + metal(i,j,k)/mu_metal
               end if
            enddo
         endif

         do i = is+1, ie+1
            if ( itmask(i) ) then
               tgas(i) = max(p2d(i)*utem/mmw(i), temstart)
               mmw(i) = d(i,j,k) / mmw(i)
            end if
         enddo

!     Correct temperature for gamma from H2

         if (ispecies .gt. 1) then
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  nH2 = 0.5_DKIND*(H2I(i,j,k) + H2II(i,j,k))
                  nother = (HeI(i,j,k) + HeII(i,j,k) +
     &                 HeIII(i,j,k))/4._DKIND +
     &                 HI(i,j,k) + HII(i,j,k) + de(i,j,k)
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
                  iter_tgas = 0
                  tgas_err = huge8
                  do while ((iter_tgas .lt. 100)
     &                 .and.(tgas_err .gt. 1.d-3))
                  tgas0 = tgas(i)
#endif
                  if (nH2/nother .gt. 1.0e-3_DKIND) then
                     x = 6100._DKIND/tgas(i) ! not quite self-consistent
                     if (x .gt. 10._DKIND) then
                        gamma2 = 0.5_DKIND*5._DKIND
                     else
                        gamma2 = 0.5_DKIND*(5._DKIND + 2._DKIND*x**2 * 
     &                       exp(x)/(exp(x)-1)**2)
                     endif
                  else
                     gamma2 = 2.5_DKIND
                  endif
                  gamma2 = 1._DKIND + (nH2 + nother)/
     &                 (nH2*gamma2 + nother/(gamma-1._DKIND))
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
                  tgas(i) = max((gamma2 - 1._DKIND)*mmw(i)*e(i,j,k)*utem
     &                        , temstart)
                  tgas_err = dabs(tgas0 - tgas(i)) / tgas0
                  iter_tgas = iter_tgas + 1
                  end do
#else
                  tgas(i) = tgas(i) * (gamma2 - 1._DKIND)/
     &                 (gamma - 1._DKIND)
#endif
               end if
            enddo
         endif

      endif

!     Calculate metallicity and H number density

      if (imetal .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
               metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar
            endif
         enddo
      endif

      do i = is+1, ie+1
         if ( itmask(i) ) then
            mynh(i) = rhoH(i) * dom
         end if
      enddo

!     If this is the first time through, just set tgasold to tgas

      if (iter .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            tgasold(i) = tgas(i)
            end if
         enddo
      endif

#ifdef GRACKLE_MD

!     Compute log densities

      logdom = log10(dom)
      do i = is+1, ie+1
         if ( itmask(i) ) then
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
            logdvdr(i) = -8.79947961814d0 + 0.5d0 * logrho(i) ! km/s / cm
            lshield_con(i) = 
     &         sqrt((gamma * pi_val * kboltz * tgas(i)) /
     &         (GravConst *  mmw(i)*mh * d(i,j,k)*dom*mh))

          endif
      enddo
#endif
     
      do i = is+1, ie+1
         if ( itmask(i) ) then

!        Compute log temperature and truncate if above/below table max/min

         logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

         endif
      enddo

!     Compute interpolation indices

      if (interp) then
         do i = is+1, ie+1
            if ( itmask(i) ) then

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
         if ( itmask(i) ) then

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
         if ( itmask(i) ) then
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

#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'pri', tgas(i), edot(i)
#endif

#ifdef UNUSED
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
#endif         

         end if
      enddo

      endif

!     --- H2 cooling ---

      if (ispecies .gt. 1) then
#define USE_CHIAKI_WISE2019
#if defined(GRACKLE_MD) && defined(USE_CHIAKI_WISE2019)

      do i = is+1, ie+1
      if ( itmask(i) ) then

      if (logH2I(i) .gt. - 5.d0) then
         lognhat = logH2I(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int8(LH2_N), 
     &     LH2_D, LH2_dD, LH2_T, LH2_dT, LH2_H, LH2_dH,
     &     int8(LH2_Size), LH2_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int8(LH2_N), 
     &     LH2_D, LH2_dD, LH2_T, LH2_dT, LH2_H, LH2_dH,
     &     int8(LH2_Size), LH2_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LH2(i) = - (L - G) / dom * H2I(i,j,k)/2.d0
         if (LH2(i) .ne. LH2(i)) LH2(i) = 0.d0
         edot(i) = edot(i) + LH2(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'H2 ', logH2I(i), logdvdr(i), logT(i), logH(i)
     &     , L * coolunit
#endif
      endif

      endif
      enddo
#else  /* USE_CHIAKI_WISE2019 */

#define USE_GLOVER_ABEL2008
#ifdef USE_GLOVER_ABEL2008
         do i = is+1, ie+1
            if ( itmask(i) ) then
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
            if ( itmask(i) ) then
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
            edot(i) = edot(i) - real(ih2co, DKIND)*fudge*H2I(i,j,k)*
     &           h2lte(i)/(1._DKIND + gphdl1/galdl(i)) / (2._DKIND*dom)

            end if
         enddo
#else

#define USE_GALLI_PALLA1999
#define NO_OPTICAL_DEPTH_FUDGE

!        Use the Galli and Palla (1999) cooling rates for molecular H.

#ifdef USE_GALLI_PALLA1999

         do i = is+1, ie+1
            if ( itmask(i) ) then
            gpldl(i) = gpldla(indixe(i)) + tdef(i)
     &         *(gpldla(indixe(i)+1) - gpldla(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then

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
            edot(i) = edot(i) - real(ih2co, DKIND)*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1._DKIND + gphdl1/gpldl(i)) / (2._DKIND*dom)

            end if
         enddo

#else /* USE_GALLI_PALLA1999 */

         do i = is+1, ie+1
            if ( itmask(i) ) then
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
            if ( itmask(i) ) then
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

            edot(i) = edot(i) - real(ih2co, DKIND)*fudge*H2I(i,j,k)*(
     &           vibh(i)/(1._DKIND+vibh(i)/max(   vibl     ,tiny)) +
     &           roth(i)/(1._DKIND+roth(i)/max(qq*rotl(i),tiny))     
     &           )/2._DKIND/dom
            end if
         enddo

#endif /* USE_GALLI_PALLA1999 */
#endif /* USE_GLOVER_ABEL2008 */

#endif /* USE_CHIAKI_WISE2019 */

c     CIE
c     cooling from H2-H2 and He-H2 collisional induced emission comes
C     with its own radiative transfer correction as discussed in
C     Ripamonti & Abel 2003
         if (iciecool.eq.1) then
            do i = is+1, ie+1
            if (itmask(i)) then
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
#define USE_YOSHIDA06
#if defined(GRACKLE_MD) && defined(USE_YOSHIDA06)
                cieY06(i) = cieY06a(indixe(i)) + tdef(i)
     &             *(cieY06a(indixe(i)+1) - cieY06a(indixe(i)))
                LCIE(i) = - cieY06(i) * (H2I(i,j,k)/2.d0)**2
                edot(i) = edot(i) + LCIE(i)
#else
                edot(i) = ciefudge*(edot(i) - 
     &                  H2I(i,j,k)*(d(i,j,k)*cieco(i)))
#endif /* USE_YOSHIDA06 */
#ifdef CHEM_DEBUG
                if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &          write(*,*) 'CIE', LCIE(i), edot(i)
#endif
              endif
            endif
            enddo
         endif

      endif

!     --- Cooling from HD ---

      if (ispecies .gt. 2) then
#define USE_CHIAKI_WISE2019
#if defined(GRACKLE_MD) && defined(USE_CHIAKI_WISE2019)
      do i = is+1, ie+1
      if ( itmask(i) ) then

      if (logHDI(i) .gt. - 9.5d0) then
         lognhat = logHDI(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int8(LHD_N), 
     &     LHD_D, LHD_dD, LHD_T, LHD_dT, LHD_H, LHD_dH,
     &     int8(LHD_Size), LHD_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int8(LHD_N), 
     &     LHD_D, LHD_dD, LHD_T, LHD_dT, LHD_H, LHD_dH,
     &     int8(LHD_Size), LHD_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LHD(i) = - (L - G) / dom * HDI(i,j,k)/3.d0
         if (LHD(i) .ne. LHD(i)) LHD(i) = 0.d0
         edot(i) = edot(i) + LHD(i) !!!! CAUTION !!!!
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'HD ', LHD(i), edot(i)
#endif
      endif

      endif
      enddo
#else
         do i = is+1, ie+1
            if ( itmask(i) ) then
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
            if ( itmask(i) ) then
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
#endif /* USE_CHIAKI_WISE2019 */ 
      endif

!     Calculate dust to gas ratio

      if (anydust .or. (igammah .gt. 0)) then
         if (idustfield .gt. 0) then
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  dust2gas(i) = dust(i,j,k) / d(i,j,k)
               endif
            enddo
         else
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  dust2gas(i) = fgr * metallicity(i)
               endif
            enddo
         endif
      endif

!     Calculate interstellar radiation field

      if (anydust .or. (igammah .gt. 1)) then
         if (iisrffield .gt. 0) then
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  myisrf(i) = isrf_habing(i,j,k)
               endif
            enddo
         else
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  myisrf(i) = isrf
               endif
            enddo
         endif
      endif

#ifdef GRACKLE_MD

!     Compute grain size increment

      if (anydust) then
         if (impop3 .gt. 0) then

            call calc_grain_size_increment_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , metal, metal_C30, metal_F13
     &              , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &              , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     &              , C30_fFeS, C30_fAl2O3
     &              , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4, C30_r0MgSiO3
     &              , C30_r0Fe3O4, C30_r0AC, C30_r0SiO2D, C30_r0MgO
     &              , C30_r0FeS, C30_r0Al2O3
     &              , C30_a0SiM, C30_a0FeM, C30_a0Mg2SiO4, C30_a0MgSiO3
     &              , C30_a0Fe3O4, C30_a0AC, C30_a0SiO2D, C30_a0MgO
     &              , C30_a0FeS, C30_a0Al2O3
     &              , C30_v0SiM, C30_v0FeM, C30_v0Mg2SiO4, C30_v0MgSiO3
     &              , C30_v0Fe3O4, C30_v0AC, C30_v0SiO2D, C30_v0MgO
     &              , C30_v0FeS, C30_v0Al2O3
     &              , C30_N
     &              , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     &              , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     &              , C30_RFeS, C30_RAl2O3
     &              , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &              , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     &              , F13_fFeS, F13_fAl2O3
     &              , F13_r0SiM, F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     &              , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0MgO
     &              , F13_r0FeS, F13_r0Al2O3
     &              , F13_a0SiM, F13_a0FeM, F13_a0Mg2SiO4, F13_a0MgSiO3
     &              , F13_a0Fe3O4, F13_a0AC, F13_a0SiO2D, F13_a0MgO
     &              , F13_a0FeS, F13_a0Al2O3
     &              , F13_v0SiM, F13_v0FeM, F13_v0Mg2SiO4, F13_v0MgSiO3
     &              , F13_v0Fe3O4, F13_v0AC, F13_v0SiO2D, F13_v0MgO
     &              , F13_v0FeS, F13_v0Al2O3
     &              , F13_N
     &              , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     &              , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     &              , F13_RFeS, F13_RAl2O3
     &              , FC30_SiM, FC30_FeM, FC30_Mg2SiO4, FC30_MgSiO3
     &              , FC30_Fe3O4, FC30_AC, FC30_SiO2D, FC30_MgO
     &              , FC30_FeS, FC30_Al2O3
     &              , FF13_SiM, FF13_FeM, FF13_Mg2SiO4, FF13_MgSiO3
     &              , FF13_Fe3O4, FF13_AC, FF13_SiO2D, FF13_MgO
     &              , FF13_FeS, FF13_Al2O3
     &            )

         endif
      endif
#endif

!     --- Gas to grain heat transfer ---

      if (anydust) then

#ifdef GRACKLE_MD
         if ( impop3 .eq. 0 ) then
            do i = is+1, ie+1
               if ( itmask(i) ) then
          
            call interpolate_2D_g(
     &        logrho(i), logT(i), int8(grain_N), grain_D, grain_dD, 
     &        grain_T, grain_dT, int8(grain_Size), 
     &        Tgrain, tdust(i))

               endif
            enddo
         endif
#else
!     Look up gas/grain heat transfer rates

         do i = is+1, ie+1
            if ( itmask(i) ) then
               gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &              *(gasgra(indixe(i)+1) -gasgra(indixe(i)))
               gasgr_tdust(i) = fgr * gasgr(i) * coolunit / mh
            endif
         enddo

!     Compute dust temperature

         call calc_tdust_1d_g(tdust, tgas, mynh, gasgr_tdust,
     &        gamma_isrfa, myisrf, itmask, comp2, in, is, ie, j, k)
#endif
          
!     Calculate dust cooling rate

#ifdef GRACKLE_MD
         do i = is+1, ie+1
            if ( itmask(i) ) then

            call interpolate_2D_g(
     &        logrho(i), logT(i), int8(grain_N), grain_D, grain_dD, 
     &        grain_T, grain_dT, int8(grain_Size), 
     &        Lgrain, log_L)
              Ldst(i) = 1.d1**(log_L)

            if (icmbTfloor .eq. 1) then
            call interpolate_2D_g(
     &        logrho(i), logTcmb(i), int8(grain_N), grain_D, grain_dD, 
     &        grain_T, grain_dT, int8(grain_Size), 
     &        Lgrain, log_G)
              Ldst(i) = Ldst(i) - 1.d1**(log_G)
            endif

            Ldst(i) = - Ldst(i) / (dom * dom)
     &               * metal_loc(i,j,k)/(0.02d0*d(i,j,k))

            endif
         enddo

         if ( impop3 .gt. 0 ) then

            call lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FC30_SiM, FC30_FeM, FC30_Mg2SiO4, FC30_MgSiO3
     &              , FC30_Fe3O4, FC30_AC, FC30_SiO2D, FC30_MgO
     &              , FC30_FeS, FC30_Al2O3
     &              , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &              , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     &              , C30_fFeS, C30_fAl2O3
     &              , C30_N, C30_Size
     &              , C30_D, C30_T, C30_dD, C30_dT
     &              , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     &              , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     &              , C30_RFeS, C30_RAl2O3
     &              , C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     &              , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     &              , C30_dRFeS, C30_dRAl2O3
     &              , LC30_SiM, LC30_FeM, LC30_Mg2SiO4, LC30_MgSiO3
     &              , LC30_Fe3O4, LC30_AC, LC30_SiO2D, LC30_MgO
     &              , LC30_FeS, LC30_Al2O3
     &              , LC30
     &            )
            call lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FF13_SiM, FF13_FeM, FF13_Mg2SiO4, FF13_MgSiO3
     &              , FF13_Fe3O4, FF13_AC, FF13_SiO2D, FF13_MgO
     &              , FF13_FeS, FF13_Al2O3
     &              , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &              , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     &              , F13_fFeS, F13_fAl2O3
     &              , F13_N, F13_Size
     &              , F13_D, F13_T, F13_dD, F13_dT
     &              , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     &              , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     &              , F13_RFeS, F13_RAl2O3
     &              , F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     &              , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     &              , F13_dRFeS, F13_dRAl2O3
     &              , LF13_SiM, LF13_FeM, LF13_Mg2SiO4, LF13_MgSiO3
     &              , LF13_Fe3O4, LF13_AC, LF13_SiO2D, LF13_MgO
     &              , LF13_FeS, LF13_Al2O3
     &              , LF13
     &            )

            if (icmbTfloor .eq. 1) then
            call lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logTcmb
     &              , FC30_SiM, FC30_FeM, FC30_Mg2SiO4, FC30_MgSiO3
     &              , FC30_Fe3O4, FC30_AC, FC30_SiO2D, FC30_MgO
     &              , FC30_FeS, FC30_Al2O3
     &              , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &              , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     &              , C30_fFeS, C30_fAl2O3
     &              , C30_N, C30_Size
     &              , C30_D, C30_T, C30_dD, C30_dT
     &              , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     &              , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     &              , C30_RFeS, C30_RAl2O3
     &              , C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     &              , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     &              , C30_dRFeS, C30_dRAl2O3
     &              , LC30_SiM, LC30_FeM, LC30_Mg2SiO4, LC30_MgSiO3
     &              , LC30_Fe3O4, LC30_AC, LC30_SiO2D, LC30_MgO
     &              , LC30_FeS, LC30_Al2O3
     &              , GC30
     &            )
            call lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logTcmb
     &              , FF13_SiM, FF13_FeM, FF13_Mg2SiO4, FF13_MgSiO3
     &              , FF13_Fe3O4, FF13_AC, FF13_SiO2D, FF13_MgO
     &              , FF13_FeS, FF13_Al2O3
     &              , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &              , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     &              , F13_fFeS, F13_fAl2O3
     &              , F13_N, F13_Size
     &              , F13_D, F13_T, F13_dD, F13_dT
     &              , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     &              , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     &              , F13_RFeS, F13_RAl2O3
     &              , F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     &              , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     &              , F13_dRFeS, F13_dRAl2O3
     &              , LF13_SiM, LF13_FeM, LF13_Mg2SiO4, LF13_MgSiO3
     &              , LF13_Fe3O4, LF13_AC, LF13_SiO2D, LF13_MgO
     &              , LF13_FeS, LF13_Al2O3
     &              , GF13
     &            )
               do i = is+1, ie+1
                  if ( itmask(i) ) then
                     LC30(i) = LC30(i) - GC30(i)
                     LF13(i) = LF13(i) - GF13(i)
                  endif
               enddo
            endif

            do i = is+1, ie+1
               if ( itmask(i) ) then
                  Ldst(i) = Ldst(i) -
     &              ( LC30(i) * metal_C30(i,j,k)
     &              + LF13(i) * metal_F13(i,j,k)
     &              )/(0.02d0*d(i,j,k)) / (dom * dom)
               endif
            enddo

         endif !! impop3

         do i = is+1, ie+1
            if ( itmask(i) ) then
            edot(i) = edot(i) + Ldst(i)
#ifdef CHEM_DEBUG
            if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &      write(*,*) 'dst', Ldst(i), edot(i)
#endif
            endif
         enddo

#else

         do i = is+1, ie+1
            if ( itmask(i) ) then
               edot(i) = edot(i) - 
     &              gasgr(i) * (tgas(i) - tdust(i)) * 
     &              dust2gas(i) * rhoH(i) * rhoH(i)
            endif
         enddo

#endif
      endif

#ifdef GRACKLE_MD

!     Compute continuum opacity

      do i = is+1, ie+1
         if ( itmask(i) ) then

            call interpolate_2D_g(
     &        logrho(i), logT(i), int8(alphap_N), alphap_D, alphap_dD, 
     &        alphap_T, alphap_dT, int8(alphap_Size), 
     &        alphap_Data, log_a)
              alpha(i) = 1.d1**log_a

         endif
      enddo

      if (anydust) then

         do i = is+1, ie+1
            if ( itmask(i) ) then

            call interpolate_2D_g(
     &        logrho(i), logT(i), int8(grain_N), grain_D, grain_dD, 
     &        grain_T, grain_dT, int8(grain_Size), 
     &        Ograin, log_a)

            alpha(i) = alpha(i) + 1.d1**log_a
     &             *metal_loc(i,j,k)/(0.02d0*d(i,j,k))

            endif
         enddo

         if ( impop3 .gt. 0 ) then

            call lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FC30_SiM, FC30_FeM, FC30_Mg2SiO4, FC30_MgSiO3
     &              , FC30_Fe3O4, FC30_AC, FC30_SiO2D, FC30_MgO
     &              , FC30_FeS, FC30_Al2O3
     &              , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &              , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     &              , C30_fFeS, C30_fAl2O3
     &              , C30_N, C30_Size
     &              , C30_D, C30_T, C30_dD, C30_dT
     &              , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     &              , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     &              , C30_RFeS, C30_RAl2O3
     &              , C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     &              , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     &              , C30_dRFeS, C30_dRAl2O3
     &              , OC30_SiM, OC30_FeM, OC30_Mg2SiO4, OC30_MgSiO3
     &              , OC30_Fe3O4, OC30_AC, OC30_SiO2D, OC30_MgO
     &              , OC30_FeS, OC30_Al2O3
     &              , OC30
     &            )
            call lookup_dust_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FF13_SiM, FF13_FeM, FF13_Mg2SiO4, FF13_MgSiO3
     &              , FF13_Fe3O4, FF13_AC, FF13_SiO2D, FF13_MgO
     &              , FF13_FeS, FF13_Al2O3
     &              , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &              , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     &              , F13_fFeS, F13_fAl2O3
     &              , F13_N, F13_Size
     &              , F13_D, F13_T, F13_dD, F13_dT
     &              , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     &              , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     &              , F13_RFeS, F13_RAl2O3
     &              , F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     &              , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     &              , F13_dRFeS, F13_dRAl2O3
     &              , OF13_SiM, OF13_FeM, OF13_Mg2SiO4, OF13_MgSiO3
     &              , OF13_Fe3O4, OF13_AC, OF13_SiO2D, OF13_MgO
     &              , OF13_FeS, OF13_Al2O3
     &              , OF13
     &            )

            do i = is+1, ie+1
               if ( itmask(i) ) then
                  alpha(i) = alpha(i) +
     &              ( OC30(i) * metal_C30(i,j,k)
     &              + OF13(i) * metal_F13(i,j,k)
     &              )/(0.02d0*d(i,j,k))
               endif
            enddo

         endif !! impop3

      endif !! anydust

      do i = is+1, ie+1
         if ( itmask(i) ) then
         tau_con(i) = alpha(i) * lshield_con(i)
         endif
      enddo
#endif

!     --- Compute (external) radiative heating terms ---
!     Photoionization heating

      if (ispecies .gt. 0) then

      if (iradshield == 0) then ! no shielding
        do i = is+1, ie+1
           if ( itmask(i) ) then
              edot(i) = edot(i) + real(ipiht, DKIND)*(
     &               piHI  *HI  (i,j,k) ! pi of HI
     &             + piHeI *HeI (i,j,k)*0.25_DKIND ! pi of HeI
     &             + piHeII*HeII(i,j,k)*0.25_DKIND ! pi of HeII
     &             )/dom
#ifdef CHEM_DEBUG
           if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &     write(*,*) 'rs ', piHI, edot(i)
#endif
           end if
        enddo

      else if (iradshield == 1) then
!
!     approximate self shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     to shield HI, while leaving HeI and HeII optically thin
!

        do i = is+1, ie+1
           if (itmask(i)) then
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
#ifdef CHEM_DEBUG
           if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &     write(*,*) 'ss ', piHI, edot(i)
#endif
           endif
        enddo

      else if (iradshield == 2)  then
!
!     Better self-shielding in HI using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     approximate self shielding in HeI and HeII
!

        do i = is+1, ie+1
           if ( itmask(i) ) then
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
#ifdef CHEM_DEBUG
           if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &     write(*,*) 'rad', piHI, edot(i)
#endif
           endif
        enddo

      else if (iradshield == 3) then
!
!     shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     in HI and HeI, but ignoring HeII heating entirely
!

        do i = is+1, ie+1
           if ( itmask(i) ) then
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
#ifdef CHEM_DEBUG
             if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &       write(*,*) 'ss2', piHI, edot(i)
#endif
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
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'clo', ispecies, edot(i)
#endif

!     Calculate electron density from mean molecular weight

         do i = is+1, ie+1
            if ( itmask(i) ) then

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
             if (itmask(i)) then
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
             if (itmask(i)) then
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
             if (itmask(i)) then
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
             if (itmask(i)) then
                edot(i) = edot(i) + gammaha_eff(i) * rhoH(i) *
     &               dom_inv * dust2gas(i) / fgr
#ifdef CHEM_DEBUG
             if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &       write(*,*) 'gmh', gammaha_eff(i), edot(i)
#endif
            endif
         enddo
      endif

!     Electron recombination onto dust grains (eqn. 9 of Wolfire 1995)

      if (idustall .gt. 0) then

          do i = is + 1, ie + 1
             if (itmask(i)) then
                regr(i) = regra(indixe(i)) + tdef(i)
     &               *(regra(indixe(i)+1) -regra(indixe(i)))
             endif
          enddo

          do i = is + 1, ie + 1
             if (itmask(i)) then
                grbeta = 0.74_DKIND / tgas(i)**0.068_DKIND
                edot(i) = edot(i) -
     &               regr(i) * (myisrf(i)*dom_inv / myde(i))**grbeta *
     &               myde(i) * rhoH(i)
#ifdef CHEM_DEBUG
             if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &       write(*,*) 'erd', regr(i), edot(i)
#endif
             endif
          enddo

      endif

!     Compton cooling or heating and X-ray compton heating

      do i = is + 1, ie + 1
         if (itmask(i)) then

            edot(i) = edot(i)

!                  Compton cooling or heating

     &           - comp1      * (tgas(i) - comp2)     * myde(i)*dom_inv

!                  X-ray compton heating

     &           - comp_xraya * (tgas(i) - comp_temp) * myde(i)*dom_inv

#ifdef CHEM_DEBUG
             if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &       write(*,*) 'com', comp1, edot(i)
#endif
         endif
      enddo
 
!     Photoheating from radiative transfer

      if (iradtrans .eq. 1) then
          do i = is + 1, ie + 1
            if (itmask(i)) then
              edot(i) = edot(i) + real(ipiht, DKIND) * photogamma(i,j,k)
     &                          / coolunit * HI(i,j,k) / dom
#ifdef CHEM_DEBUG
              if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &        write(*,*) 'pg ', photogamma(i,j,k), edot(i)
#endif

#ifdef UNUSED
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
#endif

            endif
          enddo
        endif

!     --- Cloudy metal cooling and heating ---

      if (imcool .eq. 1) then
#ifdef GRACKLE_MD
      if (imchem .eq. 1) then

!     --- C/O fine-structure, metal molecular rotational cooling for low temperatures ---

         itmask_tmp  = itmask
         itmask_line = itmask
         do i = is+1, ie+1
            if ( tgas(i) .lt. 3.d3 ) then
               itmask(i)      = .false.
            else
               itmask_line(i) = .false.
            endif
         enddo

         do i = is+1, ie+1
            if ( itmask_line(i) ) then

!     C/O fine-structure cooling

      if (logCI(i) .gt. - 8.d0) then
         lognhat = logCI(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int8(LCI_N), 
     &     LCI_D, LCI_dD, LCI_T, LCI_dT, LCI_H, LCI_dH,
     &     int8(LCI_Size), LCI_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int8(LCI_N), 
     &     LCI_D, LCI_dD, LCI_T, LCI_dT, LCI_H, LCI_dH,
     &     int8(LCI_Size), LCI_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LCI(i) = - (L - G) / dom * CI(i,j,k)/12.d0
         if (LCI(i) .ne. LCI(i)) LCI(i) = 0.d0
         edot(i) = edot(i) + LCI(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'CI ', LCI(i), edot(i)
#endif
      endif

      if (logCII(i) .gt. - 8.d0) then
         lognhat = logCII(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int8(LCII_N), 
     &     LCII_D, LCII_dD, LCII_T, LCII_dT, LCII_H, LCII_dH,
     &     int8(LCII_Size), LCII_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int8(LCII_N), 
     &     LCII_D, LCII_dD, LCII_T, LCII_dT, LCII_H, LCII_dH,
     &     int8(LCII_Size), LCII_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LCII(i) = - (L - G) / dom * CII(i,j,k)/12.d0
         if (LCII(i) .ne. LCII(i)) LCII(i) = 0.d0
         edot(i) = edot(i) + LCII(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'CII', LCII(i), edot(i)
#endif
      endif

      if (logOI(i) .gt. - 8.d0) then
         lognhat = logOI(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH(i), int8(LOI_N), 
     &     LOI_D, LOI_dD, LOI_T, LOI_dT, LOI_H, LOI_dH,
     &     int8(LOI_Size), LOI_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH(i), int8(LOI_N), 
     &     LOI_D, LOI_dD, LOI_T, LOI_dT, LOI_H, LOI_dH,
     &     int8(LOI_Size), LOI_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LOI(i) = - (L - G) / dom * OI(i,j,k)/16.d0
         if (LOI(i) .ne. LOI(i)) LOI(i) = 0.d0
         edot(i) = edot(i) + LOI(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'OI ', LOI(i), edot(i)
#endif
      endif

            endif
         enddo

      endif ! imchem
#endif

         if (clnew .eq. 1) then

            iZscale = 1
            call cool1d_cloudy_g(d, rhoH, metallicity,
     &           in, jn, kn, is, ie, j, k,
     &           logtem, edot, comp2, dom, zr,
     &           icmbTfloor, iClHeat, iZscale,
     &           metGridRank, metGridDim,
     &           metPar1, metPar2, metPar3,
     &           metDataSize, metCooling, metHeating,
     &           itmask)

         else

            call cool1D_cloudy_old_tables_g(
     &           d, de, rhoH, metallicity,
     &           in, jn, kn, is, ie, j, k,
     &           logtem, edot, comp2, ispecies, dom, zr,
     &           icmbTfloor, iClHeat, 
     &           clEleFra, metGridRank, metGridDim,
     &           metPar1, metPar2, metPar3, metPar4, metPar5,
     &           metDataSize, metCooling, metHeating, 
     &           itmask)

         endif

#ifdef GRACKLE_MD
      if (imchem .eq. 1) then

         itmask = itmask_tmp


         do i = is+1, ie+1
            if ( itmask(i) ) then

!     metal molecules rotational cooling

      if (logCO(i) .gt. - 8.d0) then
         lognhat = logCO(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH2(i), int8(LCO_N), 
     &     LCO_D, LCO_dD, LCO_T, LCO_dT, LCO_H, LCO_dH,
     &     int8(LCO_Size), LCO_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH2(i), int8(LCO_N), 
     &     LCO_D, LCO_dD, LCO_T, LCO_dT, LCO_H, LCO_dH,
     &     int8(LCO_Size), LCO_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LCO(i) = - (L - G) / dom * CO(i,j,k)/28.d0
         if (LCO(i) .ne. LCO(i)) LCO(i) = 0.d0
         edot(i) = edot(i) + LCO(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'CO ', LCO(i), edot(i)
#endif
      endif

      if (logOH(i) .gt. - 13.d0) then
         lognhat = logOH(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH2(i), int8(LOH_N), 
     &     LOH_D, LOH_dD, LOH_T, LOH_dT, LOH_H, LOH_dH,
     &     int8(LOH_Size), LOH_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH2(i), int8(LOH_N), 
     &     LOH_D, LOH_dD, LOH_T, LOH_dT, LOH_H, LOH_dH,
     &     int8(LOH_Size), LOH_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LOH(i) = - (L - G) / dom * OH(i,j,k)/17.d0
         if (LOH(i) .ne. LOH(i)) LOH(i) = 0.d0
         edot(i) = edot(i) + LOH(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'OH ', LOH(i), edot(i)
#endif
      endif

      if (logH2O(i) .gt. - 8.d0) then
         lognhat = logH2O(i) - logdvdr(i)

         call interpolate_3D_g(
     &     lognhat, logT(i), logH2(i), int8(LH2O_N), 
     &     LH2O_D, LH2O_dD, LH2O_T, LH2O_dT, LH2O_H, LH2O_dH,
     &     int8(LH2O_Size), LH2O_L, log_Linv)
         L = 1.d1**(-log_Linv)

         if (icmbTfloor .eq. 1) then
         call interpolate_3D_g(
     &     lognhat, logTcmb(i), logH2(i), int8(LH2O_N), 
     &     LH2O_D, LH2O_dD, LH2O_T, LH2O_dT, LH2O_H, LH2O_dH,
     &     int8(LH2O_Size), LH2O_L, log_Ginv)
         G = 1.d1**(-log_Ginv)
         else
         G = tiny8
         endif

         LH2O(i) = - (L - G) / dom * H2O(i,j,k)/18.d0
         if (LH2O(i) .ne. LH2O(i)) LH2O(i) = 0.d0
         edot(i) = edot(i) + LH2O(i)
#ifdef CHEM_DEBUG
         if((i.eq.iref).and.(j.eq.jref).and.(k.eq.kref))
     &   write(*,*) 'H2O', LH2O(i), edot(i)
#endif
      endif

            endif
         enddo

      endif ! imchem

#endif
      endif ! imcool

#ifdef GRACKLE_MD

!     Continuum opacity

      do i = is+1, ie+1
         if ( itmask(i) ) then
            if ( tau_con(i) .gt. 1.d0 ) then
!              if ( tau_con(i) .lt. 1.d2 ) then
!!                edot(i) = edot(i) * exp(-tau_con(i))
                  edot(i) = edot(i) / (tau_con(i)*tau_con(i))
!              else
!                 edot(i) = 0.d0
!              endif
            endif
         end if
      enddo
#endif

!     Add user-provided volumetric and/or specific heating terms

      if (iVheat .eq. 1) then

      do i = is+1, ie+1
         if ( itmask(i) ) then
            edot(i) = edot(i) + Vheat(i,j,k) / coolunit / dom**2
         end if
      enddo

      endif

      if (iMheat .eq. 1) then

      do i = is+1, ie+1
         if ( itmask(i) ) then
            edot(i) = edot(i) + Mheat(i,j,k) * d(i,j,k) * mh
     &          / coolunit / dom
         end if
      enddo

      endif

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) ) then
         tgasold(i) = tgas(i)
         end if
      enddo

!#ifdef CHEM_DEBUG
!     output cooling rate for debug
!!!   write(*,*) coolunit*dom/mh = vbase1**2 / tbase1
!!!                              = 2.5444239786382621E-004
      do i = is+1, ie+1
         if ( itmask(i) ) then

!        if((i.eq.1).and.(j.eq.1).and.(k.eq.1)) then
!     write(*,*) 0.75*d(i,j,k)*dom, tgas(i)
!    &   , LH2 (i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LHD (i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LCI (i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LCII(i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LOI (i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LCO (i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LOH (i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LH2O(i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , Ldst(i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , LCIE(i)/d(i,j,k)*2.5444239786382621E-004 * exp(-tau_con(i))
!    &   , tau_con(i)
!        endif

         endif
      enddo
!#endif

      return
      end
