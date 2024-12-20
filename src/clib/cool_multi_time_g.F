!=======================================================================
!////////////////////  SUBROUTINE COOL_MULTI_TIME_G  \\\\\\\\\\\\\\\\\\\

      subroutine cool_multi_time_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                  cooltime,
     &                in, jn, kn, nratec, iexpand, 
     &                ispecies, imetal, imcool, idust,
     &                idustall, idustfield, idustrec, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh, z_solar, fgr,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha, isrf,
     &                regra, gamma_isrfa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal, dust,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                h2ltea, gasgra,
     &                iradshield, avgsighi, avgsighei, avgsigheii,
     &                k24, k26, iradtrans, photogamma,
     &                ih2optical, iciecool, ih2cr, ihdcr, ciecoa,
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat,
     &                iTfloor, Tfloor_scalar, Tfloor
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
     &              , gasgr2a, gamma_isrf2a,
     &                iisrffield, isrf_habing)

!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
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
#ifdef _OPENMP
#include "omp_lib.h"
#endif

!  Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, nratec,
     &        iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        imcool, idust, idustall,
     &        idustfield, idustrec, igammah, ih2optical, iciecool,
     &        ih2cr, ihdcr,
     &        clnew, iVheat, iMheat, iTfloor, iradtrans, iradshield,
     &        iisrffield
     &       ,imchem, igrgr, ipcont

      real*8  aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma, fh, z_solar, fgr, clEleFra, Tfloor_scalar,
     &        tmcool
      R_PREC  d(in,jn,kn),   e(in,jn,kn),
     &        u(in,jn,kn),   v(in,jn,kn),    w(in,jn,kn),
     &        de(in,jn,kn),  HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn),  H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn),  HDI(in,jn,kn),
     &        metal(in,jn,kn), dust(in,jn,kn),
     &        Vheat(in,jn,kn), Mheat(in,jn,kn), Tfloor(in,jn,kn),
     &        isrf_habing(in,jn,kn),
     &        cooltime(in,jn,kn)
      real*8  photogamma(in,jn,kn)
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
     &        avgsighi, avgsighei, avgsigheii,
     &        k24, k26
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
       real*8  cieY06a(nratec)
       integer*8 LH2_N(3), LH2_Size
       real*8  LH2_D(LH2_N(1)), LH2_T(LH2_N(2)), LH2_H(LH2_N(3))
     &       , LH2_dD, LH2_dT, LH2_dH, LH2_L(LH2_Size)
       integer*8 LHD_N(3), LHD_Size
       real*8  LHD_D(LHD_N(1)), LHD_T(LHD_N(2)), LHD_H(LHD_N(3))
     &       , LHD_dD, LHD_dT, LHD_dH, LHD_L(LHD_Size)
       integer*8 LCI_N(3), LCI_Size
       real*8  LCI_D(LCI_N(1)), LCI_T(LCI_N(2)), LCI_H(LCI_N(3))
     &       , LCI_dD, LCI_dT, LCI_dH, LCI_L(LCI_Size)
       integer*8 LCII_N(3), LCII_Size
       real*8  LCII_D(LCII_N(1)), LCII_T(LCII_N(2)), LCII_H(LCII_N(3))
     &       , LCII_dD, LCII_dT, LCII_dH, LCII_L(LCII_Size)
       integer*8 LOI_N(3), LOI_Size
       real*8  LOI_D(LOI_N(1)), LOI_T(LOI_N(2)), LOI_H(LOI_N(3))
     &       , LOI_dD, LOI_dT, LOI_dH, LOI_L(LOI_Size)
       integer*8 LCO_N(3), LCO_Size
       real*8  LCO_D(LCO_N(1)), LCO_T(LCO_N(2)), LCO_H(LCO_N(3))
     &       , LCO_dD, LCO_dT, LCO_dH, LCO_L(LCO_Size)
       integer*8 LOH_N(3), LOH_Size
       real*8  LOH_D(LOH_N(1)), LOH_T(LOH_N(2)), LOH_H(LOH_N(3))
     &       , LOH_dD, LOH_dT, LOH_dH, LOH_L(LOH_Size)
       integer*8 LH2O_N(3), LH2O_Size
       real*8  LH2O_D(LH2O_N(1)), LH2O_T(LH2O_N(2)), LH2O_H(LH2O_N(3))
     &       , LH2O_dD, LH2O_dT, LH2O_dH, LH2O_L(LH2O_Size)
       integer*8 alphap_N(2), alphap_Size
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
!     grain temperature
      real*8  tSiM(in), tFeM(in), tMg2SiO4(in)
     &      , tMgSiO3(in), tFe3O4(in), tAC(in)
     &      , tSiO2D(in), tMgO(in), tFeS(in)
     &      , tAl2O3(in)
     &      , treforg(in), tvolorg(in), tH2Oice(in)
      real*8  gasgr2a(nratec), gamma_isrf2a

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

!  Locals

      integer i, j, k
      integer t, dj, dk
      real*8 comp1, comp2, energy
      R_PREC factor

!  Slice locals
 
      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &       tgas(in), tgasold(in), mmw(in),
     &       tdust(in), metallicity(in), dust2gas(in), rhoH(in),
     &       mynh(in), myde(in), gammaha_eff(in),
     &       gasgr_tdust(in), regr(in)
      real*8 edot(in)

!  Cooling/heating slice locals

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), cieco(in),
     &       hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in)
      integer dummy_iter_arg
          
!  Iteration mask for multi_cool

      MASK_TYPE itmask(in), itmask_metal(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
      dk = ke - ks + 1
      dj = je - js + 1

!     Convert densities from comoving to 'proper'

      if (iexpand .eq. 1) then

         factor = real(aye**(-3), RKIND)

         call scale_fields_g(
     &     ispecies, imetal, idustfield, imchem,
     &     imabund, idspecies, immulti, igrgr, idsub, factor,
     &     is, ie, js, je, ks, ke, in, jn, kn,
     &     d, de, HI, HII, HeI, HeII, HeIII,
     &     HM, H2I, H2II, DI, DII, HDI, DM, HDII, HeHII,
     &     metal, dust,
     &     CI, CII, CO, CO2,
     &     OI, OH, H2O, O2,
     &     SiI, SiOI, SiO2I,
     &     CH, CH2, COII, OII,
     &     OHII, H2OII, H3OII, O2II,
     &     Mg, Al, S, Fe,
     &     SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &     AC, SiO2D, MgO, FeS, Al2O3,
     &     reforg, volorg, H2Oice,
     &     metal_loc, metal_C13, metal_C20, metal_C25, metal_C30,
     &     metal_F13, metal_F15, metal_F50, metal_F80,
     &     metal_P170, metal_P200, metal_Y19)
 
      endif

!     Loop over slices (in the k-direction)

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   comp1, comp2, energy,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef, p2d,
!$omp&   tgas, tgasold,
!$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
!$omp&   mynh, myde, gammaha_eff, gasgr_tdust, regr, edot,
!$omp&   ceHI, ceHeI, ceHeII,
!$omp&   ciHI, ciHeI, ciHeIS, ciHeII,
!$omp&   reHII, reHeII1, reHeII2, reHeIII,
!$omp&   brem, cieco,
!$omp&   hyd01k, h2k01, vibh, roth, rotl,
!$omp&   gpldl, gphdl, hdlte, hdlow,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

         do i = is+1, ie+1
            itmask(i) = MASK_TRUE
         end do

!        Compute the cooling rate
         dummy_iter_arg=1

         call cool1d_multi_g(     
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool,
     &                idust, idustall, idustfield, idustrec,
     &                idim, is, ie, j, k, ih2co, ipiht,
     &                dummy_iter_arg, igammah,
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
     &                iisrffield, isrf_habing, itmask, itmask_metal
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

!        Compute the cooling time on the slice
!          (the gamma used here is the same as used to calculate the pressure
!           in cool1d_multi_g)

         do i = is+1, ie+1
            energy = max(p2d(i)/(gamma-1._DKIND), tiny)   
            cooltime(i,j,k) = real(energy/edot(i), RKIND)
         enddo

      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

!     Convert densities back to comoving from 'proper'

      if (iexpand .eq. 1) then

         factor = real(aye**3, RKIND)

         call scale_fields_g(
     &     ispecies, imetal, idustfield, imchem,
     &     imabund, idspecies, immulti, igrgr, idsub, factor,
     &     is, ie, js, je, ks, ke, in, jn, kn,
     &     d, de, HI, HII, HeI, HeII, HeIII,
     &     HM, H2I, H2II, DI, DII, HDI, DM, HDII, HeHII,
     &     metal, dust,
     &     CI, CII, CO, CO2,
     &     OI, OH, H2O, O2,
     &     SiI, SiOI, SiO2I,
     &     CH, CH2, COII, OII,
     &     OHII, H2OII, H3OII, O2II,
     &     Mg, Al, S, Fe,
     &     SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &     AC, SiO2D, MgO, FeS, Al2O3,
     &     reforg, volorg, H2Oice,
     &     metal_loc, metal_C13, metal_C20, metal_C25, metal_C30,
     &     metal_F13, metal_F15, metal_F50, metal_F80,
     &     metal_P170, metal_P200, metal_Y19)

      endif

      return
      end

