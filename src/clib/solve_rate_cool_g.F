#include "phys_const.def"
#define CHEM_DEBUG
#define CORRECTED_BY_CHIAKI2020FEB
#define CORRECTED_BY_CHIAKI2020MAY

!=======================================================================
!///////////////////  SUBROUTINE SOLVE_RATE_COOL_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_rate_cool_g(icool, d, e, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, 
     &                ispecies, imetal, imcool, idust,
     &                idustall, idustfield, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                dx, dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh, dtoh, z_solar, fgr,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k22a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                k57a, k58a,
     &                ndratec, dtemstart, dtemend, h2dusta,
     &                ncrna, ncrd1a, ncrd2a,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha, isrf,
     &                regra, gamma_isrfa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal, dust,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela, 
     &                h2ltea, gasgra, iH2shield, iradshield,
     &                avgsighi, avgsighei, avgsigheii,
     &                iradtrans, iradcoupled, iradstep,
     &                irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                photogamma, xH2shield,
     &                ierr,
     &                ih2optical, iciecool, ithreebody, ciecoa, 
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat,
     &                iisrffield, isrf_habing
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
     &              , loc_XC, loc_XO, loc_XMg, loc_XAl, loc_XSi
     &              , loc_XS, loc_XFe
     &              , loc_fC, loc_fO, loc_fMg, loc_fAl, loc_fSi
     &              , loc_fS, loc_fFe
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
     &              , Hgrain, Tgrain, Ograin, Lgrain
     &              , impop3
     &              , metal_loc, metal_C30, metal_F13
     &              , C30_XC , C30_XO, C30_XMg, C30_XAl
     &              , C30_XSi, C30_XS, C30_XFe
     &              , C30_fC , C30_fO, C30_fMg, C30_fAl
     &              , C30_fSi, C30_fS, C30_fFe
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
     &              , HC30_SiM, HC30_FeM, HC30_Mg2SiO4, HC30_MgSiO3
     &              , HC30_Fe3O4, HC30_AC, HC30_SiO2D, HC30_MgO
     &              , HC30_FeS, HC30_Al2O3
     &              , OC30_SiM, OC30_FeM, OC30_Mg2SiO4, OC30_MgSiO3
     &              , OC30_Fe3O4, OC30_AC, OC30_SiO2D, OC30_MgO
     &              , OC30_FeS, OC30_Al2O3
     &              , LC30_SiM, LC30_FeM, LC30_Mg2SiO4, LC30_MgSiO3
     &              , LC30_Fe3O4, LC30_AC, LC30_SiO2D, LC30_MgO
     &              , LC30_FeS, LC30_Al2O3
     &              , KC30_SiM, KC30_FeM, KC30_Mg2SiO4, KC30_MgSiO3
     &              , KC30_Fe3O4, KC30_AC, KC30_SiO2D, KC30_MgO
     &              , KC30_FeS, KC30_Al2O3
     &              , F13_XC , F13_XO, F13_XMg, F13_XAl
     &              , F13_XSi, F13_XS, F13_XFe
     &              , F13_fC , F13_fO, F13_fMg, F13_fAl
     &              , F13_fSi, F13_fS, F13_fFe
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
     &              , HF13_SiM, HF13_FeM, HF13_Mg2SiO4, HF13_MgSiO3
     &              , HF13_Fe3O4, HF13_AC, HF13_SiO2D, HF13_MgO
     &              , HF13_FeS, HF13_Al2O3
     &              , OF13_SiM, OF13_FeM, OF13_Mg2SiO4, OF13_MgSiO3
     &              , OF13_Fe3O4, OF13_AC, OF13_SiO2D, OF13_MgO
     &              , OF13_FeS, OF13_Al2O3
     &              , LF13_SiM, LF13_FeM, LF13_Mg2SiO4, LF13_MgSiO3
     &              , LF13_Fe3O4, LF13_AC, LF13_SiO2D, LF13_MgO
     &              , LF13_FeS, LF13_Al2O3
     &              , KF13_SiM, KF13_FeM, KF13_Mg2SiO4, KF13_MgSiO3
     &              , KF13_Fe3O4, KF13_AC, KF13_SiO2D, KF13_MgO
     &              , KF13_FeS, KF13_Al2O3
#endif
     &          )

!
!  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
!  modified2:  October, 1996 by GB; adapted to AMR
!  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
!  modified4:  June,    2005 by GB to solve rate & cool at same time
!  modified5:  April,   2009 by JHW to include radiative transfer
!  modified6:  September, 2009 by BDS to include cloudy cooling
!
!  PURPOSE:
!    Solve the multi-species rate and cool equations.
!
!  INPUTS:
!    icool    - flag to update energy from radiative cooling
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!    HI,HII   - H density fields (neutral & ionized)
!    HeI/II/III - He density fields
!    DI/II    - D density fields (neutral & ionized)
!    HDI      - neutral HD molecule density field
!    HM       - H- density field
!    H2I      - H_2 (molecular H) density field
!    H2II     - H_2+ density field
!    metal    - metal density field
!    dust     - dust density field
!    kph*     - photoionization fields
!    gamma*   - photoheating fields
!
!    is,ie    - start and end indices of active region (zero based)
!    iexpand  - comoving coordinates flag (0 = off, 1 = on)
!    idim     - dimensionality (rank) of problem
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!    imcool   - flag if there is metal cooling
!    idust    - flag for H2 formation on dust grains
!    idustall - flag for dust (0 - none, 1 - heating/cooling + H2 form.)
!    idustfield - flag if a dust density field is present
!    iisrffield - flag if a field for the interstellar radiation field is present
!    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
!    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
!    iH2shield - flag for approximate self-shielding of H2 (Wolcott-Green+ 2011)
!    iradshield - flag for approximate self-shielding of UV background
!    avgsighi   - spectrum averaged ionization crs for HI for use with shielding
!    avgsighei  - spectrum averaged ionization crs for HeI for use with shielding
!    avgsigheii - spectrum averaged ionization crs for HeII for use with shielding
!    iradtrans - flag to include radiative transfer (1 = on, 0 = off)
!    iradcoupled - flag to indicate coupled radiative transfer
!    iradstep  - flag to indicate intermediate coupled radiative transfer timestep
!    irt_honly - flag to indicate applying RT ionization and heating to HI only

!    fh       - Hydrogen mass fraction (typically 0.76)
!    dtoh     - Deuterium to H mass ratio
!    z_solar  - Solar metal mass fraction
!    fgr      - the local dust to gas ratio (by mass)
!    dt       - timestep to integrate over
!    aye      - expansion factor (in code units)
!
!    utim     - time units (i.e. code units to CGS conversion factor)
!    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
!    urho     - density units
!    uxyz     - length units
!    utem     - temperature(-like) units
!
!    temstart, temend - start and end of temperature range for rate table
!    nratec   - dimensions of chemical rate arrays (functions of temperature)
!    dtemstart, dtemend - start and end of dust temperature range
!    ndratec  - extra dimension for H2 formation on dust rate (dust temperature)
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    priGridRank - rank of cloudy primordial cooling data grid
!    priGridDim  - array containing dimensions of cloudy primordial data
!    priPar1, priPar2, priPar3 - arrays containing primordial grid parameter values
!    priDataSize - total size of flattened 1D primordial cooling data array
!    priCooling  - primordial cooling data
!    priHeating  - primordial heating data
!    priMMW      - primordial mmw data
!    metGridRank - rank of cloudy metal cooling data grid
!    metGridDim  - array containing dimensions of cloudy metal data
!    metPar1, metPar2, metPar3 - arrays containing metal grid parameter values
!    metDataSize - total size of flattened 1D metal cooling data array
!    metCooling  - metal cooling data
!    metHeating  - metal heating data
!    iVheat      - flag for using volumetric heating rate
!    iMheat      - flag for using specific heating rate
!    Vheat       - array of volumetric heating rates
!    Mheat       - array of specific heating rates
!
!  OUTPUTS:
!    update chemical rate densities (HI, HII, etc)
!
!  PARAMETERS:
!    itmax   - maximum allowed sub-cycle iterations
!    mh      - H mass in cgs units
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"
#ifdef _OPENMP
#include "omp_lib.h"
#endif

!  General Arguments

      integer icool, in, jn, kn, is, js, ks, ie, je, ke, nratec, 
     &        iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        ierr, imcool, idust, idustall, idustfield,
     &        igammah, ih2optical, iciecool, ithreebody,
     &        ndratec, clnew, iVheat, iMheat, iH2shield, iradshield,
     &        iradtrans, iradcoupled, iradstep, irt_honly,
     &        iisrffield
#ifdef GRACKLE_MD
     &       ,imchem, igrgr
#endif

      real*8  dx, dt, aye, temstart, temend, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh, z_solar, 
     &        fgr, dtemstart, dtemend, clEleFra

!  Density, energy and velocity fields fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        d(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn), dust(in,jn,kn),
     &        Vheat(in,jn,kn), Mheat(in,jn,kn)
#ifdef GRACKLE_MD
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
      real*8  loc_XC , loc_XO, loc_XMg, loc_XAl
     &      , loc_XSi, loc_XS, loc_XFe
      real*8  loc_fC , loc_fO, loc_fMg, loc_fAl
     &      , loc_fSi, loc_fS, loc_fFe
       R_PREC DM(in,jn,kn), HDII(in,jn,kn), HeHII(in,jn,kn)
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
#endif

!  Radiative transfer fields

      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn), photogamma(in,jn,kn)

!  H2 self-shielding length-scale field

      R_PREC  xH2shield(in,jn,kn)

!  Interstellar radiation field for dust heating

      R_PREC  isrf_habing(in,jn,kn)

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
#ifdef GRACKLE_MD
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
      integer impop3
      real*8  C30_XC , C30_XO, C30_XMg, C30_XAl
     &      , C30_XSi, C30_XS, C30_XFe
      real*8  C30_fC , C30_fO, C30_fMg, C30_fAl
     &      , C30_fSi, C30_fS, C30_fFe
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
      real*8  HC30_SiM(C30_Size), HC30_FeM(C30_Size)
     &      , HC30_Mg2SiO4(C30_Size), HC30_MgSiO3(C30_Size)
     &      , HC30_Fe3O4(C30_Size), HC30_AC(C30_Size)
     &      , HC30_SiO2D(C30_Size), HC30_MgO(C30_Size)
     &      , HC30_FeS(C30_Size), HC30_Al2O3(C30_Size)
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
      real*8  KC30_SiM(C30_Size), KC30_FeM(C30_Size)
     &      , KC30_Mg2SiO4(C30_Size), KC30_MgSiO3(C30_Size)
     &      , KC30_Fe3O4(C30_Size), KC30_AC(C30_Size)
     &      , KC30_SiO2D(C30_Size), KC30_MgO(C30_Size)
     &      , KC30_FeS(C30_Size), KC30_Al2O3(C30_Size)
      real*8  F13_XC , F13_XO, F13_XMg, F13_XAl
     &      , F13_XSi, F13_XS, F13_XFe
      real*8  F13_fC , F13_fO, F13_fMg, F13_fAl
     &      , F13_fSi, F13_fS, F13_fFe
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
      real*8  HF13_SiM(F13_Size), HF13_FeM(F13_Size)
     &      , HF13_Mg2SiO4(F13_Size), HF13_MgSiO3(F13_Size)
     &      , HF13_Fe3O4(F13_Size), HF13_AC(F13_Size)
     &      , HF13_SiO2D(F13_Size), HF13_MgO(F13_Size)
     &      , HF13_FeS(F13_Size), HF13_Al2O3(F13_Size)
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
      real*8  KF13_SiM(F13_Size), KF13_FeM(F13_Size)
     &      , KF13_Mg2SiO4(F13_Size), KF13_MgSiO3(F13_Size)
     &      , KF13_Fe3O4(F13_Size), KF13_AC(F13_Size)
     &      , KF13_SiO2D(F13_Size), KF13_MgO(F13_Size)
     &      , KF13_FeS(F13_Size), KF13_Al2O3(F13_Size)
#endif

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
#ifdef GRACKLE_MD
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
#endif

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
     &       eqt1, eqt2, eqtdef, dheq, heq, dlogtem, dx_cgs,
     &       c_ljeans

!  row temporaries

      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in),
     &       dtit(in), ttot(in), p2d(in), tgas(in), tgasold(in),
     &       tdust(in), metallicity(in), dust2gas(in),
     &       rhoH(in), mmw(in), mynh(in), myde(in), gammaha_eff(in),
     &       gasgr_tdust(in), regr(in), olddtit

!  Rate equation row temporaries

      real*8 HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &       HMp(in), H2Ip(in), H2IIp(in),
     &       dep(in), dedot(in),HIdot(in), dedot_prev(in),
     &       DIp(in), DIIp(in), HDIp(in), HIdot_prev(in),
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k28shield(in), k29shield(in), k30shield(in),
     &       k31shield(in),
     &       k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in),
     &       k13dd(in, 14), h2dust(in),
     &       ncrn(in), ncrd1(in), ncrd2(in)
#ifdef GRACKLE_MD
      real*8  DMp(in)   , HDIIp(in)   , HeHIIp(in)
     &      , CIp(in)   , CIIp(in)    , COp(in)
     &      , CO2p(in)  , OIp(in)     , OHp(in)
     &      , H2Op(in)  , O2p(in)     , SiIp(in)
     &      , SiOIp(in) , SiO2Ip(in)  , CHp(in)
     &      , CH2p(in)  , COIIp(in)   , OIIp(in)
     &      , OHIIp(in) , H2OIIp(in)  , H3OIIp(in)
     &      , O2IIp(in) , Mgp(in)     , Alp(in)
     &      , Sp(in)    , Fep(in)     , SiMp(in)
     &      , FeMp(in)  , Mg2SiO4p(in), MgSiO3p(in)
     &      , Fe3O4p(in), ACp(in)     , SiO2Dp(in)
     &      , MgOp(in)  , FeSp(in)    , Al2O3p(in)

      real*8 k125(in), k129(in), k130(in), k131(in), k132(in)
     &     , k133(in), k134(in), k135(in), k136(in), k137(in)
     &     , k148(in), k149(in), k150(in), k151(in), k152(in)
     &     , k153(in)
     &     , kz15(in), kz16(in), kz17(in), kz18(in), kz19(in)
     &     , kz20(in), kz21(in), kz22(in), kz23(in), kz24(in)
     &     , kz25(in), kz26(in), kz27(in), kz28(in), kz29(in)
     &     , kz30(in), kz31(in), kz32(in), kz33(in), kz34(in)
     &     , kz35(in), kz36(in), kz37(in), kz38(in), kz39(in)
     &     , kz40(in), kz41(in), kz42(in), kz43(in), kz44(in)
     &     , kz45(in), kz46(in), kz47(in), kz48(in), kz49(in)
     &     , kz50(in), kz51(in), kz52(in), kz53(in), kz54(in)
      real*8  kd0(in), kd1(in), kd2(in), kd3(in), kd4(in)
     &      , kd5(in), kd6(in), kd7(in), kd8(in), kd9(in)
#endif

!  Cooling/heating row locals

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), edot(in)
      real*8 hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in), cieco(in)

!  Iteration mask

      logical itmask(in), anydust
#ifdef SWITCH_TO_NR_SCHEME
      logical itmask_tmp(in), itmask_nr(in)
      integer itr, imp_eng, itr_time
      integer nsp, isp, jsp, id
      real*8  dspj, err, err_max
#ifdef CHEM_DEBUG
      integer ierr_max
#endif
#ifdef GRACKLE_MD
      integer,parameter   :: i_eng = 49
#else
      integer,parameter   :: i_eng = 13
#endif
      real*8  dsp(i_eng), dsp0(i_eng), dsp1(i_eng)
     &      , dspdot(i_eng), dspdot1(i_eng)
     &      , ddsp(i_eng), der(i_eng,i_eng)
      integer,allocatable :: idsp(:)
      real*8, allocatable :: mtrx(:,:), vec(:)
#ifdef CHEM_DEBUG
!     real*8, allocatable :: dummtrx(:,:), dumvec(:)
!     real*8 vecDtot
#endif
      real*8, parameter   :: eps = 1.d-4
#endif
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set error indicator

      ierr = 0

!     Set flag for dust-related options

      anydust = (idust .gt. 0) .or. (idustall .gt. 0)
      
!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      uvel     = uxyz / utim
#ifdef CORRECTED_BY_CHIAKI2020MAY
      chunit = (1.60218e-12_DKIND)/((xbase1/tbase1)**2*mh)   ! 1 eV per REACTION
#else
c     chunit = (7.17775e-12_DKIND)/(2._DKIND*uvel*uvel*mh)   ! 4.5 eV per H2 formed
      chunit = (1.60218e-12_DKIND)/(2._DKIND*uvel*uvel*mh)   ! 1 eV per H2 formed
#endif

      dx_cgs = dx * xbase1
      c_ljeans = sqrt((gamma * pi * kboltz) /
     &     (GravConst * mh * dbase1))

      dlogtem = (log(temend) - log(temstart))/real(nratec-1, DKIND)

#ifdef CORRECTED_BY_CHIAKI2020FEB
      if (ispecies .gt. 0) then

!     Correct the species to ensure consistency (i.e. type conservation)

#define ABUNDANCE_CORRECTION
#ifdef ABUNDANCE_CORRECTION
      call make_consistent_g(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, 
     &                     d, is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, fh, dtoh
#ifdef GRACKLE_MD
     &                   , imchem, igrgr, dom
     &                   , DM, HDII, HeHII
     &                   , CI, CII, CO, CO2
     &                   , OI, OH, H2O, O2
     &                   , SiI, SiOI, SiO2I
     &                   , CH, CH2, COII, OII
     &                   , OHII, H2OII, H3OII, O2II
     &                   , Mg, Al, S, Fe
     &                   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                   , AC, SiO2D, MgO, FeS, Al2O3
     &                   , impop3
     &                   , metal_loc, metal_C30, metal_F13
     &                   , loc_XC, loc_XO, loc_XMg, loc_XAl, loc_XSi
     &                   , loc_XS, loc_XFe
     &                   , loc_fC, loc_fO, loc_fMg, loc_fAl, loc_fSi
     &                   , loc_fS, loc_fFe
     &                   , C30_XC, C30_XO, C30_XMg, C30_XAl, C30_XSi
     &                   , C30_XS, C30_XFe
     &                   , C30_fC, C30_fO, C30_fMg, C30_fAl, C30_fSi
     &                   , C30_fS, C30_fFe
     &                   , F13_XC, F13_XO, F13_XMg, F13_XAl, F13_XSi
     &                   , F13_XS, F13_XFe
     &                   , F13_fC, F13_fO, F13_fMg, F13_fAl, F13_fSi
     &                   , F13_fS, F13_fFe
#endif
     &                      )
#endif

      endif
#endif

!  Convert densities from comoving to proper

      if (iexpand .eq. 1) then

         call scale_fields_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal, dust,
#ifdef GRACKLE_MD
     &                  imchem, igrgr,
     &                  DM, HDII, HeHII,
     &                  CI, CII, CO, CO2,
     &                  OI, OH, H2O, O2,
     &                  SiI, SiOI, SiO2I,
     &                  CH, CH2, COII, OII,
     &                  OHII, H2OII, H3OII, O2II,
     &                  Mg, Al, S, Fe,
     &                  SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &                  AC, SiO2D, MgO, FeS, Al2O3,
     &                  impop3,
     &                  metal_loc, metal_C30, metal_F13,
#endif
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, idustfield,
     &                  aye**(-3))

      endif

!!!         write(*,*) 'CHECK1', d(1,1,1), e(1,1,1), tbase1, xbase1
#ifdef ABUNDANCE_CORRECTION
      call ceiling_species_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal,
#ifdef GRACKLE_MD
     &                     imchem, igrgr,
     &                     DM, HDII, HeHII,
     &                     CI, CII, CO, CO2,
     &                     OI, OH, H2O, O2,
     &                     SiI, SiOI, SiO2I,
     &                     CH, CH2, COII, OII,
     &                     OHII, H2OII, H3OII, O2II,
     &                     Mg, Al, S, Fe,
     &                     SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &                     AC, SiO2D, MgO, FeS, Al2O3,
     &                     impop3,
     &                     metal_loc, metal_C30, metal_F13,
#endif
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal)
#endif
!!!         write(*,*) 'CHECK2', d(1,1,1), e(1,1,1)


!  Loop over zones, and do an entire i-column in one go
      dk = ke - ks + 1
      dj = je - js + 1

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
! ierr is declared as shared and should be modified with atomic operation
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k, iter,
!$omp&   ttmin, energy, comp1, comp2,
!$omp&   heq1, heq2, eqk221, eqk222, eqk131, eqk132,
!$omp&   eqt1, eqt2, eqtdef, dheq, heq,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef,
!$omp&   dtit, ttot, p2d, tgas, tgasold,
!$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
!$omp&   mynh, myde, gammaha_eff, gasgr_tdust, regr, olddtit,
!$omp&   HIp, HIIp, HeIp, HeIIp, HeIIIp,
!$omp&   HMp, H2Ip, H2IIp,
!$omp&   dep, dedot,HIdot, dedot_prev,
!$omp&   DIp, DIIp, HDIp, HIdot_prev,
!$omp&   k24shield, k25shield, k26shield,
!$omp&   k28shield, k29shield, k30shield,
!$omp&   k31shield,
!$omp&   k1 , k2 , k3 , k4 , k5,
!$omp&   k6 , k7 , k8 , k9 , k10,
!$omp&   k11, k12, k13, k14, k15,
!$omp&   k16, k17, k18, k19, k22,
!$omp&   k50, k51, k52, k53, k54,
!$omp&   k55, k56, k57, k58, k13dd, h2dust,
!$omp&   ncrn, ncrd1, ncrd2,
!$omp&   ceHI, ceHeI, ceHeII,
!$omp&   ciHI, ciHeI, ciHeIS, ciHeII,
!$omp&   reHII, reHeII1, reHeII2, reHeIII,
!$omp&   brem, edot,
!$omp&   hyd01k, h2k01, vibh, roth, rotl,
!$omp&   gpldl, gphdl, hdlte, hdlow, cieco,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

!       tolerance = 1.0e-06_DKIND * dt

!     Initialize iteration mask to true for all cells.

         do i = is+1, ie+1
            itmask(i) = .true.
         enddo

!      If we are using coupled radiation with intermediate stepping,
!      set iteration mask to include only cells with radiation in the
!      intermediate coupled chemistry / energy step 
         if (iradtrans .eq. 1) then
            if (iradcoupled .eq. 1 .and. iradstep .eq. 1) then
               do i = is+1, ie+1
                  if (kphHI(i,j,k) .gt. 0) then
                      itmask(i) = .true.
                  else
                      itmask(i) = .false.
                  endif
               enddo
            endif

!      Normal rate solver, but don't double count cells with radiation
            if (iradcoupled .eq. 1 .and. iradstep .eq. 0) then
               do i = is+1, ie + 1
                  if (kphHI(i,j,k) .gt. 0) then
                     itmask(i) = .false.
                  else
                     itmask(i) = .true.
                  endif
               enddo
            endif
         endif ! end rad trans check (divergent from original code)
          


!        Set time elapsed to zero for each cell in 1D section

         do i = is+1, ie+1
            ttot(i) = 0._DKIND
         enddo

!        ------------------ Loop over subcycles ----------------

         do iter = 1, itmax

            do i = is+1, ie+1
               if (itmask(i)) then
                  dtit(i) = huge8
               endif
            enddo

!           Compute the cooling rate, tgas, tdust, and metallicity for this row

            call cool1d_multi_g(
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
     &          )


            if (ispecies .gt. 0) then

!        Look-up rates as a function of temperature for 1D set of zones
!         (maybe should add itmask to this call)

            call lookup_cool_rates1d_g(temstart, temend, nratec, j, k,
     &               is, ie, ithreebody,
     &               in, jn, kn, ispecies, anydust,
     &               iH2shield, iradshield,
     &               tgas, mmw, d, HI, HII, HeI, HeII, HeIII, 
     &               HM, H2I, H2II, DI, DII, HDI,
     &               tdust, dust2gas,
     &               k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &               k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &               k17a, k18a, k19a, k22a,
     &               k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &               k57a, k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &               ncrna, ncrd1a, ncrd2a, 
     &               avgsighi, avgsighei, avgsigheii, piHI, piHeI,
     &               k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &               k11, k12, k13, k14, k15, k16, k17, k18,
     &               k19, k22, k24, k25, k26, k28, k29, k30, k31,
     &               k50, k51, k52, k53, k54, k55, k56, k57,
     &               k58, k13dd, k24shield, k25shield, k26shield,
     &               k28shield, k29shield, k30shield,
     &               k31shield, h2dust, ncrn, ncrd1, ncrd2, 
     &               t1, t2, tdef, logtem, indixe, 
     &               dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
     &               iradtrans, kdissH2I, xH2shield, itmask
#ifdef GRACKLE_MD
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
     &              , HC30_SiM, HC30_FeM, HC30_Mg2SiO4, HC30_MgSiO3
     &              , HC30_Fe3O4, HC30_AC, HC30_SiO2D, HC30_MgO
     &              , HC30_FeS, HC30_Al2O3
     &              , KC30_SiM, KC30_FeM, KC30_Mg2SiO4, KC30_MgSiO3
     &              , KC30_Fe3O4, KC30_AC, KC30_SiO2D, KC30_MgO
     &              , KC30_FeS, KC30_Al2O3
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
     &              , HF13_SiM, HF13_FeM, HF13_Mg2SiO4, HF13_MgSiO3
     &              , HF13_Fe3O4, HF13_AC, HF13_SiO2D, HF13_MgO
     &              , HF13_FeS, HF13_Al2O3
     &              , KF13_SiM, KF13_FeM, KF13_Mg2SiO4, KF13_MgSiO3
     &              , KF13_Fe3O4, KF13_AC, KF13_SiO2D, KF13_MgO
     &              , KF13_FeS, KF13_Al2O3
     &              , kd0, kd1, kd2, kd3, kd4, kd5, kd6, kd7, kd8, kd9
#endif
     &      )

!           Compute dedot and HIdot, the rates of change of de and HI
!             (should add itmask to this call)

            call rate_timestep_g(
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
     &                     itmask, edot, chunit, dom)

#ifdef SWITCH_TO_NR_SCHEME
!           move itmask temporary array
!           then split cells with low densities
!              => Gauss-Seidel scheme
!                        and with high densities
!              => Newton-Raphson scheme

            itmask_tmp = itmask
            itmask_nr  = itmask
            do i = is+1, ie+1
            if ( itmask_tmp(i) ) then

!              if (d(i,j,k)*dom .gt. 1.e6_DKIND) then
!!             if (d(i,j,k)*dom .gt.
!!   &          min(1.e6_DKIND/(metal(i,j,k)/d(i,j,k)/0.02d-4)**2
!!   &             ,1.e6_DKIND)) then
               if (
     &         ( ( metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k) )
     &      .and.( d(i,j,k)*dom .gt. 1.e6_DKIND) ) .or.
     &         ( ( metal(i,j,k) .le. 1.e-9_DKIND * d(i,j,k) )
     &      .and.( d(i,j,k)*dom .gt. 1.e8_DKIND) )
     &            ) then
                  itmask(i) = .false.
               else
                  itmask_nr(i) = .false.
               endif

            endif
            enddo
#endif

!           Find timestep that keeps relative chemical changes below 10%

            do i = is+1, ie+1
               if (itmask(i)) then
!              Bound from below to prevent numerical errors
               
	       if (abs(dedot(i)) .lt. tiny8) 
     &             dedot(i) = min(tiny,de(i,j,k))
	       if (abs(HIdot(i)) .lt. tiny8)
     &             HIdot(i) = min(tiny,HI(i,j,k))

!              If the net rate is almost perfectly balanced then set
!                  it to zero (since it is zero to available precision)

               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)),
     &                 abs(k2(i)*HII(i,j,k)*de(i,j,k)))/
     &             max(abs(dedot(i)),abs(HIdot(i))) .gt.
     &              1.0e6_DKIND) then
                  dedot(i) = tiny8
                  HIdot(i) = tiny8
               endif

!              If the iteration count is high then take the smaller of
!                the calculated dedot and last time step's actual dedot.
!                This is intended to get around the problem of a low
!                electron or HI fraction which is in equilibrium with high
!                individual terms (which all nearly cancel).

               if (iter .gt. 50) then
                  dedot(i) = min(abs(dedot(i)), abs(dedot_prev(i)))
                  HIdot(i) = min(abs(HIdot(i)), abs(HIdot_prev(i)))
               endif

!              compute minimum rate timestep

               olddtit = dtit(i)
               dtit(i) = min(abs(0.1_DKIND*de(i,j,k)/dedot(i)), 
     &              abs(0.1_DKIND*HI(i,j,k)/HIdot(i)),
     &              dt-ttot(i), 0.5_DKIND*dt)

               if (d(i,j,k)*dom .gt. 1e8_DKIND .and. 
     &              edot(i) .gt. 0._DKIND) then
!              Equilibrium value for H is:
!              H = (-1._DKIND / (4*k22)) * (k13 - sqrt(8 k13 k22 rho + k13^2))
!              We now want this to change by 10% or less, but we're only
!              differentiating by dT.  We have de/dt.  We need dT/de.
!              T = (g-1)*p2d*utem/N; tgas == (g-1)(p2d*utem/N)
!              dH_eq / dt = (dH_eq/dT) * (dT/de) * (de/dt)
!              dH_eq / dT (see above; we can calculate the derivative here)
!              dT / de = utem * (gamma - 1._DKIND) / N == tgas / p2d
!              de / dt = edot
!              Now we use our estimate of dT/de to get the estimated
!              difference in the equilibrium
                eqt2 = min(log(tgas(i)) + 0.1_DKIND*dlogtem, t2(i))
                eqtdef = (eqt2 - t1(i))/(t2(i) - t1(i))
                eqk222 = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk132 = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq2 = (-1._DKIND / (4._DKIND*eqk222)) * (eqk132-
     &               sqrt(8._DKIND*eqk132*eqk222*
     &               fh*d(i,j,k)+eqk132**2._DKIND))

                eqt1 = max(log(tgas(i)) - 0.1_DKIND*dlogtem, t1(i))
                eqtdef = (eqt1 - t1(i))/(t2(i) - t1(i))
                eqk221 = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk131 = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq1 = (-1._DKIND / (4._DKIND*eqk221)) * (eqk131-
     &               sqrt(8._DKIND*eqk131*eqk221*
     &               fh*d(i,j,k)+eqk131**2._DKIND))

                dheq = (abs(heq2-heq1)/(exp(eqt2) - exp(eqt1)))
     &               * (tgas(i)/p2d(i)) * edot(i)
                heq = (-1._DKIND / (4._DKIND*k22(i))) * (k13(i)-
     &               sqrt(8._DKIND*k13(i)*k22(i)*
     &               fh*d(i,j,k)+k13(i)**2._DKIND))
                !write(0,*) heq2, heq1, eqt2, eqt1, tgas(i), p2d(i), 
!     &                     edot(i)
                if (d(i,j,k)*dom.gt.1e18_DKIND.and.i.eq.4) then
#ifdef _OPENMP
!$omp critical
#endif
                  write(0,*) HI(i,j,k)/heq, edot(i), tgas(i)
#ifdef _OPENMP
!$omp end critical
#endif
                endif
                dtit(i) = min(dtit(i), 0.1_DKIND*heq/dheq)
              endif
              if (iter.gt.10_DKIND) then
                 dtit(i) = min(olddtit*1.5_DKIND, dtit(i))
              endif

#define DONT_WRITE_COOLING_DEBUG
#ifdef WRITE_COOLING_DEBUG
!              Output some debugging information if required

               if (dtit(i)/dt .lt. 1.0e-2_DKIND .and.
     &             iter .gt. 800 .and.
     &             abs((dt-ttot(i))/dt) .gt. 1.0e-3_DKIND) then
#ifdef _OPENMP
!$omp critical
#endif
                  write(4,1000) iter,i,j,k,dtit(i),
     &              ttot(i),dt,de(i,j,k),dedot(i),HI(i,j,k),HIdot(i),
     &              tgas(i), dedot_prev(i), HIdot_prev(i)
                  write(4,1100) HI(i,j,k),HII(i,j,k),
     &              HeI(i,j,k),HeII(i,j,k),HeIII(i,j,k),
     &              HM(i,j,k),H2I(i,j,k),H2II(i,j,k),de(i,j,k)
                  write(4,1100)
     &               -      k1(i) *de(i,j,k)    *HI(i,j,k)  ,
     &               -      k7(i) *de(i,j,k)    *HI(i,j,k),
     &               -      k8(i) *HM(i,j,k)    *HI(i,j,k),
     &               -      k9(i) *HII(i,j,k)   *HI(i,j,k),
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2._DKIND,
     &               - 2._DKIND*k22(i)*HI(i,j,k)**2 *HI(i,j,k),
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) ,
     &               + 2._DKIND*k13(i)*HI(i,j,k)*H2I(i,j,k)/2._DKIND,
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._DKIND,
     &               + 2._DKIND*k12(i)*de(i,j,k)*H2I(i,j,k)/2._DKIND,
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k),
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k),
     &               + 2._DKIND*k16(i)*HM(i,j,k)    *HII(i,j,k),
     &               + 2._DKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._DKIND,
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._DKIND,
     &               -      k57(i)*HI(i,j,k)    *HI(i,j,k),
     &               -      k58(i)*HI(i,j,k)    *HeI(i,j,k)/4._DKIND
#ifdef _OPENMP
!$omp end critical
#endif
               endif
 1000          format(i5,3(i3,1x),1p,11(e11.3))
 1100          format(1p,22(e11.3))
#endif /* WRITE_COOLING_DEBUG */
#ifdef SWITCH_TO_NR_SCHEME
            else if (itmask_nr(i)) then
               dtit(i) = dt - ttot(i);
#endif
            else               ! itmask
               dtit(i) = dt;
            endif
            enddo               ! end loop over i

            endif ! end if (ispecies .gt. 0)

!           Compute maximum timestep for cooling/heating

            do i = is+1, ie+1
               if (itmask(i)) then
!              Set energy per unit volume of this cell based in the pressure
!              (the gamma used here is the right one even for H2 since p2d 
!               is calculated with this gamma).

               energy = max(p2d(i)/(gamma-1._DKIND), tiny8)

!              If the temperature is at the bottom of the temperature look-up 
!              table and edot < 0, then shut off the cooling.

               if (tgas(i) .le. 1.01_DKIND*temstart .and.
     &              edot(i) .lt. 0._DKIND) 
     &              edot(i) = tiny8
	       if (abs(edot(i)) .lt. tiny8) edot(i) = tiny8

!              Compute timestep for 10% change

                  dtit(i) = min(real(abs(0.1_DKIND*
     &              energy/edot(i)), DKIND), 
     &              dt-ttot(i), dtit(i))

               if (dtit(i) .ne. dtit(i)) then
#ifdef _OPENMP
!$omp critical
#endif
                 write(6,*) 'HUGE dtit :: ', energy, edot(i), dtit(i),
     &                      dt, ttot(i), abs(0.1_DKIND*energy/edot(i)), 
     &                      real(abs(0.1_DKIND*energy/edot(i)), DKIND)
#ifdef _OPENMP
!$omp end critical
#endif
               endif

#ifdef WRITE_COOLING_DEBUG
!              If the timestep is too small, then output some debugging info

               if (((dtit(i)/dt .lt. 1.0e-2_DKIND .and. iter .gt. 800) 
     &               .or. iter .gt. itmax-100) .and.
     &              abs((dt-ttot(i))/dt) .gt. 1.0e-3_DKIND) then
#ifdef _OPENMP
!$omp critical
#endif
                 write(3,2000) i,j,k,iter,e(i,j,k),edot(i),tgas(i),
     &             energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k),dtit(i)
#ifdef _OPENMP
!$omp end critical
#endif
               endif
 2000          format(4(i4,1x),1p,10(e14.3))
#endif /* WRITE_COOLING_DEBUG */
            endif   ! itmask
            enddo   ! end loop over i

!           Update total and gas energy
!!!         if((j.eq.1) .and. (k.eq.1))
!!!  &      write(*,*) 'CHECK4', d(1,1,1), e(1,1,1), edot(1)
            if (icool .eq. 1) then
            do i = is+1, ie+1
               if (itmask(i)) then
               e(i,j,k)  = e(i,j,k) +
     &                 real(edot(i)/d(i,j,k)*dtit(i), RKIND)
#ifdef WRITE_COOLING_DEBUG
               if (e(i,j,k) .ne. e(i,j,k)) then
#ifdef _OPENMP
!$omp critical
#endif
                    write(3,*) edot(i),d(i,j,k),dtit(i)
#ifdef _OPENMP
!$omp end critical
#endif
               endif
#endif /* WRITE_COOLING_DEBUG */

            endif               ! itmask
            enddo
            endif
!!          if((j.eq.10) .and. (k.eq.7)) then
!!             i = 10 
!!             write(*,*) 'CHECK5', d(i,j,k)*dom, tgas(i)
!!   &             , de(i,j,k)/d(i,j,k)
!!   &             , HI(i,j,k)/d(i,j,k)
!!   &             , HII(i,j,k)/d(i,j,k)
!!          endif

            if (ispecies .gt. 0) then

!           Solve rate equations with one linearly implicit Gauss-Seidel 
!           sweep of a backward Euler method ---

            call step_rate_g(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k,
     &                     ispecies, anydust,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58, 
     &                     h2dust, rhoH,
     &                     k24shield, k25shield, k26shield, 
     &                     k28shield, k29shield, k30shield, k31shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, 
     &                     kphHI, kphHeI, kphHeII,
     &                     itmask
#ifdef GRACKLE_MD
     &                   , DM, HDII, HeHII, imetal, metal, imchem, igrgr
     &                   , CI, CII, CO, CO2
     &                   , OI, OH, H2O, O2
     &                   , SiI, SiOI, SiO2I
     &                   , CH, CH2, COII, OII
     &                   , OHII, H2OII, H3OII, O2II
     &                   , Mg, Al, S, Fe
     &                   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                   , AC, SiO2D, MgO, FeS, Al2O3
     &                   , k125, k129, k130, k131, k132
     &                   , k133, k134, k135, k136, k137
     &                   , k148, k149, k150, k151, k152
     &                   , k153
     &                   , kz15 , kz16 , kz17 , kz18 , kz19 
     &                   , kz20 , kz21 , kz22 , kz23 , kz24 
     &                   , kz25 , kz26 , kz27 , kz28 , kz29 
     &                   , kz30 , kz31 , kz32 , kz33 , kz34 
     &                   , kz35 , kz36 , kz37 , kz38 , kz39 
     &                   , kz40 , kz41 , kz42 , kz43 , kz44 
     &                   , kz45 , kz46 , kz47 , kz48 , kz49 
     &                   , kz50 , kz51 , kz52 , kz53 , kz54 
     &                   , DMp, HDIIp, HeHIIp
     &                   , CIp, CIIp, COp, CO2p
     &                   , OIp, OHp, H2Op, O2p
     &                   , SiIp, SiOIp, SiO2Ip
     &                   , CHp, CH2p, COIIp, OIIp
     &                   , OHIIp, H2OIIp, H3OIIp, O2IIp
     &                   , Mgp, Alp, Sp, Fep
     &                   , SiMp, FeMp, Mg2SiO4p, MgSiO3p, Fe3O4p
     &                   , ACp, SiO2Dp, MgOp, FeSp, Al2O3p
     &                   , kd0, kd1, kd2, kd3, kd4
     &                   , kd5, kd6, kd7, kd8, kd9
#endif
     &           )
!!!         if((j.eq.1) .and. (k.eq.1))
!!!  &      write(*,*) 'CHECKS', d(1,1,1), e(1,1,1)

#ifdef SWITCH_TO_NR_SCHEME
            do i = is+1, ie+1
            if (itmask_nr(i)) then

!!             if((i.eq.12).and.(j.eq.6).and.(k.eq.6))
!              if( ((j-1) + (k-1)*jn)*in+i .eq. 1410 )
!    &            write(*,*) 'CHECKa', d(i,j,k)*dom, tgas(i)
!    &              , de (i,j,k)/d(i,j,k)
!    &              , HI (i,j,k)/d(i,j,k)
!    &              , H2I(i,j,k)/d(i,j,k)

               if (   (d(i,j,k)*dom .gt. 1.d7)
     &           .and.(tgas(i) .gt. 1650.d0) ) then
                  imp_eng = 1
               else
                  imp_eng = 0
               endif

!           If density and temperature are low, update total and gas energy

CCCC           if (icool .eq. 1) then
CCCC              if (imp_eng .eq. 0) then
CCCC                 e(i,j,k)  = e(i,j,k) +
CCCC &                    real(edot(i)/d(i,j,k)*dtit(i), RKIND)
CCCC              endif
CCCC           endif

!              initialize arrays
               if (ispecies .gt. 0) nsp = 6
               if (ispecies .gt. 1) nsp = nsp + 3
               if (ispecies .gt. 2) nsp = nsp + 3
#ifdef GRACKLE_MD
               if (ispecies .gt. 3) nsp = nsp + 3
               if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
                 if (imchem .eq. 1) nsp = nsp + 19
                 if (igrgr  .eq. 1) nsp = nsp + 14
               endif
#endif
               nsp = nsp + imp_eng
               allocate(idsp(nsp))
               allocate(mtrx(nsp,nsp))
               allocate(vec(nsp))
#ifdef CHEM_DEBUG
!              allocate(dummtrx(nsp,nsp))
!              allocate(dumvec(nsp))
#endif

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
#ifdef GRACKLE_MD
               if ( ispecies .gt. 3 ) then
                  dsp(13) = DM(i,j,k)
                  dsp(14) = HDII(i,j,k)
                  dsp(15) = HeHII(i,j,k)
               endif
               if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
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
               endif
               if ( igrgr .eq. 1 ) then
                  dsp(35) = Mg(i,j,k)
                  dsp(36) = Al(i,j,k)
                  dsp(37) = S(i,j,k)
                  dsp(38) = Fe(i,j,k)
                  dsp(39) = SiM(i,j,k)
                  dsp(40) = FeM(i,j,k)
                  dsp(41) = Mg2SiO4(i,j,k)
                  dsp(42) = MgSiO3(i,j,k)
                  dsp(43) = Fe3O4(i,j,k)
                  dsp(44) = AC(i,j,k)
                  dsp(45) = SiO2D(i,j,k)
                  dsp(46) = MgO(i,j,k)
                  dsp(47) = FeS(i,j,k)
                  dsp(48) = Al2O3(i,j,k)
               endif
               endif
#endif
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
#ifdef GRACKLE_MD
               if (ispecies .gt. 3) then
                  do isp = 13, 15
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
               if (imchem .eq. 1) then
                  do isp = 16, 34
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               if (igrgr .eq. 1) then
                  do isp = 35, 48
                     id = id + 1
                     idsp(id) = isp
                  enddo
               endif
               endif
#endif
               if ( imp_eng .eq.1 ) then
                  id = id + 1
                  idsp(id) = i_eng
               endif

               dsp0 = dsp
               ddsp = 0.d0

!              iteration to timestep
               ierr=1
               itr_time=0
               do while ((ierr.eq.1))

                  ! if timestep becomes too short, do not update

                  if(itr_time .gt. 5) then
                     write(6,*) 'DO NOT CONVERGE '
                     write(6,*) i,j,k, d(i,j,k)*dom, tgas(i)
                     write(6,*) aye, utem, uxyz, uaye, urho, utim
                     write(6,*) dt, d(i,j,k), e(i,j,k)
                     if ( ispecies .gt. 0 ) then
                        write(6,*) de(i,j,k), HI(i,j,k), HII(i,j,k)
     &                           , HeI(i,j,k), HeII(i,j,k), HeIII(i,j,k)
                     endif
                     if ( ispecies .gt. 1 ) then
                        write(6,*) HM(i,j,k), H2I(i,j,k), H2II(i,j,k)
                     endif
                     if ( ispecies .gt. 2 ) then
                        write(6,*) DI(i,j,k), DII(i,j,k), HDI(i,j,k)
                     endif
#ifdef GRACKLE_MD
                     if ( ispecies .gt. 3 ) then
                        write(6,*) DM(i,j,k), HDII(i,j,k), HeHII(i,j,k)
                     endif
                     if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
                     if ( imchem .eq. 1 ) then
                        write(6,*) CI(i,j,k), CII(i,j,k), CO(i,j,k)
     &                           , CO2(i,j,k)
                        write(6,*) OI(i,j,k), OH(i,j,k), H2O(i,j,k)
     &                           , O2(i,j,k)
                        write(6,*) SiI(i,j,k), SiOI(i,j,k), SiO2I(i,j,k)
                        write(6,*) CH(i,j,k), CH2(i,j,k), COII(i,j,k)
                        write(6,*) OII(i,j,k), OHII(i,j,k), H2OII(i,j,k)
     &                           , H3OII(i,j,k), O2II(i,j,k)
                     endif
                     if ( igrgr .eq. 1 ) then
                        write(6,*) Mg(i,j,k), Al(i,j,k)
     &                           , S(i,j,k), Fe(i,j,k)
                        write(6,*) SiM(i,j,k), FeM(i,j,k)
     &                           , Mg2SiO4(i,j,k), MgSiO3(i,j,k)
     &                           , Fe3O4(i,j,k)
                        write(6,*) AC(i,j,k), SiO2D(i,j,k), MgO(i,j,k)
     &                           , FeS(i,j,k), Al2O3(i,j,k)
                     endif
                     endif
#endif
                     if ( imetal .eq. 1 ) then
                        write(6,*) metal(i,j,k)
                     endif
                     if (iradtrans .eq. 1) then
                        write(6,*) kphHI(i,j,k), kphHeI(i,j,k)
     &                           , kphHeII(i,j,k), kdissH2I(i,j,k)
     &                           , photogamma(i,j,k)
                     endif

                     stop
                     dtit(i) = huge8
                     ttot(i) = huge8
                     go to 9997
                  endif

                  dsp = dsp0
                  ddsp = 0.d0

!!             if((i.eq.12).and.(j.eq.6).and.(k.eq.6))
!              if( ((j-1) + (k-1)*jn)*in+i .eq. 1410 )
!    &            write(*,*) 'CHECKb', d(i,j,k)*dom, tgas(i)
!    &              , dsp( 1)/d(i,j,k)
!    &              , dsp( 2)/d(i,j,k)
!    &              , dsp( 8)/d(i,j,k)

!                 iteration to solve ODEs
                  err_max=1.d2
                  itr=0
                  do while ((itr.lt.20) .and. (err_max.gt.1.d-8))

!!                  write(*,*) 'A', 0.75*d(i,j,k)*dom, itr_time, itr
!!                  write(*,*) 'a', idsp, dtit(i)/dt
!!                  write(*,*) 'a', dsp/d(i,j,k)

            call lookup_cool_rates0d(itr_time+itr, dtit(i),
     &   d(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), 
     &   nsp, dsp, dspdot, nratec, 
     &   iexpand, ispecies, imetal, imcool,
     &   idust, idustall, idustfield,
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
     &   ih2optical, iciecool, ciecoa, cieco(i),
     &   icmbTfloor, iClHeat, clEleFra,
     &   priGridRank, priGridDim,
     &   priPar1, priPar2, priPar3, priPar4, priPar5,
     &   priDataSize, priCooling, priHeating, priMMW,
     &   metGridRank, metGridDim,
     &   metPar1, metPar2, metPar3, metPar4, metPar5,
     &   metDataSize, metCooling, metHeating, clnew,
     &   iVheat, iMheat, Vheat(i,j,k), Mheat(i,j,k),
     &   iisrffield, isrf_habing,
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
     &   photogamma(i,j,k), xH2shield(i,j,k), chunit, itmask_nr(i)
#ifdef GRACKLE_MD
     & , imchem, igrgr
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
     & , grain_N, grain_Size
     & , grain_D, grain_T, grain_dD, grain_dT
     & , Hgrain, Tgrain, Ograin, Lgrain
     & , impop3
     & , metal_loc(i,j,k), metal_C30(i,j,k), metal_F13(i,j,k)
     & , C30_XC , C30_XO, C30_XMg, C30_XAl
     & , C30_XSi, C30_XS, C30_XFe
     & , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     & , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     & , C30_fFeS, C30_fAl2O3
     & , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4, C30_r0MgSiO3
     & , C30_r0Fe3O4, C30_r0AC, C30_r0SiO2D, C30_r0MgO
     & , C30_r0FeS, C30_r0Al2O3
     & , C30_a0SiM, C30_a0FeM, C30_a0Mg2SiO4, C30_a0MgSiO3
     & , C30_a0Fe3O4, C30_a0AC, C30_a0SiO2D, C30_a0MgO
     & , C30_a0FeS, C30_a0Al2O3
     & , C30_v0SiM, C30_v0FeM, C30_v0Mg2SiO4, C30_v0MgSiO3
     & , C30_v0Fe3O4, C30_v0AC, C30_v0SiO2D, C30_v0MgO
     & , C30_v0FeS, C30_v0Al2O3
     & , C30_N, C30_Size
     & , C30_D, C30_T, C30_dD, C30_dT
     & , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     & , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     & , C30_RFeS, C30_RAl2O3
     & , C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     & , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     & , C30_dRFeS, C30_dRAl2O3
     & , HC30_SiM, HC30_FeM, HC30_Mg2SiO4, HC30_MgSiO3
     & , HC30_Fe3O4, HC30_AC, HC30_SiO2D, HC30_MgO
     & , HC30_FeS, HC30_Al2O3
     & , OC30_SiM, OC30_FeM, OC30_Mg2SiO4, OC30_MgSiO3
     & , OC30_Fe3O4, OC30_AC, OC30_SiO2D, OC30_MgO
     & , OC30_FeS, OC30_Al2O3
     & , LC30_SiM, LC30_FeM, LC30_Mg2SiO4, LC30_MgSiO3
     & , LC30_Fe3O4, LC30_AC, LC30_SiO2D, LC30_MgO
     & , LC30_FeS, LC30_Al2O3
     & , KC30_SiM, KC30_FeM, KC30_Mg2SiO4, KC30_MgSiO3
     & , KC30_Fe3O4, KC30_AC, KC30_SiO2D, KC30_MgO
     & , KC30_FeS, KC30_Al2O3
     & , F13_XC , F13_XO, F13_XMg, F13_XAl
     & , F13_XSi, F13_XS, F13_XFe
     & , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     & , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     & , F13_fFeS, F13_fAl2O3
     & , F13_r0SiM, F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     & , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0MgO
     & , F13_r0FeS, F13_r0Al2O3
     & , F13_a0SiM, F13_a0FeM, F13_a0Mg2SiO4, F13_a0MgSiO3
     & , F13_a0Fe3O4, F13_a0AC, F13_a0SiO2D, F13_a0MgO
     & , F13_a0FeS, F13_a0Al2O3
     & , F13_v0SiM, F13_v0FeM, F13_v0Mg2SiO4, F13_v0MgSiO3
     & , F13_v0Fe3O4, F13_v0AC, F13_v0SiO2D, F13_v0MgO
     & , F13_v0FeS, F13_v0Al2O3
     & , F13_N, F13_Size
     & , F13_D, F13_T, F13_dD, F13_dT
     & , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     & , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     & , F13_RFeS, F13_RAl2O3
     & , F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     & , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     & , F13_dRFeS, F13_dRAl2O3
     & , HF13_SiM, HF13_FeM, HF13_Mg2SiO4, HF13_MgSiO3
     & , HF13_Fe3O4, HF13_AC, HF13_SiO2D, HF13_MgO
     & , HF13_FeS, HF13_Al2O3
     & , OF13_SiM, OF13_FeM, OF13_Mg2SiO4, OF13_MgSiO3
     & , OF13_Fe3O4, OF13_AC, OF13_SiO2D, OF13_MgO
     & , OF13_FeS, OF13_Al2O3
     & , LF13_SiM, LF13_FeM, LF13_Mg2SiO4, LF13_MgSiO3
     & , LF13_Fe3O4, LF13_AC, LF13_SiO2D, LF13_MgO
     & , LF13_FeS, LF13_Al2O3
     & , KF13_SiM, KF13_FeM, KF13_Mg2SiO4, KF13_MgSiO3
     & , KF13_Fe3O4, KF13_AC, KF13_SiO2D, KF13_MgO
     & , KF13_FeS, KF13_Al2O3
#endif
     & , imp_eng)

!!!            if((i.eq.10).and.(j.eq.10).and.(k.eq.7)) then
!!!               write(*,*) 'CHECKc', d(i,j,k)*dom, tgas(i)
!!!  &              , dspdot( 1)/d(i,j,k)*dtit(i)
!!!  &              , dspdot( 2)/d(i,j,k)*dtit(i)
!!!  &              , dspdot( 3)/d(i,j,k)*dtit(i)
!!!               write(*,*) 'CHECKd', d(i,j,k)*dom, tgas(i)
!!!  &              , (dsp( 1)+dspdot( 1)*dtit(i))/d(i,j,k)
!!!  &              , (dsp( 2)+dspdot( 2)*dtit(i))/d(i,j,k)
!!!  &              , (dsp( 3)+dspdot( 3)*dtit(i))/d(i,j,k)
!!!             endif

                    do jsp = 1, nsp
                      dspj = eps * dsp(idsp(jsp))
                      do isp = 1, nsp
                        if(isp .eq. jsp) then
                          dsp1(idsp(isp)) = dsp(idsp(isp)) + dspj
                        else
                          dsp1(idsp(isp)) = dsp(idsp(isp))
                        endif
                      enddo

            call lookup_cool_rates0d(1, dtit(i),
     &   d(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), 
     &   nsp, dsp1, dspdot1, nratec, 
     &   iexpand, ispecies, imetal, imcool,
     &   idust, idustall, idustfield,
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
     &   ih2optical, iciecool, ciecoa, cieco(i),
     &   icmbTfloor, iClHeat, clEleFra,
     &   priGridRank, priGridDim,
     &   priPar1, priPar2, priPar3, priPar4, priPar5,
     &   priDataSize, priCooling, priHeating, priMMW,
     &   metGridRank, metGridDim,
     &   metPar1, metPar2, metPar3, metPar4, metPar5,
     &   metDataSize, metCooling, metHeating, clnew,
     &   iVheat, iMheat, Vheat(i,j,k), Mheat(i,j,k),
     &   iisrffield, isrf_habing,
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
     &   photogamma(i,j,k), xH2shield(i,j,k), chunit, itmask_nr(i)
#ifdef GRACKLE_MD
     & , imchem, igrgr
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
     & , grain_N, grain_Size
     & , grain_D, grain_T, grain_dD, grain_dT
     & , Hgrain, Tgrain, Ograin, Lgrain
     & , impop3
     & , metal_loc(i,j,k), metal_C30(i,j,k), metal_F13(i,j,k)
     & , C30_XC , C30_XO, C30_XMg, C30_XAl
     & , C30_XSi, C30_XS, C30_XFe
     & , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     & , C30_fFe3O4, C30_fAC, C30_fSiO2D, C30_fMgO
     & , C30_fFeS, C30_fAl2O3
     & , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4, C30_r0MgSiO3
     & , C30_r0Fe3O4, C30_r0AC, C30_r0SiO2D, C30_r0MgO
     & , C30_r0FeS, C30_r0Al2O3
     & , C30_a0SiM, C30_a0FeM, C30_a0Mg2SiO4, C30_a0MgSiO3
     & , C30_a0Fe3O4, C30_a0AC, C30_a0SiO2D, C30_a0MgO
     & , C30_a0FeS, C30_a0Al2O3
     & , C30_v0SiM, C30_v0FeM, C30_v0Mg2SiO4, C30_v0MgSiO3
     & , C30_v0Fe3O4, C30_v0AC, C30_v0SiO2D, C30_v0MgO
     & , C30_v0FeS, C30_v0Al2O3
     & , C30_N, C30_Size
     & , C30_D, C30_T, C30_dD, C30_dT
     & , C30_RSiM, C30_RFeM, C30_RMg2SiO4, C30_RMgSiO3
     & , C30_RFe3O4, C30_RAC, C30_RSiO2D, C30_RMgO
     & , C30_RFeS, C30_RAl2O3
     & , C30_dRSiM, C30_dRFeM, C30_dRMg2SiO4, C30_dRMgSiO3
     & , C30_dRFe3O4, C30_dRAC, C30_dRSiO2D, C30_dRMgO
     & , C30_dRFeS, C30_dRAl2O3
     & , HC30_SiM, HC30_FeM, HC30_Mg2SiO4, HC30_MgSiO3
     & , HC30_Fe3O4, HC30_AC, HC30_SiO2D, HC30_MgO
     & , HC30_FeS, HC30_Al2O3
     & , OC30_SiM, OC30_FeM, OC30_Mg2SiO4, OC30_MgSiO3
     & , OC30_Fe3O4, OC30_AC, OC30_SiO2D, OC30_MgO
     & , OC30_FeS, OC30_Al2O3
     & , LC30_SiM, LC30_FeM, LC30_Mg2SiO4, LC30_MgSiO3
     & , LC30_Fe3O4, LC30_AC, LC30_SiO2D, LC30_MgO
     & , LC30_FeS, LC30_Al2O3
     & , KC30_SiM, KC30_FeM, KC30_Mg2SiO4, KC30_MgSiO3
     & , KC30_Fe3O4, KC30_AC, KC30_SiO2D, KC30_MgO
     & , KC30_FeS, KC30_Al2O3
     & , F13_XC , F13_XO, F13_XMg, F13_XAl
     & , F13_XSi, F13_XS, F13_XFe
     & , F13_fSiM, F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     & , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fMgO
     & , F13_fFeS, F13_fAl2O3
     & , F13_r0SiM, F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     & , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0MgO
     & , F13_r0FeS, F13_r0Al2O3
     & , F13_a0SiM, F13_a0FeM, F13_a0Mg2SiO4, F13_a0MgSiO3
     & , F13_a0Fe3O4, F13_a0AC, F13_a0SiO2D, F13_a0MgO
     & , F13_a0FeS, F13_a0Al2O3
     & , F13_v0SiM, F13_v0FeM, F13_v0Mg2SiO4, F13_v0MgSiO3
     & , F13_v0Fe3O4, F13_v0AC, F13_v0SiO2D, F13_v0MgO
     & , F13_v0FeS, F13_v0Al2O3
     & , F13_N, F13_Size
     & , F13_D, F13_T, F13_dD, F13_dT
     & , F13_RSiM, F13_RFeM, F13_RMg2SiO4, F13_RMgSiO3
     & , F13_RFe3O4, F13_RAC, F13_RSiO2D, F13_RMgO
     & , F13_RFeS, F13_RAl2O3
     & , F13_dRSiM, F13_dRFeM, F13_dRMg2SiO4, F13_dRMgSiO3
     & , F13_dRFe3O4, F13_dRAC, F13_dRSiO2D, F13_dRMgO
     & , F13_dRFeS, F13_dRAl2O3
     & , HF13_SiM, HF13_FeM, HF13_Mg2SiO4, HF13_MgSiO3
     & , HF13_Fe3O4, HF13_AC, HF13_SiO2D, HF13_MgO
     & , HF13_FeS, HF13_Al2O3
     & , OF13_SiM, OF13_FeM, OF13_Mg2SiO4, OF13_MgSiO3
     & , OF13_Fe3O4, OF13_AC, OF13_SiO2D, OF13_MgO
     & , OF13_FeS, OF13_Al2O3
     & , LF13_SiM, LF13_FeM, LF13_Mg2SiO4, LF13_MgSiO3
     & , LF13_Fe3O4, LF13_AC, LF13_SiO2D, LF13_MgO
     & , LF13_FeS, LF13_Al2O3
     & , KF13_SiM, KF13_FeM, KF13_Mg2SiO4, KF13_MgSiO3
     & , KF13_Fe3O4, KF13_AC, KF13_SiO2D, KF13_MgO
     & , KF13_FeS, KF13_Al2O3
#endif
     & , imp_eng)

!!!                   write(*,*) 'g', d(i,j,k)*dom, tgas, dspdot1(1:15)

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

!!                  write(*,*) 'h'
!!                  do jsp = 1, 13
!!                    write(*,*) der(:,jsp)
!!                  enddo

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

#ifdef CHEM_DEBUG
!                   if((i.eq.1).and.(j.eq.1).and.(k.eq.1)) then
!                !!!   write(*,'(A6,16E10.2)') "CHECKd", vec(idsp(:))
!                !!!   dummtrx = mtrx
!                !!!   dumvec  = vec
!                      vecDtot =
!    &                     vec(10)/0.75/d(i,j,k)/ 2._DKIND
!    &                   + vec(11)/0.75/d(i,j,k)/ 2._DKIND
!    &                   + vec(12)/0.75/d(i,j,k)/ 3._DKIND
!    &                   + vec(13)/0.75/d(i,j,k)/ 2._DKIND
!    &                   + vec(14)/0.75/d(i,j,k)/ 3._DKIND
!                   endif
#endif

                    call gaussj_g(nsp, mtrx, vec, ierr)
                    if(ierr .eq. 1) then
                      go to 9998
                    endif

#ifdef CHEM_DEBUG
!                   if((i.eq.1).and.(j.eq.1).and.(k.eq.1)) then
!!!                    do isp = 1, nsp
!!!                      do jsp = 1, nsp
!!!                        dumvec(isp) = dumvec(isp)
!!!  &                       - dummtrx(isp,jsp)*vec(jsp)
!!!                      enddo
!!!                    enddo
!!!                    write(*,'(17E10.2)') 0.75*d(i,j,k)*dom
!!!  &                     , dumvec(idsp(:))
!                     if((itr_time.eq.0).and.(itr.eq.0)) 
!    &                 write(*,*) 0.75*d(i,j,k)*dom
!    &                , dspdot(10)/0.75/d(i,j,k)/ 2._DKIND*dtit(i)
!    &                + dspdot(11)/0.75/d(i,j,k)/ 2._DKIND*dtit(i)
!    &                + dspdot(12)/0.75/d(i,j,k)/ 3._DKIND*dtit(i)
!    &                + dspdot(13)/0.75/d(i,j,k)/ 2._DKIND*dtit(i)
!    &                + dspdot(14)/0.75/d(i,j,k)/ 3._DKIND*dtit(i)
!    &                   , vecDtot
!    &                   , vec(10)/0.75/d(i,j,k)/ 2._DKIND
!    &                   + vec(11)/0.75/d(i,j,k)/ 2._DKIND
!    &                   + vec(12)/0.75/d(i,j,k)/ 3._DKIND
!    &                   + vec(13)/0.75/d(i,j,k)/ 2._DKIND
!    &                   + vec(14)/0.75/d(i,j,k)/ 3._DKIND
!                   endif
#endif

                    do isp = 1, nsp
                      ddsp(idsp(isp)) = ddsp(idsp(isp)) + vec(isp)
                      dsp(idsp(isp))  = dsp(idsp(isp))  + vec(isp)
                    enddo

!!             if((i.eq.12).and.(j.eq.6).and.(k.eq.6))
!              if( ((j-1) + (k-1)*jn)*in+i .eq. 1410 )
!    &            write(*,*) 'CHECKf', d(i,j,k)*dom, tgas(i)
!    &              , dsp( 1)/d(i,j,k)
!    &              , dsp( 2)/d(i,j,k)
!    &              , dsp( 8)/d(i,j,k)

                    err_max = 0.d0
#ifdef CHEM_DEBUG
                    ierr_max = 0
#endif
                    do isp = 1, nsp
                      if(dsp(idsp(isp)) .gt. tiny8) then
                        err = dabs(vec(isp) / dsp(idsp(isp)))
                      else
                        err = 0.d0
                      endif
                      if(err .gt. err_max) then
                        err_max = err
#ifdef CHEM_DEBUG
                        ierr_max = idsp(isp)
#endif
                      endif
                    enddo
#ifdef CHEM_DEBUG
!                   write(*,*) 0.75*d(i,j,k)*dom, tgas(i)
!    &                   , itr_time, itr
!    &                   , ierr_max
!    &                   , dsp(ierr_max)/d(i,j,k)/0.75_DKIND
!!   &                   , (dsp(idsp(isp))/d(i,j,k), isp=1,nsp)
!!   &                   , dsp(2)/d(i,j,k)/0.75_DKIND
!!   &                   , dsp(8)/d(i,j,k)/0.75_DKIND/2._DKIND
#endif
                    itr=itr+1
                  enddo

9998              continue

#ifdef UNUSED
                  if ( itr .ge. 20 ) then
                    ierr = 1
                  endif
                  do isp = 1, nsp
                    if (  ( dsp(idsp(isp)) .lt. 0.d0 ) 
     &               .or. ( dsp(idsp(isp)) .ne. dsp(idsp(isp)) ) ) then
                      ierr = 1
                    endif
                  enddo
#endif
                  if(ierr .eq. 1) then
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
#ifdef GRACKLE_MD
               if ( ispecies .gt. 3 ) then
                   DM(i,j,k)      = dsp(13)
                   HDII(i,j,k)    = dsp(14)
                   HeHII(i,j,k)   = dsp(15)
               endif
               if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
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
               endif
               if ( igrgr .eq. 1 ) then
                   Mg(i,j,k)      = dsp(35)
                   Al(i,j,k)      = dsp(36)
                   S(i,j,k)       = dsp(37)
                   Fe(i,j,k)      = dsp(38)
                   SiM(i,j,k)     = dsp(39)
                   FeM(i,j,k)     = dsp(40)
                   Mg2SiO4(i,j,k) = dsp(41)
                   MgSiO3(i,j,k)  = dsp(42)
                   Fe3O4(i,j,k)   = dsp(43)
                   AC(i,j,k)      = dsp(44)
                   SiO2D(i,j,k)   = dsp(45)
                   MgO(i,j,k)     = dsp(46)
                   FeS(i,j,k)     = dsp(47)
                   Al2O3(i,j,k)   = dsp(48)
               endif
               endif
#endif
               e(i,j,k)     = dsp(i_eng)

               if (icool .eq. 1) then
                  if (imp_eng .eq. 0) then
                     e(i,j,k)  = e(i,j,k) +
     &                    real(edot(i)/d(i,j,k)*dtit(i), RKIND)
                  endif
               endif

9997           continue

               deallocate(idsp)
               deallocate(vec)
               deallocate(mtrx)
#ifdef CHEM_DEBUG
!              deallocate(dumvec)
!              deallocate(dummtrx)
#endif

!!             if((i.eq.12).and.(j.eq.6).and.(k.eq.6))
!              if( ((j-1) + (k-1)*jn)*in+i .eq. 1410 )
!    &            write(*,*) 'CHECKg', d(i,j,k)*dom, e(i,j,k)
!    &              ,  de(i,j,k)/d(i,j,k)
!    &              ,  HI(i,j,k)/d(i,j,k)
!    &              , H2I(i,j,k)/d(i,j,k)

            endif
            enddo
#endif
            endif ! if (ispecies .gt. 0) then
!!!         if((j.eq.1) .and. (k.eq.1))
!!!  &      write(*,*) 'CHECK6', d(1,1,1), e(1,1,1)

#ifdef SWITCH_TO_NR_SCHEME
!           return itmask
            do i = is+1, ie+1
               itmask(i) = itmask_tmp(i)
            enddo
!!!         write(*,*) 'j', de, HI, HII
#endif

!           Add the timestep to the elapsed time for each cell and find
!            minimum elapsed time step in this row

            ttmin = huge8
            do i = is+1, ie+1
               ttot(i) = min(ttot(i) + dtit(i), dt)
               if (abs(dt-ttot(i)) .lt.
     &              tolerance*dt) itmask(i) = .false.
               if (ttot(i).lt.ttmin) ttmin = ttot(i)
            enddo
!!!         write(*,*) 'k', de, HI, HII

!           If all cells are done (on this slice), then exit

            if (abs(dt-ttmin) .lt. tolerance*dt) go to 9999

!           Next subcycle iteration

         enddo ! end of 'do iter = 1, itmax'

 9999    continue
!!          write(*,*) 'l', de, HI, HII
!        do i = is+1, ie+1
!!!         if(tgas(i) .lt. 1.d1) then
!!!           write(*,*) 'LOWTEMP'
!           if(d(i,j,k)*dom .gt. 5.d5) then
!             write(*,*) 'HIGHDENSa'
!!!         if((d(i,j,k)*dom .gt. 5.d5)
!!!  &   .and. (tgas(i) .lt. 2.d2) ) then
!!!           write(*,*) 'STRANGE'
!           if( ((j-1) + (k-1)*jn)*in+i .eq. 1410 ) then
!             write(*,*) 'HIGHDENSa'
!    &         ,((j-1) + (k-1)*jn)*in+i
!    &         , d(i,j,k)*dom, e(i,j,k)
!    &         , kphHI(i,j,k)
!    &         , kdissH2I(i,j,k)
!           endif
!        enddo

!       Abort if iteration count exceeds maximum

         if (iter .gt. itmax) then
#ifdef _OPENMP
!$omp critical
#endif
	    write(0,*) 'inside if statement solve rate cool:',is,ie
            write(6,*) 'MULTI_COOL iter > ',itmax,' at j,k =',j,k
            write(0,*) 'FATAL error (2) in MULTI_COOL'
            write(0,'(" dt = ",1pe10.3," ttmin = ",1pe10.3)') dt, ttmin
            write(0,'((16(1pe8.1)))') (dtit(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (ttot(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (edot(i),i=is+1,ie+1)
            write(0,'((16(l3)))') (itmask(i),i=is+1,ie+1)
#ifdef _OPENMP
!$omp end critical
#endif
c            WARNING_MESSAGE
         endif

         if (iter .gt. itmax/2) then
#ifdef _OPENMP
!$omp critical
#endif
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
#ifdef _OPENMP
!$omp end critical
#endif
         end if
!     
!     Next j,k
!     
      enddo ! end of 'do t = 0, dk*dj-1'
#ifdef _OPENMP
!$omp end parallel do
#endif

!     Convert densities back to comoving from proper

      if (iexpand .eq. 1) then

         call scale_fields_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal, dust,
#ifdef GRACKLE_MD
     &                  imchem, igrgr,
     &                  DM, HDII, HeHII,
     &                  CI, CII, CO, CO2,
     &                  OI, OH, H2O, O2,
     &                  SiI, SiOI, SiO2I,
     &                  CH, CH2, COII, OII,
     &                  OHII, H2OII, H3OII, O2II,
     &                  Mg, Al, S, Fe,
     &                  SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &                  AC, SiO2D, MgO, FeS, Al2O3,
     &                  impop3,
     &                  metal_loc, metal_C30, metal_F13,
#endif
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, idustfield,
     &                  aye**3)

      endif

      if (ispecies .gt. 0) then

!     Correct the species to ensure consistency (i.e. type conservation)

!!!         write(*,*) 'CHECK7', d(1,1,1), e(1,1,1)
#ifdef ABUNDANCE_CORRECTION
      call make_consistent_g(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, 
     &                     d, is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, fh, dtoh
#ifdef GRACKLE_MD
     &                   , imchem, igrgr, dom
     &                   , DM, HDII, HeHII
     &                   , CI, CII, CO, CO2
     &                   , OI, OH, H2O, O2
     &                   , SiI, SiOI, SiO2I
     &                   , CH, CH2, COII, OII
     &                   , OHII, H2OII, H3OII, O2II
     &                   , Mg, Al, S, Fe
     &                   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                   , AC, SiO2D, MgO, FeS, Al2O3
     &                   , impop3
     &                   , metal_loc, metal_C30, metal_F13
     &                   , loc_XC, loc_XO, loc_XMg, loc_XAl, loc_XSi
     &                   , loc_XS, loc_XFe
     &                   , loc_fC, loc_fO, loc_fMg, loc_fAl, loc_fSi
     &                   , loc_fS, loc_fFe
     &                   , C30_XC, C30_XO, C30_XMg, C30_XAl, C30_XSi
     &                   , C30_XS, C30_XFe
     &                   , C30_fC, C30_fO, C30_fMg, C30_fAl, C30_fSi
     &                   , C30_fS, C30_fFe
     &                   , F13_XC, F13_XO, F13_XMg, F13_XAl, F13_XSi
     &                   , F13_XS, F13_XFe
     &                   , F13_fC, F13_fO, F13_fMg, F13_fAl, F13_fSi
     &                   , F13_fS, F13_fFe
#endif
     &                      )
#endif

      endif
!!!         write(*,*) 'CHECK8', d(1,1,1), e(1,1,1)
            flush(6)
      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal, dust,
#ifdef GRACKLE_MD
     &                        imchem, igrgr,
     &                        DM, HDII, HeHII,
     &                        CI, CII, CO, CO2,
     &                        OI, OH, H2O, O2,
     &                        SiI, SiOI, SiO2I,
     &                        CH, CH2, COII, OII,
     &                        OHII, H2OII, H3OII, O2II,
     &                        Mg, Al, S, Fe,
     &                        SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &                        AC, SiO2D, MgO, FeS, Al2O3,
     &                        impop3,
     &                        metal_loc, metal_C30, metal_F13,
#endif
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, idustfield,
     &                        factor)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal,
     &        idustfield
      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        d(in,jn,kn),   metal(in,jn,kn), dust(in,jn,kn)
#ifdef GRACKLE_MD
       integer imchem, igrgr
       R_PREC DM(in,jn,kn), HDII(in,jn,kn), HeHII(in,jn,kn)
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
      integer impop3
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
#endif
      real*8 factor

!     locals

      integer i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*factor
            enddo
         enddo
      enddo

      if (ispecies .gt. 0) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  de(i,j,k)    = de(i,j,k)*factor
                  HI(i,j,k)    = HI(i,j,k)*factor
                  HII(i,j,k)   = HII(i,j,k)*factor
                  HeI(i,j,k)   = HeI(i,j,k)*factor
                  HeII(i,j,k)  = HeII(i,j,k)*factor
                  HeIII(i,j,k) = HeIII(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
#ifdef GRACKLE_MD
      if (ispecies .gt. 3) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DM(i,j,k)      = DM(i,j,k)*factor
                  HDII(i,j,k)    = HDII(i,j,k)*factor
                  HeHII(i,j,k)   = HeHII(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
#endif
      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*factor
                     metal_loc(i,j,k) = metal_loc(i,j,k)*factor
                  if (impop3 .gt. 0) then
                     metal_C30(i,j,k) = metal_C30(i,j,k)*factor
                     metal_F13(i,j,k) = metal_F13(i,j,k)*factor
                  endif
#ifdef GRACKLE_MD
                  if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
                  if (imchem .eq. 1) then
                     CI(i,j,k)      = CI(i,j,k)*factor
                     CII(i,j,k)     = CII(i,j,k)*factor
                     CO(i,j,k)      = CO(i,j,k)*factor
                     CO2(i,j,k)     = CO2(i,j,k)*factor
                     OI(i,j,k)      = OI(i,j,k)*factor
                     OH(i,j,k)      = OH(i,j,k)*factor
                     H2O(i,j,k)     = H2O(i,j,k)*factor
                     O2(i,j,k)      = O2(i,j,k)*factor
                     SiI(i,j,k)     = SiI(i,j,k)*factor
                     SiOI(i,j,k)    = SiOI(i,j,k)*factor
                     SiO2I(i,j,k)   = SiO2I(i,j,k)*factor
                     CH(i,j,k)      = CH(i,j,k)*factor
                     CH2(i,j,k)     = CH2(i,j,k)*factor
                     COII(i,j,k)    = COII(i,j,k)*factor
                     OII(i,j,k)     = OII(i,j,k)*factor
                     OHII(i,j,k)    = OHII(i,j,k)*factor
                     H2OII(i,j,k)   = H2OII(i,j,k)*factor
                     H3OII(i,j,k)   = H3OII(i,j,k)*factor
                     O2II(i,j,k)    = O2II(i,j,k)*factor
                  endif
                  if (igrgr .eq. 1) then
                     Mg(i,j,k)      = Mg(i,j,k)*factor
                     Al(i,j,k)      = Al(i,j,k)*factor
                     S(i,j,k)       = S(i,j,k)*factor
                     Fe(i,j,k)      = Fe(i,j,k)*factor
                     SiM(i,j,k)     = SiM(i,j,k)*factor
                     FeM(i,j,k)     = FeM(i,j,k)*factor
                     Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*factor
                     MgSiO3(i,j,k)  = MgSiO3(i,j,k)*factor
                     Fe3O4(i,j,k)   = Fe3O4(i,j,k)*factor
                     AC(i,j,k)      = AC(i,j,k)*factor
                     SiO2D(i,j,k)   = SiO2D(i,j,k)*factor
                     MgO(i,j,k)     = MgO(i,j,k)*factor
                     FeS(i,j,k)     = FeS(i,j,k)*factor
                     Al2O3(i,j,k)   = Al2O3(i,j,k)*factor
                  endif
                  endif
#endif
               enddo
            enddo
         enddo
      endif
      if (idustfield .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  dust(i,j,k) = dust(i,j,k)*factor
               enddo
            enddo
         enddo
      endif

      return
      end

c -----------------------------------------------------------
!   This routine ensures that the species aren't below tiny.

      subroutine ceiling_species_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                           HM, H2I, H2II, DI, DII, HDI, metal,
#ifdef GRACKLE_MD
     &                           imchem, igrgr,
     &                           DM, HDII, HeHII,
     &                           CI, CII, CO, CO2,
     &                           OI, OH, H2O, O2,
     &                           SiI, SiOI, SiO2I,
     &                           CH, CH2, COII, OII,
     &                           OHII, H2OII, H3OII, O2II,
     &                           Mg, Al, S, Fe,
     &                           SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4,
     &                           AC, SiO2D, MgO, FeS, Al2O3,
     &                           impop3,
     &                           metal_loc, metal_C30, metal_F13,
#endif
     &                           is, ie, js, je, ks, ke,
     &                           in, jn, kn, ispecies, imetal)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC  d(in,jn,kn),
     &        de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
#ifdef GRACKLE_MD
       integer imchem, igrgr
       R_PREC DM(in,jn,kn), HDII(in,jn,kn), HeHII(in,jn,kn)
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
      integer impop3
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
#endif

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
#ifdef GRACKLE_MD
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
#endif
      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = min(max(metal(i,j,k), tiny),
     &                 0.9_RKIND*d(i,j,k))
                     metal_loc(i,j,k) = min(max(metal_loc(i,j,k), tiny),
     &                 0.9_RKIND*d(i,j,k))
                  if( impop3 .gt. 0 ) then
                     metal_C30(i,j,k) = min(max(metal_C30(i,j,k), tiny),
     &                 0.9_RKIND*d(i,j,k))
                     metal_F13(i,j,k) = min(max(metal_F13(i,j,k), tiny),
     &                 0.9_RKIND*d(i,j,k))
                  endif
#ifdef GRACKLE_MD
                  if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
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
                  endif
                  if (igrgr .eq. 1) then
                     Mg(i,j,k)      = max(Mg(i,j,k), tiny)
                     Al(i,j,k)      = max(Al(i,j,k), tiny)
                     S(i,j,k)       = max(S(i,j,k), tiny)
                     Fe(i,j,k)      = max(Fe(i,j,k), tiny)
                     SiM(i,j,k)     = max(SiM(i,j,k), tiny)
                     FeM(i,j,k)     = max(FeM(i,j,k), tiny)
                     Mg2SiO4(i,j,k) = max(Mg2SiO4(i,j,k), tiny)
                     MgSiO3(i,j,k)  = max(MgSiO3(i,j,k), tiny)
                     Fe3O4(i,j,k)   = max(Fe3O4(i,j,k), tiny)
                     AC(i,j,k)      = max(AC(i,j,k), tiny)
                     SiO2D(i,j,k)   = max(SiO2D(i,j,k), tiny)
                     MgO(i,j,k)     = max(MgO(i,j,k), tiny)
                     FeS(i,j,k)     = max(FeS(i,j,k), tiny)
                     Al2O3(i,j,k)   = max(Al2O3(i,j,k), tiny)
                  endif
                  endif
#endif
               enddo
            enddo
         enddo
      endif

      return
      end



! -----------------------------------------------------------
! This routine uses the temperature to look up the chemical
!   rates which are tabulated in a log table as a function
!   of temperature.

      subroutine lookup_cool_rates1d_g(temstart, temend, nratec, j, k,
     &                is, ie, ithreebody, in, jn, kn,
     &                ispecies, anydust, iH2shield, iradshield,
     &                tgas1d, mmw, d, HI, HII, HeI, HeII, HeIII,
     &                HM, H2I, H2II, DI, DII, HDI,
     &                tdust, dust2gas,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a, k57a,
     &                k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &                ncrna, ncrd1a, ncrd2a,
     &                avgsighi, avgsighei, avgsigheii, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26, k28, k29, k30, k31,
     &                k50, k51, k52, k53, k54, k55, k56, k57,
     &                k58, k13dd, k24shield, k25shield, k26shield,
     &                k28shield, k29shield, k30shield, k31shield,
     &                h2dust, ncrn, ncrd1, ncrd2,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, xbase1, dx_cgs, c_ljeans,
     &                iradtrans, kdissH2I, xH2shield, itmask
#ifdef GRACKLE_MD
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
     &              , HC30_SiM, HC30_FeM, HC30_Mg2SiO4, HC30_MgSiO3
     &              , HC30_Fe3O4, HC30_AC, HC30_SiO2D, HC30_MgO
     &              , HC30_FeS, HC30_Al2O3
     &              , KC30_SiM, KC30_FeM, KC30_Mg2SiO4, KC30_MgSiO3
     &              , KC30_Fe3O4, KC30_AC, KC30_SiO2D, KC30_MgO
     &              , KC30_FeS, KC30_Al2O3
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
     &              , HF13_SiM, HF13_FeM, HF13_Mg2SiO4, HF13_MgSiO3
     &              , HF13_Fe3O4, HF13_AC, HF13_SiO2D, HF13_MgO
     &              , HF13_FeS, HF13_Al2O3
     &              , KF13_SiM, KF13_FeM, KF13_Mg2SiO4, KF13_MgSiO3
     &              , KF13_Fe3O4, KF13_AC, KF13_SiO2D, KF13_MgO
     &              , KF13_FeS, KF13_Al2O3
     &              , kd0, kd1, kd2, kd3, kd4, kd5, kd6, kd7, kd8, kd9
#endif
     &      )
! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer is, ie, nratec,
     &        in, jn, kn, ispecies, ithreebody, j, k,
     &        ndratec, iH2shield, iradshield, iradtrans
      real*8 temstart, temend, tgas1d(in), mmw(in), dom,
     &       dtemstart, dtemend
      real*8 coolunit, tbase1, xbase1, dx_cgs, c_ljeans
      logical itmask(in), anydust

!     Chemistry rates as a function of temperature

      real*8 k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &       k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &       k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &       k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &       k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &       k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &       k54a(nratec), k55a(nratec), k56a(nratec), k57a(nratec),
     &       k58a(nratec), k13dda(nratec, 14), h2dusta(nratec, ndratec),
     &       ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &       k24, k25, k26, k28, k29, k30, k31,
     &       piHI, piHeI,
     &       avgsighi, avgsighei, avgsigheii
#ifdef GRACKLE_MD
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
#endif

!     Density fields

      R_PREC  d(in,jn,kn), HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn), H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn), DII(in,jn,kn), HDI(in,jn,kn)	
#ifdef GRACKLE_MD
       R_PREC metal(in,jn,kn)
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
       R_PREC DM(in,jn,kn)   , HDII(in,jn,kn)   , HeHII(in,jn,kn)
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
       integer imetal, imchem, igrgr
#endif

!     Radiation fields

      R_PREC kdissH2I(in,jn,kn)

!     H2 self-shielding length-scale field

      R_PREC  xH2shield(in,jn,kn)

!     Returned rate values

      real*8 k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in),
     &       k13dd(in, 14), h2dust(in), 
     &       ncrn(in), ncrd1(in), ncrd2(in),
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k28shield(in), k29shield(in), k30shield(in),
     &       k31shield(in)
#ifdef GRACKLE_MD
      real*8 k125(in), k129(in), k130(in), k131(in), k132(in)
     &     , k133(in), k134(in), k135(in), k136(in), k137(in)
     &     , k148(in), k149(in), k150(in), k151(in), k152(in)
     &     , k153(in)
     &     , kz15(in), kz16(in), kz17(in), kz18(in), kz19(in)
     &     , kz20(in), kz21(in), kz22(in), kz23(in), kz24(in)
     &     , kz25(in), kz26(in), kz27(in), kz28(in), kz29(in)
     &     , kz30(in), kz31(in), kz32(in), kz33(in), kz34(in)
     &     , kz35(in), kz36(in), kz37(in), kz38(in), kz39(in)
     &     , kz40(in), kz41(in), kz42(in), kz43(in), kz44(in)
     &     , kz45(in), kz46(in), kz47(in), kz48(in), kz49(in)
     &     , kz50(in), kz51(in), kz52(in), kz53(in), kz54(in)
#endif
!     1D temporaries (passed in)

      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in),
     &       tdust(in), dust2gas(in)

!     1D temporaries (not passed in)

      integer*8 d_indixe(in)
      real*8 d_t1(in), d_t2(in), d_logtem(in), d_tdef(in),
     &       dusti1(in), dusti2(in), divrhoa(6),
     &       f_shield_H(in), f_shield_He(in)

!     Parameters

      real*8 everg, e24, e26
      parameter(everg = ev2erg, e24 = 13.6_DKIND,
     &     e26 = 24.6_DKIND)

!     locals

      integer i, n1
      real*8 factor, x, logtem0, logtem9, dlogtem, nh,
     &       d_logtem0, d_logtem9, d_dlogtem, divrho, N_H2,
     &       f_shield, b_doppler, l_H2shield
      real*8 k13_CID, k13_DT
#ifdef GRACKLE_MD
      real*8 k13ind
#endif
#ifdef GRACKLE_MD
      real*8 logT(in), logrho(in)
       integer grain_N(2), grain_Size
       real*8  grain_D(grain_N(1)), grain_T(grain_N(2))
     &       , grain_dD, grain_dT, Hgrain(grain_Size)
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
      real*8  HC30_SiM(C30_Size), HC30_FeM(C30_Size)
     &      , HC30_Mg2SiO4(C30_Size), HC30_MgSiO3(C30_Size)
     &      , HC30_Fe3O4(C30_Size), HC30_AC(C30_Size)
     &      , HC30_SiO2D(C30_Size), HC30_MgO(C30_Size)
     &      , HC30_FeS(C30_Size), HC30_Al2O3(C30_Size)
      real*8  KC30_SiM(C30_Size), KC30_FeM(C30_Size)
     &      , KC30_Mg2SiO4(C30_Size), KC30_MgSiO3(C30_Size)
     &      , KC30_Fe3O4(C30_Size), KC30_AC(C30_Size)
     &      , KC30_SiO2D(C30_Size), KC30_MgO(C30_Size)
     &      , KC30_FeS(C30_Size), KC30_Al2O3(C30_Size)
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
      real*8  HF13_SiM(F13_Size), HF13_FeM(F13_Size)
     &      , HF13_Mg2SiO4(F13_Size), HF13_MgSiO3(F13_Size)
     &      , HF13_Fe3O4(F13_Size), HF13_AC(F13_Size)
     &      , HF13_SiO2D(F13_Size), HF13_MgO(F13_Size)
     &      , HF13_FeS(F13_Size), HF13_Al2O3(F13_Size)
      real*8  KF13_SiM(F13_Size), KF13_FeM(F13_Size)
     &      , KF13_Mg2SiO4(F13_Size), KF13_MgSiO3(F13_Size)
     &      , KF13_Fe3O4(F13_Size), KF13_AC(F13_Size)
     &      , KF13_SiO2D(F13_Size), KF13_MgO(F13_Size)
     &      , KF13_FeS(F13_Size), KF13_Al2O3(F13_Size)
      real*8  log_kh2, log_kgg
!     grain growth
      real*8  FC30_SiM(in), FC30_FeM(in), FC30_Mg2SiO4(in)
     &      , FC30_MgSiO3(in), FC30_Fe3O4(in), FC30_AC(in)
     &      , FC30_SiO2D(in), FC30_MgO(in), FC30_FeS(in)
     &      , FC30_Al2O3(in)
      real*8  FF13_SiM(in), FF13_FeM(in), FF13_Mg2SiO4(in)
     &      , FF13_MgSiO3(in), FF13_Fe3O4(in), FF13_AC(in)
     &      , FF13_SiO2D(in), FF13_MgO(in), FF13_FeS(in)
     &      , FF13_Al2O3(in)
      real*8  kd0(in), kd1(in), kd2(in), kd3(in), kd4(in)
     &      , kd5(in), kd6(in), kd7(in), kd8(in), kd9(in)
      real*8 fh, mh
      parameter (mh = mass_h)      !DPC
      real*8 HC30(in)
     &     , kd0_C30(in), kd1_C30(in), kd2_C30(in), kd3_C30(in)
     &     , kd4_C30(in), kd5_C30(in), kd6_C30(in), kd7_C30(in)
     &     , kd8_C30(in), kd9_C30(in)
      real*8 HF13(in)
     &     , kd0_F13(in), kd1_F13(in), kd2_F13(in), kd3_F13(in)
     &     , kd4_F13(in), kd5_F13(in), kd6_F13(in), kd7_F13(in)
     &     , kd8_F13(in), kd9_F13(in)
#endif

      real*8 nSSh, nratio
!!    if(itmask(1))
!!   &write(*,*) 'b', 0.d0, HI, HII, HeI, HeII, HeIII,
!!   &           HM, H2I, H2II, DI, DII, HDI

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/real(nratec-1, DKIND)

      do i = is+1, ie+1
         if (itmask(i)) then
!        Compute temp-centered temperature (and log)

!        logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
         logtem(i) = log(tgas1d(i))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Find index into tble and precompute interpolation values

         indixe(i) = min(nratec-1,
     &        max(1,int((logtem(i)-logtem0)/dlogtem, DIKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!        Do linear table lookup (in log temperature)

         k1(i) = k1a(indixe(i)) +
     &           (k1a(indixe(i)+1) -k1a(indixe(i)))*tdef(i)
         k2(i) = k2a(indixe(i)) +
     &           (k2a(indixe(i)+1) -k2a(indixe(i)))*tdef(i)
         k3(i) = k3a(indixe(i)) +
     &           (k3a(indixe(i)+1) -k3a(indixe(i)))*tdef(i)
         k4(i) = k4a(indixe(i)) +
     &           (k4a(indixe(i)+1) -k4a(indixe(i)))*tdef(i)
         k5(i) = k5a(indixe(i)) +
     &           (k5a(indixe(i)+1) -k5a(indixe(i)))*tdef(i)
         k6(i) = k6a(indixe(i)) +
     &           (k6a(indixe(i)+1) -k6a(indixe(i)))*tdef(i)
         k57(i) = k57a(indixe(i)) +
     &            (k57a(indixe(i)+1) -k57a(indixe(i)))*tdef(i)
         k58(i) = k58a(indixe(i)) +
     &            (k58a(indixe(i)+1) -k58a(indixe(i)))*tdef(i)
      endif
      enddo

!     Look-up for 9-species model

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k7(i) = k7a(indixe(i)) +
     &            (k7a(indixe(i)+1) -k7a(indixe(i)))*tdef(i)
            k8(i) = k8a(indixe(i)) +
     &            (k8a(indixe(i)+1) -k8a(indixe(i)))*tdef(i)
            k9(i) = k9a(indixe(i)) +
     &            (k9a(indixe(i)+1) -k9a(indixe(i)))*tdef(i)
            k10(i) = k10a(indixe(i)) +
     &            (k10a(indixe(i)+1) -k10a(indixe(i)))*tdef(i)
            k11(i) = k11a(indixe(i)) +
     &            (k11a(indixe(i)+1) -k11a(indixe(i)))*tdef(i)
            k12(i) = k12a(indixe(i)) +
     &            (k12a(indixe(i)+1) -k12a(indixe(i)))*tdef(i)
            k13(i) = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*tdef(i)
            k14(i) = k14a(indixe(i)) +
     &            (k14a(indixe(i)+1) -k14a(indixe(i)))*tdef(i)
            k15(i) = k15a(indixe(i)) +
     &            (k15a(indixe(i)+1) -k15a(indixe(i)))*tdef(i)
            k16(i) = k16a(indixe(i)) +
     &            (k16a(indixe(i)+1) -k16a(indixe(i)))*tdef(i)
            k17(i) = k17a(indixe(i)) +
     &            (k17a(indixe(i)+1) -k17a(indixe(i)))*tdef(i)
            k18(i) = k18a(indixe(i)) +
     &            (k18a(indixe(i)+1) -k18a(indixe(i)))*tdef(i)
            k19(i) = k19a(indixe(i)) +
     &            (k19a(indixe(i)+1) -k19a(indixe(i)))*tdef(i)
            k22(i) = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*tdef(i)

!     H2 formation heating terms.

            ncrn(i) = ncrna(indixe(i)) +
     &           (ncrna(indixe(i)+1) -ncrna(indixe(i)))*tdef(i)
            ncrd1(i) = ncrd1a(indixe(i)) +
     &           (ncrd1a(indixe(i)+1) -ncrd1a(indixe(i)))*tdef(i)
            ncrd2(i) = ncrd2a(indixe(i)) +
     &           (ncrd2a(indixe(i)+1) -ncrd2a(indixe(i)))*tdef(i)

         endif
         enddo

         do n1 = 1, 14
            do i = is+1, ie+1
               if (itmask(i)) then
               k13dd(i,n1) = k13dda(indixe(i),n1) +
     &             (k13dda(indixe(i)+1,n1) - 
     &               k13dda(indixe(i)  ,n1) )*tdef(i)
            endif
            enddo
         enddo         

      endif

!     Look-up for 12-species model

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k50(i) = k50a(indixe(i)) +
     &            (k50a(indixe(i)+1) -k50a(indixe(i)))*tdef(i)
            k51(i) = k51a(indixe(i)) +
     &            (k51a(indixe(i)+1) -k51a(indixe(i)))*tdef(i)
            k52(i) = k52a(indixe(i)) +
     &            (k52a(indixe(i)+1) -k52a(indixe(i)))*tdef(i)
            k53(i) = k53a(indixe(i)) +
     &            (k53a(indixe(i)+1) -k53a(indixe(i)))*tdef(i)
            k54(i) = k54a(indixe(i)) +
     &            (k54a(indixe(i)+1) -k54a(indixe(i)))*tdef(i)
            k55(i) = k55a(indixe(i)) +
     &            (k55a(indixe(i)+1) -k55a(indixe(i)))*tdef(i)
            k56(i) = k56a(indixe(i)) +
     &            (k56a(indixe(i)+1) -k56a(indixe(i)))*tdef(i)
         endif
         enddo
      endif
#ifdef GRACKLE_MD
!     Look-up for 15-species model

      if (ispecies .gt. 3) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k125(i) = k125a(indixe(i)) +
     &            (k125a(indixe(i)+1) -k125a(indixe(i)))*tdef(i)
            k129(i) = k129a(indixe(i)) +
     &            (k129a(indixe(i)+1) -k129a(indixe(i)))*tdef(i)
            k130(i) = k130a(indixe(i)) +
     &            (k130a(indixe(i)+1) -k130a(indixe(i)))*tdef(i)
            k131(i) = k131a(indixe(i)) +
     &            (k131a(indixe(i)+1) -k131a(indixe(i)))*tdef(i)
            k132(i) = k132a(indixe(i)) +
     &            (k132a(indixe(i)+1) -k132a(indixe(i)))*tdef(i)
            k133(i) = k133a(indixe(i)) +
     &            (k133a(indixe(i)+1) -k133a(indixe(i)))*tdef(i)
            k134(i) = k134a(indixe(i)) +
     &            (k134a(indixe(i)+1) -k134a(indixe(i)))*tdef(i)
            k135(i) = k135a(indixe(i)) +
     &            (k135a(indixe(i)+1) -k135a(indixe(i)))*tdef(i)
            k136(i) = k136a(indixe(i)) +
     &            (k136a(indixe(i)+1) -k136a(indixe(i)))*tdef(i)
            k137(i) = k137a(indixe(i)) +
     &            (k137a(indixe(i)+1) -k137a(indixe(i)))*tdef(i)
            k148(i) = k148a(indixe(i)) +
     &            (k148a(indixe(i)+1) -k148a(indixe(i)))*tdef(i)
            k149(i) = k149a(indixe(i)) +
     &            (k149a(indixe(i)+1) -k149a(indixe(i)))*tdef(i)
            k150(i) = k150a(indixe(i)) +
     &            (k150a(indixe(i)+1) -k150a(indixe(i)))*tdef(i)
            k151(i) = k151a(indixe(i)) +
     &            (k151a(indixe(i)+1) -k151a(indixe(i)))*tdef(i)
            k152(i) = k152a(indixe(i)) +
     &            (k152a(indixe(i)+1) -k152a(indixe(i)))*tdef(i)
            k153(i) = k153a(indixe(i)) +
     &            (k153a(indixe(i)+1) -k153a(indixe(i)))*tdef(i)
         endif
         enddo
      endif

!     Look-up for metal species model

      if (imchem .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            kz15(i) = kz15a(indixe(i)) +
     &            (kz15a(indixe(i)+1) -kz15a(indixe(i)))*tdef(i)
            kz16(i) = kz16a(indixe(i)) +
     &            (kz16a(indixe(i)+1) -kz16a(indixe(i)))*tdef(i)
            kz17(i) = kz17a(indixe(i)) +
     &            (kz17a(indixe(i)+1) -kz17a(indixe(i)))*tdef(i)
            kz18(i) = kz18a(indixe(i)) +
     &            (kz18a(indixe(i)+1) -kz18a(indixe(i)))*tdef(i)
            kz19(i) = kz19a(indixe(i)) +
     &            (kz19a(indixe(i)+1) -kz19a(indixe(i)))*tdef(i)
            kz20(i) = kz20a(indixe(i)) +
     &            (kz20a(indixe(i)+1) -kz20a(indixe(i)))*tdef(i)
            kz21(i) = kz21a(indixe(i)) +
     &            (kz21a(indixe(i)+1) -kz21a(indixe(i)))*tdef(i)
            kz22(i) = kz22a(indixe(i)) +
     &            (kz22a(indixe(i)+1) -kz22a(indixe(i)))*tdef(i)
            kz23(i) = kz23a(indixe(i)) +
     &            (kz23a(indixe(i)+1) -kz23a(indixe(i)))*tdef(i)
            kz24(i) = kz24a(indixe(i)) +
     &            (kz24a(indixe(i)+1) -kz24a(indixe(i)))*tdef(i)
            kz25(i) = kz25a(indixe(i)) +
     &            (kz25a(indixe(i)+1) -kz25a(indixe(i)))*tdef(i)
            kz26(i) = kz26a(indixe(i)) +
     &            (kz26a(indixe(i)+1) -kz26a(indixe(i)))*tdef(i)
            kz27(i) = kz27a(indixe(i)) +
     &            (kz27a(indixe(i)+1) -kz27a(indixe(i)))*tdef(i)
            kz28(i) = kz28a(indixe(i)) +
     &            (kz28a(indixe(i)+1) -kz28a(indixe(i)))*tdef(i)
            kz29(i) = kz29a(indixe(i)) +
     &            (kz29a(indixe(i)+1) -kz29a(indixe(i)))*tdef(i)
            kz30(i) = kz30a(indixe(i)) +
     &            (kz30a(indixe(i)+1) -kz30a(indixe(i)))*tdef(i)
            kz31(i) = kz31a(indixe(i)) +
     &            (kz31a(indixe(i)+1) -kz31a(indixe(i)))*tdef(i)
            kz32(i) = kz32a(indixe(i)) +
     &            (kz32a(indixe(i)+1) -kz32a(indixe(i)))*tdef(i)
            kz33(i) = kz33a(indixe(i)) +
     &            (kz33a(indixe(i)+1) -kz33a(indixe(i)))*tdef(i)
            kz34(i) = kz34a(indixe(i)) +
     &            (kz34a(indixe(i)+1) -kz34a(indixe(i)))*tdef(i)
            kz35(i) = kz35a(indixe(i)) +
     &            (kz35a(indixe(i)+1) -kz35a(indixe(i)))*tdef(i)
            kz36(i) = kz36a(indixe(i)) +
     &            (kz36a(indixe(i)+1) -kz36a(indixe(i)))*tdef(i)
            kz37(i) = kz37a(indixe(i)) +
     &            (kz37a(indixe(i)+1) -kz37a(indixe(i)))*tdef(i)
            kz38(i) = kz38a(indixe(i)) +
     &            (kz38a(indixe(i)+1) -kz38a(indixe(i)))*tdef(i)
            kz39(i) = kz39a(indixe(i)) +
     &            (kz39a(indixe(i)+1) -kz39a(indixe(i)))*tdef(i)
            kz40(i) = kz40a(indixe(i)) +
     &            (kz40a(indixe(i)+1) -kz40a(indixe(i)))*tdef(i)
            kz41(i) = kz41a(indixe(i)) +
     &            (kz41a(indixe(i)+1) -kz41a(indixe(i)))*tdef(i)
            kz42(i) = kz42a(indixe(i)) +
     &            (kz42a(indixe(i)+1) -kz42a(indixe(i)))*tdef(i)
            kz43(i) = kz43a(indixe(i)) +
     &            (kz43a(indixe(i)+1) -kz43a(indixe(i)))*tdef(i)
            kz44(i) = kz44a(indixe(i)) +
     &            (kz44a(indixe(i)+1) -kz44a(indixe(i)))*tdef(i)
            kz45(i) = kz45a(indixe(i)) +
     &            (kz45a(indixe(i)+1) -kz45a(indixe(i)))*tdef(i)
            kz46(i) = kz46a(indixe(i)) +
     &            (kz46a(indixe(i)+1) -kz46a(indixe(i)))*tdef(i)
            kz47(i) = kz47a(indixe(i)) +
     &            (kz47a(indixe(i)+1) -kz47a(indixe(i)))*tdef(i)
            kz48(i) = kz48a(indixe(i)) +
     &            (kz48a(indixe(i)+1) -kz48a(indixe(i)))*tdef(i)
            kz49(i) = kz49a(indixe(i)) +
     &            (kz49a(indixe(i)+1) -kz49a(indixe(i)))*tdef(i)
            kz50(i) = kz50a(indixe(i)) +
     &            (kz50a(indixe(i)+1) -kz50a(indixe(i)))*tdef(i)
            kz51(i) = kz51a(indixe(i)) +
     &            (kz51a(indixe(i)+1) -kz51a(indixe(i)))*tdef(i)
            kz52(i) = kz52a(indixe(i)) +
     &            (kz52a(indixe(i)+1) -kz52a(indixe(i)))*tdef(i)
            kz53(i) = kz53a(indixe(i)) +
     &            (kz53a(indixe(i)+1) -kz53a(indixe(i)))*tdef(i)
            kz54(i) = kz54a(indixe(i)) +
     &            (kz54a(indixe(i)+1) -kz54a(indixe(i)))*tdef(i)
#define CHEM_TEMPO
#ifdef CHEM_TEMPO
!!!!!!!!!! He-baring
!!!!!!!!!!  k3(i) = 0.d0
!!!!!!!!!!  k4(i) = 0.d0
!!!!!!!!!!  k5(i) = 0.d0
!!!!!!!!!!  k6(i) = 0.d0
!!!!!!!!!!  k58(i) = 0.d0
!!!!!!!!!!  k148(i) = 0.d0
!!!!!!!!!!  k149(i) = 0.d0
!!!!!!!!!!  k150(i) = 0.d0
!!!!!!!!!!  k151(i) = 0.d0
!!!!!!!!!!  k152(i) = 0.d0
!!!!!!!!!!  k153(i) = 0.d0

!!!!!!!!!! C-baring
!!!!        kz15(i) = 0.d0
!!!!        kz16(i) = 0.d0
!!!!        kz20(i) = 0.d0
!!!!        kz23(i) = 0.d0
!!          kz26(i) = 0.d0
!!!!        kz27(i) = 0.d0
!!          kz28(i) = 0.d0
!!          kz29(i) = 0.d0
!!          kz32(i) = 0.d0
!!          kz37(i) = 0.d0
!!          kz38(i) = 0.d0!!!!!
!!          kz43(i) = 0.d0
!!!!        kz44(i) = 0.d0
!!!!        kz51(i) = 0.d0

!!!!!!!!!! H3OII
!!          kz42(i) = 0.d0
!!          kz48(i) = 0.d0
!!          kz49(i) = 0.d0
#endif
         endif
         enddo
      endif
#endif

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

!     Look-up for H2 formation on dust
      if (anydust) then
#ifdef GRACKLE_MD
         do i = is+1, ie+1
            if ( itmask(i) ) then
               logT(i)   = log10(tgas1d(i))
               logrho(i) = log10(d(i,j,k) * dom*mh)
            endif
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then

            call interpolate_2D_g(
     &        logrho(i), logT(i), int8(grain_N), grain_D, grain_dD, 
     &        grain_T, grain_dT, int8(grain_Size), 
     &        Hgrain, log_kh2)

              h2dust(i) = 1.d1**(log_kh2)/ dom
     &               * metal_loc(i,j,k)/(0.02d0*d(i,j,k))

            endif
         enddo

         if(impop3 .gt. 0) then

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
     &              , HC30_SiM, HC30_FeM, HC30_Mg2SiO4, HC30_MgSiO3
     &              , HC30_Fe3O4, HC30_AC, HC30_SiO2D, HC30_MgO
     &              , HC30_FeS, HC30_Al2O3
     &              , HC30
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
     &              , HF13_SiM, HF13_FeM, HF13_Mg2SiO4, HF13_MgSiO3
     &              , HF13_Fe3O4, HF13_AC, HF13_SiO2D, HF13_MgO
     &              , HF13_FeS, HF13_Al2O3
     &              , HF13
     &            )

            do i = is+1, ie+1
               if ( itmask(i) ) then
                  h2dust(i) = h2dust(i) +
     &              ( HC30(i) * metal_C30(i,j,k)
     &              + HF13(i) * metal_F13(i,j,k)
     &              )/(0.02d0*d(i,j,k)) / dom
               endif
            enddo

         endif !! impop3

#else
         d_logtem0 = log(dtemstart)
         d_logtem9 = log(dtemend)
         d_dlogtem = (log(dtemend) - log(dtemstart))/
     &        real(ndratec-1, DKIND)

         do i = is+1, ie+1
            if (itmask(i)) then

!              Assume dust melting at T > 1500 K

               if (tdust(i) .gt. dtemend) then
                  h2dust(i) = tiny8
               else

!                 Get log dust temperature

                  d_logtem(i) = log(tdust(i))
                  d_logtem(i) = max(d_logtem(i), d_logtem0)
                  d_logtem(i) = min(d_logtem(i), d_logtem9)

!                 Find index into table and precompute interpolation values

                  d_indixe(i) = min(ndratec-1,
     &                 max(1,
     &                 int((d_logtem(i)-d_logtem0)/d_dlogtem,
     &                 DIKIND)+1))
                  d_t1(i) = (d_logtem0 + (d_indixe(i) - 1)*d_dlogtem)
                  d_t2(i) = (d_logtem0 + (d_indixe(i)    )*d_dlogtem)
                  d_tdef(i) = (d_logtem(i) - d_t1(i)) / 
     &                 (d_t2(i) - d_t1(i))

!                 Get rate from 2D interpolation

                  dusti1(i) = h2dusta(indixe(i), d_indixe(i)) +
     &                 (h2dusta(indixe(i)+1, d_indixe(i)) - 
     &                 h2dusta(indixe(i),   d_indixe(i)))*tdef(i)
                  dusti2(i) = h2dusta(indixe(i), d_indixe(i)+1) +
     &                 (h2dusta(indixe(i)+1, d_indixe(i)+1) - 
     &                 h2dusta(indixe(i),   d_indixe(i)+1))*tdef(i)
                  h2dust(i) = dusti1(i) + 
     &                 (dusti2(i) - dusti1(i))*d_tdef(i)

!                 Multiply by dust to gas ratio

                  h2dust(i) = h2dust(i) * dust2gas(i)

               endif
            endif
         enddo
#endif

#ifdef GRACKLE_MD
!     Compute grain growth rate

         if(igrgr .eq. 1) then
         if(impop3 .gt. 0) then

            call lookup_grain_growth_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FC30_SiM, FC30_FeM, FC30_Mg2SiO4, FC30_MgSiO3
     &              , FC30_Fe3O4, FC30_AC, FC30_SiO2D, FC30_MgO
     &              , FC30_FeS, FC30_Al2O3
     &              , metal_C30
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
     &              , KC30_SiM, KC30_FeM, KC30_Mg2SiO4, KC30_MgSiO3
     &              , KC30_Fe3O4, KC30_AC, KC30_SiO2D, KC30_MgO
     &              , KC30_FeS, KC30_Al2O3
     &              , kd0_C30, kd1_C30, kd2_C30, kd3_C30, kd4_C30
     &              , kd5_C30, kd6_C30, kd7_C30, kd8_C30, kd9_C30
     &            )
            call lookup_grain_growth_rates1d(
     &                itmask
     &              , in, jn, kn, is, ie, j, k
     &              , logrho, logT
     &              , FF13_SiM, FF13_FeM, FF13_Mg2SiO4, FF13_MgSiO3
     &              , FF13_Fe3O4, FF13_AC, FF13_SiO2D, FF13_MgO
     &              , FF13_FeS, FF13_Al2O3
     &              , metal_F13
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
     &              , KF13_SiM, KF13_FeM, KF13_Mg2SiO4, KF13_MgSiO3
     &              , KF13_Fe3O4, KF13_AC, KF13_SiO2D, KF13_MgO
     &              , KF13_FeS, KF13_Al2O3
     &              , kd0_F13, kd1_F13, kd2_F13, kd3_F13, kd4_F13
     &              , kd5_F13, kd6_F13, kd7_F13, kd8_F13, kd9_F13
     &            )

            do i = is+1, ie+1
               if ( itmask(i) ) then

                  kd0(i) = (kd0_C30(i) * metal_C30(i,j,k)
     &                    + kd0_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * SiI(i,j,k) / 28._DKIND**1.5_DKIND

                  kd1(i) = (kd1_C30(i) * metal_C30(i,j,k)
     &                    + kd1_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * Fe(i,j,k) / 56._DKIND**1.5_DKIND

                  kd2(i) = (kd2_C30(i) * metal_C30(i,j,k)
     &                    + kd2_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * min( Mg  (i,j,k) / 24._DKIND**1.5_DKIND
     &                       / 2._DKIND
     &                    , SiOI(i,j,k) / 44._DKIND**1.5_DKIND
     &                    , H2O (i,j,k) / 18._DKIND**1.5_DKIND
     &                       / 3._DKIND )

                  kd3(i) = (kd3_C30(i) * metal_C30(i,j,k)
     &                    + kd3_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * min( Mg  (i,j,k) / 24._DKIND**1.5_DKIND
     &                    , SiOI(i,j,k) / 44._DKIND**1.5_DKIND
     &                    , H2O (i,j,k) / 18._DKIND**1.5_DKIND
     &                       / 2._DKIND )

                  kd4(i) = (kd4_C30(i) * metal_C30(i,j,k)
     &                    + kd4_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * min( Fe  (i,j,k) / 56._DKIND**1.5_DKIND
     &                       / 3._DKIND
     &                    , H2O (i,j,k) / 18._DKIND**1.5_DKIND
     &                       / 4._DKIND )

                  kd5(i) = (kd5_C30(i) * metal_C30(i,j,k)
     &                    + kd5_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * CI(i,j,k) / 12._DKIND**1.5_DKIND

                  kd6(i) = (kd6_C30(i) * metal_C30(i,j,k)
     &                    + kd6_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * SiO2I(i,j,k) / 60._DKIND**1.5_DKIND

                  kd7(i) = (kd7_C30(i) * metal_C30(i,j,k)
     &                    + kd7_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * min( Mg  (i,j,k) / 24._DKIND**1.5_DKIND
     &                    , H2O (i,j,k) / 18._DKIND**1.5_DKIND )

                  kd8(i) = (kd8_C30(i) * metal_C30(i,j,k)
     &                    + kd8_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * min( S (i,j,k) / 32._DKIND**1.5_DKIND
     &                    , Fe(i,j,k) / 56._DKIND**1.5_DKIND )

                  kd9(i) = (kd9_C30(i) * metal_C30(i,j,k)
     &                    + kd9_F13(i) * metal_F13(i,j,k)
     &               )/(0.02d0*d(i,j,k)) / dom
     &               * min( Al (i,j,k) / 27._DKIND**1.5_DKIND
     &                       / 2._DKIND
     &                    , H2O(i,j,k) / 18._DKIND**1.5_DKIND
     &                       / 3._DKIND )

               endif
            enddo

         endif !! impop3
         endif
#endif
      endif
!        Include approximate self-shielding factors if requested

      do i = is+1, ie+1
         if (itmask(i)) then
            k24shield(i) = k24
            k25shield(i) = k25
            k26shield(i) = k26
            k28shield(i) = k28
            k29shield(i) = k29
            k30shield(i) = k30
         endif
      enddo

!
!     H2 self-shielding (Sobolev-like, spherically averaged, Wolcott-Green+ 2011)
!

      if (ispecies .gt. 1) then
      if (iradtrans == 0) then
         do i = is+1, ie+1
            if (itmask(i)) then
               k31shield(i) = k31
            endif 
         enddo
      else
         do i = is+1, ie+1
            if (itmask(i)) then
               k31shield(i) = k31 + kdissH2I(i,j,k)
            endif
         enddo

      endif

      if (iH2shield .gt. 0) then
         do i = is+1, ie+1
            if (itmask(i)) then

!              Calculate a Sobolev-like length assuming a 3D grid.
               if (iH2shield == 1) then

               divrhoa(1) = d(i+1, j  , k  ) - d(i,j,k)
               divrhoa(2) = d(i-1, j  , k  ) - d(i,j,k)
               divrhoa(3) = d(i  , j+1, k  ) - d(i,j,k)
               divrhoa(4) = d(i  , j-1, k  ) - d(i,j,k)
               divrhoa(5) = d(i  , j  , k+1) - d(i,j,k)
               divrhoa(6) = d(i  , j  , k-1) - d(i,j,k)
               divrho = tiny
!              Exclude directions with (drho/ds > 0)
               do n1 = 1, 6
                  if (divrhoa(n1) .lt. 0._DKIND) then
                     divrho = divrho + divrhoa(n1)
                  endif
               enddo
!              (rho / divrho) is the Sobolev-like length in cell widths
               l_H2shield = min(dx_cgs * d(i,j,k) / abs(divrho), xbase1)

!              User-supplied length-scale field.
               else if (iH2shield == 2) then
                  l_H2shield = xH2shield(i,j,k) * xbase1

!              Jeans Length
               else if (iH2shield == 3) then
                  l_H2shield = c_ljeans *
     &                 sqrt(tgas1d(i) / (d(i,j,k) * mmw(i)))

               else
                  l_H2shield = 0._RKIND
               endif

               N_H2 = dom*H2I(i,j,k) * l_H2shield

               x = 2.0E-15_DKIND * N_H2
               b_doppler = 1E-5_DKIND *
     &                 sqrt(2._DKIND * kboltz *
     &                      tgas1d(i) / mass_h)
               f_shield = 0.965_DKIND / (1._DKIND + x/b_doppler)**1.1+
     &           0.035_DKIND * exp(-8.5E-4_DKIND * sqrt(1._DKIND +x))/
     &           sqrt(1._DKIND + x)
               k31shield(i) = f_shield * k31shield(i)
            endif
         enddo
      endif
      endif

      if (iradshield > 0) then
!     Compute shielding factors
        do i = is+1, ie+1
          if (itmask(i)) then

!         Compute shielding factor for H
            nSSh = 6.73e-3_DKIND *
     &           (avgsighi /2.49e-18_DKIND)**(-2._DKIND/3._DKIND) *
     &           (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND) *
     &           (k24/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)

!           Compute the total Hydrogen number density
            nratio = (HI(i,j,k) + HII(i,j,k))
            if (ispecies .gt. 1) then
              nratio = nratio +
     &                          HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)

              if (ispecies .gt. 2) then
                nratio = nratio +
     &                       0.5_DKIND*(DI(i,j,k) + DII(i,j,k)) +
     &                       2.0_DKIND*HDI(i,j,k)/3.0_DKIND
              endif
            endif

            nratio = nratio*dom/nSSh

            f_shield_H(i) = (0.98_DKIND*
     &           (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &            0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))

!       Compute shielding factor for He

            nSSh = 6.73e-3_DKIND *
     &           (avgsighei/2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &           (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &           (k26/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)

            nratio = 0.25_DKIND*
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))*dom/nSSh

            f_shield_He(i) = (0.98_DKIND*
     &           (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &            0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))

          endif
        enddo
      endif

      if (iradshield == 1) then
!
!     approximate self shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     to shield HI, while leaving HeI and HeII optically thin
!
!       Attenuate radiation rates for direct H2 ionization (15.4 eV)
!       using same scaling. (rate k29)
!
        do i = is+1, ie+1
          if (itmask(i)) then

            if (k24 .lt. tiny8) then
              k24shield(i) = 0._DKIND
            else
              k24shield(i) = k24shield(i)*f_shield_H(i)
            endif

!     Scale H2 direct ionization radiation
            if (k29 .lt. tiny8) then
              k29shield(i) = 0._DKIND
            else
              k29shield(i) = k29shield(i)*f_shield_H(i)
            endif

            k25shield(i) = k25
            k26shield(i) = k26
          endif
        enddo

      else if (iradshield == 2) then
!
!     Better self-shielding in HI using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     approximate self shielding in HeI and HeII
!
!       Attenuate radiation rates for direct H2 ionization (15.4 eV)
!       using same scaling as HI. (rate k29)
!
!       Attenuate radiation rates for H2+ dissociation (30 eV)
!       using same scaling as HeII. (rate k28 and k30)
!

        do i = is+1, ie+1
          if (itmask(i)) then

            if (k24 .lt. tiny8) then
               k24shield(i) = 0._DKIND
            else
               k24shield(i) = k24shield(i)*f_shield_H(i)
             endif

!     Scale H2 direct ionization radiation
            if (k29 .lt. tiny8) then
              k29shield(i) = 0._DKIND
            else
              k29shield(i) = k29shield(i)*f_shield_H(i)
            endif

!
!     Apply same equations to HeI (assumes HeI closely follows HI)
!

            if (k26 .lt. tiny8) then
               k26shield(i) = 0._DKIND
            else
               k26shield(i) = k26shield(i)*f_shield_He(i)
            endif

!     Scale H2+ dissociation radiation
            if (k28 .lt. tiny8) then
                k28shield(i) = 0.0_DKIND
            else
                k28shield(i) = k28shield(i)*f_shield_He(i)
            endif

            if (k30 .lt. tiny8) then                 
                k30shield(i) = 0.0_DKIND
            else
                k30shield(i) = k30shield(i)*f_shield_He(i)
            endif

            k25shield(i) = k25
          endif
        enddo

      else if (iradshield == 3) then
!
!     shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     in HI and HeI, but ignoring HeII heating entirely
!
        do i = is+1, ie+1
          if (itmask(i)) then

            if (k24 .lt. tiny8) then
               k24shield(i) = 0._DKIND
            else
               k24shield(i)=k24shield(i)*f_shield_H(i)
            endif

!     Scale H2 direct ionization radiation
            if (k29 .lt. tiny8) then
              k29shield(i) = 0._DKIND
            else
              k29shield(i) = k29shield(i)*f_shield_H(i)
            endif

!
!     Apply same equations to HeI (assumes HeI closely follows HI)
!

            if (k26 .lt. tiny8) then
              k26shield(i) = 0._DKIND
            else
              k26shield(i) = k26shield(i)*f_shield_He(i)
            endif

!     Scale H2+ dissociation radiation
            if (k28 .lt. tiny8) then                 
              k28shield(i) = 0.0_DKIND
            else
              k28shield(i) = k28shield(i)*f_shield_He(i)
            endif

            if (k30 .lt. tiny8) then
              k30shield(i) = 0.0_DKIND
            else
              k30shield(i) = k30shield(i)*f_shield_He(i)
            endif

            k25shield(i) = 0.0_DKIND

          endif
        enddo

      endif

#ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
!        If using a high-energy radiation field, then account for
!          effects of secondary electrons (Shull * Steenberg 1985)
!          (see calc_rate.src)

      do i = is+1, ie+1
         if (itmask(i)) then
            x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4_DKIND)
            factor = 0.3908_DKIND*(1._DKIND -
     &           x**0.4092_DKIND)**1.7592_DKIND
            k24shield(i) = k24shield(i) + 
     &           factor*(piHI + 0.08_DKIND*piHeI)/
     &           (e24*everg) *coolunit*tbase1
            factor = 0.0554_DKIND*
     &           (1._DKIND - x**0.4614_DKIND)**1.6660_DKIND
            k26shield(i) = k26shield(i) + 
     &           factor*(piHI/0.08_DKIND + piHeI)/
     &           (e26*everg) *coolunit*tbase1
         endif
      enddo
#endif


!           If using H2, and using the density-dependent collisional
!             H2 dissociation rate, then replace the the density-independant
!                k13 rate with the new one.
!         May/00: there appears to be a problem with the density-dependent
!             collisional rates.  Currently turned off until further notice.

#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
            if (ispecies .gt. 1 .and. ithreebody .eq. 0) then
               do i = is+1, ie+1
                  if (itmask(i)) then
                  nh = min(HI(i,j,k)*dom, 1.0e9_DKIND)
                  k13(i) = tiny8
                  if (tgas1d(i) .ge. 500._DKIND .and.
     &                tgas1d(i) .lt. 1.0e6_DKIND) then
c Direct collisional dissociation
                     k13_CID = k13dd(i,1)-k13dd(i,2)/
     &                          (1._DKIND+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1._DKIND+(nh/k13dd(i,6))**k13dd(i,7))
                     k13_CID = max(10._DKIND**k13_CID, tiny8)
c Dissociative tunnelling
                     k13_DT  = k13dd(i,8)-k13dd(i,9)/
     &                          (1._DKIND+(nh/k13dd(i,12))**k13dd(i,14))
     &                     + k13dd(i,10)-k13dd(i,11)/
     &                          (1._DKIND+(nh/k13dd(i,13))**k13dd(i,14))
                     k13_DT  = max(10._DKIND**k13_DT, tiny8)
c
                     k13(i)  = k13_DT + k13_CID
                  endif
               endif
               enddo
            endif
!!#define USE_PALLA_SALPETER_STAHLER1983
!!#if defined(GRACKLE_MD) && defined(USE_PALLA_SALPETER_STAHLER1983)
!!          if (ispecies .gt. 1 .and. ithreebody .eq. 1) then
!!             do i = is+1, ie+1
!!                if (itmask(i)) then
!!                nh = (HI(i,j,k) + H2I(i,j,k)/2._DKIND)*dom
!!                k13ind = 1._DKIND / (1._DKIND + nh / k13dd(i,3))
!!                k13(i) = 10._DKIND**(
!!   &                     (1._DKIND-k13ind) * k13dd(i,2)
!!   &                             + k13ind  * k13dd(i,1) )
!!                endif
!!             enddo
!!          endif
!!#endif
#endif /*  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE */
!!    if(itmask(1))
!!   &   write(*,'(A1,36E9.1)') 'c',
!!   &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
!!   &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
!!   &                     k24, k25, k26,      k28, k29, k30,
!!   &                     k50, k51, k52, k53, k54, k55, k56, k57, k58
!!!                                                       1    1
!!    if(itmask(1))
!!   &   write(*,*) 'c', k13, k22

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
     &                     itmask, edot, chunit, dom)

! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     arguments

      integer ispecies, is, ie, j, k, in, jn, kn,
     &        iradtrans, irt_honly
      real*8 dedot(in), HIdot(in), dom
      real*8 edot(in)
      logical itmask(in), anydust

!     Density fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)

!      Radiative Transfer Fields
      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn)

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

      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then
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
            if (itmask(i)) then
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

            if (anydust) then
#ifdef GRACKLE_MD
               HIdot(i) = HIdot(i) 
     &              - 2._DKIND * h2dust(i)
#else
               HIdot(i) = HIdot(i) 
     &              - 2._DKIND * h2dust(i) * rhoH(i)
#endif
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

!     H2 formation heating

!     Equation 23 from Omukai (2000)
            h2heatfac(i) = (1._DKIND + (ncrn(i) / (dom *
     &           (HI(i,j,k) * ncrd1(i) +
     &           H2I(i,j,k) * 0.5_DKIND * ncrd2(i)))))**(-1._DKIND)

            ! We only want to apply this if the formation dominates, but we
            ! need to apply it outside the delta calculation.
#ifdef CORRECTED_BY_CHIAKI2020FEB
            H2delta(i) = 
     &          HI(i,j,k) *
     &           ( (3.53_DKIND * k8 (i) * HM(i,j,k)
     &            + 4.48_DKIND * k22(i) * HI(i,j,k)**2._DKIND)
     &            * h2heatfac(i)
     &            - 4.48_DKIND * k13(i) * H2I(i,j,k)/2._DKIND)
#else
            if(H2delta(i).gt.0._DKIND) then
              H2delta(i) = H2delta(i) * h2heatfac(i)
            endif
#endif

            if (anydust) then
#ifdef GRACKLE_MD
               H2delta(i) = H2delta(i) + 
     &              h2dust(i) * HI(i,j,k) * 
     &              (0.2_DKIND + 4.2_DKIND * h2heatfac(i))
#else
               H2delta(i) = H2delta(i) + 
     &              h2dust(i) * HI(i,j,k) * rhoH(i) * 
     &              (0.2_DKIND + 4.2_DKIND * h2heatfac(i))
#endif
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
#ifdef CHEM_DEBUG
CCCCC chunit = (1.60218e-12_DKIND)/(2._DKIND*(xbase1/tbase1)**2*mh)   ! 1 eV per H2 formed
CCCCC    dom = 1.333333333333333
CCCCC    tbase1 = 1406663214374408.5
CCCCC    xbase1 = 8.4155002150553872E+020
CCCCC    coolunit*dom/mh = vbase**2 / tbase1 =  2.5444239786382621E-004
CCCCC coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
CCCCC          = vbase1**2 * mh / dom / tbase1
!        if((i.eq.1).and.(j.eq.1).and.(k.eq.1)) then
!!!         write(*,*)
!!!  &     d(i,j,k)
!!!  &   , k22(i) * 0.75**3 * dom / 1406663214374408.5
!!!  &   , HI(i,j,k) / 0.75
!!!  &   , k22(i) * HI(i,j,k)**3 / 0.75 /1406663214374408.5
!!          write(*,*)
!!   &     d(i,j,k)
!!   &   , 4.48_DKIND*(1.60218e-12_DKIND)
!!   &   , k22(i) * HI(i,j,k)**3 / 0.75 /1406663214374408.5
!!   &   , h2heatfac(i)
!!   &   , 4.48_DKIND*(1.60218e-12_DKIND) ! erg
!!   &   * k22(i) * HI(i,j,k)**3 / 0.75 /1406663214374408.5 ! /cc/s
!!   &   * h2heatfac(i)
!!   &   / (d(i,j,k) / 0.75 * 1.6726d-24) ! / (g/cc)
!!   &   , HI(i,j,k) *
!!   &      (  4.48_DKIND * k22(i) * HI(i,j,k)**2._DKIND )
!!   &       * h2heatfac(i)
!!   &       /d(i,j,k)*1.60218e-12_DKIND
!!   &        / 2._DKIND / 1.6726d-24 / 1406663214374408.5
!           write(*,*)
!    &     d(i,j,k)
!    &   , h2heatfac(i)
!    &   , HI(i,j,k) *
!    &      (  3.53_DKIND * k8 (i) * HM(i,j,k) )
!    &       * h2heatfac(i)
!    &       *chunit/d(i,j,k)*2.5444239786382621E-004
!    &   , HI(i,j,k) *
!    &      (  4.48_DKIND * k22(i) * HI(i,j,k)**2._DKIND )
!    &       * h2heatfac(i)
!    &       *chunit/d(i,j,k)*2.5444239786382621E-004
!    &   , HI(i,j,k) *
!    &      (- 4.48_DKIND * k13(i) * H2I(i,j,k)/2._DKIND )
!    &       *chunit/d(i,j,k)*2.5444239786382621E-004
!    &   , h2dust(i) * HI(i,j,k) * 
!    &     (0.2_DKIND + 4.2_DKIND * h2heatfac(i))
!    &       *chunit/d(i,j,k)*2.5444239786382621E-004
!    &   , chunit*H2delta(i)/d(i,j,k)*2.5444239786382621E-004
!        endif
#endif
         endif                  ! itmask
         enddo
      endif

!     Add photo-ionization rates if needed

      if (iradtrans .eq. 1) then
         if (irt_honly .eq. 0) then
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     &                 + kphHeI(i,j,k) * HeI(i,j,k) / 4._DKIND
     &                 + kphHeII(i,j,k) * HeII(i,j,k) / 4._DKIND
               endif
            enddo
         else
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
               endif
            enddo
         endif
      endif

      


      return
      end


! -----------------------------------------------------------
!  This routine uses one linearly implicit Gauss-Seidel sweep of 
!   a backward-Euler time integrator to advance the rate equations 
!   by one (sub-)cycle (dtit).

      subroutine step_rate_g(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies, anydust,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58,
     &                     h2dust, rhoH,
     &                     k24shield, k25shield, k26shield, 
     &                     k28shield, k29shield, k30shield, k31shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly,
     &                     kphHI, kphHeI, kphHeII,
     &                     itmask
#ifdef GRACKLE_MD
     &                   , DM, HDII, HeHII, imetal, metal, imchem, igrgr
     &                   , CI, CII, CO, CO2
     &                   , OI, OH, H2O, O2
     &                   , SiI, SiOI, SiO2I
     &                   , CH, CH2, COII, OII
     &                   , OHII, H2OII, H3OII, O2II
     &                   , Mg, Al, S, Fe
     &                   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                   , AC, SiO2D, MgO, FeS, Al2O3
     &                   , k125, k129, k130, k131, k132
     &                   , k133, k134, k135, k136, k137
     &                   , k148, k149, k150, k151, k152
     &                   , k153
     &                   , kz15 , kz16 , kz17 , kz18 , kz19 
     &                   , kz20 , kz21 , kz22 , kz23 , kz24 
     &                   , kz25 , kz26 , kz27 , kz28 , kz29 
     &                   , kz30 , kz31 , kz32 , kz33 , kz34 
     &                   , kz35 , kz36 , kz37 , kz38 , kz39 
     &                   , kz40 , kz41 , kz42 , kz43 , kz44 
     &                   , kz45 , kz46 , kz47 , kz48 , kz49 
     &                   , kz50 , kz51 , kz52 , kz53 , kz54 
     &                   , DMp, HDIIp, HeHIIp
     &                   , CIp, CIIp, COp, CO2p
     &                   , OIp, OHp, H2Op, O2p
     &                   , SiIp, SiOIp, SiO2Ip
     &                   , CHp, CH2p, COIIp, OIIp
     &                   , OHIIp, H2OIIp, H3OIIp, O2IIp
     &                   , Mgp, Alp, Sp, Fep
     &                   , SiMp, FeMp, Mg2SiO4p, MgSiO3p, Fe3O4p
     &                   , ACp, SiO2Dp, MgOp, FeSp, Al2O3p
     &                   , kd0, kd1, kd2, kd3, kd4
     &                   , kd5, kd6, kd7, kd8, kd9
#endif
     &           )
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     arguments

      integer ispecies, in, jn, kn, is, ie, j, k,
     &        iradtrans, irt_honly
      real*8 dtit(in), dedot_prev(in), HIdot_prev(in)
      logical itmask(in), anydust

!     Density fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
#ifdef GRACKLE_MD
       integer imetal, imchem, igrgr
       R_PREC metal(in,jn,kn)
       R_PREC DM(in,jn,kn)   , HDII(in,jn,kn)   , HeHII(in,jn,kn)
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
#endif

!     Radiation Fields
      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn)


!     Rate values

      real*8 k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in), 
     &       h2dust(in), rhoH(in),
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k28shield(in), k29shield(in), k30shield(in),
     &       k31shield(in),
     &       k24, k25, k26, k27, k28, k29, k30
#ifdef GRACKLE_MD
      real*8 k125(in), k129(in), k130(in), k131(in), k132(in)
     &     , k133(in), k134(in), k135(in), k136(in), k137(in)
     &     , k148(in), k149(in), k150(in), k151(in), k152(in)
     &     , k153(in)
     &     , kz15(in), kz16(in), kz17(in), kz18(in), kz19(in)
     &     , kz20(in), kz21(in), kz22(in), kz23(in), kz24(in)
     &     , kz25(in), kz26(in), kz27(in), kz28(in), kz29(in)
     &     , kz30(in), kz31(in), kz32(in), kz33(in), kz34(in)
     &     , kz35(in), kz36(in), kz37(in), kz38(in), kz39(in)
     &     , kz40(in), kz41(in), kz42(in), kz43(in), kz44(in)
     &     , kz45(in), kz46(in), kz47(in), kz48(in), kz49(in)
     &     , kz50(in), kz51(in), kz52(in), kz53(in), kz54(in)
#endif

!     temporaries (passed in)

      real*8 HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &       HMp(in), H2Ip(in), H2IIp(in), dep(in),
     &       DIp(in), DIIp(in), HDIp(in)
#ifdef GRACKLE_MD
      real*8 DMp(in)   , HDIIp(in)   , HeHIIp(in)
     &     , CIp(in)   , CIIp(in)    , COp(in)
     &     , CO2p(in)  , OIp(in)     , OHp(in)
     &     , H2Op(in)  , O2p(in)     , SiIp(in)
     &     , SiOIp(in) , SiO2Ip(in)  , CHp(in)
     &     , CH2p(in)  , COIIp(in)   , OIIp(in)
     &     , OHIIp(in) , H2OIIp(in)  , H3OIIp(in)
     &     , O2IIp(in) , Mgp(in)     , Alp(in)
     &     , Sp(in)    , Fep(in)     , SiMp(in)
     &     , FeMp(in)  , Mg2SiO4p(in), MgSiO3p(in)
     &     , Fe3O4p(in), ACp(in)     , SiO2Dp(in)
     &     , MgOp(in)  , FeSp(in)    , Al2O3p(in)
      real*8  kd0(in), kd1(in), kd2(in), kd3(in), kd4(in)
     &      , kd5(in), kd6(in), kd7(in), kd8(in), kd9(in)
#endif

!     locals

      integer i
      real*8 scoef, acoef

!   A) the 6-species integrator
!      
      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI

            scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
            acoef  = k1(i)*de(i,j,k)
     &             + k57(i)*HI(i,j,k)
     &             + k58(i)*HeI(i,j,k)/4._DKIND
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
            HIp(i)  = (scoef*dtit(i) + HI(i,j,k))/
     &           (1._DKIND + acoef*dtit(i))
            if (HIp(i) .ne. HIp(i)) then
#ifdef _OPENMP
!$omp critical
#endif
               write(*,*) 'HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), kphHI(i,j,k),
     $              scoef, acoef, dtit(i)
#ifdef _OPENMP
!$omp end critical
#endif
c               ERROR_MESSAGE
            endif

!        2) HII
c 
            scoef  = k1(i)*HIp(i)*de(i,j,k)
     &             + k57(i)*HIp(i)*HIp(i)
     &             + k58(i)*HIp(i)*HeI(i,j,k)/4._DKIND
     &             + k24shield(i)*HIp(i)
            if (iradtrans .eq. 1) 
     &          scoef = scoef + kphHI(i,j,k)*HIp(i)
            acoef  = k2(i)*de (i,j,k)
            HIIp(i) = (scoef*dtit(i) + HII(i,j,k))/
     &           (1._DKIND +acoef*dtit(i))
!
            if (HIIp(i) .le. 0._DKIND) then  !#####
#ifdef _OPENMP
!$omp critical
#endif
               write(*,*) 'negative HIIp! :: ', i, j, k, HIIp(i), 
     $              scoef, dtit(i), HII(i,j,k), acoef,
     $              k2(i), de(i,j,k),
     $              kphHI(i,j,k), HIp(i),
     $              k24shield(i)
#ifdef _OPENMP
!$omp end critical
#endif
            endif

!        3) Electron density

            scoef = 0._DKIND
     &                 + k57(i)*HIp(i)*HIp(i)
     &                 + k58(i)*HIp(i)*HeI(i,j,k)/4._DKIND
     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4._DKIND
     &                 + k26shield(i)*HeI(i,j,k)/4._DKIND

            if ( (iradtrans .eq. 1) .and. ( irt_honly .eq. 0) )
     &          scoef = scoef + kphHI(i,j,k) * HI(i,j,k)
     &                + kphHeI(i,j,k)  * HeI(i,j,k)  / 4._DKIND
     &                + kphHeII(i,j,k) * HeII(i,j,k) / 4._DKIND
            if ( (iradtrans .eq. 1) .and. ( irt_honly .eq. 1) )
     &          scoef = scoef + kphHI(i,j,k) * HI(i,j,k)



            acoef = -(k1(i)*HI(i,j,k)      - k2(i)*HII(i,j,k)
     &              + k3(i)*HeI(i,j,k)/4._DKIND -
     &           k6(i)*HeIII(i,j,k)/4._DKIND
     &              + k5(i)*HeII(i,j,k)/4._DKIND -
     &           k4(i)*HeII(i,j,k)/4._DKIND)
            dep(i)   = (scoef*dtit(i) + de(i,j,k))
     &                     / (1._DKIND + acoef*dtit(i))

         endif                  ! itmask
         enddo

      endif                     ! (ispecies .eq. 1)

!  --- (B) Do helium chemistry in any case: (for all ispecies values) ---

      do i = is+1, ie+1
         if (itmask(i)) then

!        4) HeI

         scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
         acoef  = k3(i)*de(i,j,k)
     &                + k26shield(i)

         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       acoef = acoef + kphHeI(i,j,k)
#ifdef GRACKLE_MD
         if (ispecies .gt. 3) then
            scoef = scoef +  4._DKIND * ( 0._DKIND
     &          + k152(i) * HeHII(i,j,k) *    HI(i,j,k) /  5._DKIND
     &          + k153(i) * HeHII(i,j,k) *    de(i,j,k) /  5._DKIND
     &          )
            acoef = acoef
     &          + k148(i) *   HII(i,j,k)
     &          + k149(i) *   HII(i,j,k)
     &          + k150(i) *  H2II(i,j,k) /  2._DKIND
         endif
#endif
         HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &              / ( 1._DKIND + acoef*dtit(i) )

!        5) HeII

         scoef  = k3(i)*HeIp(i)*de(i,j,k)
     &          + k6(i)*HeIII(i,j,k)*de(i,j,k)
     &          + k26shield(i)*HeIp(i)
     
         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       scoef = scoef + kphHeI(i,j,k)*HeIp(i)

         acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)
     &          + k25shield(i)
     
         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       acoef = acoef + kphHeII(i,j,k)
#ifdef GRACKLE_MD
         if (ispecies .gt. 3) then
            acoef = acoef
     &          + k151(i) *    HI(i,j,k)
         endif
#endif
         HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &              / ( 1._DKIND + acoef*dtit(i) )

!       6) HeIII

         scoef   = k5(i)*HeIIp(i)*de(i,j,k)
     &           + k25shield(i)*HeIIp(i)
         if ((iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       scoef = scoef + kphHeII(i,j,k) * HeIIp(i)
         acoef   = k6(i)*de(i,j,k)
         HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                / ( 1._DKIND + acoef*dtit(i) )

      endif                     ! itmask
      enddo

c --- (C) Now do extra 3-species for molecular hydrogen ---

      if (ispecies .gt. 1) then

!        First, do HI/HII with molecular hydrogen terms

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI
!     
            scoef  =      k2(i) * HII(i,j,k) * de(i,j,k) 
     &             + 2._DKIND*k13(i)* HI(i,j,k)  * H2I(i,j,k)/2._DKIND
     &             +      k11(i)* HII(i,j,k) * H2I(i,j,k)/2._DKIND
     &             + 2._DKIND*k12(i)* de(i,j,k)  * H2I(i,j,k)/2._DKIND
     &             +      k14(i)* HM(i,j,k)  * de(i,j,k)
     &             +      k15(i)* HM(i,j,k)  * HI(i,j,k)
     &             + 2._DKIND*k16(i)* HM(i,j,k)  * HII(i,j,k)
     &             + 2._DKIND*k18(i)* H2II(i,j,k)* de(i,j,k)/2._DKIND
     &             +      k19(i)* H2II(i,j,k)* HM(i,j,k)/2._DKIND
     &             + 2._DKIND*k31shield(i)   * H2I(i,j,k)/2._DKIND

            acoef  =      k1(i) * de(i,j,k)
     &             +      k7(i) * de(i,j,k)  
     &             +      k8(i) * HM(i,j,k)
     &             +      k9(i) * HII(i,j,k)
     &             +      k10(i)* H2II(i,j,k)/2._DKIND
     &             + 2._DKIND*k22(i)* HI(i,j,k)**2
     &             +      k57(i)* HI(i,j,k)
     &             +      k58(i)* HeI(i,j,k)/4._DKIND
     &             + k24shield(i)
!!!         if((i.eq.10).and.(j.eq.10).and.(k.eq.7))
!!!  &        write(*,*) 'HI a', scoef, acoef

            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
!!!         if((i.eq.10).and.(j.eq.10).and.(k.eq.7))
!!!  &        write(*,*) 'HI b', scoef, acoef

            if (anydust) then
#ifdef GRACKLE_MD
               acoef = acoef + 2._DKIND * h2dust(i)
#else
               acoef = acoef + 2._DKIND * h2dust(i) * rhoH(i)
#endif
            endif
!!!         if((i.eq.10).and.(j.eq.10).and.(k.eq.7))
!!!  &        write(*,*) 'HI c', scoef, acoef
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
            if (ispecies .gt. 2) then
               scoef = scoef 
     &               + k50(i) * HII(i,j,k) * DI(i,j,k)  / 2._DKIND
     &               + k54(i) * H2I(i,j,k) * DI(i,j,k)  / 4._DKIND
               acoef = acoef
     &               + k51(i) * DII(i,j,k) / 2._DKIND
     &               + k55(i) * HDI(i,j,k) / 3._DKIND
            endif
#endif
!!!         if((i.eq.10).and.(j.eq.10).and.(k.eq.7))
!!!  &        write(*,*) 'HI d', scoef, acoef
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k131(i) *  HDII(i,j,k) *    de(i,j,k) /  3._DKIND
     &             + k134(i) *   HII(i,j,k) *    DM(i,j,k) /  2._DKIND
     &             + k135(i) *    HM(i,j,k) *    DI(i,j,k) /  2._DKIND
     &             + k150(i) *   HeI(i,j,k) *  H2II(i,j,k) /  8._DKIND
     &             + k153(i) * HeHII(i,j,k) *    de(i,j,k) /  5._DKIND
               acoef = acoef
     &             + k125(i) *  HDII(i,j,k) /  3._DKIND
     &             + k130(i) *   DII(i,j,k) /  2._DKIND
     &             + k136(i) *    DM(i,j,k) /  2._DKIND
     &             + k137(i) *    DM(i,j,k) /  2._DKIND
     &             + k151(i) *  HeII(i,j,k) /  4._DKIND
     &             + k152(i) * HeHII(i,j,k) /  5._DKIND
            endif
!!!         if((i.eq.10).and.(j.eq.10).and.(k.eq.7))
!!!  &        write(*,*) 'HI e', scoef, acoef
            if (   (imchem .eq. 1) 
     &       .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) ) then
               scoef = scoef
     &             + kz20(i) *    CI(i,j,k) *   H2I(i,j,k) / 24._DKIND
     &             + kz21(i) *    OI(i,j,k) *   H2I(i,j,k) / 32._DKIND
     &             + kz22(i) *   HII(i,j,k) *    OI(i,j,k) / 16._DKIND
     &             + kz23(i) *   H2I(i,j,k) *    CH(i,j,k) / 26._DKIND
     &             + kz24(i) *   H2I(i,j,k) *    OH(i,j,k) / 34._DKIND
     &             + kz26(i) *    OH(i,j,k) *    CO(i,j,k) / 476._DKIND
     &             + kz28(i) *    CI(i,j,k) *    OH(i,j,k) / 204._DKIND
     &             + kz32(i) *    OI(i,j,k) *    CH(i,j,k) / 208._DKIND
     &             + kz33(i) *    OI(i,j,k) *    OH(i,j,k) / 272._DKIND
     &             + kz34(i) *   HII(i,j,k) *    OH(i,j,k) / 17._DKIND
     &             + kz35(i) *   HII(i,j,k) *   H2O(i,j,k) / 18._DKIND
     &             + kz36(i) *   HII(i,j,k) *    O2(i,j,k) / 32._DKIND
     &             + kz37(i) *   CII(i,j,k) *    OH(i,j,k) / 204._DKIND
     &             + kz40(i) *   OII(i,j,k) *   H2I(i,j,k) / 32._DKIND
     &             + kz41(i) *  OHII(i,j,k) *   H2I(i,j,k) / 34._DKIND
     &             + kz42(i) * H2OII(i,j,k) *   H2I(i,j,k) / 36._DKIND
     &             + kz46(i) * H2OII(i,j,k) *    de(i,j,k) / 18._DKIND
     &             + kz48(i) * H3OII(i,j,k) *    de(i,j,k) / 19._DKIND
     &             + kz49(i) * H3OII(i,j,k) *    de(i,j,k) / 9.5_DKIND
     &             + kz52(i) *   SiI(i,j,k) *    OH(i,j,k) / 476._DKIND
     &             + kz54(i) *  SiOI(i,j,k) *    OH(i,j,k) / 748._DKIND
               acoef = acoef
     &             + kz15(i) *    CH(i,j,k) / 13._DKIND
     &             + kz16(i) *   CH2(i,j,k) / 14._DKIND
     &             + kz17(i) *    OH(i,j,k) / 17._DKIND
     &             + kz18(i) *   H2O(i,j,k) / 18._DKIND
     &             + kz19(i) *    O2(i,j,k) / 32._DKIND
     &             + kz27(i) *    CI(i,j,k) / 12._DKIND
     &             + kz30(i) *    OI(i,j,k) / 16._DKIND
     &             + kz39(i) *   OII(i,j,k) / 16._DKIND
     &             + kz43(i) *  COII(i,j,k) / 28._DKIND
            endif
!!!         if((i.eq.10).and.(j.eq.10).and.(k.eq.7))
!!!  &        write(*,*) 'HI f', scoef, acoef
#endif 
            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1. + acoef*dtit(i) )
            if (HIp(i) .ne. HIp(i)) then
#ifdef _OPENMP
!$omp critical
#endif
               write(*,*) 'HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), H2I(i,j,k),
     $              kphHI(i,j,k)
#ifdef CORRECTED_BY_CHIAKI2020MAY
        write(*,*) 'p1', HeI(i,j,k), HeII(i,j,k), HeIII(i,j,k)
        write(*,*) 'p2', HM(i,j,k), H2II(i,j,k)
        if (ispecies .gt. 2) then
        write(*,*) 'p3', DI(i,j,k), DII(i,j,k), HDI(i,j,k)
        endif
        if (ispecies .gt. 3) then
        write(*,*) 'p4', DM(i,j,k), HDII(i,j,k), HeHII(i,j,k)
        endif
        if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
        if (imchem .eq. 1) then
        write(*,*) 'm1', CI(i,j,k), CII(i,j,k), CO(i,j,k), CO2(i,j,k)
        write(*,*) 'm2', OI(i,j,k), OH(i,j,k), H2O(i,j,k), O2(i,j,k)
        write(*,*) 'm3', SiI(i,j,k), SiOI(i,j,k), SiO2I(i,j,k)
        write(*,*) 'm4', CH(i,j,k), CH2(i,j,k), COII(i,j,k), OII(i,j,k)
        write(*,*) 'm5', OHII(i,j,k), H2OII(i,j,k)
     &           , H3OII(i,j,k), O2II(i,j,k)
        endif
        if (igrgr .eq. 1) then
        write(*,*) 'm6', Mg(i,j,k), Al(i,j,k), S(i,j,k), Fe(i,j,k)
        write(*,*) 'd1', SiM(i,j,k), FeM(i,j,k), Mg2SiO4(i,j,k)
     &           , MgSiO3(i,j,k), Fe3O4(i,j,k)
        write(*,*) 'd2', AC(i,j,k), SiO2D(i,j,k), MgO(i,j,k)
     &           , FeS(i,j,k), Al2O3(i,j,k)
        endif
        endif
#endif
#ifdef _OPENMP
!$omp end critical
#endif
            endif

!          2) HII

            scoef  =    k1(i)  * HI(i,j,k) * de(i,j,k)
     &             +    k10(i) * H2II(i,j,k)*HI(i,j,k)/2._DKIND
     &             +    k57(i) * HI(i,j,k) * HI(i,j,k)
     &             +    k58(i) * HI(i,j,k) * HeI(i,j,k)/4._DKIND
     &             + k24shield(i)*HI(i,j,k)

            if (iradtrans .eq. 1) 
     &          scoef = scoef + kphHI(i,j,k) * HI(i,j,k)

            acoef  =    k2(i)  * de(i,j,k)
     &             +    k9(i)  * HI(i,j,k)
     &             +    k11(i) * H2I(i,j,k)/2._DKIND
     &             +    k16(i) * HM(i,j,k)
     &             +    k17(i) * HM(i,j,k)
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
            if (ispecies .gt. 2) then
               scoef = scoef
     &               + k51(i) * HI (i,j,k) * DII(i,j,k) / 2._DKIND
     &               + k52(i) * H2I(i,j,k) * DII(i,j,k) / 4._DKIND
               acoef = acoef
     &               + k50(i) * DI (i,j,k) / 2._DKIND
     &               + k53(i) * HDI(i,j,k) / 3._DKIND
            endif
#endif
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k125(i) *  HDII(i,j,k) *    HI(i,j,k) /  3._DKIND
               acoef = acoef
     &             + k129(i) *    DI(i,j,k) /  2._DKIND
     &             + k134(i) *    DM(i,j,k) /  2._DKIND
     &             + k148(i) *   HeI(i,j,k) /  4._DKIND
     &             + k149(i) *   HeI(i,j,k) /  4._DKIND
            endif
            if (   (imchem .eq. 1) 
     &       .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) ) then
               scoef = scoef
     &             + kz39(i) *   OII(i,j,k) *    HI(i,j,k) / 16._DKIND
     &             + kz43(i) *  COII(i,j,k) *    HI(i,j,k) / 28._DKIND
               acoef = acoef
     &             + kz22(i) *    OI(i,j,k) / 16._DKIND
     &             + kz34(i) *    OH(i,j,k) / 17._DKIND
     &             + kz35(i) *   H2O(i,j,k) / 18._DKIND
     &             + kz36(i) *    O2(i,j,k) / 32._DKIND
            endif
#endif
            HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1._DKIND + acoef*dtit(i) )
!     
!          3) electrons:

            scoef =   k8(i) * HM(i,j,k) * HI(i,j,k)
     &             +  k15(i)* HM(i,j,k) * HI(i,j,k)
     &             +  k17(i)* HM(i,j,k) * HII(i,j,k)
     &             +  k57(i)* HI(i,j,k) * HI(i,j,k)
     &             +  k58(i)* HI(i,j,k) * HeI(i,j,k)/4._DKIND
!                  
     &             + k24shield(i)*HIp(i)
     &             + k25shield(i)*HeIIp(i)/4._DKIND
     &             + k26shield(i)*HeIp(i)/4._DKIND

            if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0) )
     &          scoef = scoef + kphHI(i,j,k) * HIp(i)
     &                + kphHeI(i,j,k)  * HeIp(i)  / 4._DKIND
     &                + kphHeII(i,j,k) * HeIIp(i) / 4._DKIND
            if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 1) )
     &          scoef = scoef + kphHI(i,j,k) * HIp(i)

            acoef = - (k1(i) *HI(i,j,k)    - k2(i)*HII(i,j,k)
     &              +  k3(i) *HeI(i,j,k)/4._DKIND -
     &           k6(i)*HeIII(i,j,k)/4._DKIND
     &              +  k5(i) *HeII(i,j,k)/4._DKIND -
     &           k4(i)*HeII(i,j,k)/4._DKIND
     &              +  k14(i)*HM(i,j,k)
     &              -  k7(i) *HI(i,j,k)
     &              -  k18(i)*H2II(i,j,k)/2._DKIND)
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
            if (ispecies .gt. 2) then
               scoef = scoef
     &               + k56(i) * DI (i,j,k) * HM(i,j,k) / 2._DKIND
               acoef = acoef
     &               - k1 (i) * DI (i,j,k) / 2._DKIND
     &               + k2 (i) * DII(i,j,k) / 2._DKIND
            endif
#endif
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k137(i) *    DM(i,j,k) *    HI(i,j,k) /  2._DKIND
               acoef = acoef
     &             + k131(i) *  HDII(i,j,k) /  3._DKIND
     &             + k132(i) *    DI(i,j,k) /  2._DKIND
     &             + k153(i) * HeHII(i,j,k) /  5._DKIND
            endif
            if (   (imchem .eq. 1) 
     &       .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) ) then
               scoef = scoef
               acoef = acoef
     &             + kz44(i) *   CII(i,j,k) / 12._DKIND
     &             + kz45(i) *   OII(i,j,k) / 16._DKIND
     &             + kz46(i) * H2OII(i,j,k) / 18._DKIND
     &             + kz47(i) * H2OII(i,j,k) / 18._DKIND
     &             + kz48(i) * H3OII(i,j,k) / 19._DKIND
     &             + kz49(i) * H3OII(i,j,k) / 19._DKIND
     &             + kz50(i) *  O2II(i,j,k) / 32._DKIND
            endif
#endif
            dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                / ( 1._DKIND + acoef*dtit(i) )

!           7) H2

            scoef = 2._DKIND*(k8(i)  * HM(i,j,k)   * HI(i,j,k)
     &            +       k10(i) * H2II(i,j,k) * HI(i,j,k)/2._DKIND
     &            +       k19(i) * H2II(i,j,k) * HM(i,j,k)/2._DKIND
     &            +       k22(i) * HI(i,j,k) * (HI(i,j,k))**2._DKIND)
            acoef = ( k13(i)*HI(i,j,k) + k11(i)*HII(i,j,k)
     &              + k12(i)*de(i,j,k) )
     &              + k29shield(i) + k31shield(i)

            if (anydust) then
#ifdef GRACKLE_MD
               scoef = scoef + 2._DKIND * h2dust(i) *
     &              HI(i,j,k)
#else
               scoef = scoef + 2._DKIND * h2dust(i) *
     &              HI(i,j,k) * rhoH(i)
#endif
            endif
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
            if (ispecies .gt. 2) then
               scoef = scoef + 2._DKIND * (
     &                 k53(i) * HDI(i,j,k) * HII(i,j,k) / 3._DKIND
     &               + k55(i) * HDI(i,j,k) * HI (i,j,k) / 3._DKIND
     &                  )
               acoef = acoef
     &               + k52(i) * DII(i,j,k) / 2._DKIND
     &               + k54(i) * DI (i,j,k) / 2._DKIND
            endif
#endif
#ifdef GRACKLE_MD
            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
            if (imchem .eq. 1) then
               scoef = scoef +  2._DKIND * ( 0._DKIND
     &             + kz15(i) *    HI(i,j,k) *    CH(i,j,k) / 13._DKIND
     &             + kz16(i) *    HI(i,j,k) *   CH2(i,j,k) / 14._DKIND
     &             + kz17(i) *    HI(i,j,k) *    OH(i,j,k) / 17._DKIND
     &             + kz18(i) *    HI(i,j,k) *   H2O(i,j,k) / 18._DKIND
     &             + kz47(i) * H2OII(i,j,k) *    de(i,j,k) / 18._DKIND
     &            )
               acoef = acoef
     &             + kz20(i) *    CI(i,j,k) / 12._DKIND
     &             + kz21(i) *    OI(i,j,k) / 16._DKIND
     &             + kz23(i) *    CH(i,j,k) / 13._DKIND
     &             + kz24(i) *    OH(i,j,k) / 17._DKIND
     &             + kz40(i) *   OII(i,j,k) / 16._DKIND
     &             + kz41(i) *  OHII(i,j,k) / 17._DKIND
     &             + kz42(i) * H2OII(i,j,k) / 18._DKIND
     &             + kz51(i) *    CI(i,j,k) / 12._DKIND
            endif
            if(igrgr .eq. 1) then
               scoef = scoef + 2._DKIND * (
     &               kd2(i) * 3._DKIND
     &             + kd3(i) * 2._DKIND
     &             + kd4(i) * 4._DKIND
     &             + kd7(i)
     &             + kd9(i) * 3._DKIND
     &           )
            endif
            endif
#endif
            H2Ip(i) = ( scoef*dtit(i) + H2I(i,j,k) )
     &                / ( 1._DKIND + acoef*dtit(i) )

!           8) H-

            scoef = k7(i) * HI(i,j,k) * de(i,j,k) 
            acoef = (k8(i)  + k15(i))  * HI(i,j,k) + 
     &              (k16(i) + k17(i))  * HII(i,j,k) +  
     &	            k14(i) * de(i,j,k) + k19(i) * H2II(i,j,k)/2.0 +
     &	            k27
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
            if (ispecies .gt. 2) then
               acoef = acoef
     &               + k56(i) * DI (i,j,k) / 2._DKIND
            endif
#endif
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               scoef = scoef
     &             + k136(i) *    DM(i,j,k) *    HI(i,j,k) /  2._DKIND
               acoef = acoef
     &             + k135(i) *    DI(i,j,k) /  2._DKIND
            endif
#endif
            HMp(i) = (scoef*dtit(i) + HM(i,j,k))
     &           / (1.0 + acoef*dtit(i))


!           9) H2+

            H2IIp(i) = 2._DKIND*( k9 (i)*HIp(i)*HIIp(i)
     &                    +   k11(i)*H2Ip(i)/2._DKIND*HIIp(i)
     &                    +   k17(i)*HMp(i)*HIIp(i)
     &                    + k29shield(i)*H2Ip(i)
     &                    )
     &                 /  ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                    + k19(i)*HMp(i)
     &                    + (k28shield(i)+k30shield(i))
     &                    )
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
             H2IIp(i) = 2._DKIND * (  k9 (i)*HIp(i)*HIIp(i)
     &                    +   k11(i)*H2Ip(i)/2._DKIND*HIIp(i)
     &                    +   k17(i)*HMp(i)*HIIp(i)
     &                    + k29shield(i)*H2Ip(i)
     &                    + k152(i)*HeHII(i,j,k)*HIp(i)/5._DKIND
     &                    )
     &                 /  ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                    + k19(i)*HMp(i)
     &                    + (k28shield(i)+k30shield(i))
     &                    + k150(i)*HeIp(i)/4._DKIND
     &                    )
            endif
#endif
         endif                  ! itmask
         enddo
!     
      endif                     ! H2

!  --- (D) Now do extra 3-species for molecular HD ---
!     
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
!     
!         1) DI
!     
            scoef =   (       k2(i) * DII(i,j,k) * de(i,j,k)
     &                 +      k51(i)* DII(i,j,k) * HI(i,j,k)
     &                 + 2._DKIND*k55(i)* HDI(i,j,k) *
     &              HI(i,j,k)/3._DKIND
     &                 )
            acoef  =    k1(i) * de(i,j,k)
     &             +    k50(i) * HII(i,j,k)
     &             +    k54(i) * H2I(i,j,k)/2._DKIND
     &             +    k56(i) * HM(i,j,k)
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               scoef = scoef +  2._DKIND * ( 0._DKIND
     &             + k131(i) *  HDII(i,j,k) *    de(i,j,k) /  3._DKIND
     &             + k133(i) *   DII(i,j,k) *    DM(i,j,k) /  2._DKIND
     &             + k134(i) *   HII(i,j,k) *    DM(i,j,k) /  2._DKIND
     &             + k136(i) *    DM(i,j,k) *    HI(i,j,k) /  2._DKIND
     &             )
               acoef = acoef
     &             + k129(i) *   HII(i,j,k)
     &             + k132(i) *    de(i,j,k)
     &             + k135(i) *    HM(i,j,k)
            endif
#endif
            DIp(i)    = ( scoef*dtit(i) + DI(i,j,k) ) / 
     &                  ( 1._DKIND + acoef*dtit(i) )

!         2) DII
c 
            scoef =   (   k1(i)  * DI(i,j,k) * de(i,j,k)
     &            +       k50(i) * HII(i,j,k)* DI(i,j,k)
     &            +  2._DKIND*k53(i) * HII(i,j,k)* HDI(i,j,k)/3._DKIND
     &            )
     &            + k24shield(i)*DI(i,j,k)
#ifdef CORRECTED_BY_CHIAKI2020MAY
            acoef = 0._DKIND
#endif
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*DI(i,j,k)
            acoef =    k2(i)  * de(i,j,k)
     &            +    k51(i) * HI(i,j,k)
     &            +    k52(i) * H2I(i,j,k)/2._DKIND
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               acoef = acoef
     &             + k130(i) *    HI(i,j,k)
     &             + k133(i) *    DM(i,j,k) /  2._DKIND
            endif
#endif
            DIIp(i)   = ( scoef*dtit(i) + DII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

!          3) HDI
c 
            scoef = 3._DKIND*(k52(i) * DII(i,j,k)* 
     &           H2I(i,j,k)/2._DKIND/2._DKIND
     &           + k54(i) * DI(i,j,k) * H2I(i,j,k)/2._DKIND/2._DKIND
#ifdef GRACKLE_MD
     &           +          k56(i) * DI(i,j,k) * HM(i,j,k)/2._DKIND
#else
     &           + 2._DKIND*k56(i) * DI(i,j,k) * HM(i,j,k)/2._DKIND
#endif
     &                 )
            acoef  =    k53(i) * HII(i,j,k)
     &             +    k55(i) * HI(i,j,k)
#ifdef GRACKLE_MD
            if (ispecies .gt. 3) then
               scoef = scoef +  3._DKIND * ( 0._DKIND
     &             + k125(i) *  HDII(i,j,k) *    HI(i,j,k) /  3._DKIND
     &             + k137(i) *    DM(i,j,k) *    HI(i,j,k) /  2._DKIND
     &             )
            endif
#endif
            HDIp(i)   = ( scoef*dtit(i) + HDI(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

         endif                  ! itmask
         enddo
      endif
#ifdef GRACKLE_MD
!  --- (D2) Now do extra 3-species for minor primordial species ---
!     
      if (ispecies .gt. 3) then
         do i = is+1, ie+1
            if (itmask(i)) then
!     
!         1) DM
!     
            scoef =
     &            k132(i) *    DI(i,j,k) *    de(i,j,k)
     &          + k135(i) *    HM(i,j,k) *    DI(i,j,k)
            acoef = 
     &            k133(i) *   DII(i,j,k) /  2._DKIND
     &          + k134(i) *   HII(i,j,k)
     &          + k136(i) *    HI(i,j,k)
     &          + k137(i) *    HI(i,j,k)

            DMp(i)    = ( scoef*dtit(i) + DM(i,j,k) ) / 
     &                  ( 1._DKIND + acoef*dtit(i) )

!         2) HDII
c 
            scoef = 3._DKIND * (
     &            k129(i) *    DI(i,j,k) *   HII(i,j,k) /  2._DKIND
     &          + k130(i) *   DII(i,j,k) *    HI(i,j,k) /  2._DKIND
     &          )
            acoef = 
     &            k125(i) *    HI(i,j,k)
     &          + k131(i) *    de(i,j,k)

            HDIIp(i)   = ( scoef*dtit(i) + HDII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

!          3) HeHII
c 
            scoef = 5._DKIND * (
     &            k148(i) *   HeI(i,j,k) *   HII(i,j,k) /  4._DKIND
     &          + k149(i) *   HeI(i,j,k) *   HII(i,j,k) /  4._DKIND
     &          + k150(i) *   HeI(i,j,k) *  H2II(i,j,k) /  8._DKIND
     &          + k151(i) *  HeII(i,j,k) *    HI(i,j,k) /  4._DKIND
     &          )
            acoef = 
     &            k152(i) *    HI(i,j,k)
     &          + k153(i) *    de(i,j,k)

            HeHIIp(i)   = ( scoef*dtit(i) + HeHII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

         endif                  ! itmask
         enddo
      endif

!  --- (D3) Now do metal species ---
!     
      if (imchem .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then

C***** CI **********
               scoef = 0._DKIND + 12._DKIND * ( 0._DKIND
     &             + kz15(i) *    HI(i,j,k) *    CH(i,j,k) / 13._DKIND
     &             + kz44(i) *   CII(i,j,k) *    de(i,j,k) / 12._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz20(i) *   H2I(i,j,k) /  2._DKIND
     &             + kz27(i) *    HI(i,j,k)
     &             + kz28(i) *    OH(i,j,k) / 17._DKIND
     &             + kz29(i) *    O2(i,j,k) / 32._DKIND
     &             + kz51(i) *   H2I(i,j,k) /  2._DKIND
               if(igrgr .eq. 1) then
                 acoef = acoef 
     &         + kd5(i) / CI(i,j,k) * 12._DKIND
               endif

            CIp(i)   = ( scoef*dtit(i) + CI(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** CII **********
               scoef = 0._DKIND + 12._DKIND * ( 0._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz37(i) *    OH(i,j,k) / 17._DKIND
     &             + kz38(i) *    O2(i,j,k) / 32._DKIND
     &             + kz44(i) *    de(i,j,k)

            CIIp(i)   = ( scoef*dtit(i) + CII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** CO **********
               scoef = 0._DKIND + 28._DKIND * ( 0._DKIND
     &             + kz28(i) *    CI(i,j,k) *    OH(i,j,k) / 204._DKIND
     &             + kz29(i) *    CI(i,j,k) *    O2(i,j,k) / 384._DKIND
     &             + kz32(i) *    OI(i,j,k) *    CH(i,j,k) / 208._DKIND
     &             + kz38(i) *   CII(i,j,k) *    O2(i,j,k) / 384._DKIND
     &             + kz43(i) *  COII(i,j,k) *    HI(i,j,k) / 28._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz26(i) *    OH(i,j,k) / 17._DKIND

            COp(i)   = ( scoef*dtit(i) + CO(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** CO2 **********
               scoef = 0._DKIND + 44._DKIND * ( 0._DKIND
     &             + kz26(i) *    OH(i,j,k) *    CO(i,j,k) / 476._DKIND
     &            )
               acoef = 0._DKIND

            CO2p(i)   = ( scoef*dtit(i) + CO2(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** OI **********
               scoef = 0._DKIND + 16._DKIND * ( 0._DKIND
     &             + kz17(i) *    HI(i,j,k) *    OH(i,j,k) / 17._DKIND
     &             + kz19(i) *    HI(i,j,k) *    O2(i,j,k) / 32._DKIND
     &             + kz25(i) *    OH(i,j,k) *    OH(i,j,k) / 289._DKIND
     &             + kz29(i) *    CI(i,j,k) *    O2(i,j,k) / 384._DKIND
     &             + kz39(i) *   OII(i,j,k) *    HI(i,j,k) / 16._DKIND
     &             + kz45(i) *   OII(i,j,k) *    de(i,j,k) / 16._DKIND
     &             + kz47(i) * H2OII(i,j,k) *    de(i,j,k) / 18._DKIND
     &             + kz50(i) *  O2II(i,j,k) *    de(i,j,k) / 16._DKIND
     &             + kz53(i) *   SiI(i,j,k) *    O2(i,j,k) / 896._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz21(i) *   H2I(i,j,k) /  2._DKIND
     &             + kz22(i) *   HII(i,j,k)
     &             + kz30(i) *    HI(i,j,k)
     &             + kz31(i) *    OI(i,j,k) / 8._DKIND
     &             + kz32(i) *    CH(i,j,k) / 13._DKIND
     &             + kz33(i) *    OH(i,j,k) / 17._DKIND

            OIp(i)   = ( scoef*dtit(i) + OI(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** OH **********
               scoef = 0._DKIND + 17._DKIND * ( 0._DKIND
     &             + kz18(i) *    HI(i,j,k) *   H2O(i,j,k) / 18._DKIND
     &             + kz19(i) *    HI(i,j,k) *    O2(i,j,k) / 32._DKIND
     &             + kz21(i) *    OI(i,j,k) *   H2I(i,j,k) / 32._DKIND
     &             + kz30(i) *    OI(i,j,k) *    HI(i,j,k) / 16._DKIND
     &             + kz46(i) * H2OII(i,j,k) *    de(i,j,k) / 18._DKIND
     &             + kz49(i) * H3OII(i,j,k) *    de(i,j,k) / 19._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz17(i) *    HI(i,j,k)
     &             + kz24(i) *   H2I(i,j,k) /  2._DKIND
     &             + kz25(i) *    OH(i,j,k) / 8.5_DKIND
     &             + kz26(i) *    CO(i,j,k) / 28._DKIND
     &             + kz28(i) *    CI(i,j,k) / 12._DKIND
     &             + kz33(i) *    OI(i,j,k) / 16._DKIND
     &             + kz34(i) *   HII(i,j,k)
     &             + kz37(i) *   CII(i,j,k) / 12._DKIND
     &             + kz52(i) *   SiI(i,j,k) / 28._DKIND
     &             + kz54(i) *  SiOI(i,j,k) / 44._DKIND

            OHp(i)   = ( scoef*dtit(i) + OH(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** H2O **********
               scoef = 0._DKIND + 18._DKIND * ( 0._DKIND
     &             + kz24(i) *   H2I(i,j,k) *    OH(i,j,k) / 34._DKIND
     &             + kz25(i) *    OH(i,j,k) *    OH(i,j,k) / 289._DKIND
     &             + kz48(i) * H3OII(i,j,k) *    de(i,j,k) / 19._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz18(i) *    HI(i,j,k)
     &             + kz35(i) *   HII(i,j,k)
               if(igrgr .eq. 1) then
                 acoef = acoef 
     &         + kd2(i) / H2O(i,j,k) * 18._DKIND * 3._DKIND
     &         + kd3(i) / H2O(i,j,k) * 18._DKIND * 2._DKIND
     &         + kd4(i) / H2O(i,j,k) * 18._DKIND * 4._DKIND
     &         + kd7(i) / H2O(i,j,k) * 18._DKIND
     &         + kd9(i) / H2O(i,j,k) * 18._DKIND * 3._DKIND
               endif

            H2Op(i)   = ( scoef*dtit(i) + H2O(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** O2 **********
               scoef = 0._DKIND + 32._DKIND * ( 0._DKIND
     &             + kz31(i) *    OI(i,j,k) *    OI(i,j,k) / 256._DKIND
     &             + kz33(i) *    OI(i,j,k) *    OH(i,j,k) / 272._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz19(i) *    HI(i,j,k)
     &             + kz29(i) *    CI(i,j,k) / 12._DKIND
     &             + kz36(i) *   HII(i,j,k)
     &             + kz38(i) *   CII(i,j,k) / 12._DKIND
     &             + kz53(i) *   SiI(i,j,k) / 28._DKIND

            O2p(i)   = ( scoef*dtit(i) + O2(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** SiI **********
               scoef = 0._DKIND + 28._DKIND * ( 0._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz52(i) *    OH(i,j,k) / 17._DKIND
     &             + kz53(i) *    O2(i,j,k) / 32._DKIND
            if (igrgr .eq. 1) then
               acoef = acoef 
     &         + kd0(i) / SiI(i,j,k) * 28._DKIND
            endif

            SiIp(i)   = ( scoef*dtit(i) + SiI(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** SiOI **********
               scoef = 0._DKIND + 44._DKIND * ( 0._DKIND
     &             + kz52(i) *   SiI(i,j,k) *    OH(i,j,k) / 476._DKIND
     &             + kz53(i) *   SiI(i,j,k) *    O2(i,j,k) / 896._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz54(i) *    OH(i,j,k) / 17._DKIND
               if(igrgr .eq. 1) then
                 acoef = acoef 
     &         + kd2(i) / SiOI(i,j,k) * 44._DKIND
     &         + kd3(i) / SiOI(i,j,k) * 44._DKIND
               endif

            SiOIp(i)   = ( scoef*dtit(i) + SiOI(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** SiO2I **********
               scoef = 0._DKIND + 60._DKIND * ( 0._DKIND
     &             + kz54(i) *  SiOI(i,j,k) *    OH(i,j,k) / 748._DKIND
     &            )
               acoef = 0._DKIND
               if(igrgr .eq. 1) then
                 acoef = acoef 
     &         + kd6(i) / SiO2I(i,j,k) * 60._DKIND
               endif

            SiO2Ip(i)   = ( scoef*dtit(i) + SiO2I(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

            ! MINOR BUT IMPORTANT SPECIES FOR MOLECULAR FORMATION
C***** CH **********
               scoef = 0._DKIND + 13._DKIND * ( 0._DKIND
     &             + kz16(i) *    HI(i,j,k) *   CH2(i,j,k) / 14._DKIND
     &             + kz20(i) *    CI(i,j,k) *   H2I(i,j,k) / 24._DKIND
     &             + kz27(i) *    CI(i,j,k) *    HI(i,j,k) / 12._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz15(i) *    HI(i,j,k)
     &             + kz23(i) *   H2I(i,j,k) /  2._DKIND
     &             + kz32(i) *    OI(i,j,k) / 16._DKIND

            CHp(i)   = ( scoef*dtit(i) + CH(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** CH2 **********
               scoef = 0._DKIND + 14._DKIND * ( 0._DKIND
     &             + kz23(i) *   H2I(i,j,k) *    CH(i,j,k) / 26._DKIND
     &             + kz51(i) *   H2I(i,j,k) *    CI(i,j,k) / 24._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz16(i) *    HI(i,j,k)

            CH2p(i)   = ( scoef*dtit(i) + CH2(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** COII **********
               scoef = 0._DKIND + 28._DKIND * ( 0._DKIND
     &             + kz37(i) *   CII(i,j,k) *    OH(i,j,k) / 204._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz43(i) *    HI(i,j,k)

            COIIp(i)   = ( scoef*dtit(i) + COII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** OII **********
               scoef = 0._DKIND + 16._DKIND * ( 0._DKIND
     &             + kz22(i) *   HII(i,j,k) *    OI(i,j,k) / 16._DKIND
     &             + kz38(i) *   CII(i,j,k) *    O2(i,j,k) / 384._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz39(i) *    HI(i,j,k)
     &             + kz40(i) *   H2I(i,j,k) /  2._DKIND
     &             + kz45(i) *    de(i,j,k)

            OIIp(i)   = ( scoef*dtit(i) + OII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** OHII **********
               scoef = 0._DKIND + 17._DKIND * ( 0._DKIND
     &             + kz34(i) *   HII(i,j,k) *    OH(i,j,k) / 17._DKIND
     &             + kz40(i) *   OII(i,j,k) *   H2I(i,j,k) / 32._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz41(i) *   H2I(i,j,k) /  2._DKIND
!!   &             + kz45(i) *    de(i,j,k)

            OHIIp(i)   = ( scoef*dtit(i) + OHII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** H2OII **********
               scoef = 0._DKIND + 18._DKIND * ( 0._DKIND
     &             + kz35(i) *   HII(i,j,k) *   H2O(i,j,k) / 18._DKIND
     &             + kz41(i) *  OHII(i,j,k) *   H2I(i,j,k) / 34._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz42(i) *   H2I(i,j,k) /  2._DKIND
     &             + kz46(i) *    de(i,j,k)
     &             + kz47(i) *    de(i,j,k)

            H2OIIp(i)   = ( scoef*dtit(i) + H2OII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** H3OII **********
               scoef = 0._DKIND + 19._DKIND * ( 0._DKIND
     &             + kz42(i) * H2OII(i,j,k) *   H2I(i,j,k) / 36._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz48(i) *    de(i,j,k)
     &             + kz49(i) *    de(i,j,k)

            H3OIIp(i)   = ( scoef*dtit(i) + H3OII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** O2II **********
               scoef = 0._DKIND + 32._DKIND * ( 0._DKIND
     &             + kz36(i) *   HII(i,j,k) *    O2(i,j,k) / 32._DKIND
     &            )
               acoef = 0._DKIND
     &             + kz50(i) *    de(i,j,k)

            O2IIp(i)   = ( scoef*dtit(i) + O2II(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


         endif
         endif                  ! itmask
         enddo
      endif

!  --- (D4) Now do dust species ---
!     
      if (igrgr .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then

C***** Mg **********
               scoef = 0._DKIND
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kd2(i) / Mg(i,j,k) * 24._DKIND * 2._DKIND
     &         + kd3(i) / Mg(i,j,k) * 24._DKIND
     &         + kd7(i) / Mg(i,j,k) * 24._DKIND

            Mgp(i)   = ( scoef*dtit(i) + Mg(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** Al **********
               scoef = 0._DKIND 
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kd9(i) / Al(i,j,k) * 27._DKIND * 2._DKIND

            Alp(i)   = ( scoef*dtit(i) + Al(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** S  **********
               scoef = 0._DKIND 
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kd8(i) / S(i,j,k) * 32._DKIND

            Sp(i)    = ( scoef*dtit(i) + S(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** Fe **********
               scoef = 0._DKIND 
               acoef = 0._DKIND
                 acoef = acoef 
     &         + kd1(i) / Fe(i,j,k) * 56._DKIND
     &         + kd4(i) / Fe(i,j,k) * 56._DKIND * 3._DKIND
     &         + kd8(i) / Fe(i,j,k) * 56._DKIND

            Fep(i)   = ( scoef*dtit(i) + Fe(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** SiM **********
               scoef = 0._DKIND
                 scoef = scoef 
     &         + kd0(i) * 28._DKIND
               acoef = 0._DKIND

            SiMp(i)   = ( scoef*dtit(i) + SiM(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** FeM **********
               scoef = 0._DKIND
                 scoef = scoef 
     &         + kd1(i) * 56._DKIND
               acoef = 0._DKIND

            FeMp(i)   = ( scoef*dtit(i) + FeM(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** Mg2SiO4 **********
               scoef = 0._DKIND
                 scoef = scoef 
     &         + kd2(i) * 140._DKIND
               acoef = 0._DKIND

            Mg2SiO4p(i)   = ( scoef*dtit(i) + Mg2SiO4(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** MgSiO3 **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd3(i) * 100._DKIND
               acoef = 0._DKIND

            MgSiO3p(i)   = ( scoef*dtit(i) + MgSiO3(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** Fe3O4 **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd4(i) * 232._DKIND
               acoef = 0._DKIND

            Fe3O4p(i)   = ( scoef*dtit(i) + Fe3O4(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** AC **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd5(i) * 12._DKIND
               acoef = 0._DKIND

            ACp(i)   = ( scoef*dtit(i) + AC(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** SiO2D **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd6(i) * 60._DKIND
               acoef = 0._DKIND

            SiO2Dp(i)   = ( scoef*dtit(i) + SiO2D(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** MgO **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd7(i) * 40._DKIND
               acoef = 0._DKIND

            MgOp(i)   = ( scoef*dtit(i) + MgO(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** FeS **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd8(i) * 88._DKIND
               acoef = 0._DKIND

            FeSp(i)   = ( scoef*dtit(i) + FeS(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


C***** Al2O3 **********
               scoef = 0._DKIND 
                 scoef = scoef 
     &         + kd9(i) * 102._DKIND
               acoef = 0._DKIND

            Al2O3p(i)   = ( scoef*dtit(i) + Al2O3(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )


         endif
         endif                  ! itmask
         enddo
      endif
#endif
!   --- (E) Set densities from 1D temps to 3D fields ---

      do i = is+1, ie+1
         if (itmask(i)) then
         HIdot_prev(i) = abs(HI(i,j,k)-HIp(i)) /
     &           max(real(dtit(i), DKIND), tiny8)
         HI(i,j,k)    = max(real(HIp(i), RKIND), tiny)
         HII(i,j,k)   = max(real(HIIp(i), RKIND), tiny)
         HeI(i,j,k)   = max(real(HeIp(i), RKIND), tiny)
         HeII(i,j,k)  = max(real(HeIIp(i), RKIND), tiny)
         HeIII(i,j,k) = max(real(HeIIIp(i), RKIND), 1e-5_RKIND*tiny)

!        de(i,j,k)    = dep(i)

!        Use charge conservation to determine electron fraction

         dedot_prev(i) = de(i,j,k)
         de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND +
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies .gt. 1) 
     &        de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2._RKIND
#ifdef GRACKLE_MD
         if (ispecies .gt. 2)
     &        de(i,j,k) = de(i,j,k) + DII(i,j,k)/2._RKIND
         if (ispecies .gt. 3)
     &        de(i,j,k) = de(i,j,k) - DM(i,j,k)/2._RKIND
     &             + HDII(i,j,k)/3._RKIND + HeHII(i,j,k)/5._RKIND
         if (   (imchem .eq. 1) 
     &    .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) )
     &        de(i,j,k) = de(i,j,k)
     &             + CII(i,j,k)/12._RKIND + COII(i,j,k)/28._RKIND
     &             + OII(i,j,k)/16._RKIND + OHII(i,j,k)/17._RKIND
     &             + H2OII(i,j,k)/18._RKIND + H3OII(i,j,k)/19._RKIND
     &             + O2II(i,j,k)/32._RKIND
#endif
         dedot_prev(i) = abs(de(i,j,k)-dedot_prev(i))/
     &        max(dtit(i),tiny8)

         if (ispecies .gt. 1) then
            HM(i,j,k)    = max(real(HMp(i), RKIND), tiny)
            H2I(i,j,k)   = max(real(H2Ip(i),RKIND), tiny)
            H2II(i,j,k)  = max(real(H2IIp(i), RKIND), tiny)
         endif

         if (ispecies .gt. 2) then
            DI(i,j,k)    = max(real(DIp(i), RKIND), tiny)
            DII(i,j,k)   = max(real(DIIp(i), RKIND), tiny)
            HDI(i,j,k)   = max(real(HDIp(i), RKIND), tiny)
         endif
#ifdef GRACKLE_MD
         if (ispecies .gt. 3) then
            DM(i,j,k)    = max(real(DMp(i), RKIND), tiny)
            HDII(i,j,k)  = max(real(HDIIp(i), RKIND), tiny)
            HeHII(i,j,k) = max(real(HeHIIp(i), RKIND), tiny)
         endif

         if (   (imchem .eq. 1) 
     &    .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) ) then ! do not remove
            CI(i,j,k)      = max(real(CIp(i)     , RKIND), tiny)
            CII(i,j,k)     = max(real(CIIp(i)    , RKIND), tiny)
            CO(i,j,k)      = max(real(COp(i)     , RKIND), tiny)
            CO2(i,j,k)     = max(real(CO2p(i)    , RKIND), tiny)
            OI(i,j,k)      = max(real(OIp(i)     , RKIND), tiny)
            OH(i,j,k)      = max(real(OHp(i)     , RKIND), tiny)
            H2O(i,j,k)     = max(real(H2Op(i)    , RKIND), tiny)
            O2(i,j,k)      = max(real(O2p(i)     , RKIND), tiny)
            SiI(i,j,k)     = max(real(SiIp(i)    , RKIND), tiny)
            SiOI(i,j,k)    = max(real(SiOIp(i)   , RKIND), tiny)
            SiO2I(i,j,k)   = max(real(SiO2Ip(i)  , RKIND), tiny)
            CH(i,j,k)      = max(real(CHp(i)     , RKIND), tiny)
            CH2(i,j,k)     = max(real(CH2p(i)    , RKIND), tiny)
            COII(i,j,k)    = max(real(COIIp(i)   , RKIND), tiny)
            OII(i,j,k)     = max(real(OIIp(i)    , RKIND), tiny)
            OHII(i,j,k)    = max(real(OHIIp(i)   , RKIND), tiny)
            H2OII(i,j,k)   = max(real(H2OIIp(i)  , RKIND), tiny)
            H3OII(i,j,k)   = max(real(H3OIIp(i)  , RKIND), tiny)
            O2II(i,j,k)    = max(real(O2IIp(i)   , RKIND), tiny)
         endif

         if (   (igrgr .eq. 1) 
     &    .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) ) then ! do not remove
            Mg(i,j,k)      = max(real(Mgp(i)     , RKIND), tiny)
            Al(i,j,k)      = max(real(Alp(i)     , RKIND), tiny)
            S(i,j,k)       = max(real(Sp(i)      , RKIND), tiny)
            Fe(i,j,k)      = max(real(Fep(i)     , RKIND), tiny)
            SiM(i,j,k)     = max(real(SiMp(i)    , RKIND), tiny)
            FeM(i,j,k)     = max(real(FeMp(i)    , RKIND), tiny)
            Mg2SiO4(i,j,k) = max(real(Mg2SiO4p(i), RKIND), tiny)
            MgSiO3(i,j,k)  = max(real(MgSiO3p(i) , RKIND), tiny)
            Fe3O4(i,j,k)   = max(real(Fe3O4p(i)  , RKIND), tiny)
            AC(i,j,k)      = max(real(ACp(i)     , RKIND), tiny)
            SiO2D(i,j,k)   = max(real(SiO2Dp(i)  , RKIND), tiny)
            MgO(i,j,k)     = max(real(MgOp(i)    , RKIND), tiny)
            FeS(i,j,k)     = max(real(FeSp(i)    , RKIND), tiny)
            Al2O3(i,j,k)   = max(real(Al2O3p(i)  , RKIND), tiny)
         endif
#endif
      endif                     ! itmask
!     

      if (HI(i,j,k) .ne. HI(i,j,k)) then
#ifdef _OPENMP
!$omp critical
#endif
         write(*,*) 'HUGE HI! :: ', i, j, k, HI(i,j,k)
#ifdef CORRECTED_BY_CHIAKI2020MAY
        write(*,*) 'p1', de(i,j,k), HII(i,j,k)
        write(*,*) 'p1', HeI(i,j,k), HeII(i,j,k), HeIII(i,j,k)
        write(*,*) 'p2', HM(i,j,k), H2I(i,j,k), H2II(i,j,k)
        if (ispecies .gt. 2) then
        write(*,*) 'p3', DI(i,j,k), DII(i,j,k), HDI(i,j,k)
        endif
        if (ispecies .gt. 3) then
        write(*,*) 'p4', DM(i,j,k), HDII(i,j,k), HeHII(i,j,k)
        endif
        if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
        if (imchem .eq. 1) then
        write(*,*) 'm1', CI(i,j,k), CII(i,j,k), CO(i,j,k), CO2(i,j,k)
        write(*,*) 'm2', OI(i,j,k), OH(i,j,k), H2O(i,j,k), O2(i,j,k)
        write(*,*) 'm3', SiI(i,j,k), SiOI(i,j,k), SiO2I(i,j,k)
        write(*,*) 'm4', CH(i,j,k), CH2(i,j,k), COII(i,j,k), OII(i,j,k)
        write(*,*) 'm5', OHII(i,j,k), H2OII(i,j,k)
     &           , H3OII(i,j,k), O2II(i,j,k)
        endif
        if (igrgr .eq. 1) then
        write(*,*) 'm6', Mg(i,j,k), Al(i,j,k), S(i,j,k), Fe(i,j,k)
        write(*,*) 'd1', SiM(i,j,k), FeM(i,j,k), Mg2SiO4(i,j,k)
     &           , MgSiO3(i,j,k), Fe3O4(i,j,k)
        write(*,*) 'd2', AC(i,j,k), SiO2D(i,j,k), MgO(i,j,k)
     &           , FeS(i,j,k), Al2O3(i,j,k)
        endif
        endif
#endif
#ifdef _OPENMP
!$omp end critical
#endif
      endif

      enddo                     ! end loop over i

      return
      end


! ------------------------------------------------------------------
!   This routine correct the highest abundence species to
!     insure conservation of particle number and charge.

      subroutine make_consistent_g(de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, fh, dtoh
#ifdef GRACKLE_MD
     &                      , imchem, igrgr, dom
     &                      , DM, HDII, HeHII
     &                      , CI, CII, CO, CO2
     &                      , OI, OH, H2O, O2
     &                      , SiI, SiOI, SiO2I
     &                      , CH, CH2, COII, OII
     &                      , OHII, H2OII, H3OII, O2II
     &                      , Mg, Al, S, Fe
     &                      , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &                      , AC, SiO2D, MgO, FeS, Al2O3
     &                      , impop3
     &                      , metal_loc, metal_C30, metal_F13
     &                      , loc_XC, loc_XO, loc_XMg, loc_XAl, loc_XSi
     &                      , loc_XS, loc_XFe
     &                      , loc_fC, loc_fO, loc_fMg, loc_fAl, loc_fSi
     &                      , loc_fS, loc_fFe
     &                      , C30_XC, C30_XO, C30_XMg, C30_XAl, C30_XSi
     &                      , C30_XS, C30_XFe
     &                      , C30_fC, C30_fO, C30_fMg, C30_fAl, C30_fSi
     &                      , C30_fS, C30_fFe
     &                      , F13_XC, F13_XO, F13_XMg, F13_XAl, F13_XSi
     &                      , F13_XS, F13_XFe
     &                      , F13_fC, F13_fO, F13_fMg, F13_fAl, F13_fSi
     &                      , F13_fS, F13_fFe
#endif
     &                            )
! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
#ifdef GRACKLE_MD
      integer imchem, igrgr
#endif
      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn), metal(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real*8 fh, dtoh
#ifdef GRACKLE_MD
      real*8 dom
       R_PREC DM(in,jn,kn)   , HDII(in,jn,kn)   , HeHII(in,jn,kn)
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
      integer impop3
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
      real*8  loc_XC , loc_XO, loc_XMg, loc_XAl
     &      , loc_XSi, loc_XS, loc_XFe
      real*8  loc_fC , loc_fO, loc_fMg, loc_fAl
     &      , loc_fSi, loc_fS, loc_fFe
      real*8  C30_XC , C30_XO, C30_XMg, C30_XAl
     &      , C30_XSi, C30_XS, C30_XFe
      real*8  C30_fC , C30_fO, C30_fMg, C30_fAl
     &      , C30_fSi, C30_fS, C30_fFe
      real*8  F13_XC , F13_XO, F13_XMg, F13_XAl
     &      , F13_XSi, F13_XS, F13_XFe
      real*8  F13_fC , F13_fO, F13_fMg, F13_fAl
     &      , F13_fSi, F13_fS, F13_fFe
#endif

!     locals

      integer i, j, k
      real*8 totalH(in), totalHe(in),
     &       totalD, metalfree(in)
      R_PREC correctH, correctHe, correctD
#ifdef GRACKLE_MD
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
      R_PREC tinymetal
      real*8 Ct(in), Ot(in), Mgt(in), Alt(in)
     &     , Sit(in), St(in), Fet(in)
      real*8 Cg(in), Og(in), Mgg(in), Alg(in)
     &     , Sig(in), Sg(in), Feg(in)
      real*8 Cd(in), Od(in), Mgd(in), Ald(in)
     &     , Sid(in), Sd(in), Fed(in)
#endif

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

#ifdef CORRECTED_BY_CHIAKI2020MAY
      if(ispecies .gt. 2) then
         do i = is+1, ie+1
            HDI(i,j,k) = abs(HDI(i,j,k))
            totalH (i) = totalH (i)
     &           + 1._DKIND/3._DKIND*HDI(i,j,k)
         enddo
      endif
#endif
#ifdef GRACKLE_MD
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

      do i = is+1, ie+1
         totalZ = metal_loc(i,j,k)
         if(impop3 .gt. 0) then
            totalZ = totalZ
     &             + metal_C30(i,j,k)
     &             + metal_F13(i,j,k)
         endif
         correctZ = metal(i,j,k) / totalZ
         metal_loc(i,j,k) = metal_loc(i,j,k) * correctZ
         if(impop3 .gt. 0) then
            metal_C30(i,j,k) = metal_C30(i,j,k) * correctZ
            metal_F13(i,j,k) = metal_F13(i,j,k) * correctZ
         endif

          Ct(i) =          loc_XC  * metal_loc(i,j,k)
          Ot(i) =          loc_XO  * metal_loc(i,j,k)
         Mgt(i) =          loc_XMg * metal_loc(i,j,k)
         Alt(i) =          loc_XAl * metal_loc(i,j,k)
         Sit(i) =          loc_XSi * metal_loc(i,j,k)
          St(i) =          loc_XS  * metal_loc(i,j,k)
         Fet(i) =          loc_XFe * metal_loc(i,j,k)

          Ct(i) =  Ct(i) + C30_XC  * metal_C30(i,j,k)
          Ot(i) =  Ot(i) + C30_XO  * metal_C30(i,j,k)
         Mgt(i) = Mgt(i) + C30_XMg * metal_C30(i,j,k)
         Alt(i) = Alt(i) + C30_XAl * metal_C30(i,j,k)
         Sit(i) = Sit(i) + C30_XSi * metal_C30(i,j,k)
          St(i) =  St(i) + C30_XS  * metal_C30(i,j,k)
         Fet(i) = Fet(i) + C30_XFe * metal_C30(i,j,k)

          Ct(i) =  Ct(i) + F13_XC  * metal_F13(i,j,k)
          Ot(i) =  Ot(i) + F13_XO  * metal_F13(i,j,k)
         Mgt(i) = Mgt(i) + F13_XMg * metal_F13(i,j,k)
         Alt(i) = Alt(i) + F13_XAl * metal_F13(i,j,k)
         Sit(i) = Sit(i) + F13_XSi * metal_F13(i,j,k)
          St(i) =  St(i) + F13_XS  * metal_F13(i,j,k)
         Fet(i) = Fet(i) + F13_XFe * metal_F13(i,j,k)

          Cg(i) =          loc_fC  * metal_loc(i,j,k)
          Og(i) =          loc_fO  * metal_loc(i,j,k)
         Mgg(i) =          loc_fMg * metal_loc(i,j,k)
         Alg(i) =          loc_fAl * metal_loc(i,j,k)
         Sig(i) =          loc_fSi * metal_loc(i,j,k)
          Sg(i) =          loc_fS  * metal_loc(i,j,k)
         Feg(i) =          loc_fFe * metal_loc(i,j,k)

          Cg(i) =  Cg(i) + C30_fC  * metal_C30(i,j,k)
          Og(i) =  Og(i) + C30_fO  * metal_C30(i,j,k)
         Mgg(i) = Mgg(i) + C30_fMg * metal_C30(i,j,k)
         Alg(i) = Alg(i) + C30_fAl * metal_C30(i,j,k)
         Sig(i) = Sig(i) + C30_fSi * metal_C30(i,j,k)
          Sg(i) =  Sg(i) + C30_fS  * metal_C30(i,j,k)
         Feg(i) = Feg(i) + C30_fFe * metal_C30(i,j,k)

          Cg(i) =  Cg(i) + F13_fC  * metal_F13(i,j,k)
          Og(i) =  Og(i) + F13_fO  * metal_F13(i,j,k)
         Mgg(i) = Mgg(i) + F13_fMg * metal_F13(i,j,k)
         Alg(i) = Alg(i) + F13_fAl * metal_F13(i,j,k)
         Sig(i) = Sig(i) + F13_fSi * metal_F13(i,j,k)
          Sg(i) =  Sg(i) + F13_fS  * metal_F13(i,j,k)
         Feg(i) = Feg(i) + F13_fFe * metal_F13(i,j,k)

          Cd(i) =  Ct(i) -  Cg(i)
          Od(i) =  Ot(i) -  Og(i)
         Mgd(i) = Mgt(i) - Mgg(i)
         Ald(i) = Alt(i) - Alg(i)
         Sid(i) = Sit(i) - Sig(i)
          Sd(i) =  St(i) -  Sg(i)
         Fed(i) = Fet(i) - Feg(i)
      enddo

      if (imchem .eq. 1) then
         do i = is+1, ie+1
         if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
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
         endif
         enddo
      endif
#endif

!     Correct densities by keeping fractions the same

      do i = is+1, ie+1
#ifdef CORRECTED_BY_CHIAKI2020MAY
         correctH = real(fh*(1._DKIND - dtoh)*metalfree(i)/totalH(i)
     &              , RKIND)
#else
         correctH = real(fh*metalfree(i)/totalH(i), RKIND)
#endif
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
#ifdef GRACKLE_MD
         if(ispecies .gt. 3) then
!!          HDII (i,j,k) = HDII (i,j,k)*correctH
            HeHII(i,j,k) = HeHII(i,j,k)*correctHe
         endif
#endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) +
     &           2._DKIND/3._DKIND*HDI(i,j,k)
#ifdef GRACKLE_MD
            if(ispecies .gt. 3) then
               DM   (i,j,k) = abs(DM   (i,j,k))
               HDII (i,j,k) = abs(HDII (i,j,k))
               totalD = totalD + DM(i,j,k) + 
     &           2._DKIND/3._DKIND*HDII(i,j,k)
            endif
#endif
            correctD = real(fh*dtoh*metalfree(i)/totalD, RKIND)
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
#ifdef GRACKLE_MD
            if(ispecies .gt. 3) then
               DM   (i,j,k) = DM   (i,j,k)*correctD
               HDII (i,j,k) = HDII (i,j,k)*correctD
            endif
#endif
         enddo
      endif

#ifdef GRACKLE_MD
!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (imchem .eq. 1) then
      do i = is+1, ie+1
      if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then

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
         if (igrgr .eq. 1) then
            Mg(i,j,k)      = abs(Mg(i,j,k)     )
            Al(i,j,k)      = abs(Al(i,j,k)     )
            S(i,j,k)       = abs(S(i,j,k)      )
            Fe(i,j,k)      = abs(Fe(i,j,k)     )
            SiM(i,j,k)     = abs(SiM(i,j,k)    )
            FeM(i,j,k)     = abs(FeM(i,j,k)    )
            Mg2SiO4(i,j,k) = abs(Mg2SiO4(i,j,k))
            MgSiO3(i,j,k)  = abs(MgSiO3(i,j,k) )
            Fe3O4(i,j,k)   = abs(Fe3O4(i,j,k)  )
            AC(i,j,k)      = abs(AC(i,j,k)     )
            SiO2D(i,j,k)   = abs(SiO2D(i,j,k)  )
            MgO(i,j,k)     = abs(MgO(i,j,k)    )
            FeS(i,j,k)     = abs(FeS(i,j,k)    )
            Al2O3(i,j,k)   = abs(Al2O3(i,j,k)  )
         endif

!!       if (d(i,j,k)*dom .lt. 1.e6_DKIND) then
         if (d(i,j,k)*dom .lt. 
     &    min(1.e6_DKIND/(metal(i,j,k)/d(i,j,k)/0.02d-4)**2
     &       ,1.e6_DKIND)) then
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
         if (igrgr .eq. 1) then
            totalOd = 64._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 48._DKIND/100._DKIND* MgSiO3(i,j,k)
     &              + 64._DKIND/232._DKIND*  Fe3O4(i,j,k)  
     &              + 32._DKIND/ 60._DKIND*  SiO2D(i,j,k)  
     &              + 16._DKIND/ 40._DKIND*    MgO(i,j,k)    
     &              + 48._DKIND/102._DKIND*  Al2O3(i,j,k)  
            correctOd = real(Od(i)/totalOd, RKIND)
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctOd
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctOd
              Fe3O4(i,j,k) =   Fe3O4(i,j,k)*correctOd
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctOd
                MgO(i,j,k) =     MgO(i,j,k)*correctOd
              Al2O3(i,j,k) =   Al2O3(i,j,k)*correctOd
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
         if (igrgr .eq. 1) then
            totalCd = AC(i,j,k)
            correctCd = real(Cd(i)/totalCd, RKIND)
            AC(i,j,k) = AC(i,j,k)*correctCd
         endif
            
            totalSig =                        SiI(i,j,k)
     &               + 28._DKIND/ 44._DKIND* SiOI(i,j,k)
     &               + 28._DKIND/ 60._DKIND*SiO2I(i,j,k)
            correctSig = real(Sig(i)/totalSig, RKIND)
              SiI(i,j,k) =   SiI(i,j,k)*correctSig
             SiOI(i,j,k) =  SiOI(i,j,k)*correctSig
            SiO2I(i,j,k) = SiO2I(i,j,k)*correctSig
         if (igrgr .eq. 1) then
            totalSid =                          SiM(i,j,k)
     &               + 28._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &               + 28._DKIND/100._DKIND* MgSiO3(i,j,k)
     &               + 28._DKIND/ 60._DKIND*  SiO2D(i,j,k)
            correctSid = real(Sid(i)/totalSid, RKIND)
                SiM(i,j,k) =     SiM(i,j,k)*correctSid
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctSid
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctSid
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctSid
         endif

         if (igrgr .eq. 1) then
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

            totalMgg =                   Mg(i,j,k)
            correctMgg = real( Mgg(i)/totalMgg, RKIND)
            Mg(i,j,k)      = Mg(i,j,k)     *correctMgg

            totalMgd = 48._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &               + 24._DKIND/100._DKIND* MgSiO3(i,j,k) 
     &               + 24._DKIND/ 40._DKIND*    MgO(i,j,k)    
            correctMgd = real( Mgd(i)/totalMgd, RKIND)
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctMgd
            MgSiO3(i,j,k)  = MgSiO3(i,j,k) *correctMgd
            MgO(i,j,k)     = MgO(i,j,k)    *correctMgd

              S(i,j,k) =                         Sg(i) * metal(i,j,k)
            FeS(i,j,k) = 88._DKIND / 32._DKIND * Sd(i) * metal(i,j,k)

               Al(i,j,k) =                         Alg(i) * metal(i,j,k)
            Al2O3(i,j,k) =  102._DKIND/54._DKIND * Ald(i) * metal(i,j,k)
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
         if (igrgr .eq. 1) then
            totalO  = totalO
     &              + 64._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 48._DKIND/100._DKIND* MgSiO3(i,j,k)
     &              + 64._DKIND/232._DKIND*  Fe3O4(i,j,k)  
     &              + 32._DKIND/ 60._DKIND*  SiO2D(i,j,k)  
     &              + 16._DKIND/ 40._DKIND*    MgO(i,j,k)    
     &              + 48._DKIND/102._DKIND*  Al2O3(i,j,k)  
         endif
         if(igrgr .eq. 0) then
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
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctO 
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctO 
              Fe3O4(i,j,k) =   Fe3O4(i,j,k)*correctO 
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctO 
                MgO(i,j,k) =     MgO(i,j,k)*correctO 
              Al2O3(i,j,k) =   Al2O3(i,j,k)*correctO 
         endif
  
            totalC  =                       CI(i,j,k)
     &              +                      CII(i,j,k)
     &              + 12._DKIND/28._DKIND*  CO(i,j,k)
     &              + 12._DKIND/44._DKIND* CO2(i,j,k)
     &              + 12._DKIND/13._DKIND*  CH(i,j,k)
     &              + 12._DKIND/14._DKIND* CH2(i,j,k)
     &              + 12._DKIND/28._DKIND*COII(i,j,k)
         if (igrgr .eq. 1) then
            totalC  = totalC
     &              +                       AC(i,j,k)
         endif
         if (igrgr .eq. 0) then
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
               CI(i,j,k) =   CI(i,j,k)*correctC
              CII(i,j,k) =  CII(i,j,k)*correctC
               CO(i,j,k) =   CO(i,j,k)*correctC
              CO2(i,j,k) =  CO2(i,j,k)*correctC
               CH(i,j,k) =   CH(i,j,k)*correctC
              CH2(i,j,k) =  CH2(i,j,k)*correctC
             COII(i,j,k) = COII(i,j,k)*correctC
               AC(i,j,k) =   AC(i,j,k)*correctC
         endif
            
            totalSi =                       SiI(i,j,k)
     &              + 28._DKIND/ 44._DKIND*SiOI(i,j,k)
     &              + 28._DKIND/ 60._DKIND*SiO2I(i,j,k)
         if (igrgr .eq. 1) then
            totalSi = totalSi
     &              +                          SiM(i,j,k)
     &              + 28._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 28._DKIND/100._DKIND* MgSiO3(i,j,k)
     &              + 28._DKIND/ 60._DKIND*  SiO2D(i,j,k)
         endif
         if (igrgr .eq. 0) then
            correctSi = real(Sig(i)/totalSi, RKIND)
                SiI(i,j,k) =   SiI(i,j,k)*correctSi
               SiOI(i,j,k) =  SiOI(i,j,k)*correctSi
              SiO2I(i,j,k) = SiO2I(i,j,k)*correctSi
         else
            correctSi = real(Sit(i)/totalSi, RKIND)
                SiI(i,j,k) =   SiI(i,j,k)*correctSi
               SiOI(i,j,k) =  SiOI(i,j,k)*correctSi
              SiO2I(i,j,k) = SiO2I(i,j,k)*correctSi
                SiM(i,j,k) =     SiM(i,j,k)*correctSi
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctSi
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctSi
              SiO2D(i,j,k) =   SiO2D(i,j,k)*correctSi
         endif

         if (igrgr .eq. 1) then
            totalFe =                         Fe(i,j,k)     
     &              +                        FeM(i,j,k)  
     &              +168._DKIND/232._DKIND*Fe3O4(i,j,k)
     &              + 56._DKIND/ 88._DKIND*  FeS(i,j,k)  
            correctFe = real( Fet(i)/totalFe, RKIND)
            Fe(i,j,k)    =    Fe(i,j,k)*correctFe
            FeM(i,j,k)   =   FeM(i,j,k)*correctFe
            Fe3O4(i,j,k) = Fe3O4(i,j,k)*correctFe
            FeS(i,j,k)   =   FeS(i,j,k)*correctFe

            totalMg =                           Mg(i,j,k)     
     &              + 48._DKIND/140._DKIND*Mg2SiO4(i,j,k)
     &              + 24._DKIND/100._DKIND* MgSiO3(i,j,k) 
     &              + 24._DKIND/ 40._DKIND*    MgO(i,j,k)    
            correctMg = real( Mgt(i)/totalMg, RKIND)
                 Mg(i,j,k) =      Mg(i,j,k)*correctMg
            Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*correctMg
             MgSiO3(i,j,k) =  MgSiO3(i,j,k)*correctMg
                MgO(i,j,k) =     MgO(i,j,k)*correctMg

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

         endif ! d(i,j,k)

            tinymetal = tiny8 * d(i,j,k)
            CI(i,j,k)      = max(CI(i,j,k), tinymetal)
            CII(i,j,k)     = max(CII(i,j,k), tinymetal)
            CO(i,j,k)      = max(CO(i,j,k), tinymetal)
            CO2(i,j,k)     = max(CO2(i,j,k), tinymetal)
            OI(i,j,k)      = max(OI(i,j,k), tinymetal)
            OH(i,j,k)      = max(OH(i,j,k), tinymetal)
            H2O(i,j,k)     = max(H2O(i,j,k), tinymetal)
            O2(i,j,k)      = max(O2(i,j,k), tinymetal)
            SiI(i,j,k)     = max(SiI(i,j,k), tinymetal)
            SiOI(i,j,k)    = max(SiOI(i,j,k), tinymetal)
            SiO2I(i,j,k)   = max(SiO2I(i,j,k), tinymetal)
            CH(i,j,k)      = max(CH(i,j,k), tinymetal)
            CH2(i,j,k)     = max(CH2(i,j,k), tinymetal)
            COII(i,j,k)    = max(COII(i,j,k), tinymetal)
            OII(i,j,k)     = max(OII(i,j,k), tinymetal)
            OHII(i,j,k)    = max(OHII(i,j,k), tinymetal)
            H2OII(i,j,k)   = max(H2OII(i,j,k), tinymetal)
            H3OII(i,j,k)   = max(H3OII(i,j,k), tinymetal)
            O2II(i,j,k)    = max(O2II(i,j,k), tinymetal)
         if (igrgr .eq. 1) then
            Mg(i,j,k)      = max(Mg(i,j,k), tinymetal)
            Al(i,j,k)      = max(Al(i,j,k), tinymetal)
            S(i,j,k)       = max(S(i,j,k), tinymetal)
            Fe(i,j,k)      = max(Fe(i,j,k), tinymetal)
            SiM(i,j,k)     = max(SiM(i,j,k), tinymetal)
            FeM(i,j,k)     = max(FeM(i,j,k), tinymetal)
            Mg2SiO4(i,j,k) = max(Mg2SiO4(i,j,k), tinymetal)
            MgSiO3(i,j,k)  = max(MgSiO3(i,j,k), tinymetal)
            Fe3O4(i,j,k)   = max(Fe3O4(i,j,k), tinymetal)
            AC(i,j,k)      = max(AC(i,j,k), tinymetal)
            SiO2D(i,j,k)   = max(SiO2D(i,j,k), tinymetal)
            MgO(i,j,k)     = max(MgO(i,j,k), tinymetal)
            FeS(i,j,k)     = max(FeS(i,j,k), tinymetal)
            Al2O3(i,j,k)   = max(Al2O3(i,j,k), tinymetal)
         endif

      endif
      enddo
      endif ! imchem
#endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND +
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &        - HM(i,j,k) + H2II(i,j,k)/2._RKIND
#ifdef GRACKLE_MD
         if (ispecies .gt. 3) de(i,j,k) = de(i,j,k)
     &        - DM   (i,j,k)/2._RKIND
     &        + HDII (i,j,k)/3._RKIND
     &        + HeHII(i,j,k)/5._RKIND
         if ( (imchem .eq. 1)
     &  .and. (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) )
     &      de(i,j,k) = de(i,j,k)
     &        + CII  (i,j,k)/12._RKIND
     &        + COII (i,j,k)/28._RKIND
     &        + OII  (i,j,k)/16._RKIND
     &        + OHII (i,j,k)/17._RKIND
     &        + H2OII(i,j,k)/18._RKIND
     &        + H3OII(i,j,k)/19._RKIND
     &        + O2II (i,j,k)/32._RKIND
#endif
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k

      return
      end
