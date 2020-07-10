#include "phys_const.def"

!     Compute grain size increment

      subroutine calc_grain_size_increment_1d(
     &                igrgr, itmask
     &              , in, jn, kn, is, ie, j, k, dom, d
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , reforg, volorg, H2Oice
     &              , metal, metal_loc, metal_C30, metal_F13
     &              , loc_fFeM, loc_fMg2SiO4, loc_fMgSiO3, loc_fFeS
     &              , loc_freforg, loc_fvolorg, loc_fH2Oice
     &              , loc_r0FeM, loc_r0Mg2SiO4, loc_r0MgSiO3, loc_r0FeS
     &              , loc_r0reforg, loc_r0volorg, loc_r0H2Oice 
     &              , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &              , C30_fAC, C30_fSiO2D, C30_fMgO
     &              , C30_fFeS, C30_fAl2O3
     &              , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4
     &              , C30_r0MgSiO3, C30_r0AC, C30_r0SiO2D
     &              , C30_r0MgO, C30_r0FeS, C30_r0Al2O3
     &              , F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &              , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fAl2O3
     &              , F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     &              , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0Al2O3
     &              , gr_N, gr_Size, gr_dT, gr_Td
     &              , loc_kpFeM, loc_kpMg2SiO4, loc_kpMgSiO3, loc_kpFeS
     &              , loc_kpreforg, loc_kpvolorg, loc_kpH2Oice
     &              , C30_kpSiM, C30_kpFeM, C30_kpMg2SiO4
     &              , C30_kpMgSiO3, C30_kpAC, C30_kpSiO2D
     &              , C30_kpMgO, C30_kpFeS, C30_kpAl2O3
     &              , F13_kpFeM, F13_kpMg2SiO4, F13_kpMgSiO3
     &              , F13_kpFe3O4, F13_kpAC, F13_kpSiO2D, F13_kpAl2O3
     &              , sgSiM, sgFeM, sgMg2SiO4, sgMgSiO3, sgFe3O4, sgAC
     &              , sgSiO2D, sgMgO, sgFeS, sgAl2O3
     &              , sgreforg, sgvolorg, sgH2Oice, sgtot
     &              , alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC
     &              , alSiO2D, alMgO, alFeS, alAl2O3
     &              , alreforg, alvolorg, alH2Oice, altot
     &            )

      implicit NONE
#include "grackle_fortran_types.def"
#include "dust_const.def"

!     in
      integer in, jn, kn, is, ie, j, k
      logical itmask(in)
      integer igrgr
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  SiM(in,jn,kn), FeM(in,jn,kn), Mg2SiO4(in,jn,kn)
     &      , MgSiO3(in,jn,kn), Fe3O4(in,jn,kn), AC(in,jn,kn)
     &      , SiO2D(in,jn,kn), MgO(in,jn,kn), FeS(in,jn,kn)
     &      , Al2O3(in,jn,kn)
     &      , reforg(in,jn,kn), volorg(in,jn,kn), H2Oice(in,jn,kn)
      R_PREC  metal(in,jn,kn)
     &      , metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
!     table
      real*8  loc_fFeM, loc_fMg2SiO4, loc_fMgSiO3, loc_fFeS
     &      , loc_freforg, loc_fvolorg, loc_fH2Oice
      real*8  loc_r0FeM(3), loc_r0Mg2SiO4(3), loc_r0MgSiO3(3)
     &      , loc_r0FeS(3), loc_r0reforg(3), loc_r0volorg(3)
     &      , loc_r0H2Oice(3)
      real*8  C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &      , C30_fAC, C30_fSiO2D, C30_fMgO
     &      , C30_fFeS, C30_fAl2O3
      real*8  C30_r0SiM(3), C30_r0FeM(3), C30_r0Mg2SiO4(3)
     &      , C30_r0MgSiO3(3), C30_r0AC(3), C30_r0SiO2D(3)
     &      , C30_r0MgO(3), C30_r0FeS(3), C30_r0Al2O3(3)
      real*8  F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &      , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fAl2O3
      real*8  F13_r0FeM(3), F13_r0Mg2SiO4(3), F13_r0MgSiO3(3)
     &      , F13_r0Fe3O4(3), F13_r0AC(3), F13_r0SiO2D(3)
     &      , F13_r0Al2O3(3)
!     opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  loc_kpFeM(gr_Size), loc_kpMg2SiO4(gr_Size)
     &      , loc_kpMgSiO3(gr_Size), loc_kpFeS(gr_Size)
     &      , loc_kpreforg(gr_Size), loc_kpvolorg(gr_Size)
     &      , loc_kpH2Oice(gr_Size)
      real*8  C30_kpSiM(gr_Size), C30_kpFeM(gr_Size)
     &      , C30_kpMg2SiO4(gr_Size), C30_kpMgSiO3(gr_Size)
     &      , C30_kpFe3O4(gr_Size), C30_kpAC(gr_Size)
     &      , C30_kpSiO2D(gr_Size), C30_kpMgO(gr_Size)
     &      , C30_kpFeS(gr_Size), C30_kpAl2O3(gr_Size)
      real*8  F13_kpFeM(gr_Size)
     &      , F13_kpMg2SiO4(gr_Size), F13_kpMgSiO3(gr_Size)
     &      , F13_kpFe3O4(gr_Size), F13_kpAC(gr_Size)
     &      , F13_kpSiO2D(gr_Size), F13_kpAl2O3(gr_Size)
!     out
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
!     array
      integer iSN, nSN, nSNmax
      parameter (nSNmax = 10)
      R_PREC  SN_metal(in, nSNmax)
      real*8  SN_fSiM(nSNmax), SN_fFeM(nSNmax)
     &      , SN_fMg2SiO4(nSNmax), SN_fMgSiO3(nSNmax)
     &      , SN_fFe3O4(nSNmax), SN_fAC(nSNmax)
     &      , SN_fSiO2D(nSNmax), SN_fMgO(nSNmax)
     &      , SN_fFeS(nSNmax), SN_fAl2O3(nSNmax)
     &      , SN_freforg(nSNmax)
     &      , SN_fvolorg(nSNmax), SN_fH2Oice(nSNmax)
      real*8  SN_r0SiM(3,nSNmax), SN_r0FeM(3,nSNmax)
     &      , SN_r0Mg2SiO4(3,nSNmax), SN_r0MgSiO3(3,nSNmax)
     &      , SN_r0Fe3O4(3,nSNmax), SN_r0AC(3,nSNmax)
     &      , SN_r0SiO2D(3,nSNmax), SN_r0MgO(3,nSNmax)
     &      , SN_r0FeS(3,nSNmax), SN_r0Al2O3(3,nSNmax)
     &      , SN_r0reforg(3,nSNmax)
     &      , SN_r0volorg(3,nSNmax), SN_r0H2Oice(3,nSNmax)
      real*8  SN_kpSiM(gr_Size,nSNmax), SN_kpFeM(gr_Size,nSNmax)
     &      , SN_kpMg2SiO4(gr_Size,nSNmax), SN_kpMgSiO3(gr_Size,nSNmax)
     &      , SN_kpFe3O4(gr_Size,nSNmax), SN_kpAC(gr_Size,nSNmax)
     &      , SN_kpSiO2D(gr_Size,nSNmax), SN_kpMgO(gr_Size,nSNmax)
     &      , SN_kpFeS(gr_Size,nSNmax), SN_kpAl2O3(gr_Size,nSNmax)
     &      , SN_kpreforg(gr_Size,nSNmax)
     &      , SN_kpvolorg(gr_Size,nSNmax), SN_kpH2Oice(gr_Size,nSNmax)
!     local
      integer i
        
      ! make arrays
      nSN = 3

      SN_metal(:,1) = metal_loc(:,j,k)
      SN_metal(:,2) = metal_C30(:,j,k)
      SN_metal(:,3) = metal_F13(:,j,k)

      SN_fSiM        (1) = 0.0
      SN_fFeM        (1) = loc_fFeM     
      SN_fMg2SiO4    (1) = loc_fMg2SiO4 
      SN_fMgSiO3     (1) = loc_fMgSiO3  
      SN_fFe3O4      (1) = 0.0
      SN_fAC         (1) = 0.0
      SN_fSiO2D      (1) = 0.0
      SN_fMgO        (1) = 0.0
      SN_fFeS        (1) = loc_fFeS     
      SN_fAl2O3      (1) = 0.0
      SN_freforg     (1) = loc_freforg  
      SN_fvolorg     (1) = loc_fvolorg  
      SN_fH2Oice     (1) = loc_fH2Oice  
      SN_r0SiM     (:,1) = 0.0
      SN_r0FeM     (:,1) = loc_r0FeM     (:)
      SN_r0Mg2SiO4 (:,1) = loc_r0Mg2SiO4 (:)
      SN_r0MgSiO3  (:,1) = loc_r0MgSiO3  (:)
      SN_r0Fe3O4   (:,1) = 0.0
      SN_r0AC      (:,1) = 0.0
      SN_r0SiO2D   (:,1) = 0.0
      SN_r0MgO     (:,1) = 0.0
      SN_r0FeS     (:,1) = loc_r0FeS     (:)
      SN_r0Al2O3   (:,1) = 0.0
      SN_r0reforg  (:,1) = loc_r0reforg  (:)
      SN_r0volorg  (:,1) = loc_r0volorg  (:)
      SN_r0H2Oice  (:,1) = loc_r0H2Oice  (:)
      SN_kpSiM     (:,1) = 0.0
      SN_kpFeM     (:,1) = loc_kpFeM     (:)
      SN_kpMg2SiO4 (:,1) = loc_kpMg2SiO4 (:)
      SN_kpMgSiO3  (:,1) = loc_kpMgSiO3  (:)
      SN_kpFe3O4   (:,1) = 0.0
      SN_kpAC      (:,1) = 0.0
      SN_kpSiO2D   (:,1) = 0.0
      SN_kpMgO     (:,1) = 0.0
      SN_kpFeS     (:,1) = loc_kpFeS     (:)
      SN_kpAl2O3   (:,1) = 0.0
      SN_kpreforg  (:,1) = loc_kpreforg  (:)
      SN_kpvolorg  (:,1) = loc_kpvolorg  (:)
      SN_kpH2Oice  (:,1) = loc_kpH2Oice  (:)

      SN_fSiM        (2) = C30_fSiM     
      SN_fFeM        (2) = C30_fFeM     
      SN_fMg2SiO4    (2) = C30_fMg2SiO4 
      SN_fMgSiO3     (2) = C30_fMgSiO3  
      SN_fFe3O4      (2) = 0.0
      SN_fAC         (2) = C30_fAC      
      SN_fSiO2D      (2) = C30_fSiO2D   
      SN_fMgO        (2) = C30_fMgO     
      SN_fFeS        (2) = C30_fFeS     
      SN_fAl2O3      (2) = C30_fAl2O3   
      SN_freforg     (2) = 0.0
      SN_fvolorg     (2) = 0.0
      SN_fH2Oice     (2) = 0.0
      SN_r0SiM     (:,2) = C30_r0SiM     (:)
      SN_r0FeM     (:,2) = C30_r0FeM     (:)
      SN_r0Mg2SiO4 (:,2) = C30_r0Mg2SiO4 (:)
      SN_r0MgSiO3  (:,2) = C30_r0MgSiO3  (:)
      SN_r0Fe3O4   (:,2) = 0.0
      SN_r0AC      (:,2) = C30_r0AC      (:)
      SN_r0SiO2D   (:,2) = C30_r0SiO2D   (:)
      SN_r0MgO     (:,2) = C30_r0MgO     (:)
      SN_r0FeS     (:,2) = C30_r0FeS     (:)
      SN_r0Al2O3   (:,2) = C30_r0Al2O3   (:)
      SN_r0reforg  (:,2) = 0.0
      SN_r0volorg  (:,2) = 0.0
      SN_r0H2Oice  (:,2) = 0.0
      SN_kpSiM     (:,2) = C30_kpSiM     (:)
      SN_kpFeM     (:,2) = C30_kpFeM     (:)
      SN_kpMg2SiO4 (:,2) = C30_kpMg2SiO4 (:)
      SN_kpMgSiO3  (:,2) = C30_kpMgSiO3  (:)
      SN_kpFe3O4   (:,2) = 0.0
      SN_kpAC      (:,2) = C30_kpAC      (:)
      SN_kpSiO2D   (:,2) = C30_kpSiO2D   (:)
      SN_kpMgO     (:,2) = C30_kpMgO     (:)
      SN_kpFeS     (:,2) = C30_kpFeS     (:)
      SN_kpAl2O3   (:,2) = C30_kpAl2O3   (:)
      SN_kpreforg  (:,2) = 0.0
      SN_kpvolorg  (:,2) = 0.0
      SN_kpH2Oice  (:,2) = 0.0

      SN_fSiM        (3) = 0.0    
      SN_fFeM        (3) = F13_fFeM     
      SN_fMg2SiO4    (3) = F13_fMg2SiO4 
      SN_fMgSiO3     (3) = F13_fMgSiO3  
      SN_fFe3O4      (3) = F13_fFe3O4   
      SN_fAC         (3) = F13_fAC      
      SN_fSiO2D      (3) = F13_fSiO2D   
      SN_fMgO        (3) = 0.0     
      SN_fFeS        (3) = 0.0     
      SN_fAl2O3      (3) = F13_fAl2O3   
      SN_freforg     (3) = 0.0
      SN_fvolorg     (3) = 0.0
      SN_fH2Oice     (3) = 0.0
      SN_r0SiM     (:,3) = 0.0
      SN_r0FeM     (:,3) = F13_r0FeM     (:)
      SN_r0Mg2SiO4 (:,3) = F13_r0Mg2SiO4 (:)
      SN_r0MgSiO3  (:,3) = F13_r0MgSiO3  (:)
      SN_r0Fe3O4   (:,3) = F13_r0Fe3O4   (:)
      SN_r0AC      (:,3) = F13_r0AC      (:)
      SN_r0SiO2D   (:,3) = F13_r0SiO2D   (:)
      SN_r0MgO     (:,3) = 0.0
      SN_r0FeS     (:,3) = 0.0
      SN_r0Al2O3   (:,3) = F13_r0Al2O3   (:)
      SN_r0reforg  (:,3) = 0.0
      SN_r0volorg  (:,3) = 0.0
      SN_r0H2Oice  (:,3) = 0.0
      SN_kpSiM     (:,3) = 0.0
      SN_kpFeM     (:,3) = F13_kpFeM     (:)
      SN_kpMg2SiO4 (:,3) = F13_kpMg2SiO4 (:)
      SN_kpMgSiO3  (:,3) = F13_kpMgSiO3  (:)
      SN_kpFe3O4   (:,3) = F13_kpFe3O4   (:)
      SN_kpAC      (:,3) = F13_kpAC      (:)
      SN_kpSiO2D   (:,3) = F13_kpSiO2D   (:)
      SN_kpMgO     (:,3) = 0.0
      SN_kpFeS     (:,3) = 0.0
      SN_kpAl2O3   (:,3) = F13_kpAl2O3   (:)
      SN_kpreforg  (:,3) = 0.0
      SN_kpvolorg  (:,3) = 0.0
      SN_kpH2Oice  (:,3) = 0.0

      !! calculate size increment

!!    write(*,*) 'SiM'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, SiM     , SN_metal, SN_fSiM     , SN_r0SiM     
     &        , sSiM     , sgSiM     , alSiM     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpSiM
     &            )

!!    write(*,*) 'FeM'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, FeM     , SN_metal, SN_fFeM     , SN_r0FeM     
     &        , sFeM     , sgFeM     , alFeM     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFeM
     &            )

!!    write(*,*) 'Mg2SiO4'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Mg2SiO4 , SN_metal, SN_fMg2SiO4 , SN_r0Mg2SiO4 
     &        , sMg2SiO4 , sgMg2SiO4 , alMg2SiO4 
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMg2SiO4
     &            )

!!    write(*,*) 'MgSiO3'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, MgSiO3  , SN_metal, SN_fMgSiO3  , SN_r0MgSiO3  
     &        , sMgSiO3  , sgMgSiO3  , alMgSiO3  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMgSiO3
     &            )

!!    write(*,*) 'Fe3O4'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Fe3O4   , SN_metal, SN_fFe3O4   , SN_r0Fe3O4   
     &        , sFe3O4   , sgFe3O4   , alFe3O4   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFe3O4
     &            )

!!    write(*,*) 'AC'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, AC      , SN_metal, SN_fAC      , SN_r0AC      
     &        , sAC      , sgAC      , alAC      
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpAC
     &            )

!!    write(*,*) 'SiO2D'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, SiO2D   , SN_metal, SN_fSiO2D   , SN_r0SiO2D   
     &        , sSiO2D   , sgSiO2D   , alSiO2D   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpSiO2D
     &            )

!!    write(*,*) 'MgO'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, MgO     , SN_metal, SN_fMgO     , SN_r0MgO     
     &        , sMgO     , sgMgO     , alMgO     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMgO
     &            )

!!    write(*,*) 'FeS'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, FeS     , SN_metal, SN_fFeS     , SN_r0FeS     
     &        , sFeS     , sgFeS     , alFeS     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFeS
     &            )

!!    write(*,*) 'Al2O3'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Al2O3   , SN_metal, SN_fAl2O3   , SN_r0Al2O3   
     &        , sAl2O3   , sgAl2O3   , alAl2O3   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpAl2O3
     &            )

!!    write(*,*) 'reforg'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, reforg  , SN_metal, SN_freforg  , SN_r0reforg  
     &        , sreforg  , sgreforg  , alreforg  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpreforg
     &            )

!!    write(*,*) 'volorg'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, volorg  , SN_metal, SN_fvolorg  , SN_r0volorg  
     &        , svolorg  , sgvolorg  , alvolorg  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpvolorg
     &            )

!!    write(*,*) 'H2Oice'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, H2Oice  , SN_metal, SN_fH2Oice  , SN_r0H2Oice  
     &        , sH2Oice  , sgH2Oice  , alH2Oice  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpH2Oice
     &            )

      do i = is+1, ie+1
         if ( itmask(i) ) then
            
            sgtot  (i) = sgSiM       (i)
     &                 + sgFeM       (i)
     &                 + sgMg2SiO4   (i)
     &                 + sgMgSiO3    (i)
     &                 + sgFe3O4     (i)
     &                 + sgAC        (i)
     &                 + sgSiO2D     (i)
     &                 + sgMgO       (i)
     &                 + sgFeS       (i)
     &                 + sgAl2O3     (i)
     &                 + sgreforg    (i)
     &                 + sgvolorg    (i)
     &                 + sgH2Oice    (i)
            
            altot(:,i) = alSiM     (:,i)
     &                 + alFeM     (:,i)
     &                 + alMg2SiO4 (:,i)
     &                 + alMgSiO3  (:,i)
     &                 + alFe3O4   (:,i)
     &                 + alAC      (:,i)
     &                 + alSiO2D   (:,i)
     &                 + alMgO     (:,i)
     &                 + alFeS     (:,i)
     &                 + alAl2O3   (:,i)
     &                 + alreforg  (:,i)
     &                 + alvolorg  (:,i)
     &                 + alH2Oice  (:,i)

         endif
      enddo

      return
      end

!=======================================================================
!//////////////  CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D \\\\\\\\\\\\\\\\\\

      subroutine calc_grain_size_increment_species_1d(
     &          igrgr, itmask
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, dsp, SN_metal, SN_fsp, SN_r0sp
     &        , ssp, sgsp, alsp
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kp0sp
     &            )


      implicit NONE
! input
      integer in, jn, kn, is, ie, j, k
      logical itmask(in)
      integer igrgr
      integer iSN, nSN, nSNmax
      parameter (nSNmax = 10)
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  dsp(in,jn,kn), SN_metal(in, nSNmax)
! table
      real*8  SN_fsp(nSNmax), SN_r0sp(3,nSNmax)
      real*8  ssp
! opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  SN_kp0sp(gr_Size,nSNmax)
! output
      real*8 sgsp(in), alsp(gr_N(2),in)
! local
      integer i
      real*8 coef0, coef1, coef2, coef3
      real*8 dsp0, SN_sgsp, SN_kpsp
      real*8 SN_dsp0(nSNmax), SN_nsp0(nSNmax)
      real*8 drsp(in)
      real*8 pi, mh
      parameter (pi = pi_val)
      parameter (mh = mass_h)
!     debug
      real*8 SN_dsp(nSNmax), SN_msp(nSNmax), dsp1
      integer iTd, iTd0

!       write(*,*) 
!    &             SN_r0sp(1,1)
!    &           , SN_r0sp(2,1)
!    &           , SN_r0sp(3,1)
!    &           , SN_r0sp(4,1)
!     do iTd = 1, gr_N(2)
!       iTd0 = (iTd-1)*gr_N(1)
!       write(*,*) gr_Td(iTd)
!    &           , SN_kp0sp(iTd0+1,1)
!    &           , SN_kp0sp(iTd0+2,1)
!    &           , SN_kp0sp(iTd0+3,1)
!    &           , SN_kp0sp(iTd0+4,1)
!     enddo

      do i = is+1, ie+1
         if ( itmask(i) ) then

         if( dsp(i,j,k) .gt. 1.d-15*d(i,j,k) ) then

            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_dsp0(iSN) = SN_fsp(iSN) * SN_metal(i, iSN)
!!                write(*,*) iSN, SN_fsp(iSN),
!!   &                             SN_metal(i, iSN) 
               endif
            enddo

            if(igrgr .eq. 0) then

            drsp(i) = 0.d0

            else

!!! TEST !!!!!! TEST !!!!!! TEST !!!!!! TEST !!!!!! TEST !!!!!! TEST !!!
!!!         dsp(i,j,k) = 3.d0 * dsp(i,j,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! calculate grain size increment from the conservation of mass

            coef0 = 0.d0
            coef1 = 0.d0
            coef2 = 0.d0
            coef3 = 0.d0
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_nsp0(iSN) = SN_dsp0(iSN) / SN_r0sp(3,iSN)
                  coef0 = coef0 +        SN_dsp0(iSN)
                  coef1 = coef1 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(2,iSN)
                  coef2 = coef2 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(1,iSN)
                  coef3 = coef3 +        SN_nsp0(iSN)
               endif
            enddo
!!          write(*,*) coef0, dsp(i,j,k)
            coef0 = coef0 - dsp(i,j,k)
!!          write(*,*) coef0, coef1, coef2, coef3

            coef0 = coef0 / coef3
            coef1 = coef1 / coef3
            coef2 = coef2 / coef3
            
            call solve_cubic_equation(coef2, coef1, coef0, drsp(i))
!!          write(*,*) coef0, coef1, coef2, coef3, drsp(i)
            drsp(i) = max(drsp(i), 0.d0)

            endif


! calculate number density (code_density / g)

            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_nsp0(iSN) = SN_dsp0(iSN)
     &             / (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
!                 if(iSN.eq.1)
!    &              write(*,*) 'aa', SN_dsp0(iSN), dsp(i,j,k)
!    &             , (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
!    &                      , SN_nsp0(iSN)
               else
                  SN_nsp0(iSN) = 0.d0
               endif
            enddo

! calculate geometrical cross-section per unit gas mass
            sgsp(i) = 0.d0
            do iSN = 1, nSN
               if( SN_fsp(iSN) .gt. 0.d0) then
                  SN_sgsp = pi *
     &               (        SN_r0sp(2,iSN)
     &               + 2.d0 * SN_r0sp(1,iSN) * drsp(i)
     &               +                         drsp(i)**2
     &               )
               else
                  SN_sgsp = 0.d0
               endif
               sgsp(i) = sgsp(i) + SN_nsp0(iSN) * SN_sgsp
!!             if(iSN.eq.1)
!!   &           write(*,*) dsp(i,j,k)/d(i,j,k)
!!   &                  , SN_dsp0(iSN)/d(i,j,k)
!!   &              , ssp, SN_r0sp(3, iSN), SN_r0sp(2,iSN) 
!!   &                  , drsp(i), sgsp(i)/d(i,j,k)
            enddo
            sgsp(i) = sgsp(i) / d(i,j,k)
    
! calculate optical cross-section per unit gas mass
            do iTd = 1, gr_N(2)
               iTd0 = (iTd-1)*gr_N(1)
               alsp(iTd,i) = 0.d0
               do iSN = 1, nSN
                  if( SN_fsp(iSN) .gt. 0.d0) then
                     SN_kpsp = 4.d0*pi/3.d0 * ssp *
     &                  (        SN_kp0sp(iTd0+4,iSN)
     &                  + 3.d0 * SN_kp0sp(iTd0+3,iSN) * drsp(i)
     &                  + 3.d0 * SN_kp0sp(iTd0+2,iSN) * drsp(i)**2
     &                  +        SN_kp0sp(iTd0+1,iSN) * drsp(i)**3
     &                  )
                  else
                     SN_kpsp = 0.d0
                  endif
                  alsp(iTd,i) = alsp(iTd,i) + SN_nsp0(iSN) * SN_kpsp
!                 if(iSN.eq.1)
!    &               write(*,*) gr_Td(iTd)
!    &                        , SN_dsp0(iSN)/d(i,j,k)
!    &             , (4.d0*pi/3.d0 * ssp * SN_kp0sp(iTd0+4,iSN))
!    &             / (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
               enddo
               alsp(iTd,i) = alsp(iTd,i) / d(i,j,k)
            enddo

! decompose grain mass density for each SN model !!!!!!!!!!!!!!!!!!!!!!!
!           dsp0 = 0.d0
!           dsp1 = 0.d0
!           do iSN = 1, nSN
!              dsp0 = dsp0 + SN_dsp0(iSN)
!              if( SN_fsp(iSN) .gt. 0.d0) then
!                 SN_dsp(iSN) = 4.d0*pi/3.d0 * ssp * (
!    &                        SN_nsp0(iSN) * SN_r0sp(3,iSN)
!    &               + 3.d0 * SN_nsp0(iSN) * SN_r0sp(2,iSN) * drsp(i)
!    &               + 3.d0 * SN_nsp0(iSN) * SN_r0sp(1,iSN) * drsp(i)**2
!    &               +        SN_nsp0(iSN)                  * drsp(i)**3
!    &               )
!                 SN_msp(iSN) = 4.d0*pi/3.d0 * ssp * (
!    &                        SN_r0sp(3,iSN)
!    &               + 3.d0 * SN_r0sp(2,iSN) * drsp(i)
!    &               + 3.d0 * SN_r0sp(1,iSN) * drsp(i)**2
!    &               +                         drsp(i)**3
!    &               )
!!                if(iSN.eq.2)
!!   &            write(*,*) iSN, (SN_dsp(iSN)/SN_msp(iSN))/SN_nsp0(iSN)
!                 dsp1 = dsp1 + SN_dsp(iSN)
!              endif
!           enddo
!           write(*,*) drsp(i), dsp1/dsp(i,j,k)

         else

            sgsp  (i) = tiny8
            alsp(:,i) = tiny8

         endif


         endif ! itmask
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
      parameter (pi = pi_val)

      q = (a*a - 3.d0*b)/9.d0
      r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c)/54.d0
      m = r*r - q*q*q
    
      if (m .lt. 0.d0) then !! three real roots 
        th = acos( r / sqrt(q*q*q) )
        root1 = -(2.d0*sqrt(q)*cos( th         /3.d0))-a/3.d0
        root2 = -(2.d0*sqrt(q)*cos((th+2.d0*pi)/3.d0))-a/3.d0
        root3 = -(2.d0*sqrt(q)*cos((th-2.d0*pi)/3.d0))-a/3.d0
        write(*,*) "three real roots", root1, root2, root3
        stop
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
