#include "phys_const.def"

!     Compute grain size increment

      subroutine calc_grain_size_increment_1d(
     &                immulti, imabund, idspecies, igrgr, itmask
     &              , in, jn, kn, is, ie, j, k, dom, d
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , reforg, volorg, H2Oice
     &              , metal, metal_loc
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
      MASK_TYPE itmask(in)
      integer immulti, imabund, idspecies, igrgr
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  SiM(in,jn,kn), FeM(in,jn,kn), Mg2SiO4(in,jn,kn)
     &      , MgSiO3(in,jn,kn), Fe3O4(in,jn,kn), AC(in,jn,kn)
     &      , SiO2D(in,jn,kn), MgO(in,jn,kn), FeS(in,jn,kn)
     &      , Al2O3(in,jn,kn)
     &      , reforg(in,jn,kn), volorg(in,jn,kn), H2Oice(in,jn,kn)
      R_PREC  metal(in,jn,kn)
     &      , metal_loc(in,jn,kn)
     &      , metal_C13(in,jn,kn), metal_C20(in,jn,kn)
     &      , metal_C25(in,jn,kn), metal_C30(in,jn,kn)
     &      , metal_F13(in,jn,kn), metal_F15(in,jn,kn)
     &      , metal_F50(in,jn,kn), metal_F80(in,jn,kn)
     &      , metal_P170(in,jn,kn), metal_P200(in,jn,kn)
     &      , metal_Y19(in,jn,kn)
!     table
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
      integer iSN, nSN, iSN0
      integer SN_i(SN0_N)
      R_PREC  SN_metal(in, SN0_N)
      real*8  SN_fSiM(SN0_N), SN_fFeM(SN0_N)
     &      , SN_fMg2SiO4(SN0_N), SN_fMgSiO3(SN0_N)
     &      , SN_fFe3O4(SN0_N), SN_fAC(SN0_N)
     &      , SN_fSiO2D(SN0_N), SN_fMgO(SN0_N)
     &      , SN_fFeS(SN0_N), SN_fAl2O3(SN0_N)
     &      , SN_freforg(SN0_N)
     &      , SN_fvolorg(SN0_N), SN_fH2Oice(SN0_N)
      real*8  SN_r0SiM(3,SN0_N), SN_r0FeM(3,SN0_N)
     &      , SN_r0Mg2SiO4(3,SN0_N), SN_r0MgSiO3(3,SN0_N)
     &      , SN_r0Fe3O4(3,SN0_N), SN_r0AC(3,SN0_N)
     &      , SN_r0SiO2D(3,SN0_N), SN_r0MgO(3,SN0_N)
     &      , SN_r0FeS(3,SN0_N), SN_r0Al2O3(3,SN0_N)
     &      , SN_r0reforg(3,SN0_N)
     &      , SN_r0volorg(3,SN0_N), SN_r0H2Oice(3,SN0_N)
      real*8  SN_kpSiM(gr_Size,SN0_N), SN_kpFeM(gr_Size,SN0_N)
     &      , SN_kpMg2SiO4(gr_Size,SN0_N), SN_kpMgSiO3(gr_Size,SN0_N)
     &      , SN_kpFe3O4(gr_Size,SN0_N), SN_kpAC(gr_Size,SN0_N)
     &      , SN_kpSiO2D(gr_Size,SN0_N), SN_kpMgO(gr_Size,SN0_N)
     &      , SN_kpFeS(gr_Size,SN0_N), SN_kpAl2O3(gr_Size,SN0_N)
     &      , SN_kpreforg(gr_Size,SN0_N)
     &      , SN_kpvolorg(gr_Size,SN0_N), SN_kpH2Oice(gr_Size,SN0_N)
!     local
      integer i, idx
!!    integer iTd, iTd0

!!    iSN0 = 2
!!    write(*,*) SN0_fMgSiO3(iSN0)
!!    write(*,*) SN0_r0MgSiO3  (:,iSN0)
!!    do iTd = 1, gr_N(2)
!!       iTd0 = (iTd - 1) * gr_N(1)
!!       write(*,*) SN0_kpMgSiO3  (iTd0+1,iSN0)
!!   &            , SN0_kpMgSiO3  (iTd0+2,iSN0)
!!   &            , SN0_kpMgSiO3  (iTd0+3,iSN0)
!!   &            , SN0_kpMgSiO3  (iTd0+4,iSN0)
!!    enddo
!!    stop
        
      ! make arrays
      if (immulti .eq. 0) then
         nSN = 1
         SN_i(1) = imabund + 1

         do i = is+1, ie+1
           SN_metal(i,1) = metal(i,j,k)
         enddo
      else
         nSN = 0
         if(maxval(metal_loc(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  1
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_loc(i,j,k)
            enddo
         endif
         if(maxval(metal_C13(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  2
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_C13(i,j,k)
            enddo
         endif
         if(maxval(metal_C20(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  3
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_C20(i,j,k)
            enddo
         endif
         if(maxval(metal_C25(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  4
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_C25(i,j,k)
            enddo
         endif
         if(maxval(metal_C30(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  5
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_C30(i,j,k)
            enddo
         endif
         if(maxval(metal_F13(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  6
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_F13(i,j,k)
            enddo
         endif
         if(maxval(metal_F15(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  7
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_F15(i,j,k)
            enddo
         endif
         if(maxval(metal_F50(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  8
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_F50(i,j,k)
            enddo
         endif
         if(maxval(metal_F80(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  9
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_F80(i,j,k)
            enddo
         endif
         if(maxval(metal_P170(:,j,k)/metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) = 10
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_P170(i,j,k)
            enddo
         endif
         if(maxval(metal_P200(:,j,k)/metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) = 11
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_P200(i,j,k)
            enddo
         endif
         if(maxval(metal_Y19(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) = 12
            do i = is+1, ie+1
              SN_metal(i,nSN) = metal_Y19(i,j,k)
            enddo
         endif
      endif
         
      do iSN = 1, nSN
         iSN0 = SN_i(iSN)
         if ( idspecies .gt. 0 ) then
            SN_fMgSiO3     (iSN) = SN0_fMgSiO3     (iSN0)
            SN_fAC         (iSN) = SN0_fAC         (iSN0)

            do idx = 1, 3
              SN_r0MgSiO3  (idx,iSN) = SN0_r0MgSiO3  (idx,iSN0)
              SN_r0AC      (idx,iSN) = SN0_r0AC      (idx,iSN0)
            enddo
            do idx = 1, gr_Size
              SN_kpMgSiO3  (idx,iSN) = SN0_kpMgSiO3  (idx,iSN0)
              SN_kpAC      (idx,iSN) = SN0_kpAC      (idx,iSN0)
            enddo
         endif
         if ( idspecies .gt. 1 ) then
            SN_fSiM        (iSN) = SN0_fSiM        (iSN0)
            SN_fFeM        (iSN) = SN0_fFeM        (iSN0)
            SN_fMg2SiO4    (iSN) = SN0_fMg2SiO4    (iSN0)
            SN_fFe3O4      (iSN) = SN0_fFe3O4      (iSN0)
            SN_fSiO2D      (iSN) = SN0_fSiO2D      (iSN0)
            SN_fMgO        (iSN) = SN0_fMgO        (iSN0)
            SN_fFeS        (iSN) = SN0_fFeS        (iSN0)
            SN_fAl2O3      (iSN) = SN0_fAl2O3      (iSN0)

            do idx = 1, 3
              SN_r0SiM     (idx,iSN) = SN0_r0SiM     (idx,iSN0)
              SN_r0FeM     (idx,iSN) = SN0_r0FeM     (idx,iSN0)
              SN_r0Mg2SiO4 (idx,iSN) = SN0_r0Mg2SiO4 (idx,iSN0)
              SN_r0Fe3O4   (idx,iSN) = SN0_r0Fe3O4   (idx,iSN0)
              SN_r0SiO2D   (idx,iSN) = SN0_r0SiO2D   (idx,iSN0)
              SN_r0MgO     (idx,iSN) = SN0_r0MgO     (idx,iSN0)
              SN_r0FeS     (idx,iSN) = SN0_r0FeS     (idx,iSN0)
              SN_r0Al2O3   (idx,iSN) = SN0_r0Al2O3   (idx,iSN0)
            enddo
            do idx = 1, gr_Size
              SN_kpSiM     (idx,iSN) = SN0_kpSiM     (idx,iSN0)
              SN_kpFeM     (idx,iSN) = SN0_kpFeM     (idx,iSN0)
              SN_kpMg2SiO4 (idx,iSN) = SN0_kpMg2SiO4 (idx,iSN0)
              SN_kpFe3O4   (idx,iSN) = SN0_kpFe3O4   (idx,iSN0)
              SN_kpSiO2D   (idx,iSN) = SN0_kpSiO2D   (idx,iSN0)
              SN_kpMgO     (idx,iSN) = SN0_kpMgO     (idx,iSN0)
              SN_kpFeS     (idx,iSN) = SN0_kpFeS     (idx,iSN0)
              SN_kpAl2O3   (idx,iSN) = SN0_kpAl2O3   (idx,iSN0)
            enddo
         endif
         if ( idspecies .gt. 2 ) then
            SN_freforg     (iSN) = SN0_freforg     (iSN0)
            SN_fvolorg     (iSN) = SN0_fvolorg     (iSN0)
            SN_fH2Oice     (iSN) = SN0_fH2Oice     (iSN0)

            do idx = 1, 3
              SN_r0reforg  (idx,iSN) = SN0_r0reforg  (idx,iSN0)
              SN_r0volorg  (idx,iSN) = SN0_r0volorg  (idx,iSN0)
              SN_r0H2Oice  (idx,iSN) = SN0_r0H2Oice  (idx,iSN0)
            enddo
            do idx = 1, gr_Size
              SN_kpreforg  (idx,iSN) = SN0_kpreforg  (idx,iSN0)
              SN_kpvolorg  (idx,iSN) = SN0_kpvolorg  (idx,iSN0)
              SN_kpH2Oice  (idx,iSN) = SN0_kpH2Oice  (idx,iSN0)
            enddo
         endif
      enddo

      !! calculate size increment

      if (idspecies .gt. 0) then
!!    write(*,*) 'MgSiO3'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, MgSiO3  , SN_metal, SN_fMgSiO3  , SN_r0MgSiO3  
     &        , sMgSiO3  , sgMgSiO3  , alMgSiO3  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMgSiO3
     &            )

!!    write(*,*) 'AC'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, AC      , SN_metal, SN_fAC      , SN_r0AC      
     &        , sAC      , sgAC      , alAC      
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpAC
     &            )
      endif

      if (idspecies .gt. 1) then
!!    write(*,*) 'SiM'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, SiM     , SN_metal, SN_fSiM     , SN_r0SiM     
     &        , sSiM     , sgSiM     , alSiM     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpSiM
     &            )

!!    write(*,*) 'FeM'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, FeM     , SN_metal, SN_fFeM     , SN_r0FeM     
     &        , sFeM     , sgFeM     , alFeM     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFeM
     &            )

!!    write(*,*) 'Mg2SiO4'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Mg2SiO4 , SN_metal, SN_fMg2SiO4 , SN_r0Mg2SiO4 
     &        , sMg2SiO4 , sgMg2SiO4 , alMg2SiO4 
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMg2SiO4
     &            )

!!    write(*,*) 'Fe3O4'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Fe3O4   , SN_metal, SN_fFe3O4   , SN_r0Fe3O4   
     &        , sFe3O4   , sgFe3O4   , alFe3O4   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFe3O4
     &            )

!!    write(*,*) 'SiO2D'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, SiO2D   , SN_metal, SN_fSiO2D   , SN_r0SiO2D   
     &        , sSiO2D   , sgSiO2D   , alSiO2D   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpSiO2D
     &            )

!!    write(*,*) 'MgO'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, MgO     , SN_metal, SN_fMgO     , SN_r0MgO     
     &        , sMgO     , sgMgO     , alMgO     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMgO
     &            )

!!    write(*,*) 'FeS'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, FeS     , SN_metal, SN_fFeS     , SN_r0FeS     
     &        , sFeS     , sgFeS     , alFeS     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFeS
     &            )

!!    write(*,*) 'Al2O3'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Al2O3   , SN_metal, SN_fAl2O3   , SN_r0Al2O3   
     &        , sAl2O3   , sgAl2O3   , alAl2O3   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpAl2O3
     &            )
      endif

      if (idspecies .gt. 2) then
!!    write(*,*) 'reforg'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, reforg  , SN_metal, SN_freforg  , SN_r0reforg  
     &        , sreforg  , sgreforg  , alreforg  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpreforg
     &            )

!!    write(*,*) 'volorg'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, volorg  , SN_metal, SN_fvolorg  , SN_r0volorg  
     &        , svolorg  , sgvolorg  , alvolorg  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpvolorg
     &            )

!!    write(*,*) 'H2Oice'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, H2Oice  , SN_metal, SN_fH2Oice  , SN_r0H2Oice  
     &        , sH2Oice  , sgH2Oice  , alH2Oice  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpH2Oice
     &            )
      endif

      do i = is+1, ie+1
      if ( itmask(i) .ne. MASK_FALSE ) then
            
         if (idspecies .gt. 0) then
            sgtot  (i) = sgMgSiO3    (i)
     &                 + sgAC        (i)
         endif
         if (idspecies .gt. 1) then
            sgtot  (i) = sgtot       (i)
     &                 + sgSiM       (i)
     &                 + sgFeM       (i)
     &                 + sgMg2SiO4   (i)
     &                 + sgFe3O4     (i)
     &                 + sgSiO2D     (i)
     &                 + sgMgO       (i)
     &                 + sgFeS       (i)
     &                 + sgAl2O3     (i)
         endif
         if (idspecies .gt. 2) then
            sgtot  (i) = sgtot(i)
     &                 + sgreforg    (i)
     &                 + sgvolorg    (i)
     &                 + sgH2Oice    (i)
         endif
            
         do idx = 1, gr_N(2)
           if (idspecies .gt. 0) then
              altot(idx,i) = alMgSiO3  (idx,i)
     &                     + alAC      (idx,i)
           endif
           if (idspecies .gt. 1) then
              altot(idx,i) = altot     (idx,i)
     &                     + alSiM     (idx,i)
     &                     + alFeM     (idx,i)
     &                     + alMg2SiO4 (idx,i)
     &                     + alFe3O4   (idx,i)
     &                     + alSiO2D   (idx,i)
     &                     + alMgO     (idx,i)
     &                     + alFeS     (idx,i)
     &                     + alAl2O3   (idx,i)
           endif
           if (idspecies .gt. 2) then
              altot(idx,i) = altot     (idx,i)
     &                     + alreforg  (idx,i)
     &                     + alvolorg  (idx,i)
     &                     + alH2Oice  (idx,i)
           endif
         enddo

      endif
      enddo

      return
      end

!=======================================================================
!//////////////  CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D \\\\\\\\\\\\\\\\\\

      subroutine calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, dsp, SN_metal, SN_fsp, SN_r0sp
     &        , ssp, sgsp, alsp
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kp0sp
     &            )

!
!  Calculates properties that are derived from the grain size increment
!  for a single grain species.
!
!  NOTE: the subroutine's name should be more descriptive
!
!  INPUTS:
!    ssp - the density of an individual dust grain species. At the time
!          of writing it has units of g/cm^3
!
!
      implicit NONE

#include "grackle_fortran_types.def"

! input
      integer in, jn, kn, is, ie, j, k
      MASK_TYPE itmask(in)
      integer igrgr
      integer iSN, nSN, SN0_N
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  dsp(in,jn,kn), SN_metal(in, SN0_N)
! table
      real*8  SN_fsp(SN0_N), SN_r0sp(3,SN0_N)
      real*8  ssp
! opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  SN_kp0sp(gr_Size,SN0_N)
! output
      real*8 sgsp(in), alsp(gr_N(2),in)
! local
      integer i
      real*8 coef0, coef1, coef2, coef3
      real*8 dsp_inject_sum
      real*8 dsp0, SN_sgsp, SN_kpsp
      real*8 SN_dsp0(SN0_N), SN_nsp0(SN0_N)
      real*8 drsp(in)
      real*8, parameter :: pi = pi_val
      real*8, parameter :: mh = mass_h
      integer iTd, iTd0
 
      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then

! Step 1: compute the total mass density of the current grain species
!         that was injected (by summing the amounts injected by each
!         injection pathway)
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_dsp0(iSN) = SN_fsp(iSN) * SN_metal(i, iSN)
               endif
            enddo

! Step 2: Compute the size increment for the grain species
!
! Let's go into detail:
! - todo: maybe we move most of this description into the docstring
!   and possibly move some of it into the narrative documentation
! - let's define some notation:
!    - for some time t, let φ(r,t) denote the grain species's
!      differential size distribution at time t, ρ(t) denote the grain
!      species's mass density, and n(t) denote it number density.
!    - we define <rᵖ>(t) as the `p`th order moment of `φ(r,t)`. In
!      other words, <rᵖ>ⱼ = ∫ rᵖ φⱼ dr, where the integral is taken
!      from 0 to infinity
!    - for additional context:
!      - φ(r,t) is normalized such that the p=0 moment is 1
!      - n(t)*φ(r,t) specifies the number density of the grain with
!        radii between r and r+dr
!      - given the mass density of a single grain, ζ, the variables
!        are related via ρ(t) = (4π/3) ζ n(t) <r³>(t).
!    - we define Φⱼ(r) as the grain species's initial differential
!      size distribution when distributed by injection pathway j.
!    - for notational convenience, we define <rᵖ>ⱼ as the `p`th order
!      moment of `Φⱼ(r)`
!    - NOTE: we don't explicitly track φ(r,t) or Φⱼ(r), instead we
!      track moments
! - The underlying model assumes that all injection pathways inject
!   the grain species at a single time t_inj.
!   - In other words:  φ(r,t_inj) = (∑ⱼ nⱼ Φⱼ(r)) / (∑ⱼ nⱼ)
!   - Consequently:
!         <rᵖ>(t_inj) = (∑ⱼ nⱼ <rᵖ>ⱼ) / (∑ⱼ nⱼ)
!         OR
!         <rᵖ>(t_inj) = (∑ⱼ <rᵖ>ⱼ ρⱼ/<r³>ⱼ) / (∑ⱼ ρⱼ/<r³>ⱼ)
! - Furthermore we assume that:
!   - grains can only be created via injection (i.e. the number
!     density of the grain species is frozen)
!   - φ(r,t) can't be deformed. However it can be translated (which
!     is exactly what happens when grains undergo growth).
! - Thus, when we model growth, φ(r,t) = φ(r - δr(t), t_inj) where
!   we call δr(t) the "size increment"
!
! Anybody reading this might notice that several limitations to this
! model. This is discussed in detail in GH Issue 444 (the description
! should be made part of the narrative docs

            if(igrgr .eq. 0) then
! the model is constructed such that δr(t) = 0 if we aren't
! modelling grain growth

            drsp(i) = 0.d0

            else

! here we calculate the current grain size increment, δr(t), using
! conservation of grain number (i.e. grains are only created via
! injection pathways)
! - recall that to model growth we effectively assume that each
!   occurrence of a grain species was injected at a single time, t_inj
! - conservation of grain number tells us that
!   <r³>(t) = <r³>(t_inj) * ρ(t) / ρ(t_inj),
!   where ρ(t_inj) = ∑ⱼ ρⱼ
! - Equation A3 of Chiaki & Wise 2019 indicates that the grain-size
!   increment is the root of the following cubic equation
!    [<r³>(t_inj) - <r³>(t)] + [3 * <r²>(t_inj)] * δr(t)
!       + [3 * <r¹>(t_inj)] * δr²(t) + δr³(t) = 0
! - Let's rewrite the expression (<r³>(t_inj) - <r³>(t)). First,
!   let's plug in the equation from conservation of grain number:
!     <r³>(t_inj) - <r³>(t) = <r³>(t_inj) * [ 1 - (ρ(t)/∑ⱼρⱼ)]
!   Now, let's plug in <r³>(t_inj) = (∑ⱼρⱼ) / (∑ⱼ ρⱼ/<r³>ⱼ), which
!   comes from the generic formula for <rᵖ>(t_inj) defined earlier:
!     <r³>(t_inj) - <r³>(t) = [(∑ⱼρⱼ) - ρ(t)] / (∑ⱼ ρⱼ/<r³>ⱼ)
! - If we plug this into the preceding cubic equation and multiply the
!   resulting equation by (∑ⱼ ρⱼ/<r³>ⱼ), we get:
!       coef0 + coef1 * δr(t) + coef2 * δr²(t) + coef3 * δr³(t) = 0
!   where:
!       - coef0 = (∑ⱼ ρⱼ) - ρ(t)
!       - coef1 = ∑ⱼ 3 <r²>ⱼ ρⱼ/<r³>ⱼ
!       - coef2 = ∑ⱼ 3 <r¹>ⱼ ρⱼ/<r³>ⱼ
!       - coef3 = ∑ⱼ ρⱼ/<r³>ⱼ

! Let's compute these coefficients
            dsp_inject_sum = 0.d0
! We'll use dsp_inject_sum to compute coef0
            coef1 = 0.d0
            coef2 = 0.d0
            coef3 = 0.d0

!  Loop over each injection pathway
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
!  Calculate 4πζnⱼ/3 = ρⱼ/<r³>ⱼ
!  -> recall: that ζ is the mass density of a single grain of the
!     current grain species (i.e. it's a constant)
!  -> TODO: its very confusing that we store the result within SN_nsp0
!     since that variable is later reused to directly hold nⱼ
                  SN_nsp0(iSN) = SN_dsp0(iSN) / SN_r0sp(3,iSN)

                  dsp_inject_sum = dsp_inject_sum + SN_dsp0(iSN)
                  coef1 = coef1 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(2,iSN)
                  coef2 = coef2 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(1,iSN)
                  coef3 = coef3 +        SN_nsp0(iSN)
               endif
            enddo

            coef0 = dsp_inject_sum - dsp(i,j,k)

!  Let's actually solve for the δr(t), the root of the cubic equation
!  - before we do that, we divide both sides of the equation by the
!    value of coef3. This ensures that the rescaled value of coef3 is
!    exactly 1 (a requirement of the solve_cubic_equation routine)
!  - todo: consider not modifying the variables in-place or performing
!    this operation within the solve_cubic_equation
            coef0 = coef0 / coef3
            coef1 = coef1 / coef3
            coef2 = coef2 / coef3

            call solve_cubic_equation(coef2, coef1, coef0, drsp(i))

            drsp(i) = max(drsp(i), 0.d0)

            endif


! Step 3: calculate number density (code_density / g)

            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_nsp0(iSN) = SN_dsp0(iSN)
     &             / (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
               else
                  SN_nsp0(iSN) = 0.d0
               endif
            enddo

! Step 4: calculate geometrical cross-section per unit gas mass
! -> units of cm^2/g
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
            enddo
            sgsp(i) = sgsp(i) / d(i,j,k)
    
! Step 5: calculate optical opacity related quantities
! -> we are effectively constructing a 1d table of values, at various
!    possible grain temperature for each injection pathway (in other
!    words, its a 2D table)
! -> I'm pretty sure that the values we are the opacity (per unit gas
!    mass)
! -> I'm pretty confident that the units are independent of code units
!    (I think the units are cm^2/g)
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
               enddo
               alsp(iTd,i) = alsp(iTd,i) / d(i,j,k)
            enddo


         endif ! itmask
      enddo

      return
      end

!=======================================================================
!////////////////  SUBROUTINE SOLVE_CUBIC_EQUATION \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_cubic_equation(a, b, c, root)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      real*8 a, b, c
      real*8 root
      real*8 root1, root2, root3

!  Locals

      real*8 q, r, m
      real*8 th
      real*8 s,t
      real*8, parameter :: pi = pi_val

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
