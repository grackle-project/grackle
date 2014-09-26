#include "fortran.def"
!=======================================================================
!////////////////////  SUBROUTINE COOL_MULTI_TIME_G  \\\\\\\\\\\\\\\\\\\

      subroutine cool_multi_time_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                  cooltime,
     &                in, jn, kn, nratec, iexpand, 
     &                ispecies, imetal, imcool, idust, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh, z_solar,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                gasgra, 
     &                ih2optical, iciecool, ciecoa, 
     &                icmbTfloor, iClHeat,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3,
     &                priDataSize, priCooling, priHeating,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3,
     &                metDataSize, metCooling, metHeating,
     &                mutaba)

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

!  Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, nratec, 
     &        iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        imcool, idust, igammah, ih2optical, iciecool
      real*8  aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma, fh, z_solar
      R_PREC  d(in,jn,kn),   e(in,jn,kn),
     &        u(in,jn,kn),   v(in,jn,kn),    w(in,jn,kn),
     &        de(in,jn,kn),  HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn),  H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn),  HDI(in,jn,kn),
     &        metal(in,jn,kn),
     &        cooltime(in,jn,kn)
      real*8  hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     &        gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), gasgra(nratec),
     &        ciecoa(nratec),
     &        ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec),
     &        mutaba(nratec),
     &        compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        gammaha

!  Cloudy cooling data

      integer icmbTfloor, iClHeat
      integer*8 priGridRank, priDataSize,
     &     metGridRank, metDataSize,
     &     priGridDim(priGridRank), metGridDim(metGridRank)
      real*8 priPar1(priGridDim(1)), priPar2(priGridDim(2)), 
     &     priPar3(priGridDim(3)),
     &     metPar1(metGridDim(1)), metPar2(metGridDim(2)), 
     &     metPar3(metGridDim(3)),
     &     priCooling(priDataSize), priHeating(priDataSize),
     &     metCooling(metDataSize), metHeating(metDataSize)

!  Parameters

      integer ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)

!  Locals

      integer i, j, k
      real*8 comp1, comp2, energy

!  Slice locals
 
      integer*8 indixe(ijk)
      real*8 t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &       tgas(ijk), tgasold(ijk),
     &       tdust(ijk), metallicity(ijk), rhoH(ijk)
      real*8 edot(ijk)

!  Cooling/heating slice locals

      real*8 ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &       ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &       reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &       brem(ijk), cieco(ijk),
     &       hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &       gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk)
          
!  Iteration mask for multi_cool

      logical itmask(ijk)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Convert densities from comoving to 'proper'

      if (iexpand .eq. 1) then

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)/real(aye**3, RKIND)
            enddo
            if (ispecies .gt. 0) then
               do i = is+1, ie+1
                  de(i,j,k)    = de(i,j,k)/real(aye**3, RKIND)
                  HI(i,j,k)    = HI(i,j,k)/real(aye**3, RKIND)
                  HII(i,j,k)   = HII(i,j,k)/real(aye**3, RKIND)
                  HeI(i,j,k)   = HeI(i,j,k)/real(aye**3, RKIND)
                  HeII(i,j,k)  = HeII(i,j,k)/real(aye**3, RKIND)
                  HeIII(i,j,k) = HeIII(i,j,k)/real(aye**3, RKIND)
               enddo
            endif
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)/real(aye**3, RKIND)
                  H2I(i,j,k)  = H2I(i,j,k)/real(aye**3, RKIND)
                  H2II(i,j,k) = H2II(i,j,k)/real(aye**3, RKIND)
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)/real(aye**3, RKIND)
                  DII(i,j,k) = DII(i,j,k)/real(aye**3, RKIND)
                  HDI(i,j,k) = HDI(i,j,k)/real(aye**3, RKIND)
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)/real(aye**3, RKIND)
               enddo
            endif
         enddo
      enddo
 
      endif

!     Loop over slices (in the k-direction)

      do k = ks+1, ke+1
       do j = js+1, je+1

         do i = is+1, ie+1
            itmask(i) = .true.
         end do

!        Compute the cooling rate

         call cool1d_multi_g(     
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool, idust, idim,
     &                is, ie, j, k, ih2co, ipiht, 1, igammah,
     &                aye, temstart, temend, z_solar,
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                gasgra, mutaba,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, tdust, metallicity, rhoH,
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat, 
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3,
     &                priDataSize, priCooling, priHeating,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3,
     &                metDataSize, metCooling, metHeating,
     &                itmask)

!        Compute the cooling time on the slice
!          (the gamma used here is the same as used to calculate the pressure
!           in cool1d_multi_g)

         do i = is+1, ie+1
            energy = max(p2d(i)/(gamma-1._DKIND), tiny)   
            cooltime(i,j,k) = real(abs(energy/edot(i)), RKIND)
         enddo

       enddo
      enddo

!     Convert densities back to comoving from 'proper'

      if (iexpand .eq. 1) then

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*real(aye**3, RKIND)
            enddo
            if (ispecies .gt. 0) then
               do i = is+1, ie+1
                  de(i,j,k)    = de(i,j,k)*real(aye**3, RKIND)
                  HI(i,j,k)    = HI(i,j,k)*real(aye**3, RKIND)
                  HII(i,j,k)   = HII(i,j,k)*real(aye**3, RKIND)
                  HeI(i,j,k)   = HeI(i,j,k)*real(aye**3, RKIND)
                  HeII(i,j,k)  = HeII(i,j,k)*real(aye**3, RKIND)
                  HeIII(i,j,k) = HeIII(i,j,k)*real(aye**3, RKIND)
               enddo
            endif
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*real(aye**3, RKIND)
                  H2I(i,j,k)  = H2I(i,j,k)*real(aye**3, RKIND)
                  H2II(i,j,k) = H2II(i,j,k)*real(aye**3, RKIND)
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*real(aye**3, RKIND)
                  DII(i,j,k) = DII(i,j,k)*real(aye**3, RKIND)
                  HDI(i,j,k) = HDI(i,j,k)*real(aye**3, RKIND)
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*real(aye**3, RKIND)
               enddo
            endif
         enddo
      enddo

      endif

      return
      end

