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
     &                h2ltea, gasgra,
     &                iradshield, avgsighi, avgsighei, avgsigheii,
     &                k24, k26, iradtrans, photogamma,
     &                ih2optical, iciecool, ciecoa, 
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat)

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
     &        imcool, idust, igammah, ih2optical, iciecool,
     &        clnew, iVheat, iMheat, iradtrans, iradshield

      real*8  aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma, fh, z_solar, clEleFra
      R_PREC  d(in,jn,kn),   e(in,jn,kn),
     &        u(in,jn,kn),   v(in,jn,kn),    w(in,jn,kn),
     &        de(in,jn,kn),  HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn),  H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn),  HDI(in,jn,kn),
     &        metal(in,jn,kn), Vheat(in,jn,kn), Mheat(in,jn,kn),
     &        cooltime(in,jn,kn)
      real*8  photogamma(in,jn,kn)
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
     &        gammaha, avgsighi, avgsighei, avgsigheii,
     &        k24, k26
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

!  Slice locals
 
      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &       tgas(in), tgasold(in),
     &       tdust(in), metallicity(in), rhoH(in)
      real*8 edot(in)

!  Cooling/heating slice locals

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), cieco(in),
     &       hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in)
          
!  Iteration mask for multi_cool

      logical itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
      dk = ke - ks + 1
      dj = je - js + 1

!     Convert densities from comoving to 'proper'

      if (iexpand .eq. 1) then

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   comp1, comp2, energy,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef, p2d,
!$omp&   tgas, tgasold,
!$omp&   tdust, metallicity, rhoH,
!$omp&   edot,
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
#ifdef _OPENMP
!$omp end parallel do
#endif
 
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
!$omp&   tdust, metallicity, rhoH,
!$omp&   edot,
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
     &                h2ltea, gasgra,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, tdust, metallicity, rhoH,
     &                iradshield, avgsighi, avgsighei, avgsigheii,
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
     &                itmask)

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

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   comp1, comp2, energy,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef, p2d,
!$omp&   tgas, tgasold,
!$omp&   tdust, metallicity, rhoH,
!$omp&   edot,
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
#ifdef _OPENMP
!$omp end parallel do
#endif

      endif

      return
      end

