#include "phys_const.def"

!=======================================================================
!/////////////////////  SUBROUTINE COOL1D_MULTI_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_multi_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool, idust, idim, 
     &                is, ie, j, k, ih2co, ipiht, iter, igammah,
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
     &                itmask)

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
     &        iexpand, ih2co, ipiht, ispecies, imcool, idust, 
     &        imetal, igammah, ih2optical, iciecool, clnew,
     &        iVheat, iMheat, iradtrans, iradshield

      real*8  aye, temstart, temend, z_solar,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma, fh, clEleFra
      R_PREC  d(in,jn,kn),    e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn), Vheat(in,jn,kn), Mheat(in,jn,kn)
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
     &        gammaha, gammaha_eff, avgsighi, avgsighei,
     &        avgsigheii, k24, k26

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
     &       nH2, nother, myde, nSSh, nratio, nssh_he, nratio_he,
     &       fSShHI, fSShHeI

!  Slice locals
 
      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &       tgas(in), tgasold(in), tdust(in), rhoH(in), nh(in),
     &       metallicity(in), edot(in)

      real*8  mmw(in)

!  Cooling/heating slice locals

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), cieco(in),
     &       hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in),
     &       gaHI(in), gaH2(in), gaHe(in), gaHp(in), gael(in),
     &       h2lte(in), galdl(in), gasgr(in), gasgr_tdust(in)

!  Iteration mask

      logical itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

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

         do i = is+1, ie+1
            if ( itmask(i) ) then
               rhoH(i) = fh * d(i,j,k)
            endif
         enddo

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
               tgas(i) = 
     &              (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4._DKIND +
     &              HI(i,j,k) + HII(i,j,k) + de(i,j,k)
               rhoH(i) = HI(i,j,k) + HII(i,j,k)
            end if
         enddo

!     (include molecular hydrogen, but ignore deuterium)

         if (ispecies .gt. 1) then
            do i = is+1, ie+1
               if ( itmask(i) ) then
                  tgas(i) = tgas(i) +
     &                 HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2._DKIND
                  rhoH(i) = rhoH(i) + H2I(i,j,k) + H2II(i,j,k)
               end if
            enddo
         endif

!     Include metal species

         if (imetal .eq. 1) then
            do i = is+1, ie+1
               tgas(i) = tgas(i) + metal(i,j,k)/mu_metal
            enddo
         endif

         do i = is+1, ie+1
            if ( itmask(i) ) then
               tgas(i) = max(p2d(i)*utem/tgas(i), temstart)
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
                  tgas(i) = tgas(i) * (gamma2 - 1._DKIND)/
     &                 (gamma - 1._DKIND)
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
            nh(i) = rhoH(i) * dom
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

      do i = is+1, ie+1
         if ( itmask(i) ) then

!        Compute log temperature and truncate if above/below table max/min

         logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

         endif
      enddo

!     --- 6 species cooling ---

      if (ispecies .gt. 0) then

      do i = is+1, ie+1
         if ( itmask(i) ) then

!        Compute index into the table and precompute parts of linear interp

         indixe(i) = min(nratec-1,
     &           max(1,int((logtem(i)-logtem0)/dlogtem, DIKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

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
         if ( tgas(i) > 2.d4 ) then
            gammaha_eff = 0 
         else
            gammaha_eff = gammaha
         endif
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

!                    Compton cooling or heating

     &             - comp1*(tgas(i)-comp2)*de(i,j,k)*dom_inv

!                    X-ray compton heating

     &             - comp_xraya*(tgas(i)-comp_temp)*de(i,j,k)*dom_inv

!                    Bremsstrahlung

     &             - brem(i)*(HII(i,j,k)+HeII(i,j,k)/4._DKIND +
     &           HeIII(i,j,k)) * de(i,j,k)

!                    Photoelectric heating by UV-irradiated dust

     &             + real(igammah, DKIND) * gammaha_eff *
     &           (HI(i,j,k) + HII(i,j,k)) * dom_inv

     &             )

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
                edot(i) = ciefudge*(edot(i) - 
     &                  H2I(i,j,k)*(d(i,j,k)*cieco(i)))
              endif
            endif
            enddo
         endif

      endif

!     --- Cooling from HD ---

      if (ispecies .gt. 2) then
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
     .              (hdlte(i)/(1._DKIND + hdlte1/hdlow1)) /
     &              (3._DKIND*dom)
            end if
         enddo
      endif

!     --- Gas to grain heat transfer ---

      if (idust .gt. 0) then

!     Look up gas/grain heat transfer rates

         do i = is+1, ie+1
            if ( itmask(i) ) then
               gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &              *(gasgra(indixe(i)+1) -gasgra(indixe(i)))
               gasgr_tdust(i) = gasgr(i) * coolunit / mh
            endif
         enddo

!     Compute dust temperature

         call calc_tdust_1d_g(tdust, tgas, nh, gasgr_tdust, 
     &           itmask, comp2, in, is, ie, j, k)

!     Calculate cooling rate

         do i = is+1, ie+1
            if ( itmask(i) ) then
               edot(i) = edot(i) - 
     &              gasgr(i) * (tgas(i) - tdust(i)) * 
     &              metallicity(i) * rhoH(i) * rhoH(i)
            endif
         enddo

      endif

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
     &            + piHeI * HeI(i,j,k)*0.25_DKIND*
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

!     Compton cooling or heating/photoelectric heating

         do i = is+1, ie+1
            if ( itmask(i) ) then
            if ( tgas(i) > 2.d4 ) then
               gammaha_eff = 0 
            else
               gammaha_eff = gammaha
            endif

!     Calculate electron density from mean molecular weight

               myde = 1 - mmw(i) * (3.0_DKIND * fh + 1.0_DKIND) /
     &              4.0_DKIND
               if (imetal .eq. 1) then
                  myde = myde - mmw(i) * metal(i,j,k) / 
     &                 (d(i,j,k) * mu_metal)
               endif
               myde = d(i,j,k) * myde / mmw(i)

               edot(i) = edot(i) 
     &              - comp1*(tgas(i)-comp2)*myde*dom_inv
     &              + real(igammah, DKIND)*gammaha_eff*rhoH(i)*dom_inv

            end if
         enddo

      endif
 
!     Photoheating from radiative transfer

      if (iradtrans .eq. 1) then
          do i = is + 1, ie + 1
            if (itmask(i)) then
              edot(i) = edot(i) + real(ipiht, DKIND) * photogamma(i,j,k)
     &                          / coolunit * HI(i,j,k) / dom
c              if (photogamma(i,j,k) .gt. 0) then
c                print*, i,j,k,edot(i), photogamma(i,j,k), dom,
c      $                 aye, utim, coolunit, d(i,j,k), hi(i,j,k), tgas(i)
c              endif
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

      endif

!     Add user-provided volumetric and/or specific heating terms

      if (iVheat .eq. 1) then

      do i = is+1, ie+1
         if ( itmask(i) ) then
            edot(i) = edot(i) + (
     &           Vheat(i,j,k) * d(i,j,k) / coolunit
     &           ) / dom
         end if
      enddo

      endif

      if (iMheat .eq. 1) then

      do i = is+1, ie+1
         if ( itmask(i) ) then
            edot(i) = edot(i) + (
     &           Mheat(i,j,k) * d(i,j,k)**2 / coolunit)
         end if
      enddo

      endif

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) ) then
         tgasold(i) = tgas(i)
         end if
      enddo

      return
      end
