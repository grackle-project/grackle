#include "phys_const.def"
#include "dust_evol.def"

!////////////////////////////////////////////////////////////
! calculate characteristic timescales used for subcycling
!////////////////////////////////////////////////////////////
      subroutine dust_tau_for_timestep(in, jn, kn, i, j, k,
     &              Density, Tgas,
     &              Metallicity1, Metallicity2, Metallicity3,
     &              Metallicity4, Metallicity5, Metallicity6,
     &              Metallicity7, Metallicity8, Metallicity9,
     &              Metallicity10, Metallicity11, SNe_ThisTimeStep,
     &              tau_accr, tau_dest, tau_sput,
     &              SolarAbundances,
     &              time_units, dens_units, len_units, aye, dtime)
      implicit NONE
#include "grackle_fortran_types.def"

      integer in, jn, kn, i, j, k
      real*8 tau_accr(in), tau_dest(in), tau_sput(in),
     &       SolarAbundances(NUM_METAL_SPECIES),
     &       time_units, dens_units, len_units, aye, dtime
      R_PREC Density(in,jn,kn), Tgas(in,jn,kn),
     &       Metallicity1(in,jn,kn), Metallicity2(in,jn,kn),
     &       Metallicity3(in,jn,kn), Metallicity4(in,jn,kn),
     &       Metallicity5(in,jn,kn), Metallicity6(in,jn,kn),
     &       Metallicity7(in,jn,kn), Metallicity8(in,jn,kn),
     &       Metallicity9(in,jn,kn), Metallicity10(in,jn,kn),
     &       Metallicity11(in,jn,kn), SNe_ThisTimeStep(in,jn,kn)
    
!   local variables
      real*8 tau_ref, Ms100, tau_accr0, dens_proper

      dens_proper = dens_units * aye**3

      Ms100 = 6800.0*SNE_COEFF
     &        *(100.0/SNE_SHOCKSPEED)*(100.0/SNE_SHOCKSPEED)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
! check: 
      tau_ref = DUST_GROWTH_TAUREF * 1e9 * SEC_PER_YEAR / time_units

!   growth
      tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF / dens_proper)
     &            / Density(i,j,k) * (DUST_GROWTH_TREF/Tgas(i,j,k))**0.5
      tau_accr(i) = 1.d+20
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(1)
     &     / Metallicity1(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(2)
     &     / Metallicity2(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(3)
     &     / Metallicity3(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(4)
     &     / Metallicity4(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(5)
     &     / Metallicity5(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(6)
     &     / Metallicity6(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(7)
     &     / Metallicity7(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(8)
     &     / Metallicity8(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(9)
     &     / Metallicity9(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(10)
     &     / Metallicity10(i,j,k)),tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(11)
     &     / Metallicity11(i,j,k)),tau_accr(i))

!   destruction by SN shocks
      if (SNe_ThisTimeStep(i,j,k).le.0.0) then
          tau_dest(i) = 1.d+20
      else
          tau_dest(i) = Density(i,j,k)/(Ms100*SNe_ThisTimeStep(i,j,k)
     &                      *DUST_DESTRUCTION_EFF) * dtime
      endif

!   destruction by thermal sputtering
      tau_sput(i) = 1.7e8 * SEC_PER_YEAR / time_units
     &     * (DUST_GRAINSIZE/0.1) * (1e-27/(dens_proper*Density(i,j,k)))
     &     * ((2e6/Tgas(i,j,k))**2.5+1.0) ! Tsai & Mathews (1995)

      return
      end ! end of dust_tau_for_timestep


!////////////////////////////////////////////////////////////
!   main routine for dust growth and destruction
!////////////////////////////////////////////////////////////
      subroutine dust_growth_and_destruction(in, jn, kn, i, j, k,
     &     Density, dust_Density, metal_Density, Tgas,
     &     Metallicity1, Metallicity2, Metallicity3,
     &     Metallicity4, Metallicity5, Metallicity6,
     &     Metallicity7, Metallicity8, Metallicity9,
     &     Metallicity10, Metallicity11,
     &     dust_Metallicity1, dust_Metallicity2, dust_Metallicity3,
     &     dust_Metallicity4, dust_Metallicity5, dust_Metallicity6,
     &     dust_Metallicity7, dust_Metallicity8, dust_Metallicity9,
     &     dust_Metallicity10, dust_Metallicity11,
     &     SNe_ThisTimeStep,
     &     tau_dest, tau_sput,
     &     SolarAbundances,
     &     time_units, dens_units, len_units, aye, dtime, metpt, dmetpt)
      implicit NONE
#include "grackle_fortran_types.def"
!   arguments
      integer in, jn, kn, i, j, k
      real*8 tau_dest(in), tau_sput(in),
     &       SolarAbundances(NUM_METAL_SPECIES), 
     &       time_units, dens_units, len_units, aye, dtime
      R_PREC Density(in,jn,kn), Tgas(in,jn,kn),
     &       Metallicity1(in,jn,kn), Metallicity2(in,jn,kn),
     &       Metallicity3(in,jn,kn), Metallicity4(in,jn,kn),
     &       Metallicity5(in,jn,kn), Metallicity6(in,jn,kn),
     &       Metallicity7(in,jn,kn), Metallicity8(in,jn,kn),
     &       Metallicity9(in,jn,kn), Metallicity10(in,jn,kn),
     &       Metallicity11(in,jn,kn),
     &       dust_Metallicity1(in,jn,kn), dust_Metallicity2(in,jn,kn),
     &       dust_Metallicity3(in,jn,kn), dust_Metallicity4(in,jn,kn),
     &       dust_Metallicity5(in,jn,kn), dust_Metallicity6(in,jn,kn),
     &       dust_Metallicity7(in,jn,kn), dust_Metallicity8(in,jn,kn),
     &       dust_Metallicity9(in,jn,kn), dust_Metallicity10(in,jn,kn),
     &       dust_Metallicity11(in,jn,kn),
     &       metal_Density(in,jn,kn), dust_Density(in,jn,kn),
     &       SNe_ThisTimeStep(in,jn,kn),
     &       metpt(NUM_METAL_SPECIES), dmetpt(NUM_METAL_SPECIES)

!   local variables
      integer n
      real*8 tau_ref, Ms100, tau_accr0, tau_accr
      R_PREC dM(NUM_METAL_SPECIES), dMs,
     &       Mmet(NUM_METAL_SPECIES),Mdust(NUM_METAL_SPECIES),
     &       dens_proper

      dens_proper = dens_units * aye**3

!     copy metallicites into a single array
      metpt(1) = Metallicity1(i,j,k)
      metpt(2) = Metallicity2(i,j,k)
      metpt(3) = Metallicity3(i,j,k)
      metpt(4) = Metallicity4(i,j,k)
      metpt(5) = Metallicity5(i,j,k)
      metpt(6) = Metallicity6(i,j,k)
      metpt(7) = Metallicity7(i,j,k)
      metpt(8) = Metallicity8(i,j,k)
      metpt(9) = Metallicity9(i,j,k)
      metpt(10) = Metallicity10(i,j,k)
      metpt(11) = Metallicity11(i,j,k)
      dmetpt(1) = dust_Metallicity1(i,j,k)
      dmetpt(2) = dust_Metallicity2(i,j,k)
      dmetpt(3) = dust_Metallicity3(i,j,k)
      dmetpt(4) = dust_Metallicity4(i,j,k)
      dmetpt(5) = dust_Metallicity5(i,j,k)
      dmetpt(6) = dust_Metallicity6(i,j,k)
      dmetpt(7) = dust_Metallicity7(i,j,k)
      dmetpt(8) = dust_Metallicity8(i,j,k)
      dmetpt(9) = dust_Metallicity9(i,j,k)
      dmetpt(10) = dust_Metallicity10(i,j,k)
      dmetpt(11) = dust_Metallicity11(i,j,k)

      Ms100 = 6800.0*SNE_COEFF
     &        *(100.0/SNE_SHOCKSPEED)*(100.0/SNE_SHOCKSPEED)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
      tau_ref = DUST_GROWTH_TAUREF * 1e9 * SEC_PER_YEAR / time_units

      do n=0+1,NUM_METAL_SPECIES 
          Mdust(n) = dmetpt(n) * Density(i,j,k)
          dM(n) = 0.
      enddo

!   metals accreted from gas to dust
      tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF/dens_proper
     &            /Density(i,j,k)) * (DUST_GROWTH_TREF/Tgas(i,j,k))**0.5
      do n=1+1,NUM_METAL_SPECIES
          Mmet(n) = metpt(n) * Density(i,j,k) + Mdust(n)
          tau_accr = tau_accr0 * (SolarAbundances(n) 
     &               / metpt(n))
          dM(n) = dM(n) + min((1 - Mdust(n) / Mmet(n)) 
     &            * (Mdust(n)/tau_accr) * dtime,
     &            Mmet(n) - Mdust(n)) ! growth capped by gas in a single cell
          dM(1) = dM(1) + dM(n)
      enddo

!   SNe shock and thermal sputtering
      if (SNe_ThisTimeStep(i,j,k) .le. 0.0) then
          dMs = 0.0
      else
          dMs = min(dust_Density(i,j,k)/tau_dest(i)*dtime,
     &          dust_Density(i,j,k))
      endif
      if (dMs .ge. dust_Density(i,j,k)) then
          write (*,*) "WARNING: dMs >= Mdust SNe shock destruction:", 
     &                SNe_ThisTimeStep(i,j,k), tau_dest(i)
      else
          dMs = dMs + dust_Density(i,j,k) / tau_sput(i) *3.0* dtime
          dMs = min(dMs,dust_Density(i,j,k))
      endif
      do n=0+1,NUM_METAL_SPECIES
          dM(n) = dM(n) - dmetpt(n) * dMs
      enddo

!   recalculate metallicity
      dust_Density(i,j,k) = dust_Density(i,j,k) + dM(1)
      metal_Density(i,j,k) = metal_Density(i,j,k) - dM(1)
      do n=0+1,NUM_METAL_SPECIES
          if (dust_Density(i,j,k) .gt. 0.0) then
              dmetpt(n) =  (dM(n) + Mdust(n))
     &        / Density(i,j,k)  ! SNe shocks on neighbors are ignored
          else
              dmetpt(n) = 0.0
          endif
          metpt(n) = (metpt(n) * Density(i,j,k)
     &                            - dM(n)) / (Density(i,j,k) - dM(1))
      enddo

!     copy metallicites back into fields arrays
      Metallicity1(i,j,k) = metpt(1)
      Metallicity2(i,j,k) = metpt(2)
      Metallicity3(i,j,k) = metpt(3)
      Metallicity4(i,j,k) = metpt(4)
      Metallicity5(i,j,k) = metpt(5)
      Metallicity6(i,j,k) = metpt(6)
      Metallicity7(i,j,k) = metpt(7)
      Metallicity8(i,j,k) = metpt(8)
      Metallicity9(i,j,k) = metpt(9)
      Metallicity10(i,j,k) = metpt(10)
      Metallicity11(i,j,k) = metpt(11)
      dust_Metallicity1(i,j,k) = dmetpt(1)
      dust_Metallicity2(i,j,k) = dmetpt(2)
      dust_Metallicity3(i,j,k) = dmetpt(3)
      dust_Metallicity4(i,j,k) = dmetpt(4)
      dust_Metallicity5(i,j,k) = dmetpt(5)
      dust_Metallicity6(i,j,k) = dmetpt(6)
      dust_Metallicity7(i,j,k) = dmetpt(7)
      dust_Metallicity8(i,j,k) = dmetpt(8)
      dust_Metallicity9(i,j,k) = dmetpt(9)
      dust_Metallicity10(i,j,k) = dmetpt(10)
      dust_Metallicity11(i,j,k) = dmetpt(11)

      return
      end ! end of dust growth and destruction module
