! Density must be in [CGS], if not, covert the unit before calling routines here
! Mass must be in code units
! Timescales are in code units

#include "dust_evol.h"
#define NUM_METAL_SPECIES 11

!////////////////////////////////////////////////////////////
! calculate characteristic timescales used for subcycling
!////////////////////////////////////////////////////////////
      subroutine dust_tau_for_timestep(in, jn, kn, i, j, k,
     &              Density, Tgas,
     &              Metallicity, SNe_ThisTimeStep,
     &              tau_accr, tau_dest, tau_sput,
     &              SolarAbundances,
     &              time_units, dens_units, len_units, dtime)
      implicit NONE

      integer in, jn, kn, i, j, k
      real*8 tau_accr(in), tau_dest(in), tau_sput(in)
     &       SolarAbundances(NUM_METAL_SPECIES),
     &       time_units, dens_units, len_units, dtime
      R_PREC Density(in,jn,kn), Tgas(in,jn,kn)
     &       dust_Density(in,jn,kn), 
     &       Metallicity(NUM_METAL_SPECIES,in,jn,kn), 
     &       SNe_ThisTimeStep(in,jn,kn)
    
!   local variables
      integer n
      real*8 tau_ref, Ms100, tau_accr0

      Ms100 = 6800.0*SNE_COEFF
     &        *(100.0/SNE_SHOCKSPEED)*(100.0/SNE_SHOCKSPEED)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
! check: 
      tau_ref = DUST_GROWTH_TAUREF * 1e9 * SEC_PER_YEAR / time_units

!   growth
      tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF / dens_units) 
     &            / Density(i,j,k) * (DUST_GROWTH_TREF/Tgas(i,j,k))**0.5
      tau_accr(i) = 1.d+20
      do n=1+1,NUM_METAL_SPECIES ! CHECK FORTRAN ARRAY INDEX, if i remember correctly it starts at 1 instead of 0
          tau_accr(i) = min(tau_accr0 * (SolarAbundances(n)
     &                  / Metallicity(n,i,j,k)),tau_accr(i,j,k))

!   destruction by SN shocks
      if (SNe_ThisTimeStep.le.0.0) then ! CHECK float type
          tau_dest(i) = 1.d+20
      else
          tau_dest(i) = Density/(Ms100*SNe_ThisTimeStep(i,j,k)
     &                      *DUST_DESTRUCTION_EFF) * dtime
      endif

!   destruction by thermal sputtering
      tau_sput(i,j,k) = 1.7e8 * SEC_PER_YEAR / time_units
     &                  * (DUST_GRAINSIZE/0.1) * (1e-27/Density(i,j,k))
     &                  * ((2e6/Tgas(i,j,k))**2.5+1.0) ! Tsai & Mathews (1995)

      return
      end ! end of dust_tau_for_timestep


!////////////////////////////////////////////////////////////
!   main routine for dust growth and destruction
!////////////////////////////////////////////////////////////
      subroutine dust_growth_and_destruction(in, jn, kn, i, j, k,
     &              Density, dust_Density, fgr, Tgas,
     &              Metallicity, dust_Metallicity, SNe_ThisTimeStep,
     &              tau_dest, tau_sput, 
     &              SolarAbundances, 
     &              time_units, dens_units, len_units dtime)
      implicit NONE
!   arguments
      integer in, jn, kn, i, j, k
      real*8 tau_dest(in,jn,kn), tau_sput(in,jn,kn)
     &       SolarAbundances(NUM_METAL_SPECIES), 
     &       time_units, dens_units, len_units, dtime
      R_PREC Density(in,jn,kn), Tgas(in,jn,kn),
     &       dust_Density(in,jn,kn), fgr(in,jn,kn)
     &       Metallicity(NUM_METAL_SPECIES,in,jn,kn), 
     &       SNe_ThisTimeStep(in,jn,kn),
     &       dust_Metallicity(NUM_METAL_SPECIES,in,jn,kn)

!   local variables
      integer n
      real*8 tau_ref, Ms100, tau_accr0, tau_accr
      R_PREC dM(NUM_METAL_SPECIES) = 0, dMs,
     &       Mmet(NUM_METAL_SPECIES) = 0 ,Mdust[NUM_METAL_SPECIES]

      Ms100 = 6800.0*SNE_COEFF
     &        *(100.0/SNE_SHOCKSPEED)*(100.0/SNE_SHOCKSPEED)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
      tau_ref = DUST_GROWTH_TAUREF * 1e9 * SEC_PER_YEAR / time_units

      do n=1+1,NUM_METAL_SPECIES ! CHECK FORTRAN ARRAY INDEX, if i remember correctly it starts at 1 instead of 0 
          Mdust(n) = dust_Metallicity(n,i,j,k) * dust_Density(i,j,k)
          dM(n) = 0.
      enddo

!   metals accreted from gas to dust
      tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF/dens_units
     &            /Density(i,j,k)) * (DUST_GROWTH_TREF/Tgas(i,j,k))**0.5
      do n=1+1,NUM_METAL_SPECIES
          Mmet(n) = Metallicity(n,i,j,k) * Density(i,j,k) + Mdust(n)
          tau_accr = tau_accr0 * (SolarAbundances(n) 
     &               / Metallicity(n,i,j,k))
          dM(n) = dM(n) + min((1 - Mdust(n) / Mmet(n)) 
     &            * (Mdust(n)/tau_accr) * dtime_sub,
     &            Mmet(n) - Mdust(n)) ! growth capped by gas in a single cell
          dM(1) = dM(1) + dM(n)
      enddo

!   SNe shock and thermal sputtering
      if (SNe_ThisTimeStep(i,j,k) .le. 0.0)
          dMs = 0.0
      else
          dMs = min(dust_Mass(i,j,k)/tau_dest(i,j,k)*dtime,
     &          dust_Mass(i,j,k))
      if (dMs .ge. dust_Mass(i,j,k)) then
          write (*,*) "WARNING: dMs >= Mdust SNe shock destruction:", 
     &                SNe_ThisTimeStep(i,j,k), tau_dest(i,j,k)
      else
          dMs = dMs + dust_Mass(i,j,k) / tau_sput(i,j,k) * 3.0 * dtime
          dMs = min(dMs,dust_Mass(i,j,k))
      endif
      do n=0+1,NUM_METAL_SPECIES
          dM(n) = dM(n) - dust_Metallicity(n,i,j,k) * dMs
      enddo

!   recalculate metallicity
      dust_Density(i,j,k) = dust_Density(i,j,k) + dM(1)
      do n=0+1,NUM_METAL_SPECIES
          if (dust_Density(i,j,k) .lt. 0.0) then 
              dust_Metallicity(n,i,j,k) =  (dM(n) + Mdust(n))
     &        / dust_Density(i,j,k)  ! SNe shocks on neighbors are ignored
          else
              dust_Metallicity(n,i,j,k) = 0.0
          endif
          Metallicity(n,i,j,k) = (Metallicity(n,i,j,k) * Density(i,j,k) 
     &                            - dM(n)) / (Density(i,j,k) - dM(1))
      enddo

      return
      end ! end of dust growth and destruction module
