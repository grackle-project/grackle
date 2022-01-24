#include "phys_const.def"
#include "dust_evol.def"

!////////////////////////////////////////////////////////////
! calculate characteristic timescales used for subcycling
!////////////////////////////////////////////////////////////
      subroutine dust_tau_for_timestep(in, jn, kn, i, j, k,
     &              Density, Tgas, Tdust, metal_Density,
     &              gas_metal1, gas_metal2, gas_metal3,
     &              gas_metal4, gas_metal5, gas_metal6,
     &              gas_metal7, gas_metal8, gas_metal9,
     &              gas_metal10, SNe_ThisTimeStep,
     &              tau_accr, tau_dest, tau_sput,
     &              z_solar, SolarAbundances,
     &              time_units, dens_units, len_units, aye, dtime)
      implicit NONE
#include "grackle_fortran_types.def"

      integer in, jn, kn, i, j, k
      real*8 tau_accr(in), tau_dest(in), tau_sput(in),
     &       z_solar, SolarAbundances(NUM_METAL_SPECIES_GRACKLE),
     &       time_units, dens_units, len_units, aye, dtime
      R_PREC Density(in,jn,kn), Tgas(in,jn,kn), Tdust(in,jn,kn), metal_Density(in,jn,kn),
     &       gas_metal1(in,jn,kn), gas_metal2(in,jn,kn),
     &       gas_metal3(in,jn,kn), gas_metal4(in,jn,kn),
     &       gas_metal5(in,jn,kn), gas_metal6(in,jn,kn),
     &       gas_metal7(in,jn,kn), gas_metal8(in,jn,kn),
     &       gas_metal9(in,jn,kn), gas_metal10(in,jn,kn),
     &       SNe_ThisTimeStep(in,jn,kn)
    
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
     &     * (Tdust(i,j,k)/Tgas(i,j,k))**0.5
      tau_accr(i) = 1.d+20
      tau_accr(i) = min(tau_accr0 * (z_solar
     &     / metal_Density(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(1)
     &     / gas_metal1(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(2)
     &     / gas_metal2(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(3)
     &     / gas_metal3(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(4)
     &     / gas_metal4(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(5)
     &     / gas_metal5(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(6)
     &     / gas_metal6(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(7)
     &     / gas_metal7(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(8)
     &     / gas_metal8(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(9)
     &     / gas_metal9(i,j,k)), tau_accr(i))
      tau_accr(i) = min(tau_accr0 * (SolarAbundances(10)
     &     / gas_metal10(i,j,k)), tau_accr(i))

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
     &     dust_Density, metal_Density, Tgas, Tdust,
     &     gas_metal1, gas_metal2, gas_metal3,
     &     gas_metal4, gas_metal5, gas_metal6,
     &     gas_metal7, gas_metal8, gas_metal9,
     &     gas_metal10,
     &     dust_metal1, dust_metal2, dust_metal3,
     &     dust_metal4, dust_metal5, dust_metal6,
     &     dust_metal7, dust_metal8, dust_metal9,
     &     dust_metal10,
     &     SNe_ThisTimeStep,
     &     tau_dest, tau_sput,
     &     z_solar, SolarAbundances,
     &     time_units, dens_units, len_units, aye, dtime, metpt, dmetpt)
      implicit NONE
#include "grackle_fortran_types.def"
!   arguments
      integer in, jn, kn, i, j, k
      real*8 tau_dest(in), tau_sput(in),
     &       z_solar, SolarAbundances(NUM_METAL_SPECIES_GRACKLE),
     &       time_units, dens_units, len_units, aye, dtime
      R_PREC Tgas(in,jn,kn), Tdust(in,jn,kn),
     &       gas_metal1(in,jn,kn), gas_metal2(in,jn,kn),
     &       gas_metal3(in,jn,kn), gas_metal4(in,jn,kn),
     &       gas_metal5(in,jn,kn), gas_metal6(in,jn,kn),
     &       gas_metal7(in,jn,kn), gas_metal8(in,jn,kn),
     &       gas_metal9(in,jn,kn), gas_metal10(in,jn,kn),
     &       dust_metal1(in,jn,kn), dust_metal2(in,jn,kn),
     &       dust_metal3(in,jn,kn), dust_metal4(in,jn,kn),
     &       dust_metal5(in,jn,kn), dust_metal6(in,jn,kn),
     &       dust_metal7(in,jn,kn), dust_metal8(in,jn,kn),
     &       dust_metal9(in,jn,kn), dust_metal10(in,jn,kn),
     &       metal_Density(in,jn,kn), dust_Density(in,jn,kn),
     &       SNe_ThisTimeStep(in,jn,kn),
     &       metpt(NUM_METAL_SPECIES_GRACKLE+1), dmetpt(NUM_METAL_SPECIES_GRACKLE+1)

!   local variables
      integer n
      character(80) :: userinp
      real*8 tau_ref, Ms100, tau_accr0, tau_accr
      R_PREC dM(NUM_METAL_SPECIES_GRACKLE+1), dMs,
     &       Mmet(NUM_METAL_SPECIES_GRACKLE+1),
     &       dens_proper

      dens_proper = dens_units * aye**3

!     copy metallicites into a single array
      metpt(1) = metal_Density(i,j,k)
      metpt(2) = gas_metal1(i,j,k)
      metpt(3) = gas_metal2(i,j,k)
      metpt(4) = gas_metal3(i,j,k)
      metpt(5) = gas_metal4(i,j,k)
      metpt(6) = gas_metal5(i,j,k)
      metpt(7) = gas_metal6(i,j,k)
      metpt(8) = gas_metal7(i,j,k)
      metpt(9) = gas_metal8(i,j,k)
      metpt(10) = gas_metal9(i,j,k)
      metpt(11) = gas_metal10(i,j,k)
      dmetpt(1) = dust_Density(i,j,k)
      dmetpt(2) = dust_metal1(i,j,k)
      dmetpt(3) = dust_metal2(i,j,k)
      dmetpt(4) = dust_metal3(i,j,k)
      dmetpt(5) = dust_metal4(i,j,k)
      dmetpt(6) = dust_metal5(i,j,k)
      dmetpt(7) = dust_metal6(i,j,k)
      dmetpt(8) = dust_metal7(i,j,k)
      dmetpt(9) = dust_metal8(i,j,k)
      dmetpt(10) = dust_metal9(i,j,k)
      dmetpt(11) = dust_metal10(i,j,k)

      Ms100 = 6800.0*SNE_COEFF
     &        *(100.0/SNE_SHOCKSPEED)*(100.0/SNE_SHOCKSPEED)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
      tau_ref = DUST_GROWTH_TAUREF * 1e9 * SEC_PER_YEAR / time_units

      do n=0+1,NUM_METAL_SPECIES_GRACKLE+1
          dM(n) = 0.
      enddo

!   metals accreted from gas to dust
      tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF/dens_proper)
     &     * (Tdust(i,j,k)/Tgas(i,j,k))**0.5
      do n=1+1,NUM_METAL_SPECIES_GRACKLE+1
          Mmet(n) = metpt(n) + dmetpt(n)
          tau_accr = tau_accr0 * SolarAbundances(n-1) / metpt(n)

          if (dmetpt(n) .ne. dmetpt(n)) then
            write(*,*) "dmetpt at index", n, 
     &                 "calculated as NaN in line 174 of dust model", dmetpt(n)
          endif 

          if (Mmet(n) .ne. Mmet(n)) then
            write(*,*) "Mmet at index", n, 
     &                 "calculated as NaN in line 174 of dust model", Mmet(n)
          endif 

          if (tau_accr .ne. tau_accr) then
            write(*,*) "tau_accr",
     &                 "calculated as NaN in line 174 of dust model", tau_accr
          endif 

          if (dtime .ne. dtime) then
            write(*,*) "dtime",
     &                 "calculated as NaN in line 174 of dust model", dtime
          endif 

          dM(n) = dM(n) + min((1 - dmetpt(n) / Mmet(n))
     &            * (dmetpt(n)/tau_accr) * dtime,
     &            Mmet(n) - dmetpt(n)) ! growth capped by gas in a single cell

          if (dM(n) .ne. dM(n)) then
              write(*,*) "dM at index", n, 
     &                   "calculated as NaN in line 174 of dust model", dM(n)
              read(*,*) userinp
          endif 

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
            if (dMs .gt. dust_Density(i,j,k)) then
                  write (*,*) "WARNING: dMs > Mdust SNe shock destruction:", 
     &                SNe_ThisTimeStep(i,j,k), tau_dest(i)
            endif
      else
          dMs = dMs + dust_Density(i,j,k) / tau_sput(i) *3.0* dtime
          dMs = min(dMs,dust_Density(i,j,k))
      endif
      do n=0+1,NUM_METAL_SPECIES_GRACKLE+1
          dM(n) = dM(n) - dmetpt(n) * dMs

          if (dM(n) .ne. dM(n)) then
            write(*,*) "dM at index", n, 
     &                 "calculated as NaN in line 202 of dust model", dM(n)
            read(*,*) userinp
        endif 
          
      enddo

!   recalculate metallicity
      do n=0+1,NUM_METAL_SPECIES_GRACKLE+1
          if (dust_Density(i,j,k) .gt. 0.0) then
              dmetpt(n) = dM(n) + dmetpt(n) ! SNe shocks on neighbors are ignored
          else
              dmetpt(n) = 0.0
          endif
          metpt(n) = metpt(n) - dM(n) 

      enddo

      

!     copy metallicites back into fields arrays
      metal_Density(i,j,k) = metpt(1)
      gas_metal1(i,j,k) = metpt(2)
      gas_metal2(i,j,k) = metpt(3)
      gas_metal3(i,j,k) = metpt(4)
      gas_metal4(i,j,k) = metpt(5)
      gas_metal5(i,j,k) = metpt(6)
      gas_metal6(i,j,k) = metpt(7)
      gas_metal7(i,j,k) = metpt(8)
      gas_metal8(i,j,k) = metpt(9)
      gas_metal9(i,j,k) = metpt(10)
      gas_metal10(i,j,k) = metpt(11)
      dust_Density(i,j,k) = dmetpt(1)
      dust_metal1(i,j,k) = dmetpt(2)
      dust_metal2(i,j,k) = dmetpt(3)
      dust_metal3(i,j,k) = dmetpt(4)
      dust_metal4(i,j,k) = dmetpt(5)
      dust_metal5(i,j,k) = dmetpt(6)
      dust_metal6(i,j,k) = dmetpt(7)
      dust_metal7(i,j,k) = dmetpt(8)
      dust_metal8(i,j,k) = dmetpt(9)
      dust_metal9(i,j,k) = dmetpt(10)
      dust_metal10(i,j,k) = dmetpt(11)

      return
      end ! end of dust growth and destruction module
