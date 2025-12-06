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
     &              dust_destruction_eff,
     &              sne_coeff, sne_shockspeed,
     &              dust_grainsize, dust_growth_densref,
     &              dust_growth_tauref,
     &              tau_accr, tau_dest, tau_sput,
     &              z_solar, SolarAbundances,
     &              time_units, dens_units, len_units, aye, dtime)
      implicit NONE
#include "grackle_fortran_types.def"

      character(80) :: userinp
      integer in, jn, kn, i, j, k
      real*8 tau_accr(in), tau_dest(in), tau_sput(in),
     &       z_solar, SolarAbundances(NUM_METAL_SPECIES_GRACKLE),
     &       dust_destruction_eff, sne_coeff, sne_shockspeed,
     &       dust_grainsize, dust_growth_densref,
     &       dust_growth_tauref,
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

      Ms100 = 6800.0*sne_coeff
     &        *(100.0/sne_shockspeed)*(100.0/sne_shockspeed)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
! check: 
      tau_ref = dust_growth_tauref * 1e9 * SEC_PER_YEAR_GRACKLE / time_units

!   growth
      tau_accr0 = tau_ref * (dust_growth_densref / dens_proper)
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

      if (tau_accr(i) .le. 0) then
            write(*,*) "-- tau_accr neagative! --", tau_accr(i)
            write(*,*) "tau_accr0=", tau_accr0
            write(*,*) "metal_density=", metal_density(i,j,k)
            write(*,*) "gas_metal1=", gas_metal1(i,j,k)
            write(*,*) "gas_metal2=", gas_metal2(i,j,k)
            write(*,*) "gas_metal3=", gas_metal3(i,j,k)
            write(*,*) "gas_metal4=", gas_metal4(i,j,k)
            write(*,*) "gas_metal5=", gas_metal5(i,j,k)
            write(*,*) "gas_metal6=", gas_metal6(i,j,k)
            write(*,*) "gas_metal7=", gas_metal7(i,j,k)
            write(*,*) "gas_metal8=", gas_metal8(i,j,k)
            write(*,*) "gas_metal9=", gas_metal9(i,j,k)
            write(*,*) "gas_metal10=", gas_metal10(i,j,k)
            write(*,*) "SolarAbundances1=", SolarAbundances(1)
            write(*,*) "SolarAbundances2=", SolarAbundances(2)
            write(*,*) "SolarAbundances3=", SolarAbundances(3)
            write(*,*) "SolarAbundances4=", SolarAbundances(4)
            write(*,*) "SolarAbundances5=", SolarAbundances(5)
            write(*,*) "SolarAbundances6=", SolarAbundances(6)
            write(*,*) "SolarAbundances7=", SolarAbundances(7)
            write(*,*) "SolarAbundances8=", SolarAbundances(8)
            write(*,*) "SolarAbundances9=", SolarAbundances(9)
            write(*,*) "SolarAbundances10=", SolarAbundances(10)
            call exit(74)
      endif

!   destruction by SN shocks
      if (SNe_ThisTimeStep(i,j,k).le.0.0) then
          tau_dest(i) = 1.d+20
      else
          tau_dest(i) = Density(i,j,k)/(Ms100*SNe_ThisTimeStep(i,j,k)
     &                      *dust_destruction_eff) * dtime
      endif

!   destruction by thermal sputtering
      tau_sput(i) = 1.7e8 * SEC_PER_YEAR_GRACKLE / time_units
     &     * (dust_grainsize/0.1) * (1e-27/(dens_proper*Density(i,j,k)))
     &     * ((2e6/Tgas(i,j,k))**2.5+1.0) ! Tsai & Mathews (1995)

      return
      end ! end of dust_tau_for_timestep


!////////////////////////////////////////////////////////////
!   main routine for dust growth and destruction
!////////////////////////////////////////////////////////////
      subroutine dust_growth_and_destruction(in, jn, kn, i, j, k,
     &     gas_Density, dust_Density, metal_Density, Tgas, Tdust,
     &     gas_metal1, gas_metal2, gas_metal3,
     &     gas_metal4, gas_metal5, gas_metal6,
     &     gas_metal7, gas_metal8, gas_metal9,
     &     gas_metal10,
     &     dust_metal1, dust_metal2, dust_metal3,
     &     dust_metal4, dust_metal5, dust_metal6,
     &     dust_metal7, dust_metal8, dust_metal9,
     &     dust_metal10,
     &     SNe_ThisTimeStep, dust_destruction_eff,
     &     sne_coeff, sne_shockspeed,
     &     dust_grainsize, dust_growth_densref,
     &     dust_growth_tauref,
     &     tau_dest, tau_sput,
     &     z_solar, SolarAbundances,
     &     time_units, dens_units, len_units, aye, dtime, metpt, dmetpt)
      implicit NONE
#include "grackle_fortran_types.def"
!   arguments
      integer in, jn, kn, i, j, k
      real*8 tau_dest(in), tau_sput(in),
     &       z_solar, SolarAbundances(NUM_METAL_SPECIES_GRACKLE),
     &       dust_destruction_eff, sne_coeff, sne_shockspeed,
     &       dust_grainsize, dust_growth_densref,
     &       dust_growth_tauref,
     &       time_units, dens_units, len_units, aye, dtime, dM_conserv
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
     &       gas_Density(in,jn,kn), SNe_ThisTimeStep(in,jn,kn),
     &       metpt(NUM_METAL_SPECIES_GRACKLE+1), dmetpt(NUM_METAL_SPECIES_GRACKLE+1)

!   local variables
      integer n
      character(80) :: userinp
      real*8 tau_ref, Ms100, tau_accr0, tau_accr
      R_PREC dM(NUM_METAL_SPECIES_GRACKLE+1), dMs,
     &       Mmet(NUM_METAL_SPECIES_GRACKLE+1),
     &       dens_proper

!   debug variables
      real*8 dM_tau_accr
      R_PREC total_density_init, total_density_final 

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

      total_density_init = metpt(1) + dmetpt(1)

      Ms100 = 6800.0*sne_coeff
     &        *(100.0/sne_shockspeed)*(100.0/sne_shockspeed)
     &        * SolarMass / (dens_units * len_units**3) ! gas mass shocked per SNe (Sedov-Taylor)
      tau_ref = dust_growth_tauref * 1e9 * SEC_PER_YEAR_GRACKLE / time_units

      do n=0+1,NUM_METAL_SPECIES_GRACKLE+1
          dM(n) = 0.
      enddo

!   metals accreted from gas to dust
      tau_accr0 = tau_ref * (dust_growth_densref/dens_proper)
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
          endif 

          dM(1) = dM(1) + dM(n)
      enddo

      dM_tau_accr = dM(1)

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
          endif
      enddo
    
!   recalculate metallicity
      do n=0+1,NUM_METAL_SPECIES_GRACKLE+1 
!        1. Place conditions on dM such that mass is conserved but dust and gas
!        mass cannot be negative
!        2. Make sure that only 90% of the gaseous metals can become dust
          dM(n) = max(-1*dmetpt(n), dM(n))
          dM(n) = min(0.9 * metpt(n), dM(n))

!         Transfer metals between gaseous and dust phase only if dust present
          dM_conserv = 0.0
          if (dust_Density(i,j,k) .gt. 0.0) then
              dmetpt(n) = dmetpt(n) + dM(n)
              metpt(n) = metpt(n) - dM(n)
          else
!             If dust density is <= 0 then makes sure metals are conserved accordingly
              dM_conserv = dmetpt(n)
              dmetpt(n) = dmetpt(n) - dM_conserv
              metpt(n) = metpt(n) + dM_conserv
          endif
!        surely this should be only happen if you add to the dust metallicity,
!        else you are removing gas metallicity and breaking conservation
      !     metpt(n) = metpt(n) - dM(n)
      enddo

!     update gas density to account for changes in metal density
      gas_Density(i,j,k) = gas_Density(i,j,k) + (metpt(1) - metal_Density(i,j,k))

      if (metpt(1) .lt. 0) then
            write(*,*) "metpt(1)=", metpt(1),
     &      "dM(1)=", dM(1)
      endif

      if (dmetpt(1) .lt. 0) then
            write(*,*) "final metal density = ", metpt(1)
            write(*,*) "inital metal density   = ", metal_Density(i,j,k)
            write(*,*) "DELTA metal density   = ", -1*dM(1)
            write(*,*) "final dust density = ", dmetpt(1)
            write(*,*) "initial dust density = ", dust_Density(i,j,k)
            write(*,*) "DELTA dust density = ", dM(1)
            write(*,*) "------------"
            write(*,*) "dM_tau_accr = ", dM_tau_accr
            call exit(21)
      endif

      total_density_final = metpt(1) + dmetpt(1)

      if (abs(total_density_final - total_density_init) .gt. 1d-8) then
            write(*,*) "Total dust + metal density not conserved!"
            write(*,*) "Initial value --> ", total_density_init
            write(*,*) "Final value --> ", total_density_final
            write(*,*) "dmetpt = ", dmetpt(1), ", metpt = ", metpt(1), ", dM = ", dM(1)
            write(*,*) "Ending run!!"
            call exit(21)
      endif

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
