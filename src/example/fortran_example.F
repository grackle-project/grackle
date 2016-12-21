      program fortran_example

      USE ISO_C_BINDING
      implicit None
#include "grackle.def"

      integer iresult, i

c     Define constants

      real*8 mh, kboltz, fH
      R_PREC tiny_number
      parameter (tiny_number = 1.0e-20_RKIND)
      parameter (mh = 1.67262171d-24)
      parameter (kboltz = 1.3806504d-16)
      parameter (fH = 0.76)

c     Initialization parameters

      real*8 initial_redshift
      character(len=80), TARGET :: filename

c     Field data arrays

      integer field_size
      parameter (field_size = 1)

      real*8 temperature_units, dt

      R_PREC, TARGET :: density(field_size), energy(field_size),
     &     x_velocity(field_size), y_velocity(field_size),
     &     z_velocity(field_size),
     &     HI_density(field_size), HII_density(field_size),
     &     HM_density(field_size),
     &     HeI_density(field_size), HeII_density(field_size),
     &     HeIII_density(field_size),
     &     H2I_density(field_size), H2II_density(field_size),
     &     DI_density(field_size), DII_density(field_size),
     &     HDI_density(field_size),
     &     e_density(field_size), metal_density(field_size),
     &     volumetric_heating_rate(field_size),
     &     specific_heating_rate(field_size),
     &     RT_HI_ionization_rate(field_size), 
     &     RT_HeI_ionization_rate(field_size),
     &     RT_HeII_ionization_rate(field_size),
     &     RT_H2_dissociation_rate(field_size),
     &     RT_heating_rate(field_size)

      R_PREC, TARGET :: cooling_time(field_size), gamma(field_size),
     &     pressure(field_size), temperature(field_Size)

c     Grid size and dimension
c     grid_start and grid_end are used to ignore ghost zones.
c     grid_dx is used in H2 self-shielding approximation only

      INTEGER, TARGET :: grid_rank, grid_dimension(3), 
     &     grid_start(3), grid_end(3)

      real*8 grid_dx

c     Define storage for grackle units, fields and parameter data

      TYPE (grackle_units) :: my_units
      TYPE (grackle_field_data) :: my_fields
      TYPE (grackle_chemistry_data) :: grackle_data

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Create a grackle chemistry object for parameters and set defaults
  
      iresult = set_default_chemistry_parameters(grackle_data)

c     Set parameters

      grackle_data%use_grackle = 1            ! chemistry on
      grackle_data%with_radiative_cooling = 1 ! cooling on
      grackle_data%primordial_chemistry = 3   ! network with H, He, D
      grackle_data%metal_cooling = 1          ! metal cooling on
      grackle_data%UVbackground = 1           ! UV background on
c     cooling data for Haardt & Madau 2012 background
      filename = "../../input/CloudyData_UVB=HM2012.h5"//C_NULL_CHAR
      grackle_data%grackle_data_file = C_LOC(filename(1:1))
      grackle_data%h2_on_dust = 0             ! no dust
      grackle_data%cmb_temperature_floor = 1  ! include CMB cooling floor
      grackle_data%Gamma = 5./3.;          ! monoatomic gas

c     Set units

      my_units%comoving_coordinates = 0
      my_units%density_units = 1.67d-24
      my_units%length_units = 1.0d0
      my_units%time_units = 1.0d12
      my_units%velocity_units = my_units%length_units / 
     &     my_units%time_units
      my_units%a_units = 1.0d0

c     Set initial expansion factor (for internal units).
c     Set expansion factor to 1 for non-cosmological simulation.
      initial_redshift = 0.;
      my_units%a_value = 1. / (1. + initial_redshift);

c     Initialize the Grackle

      write(6,*) "primordial_chemistry:",
     &     grackle_data%primordial_chemistry
      write(6,*) "metal_cooling:",
     &     grackle_data%metal_cooling
      iresult = initialize_chemistry_data(my_units)

c     Set field arrays

c     If grid rank is less than 3, set the other dimensions, 
c     start indices, and end indices to 0.
      grid_rank = 3
      do i = 1, grid_rank
         grid_dimension(i) = 1
         grid_start(i) = 0
         grid_end(i) = 0
      enddo
      grid_dx = 0.0
      grid_dimension(1) = field_size
c     0-based
      grid_end(1) = field_size - 1

c     set temperature units
      temperature_units =  mh *
     &     (my_units%a_units * my_units%length_units / 
     &      my_units%time_units)**2 / kboltz

      do i = 1,field_size
         density(i) = 1.0
         HI_density(i) = fH * density(i)
         HII_density(i) = tiny_number * density(i)
         HM_density(i) = tiny_number * density(i)
         HeI_density(i) = (1.0 - fH) * density(i)
         HeII_density(i) = tiny_number * density(i)
         HeIII_density(i) = tiny_number * density(i)
         H2I_density(i) = tiny_number * density(i)
         H2II_density(i) = tiny_number * density(i)
         DI_density(i) = 2.0 * 3.4e-5 * density(i)
         DII_density(i) = tiny_number * density(i)
         HDI_density(i) = tiny_number * density(i)
         e_density(i) = tiny_number * density(i)
c        solar metallicity
         metal_density(i) = grackle_data%SolarMetalFractionByMass *
     &        density(i)

         x_velocity(i) = 0.0
         y_velocity(i) = 0.0
         z_velocity(i) = 0.0

c        initilize internal energy (here 1000 K for no reason)
         energy(i) = 1000. / temperature_units

         volumetric_heating_rate(i) = 0.0
         specific_heating_rate(i) = 0.0
         RT_HI_ionization_rate(i) = 0.0
         RT_HeI_ionization_rate(i) = 0.0
         RT_HeII_ionization_rate(i) = 0.0
         RT_H2_dissociation_rate(i) = 0.0
         RT_heating_rate(i) = 0.0
      enddo
c
c     Fill in structure to be passed to Grackle
c
      my_fields%grid_rank = 1
      my_fields%grid_dimension = C_LOC(grid_dimension)
      my_fields%grid_start = C_LOC(grid_start)
      my_fields%grid_end = C_LOC(grid_end)
      my_fields%grid_dx  = grid_dx

      my_fields%density = C_LOC(density)
      my_fields%HI_density = C_LOC(HI_density)
      my_fields%HII_density = C_LOC(HII_density)
      my_fields%HM_density = C_LOC(HM_density)
      my_fields%HeI_density = C_LOC(HeI_density)
      my_fields%HeII_density = C_LOC(HeII_density)
      my_fields%HeIII_density = C_LOC(HeIII_density)
      my_fields%H2I_density = C_LOC(H2I_density)
      my_fields%H2II_density = C_LOC(H2II_density)
      my_fields%DI_density = C_LOC(DI_density)
      my_fields%DII_density = C_LOC(DII_density)
      my_fields%HDI_density = C_LOC(HDI_density)
      my_fields%e_density = C_LOC(e_density)
      my_fields%metal_density = C_LOC(metal_density)
      my_fields%internal_energy = C_LOC(energy)
      my_fields%x_velocity = C_LOC(x_velocity)
      my_fields%y_velocity = C_LOC(y_velocity)
      my_fields%z_velocity = C_LOC(z_velocity)
      my_fields%volumetric_heating_rate = 
     &                                C_LOC(volumetric_heating_rate)
      my_fields%specific_heating_rate = C_LOC(specific_heating_rate)
      my_fields%RT_HI_ionization_rate = C_LOC(RT_HI_ionization_rate)
      my_fields%RT_HeI_ionization_rate = C_LOC(RT_HeI_ionization_rate)
      my_fields%RT_HeII_ionization_rate = 
     &                                  C_LOC(RT_HeII_ionization_rate)
      my_fields%RT_H2_dissociation_rate = C_LOC(RT_H2_dissociation_rate)
      my_fields%RT_heating_rate = C_LOC(RT_heating_rate)
c
c     Calling the chemistry solver
c     These routines can now be called during the simulation.
      
c     Evolving the chemistry.

      dt = 3.15e7 * 1e6 / my_units%time_units    ! some timestep
      iresult = solve_chemistry(my_units, my_fields, dt)

c     Calculate cooling time.

      iresult = calculate_cooling_time(my_units, my_fields,cooling_time)
      write(*,*) "Cooling time = ", (cooling_time(1) * 
     &     my_units%time_units), "s."

c     Calculate temperature.

      iresult = calculate_temperature(my_units, my_fields, temperature)
      write(*,*) "Temperature = ", temperature(1), "K."

c     Calcualte pressure.

      iresult = calculate_pressure(my_units, my_fields, pressure)
      write(*,*) "Pressure = ", pressure(1)

c     Calculate gamma.

      iresult = calculate_gamma(my_units, my_fields, gamma)
      write(*,*) "Gamma = ", gamma(1)

      end program fortran_example
