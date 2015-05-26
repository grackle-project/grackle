      program fortran_example

      implicit NONE
#include "grackle_fortran_types.def"

      integer iresult, i

      integer initialize_grackle,
     &     solve_chemistry_table,
     &     calculate_cooling_time_table,
     &     calculate_temperature_table,
     &     calculate_pressure_table

      real*8 mh, kboltz, fH
      R_PREC tiny_number
      parameter (tiny_number = 1.0e-20_RKIND)
      parameter (mh = 1.67262171d-24)
      parameter (kboltz = 1.3806504d-16)
      parameter (fH = 0.76)

c     Initialization parameters

      integer comoving_coordinates, use_grackle,
     &     with_radiative_cooling, primordial_chemistry,
     &     metal_cooling, UVbackground, h2_on_dust, 
     &     cmb_temperature_floor

      real*8 density_units, length_units, time_units, velocity_units,
     &     a_units, a_value, initial_redshift, my_gamma

      character*36 grackle_data_file

c     Field data arrays

      integer field_size
      parameter (field_size = 1)

      real*8 temperature_units, dt

      R_PREC density(field_size), energy(field_size),
     &     x_velocity(field_size), y_velocity(field_size),
     &     z_velocity(field_size), metal_density(field_size)

      R_PREC cooling_time(field_size), gamma(field_size),
     &     pressure(field_size), temperature(field_Size)

c     Grid size and dimension
c     grid_start and grid_end are used to ignore ghost zones.

      integer grid_rank, grid_dimension(3), 
     &     grid_start(3), grid_end(3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Set parameters

c     chemistry on
      use_grackle = 1
c     cooling on
      with_radiative_cooling = 1
c     fully tabular cooling
      primordial_chemistry = 0
c     metal cooling on
      metal_cooling = 1
c     UV background on
      UVbackground = 1
c     cooling data for Haardt & Madau 2012 background
      grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"
c     no dust
      h2_on_dust = 0
c     include CMB cooling floor
      cmb_temperature_floor = 1
c     set adiabatic index
      my_gamma = 5./3.

c     Set units

      comoving_coordinates = 0
      density_units = 1.67d-24
      length_units = 1.0d0
      time_units = 1.0d12
      velocity_units = length_units / time_units
      a_units = 1.0d0

c     Set initial expansion factor (for internal units).
c     Set expansion factor to 1 for non-cosmological simulation.
      initial_redshift = 0.;
      a_value = 1. / (1. + initial_redshift);

c     Initialize the Grackle

      iresult = initialize_grackle(
     &     comoving_coordinates,
     &     density_units, length_units,
     &     time_units, velocity_units,
     &     a_units, a_value,
     &     use_grackle, with_radiative_cooling,
     &     grackle_data_file,
     &     primordial_chemistry, metal_cooling,
     &     UVbackground, h2_on_dust,
     &     cmb_temperature_floor, my_gamma)

c     Set field arrays

c     If grid rank is less than 3, set the other dimensions, 
c     start indices, and end indices to 0.
      grid_rank = 3
      do i = 1, grid_rank
         grid_dimension(i) = 1
         grid_start(i) = 0
         grid_end(i) = 0
      enddo
      grid_dimension(1) = field_size
c     0-based
      grid_end(1) = field_size - 1

c     set temperature units
      temperature_units =  mh *
     &     (a_units * length_units / time_units)**2 / kboltz

      do i = 1,field_size
         density(i) = 1.0
c        solar metallicity
         metal_density(i) = 0.0204 * density(i)

         x_velocity(i) = 0.0
         y_velocity(i) = 0.0
         z_velocity(i) = 0.0

c        initilize internal energy (here 1000 K for no reason)
         energy(i) = 1000. / temperature_units
      enddo

c
c     Calling the chemistry solver
c     These routines can now be called during the simulation.
c

c     Evolving the chemistry.

c     some timestep
      dt = 3.15e7 * 1e6 / time_units

      iresult = solve_chemistry_table(
     &     comoving_coordinates,
     &     density_units, length_units,
     &     time_units, velocity_units,
     &     a_units, a_value, dt,
     &     grid_rank, grid_dimension,
     &     grid_start, grid_end,
     &     density, energy,
     &     x_velocity, y_velocity, z_velocity,
     &     metal_density)

c     Calculate cooling time.

      iresult = calculate_cooling_time_table(
     &     comoving_coordinates,
     &     density_units, length_units,
     &     time_units, velocity_units,
     &     a_units, a_value,
     &     grid_rank, grid_dimension,
     &     grid_start, grid_end,
     &     density, energy,
     &     x_velocity, y_velocity, z_velocity,
     &     metal_density, cooling_time)

      write(*,*) "Cooling time = ", (cooling_time(1) * time_units), "s."

c     Calculate temperature.

      iresult = calculate_temperature_table(
     &     comoving_coordinates,
     &     density_units, length_units,
     &     time_units, velocity_units,
     &     a_units, a_value,
     &     grid_rank, grid_dimension,
     &     grid_start, grid_end,
     &     density, energy,
     &     metal_density, temperature)

      write(*,*) "Temperature = ", temperature(1), "K."

c     Calcualte pressure.

      iresult = calculate_pressure_table(
     &     comoving_coordinates,
     &     density_units, length_units,
     &     time_units, velocity_units,
     &     a_units, a_value,
     &     grid_rank, grid_dimension,
     &     grid_start, grid_end,
     &     density, energy,
     &     pressure)

      write(*,*) "Pressure = ", pressure(1)

      end program fortran_example
