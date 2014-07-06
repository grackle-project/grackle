      program fortran_example

      integer iresult

c     Initialization parameters

      integer*8 comoving_coordinates, use_grackle,
     &     with_radiative_cooling, primordial_chemistry,
     &     metal_cooling, h2_on_dust, cmb_temperature_floor

      real*8 density_units, length_units, time_units, velocity_units,
     &     a_units, a_value

      character*36 grackle_data_file

c     Parameters

c     chemistry on
      use_grackle = 1
c     cooling on
      with_radiative_cooling = 1
c     molecular network with H, He, D
      primordial_chemistry = 3
c     metal cooling on
      metal_cooling = 1
c     UV background on
      UVbackground = 1
c     cooling data for Haardt & Madau 2012 background
      grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"

c     Units

      comoving_coordinates = 1
      density_units = 1.67d-24
      length_units = 1.0
      time_units = 1.0d12
      velocity_units = length_units / time_units
      a_units = 1.0

      iresult = initialize_grackle(
     &     comoving_coordinates,
     &     density_units, length_units,
     &     time_units, velocity_units,
     &     a_units, a_value,
     &     use_grackle, with_radiative_cooling,
     &     grackle_data_file,
     &     primordial_chemistry, metal_cooling,
     &     h2_on_dust, cmb_temperature_floor)

      end program fortran_example
