c     Problems relating to grackle_data%grackle_data_file arise when
c     this program is compiled using certain versions of gfortran
c     (including versions 9.2.1 & 9.4.0) with an optimization level of
c     -O2 or higher.
c
c     In these cases, the compiler appears to make an aggressive
c     optimization in which the program deallocates the `filename`
c     variable before calls to certain grackle routines. This is
c     problematic because grackle_data%grackle_data_file points to the
c     address of the data held in `filename`

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

      real*8 temperature_units, pressure_units, dt

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
     &     dust_density(field_size),
     &     volumetric_heating_rate(field_size),
     &     specific_heating_rate(field_size),
     &     RT_HI_ionization_rate(field_size), 
     &     RT_HeI_ionization_rate(field_size),
     &     RT_HeII_ionization_rate(field_size),
     &     RT_H2_dissociation_rate(field_size),
     &     RT_heating_rate(field_size),
     &     isrf_habing(field_size)

      R_PREC, TARGET :: cooling_time(field_size), gamma(field_size),
     &     pressure(field_size), temperature(field_Size),
     &     dust_temperature(field_Size)

c     Grid size and dimension
c     grid_start and grid_end are used to ignore ghost zones.
c     grid_dx is used in H2 self-shielding approximation only

      INTEGER, TARGET :: grid_dimension(3),
     &     grid_start(3), grid_end(3)

c     Define storage for grackle units, fields and parameter data

      TYPE (grackle_units) :: my_units
      TYPE (grackle_field_data) :: my_fields
      TYPE (grackle_chemistry_data) :: grackle_data

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Set initial redshift (for internal units).
      initial_redshift = 0.

c     First, set up the units system.
c     These are conversions from code units to cgs.

      my_units%comoving_coordinates = 0
      my_units%density_units = 1.67d-24
      my_units%length_units = 1.0d0
      my_units%time_units = 1.0d12
      my_units%a_units = 1.0d0
c     Set expansion factor to 1 for non-cosmological simulation.
      my_units%a_value = 1. / (1. + initial_redshift) / my_units%a_units

c     The following shouldn't be necessary (it's just here for
c     historical reasons)
      call set_velocity_units(my_units)

c     Create a grackle chemistry object for parameters and set defaults
  
      iresult = set_default_chemistry_parameters(grackle_data)
      if (iresult .ne. 1) then
        write(0,*) "Error in set_default_chemistry_parameters."
        STOP 1 ! the standard recommends that this is the return code
      endif

c     Set parameters

      grackle_data%use_grackle = 1            ! chemistry on
      grackle_data%with_radiative_cooling = 1 ! cooling on
      grackle_data%primordial_chemistry = 3   ! network with H, He, D
      grackle_data%metal_cooling = 1          ! metal cooling on
      grackle_data%UVbackground = 1           ! UV background on
      grackle_data%dust_chemistry = 1         ! dust processes
      grackle_data%use_dust_density_field = 1 ! follow dust density field
      grackle_data%use_isrf_field = 1         ! follow interstellar radiation field
      filename = "../../input/CloudyData_UVB=HM2012.h5"//C_NULL_CHAR
      grackle_data%grackle_data_file = C_LOC(filename(1:1))


c     Initialize the Grackle

      write(6,*) "primordial_chemistry:",
     &     grackle_data%primordial_chemistry
      write(6,*) "metal_cooling:",
     &     grackle_data%metal_cooling

c     Finally, initialize the chemistry object
      iresult = initialize_chemistry_data(my_units)
      if (iresult .ne. 1) then
        write(0,*) "Error in initialize_chemistry_data."
        STOP 1 ! the standard recommends that this is the return code
      endif

c     Prepare the struct for storing grackle field data
      iresult = gr_initialize_field_data(my_fields)
      if (iresult .ne. 1) then
        write(0,*) "Error in gr_initialize_field_data."
        STOP 1 ! the standard recommends that this is the return code
      endif


c     If grid rank is less than 3, set the other dimensions, 
c     start indices, and end indices to 0.
      my_fields%grid_rank = 3
      my_fields%grid_dimension = C_LOC(grid_dimension)
      my_fields%grid_start = C_LOC(grid_start)
      my_fields%grid_end = C_LOC(grid_end)
c     The following is only used for H2 self-shielding approximation
c     (a value of -1 indicates it is unset)
      my_fields%grid_dx  = -1.0

      do i = 1, my_fields%grid_rank
         grid_dimension(i) = 1
         grid_start(i) = 0
         grid_end(i) = 0
      enddo
      grid_dimension(1) = field_size
      grid_end(1) = field_size - 1 ! <- 0-based index

      my_fields%density = C_LOC(density)
      my_fields%internal_energy = C_LOC(energy)
      my_fields%x_velocity = C_LOC(x_velocity)
      my_fields%y_velocity = C_LOC(y_velocity)
      my_fields%z_velocity = C_LOC(z_velocity)
c     for primordial_chemistry >= 1
      my_fields%HI_density = C_LOC(HI_density)
      my_fields%HII_density = C_LOC(HII_density)
      my_fields%HeI_density = C_LOC(HeI_density)
      my_fields%HeII_density = C_LOC(HeII_density)
      my_fields%HeIII_density = C_LOC(HeIII_density)
      my_fields%e_density = C_LOC(e_density)
c     for primordial_chemistry >= 2
      my_fields%HM_density = C_LOC(HM_density)
      my_fields%H2I_density = C_LOC(H2I_density)
      my_fields%H2II_density = C_LOC(H2II_density)
c     for primordial_chemistry >= 3
      my_fields%DI_density = C_LOC(DI_density)
      my_fields%DII_density = C_LOC(DII_density)
      my_fields%HDI_density = C_LOC(HDI_density)
c     for metal_cooling = 1
      my_fields%metal_density = C_LOC(metal_density)
c     for use_dust_density_field = 1
      my_fields%dust_density = C_LOC(dust_density)

c     volumetric heating rate (provide in units [erg s^-1 cm^-3])
      my_fields%volumetric_heating_rate =
     &                                C_LOC(volumetric_heating_rate)
c     specific heating rate (provide in units [egs s^-1 g^-1])
      my_fields%specific_heating_rate = C_LOC(specific_heating_rate)

c     radiative transfer ionization / dissociation rate fields
c     (provide in units of [1/s])
      my_fields%RT_HI_ionization_rate = C_LOC(RT_HI_ionization_rate)
      my_fields%RT_HeI_ionization_rate = C_LOC(RT_HeI_ionization_rate)
      my_fields%RT_HeII_ionization_rate =
     &                                  C_LOC(RT_HeII_ionization_rate)
      my_fields%RT_H2_dissociation_rate = C_LOC(RT_H2_dissociation_rate)

c     radiative transfer heating rate field
c     (provide in units of [erg s^-1 cm^-3])
      my_fields%RT_heating_rate = C_LOC(RT_heating_rate)

c     interstellar radiation field strength
      my_fields%isrf_habing = C_LOC(isrf_habing)

c     get the temperature_units
      temperature_units = get_temperature_units(my_units)

      do i = 1,field_size
         density(i) = 1.0
         HI_density(i) = grackle_data%HydrogenFractionByMass *
     &        density(i)
         HII_density(i) = tiny_number * density(i)
         HM_density(i) = tiny_number * density(i)
         HeI_density(i) = (1.0 - grackle_data%HydrogenFractionByMass) *
     &        density(i)
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
c        local dust to gas ratio
         dust_density(i) = grackle_data%local_dust_to_gas_ratio *
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

         isrf_habing(i) = grackle_data%interstellar_radiation_field
      enddo

c
c     Calling the chemistry solver
c     These routines can now be called during the simulation.
      
c     Evolving the chemistry.

      dt = 3.15e7 * 1e6 / my_units%time_units    ! some timestep
      iresult = solve_chemistry(my_units, my_fields, dt)
      if (iresult .ne. 1) then
        write(0,*) "Error in solve_chemistry."
        STOP 1 ! the standard recommends that this is the return code
      endif


c     Now we compute derived quantities and print the values
c     -> the Fortran formatting syntax is a little strange. It's unclear
c        why the G formatter insists on putting a zero to the left of
c        the decimal point for really small and really large numbers

c     Calculate cooling time.

      iresult = calculate_cooling_time(my_units, my_fields,cooling_time)
      if (iresult .ne. 1) then
        write(0,*) "Error in calculate_cooling_time."
        STOP 1 ! the standard recommends that this is the return code
      endif
      write(6,'("cooling_time = ",G25.17E3," s")')
     &     cooling_time(1) * my_units%time_units

c     Calculate temperature.

      iresult = calculate_temperature(my_units, my_fields, temperature)
      if (iresult .ne. 1) then
        write(0,*) "Error in calculate_temperature."
        STOP 1 ! the standard recommends that this is the return code
      endif
      write(6,'("temperature = ",G25.17E3," K")') temperature(1)

c     Calcualte pressure.

      pressure_units = my_units%density_units *
     &     my_units%velocity_units**2
      iresult = calculate_pressure(my_units, my_fields, pressure)
      if (iresult .ne. 1) then
        write(0,*) "Error in calculate_pressure."
        STOP 1 ! the standard recommends that this is the return code
      endif
      write(6,'("pressure = ",G25.17E3," dyne/cm^2")')
     &     pressure(1)*pressure_units

c     Calculate gamma.

      iresult = calculate_gamma(my_units, my_fields, gamma)
      if (iresult .ne. 1) then
        write(0,*) "Error in calculate_gamma."
        STOP 1 ! the standard recommends that this is the return code
      endif
      write(6,'("gamma = ",G25.17E3)') gamma(1)

c     Calculate dust temperature.

      iresult = calculate_dust_temperature(my_units, my_fields,
     &     dust_temperature)
      if (iresult .ne. 1) then
        write(0,*) "Error in calculate_dust_temperature."
        STOP 1 ! the standard recommends that this is the return code
      endif
      write(6,'("dust_temperature = ",G25.17E3," K")')
     &     dust_temperature(1)


      end program fortran_example
