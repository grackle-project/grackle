.. _integration:

Adding Grackle to Your Simulation Code
======================================

Example Executable
------------------

The grackle source code contains a C++ example that links against the 
grackle library.  The example is located in the **src/example** directory 
and is called **example.C**.  If you have already installed the grackle 
library, you should be able to build the example by simply typing *make*.  
To run the example, copy the file, **CloudyData_UVB=HM2012.h5**, from the 
**input** directory to the same directory as the example executable.  
Also, make sure to add the path to the directory containing the installed 
**libgrackle.so** to your LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH on Mac).

Header Files
------------

Four source files are installed with the grackle library.  They are:

    * **grackle.h** - the primary header file, containing declarations for all the available functions and data structures.  This is the only header file that needs to be included.

    * **grackle_macros.h** - this contains basic variable type definitions.

    * **chemistry_data.h** - this defines the primary data structure which all run time parameters as well as the chemistry, cooling, and UV background data.

    * **code_units.h** - this defines the structure containing conversions from code units to CGS.

The only source file that needs to be included in your simulation code is 
**grackle.h**.

Data Types
----------

The grackle library provied two variable sized data types, one for integers 
and one for floating point variables.  With **grackle.h** included, both of 
these data types are available.

    * *gr_int* - the integer data type.  This type is a 32 bit integer (int) if compiled with *integers-32* and a 64 bit integer (long int) if compiled with *integers-64*.

    * *gr_float* - the floating point data type.  This type is a 32 bit float (float) if compiled with *precisions-32* and a 64 bit float (double) if compiled with *precision-64*.

Code Units
----------

The *code_units* structure contains conversions from code units to CGS.  
If the simulation is cosmological, these conversion factors will convert from 
*comoving code units* to *proper CGS units*.  These must be set and updated 
manually.

.. code-block:: c++

  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24; // 1 m_H/cc
  my_units.length_units = 3.086e21;  // 1 kpc
  my_units.time_units = 3.15569e13;  // 1 Myr
  my_units.a_units = 1.0;            // units for the expansion factor

Units for the expansion factor are optional.  If used, they are typically 
defined as 1 / (1 + z\ :sub:`i`), where z\ :sub:`i` is the initial redshift 
of the simulation.  Following these definitions, the velocity and temperature 
units are defined as:

    * *velocity_units* = *a_units* x *length_units* / *time_units*
    * *temperature_units* = m\ :sub:`H` * *velocity_units* \ :sup:`2` / k,

where k is the Boltzmann constant.  Temperature units convert from energy to 
temperature in K.

Chemistry Data
--------------

The *chemistry_data* structure contains all of the parameters for controlling 
the behavior of the chemistry and cooling solver.  It also contains all of the 
actual chemistry and cooling rate data.  The routine, 
*set_default_chemistry_parameters* creates the *chemistry_data* structure 
with the default settings and returns it.  The parameters can then be set to 
their desired values.  See :ref:`parameters` for a full list of the available 
parameters.

.. code-block:: c++

  chemistry_data my_chemistry = set_default_chemistry_parameters();
  // Set parameter values for chemistry.
  my_chemistry.use_chemistry = 1;          // chemistry on
  my_chemistry.with_radiative_cooling = 1; // cooling on
  my_chemistry.primordial_chemistry = 3;   // molecular network with H, He, D
  my_chemistry.metal_cooling = 1;          // metal cooling on
  my_chemistry.UVbackground = 1;           // UV background on
  my_chemistry.include_metal_heating = 1;  // heating on metals by UVB on
  my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"; // data file

Once the desired parameters have been set, the chemistry and cooling rates 
must be initialized with the *initialize_chemistry_data*.  This function 
also requires the initial value of the expansion factor for setting internal 
units.  If the simulation is not cosmological, the expansion factor should be 
set to 1.  The initializing function will return an integer indicating success 
(1) or failure (0).

.. code-block:: c++

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  gr_float initial_redshift = 100.;
  gr_float a_value = 1. / (1. + initial_redshift);

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(my_chemistry, my_units, a_value) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return 0;
  }

The *chemistry_data* structure is now ready to be used.

Creating the Necessary Fields
-----------------------------

With the *code_units* and *chemistry_data* structures ready, the only thing 
left is to create the arrays to carry the species densities.  Pointers for all 
fields must be created, but the arrays only need to be allocated if the fields 
are going to be used by the chemistry network.  Variables containing the 
dimensionality of the data, the active dimensions (not including the ghost 
zones), and the starting and ending indices for each dimensions must also be 
created.

.. code-block:: c++

  // Allocate field arrays.
  gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity,
    *HI_density, *HII_density, *HM_density,
    *HeI_density, *HeII_density, *HeIII_density,
    *H2I_density, *H2II_density,
    *DI_density, *DII_density, *HDI_density,
    *e_density, *metal_density;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  gr_int field_size = 10;
  gr_int grid_rank = 3;
  // If grid rank is less than 3, set the other dimensions, 
  // start indices, and end indices to 0.
  gr_int grid_dimension[3], grid_start[3], grid_end[3];
  for (int i = 0;i < 3;i++) {
    grid_dimension[i] = 0; // the active dimension not including ghost zones.
    grid_start[i] = 0;
    grid_end[i] = 0;
  }
  grid_dimension[0] = field_size;
  grid_end[0] = field_size - 1;

  density       = new gr_float[field_size];
  energy        = new gr_float[field_size];
  x_velocity    = new gr_float[field_size];
  y_velocity    = new gr_float[field_size];
  z_velocity    = new gr_float[field_size];
  // for primordial_chemistry >= 1
  HI_density    = new gr_float[field_size];
  HII_density   = new gr_float[field_size];
  HeI_density   = new gr_float[field_size];
  HeII_density  = new gr_float[field_size];
  HeIII_density = new gr_float[field_size];
  e_density     = new gr_float[field_size];
  // for primordial_chemistry >= 2
  HM_density    = new gr_float[field_size];
  H2I_density   = new gr_float[field_size];
  H2II_density  = new gr_float[field_size];
  // for primordial_chemistry >= 3
  DI_density    = new gr_float[field_size];
  DII_density   = new gr_float[field_size];
  HDI_density   = new gr_float[field_size];
  // for metal_cooling = 1
  metal_density = new gr_float[field_size];

Calling the Available Functions
-------------------------------

There are five functions available, one to solve the chemistry and cooling 
and four others to calculate the cooling time, temperature, pressure, and the 
ratio of the specific heats (gamma).  The arguments required are the 
*code_units* and *chemistry_data* structures, the field size and dimension 
variables, and the field arrays themselves.  In some cases, the current value 
of the expansion factor must also be given and for the chemistry solving 
routine, a timestep must be given.  For the four field calculator routines, 
the array to be filled with the field values must be created and passed as an 
argument as well.

Before these are called, the UV background 
rates must be updated.  This need only be done once per timestep.

Updating the UV Background
++++++++++++++++++++++++++

.. code-block:: c++

  if (update_UVbackground_rates(my_chemistry, 
                                my_units, a_value) == 0) {
    fprintf(stderr, "Error in update_UBbackground_rates.\n");
    return 0;
  }

Solve the Chemistry and Cooling
+++++++++++++++++++++++++++++++

.. code-block:: c++

  // some timestep (one million years)
  gr_float dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(my_chemistry, my_units,
                      a_value, dt,
                      grid_rank, grid_dimension,
                      grid_start, grid_end,
                      density, energy,
                      x_velocity, y_velocity, z_velocity,
                      HI_density, HII_density, HM_density,
                      HeI_density, HeII_density, HeIII_density,
                      H2I_density, H2II_density,
                      DI_density, DII_density, HDI_density,
                      e_density, metal_density) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return 0;
  }

Calculating the Cooling Time
++++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *cooling_time;
  cooling_time = new gr_float[field_size];
  if (calculate_cooling_time(my_chemistry, my_units,
                             a_value,
                             grid_rank, grid_dimension,
                             grid_start, grid_end,
                             density, energy,
                             x_velocity, y_velocity, z_velocity,
                             HI_density, HII_density, HM_density,
                             HeI_density, HeII_density, HeIII_density,
                             H2I_density, H2II_density,
                             DI_density, DII_density, HDI_density,
                             e_density, metal_density, 
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return 0;
  }

Calculating the Temperature Field
+++++++++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *temperature;
  temperature = new gr_float[field_size];
  if (calculate_temperature(my_chemistry, my_units,
                            grid_rank, grid_dimension,
                            density, energy,
                            HI_density, HII_density, HM_density,
                            HeI_density, HeII_density, HeIII_density,
                            H2I_density, H2II_density,
                            DI_density, DII_density, HDI_density,
                            e_density, metal_density, 
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return 0;
  }

Calculating the Pressure Field
++++++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *pressure;
  pressure = new gr_float[field_size];
  if (calculate_pressure(my_chemistry, my_units,
                         grid_rank, grid_dimension,
                         density, energy,
                         HI_density, HII_density, HM_density,
                         HeI_density, HeII_density, HeIII_density,
                         H2I_density, H2II_density,
                         DI_density, DII_density, HDI_density,
                         e_density, metal_density,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return 0;
  }

Calculating the Gamma Field
+++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *gamma;
  gamma = new gr_float[field_size];
  if (calculate_gamma(my_chemistry, my_units,
                      grid_rank, grid_dimension,
                      density, energy,
                      HI_density, HII_density, HM_density,
                      HeI_density, HeII_density, HeIII_density,
                      H2I_density, H2II_density,
                      DI_density, DII_density, HDI_density,
                      e_density, metal_density,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return 0;
  }
