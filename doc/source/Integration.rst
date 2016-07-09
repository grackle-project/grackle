.. _integration:

Adding Grackle to Your Simulation Code
======================================

The majority of this document follows the implementation of Grackle in
a C++ simulation code.  For more information on adding Grackle to a
Fortran code, see :ref:`fortran`.  Full implementation examples for
C, C++, and Fortran are also available in the Grackle source.  See
:ref:`examples` for more information.  For a list of all available
functions, see the :ref:`reference`.

.. _examples:

Example Executables
-------------------

The grackle source code contains examples for C, C++, and Fortran codes.  
They are located in the **src/example** directory and provide examples
of calling all of grackle's functions.

    * **c_example.c** - C example

    * **cxx_example.C** - C++ example

    * **cxx_omp_example.C** - C++ example using OpenMP

    * **fortran_example.F** - Fortran example

Once you have already installed the grackle library, you can build the examples 
by typing *make* and the name of the file without extension.  For example, to 
build the C++ example, type:

.. code-block:: bash

  $ make cxx_example

To run the example, make sure to add the path to the directory containing 
the installed **libgrackle.so** to your LD_LIBRARY_PATH (or 
DYLD_LIBRARY_PATH on Mac).

Header Files
------------

Seven header files are installed with the grackle library.  They are:

    * **grackle.h** - the primary header file, containing declarations for all
      the available functions and data structures.  This is the only header
      file that needs to be included for C and C++ codes.

    * **grackle.def** - the header file to be used in Fortran codes.  Only
      this file needs to be included.

    * **grackle_types.h** - defines the variable type :c:type:`gr_float`, the
      field structure :c:type:`grackle_field_data`, and the units structure
      :c:type:`code_units`.

    * **grackle_chemistry_data.h** - defines the :c:type:`chemistry_data`
      structure, which stores all Grackle run-time parameters and the
      :c:type:`chemistry_data_storage` structure, which stores all chemistry
      and cooling rate data.

    * **grackle_fortran_types.def** - similar to **grackle_types.h**, but used
      with Fortran codes.  This defines the variable type :c:type:`R_PREC` as
      either real\*4 or real\*8.

    * **grackle_fortran_interface.def** - defines the Fortran interface,
      including the Fortran analogs of :c:type:`grackle_field_data`,
      :c:type:`code_units`, and :c:type:`grackle_chemistry_data`.

    * **grackle_macros.h** - contains some macros used internally.

For C and C++ codes, the only source file that needs to be included in your
simulation code is **grackle.h**.  For Fortran, use **grackle.def**.  Since
Grackle is written in C, including **grackle.h** in a C++ code requires the
*extern "C"* directive.

.. code-block:: c++

  extern "C" {
  #include <grackle.h>
  }

Data Types
----------

The grackle library provides a configurable variable type to control the 
precision of the baryon fields passed to the grackle functions.  For C and 
C++ codes, this is :c:type:`gr_float`.  For Fortran codes, this is
:c:type:`R_PREC`.  The precision of these types can be configured with the      
*precision* compile option.  Compile with *precision-32* to make
:c:type:`gr_float` and :c:type:`R_PREC` a 4 byte float (*float* for C/C++
and *real\*4* for Fortran).  Compile with *precision-64* to make
:c:type:`gr_float` and :c:type:`R_PREC` an 8 byte float (*double* for C/C++
and *real\*8* for Fortran).

.. c:type:: gr_float

   Floating point type used for the baryon fields.  This is of type *float*
   if compiled with *precision-32* and type double if compiled with
   *precision-64*.

.. c:type:: R_PREC

   The Fortran analog of :c:type:`gr_float`.  This is of type *real\*4* if
   compiled with *precision-32* and type *real\*8* if compiled with
   *precision-64*.

Enabling Output
---------------

By default, grackle will not print anything but error messages.  However,
a short summary of the running configuration can be printed by setting
``grackle_verbose`` to 1.  In a parallel code, it is recommended that
output only be enabled for the root process.

.. code-block:: c++

   // Enable output
   grackle_verbose = 1;

Code Units
----------

**It is strongly recommended to use comoving coordinates with any
cosmological simulation.**  The :c:data:`code_units` structure contains
conversions from code units to CGS.  If :c:data:`comoving_coordinates` is set to
0, it is assumed that the fields passed into the solver are in the
proper frame.  All of the units (density, length, time, velocity, and
expansion factor) must be set.  When using the proper frame, :c:data:`a_units`
(units for the expansion factor) must be set to 1.0.

.. c:type:: code_units

   This structure contains the following members.

.. c:var:: int comoving_coordinates

   If set to 1, the incoming field data is assumed to be in the comoving
   frame.  If set to 0, the incoming field data is assumed to be in the
   proper frame.

.. c:var:: double density_units

   Conversion factor to be multiplied by density fields to return
   densities in proper g/cm\ :sup:`3`\.

.. c:var:: double length_units

   Conversion factor to be multiplied by length variables to return
   lengths in proper cm.

.. c:var:: double time_units

   Conversion factor to be multiplied by time variables to return
   times in s.

.. c:var:: double velocity_units

   Conversion factor to be multiplied by velocities to return proper cm/s.

.. c:var:: double a_units

   Conversion factor to be multiplied by the expansion factor such that
   a\ :sub:`true`\  = a\ :sub:`code`\ * :c:data:`a_units`.

.. code-block:: c++

  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24; // 1 m_H/cc
  my_units.length_units = 3.086e21;  // 1 kpc
  my_units.time_units = 3.15569e13;  // 1 Myr
  my_units.velocity_units = my_units.length_units / my_units.time_units;
  my_units.a_units = 1.0;            // units for the expansion factor

If :c:data:`comoving_coordinates` is set to 1, it is assumed that the fields being 
passed to the solver are in the comoving frame.  Hence, the units must 
convert from code units in the **comoving** frame to CGS in the **proper** 
frame.  

.. note:: With :c:data:`comoving_coordinate` set to 1, velocity units need to be
   defined in the following way.

.. code-block:: c++

  my_units.velocity_units = my_units.a_units * 
    (my_units.length_units / a_value) / my_units.time_units; // since u = a * dx/dt

For an example of using comoving units, see the units system in the 
`Enzo <http://enzo-project.org/>`_ code.  For cosmological simualations, a 
comoving unit system is preferred, though not required, since it allows the 
densities to stay close to 1.0.

Chemistry Data
--------------

The main Grackle header file contains a structure of type :c:type:`chemistry_data` 
called ``grackle_data``, which 
contains all of the parameters that control the behavior of the solver as well as 
all of the actual chemistry and cooling rate data.  The routine, 
:c:func:`set_default_chemistry_parameters` is responsible for the initial setup of this 
structure and for setting of all the default parameter values.  The parameters can 
then be set to their desired values.  See :ref:`parameters` for a full list of the 
available parameters.  The function will return an integer indicating success 
(1) or failure (0).

.. c:type:: chemistry_data

   This structure holds all grackle run time parameter and all chemistry and
   cooling data arrays.

.. code-block:: c++

  if (set_default_chemistry_parameters() == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
  }

  // Set parameter values for chemistry.
  grackle_data.use_grackle = 1;            // chemistry on
  grackle_data.with_radiative_cooling = 1; // cooling on
  grackle_data.primordial_chemistry = 3;   // molecular network with H, He, D
  grackle_data.metal_cooling = 1;          // metal cooling on
  grackle_data.UVbackground = 1;           // UV background on
  grackle_data.grackle_data_file = "CloudyData_UVB=HM2012.h5"; // data file

Once the desired parameters have been set, the chemistry and cooling rates 
must be initialized with the :c:func:`initialize_chemistry_data`.  This function 
also requires the initial value of the expansion factor for setting internal 
units.  If the simulation is not cosmological, the expansion factor should be 
set to 1.  The initializing function will return an integer indicating success 
(1) or failure (0).

.. code-block:: c++

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  double initial_redshift = 100.;
  double a_value = 1. / (1. + initial_redshift) / my_units.a_units;

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units, a_value) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return 0;
  }

The Grackle is now ready to be used.

.. _openmp:

Running with OpenMP
-------------------

As of version 2.2, Grackle can be run with OpenMP parallelism.  To do this,
the library must first be compiled with OpenMP support enabled by issuing the
command, "make omp-on", before compiling.  See :ref:`compiler-settings` for
more information on how to change settings.

For an example of how to compile your code with OpenMP, see the
**cxx_table_example.C** code example (:ref:`examples`).  Once your code has
been compiled with OpenMP enabled, the number of threads used can be controlled
by setting the :c:data:`omp_nthreads` parameter, stored in the ``grackle_data``
struct.

.. code-block:: c++

   // 8 threads per process
   grackle_data.omp_nthreads = 8;

If not set, this parameter will be set to the maximum number of threads
possible, as determined by the system or as configured by setting the
``OMP_NUM_THREADS`` environment variable.

Creating the Necessary Fields
-----------------------------

With the :c:data:`code_units` and :c:data:`chemistry_data` structures ready, the only thing 
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
  int field_size = 10;
  int grid_rank = 3;
  // If grid rank is less than 3, set the other dimensions to 1 and  
  // start indices and end indices to 0.
  int grid_dimension[3], grid_start[3], grid_end[3];
  for (int i = 0;i < 3;i++) {
    grid_dimension[i] = 1; // the active dimension not including ghost zones.
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

.. note:: The electron mass density should be scaled by the ratio of the
   proton mass to the electron mass such that the electron density in the
   code is the electron number density times the **proton** mass.

Calling the Available Functions
-------------------------------

There are five functions available, one to solve the chemistry and cooling 
and four others to calculate the cooling time, temperature, pressure, and the 
ratio of the specific heats (gamma).  The arguments required are the 
:c:data:`code_units` structure, the value of the expansion factor, the field size and 
dimension variables, and the field arrays themselves.  For the chemistry solving 
routine, a timestep must also be given.  For the four field calculator routines, 
the array to be filled with the field values must be created and passed as an 
argument as well.

Solve the Chemistry and Cooling
+++++++++++++++++++++++++++++++

.. code-block:: c++

  // some timestep (one million years)
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&my_units, a_value, dt,
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
  if (calculate_cooling_time(&my_units, a_value,
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
  if (calculate_temperature(&my_units, a_value,
                            grid_rank, grid_dimension,
                            grid_start, grid_end,
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
  if (calculate_pressure(&my_units, a_value,
                         grid_rank, grid_dimension,
                         grid_start, grid_end,
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
  if (calculate_gamma(&my_units, a_value,
                      grid_rank, grid_dimension,
                      grid_start, grid_end,
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

.. _fortran:

The Fortran Interface
---------------------
