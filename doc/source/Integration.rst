
.. _integration:

Adding Grackle to Your Simulation Code
======================================

The majority of this document follows the implementation of Grackle in
a C++ simulation code.  Full implementation examples for
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

.. c:var:: double a_value

   The current value of the expansion factor in units of :c:data:`a_units`.
   The conversion from redshift to expansion factor in code units is given
   by :c:data:`a_value` = 1 / (1 + z) / :c:data:`a_units`.  If the
   simulation is not cosmological, :c:data:`a_value` should be set to 1.
   Note, if :c:data:`a_value` is set to something other than 1 in a
   non-cosmological simulation, all redshift dependent chemistry and
   cooling terms will be set corresponding to the redshift given.

.. code-block:: c++

  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24; // 1 m_H/cc
  my_units.length_units = 3.086e21;  // 1 kpc
  my_units.time_units = 3.15569e13;  // 1 Myr
  my_units.velocity_units = my_units.length_units / my_units.time_units;
  my_units.a_units = 1.0;            // units for the expansion factor
  my_units.a_value = 1. / (1. + current_redshift) / my_units.a_units;

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
`Enzo <http://enzo-project.org/>`_ code.  For cosmological simulations, a
comoving unit system is preferred, though not required, since it allows the 
densities to stay close to 1.0.

Chemistry Data
--------------

The main Grackle header file contains a structure of type
:c:type:`chemistry_data` called ``grackle_data``, which contains all of the
parameters that control the behavior of the solver.  The routine,
:c:func:`set_default_chemistry_parameters` is responsible for the initial setup
of this structure and for setting of all the default parameter values.  This
function must be handed a pointer to an instance of :c:type:`chemistry_data`,
which will then be attached to ``grackle_data``.  The function will return an
integer indicating success (1) or failure (0).  After this, parameters can then
be set to their desired values by accessing ``grackle_data``.  See
:ref:`parameters` for a full list of the available parameters.

.. c:type:: chemistry_data

   This structure holds all grackle run-time parameters, which are listed in
   :ref:`parameters`.

.. c:type:: chemistry_data_storage

   This structure holds all chemistry and cooling rate arrays.  All functions
   described here make use of an internally stored instance of this type.
   The user will not normally encounter this data type, except when using the
   :ref:`internal_functions`.

.. code-block:: c++

  chemistry_data *my_grackle_data;
  my_grackle_data = new chemistry_data;
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
  }

  // Set parameter values for chemistry.
  // Now access the global copy of the chemistry_data struct (grackle_data).
  grackle_data->use_grackle = 1;            // chemistry on
  grackle_data->with_radiative_cooling = 1; // cooling on
  grackle_data->primordial_chemistry = 3;   // molecular network with H, He, D
  grackle_data->metal_cooling = 1;          // metal cooling on
  grackle_data->UVbackground = 1;           // UV background on
  grackle_data->grackle_data_file = "CloudyData_UVB=HM2012.h5"; // data file

Once the desired parameters have been set, the chemistry and cooling rates 
must be initialized by calling :c:func:`initialize_chemistry_data` with a
pointer to the :c:data:`code_units` struct created earlier.  This function
will return an integer indicating success (1) or failure (0).

.. code-block:: c++

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units) == 0) {
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
**cxx_omp_example.C** code example (:ref:`examples`).  Once your code has
been compiled with OpenMP enabled, the number of threads used can be controlled
by setting the :c:data:`omp_nthreads` parameter, stored in the ``grackle_data``
struct.

.. code-block:: c++

   // 8 threads per process
   grackle_data->omp_nthreads = 8;

If not set, this parameter will be set to the maximum number of threads
possible, as determined by the system or as configured by setting the
``OMP_NUM_THREADS`` environment variable.

Creating the Necessary Fields
-----------------------------

As of version 3.0, the various density and energy fields are passed to
Grackle's functions using a struct of type :c:data:`grackle_field_data`.
The struct contains information about the size and shape of the field arrays
and pointers to all field arrays.

.. c:type:: grackle_field_data

   This structure is used to pass field data to Grackle's functions.  It
   contains the following members:

.. c:var:: int grid_rank

   The active dimensions (not including ignored boundary zones) of the field
   arrays.

.. c:var:: int* grid_dimension

   This should point to an array of size :c:data:`grid_rank`.  This stores
   the size of the field arrays in each dimension.

.. c:var:: int* grid_start

   This should point to an array of size :c:data:`grid_rank`.  This stores
   the starting value in each dimension for the field data.  This can be
   used to ignore boundary cells in grid data.

.. c:var:: int* grid_end

   This should point to an array of size :c:data:`grid_rank`.  This stores
   the end value in each dimension for the field data.  This can be used
   to ignore boundary cells in grid data.

.. c:var:: gr_float* grid_dx

   This is the grid cell width in c:data:`length_units`. This is currently
   used only in computing approximate H2 self-shielding when H2 is tracked
   (c:data:`primordial_chemistry` >= 2) and c:data:`H2_self_shielding` is 
   set to 1.

.. c:var:: gr_float* density

   Pointer to the density field array.

.. c:var:: gr_float* HI_density

   Pointer to the HI density field array.  Used when
   :c:data:`primordial_chemistry` is set to 1, 2, or 3.

.. c:var:: gr_float* HII_density

   Pointer to the HII density field array.  Used when
   :c:data:`primordial_chemistry` is set to 1, 2, or 3.

.. c:var:: gr_float* HM_density

   Pointer to the H\ :sup:`-`\  density field array.  Used when
   :c:data:`primordial_chemistry` is set to 2 or 3.

.. c:var:: gr_float* HeI_density

   Pointer to the HeI density field array.  Used when
   :c:data:`primordial_chemistry` is set to 1, 2, or 3.

.. c:var:: gr_float* HeII_density

   Pointer to the HeII density field array.  Used when
   :c:data:`primordial_chemistry` is set to 1, 2, or 3.

.. c:var:: gr_float* HeIII_density

   Pointer to the HeIII density field array.  Used when
   :c:data:`primordial_chemistry` is set to 1, 2, or 3.

.. c:var:: gr_float* H2I_density

   Pointer to the H\ :sub:`2`\  density field array.  Used when
   :c:data:`primordial_chemistry` is set to 2 or 3.

.. c:var:: gr_float* H2II_density

   Pointer to the H\ :sub:`2`\ \ :sup:`+`\  density field
   array.  Used when :c:data:`primordial_chemistry` is set to
   2 or 3.

.. c:var:: gr_float* DI_density

   Pointer to the DI density field array.  Used when
   :c:data:`primordial_chemistry` is set to 3.

.. c:var:: gr_float* DII_density

   Pointer to the DII density field array.  Used when
   :c:data:`primordial_chemistry` is set to 3.

.. c:var:: gr_float* HDI_density

   Pointer to the HD density field array.  Used when
   :c:data:`primordial_chemistry` is set to 3.

.. c:var:: gr_float* e_density

   Pointer to the electron density field array.  Used when
   :c:data:`primordial_chemistry` is set to 1, 2, or 3.  Note,
   the electron mass density should be scaled by the ratio of the
   proton mass to the electron mass such that the electron density
   in the code is the electron number density times the **proton**
   mass.

.. c:var:: gr_float* metal_density

   Pointer to the metal density field array.  Used when
   :c:data:`metal_cooling` is set to 1.

.. c:var:: gr_float* internal_energy

   Pointer to the internal energy field array.

.. c:var:: gr_float* x_velocity

   Pointer to the x-velocity field array.  Currently not used.

.. c:var:: gr_float* y_velocity

   Pointer to the y-velocity field array.  Currently not used.

.. c:var:: gr_float* z_velocity

   Pointer to the z-velocity field array.  Currently not used.

.. c:var:: gr_float* volumetric_heating_rate

   Pointer to values containing volumetric heating rates.  Rates
   should be in units of erg/s/cm\ :sup:`3`\.  Used when
   :c:data:`use_volumetric_heating_rate` is set to 1.

.. c:var:: gr_float* specific_heating_rate

   Pointer to values containing specific heating rates.  Rates
   should be in units of erg/s/g.  Used when
   :c:data:`use_specific_heating_rate` is set to 1.

.. c:var:: gr_float *RT_heating_rate

   Pointer to the radiation transfer heating rate field.  Rates
   should be in units of erg/s/cm\ :sup:`3`\.  Used when
   :c:data:`use_radiative_transfer` is set to 1.

.. c:var:: gr_float *RT_HI_ionization_rate

   Pointer to the HI photo-ionization rate field used with
   radiative transfer.  Rates should be in units of
   1/:c:data:`time_units`.  Used when
   :c:data:`use_radiative_transfer` is set to 1.

.. c:var:: gr_float *RT_HeI_ionization_rate

   Pointer to the HeI photo-ionization rate field used with
   radiative transfer.  Rates should be in units of
   1/:c:data:`time_units`.  Used when
   :c:data:`use_radiative_transfer` is set to 1.

.. c:var:: gr_float *RT_HeII_ionization_rate

   Pointer to the HeII photo-ionization rate field used with
   radiative transfer.  Rates should be in units of
   1/:c:data:`time_units`.  Used when
   :c:data:`use_radiative_transfer` is set to 1.

.. c:var:: gr_float *RT_H2_dissociation_rate

   Pointer to the H\ :sub:`2`\  photo-dissociation rate field
   used with radiative transfer.  Rates should be in units of
   1/:c:data:`time_units`.  Used when
   :c:data:`use_radiative_transfer` is set to 1 and
   :c:data:`primordial_chemistry` is either 2 or 3.

It is not necessary to attach a pointer to any field that you do
not intend to use.

.. code-block:: c++

  // Create struct for storing grackle field data
  grackle_field_data my_fields;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = new int[3];
  my_fields.grid_start = new int[3];
  my_fields.grid_end = new int[3];
  my_fields.grid_dx  = 1.0; // only matters if H2 self-shielding is used
  for (int i = 0;i < 3;i++) {
    my_fields.grid_dimension[i] = 1;
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;

  // Set field arrays.
  my_fields.density         = new gr_float[field_size];
  my_fields.internal_energy = new gr_float[field_size];
  my_fields.x_velocity      = new gr_float[field_size];
  my_fields.y_velocity      = new gr_float[field_size];
  my_fields.z_velocity      = new gr_float[field_size];
  // for primordial_chemistry >= 1
  my_fields.HI_density      = new gr_float[field_size];
  my_fields.HII_density     = new gr_float[field_size];
  my_fields.HeI_density     = new gr_float[field_size];
  my_fields.HeII_density    = new gr_float[field_size];
  my_fields.HeIII_density   = new gr_float[field_size];
  my_fields.e_density       = new gr_float[field_size];
  // for primordial_chemistry >= 2
  my_fields.HM_density      = new gr_float[field_size];
  my_fields.H2I_density     = new gr_float[field_size];
  my_fields.H2II_density    = new gr_float[field_size];
  // for primordial_chemistry >= 3
  my_fields.DI_density      = new gr_float[field_size];
  my_fields.DII_density     = new gr_float[field_size];
  my_fields.HDI_density     = new gr_float[field_size];
  // for metal_cooling = 1
  my_fields.metal_density   = new gr_float[field_size];
  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = new gr_float[field_size];
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = new gr_float[field_size];
  // heating rate from radiative transfer calculations (provide in units [erg s^-1 cm^-3]
  my_fields.RT_heating_rate = new gr_float[field_size];
  // HI ionization rate from radiative transfer calculations (provide in units of [ 1/time_units ]
  my_fields.RT_HI_ionization_rate = new gr_float[field_size];
  // HeI ionization rate from radiative transfer calculations (provide in units of [1/time_units]
  my_fields.RT_HeI_ionization_rate = new gr_float[field_size];
  // HeII ionization rate from radiative transfer calculations (provide in units of [1/time_units]
  my_fields.RT_HeII_ionization_rate = new gr_float[field_size];
  // H2 dissociation rate from radiative transfer calculations (provide in units of [1/time_units]
  my_fields.RT_H2_dissociation_rate = new gr_float[field_size];


.. note:: The electron mass density should be scaled by the ratio of the
   proton mass to the electron mass such that the electron density in the
   code is the electron number density times the **proton** mass.

.. _functions:

Calling the Available Functions
-------------------------------

There are five functions available, one to solve the chemistry and cooling 
and four others to calculate the cooling time, temperature, pressure, and the 
ratio of the specific heats (gamma).  The arguments required are the 
:c:data:`code_units` structure and the :c:data:`grackle_field_data` struct.
For the chemistry solving routine, a timestep must also be given.  For the
four field calculator routines, the array to be filled with the field values
must be created and passed as an argument as well.

Solve the Chemistry and Cooling
+++++++++++++++++++++++++++++++

.. code-block:: c++

  // some timestep (one million years)
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return 0;
  }

Calculating the Cooling Time
++++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *cooling_time;
  cooling_time = new gr_float[field_size];
  if (calculate_cooling_time(&my_units, &my_fields,
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return 0;
  }

Calculating the Temperature Field
+++++++++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *temperature;
  temperature = new gr_float[field_size];
  if (calculate_temperature(&my_units, &my_fields,
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }

Calculating the Pressure Field
++++++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *pressure;
  pressure = new gr_float[field_size];
  if (calculate_pressure(&my_units, &my_fields,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }

Calculating the Gamma Field
+++++++++++++++++++++++++++

.. code-block:: c++

  gr_float *gamma;
  gamma = new gr_float[field_size];
  if (calculate_gamma(&my_units, &my_fields,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }
