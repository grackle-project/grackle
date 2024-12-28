
.. _integration:

Interacting with Grackle in Your Simulation Code
================================================

The majority of this document follows the implementation of Grackle in
a C++ simulation code.  Full implementation examples for
C, C++, and Fortran are also available in the Grackle source.  See
:ref:`examples` for more information.  For a list of all available
functions, see the :ref:`reference`.

Grackle currently supports two APIs.
The :ref:`Primary API <primary_functions>`, manages some of Grackle's data structures (e.g. :c:data:`chemistry_data` and :c:data:`chemistry_data_storage`) for a downstream simulation code by making use of global variables.
In contrast, the :ref:`Local API <local_functions>`, requires that the downstream application explicitly manage pointers to these same data-structures and requires that the pointers are provided as arguments to each function.
The latter API is explicitly thread-safe as it involves no global data.

.. _examples:

Example Executables
-------------------

The grackle source code contains examples for C, C++, and Fortran codes.  
They are located in the **src/example** directory and provide examples
of calling all of grackle's functions.

    * **c_example.c** - C example

    * **c_local_example.c** - C example using only :ref:`local_functions`

    * **cxx_example.C** - C++ example

    * **cxx_omp_example.C** - C++ example using OpenMP

    * **fortran_example.F** - Fortran example

The instructions for building and executing the examples vary based on the build-system:

.. tabs::

   .. tab:: Classic Build System

      Once you have already installed the grackle library, you can build the examples by typing *make* and the name of the file without extension.
      For example, to build the C++ example, type:

      .. code-block:: shell-session

         $ make cxx_example

      To run the example, make sure to add the path to the directory containing 
      the installed **libgrackle.so** to your LD_LIBRARY_PATH (or 
      DYLD_LIBRARY_PATH on Mac).


   .. tab:: CMake Build System

      By default, the examples are automatically built with the rest of Grackle.
      The compiled example binaries can be found within *<build-dir>/example*, where *<build-dir>* is the arbitrary build-directory that you need to specify when compiling Grackle.

      It's important that *<build-dir>* is a top-level directory in the grackle repository (e.g. something like *my-build* is fine, but choices like *../my-grackle-build* and *my_builds/my-first-build* are problematic).
      If this isn't the case, then the examples won't be able to locate the input data files.

      You don't need to worry about using LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH on Mac) to run these examples with this build-system.

.. important::

   The examples make certain assumptions about the location of the input files.
   To ensure that the input files can be found, you should execute each example-binary from the same directory where the example binary is produced.

Header Files
------------

Based on the language that your simulation code is written in, you will interact with one of the 2 header files that is shipped with Grackle:

    * **grackle.h** - the primary header file, containing declarations for all
      the available functions and data structures.  This is the only header
      file that needs to be included for C and C++ codes.

    * **grackle.def** - the header file to be used in Fortran codes.  Only
      this file needs to be included.

.. note::

   While other header files are also shipped with Grackle we consider them to be implementation details.
   With that said, we currently maintain the **grackle_types.h**, **grackle_chemistry_data.h**, **grackle_fortran_types.def**, **grackle_chemistry_data.h**, **grackle_fortran_types.def**, **grackle_fortran_interface.def**, and **grackle_macros.h** because they were previously mentioned in the documentation.
   **However, their usage is now deprecated.**
   We will only maintain these for a couple releases after 3.3.
   Afterwards, we may alter their contents or delete them entirely.

.. note::

   Earlier versions of Grackle required C++ codes to enclose the ``#include <grackle.h>`` include-directive within a C "language-linkage block" (the block starts with ``extern "C" {``).
   C++ codes should now directly include the header **(without the block)**.
   The headers internally use a standard idiom to properly handle this case.

Data Types
----------

The grackle library provides special configurable data types that are used to specify the floating-point precision of the baryon fields passed to the grackle functions.
For C and C++ codes, this is the :c:type:`gr_float` type.
For Fortran codes, this is the :c:type:`R_PREC` type.
Based on how you configure Grackle during compilation, both types will either alias the respective langauges' 4-byte (32-bit/single-precision) or 8-byte (64-bit/double precision) floating-point type.
When using the :ref:`classic build system <classic_build>`, this is controlled by assigning the ``CONFIG_PRECISION`` setting to *precision-32* or *precision-64*.
When using the :ref:`CMake build system <cmake_build>`, this is controlled by setting the ``GRACKLE_USE_DOUBLE`` cmake-variable to ``OFF`` or ``ON``.

.. c:type:: gr_float

   Floating point type used for the baryon fields in C and C++.
   When Grackle is compiled with *precision-32* (under the :ref:`classic build system <classic_build>`) or ``GRACKLE_USE_DOUBLE=OFF`` (under the :ref:`CMake build system <cmake_build>`), this aliases the :c:type:`!float` type.
   When Grackle is compiled with *precision-64* or ``GRACKLE_USE_DOUBLE=ON``, this aliases the :c:type:`!double` type.

.. c:type:: R_PREC

   The Fortran analog of :c:type:`gr_float`.
   When Grackle is compiled with *precision-32* (under the :ref:`classic build system <classic_build>`) or ``GRACKLE_USE_DOUBLE=OFF`` (under the :ref:`CMake build system <cmake_build>`), this aliases the *real\*4* type.
   When Grackle is compiled with *precision-64* or ``GRACKLE_USE_DOUBLE=ON``, this aliases the *real\*8* type.

Enabling Output
---------------

By default, grackle will not print anything but error messages.  However,
a short summary of the running configuration can be printed by setting the global
``grackle_verbose`` variable to 1.  In a parallel code, it is recommended that
output only be enabled for the root process.

.. code-block:: c++

   // Enable output
   grackle_verbose = 1;

.. _code-units:

Code Units
----------

Many of the calculations involved in chemical reactions and radiative
cooling include multiplications by density squared or even density
cubed. With typical gas densities relevant to galaxy formation being
of the order of one hydrogren atom per cubic centimeter (~10\
:sup:`-24` g/cm\ :sup:`3`, give or take a few orders of
magnitude), it is easy to end up with significant roundoff or
underflow errors when quantities are stored in CGS units.

The :c:data:`code_units` structure contains conversions from code
units to CGS such that a value passed to Grackle multiplied by the
appropriate code unit gives that value in CGS units. Units for
density, length, time, and the expansion factor must be set
manually. Units for velocity are then set by calling
:c:data:`set_velocity_units`. When using the proper frame (i.e.,
setting :c:data:`comoving_coordinates` to 0), :c:data:`a_units` (units
for the expansion factor) must be set to 1.0. See below for
recommendations on choosing appropriate units.

.. c:type:: code_units

   This structure contains the following members.

.. c:var:: int comoving_coordinates

   If set to 1, the incoming field data is assumed to be in the comoving
   frame. If set to 0, the incoming field data is assumed to be in the
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
   This should be set units the :c:data:`set_velocity_units` function. Note,
   units of specific energy (i.e., conversion to erg/g) are then defined
   as :c:data:`velocity_units`\ :sup:`2` (velocity units squared).

.. c:var:: double a_units

   Conversion factor to be multiplied by the expansion factor such that
   a\ :sub:`true`\  = a\ :sub:`code`\ * :c:data:`a_units`. When using
   proper coordinates, :c:data:`a_units` must be set to 1.

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
  my_units.a_units = 1.0;            // units for the expansion factor
  my_units.a_value = 1. / (1. + current_redshift) / my_units.a_units;
  // set velocity units
  set_velocity_units(&my_units);

Choosing Appropriate Units
^^^^^^^^^^^^^^^^^^^^^^^^^^

The main consideration when setting code units is to keep density,
length, and time values close to 1. Reasonable values for density,
length, and time units are the hydrogen mass in g, 1 kpc to 1 Mpc in
cm, and 1 Myr to 1 Gyr in s.

.. _comoving_coordinates:

Comoving Coordinates
^^^^^^^^^^^^^^^^^^^^

For cosmological simulations, a comoving unit system is preferred,
though not required, since it allows the densities to stay close to 1
as the universe expands. If :c:data:`comoving_coordinates` is set to
1, it is assumed that the fields being passed to the solver are in the
comoving frame. Hence, the units must convert from code units in the
**comoving** frame to CGS in the **proper** frame. If
:c:data:`comoving_coordinates` is set to 0, it is assumed that the
fields passed into the solver are in the proper frame. For an example
of using comoving units, see the `cosmological unit system
<https://github.com/enzo-project/enzo-dev/blob/main/src/enzo/CosmologyGetUnits.C>`__
in the `Enzo <http://enzo-project.org/>`_ code.

As the unit system is designed to convert from the comoving to the
proper frame, some of the values in the :c:data:`code_units` struct
are expected to change with expansion factor (or redshift) while some
others should remain constant. Units that should remain constant
include :c:data:`time_units` and :c:data:`a_units`. Units that should
vary are :c:data:`a_value` (obviously), :c:data:`length_units`, and
:c:data:`density_units`. Moving forward in time,
:c:data:`length_units` should be increasing proportional to
:c:data:`a_value` and :c:data:`density_units` should be decreasing as
:c:data:`a_value`:sup:`-3`.

There are two important corollaries of the above behavior. First, the
:c:data:`velocity_units` should remain constant. In comoving
coordinates, velocity units are given by

.. math::

   VU = \frac{LU}{a\ TU},

where VU is :c:data:`velocity_units`, a is :c:data:`a_value`, and TU
is :c:data:`time_units`. Second, the internal unit for the cooling
rate (equivalent to [erg s\ :sup:`-1` cm\ :sup:`+3`]) should remain
constant. The cooling unit (CU) is given by

.. math::

   CU = \frac{VU^2\ m_H^2}{DU\ a^3\ TU},

where DU is :c:data:`density_units` and m\ :sub:`H` is the hydrogen
mass. The above definitions also hold for proper coordinates by
setting a to 1.

As a convenience, we provide the :c:func:`gr_query_units` function to let the user query the code units at arbitrary values of the expansion factor.
This is intended to be used when implementing Grackle support into your simulation (and debugging) in order to let you manually confirm that the code units you are passing in are correct (internally, Grackle simply assumes that the caller is giving the right values).
The values returned by this function are based on the the initial values stored in the :c:type:`code_units` data structure (that are passed into grackle while :ref:`setting up chemistry data and the associated storage <setup_data-storage>`).


.. _setup_data-storage:

Chemistry Data
--------------

Grackle's behavior is controlled by the :c:type:`chemistry_data` type.
It is a structure that stores all of the relevant run-time parameters.
After an instance of the type is configured, it is used to initialize all of the chemistry and cooling rate arrays (stored by the separate :c:type:`chemistry_data_storage` type).


.. c:type:: chemistry_data

   This structure holds all grackle run-time parameters, which are listed in
   :ref:`parameters`.

.. c:type:: chemistry_data_storage

   This structure holds all chemistry and cooling rate arrays. The user will
   not normally need to work directly with its internals. Users calling the 
   :ref:`primary_functions` will not encounter it (because the Primary 
   Functions make use of an internally stored instance of this type).
   Users implementing the :ref:`local_functions` will have to store a pointer
   to one of these.

There is a 2 step procedure for setting up Grackle's chemistry data (the precise details vary based on the choosen API flavor):


1. Configure the solver's run-time parameters.
   There are 3 parts to this step.
   First, allocate memory for a :c:type:`chemistry_data` instance.
   Next, initialize the instance by calling a function to perform the initial setup of the instance and store all the default parameter values.
   Finally, assign the desired values to the parameters tracked by the instance.
   See :ref:`parameters` for a full list of the available parameters.


.. tabs::

   .. group-tab:: Primary Functions

      In this case, the :c:func:`set_default_chemistry_parameters` routine for initializing the structure holding the parameters. 
      The function must be handed a pointer to an instance of :c:type:`chemistry_data`.
      That pointer is attached to the global :c:var:`!grackle_data` variable and the memory is properly initialized.
      The function returns an integer indicating success (1) or failure (0). 

      The convention is to assign desired parameter values by accessing :c:var:`!grackle_data`.

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

   .. group-tab:: Local Functions

      In this scenario, the :c:func:`local_initialize_chemistry_parameters` routine is used to initialize the previously allocated memory of a :c:type:`chemistry_data` instance.
      The function returns an integer indicating success (1) or failure (0).

      Unlike the other scenario, no global variable is initialized; instead the application is responsible for keeping track of the allocated object. 

      .. code-block:: c++

         chemistry_data *my_grackle_data;
         my_grackle_data = new chemistry_data;
         if (local_initialize_chemistry_parameters(my_grackle_data) == 0) {
           fprintf(stderr, "Error in local_initialize_chemistry_parameters.\n");
         }

         // Set parameter values for chemistry.
         my_grackle_data->use_grackle = 1;            // chemistry on
         my_grackle_data->with_radiative_cooling = 1; // cooling on
         my_grackle_data->primordial_chemistry = 3;   // molecular network with H, He, D
         my_grackle_data->metal_cooling = 1;          // metal cooling on
         my_grackle_data->UVbackground = 1;           // UV background on
         my_grackle_data->grackle_data_file = "CloudyData_UVB=HM2012.h5"; // data file


2. Once the desired parameters have been set, the chemistry and cooling rates must be initialized by calling with a
pointer to the :c:data:`code_units` struct created earlier.  This function
will return an integer indicating success (1) or failure (0).

.. tabs::

   .. group-tab:: Primary Functions
      
      In this case, the initialization function is called :c:func:`initialize_chemistry_data`.
      This function internally allocates memory for the type that holds the storage data, :c:type:`chemistry_data_storage`, and stores the address in a variable that is tracked behind the scenes.
      This storage data is automatically passed to other Primary Functions.

      .. code-block:: c++
      
         // Finally, initialize the chemistry object.
         if (initialize_chemistry_data(&my_units) == 0) {
           fprintf(stderr, "Error in initialize_chemistry_data.\n");
           return 0;
         }


   .. group-tab:: Local Functions

      In this case, the initialization function is called :c:func:`local_initialize_chemistry_data`.
      The application is responsible for passing the pointer to the (previously configured) :c:func:`chemistry_data` instance into the function.
      The application is also responsible for allocating the memory for the :c:type:`chemistry_data_storage` instance that the function initializes.
      **Note:** the initialization function allocates additional memory that is tracked within the :c:type:`chemistry_data_storage` for holding chemistry and cooling rates.

      Unlike the other scenario, the downstream application is responsible for tracking the initialized object (so that it can be passed to other functions).

      .. code-block:: c++

         // Allocate the chemistry_data_storage type
         chemistry_data_storage *my_grackle_rates = new chemistry_data_storage;

         // Finally, initialize the chemistry object.
         if (local_initialize_chemistry_data(my_chemistry_data, my_grackle_rates,
                                             &my_units) == 0) {
           fprintf(stderr, "Error in local_initialize_chemistry_data.\n");
           return 0;
         }

The Grackle is now ready to be used.

As an aside, see :ref:`dynamic-api` for a description of an alternative approach for configuring a :c:type:`chemistry_data` struct. This other approach may provide additional compatability with multiple versions of Grackle, and in some cases may facillitate less-verbose, easier-to-maintain code.

.. _openmp:

Running with OpenMP
-------------------

As of version 2.2, Grackle can be run with OpenMP parallelism.
To do this, the library must be compiled with OpenMP support.
When compiling Grackle

  * with the :ref:`classic build system <classic_build>`, you should execute ``make omp-on`` before compiling (more information about modifying settings in this system are provided :ref:`here <compiler-settings>`)

  * with the :ref:`CMake build system <cmake_build>`, you should configure the build by assigning the ``GRACKLE_USE_OPENMP`` cmake-variable a value of ``ON`` (more information about modifying settings in this system are provided :ref:`here <how_to_configure>`)

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

.. _density-note:

Note on Density Fields
^^^^^^^^^^^^^^^^^^^^^^

All density fields provided to Grackle should be mass densities, i.e.,
the number density of a given species multiplied by its mass. The
units should be such that the field value multiplied by
:c:data:`density_units` results in a value with units of g/cm\
:sup:`3`. See :ref:`code-units` for further discussion of Grackle unit
systems.

.. _e-density-note:

Note on the Electron Density Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For largely historical reasons, the value of the electron density
(provided in :c:data:`e_density`) should be the true value of the
electron mass density multiplied by the ratio of the proton mass to
the electron mass, i.e., :c:data:`e_density` = :math:`{\rho}`\ :sub:`e` * m\
:sub:`p` / m\ :sub:`e`, where :math:`{\rho}`\ :sub:`e` is the true
electron mass density in :c:data:`density_units` (see :ref:`density-note`).

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

.. c:var:: gr_float grid_dx

   This is the grid cell width in :c:data:`length_units`. This is currently
   used only in computing approximate H2 self-shielding when H2 is tracked
   (:c:data:`primordial_chemistry` >= 2) and :c:data:`H2_self_shielding` is
   set to 1. If this can't be assigned a meaningful value (e.g. the field data
   does not organized on a grid or the grid cells aren't perfect cubes), we
   recommend assigning it a value of -1 (so that error-handling works properly)

.. c:var:: gr_float* density

   Pointer to the gas density field array.

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
   proton mass to the electron mass. See :ref:`e-density-note` for
   more information.

.. c:var:: gr_float* metal_density

   Pointer to the metal density field array.  Used when
   :c:data:`metal_cooling` is set to 1.

.. c:var:: gr_float* dust_density

   Pointer to the dust density field array.  Used when
   :c:data:`use_dust_density_field` is set to 1.

.. c:var:: gr_float* internal_energy

   Pointer to the internal energy field array. Internal energies should be
   in units of :c:data:`velocity_units`\ :sup:`2` (velocity units squared).
   This can be converted to and from a temperature by using the
   :c:data:`get_temperature_units` function.

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

.. c:var:: gr_float* temperature_floor

   Pointer to values containing a temperature floor for each element
   in units of K. No chemistry or cooling calculations will be
   performed on an element with a temperature at or below the
   specified value. Used when :c:data:`use_temperature_floor` is
   set to 2.

.. c:var:: gr_float *RT_heating_rate

   Pointer to the radiation transfer heating rate field.  Rates
   should be in units of (erg/s/cm\ :sup:`3`\) / n\ :sub:`HI`\, where
   n\ :sub:`HI`\  is the neutral hydrogen number density. Heating rates
   for additional species are currently not yet supported.
   Used when :c:data:`use_radiative_transfer` is set to 1.

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

.. c:var:: gr_float *H2_self_shielding_length

   Pointer to a field containing lengths to be used for
   calculating molecular hydrogen column denisty for
   H2\ :sub:`2`\ self-shielding.  Used when
   :c:data:`H2_self_shielding` is set to 2.  Field data
   should be in :c:data:`length_units`.

.. c:var:: gr_float *H2_custom_shielding_factor

   Pointer to a field containing attenuation factors to 
   be multiplied with the H\ :sub:`2`\ dissociation rate.
   Used when the :c:data:`H2_custom_shielding` flag is set.

.. c:var:: gr_float *isrf_habing

   Pointer to a field containing values of the strength
   of the insterstellar radiation field used in the
   calculation of dust heating. This is used when
   :c:data:`use_isrf_field` is set to 1. The units
   of this field should be the same as those of the
   :c:data:`interstellar_radiation_field` parameter.

It is not necessary to attach a pointer to any field that you do
not intend to use.

.. code-block:: c++

  // Create struct for storing grackle field data
  grackle_field_data my_fields;

  // initialize members of my_fields to sensible defaults
  gr_initialize_field_data(&my_fields);

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = new int[3];
  my_fields.grid_start = new int[3];
  my_fields.grid_end = new int[3];
  my_fields.grid_dx  = 1.0; // only matters if H2 self-shielding is used
                            // we recommend assigning it a value of -1 if your
                            // simulation doesn't have a meaningful value for
                            // it (so that error-handling works properly)
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

There are six functions available, one to solve the chemistry and cooling
and five others to calculate the cooling time, temperature, pressure,
ratio of the specific heats (gamma), and dust temperature. The
arguments required are the :c:data:`code_units` structure and the
:c:data:`grackle_field_data` struct. For the chemistry solving
routine, a timestep must also be given. For the four field calculator
routines, the array to be filled with the field values must be created
and passed as an argument as well.

The examples below make use of Grackle's :ref:`primary_functions`, where
the parameters and rate data are stored in instances of the
:c:data:`chemistry_data` and :c:data:`chemistry_data_storage` structs
declared in **grackle.h**.  Alternatively, a set of :ref:`local_functions`
require these structs to be provided as arguments, allowing for explicitly
thread-safe code.

Solve the Chemistry and Cooling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

  // some timestep (one million years)
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return 0;
  }

Calculating the Cooling Time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

  gr_float *cooling_time;
  cooling_time = new gr_float[field_size];
  if (calculate_cooling_time(&my_units, &my_fields,
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return 0;
  }

Calculating the Temperature Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

  gr_float *temperature;
  temperature = new gr_float[field_size];
  if (calculate_temperature(&my_units, &my_fields,
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }

Calculating the Pressure Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

  gr_float *pressure;
  pressure = new gr_float[field_size];
  if (calculate_pressure(&my_units, &my_fields,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }

Calculating the Gamma Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

  gr_float *gamma;
  gamma = new gr_float[field_size];
  if (calculate_gamma(&my_units, &my_fields,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }

Calculating the Dust Temperature Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

  gr_float *dust_temperature;
  dust_temperature = new gr_float[field_size];
  if (calculate_dust_temperature(&my_units, &my_fields,
                            dust_temperature) == 0) {
    fprintf(stderr, "Error in calculate_dust_temperature.\n");
    return EXIT_FAILURE;
  }

Clearing the memory
-------------------


When using the Grackle's :ref:`primary_functions`, global structures are used and therefore the global structure ``grackle_rates`` needs to be released with

.. code-block:: c++

  free_chemistry_data();

When using the :ref:`local_functions`, you should call the counter-part function, :c:func:`local_free_chemistry_data`, to free memory used for storing chemistry and cooling rates.

.. _query-version:

Querying library version information
------------------------------------

A struct of type :c:data:`grackle_version` is used to hold version
information about the version of Grackle that is being used. The
struct contains information about the version number and particular
git revision.

.. c:type:: grackle_version

   This structure is used to organize version information for the
   library.

.. c:var:: const char* version

   Specifies the version of the library using this template:
   ``<MAJOR>.<MINOR>(.<MICRO>)(.dev<DEV_NUM>)``. In this template
   ``<MAJOR>``, ``<MINOR>``, and ``<MICRO>`` correspond to a major,
   minor, and micro version numbers (the micro version number is
   omitted if it's zero). The final section can specify a
   development version. More details about this version number can be
   found :ref:`here <primary-version-number>`.

.. c:var:: const char* branch

   Specifies the name of the git branch that the library was compiled
   from.

.. c:var:: const char* revision

   Specifies the hash identifying the git commit that the library was
   compiled from.

The :c:func:`get_grackle_version` function is used to retrieve a
properly intialized :c:data:`grackle_version` object. The following
code snippet illustrates how one might query and print this
information:

.. code-block:: c++

  grackle_version gversion = get_grackle_version();
  printf ("The Grackle Version: %s\n", gversion.version);
  printf ("Git Branch:   %s\n", gversion.branch);
  printf ("Git Revision: %s\n", gversion.revision);

.. _dynamic-api:

Dynamic configuration of Chemistry Data
---------------------------------------

The functions providing dynamic access to the fields of :c:data:`chemistry_data` are useful for maintaining backwards compatibility with older versions of Grackle (that also provide this API) as new fields get added to :c:data:`chemistry_data`. This is exemplified in the following scenario.

Suppose Grackle is updated to have a new heating/cooling mechanism, and to allow users to control that mechanism two new fields are added to :c:data:`chemistry_data`:
   * an ``int`` field called ``use_fancy_feature``
   * a ``double`` field called ``fancy_feature_param``

Now suppose a downstream simulation code, written in ``c`` or ``c++``, wanted to support configuration of this feature. In this scenario, imagine that we have a pointer to a :c:data:`chemistry_data` structure called ``my_grackle_data``.

The obvious way to configure this feature is to include the following snippet in the simulation code:

.. code-block:: c++

  if (configure_fancy_feature) {
    my_grackle_data->use_fancy_feature = 1;
    my_grackle_data->fancy_feature_param = 5.0; // arbitrary value
  }

However, inclusion of the above snippet will prevent the simulation code from compiling if the user has a version of Grackle installed in which :c:data:`chemistry_data` does not have the ``use_fancy_feature`` and ``fancy_feature_param`` fields. Consequently, such users will have to update Grackle.

  * This can be inconvenient when a user has no interest in using this new feature, but needs an unrelated feature/bugfix introduced to the code in a subsequent changeset

  * This is especially inconvenient if a user is prototying a new feature in a custom Grackle branch in which the :c:data:`chemistry_data` struct is missing these fields.

The following snippet shows how the dynamic access API can be used in the same way for versions of Grackle that include these parameters, and don't set the features in cases 

.. code-block:: c++

  if (configure_fancy_feature) {
    int* use_fancy_feature = local_chemistry_data_access_int(
      my_grackle_data, "use_fancy_feature"
    );
    double* fancy_feature_param = local_chemistry_data_access_double(
      my_grackle_data, "fancy_feature_param"
    );

    if ((use_fancy_feature == NULL) || (fancy_feature_param == NULL)){
      fprintf(stderr, "Update grackle version to use fancy feature\n");
    } else {
      *use_fancy_feature = 1;
      *fancy_feature_param = 5.0;
    }
  }

There are a few points worth noting:

  * As the above snippets show, the dynamic api clearly produces more verbose code when configuring :c:data:`chemistry_data` field-by-field. However, in codes where users configure Grackle by specifying the name of fields in the :c:data:`chemistry_data` struct and the associated values in a parameter file, the dynamic API can facillitate MUCH less verbose code. Under certain implementations, it may not even be necessary to modify a simulation code to support newly-introduced grackle parameters.

  * The dynamic API is slower than configuring :c:data:`chemistry_data` in the classic approach. However, this shouldn't be an issue since :c:data:`chemistry_data` is usually just configured once when the simulation code starts up.

  * The highlighted functions can also be used in tandem with other functions described in :ref:`dynamic_api_functions` to simplify (de)serialization of :c:data:`chemistry_data`.

  * For completeness, the dynamic API also provides an analogous function for configuring string parameters.
