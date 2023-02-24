.. _reference:

API Reference
=============

Grackle has three versions of most functions.

1. The :ref:`primary_functions`, discussed in :ref:`functions`, make
   use of internally stored instances of the :c:data:`chemistry_data`
   and :c:data:`chemistry_data_storage` structs declared in **grackle.h**.

2. :ref:`local_functions` require pointers to :c:data:`chemistry_data`
   and :c:data:`chemistry_data_storage` instances to be provided as
   arguments. These are explicity thread-safe as they use no global data.

3. (Deprecated) :ref:`internal_functions` take pointers to individual
   field arrays instead of using the :c:data:`grackle_field_data`
   struct. These functions have been deprecated and will be removed in
   versions of Grackle later than 3.2.

.. _primary_functions:

Primary Functions
-----------------

.. c:function:: int set_default_chemistry_parameters(chemistry_data *my_grackle_data);

   Initializes the ``grackle_data`` data structure. This must be called
   before run-time parameters can be set.

   :param chemistry_data* my_grackle_data: run-time parameters
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int initialize_chemistry_data(code_units *my_units);

   Loads all chemistry and cooling data, given the set run-time parameters.
   This can only be called after :c:func:`set_default_chemistry_parameters`.

   :param code_units* my_units: code units conversions
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: void set_velocity_units(code_units *my_units);

   Sets the :c:data:`velocity_units` value of the input ``my_units``
   :c:data:`code_units` struct. For proper coordinates, velocity units are equal
   to :c:data:`length_units` / :c:data:`time_units`. For comoving coordinates,
   velocity units are equal to (:c:data:`length_units` / :c:data:`a_value`) /
   :c:data:`time_units`.

   :param code_units* my_units: code units conversions

.. c:function:: double get_velocity_units(code_units *my_units);

   Returns the appropriate value for velocity units given the values of
   :c:data:`length_units`, :c:data:`a_value`, and :c:data:`time_units`
   in the input ``my_units`` :c:data:`code_units` struct. For proper coordinates,
   velocity units are equal to :c:data:`length_units` / :c:data:`time_units`.
   For comoving coordinates, velocity units are equal to (:c:data:`length_units`
   / :c:data:`a_value`) / :c:data:`time_units`. Note, this function only returns
   a value, but does not set it in the struct. To set the value in the struct, use
   :c:data:`set_velocity_units`.

   :param code_units* my_units: code units conversions
   :rtype: double
   :returns: velocity_units

.. c:function:: double get_temperature_units(code_units *my_units);

   Returns the conversion factor between specific internal energy and temperature
   assuming gamma (the adiabatic index) = 1, such that temperature in K is equal to
   :c:data:`internal_energy` * ``temperature_units``. This unit conversion is
   defined as m\ :sub:`H` * :c:data:`velocity_units`\ :sup:`2` / k\ :sub:`b`,
   where m\ :sub:`H` is the Hydrogen mass and k\ :sub:`b` is the Boltzmann constant.

   :param code_units* my_units: code units conversions
   :rtype: double
   :returns: temperature_units

.. c:function:: int solve_chemistry(code_units *my_units, grackle_field_data *my_fields, double dt_value);

   Evolves the species densities and internal energies over a given timestep 
   by solving the chemistry and cooling rate equations.

   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param double dt_value: the integration timestep in code units
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_cooling_time(code_units *my_units, grackle_field_data *my_fields, gr_float *cooling_time);

   Calculates the instantaneous cooling time.

   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* cooling_time: array which will be filled with the calculated cooling time values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_gamma(code_units *my_units, grackle_field_data *my_fields, gr_float *my_gamma);

   Calculates the effective adiabatic index. This is only useful with
   :c:data:`primordial_chemistry` >= 2 as the only thing that alters gamma from the single
   value is H\ :sub:`2`.

   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* my_gamma: array which will be filled with the calculated gamma values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_pressure(code_units *my_units, grackle_field_data *my_fields, gr_float *pressure);

   Calculates the gas pressure.

   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_temperature(code_units *my_units, grackle_field_data *my_fields, gr_float *temperature);

   Calculates the gas temperature.

   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* temperature: array which will be filled with the calculated temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_dust_temperature(code_units *my_units, grackle_field_data *my_fields, gr_float *dust_temperature);

   Calculates the dust temperature. The dust temperature calculation is
   modified from its original version (Section 4.3 of `Smith et al. 2017
   <http://ui.adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`__) to also
   include the heating of dust grains by the interstellar radiation field
   following equation B15 of `Krumholz (2014)
   <https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.1662K/abstract>`__.

   Using this function requires :c:data:`dust_chemistry` > 0 or :c:data:`h2_on_dust` > 0.

   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* dust_temperature: array which will be filled with the calculated dust temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: grackle_version get_grackle_version();

   Constructs and returns a :c:type:`grackle_version` struct that
   encodes the version information for the library.

   :rtype: `grackle_version`

.. _local_functions:

Local Functions
---------------

These can be used to create explicitly thread-safe code or to call
the various functions with different parameter values within a
single code.

.. _local_setup_data-storage:

Initializing/Configuring Chemistry Parameters and Storage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

..
   COMMENT BLOCK
   This section helps simplify the description of the local function
   (it provides a place for us to point people to when we talk about
   expectations for how the arguments are initialized). Furthermore,
   there currently isn't any other place where we currently provide
   this sort of overview of using the local functions.

   In reality this content probably fits better on the page about
   "Adding Grackle to Your Simulation Code". However, this duplicates
   a bunch of content related to using the global Primary functions.
   Maybe we could use the sphinx_tabs extension...


The approach for initializing/configuring the :c:data:`chemistry_data` and :c:data:`chemistry_data_storage` structs for use with the :ref:`local_functions` differs to some degree from the :ref:`previously described approach <setup_data-storage>` that is used with :ref:`primary_functions`.

We highlight the steps down below.
For the sake of argument let's imagine that the user is storing the chemistry data in a variable called ``my_chemistry``.

1. First, the user should allocate ``my_chemistry`` and initialize the stored parameters with :c:func:`local_initialize_chemistry_parameters`.

   .. code-block:: c++

      chemistry_data * my_chemistry_data = new chemistry_data;
      if (local_initialize_chemistry_parameters(my_chemistry_data) == 0) {
         fprintf(stderr, "Error in local_initialize_chemistry_parameters.\n");
      }

2. Next, the user can configure the stored parameters.
   They can do this by directly modifying the stored parameters (e.g. ``my_chemistry->use_grackle = 1``) or by using the :ref:`dynamic-api`.

After the user has finished initializing ``my_chemistry``, and have configured an instance of :c:data:`code_units` (more detail provided :ref:`here <code-units>`), they can initialize an instance of :c:data:`chemistry_data_storage` with 
:c:func:`local_initialize_chemistry_data`:

.. code-block:: c++

  chemistry_data_storage* my_rates = new chemistry_data_storage;
  if (local_initialize_chemistry_data(my_chemistry, my_rates, &my_units) == 0) {
    fprintf(stderr, "Error in local_initialize_chemistry_data.\n");
    return 0;
  }

Configuration/Cleanup Functions
+++++++++++++++++++++++++++++++

.. c:function:: int local_initialize_chemistry_parameters \
                (chemistry_data *my_chemistry);

   Initializes the parameters stored in the :c:type:`chemistry_data` data structure to their default values.
   This should be called before run-time parameters are set.

   This is the "local" counterpart to :c:func:`set_default_chemistry_parameters`.

   :param chemistry_data \*my_chemistry: run-time parameters
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int local_initialize_chemistry_data \
                (chemistry_data *my_chemistry, \
                chemistry_data_storage *my_rates, \
                code_units *my_units);

   Allocates storage for and initializes the values of all relevant chemistry and cooling rate data.
   This data is stored within the provided :c:data:`chemistry_data_storage` struct.
   This is the "local" counterpart to :c:func:`initialize_chemistry_data`.

   This function should only be called after the user has finished configuring both ``my_chemistry`` and ``my_units``.
   This function assumes that none of ``my_rates``'s members of pointer type hold valid memory addresses (i.e. where applicable, the function allocates fresh storage and makes no attempts to deallocate/reuse storage).

   After calling this function, the user should avoid modifying any of the fields of ``my_chemistry`` or ``my_units``.
   The **only** exception to this rule is that ``my_units``'s :c:var:`a_value` field may be mutated when using comoving units (i.e. ``my_units``'s :c:var:`comoving_coordinates` field is ``1``).

   To deallocate any storage allocated by this function, use :c:func:`free_chemistry_data`.

   :param chemistry_data \*my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage \*my_rates: chemistry and cooling rate data structure
   :param code_units \*my_units: code units conversions
   :rtype: int
   :returns: 1 (success) or 0 (failure)

   .. note::
      In addition to modifying the contents of ``my_rates``, this function may also mutate the values stored in ``my_chemistry`` to set them to "more sensible" values (based on other values stored in ``my_chemistry``).


.. c:function:: int free_chemistry_data \
                (chemistry_data *my_chemistry, \
                chemistry_data_storage *my_rates);

   Deallocates all data held by the members of ``my_rates`` allocated during its initialization in :c:func:`local_initialize_chemistry_data` (or :c:func:`initialize_chemistry_data`).
   Issues may arise if ``my_chemistry`` was mutated between the initialization of ``my_rates`` and the call to this function.

   :param chemistry_data \*my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage \*my_rates: previously initialized chemistry and cooling rate data structure
   :rtype: int
   :returns: 1 (success) or 0 (failure)
   

Chemistry Functions
+++++++++++++++++++

.. c:function:: int local_solve_chemistry(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, grackle_field_data *my_fields, double dt_value);

   Evolves the species densities and internal energies over a given timestep
   by solving the chemistry and cooling rate equations.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param double dt_value: the integration timestep in code units
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int local_calculate_cooling_time(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, grackle_field_data *my_fields, gr_float *cooling_time);

   Calculates the instantaneous cooling time.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* cooling_time: array which will be filled with the calculated cooling time values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int local_calculate_gamma(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, grackle_field_data *my_fields, gr_float *my_gamma);

   Calculates the effective adiabatic index. This is only useful with
   :c:data:`primordial_chemistry` >= 2 as the only thing that alters gamma from the single
   value is H\ :sub:`2`.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* my_gamma: array which will be filled with the calculated gamma values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int local_calculate_pressure(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, grackle_field_data *my_fields, gr_float *pressure);

   Calculates the gas pressure.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int local_calculate_temperature(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, grackle_field_data *my_fields, gr_float *temperature);

   Calculates the gas temperature.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* temperature: array which will be filled with the calculated temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int local_calculate_dust_temperature(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, grackle_field_data *my_fields, gr_float *dust_temperature);

   Calculates the dust temperature.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param grackle_field_data* my_fields: field data storage
   :param gr_float* dust_temperature: array which will be filled with the calculated dust temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. _dynamic_api_functions:

Dynamic Configuration Functions
+++++++++++++++++++++++++++++++

The following functions are used to provide dynamic access to members of the :c:data:`chemistry_data` struct. They will return ``NULL`` when ``my_chemistry`` is ``NULL``, ``param_name`` isn't a known parameter, or the ``param_name`` is not associated with the type mentioned in the function name.

.. c:function:: int* local_chemistry_data_access_int(chemistry_data *my_chemistry, const char* param_name);

   Returns the pointer to the member of ``my_chemistry`` associated with ``param_name``.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param const char* param_name: the name of the parameter to access.
   :rtype: int*

.. c:function:: double* local_chemistry_data_access_double(chemistry_data *my_chemistry, const char* param_name);

   Returns the pointer to the member of ``my_chemistry`` associated with ``param_name``.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param const char* param_name: the name of the parameter to access.
   :rtype: double*

.. c:function:: char** local_chemistry_data_access_string(chemistry_data *my_chemistry, const char* param_name);

   Returns the pointer to the member of ``my_chemistry`` associated with ``param_name``.

   :param chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param const char* param_name: the name of the parameter to access.
   :rtype: char**

The following functions are used to query the name of the ith field of the :c:data:`chemistry_data` struct of a particular type.

.. c:function:: const char* param_name_int(size_t i);

   Query the name of the ith ``int`` field from :c:data:`chemistry_data`.

   .. warning:: The order of parameters may change between different versions of Grackle.

   :param size_t i: The index of the accessed parameter
   :rtype: const char*
   :returns: Pointer to the string-literal specifying the name. This is ``NULL``, if :c:data:`chemistry_data` has ``i`` or fewer ``int`` members
   
.. c:function:: const char* param_name_double(size_t i);

   Query the name of the ith ``double`` field from :c:data:`chemistry_data`.

   .. warning:: The order of parameters may change between different versions of Grackle.

   :param size_t i: The index of the accessed parameter
   :rtype: const char*
   :returns: Pointer to the string-literal specifying the name. This is ``NULL``, if :c:data:`chemistry_data` has ``i`` or fewer ``double`` members.

.. c:function:: const char* param_name_string(size_t i);

   Query the name of the ith ``string`` field from :c:data:`chemistry_data`.

   .. warning:: The order of parameters may change between different versions of Grackle.

   :param size_t i: The index of the accessed parameter
   :rtype: const char*
   :returns: Pointer to the string-literal specifying the name. This is ``NULL``, if :c:data:`chemistry_data` has ``i`` or fewer ``string`` members.

.. _internal_functions:

Internal Functions
------------------

These functions are mostly for internal use. The initialization functions
described here can be used in conjunction with the :ref:`local_functions`.

.. c:function:: chemistry_data _set_default_chemistry_parameters(void);

   This function has been deprecated and will be removed in versions
   of Grackle later than 3.2.
   Please use :c:func:`set_default_chemistry_parameters` or :c:func:`local_initialize_chemistry_parameters`.

   Initializes and returns :c:type:`chemistry_data` data structure. This must be
   called before run-time parameters can be set.

   :returns: data structure containing all run-time parameters
   :rtype: :c:type:`chemistry_data`

.. c:function:: int _initialize_chemistry_data(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units);

   Creates all chemistry and cooling rate data and stores within the provided :c:data:`chemistry_data_storage` struct.
   This can only be called after :c:func:`_set_default_chemistry_parameters`.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _solve_chemistry(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, double dt_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density);

   This function has been deprecated and will be removed in versions
   of Grackle later than 3.2. Please use solve_chemistry or local_solve_chemistry.

   Evolves the species densities and internal energies over a given timestep
   by solving the chemistry and cooling rate equations.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param double dt_value: the integration timestep in code units
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* x_velocity: array containing the x velocity values in code units
   :param gr_float* y_velocity: array containing the y velocity values in code units
   :param gr_float* z_velocity: array containing the z velocity values in code units
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array. Used with :c:data:`metal_cooling` = 1.
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_cooling_time(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *cooling_time);

   This function has been deprecated and will be removed in versions
   of Grackle later than 3.2. Please use calculate_cooling_time or local_calculate_cooling_time.

   Calculates the instantaneous cooling time.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* x_velocity: array containing the x velocity values in code units
   :param gr_float* y_velocity: array containing the y velocity values in code units
   :param gr_float* z_velocity: array containing the z velocity values in code units
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array. Used with :c:data:`metal_cooling` = 1.
   :param gr_float* cooling_time: array which will be filled with the calculated cooling time values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_gamma(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *my_gamma);

   This function has been deprecated and will be removed in versions
   of Grackle later than 3.2. Please use calculate_gamma or local_calculate_gamma.

   Calculates the effective adiabatic index. This is only useful with
   :c:data:`primordial_chemistry` >= 2 as the only thing that alters gamma from the single 
   value is H\ :sub:`2`.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array. Used with :c:data:`metal_cooling` = 1.
   :param gr_float* my_gamma: array which will be filled with the calculated gamma values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_pressure(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *pressure);

   This function has been deprecated and will be removed in versions
   of Grackle later than 3.2. Please use calculate_pressure or local_calculate_pressure.

   Calculates the gas pressure.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array. Used with :c:data:`metal_cooling` = 1.
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_temperature(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *temperature);

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields. This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array. Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass. Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array. Used with :c:data:`metal_cooling` = 1.
   :param gr_float* temperature: array which will be filled with the calculated temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

   This function has been deprecated and will be removed in versions
   of Grackle later than 3.2. Please use calculate_temperature or local_calculate_temperature.

   Calculates the gas temperature.
