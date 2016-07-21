.. _reference:

API Reference
=============

Grackle has two versions of most functions.  The :ref:`primary_functions`,
discussed in :ref:`functions`, make use of the :c:data:`grackle_field_data`
struct and internally stored instances :c:data:`chemistry_data` and
:c:data:`chemistry_data_storage` structs.  A set of :ref:`internal_functions`
also exist and require the user to provide their own instances of the
:c:data:`chemistry_data` and :c:data:`chemistry_data_storage`.

.. _primary_functions:

Primary Functions
-----------------

.. c:function:: int set_default_chemistry_parameters(chemistry_data *my_grackle_data);

   Initializes the ``grackle_data`` data structure.  This must be called 
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

   Calculates the effective adiabatic index.  This is only useful with
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

.. _internal_functions:

Internal Functions
------------------

These functions are mostly for internal use, but can also be used to call the various
functions with different parameter values within a single code.

.. c:function:: chemistry_data _set_default_chemistry_parameters(void);

   Initializes and returns :c:type:`chemistry_data` data structure.  This must be
   called before run-time parameters can be set.

   :returns: data structure containing all run-time parameters and all chemistry and cooling data arrays
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

   Evolves the species densities and internal energies over a given timestep
   by solving the chemistry and cooling rate equations.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param double dt_value: the integration timestep in code units
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* x_velocity: array containing the x velocity values in code units
   :param gr_float* y_velocity: array containing the y velocity values in code units
   :param gr_float* z_velocity: array containing the z velocity values in code units
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array.  Used with :c:data:`primordial_chemistry` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass.  Used with :c:data:`primordial_chemistry` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array.  Used with :c:data:`metal_cooling` = 1.
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_cooling_time(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *cooling_time);

   Calculates the instantaneous cooling time.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* x_velocity, y_velocity, z_velocity: arrays containing the x, y, and z velocity values in code units
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* cooling_time: array which will be filled with the calculated cooling time values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_gamma(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *my_gamma);

   Calculates the effective adiabatic index.  This is only useful with 
   :c:data:`primordial_chemistry` >= 2 as the only thing that alters gamma from the single 
   value is H\ :sub:`2`.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* my_gamma: array which will be filled with the calculated gamma values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_pressure(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *pressure);

   Calculates the gas pressure.

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int _calculate_temperature(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *temperature);

   :param chemistry_data* my_chemistry: the structure returned by :c:func:`_set_default_chemistry_parameters`
   :param chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param code_units* my_units: code units conversions
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* temperature: array which will be filled with the calculated temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

   Calculates the gas temperature.
