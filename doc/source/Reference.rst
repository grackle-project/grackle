.. _reference:

API Reference
=============

The Grackle has a few different versions of the various functions for solving the 
chemistry and cooling and calculating related fields.  One set of functions 
requires the user to work with C structs, while the other does not.  The set that 
does not use structs is simpler to implement in Fortran codes.  Both of these rely 
internally on a *chemistry_data* type struct called ``grackle_data``, which exists 
in the grackle namespace.  A third set of functions also exists that requires the 
user to hold and pass their own *chemistry_data* struct.

Functions using structs (best for C and C++)
--------------------------------------------

These functions require the user to directly access the ``grackle_data`` 
data structure to set parameters and to creata a ``code_units`` struct 
to control the unit system.

.. c:function:: int set_default_chemistry_parameters();

   Initializes the ``grackle_data`` data structure.  This must be called 
   before run time parameters can be set.

   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int initialize_chemistry_data(code_units *my_units, double a_value);

   Loads all chemistry and cooling data, given the set run time parameters.  
   This can only be called after :c:func:`set_default_chemistry_parameters`.  

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int solve_chemistry(code_units *my_units, double a_value, double dt_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density);

   Evolves the species densities and internal energies over a given timestep 
   by solving the chemistry and cooling rate equations.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
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
   :param gr_float* HI_density: array containing the HI densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 1.
   :param gr_float* HII_density: array containing the HII densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 1.
   :param gr_float* HM_density: array containing the H\ :sup:`-`\  densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 2.
   :param gr_float* HeI_density: array containing the HeI densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 1.
   :param gr_float* HeII_density: array containing the HeII densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 1.
   :param gr_float* HeIII_density: array containing the HeIII densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 1.
   :param gr_float* H2I_density: array containing the H\ :sub:`2`:\  densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 2.
   :param gr_float* H2II_density: array containing the H\ :sub:`2`:sup:`+`\ densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` >= 2.
   :param gr_float* DI_density: array containing the DI (deuterium) densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` = 3.
   :param gr_float* DII_density: array containing the DII densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` = 3.
   :param gr_float* HDI_density: array containing the HD densities in code units equivalent those of the density array.  Used with ``primordial_chemistry`` = 3.
   :param gr_float* e_density: array containing the e\ :sup:`-`\  densities in code units equivalent those of the density array but normalized to the ratio of the proton to electron mass.  Used with ``primordial_chemistry`` >= 1.
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array.  Used with ``metal_cooling`` = 1.
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_cooling_time(code_units *my_units, double a_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *cooling_time);

   Calculates the instantaneous cooling time.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
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

.. c:function:: int calculate_gamma(code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *my_gamma);

   Calculates the effective adiabatic index.  This is only useful with 
   ``primordial_chemistry`` >= 2 as the only thing that alters gamma from the single 
   value is H\ :sub:`2`.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* my_gamma: array which will be filled with the calculated gamma values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_pressure(code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *pressure);

   Calculates the gas pressure.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_temperature(code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *temperature);

   Calculates the gas temperature.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* HI_density, HII_density, HM_density, HeI_density, HeII_density, HeIII_density, H2I_density, H2II_density, DI_density, DII_density, HDI_density, e_density, metal_density: arrays containing the species densities in code units equivalent those of the density array
   :param gr_float* temperature: array which will be filled with the calculated temperature values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

Tabular-Only Functions
^^^^^^^^^^^^^^^^^^^^^^

These are slimmed down functions that require ``primordial_chemistry`` = 0 
and use only the tabulated cooling rates (no chemistry).

.. c:function:: int solve_chemistry_table(code_units *my_units, double a_value, double dt_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density);

   Evolves the internal energies over a given timestep by solving the cooling 
   rate equations.  This version allows only for the use of the tabulated 
   cooling functions.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
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
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array.  Used with ``metal_cooling`` = 1.
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_cooling_time_table(code_units *my_units, double a_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density, gr_float *cooling_time);

   Calculates the instantaneous cooling time.  This version allows only for the 
   use of the tabulated cooling functions.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param int* grid_start: array holding the starting indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones
   :param int* grid_end: array holding the ending indices in each dimension of the active portion of the baryon fields.  This is used to ignore ghost zones.
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* x_velocity, y_velocity, z_velocity: arrays containing the x, y, and z velocity values in code units
   :param gr_float* metal_density: array containing the metal densities in code units equivalent those of the density array.  Used with ``metal_cooling`` = 1.
   :param gr_float* cooling_time: array which will be filled with the calculated cooling time values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_pressure_table(code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *pressure);

   Calculates the gas pressure.  This version allows only for the use of the 
   tabulated cooling functions.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

.. c:function:: int calculate_temperature_table(code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *metal_density, gr_float *temperature);

   Calculates the gas temperature.  This version allows only for the use of 
   the tabulated cooling functions.

   :param code_units* my_units: code units conversions
   :param double a_value: the expansion factor in code units (a_code = a / a_units)
   :param int grid_rank: the dimensionality of the grid
   :param int* grid_dimension: array holding the size of the baryon field in each dimension
   :param gr_float* density: array containing the density values in code units
   :param gr_float* internal_energy: array containing the specific internal energy values in code units corresponding to *erg/g*
   :param gr_float* pressure: array which will be filled with the calculated pressure values
   :rtype: int
   :returns: 1 (success) or 0 (failure)

Functions without structs (best for Fortran)
--------------------------------------------

.. c:function:: int initialize_grackle(int comoving_coordinates, double density_units, double length_units, double time_units, double velocity_units, double a_units, double a_value, int use_grackle, int with_radiative_cooling, char *grackle_data_file, int primordial_chemistry, int metal_cooling, int UVbackground, int h2_on_dust, int cmb_temperature_floor, double gamma);

.. c:function:: int solve_chemistry_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, double *a_value, double *dt_value, int *grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density);

.. c:function:: int calculate_cooling_time_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, double *a_value, int *grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *cooling_time);

.. c:function:: int calculate_gamma_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, int *grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *my_gamma);

.. c:function:: int calculate_pressure_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, int *grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *pressure);

.. c:function:: int calculate_temperature_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, int *grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *temperature);

Tabular-Only Functions
^^^^^^^^^^^^^^^^^^^^^^

These are slimmed down functions that require ``primordial_chemistry`` = 0 
and use only the tabulated cooling rates (no chemistry).

.. c:function:: int solve_chemistry_table_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, double *a_value, double *dt_value, int *grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density);

.. c:function:: int calculate_cooling_time_table_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, double *a_value, int *grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density, gr_float *cooling_time);

.. c:function:: int calculate_pressure_table_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, int *grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *pressure);

.. c:function:: int calculate_temperature_table_(int *comoving_coordinates, double *density_units, double *length_units, double *time_units, double *velocity_units, double *a_units, int *grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *metal_density, gr_float *temperature);

Internal Functions
------------------

.. c:function:: chemistry_data _set_default_chemistry_parameters();

.. c:function:: int _initialize_chemistry_data(chemistry_data *my_chemistry, code_units *my_units, double a_value);

.. c:function:: int _solve_chemistry(chemistry_data *my_chemistry, code_units *my_units, double a_value, double dt_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density);

.. c:function:: int _calculate_cooling_time(chemistry_data *my_chemistry, code_units *my_units, double a_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *cooling_time);

.. c:function:: int _calculate_gamma(chemistry_data *my_chemistry, code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *my_gamma);

.. c:function:: int _calculate_pressure(chemistry_data *my_chemistry, code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *pressure);

.. c:function:: int _calculate_temperature(chemistry_data *my_chemistry, code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *temperature);

Tabular-Only Functions
^^^^^^^^^^^^^^^^^^^^^^

These are slimmed down functions that require ``primordial_chemistry`` = 0 
and use only the tabulated cooling rates (no chemistry).

.. c:function:: int _solve_chemistry_table(chemistry_data *my_chemistry, code_units *my_units, double a_value, double dt_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density);

.. c:function:: int _calculate_cooling_time_table(chemistry_data *my_chemistry, code_units *my_units, double a_value, int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, gr_float *density, gr_float *internal_energy, gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density, gr_float *cooling_time);

.. c:function:: int _calculate_pressure_table(chemistry_data *my_chemistry, code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *pressure);

.. c:function:: int _calculate_temperature_table(chemistry_data *my_chemistry, code_units *my_units, int grid_rank, int *grid_dimension, gr_float *density, gr_float *internal_energy, gr_float *metal_density, gr_float *temperature);
