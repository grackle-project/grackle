.. _rate-functions:

.. role:: c_inline(code)
   :language: c

Rate Functions
=========================

Grackle supports the calculation of individual rate coefficients at specific temperatures.
When running Grackle's chemistry solving schemes, these coefficients are computed internally,
at the temperatures required by the solver. However, there may be cases where obtaining a
specific rate coefficient at a given temperature is desirable -- with Grackle's initialisation
routines ported from Fortran to C, this is now possible. Documented here is an overview of
how to use these functions, alongside relevant examples.

All rate coefficients are designated an individual function for their calculation; their definitions
are located within ``rate_coefficients.c`` and their prototypes in ``rate_coefficients.h``, both found 
within the ``grackle/src/clib/`` directory. The vast majority of the rate coefficient's functions are
of the same general form, however there are exceptions which will be discussed individually.

The General Rate Function
---------------------------

Structure
^^^^^^^^^^

As mentioned previously, most rate functions have the same structure and can therefore be called
identically:

.. code-block:: c

    double {RATE_NAME}_rate(double T, double units, chemistry_data *my_chemistry);

where {RATE_NAME} is the name of the coefficient you wish to calculate.

Inputs
""""""""

.. c:var:: double T

    Gas temperature at which you would like the coefficient to be calculated.

.. c:var:: double units

    The calculated rate will be divided through by this number; thereby facilitating
    the use of arbitrary units. If set to 1, the rate will be calculated in CGS units
    -- the units in which the rate equation is defined.

.. c:var:: chemistry_data *my_chemistry

    Pointer to the chemistry_data struct containing the parameters for your calculations.

Outputs
"""""""""

.. c:var:: double rate

    Rate coefficient calculated for the specified input parameters.

Examples
^^^^^^^^^^

Example 1
""""""""""

Suppose we wish to print the value of the ``k1`` rate coefficient at a temperature of 1000 K.
Given we have already configured our chemistry parameters within a ``chemistry_data`` struct
named ``my_chemistry``, we can obtain the result by simply calling the function as follows:

.. code-block:: c 

    #include "grackle_rate_functions.h"

    double result = k1_rate(1e3, 1., &my_chemistry);
    printf("k1 %e", result);

where our result is in cgs units.

Example 2
""""""""""

The rate coefficient ``reHII`` can be calculated via two different methods depending on the
status of the ``CaseBRecombination`` parameter within the ``chemistry_data`` struct (named
``my_chemistry`` in this example), which can take values of either 0 or 1. Suppose we wish to compare
the results of each calculation method over an arbitrary temperature range, we can achieve this by
the following:

.. code-block:: c

    #include "grackle_rate_functions.h"

    /// Define temperature range to calculate coefficients over
    double tempStart = 1e1;
    double tempEnd = 1e8;
    double numTemps = 1e3;
    double tempSpacing = (tempEnd - tempStart) / numTemps;

    // Create arrays for results storage
    double caseAResults[(int) numTemps];
    double caseBResults[(int) numTemps];

    // Set value of my_chemistry.CaseBRecombination
    for (int caseB = 0; caseB < 2; caseB++) {
        my_chemistry.CaseBRecombination = caseB;
        // Iterate over temperature range
        for (int i = 0; i < numTemps; i++) {
            double temp = tempStart + i*tempSpacing;
            // Store results in appropriate array
            if (caseB == 0) {
                caseAResults[i] = reHII_rate(temp, 1., &my_chemistry);
            } else {
                caseBResults[i] = reHII_rate(temp, 1., &my_chemistry);
            }
        }
    }

where we have created an array of reHII coefficients for both settings of ``chemistry_data.CaseBRecombination``
over the same temperature range.

The k13dd Rate Function
-------------------------

Structure
^^^^^^^^^^

The ``k13dd`` rate function, which describes the density-dependent dissociation of molecular hydrogen, is similar
in form to the general rate functions, the only difference being its additional input parameter. This is a pointer
to an array of length ``14 * sizeof(double)``, which will hold the outputs of the function. The function
always calculates fourteen rate parameters, the first seven of which correspond to direct collisional dissociation,
whilst the latter seven correspond to dissociative tunneling -- please see
`Martin, Schwarz & Mandy, 1996 <http://adsabs.harvard.edu/pdf/1996ApJ...461..265M>`_ for further details on how 
these are calculated. The structure of the function is then:

.. code-block:: c

    void k13dd_rate(double T, double units, double *results_array, chemistry_data *my_chemistry);

Inputs
""""""""

.. c:var:: double T

    Gas temperature at which you would like the coefficient to be calculated.

.. c:var:: double units

    The calculated rate will be divided through by this number; thereby facilitating
    the use of arbitrary units. If set to 1, the rate will be calculated in CGS units
    -- the units in which the rate equation is defined.

.. c:var:: double *results_array

    Pointer to array of length :c_inline:`14 * sizeof(double)` in which the calculated rate coefficients will
    be stored.

.. c:var:: chemistry_data *my_chemistry

    Pointer to the chemistry_data struct containing the parameters for your calculations.

Outputs
"""""""""

.. c:var:: None

    Results are stored within results_array, function itself is void.

Examples
^^^^^^^^^^

Example 1
""""""""""

Suppose we would like to print the rate coefficients for the dissociation of molecular hydrogen via the 
tunneling process at a temperature of 1e5 K. Given we have already configured our chemistry parameters
within a ``chemistry_data`` struct named ``my_chemistry``, we can obtain the coefficients by the following:

.. code-block:: c

    #include "grackle_rate_functions.h"

    // Create an array of the correct size for result storage.
    double results[14];

    // Call the function at the desired temperature, getting results in cgs units.
    k13dd_rate(1e5, 1., &results, &my_chemistry);

    // Print the results corresponding to dissociative tunneling.
    for (int i = 7; i < 14; i++) {
        printf("k13dd %e", results[i]);
    }

The h2dust Rate Function
-------------------------

Structure
^^^^^^^^^^

The ``h2dust`` rate function, which describes the formation of molecular hydrogen on dust grains, is similar
in form to the general rate functions, the only difference being its additional input parameter; a
``double`` which represents the dust temperature. The function returns a ``double`` just as the
general rate function, its structure is then:

.. code-block:: c

    double h2dust_rate(double T, double T_dust, double units, chemistry_data *my_chemistry);

Inputs
""""""""

.. c:var:: double T

    Gas temperature at which you would like the coefficient to be calculated.

.. c:var:: double T_dust

    Dust temperature at which you would like the coefficient to be calculated.

.. c:var:: double units

    The calculated rate will be divided through by this number; thereby facilitating
    the use of arbitrary units. If set to 1, the rate will be calculated in CGS units
    -- the units in which the rate equation is defined.

.. c:var:: chemistry_data *my_chemistry

    Pointer to the chemistry_data struct containing the parameters for your calculations.

Outputs
"""""""""

.. c:var:: double rate

    The rate coefficient for the h2dust reaction at the specified input parameters.

Examples
^^^^^^^^^^

Example 1
""""""""""

Suppose we would like to calculate the ``h2dust`` rate coefficients for a gas temperature of 1e4 K, with a 
varying dust temperature. Given we have already configured our chemistry parameters within a ``chemistry_data``
struct named ``my_chemistry``, we can obtain the coefficients by the following:

.. code-block:: c

    #include "grackle_rate_functions.h"

    // Define dust temperature range to calculate coefficients over.
    double tempStart_dust = 1e1;
    double tempEnd_dust = 1e6;
    double numTemps_dust = 1e3;
    double tempSpacing_dust = (tempEnd_dust - tempStart_dust) / numTemps_dust;

    // Create array for results storage.
    double h2dust_results[(int) numTemps_dust];

    // Loop over dust temperatures.
    for (int i=0; i < numTemps_dust; i++){
        double temp_dust = tempStart_dust + i*tempSpacing_dust;
        h2dust_results[i] = h2dust_rate(1e4, temp_dust, 1., &my_chemistry);
    }
    
The Scalar Rate Functions
---------------------------

Structure
^^^^^^^^^^
The scalar rate functions (``comp``, ``gammah``, ``gamma_isrf``) are simpler than the general rate functions
due to their temperature independence. They require only two inputs and return a single ``double``,
their structure is as follows:

.. code-block:: c

    double {SCALAR_NAME}_rate(double units, chemistry_data *my_chemistry);

where {SCALAR_NAME} is the name of the scalar rate coefficient you wish to calculate. These are
called in the same way as the general rate functions, ignoring the temperature dependancy -- 
please see their documentation for basic examples.

Inputs
""""""""

.. c:var:: double units

    The calculated rate will be divided through by this number; thereby facilitating
    the use of arbitrary units. If set to 1, the rate will be calculated in CGS units
    -- the units in which the rate equation is defined.

.. c:var:: chemistry_data *my_chemistry

    Pointer to the chemistry_data struct containing the parameters for your calculations.

Outputs
"""""""""

.. c:var:: double rate

    The rate coefficient for the specified chemistry parameters.

