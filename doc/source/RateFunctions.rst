.. _rate-functions:

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

    {RATE_NAME}_rate(double T, double units, chemistry_data *my_chemistry);

where {RATE_NAME} is the name of the coefficient you wish to calculate.

Inputs
""""""""

.. c:var:: double T

    Temperature at which you would like the coefficient to be calculated.

.. c:var:: double units

    Unit conversion factor -- will return results in cgs units when set to 1.

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

Suppose we wish to check the value of the ``k1`` rate coefficient at a temperature of 1000 K.
Given we have already configured our chemistry parameters within a ``chemistry_data`` struct
named ``my_chemistry``, we can obtain the result by simply calling the function as follows:

.. code-block:: c 

    #include "rate_coefficients.h"

    double result = k1_rate(1000., 1., *my_chemistry);

where our result will be returned in cgs units.

Example 2
""""""""""

The rate coefficient ``reHII`` can be calculated via two different methods depending on the
status of the ``CaseBRecombination`` parameter within the ``chemistry_data`` struct (named
``my_chemistry`` in this example), which can take values of either 0 or 1. Suppose we wish to compare
the results of each calculation method over an arbitrary temperature range, we can achieve this by
the following:

.. code-block:: c

    #include "rate_coefficients.h"

    // Define temperature range to calculate coefficients over
    double tempStart = 10;
    double tempEnd = 1e8:
    double numTemps = 1e3;
    double tempSpacing = (tempEnd - tempStart) / numTemps;

    // Create arrays for results storage
    double caseAResults[numTemps];
    double caseBResults[numTemps];

    // Set value of my_chemistry.CaseBRecombination
    for (int caseB = 0; caseB < 2; caseB++) {
        my_chemistry.CaseBRecombination = caseB;
        // Iterate over temperature range
        for (int i = 0; i < numTemps; i++) {
            double temp = tempStart + i*tempSpacing;
            // Store results in appropriate array
            if (caseB == 0) {
                caseAResults[i] = reHII_rate(temp, 1., *my_chemistry);
            } else {
                caseBResults[i] = reHII_rate(temp, 1., *my_chemistry);
            }
        }
    }

where we have created an array of reHII coefficients for both settings of ``chemistry_data.CaseBRecombination``
over the same temperature range.




