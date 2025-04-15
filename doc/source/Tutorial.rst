Tutorial
========

Overview
--------
In this tutorial we walk through a minimal example of using Grackle to perform radiative cooling.
We provide a more detailed guide for using Grackle's API :ref:`here <integration>`.

Background
----------

It's instructive to understand the Grackle library universally adopts the convention:

  **fluid state is primarily represented by specific internal energy and mass density.**

These quantities are always used to describe the fluid's state when calling Grackle, Grackle uses them to derive other properties (e.g. temperature/pressure) and they are evolved by Grackle (when evolving chemistry/cooling).

Specific internal energy, :math:`e`, is a temperature-like quantity that characterizes the fluid's thermodynamic state.
For a calorically perfect ideal gas,\ [#perfect-ideal]_  :math:`e = k_B T/ ((\gamma -1) m_H \mu)`, where :math:`\gamma` is adiabatic index, :math:`m_H` is Hydrogen mass, :math:`k_B` is the Boltzmann constant, :math:`\mu` is the (composition-dependent) mean molecular weight.

The total mass density, :math:`\rho`, describes the "amount" of a fluid.
Other mass densities can be used to describe the fluid's composition (not relevant to the current example).
Note: Grackle never evolves the total mass density.

The Scenario
------------

We calculate the specific internal energy of a gas parcel after it has undergone radiative cooling for :math:`10^13\ {\rm s}` (or :math:`{\approx}0.32` Myr).
The mass density is :math:`\rho = 5 \times 10^{-24}\ {\rm g}\ {\rm cm}^{-3}` (or :math:`{\approx}3` Hydrogen masses per cubic centimeter), initial specific internal energy is :math:`e_0 = 2\times 10^{12}\ {\rm cm}^2\ {\rm s}^{-2}`, and :math:`\gamma = 5/3`.

We make 3 simplifying assumptions about the parcel's composition:\ [#composition-assumptions]_ (i) there are no metals (i.e. the heaviest particle is He), (ii) the parcel is always in ionization equilibrium, and (iii) the Hydrogen to Helium ratio is described :any:`here <HydrogenFractionByMass>`.

Under these assumptions, the initial temperature is :math:`T\approx 1.5\times 10^4\ {\rm K}`.

The Example Program
-------------------

The following example program performs the calculation described just above.
The program is intended to be compiled with a C99 compiler (or newer).
To keep the program short, the error-handling is fairly crude.

.. gr-include-snippet:: ../../src/example/c99_simple_example.c
   :language: c

This example illustrates some essential aspects of using Grackle.
We highlight these aspects down blow.

Setting Up the Scenario
^^^^^^^^^^^^^^^^^^^^^^^

At the start of :c:func:`!main`, we define some variables that store inputs for our calculation.
The details are influenced by conventions of simulation codes.

.. literalinclude:: ../../src/example/c99_simple_example.c
   :language: c
   :start-after: //@![[[ANCHOR:ScenarioSetup]]]
   :end-before: //@![[[ANCHOR:SanityChecks]]]

On the first line, we define variables, (:c:var:`!densU`, :c:var:`!lenU`, and :c:var:`!timeU`) that describe a custom unit system [#about-code-units]_ that our sample-program uses to encode the input quantities.
We define a custom unit-system because it's standard convention *and* Grackle assumes that input values are defined with an :ref:`appropriate unit-system <choosing_appropriate_units>`.\ [#require-reasonable-units]_
Our unit system use rounded numbers roughly corresponding to a density unit of 1 Hydrogen masses per cubic centimeter, a length unit of kpc, and a time unit of Mpc.

On the second line, we define :c:var:`!energyU` to hold the conversion factor for specific energy units.
(**BE AWARE:** Grackle internally relates specific energy-units, length-units and time-units somewhat differently when using comoving coordinates).

Next, we define variables that hold data describing fluid properties.
We crudely approximate the standard pattern (sometimes called `struct of arrays <https://en.wikipedia.org/wiki/AoS_and_SoA>`_ or `parallel arrays <https://en.wikipedia.org/wiki/Parallel_array>_`) from grid-codes: each quantity is tracked in a distinct 1D buffer (usually called a field).
Each location in a field typically maps to 3D indices (corresponding to 3D space).
Since our example-program only considers a single parcel of gas, we store data in buffers corresponding to a ``(1,1,1)`` field.

Finally, the :c:var:`!dt` variable holds the timestep we will evolve the timestep over.

Sanity Checks
^^^^^^^^^^^^^

The next chunk of logic performs some quick sanity checks.
It confirms that Grackle configured with the correct precision and that it is properly linked.

.. literalinclude:: ../../src/example/c99_simple_example.c
   :language: c
   :start-after: //@![[[ANCHOR:SanityChecks]]]
   :end-before: //@![[[ANCHOR:SetupChemistrySolver]]]

Configuring the Solver
^^^^^^^^^^^^^^^^^^^^^^
Next, we move onto configuring the Grackle solver.

.. literalinclude:: ../../src/example/c99_simple_example.c
   :language: c
   :start-after: //@![[[ANCHOR:SetupChemistrySolver]]]
   :end-before: //@![[[ANCHOR:SetupGrackleFieldData]]]

The :c:type:`chemistry_data` type encodes runtime parameter, :c:type:`code_units` is used to describe the chosen unit-system, and :c:type:`chemistry_data_storage` is used to manage internal data.

.. note::

   As long as you don't mutate these data-structures after you configure them, they can be used in as many operations as you want.
   If you want to mutate anything, you must generally configure fresh data structures.
   :ref:`The linked page <integration>` describes a few exceptions to this rule (they're mostly related to cosmology).

Preparing to communicate Field Data to Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At this point, we are ready to setup an instance of :c:type:`grackle_field_data`, which is used to communicate the field data to the Grackle solver.

.. literalinclude:: ../../src/example/c99_simple_example.c
   :language: c
   :start-after: //@![[[ANCHOR:SetupGrackleFieldData]]]
   :end-before: //@![[[ANCHOR:ApplySolver]]]

In this snippet, we have configured :c:var:`grid_start` and :c:var:`grid_end` such that the solver operates on every fluid element (yes, exactly 1 fluid element).

Calling the Solver
^^^^^^^^^^^^^^^^^^

After all that, we can actually run the calculation and print the result

.. literalinclude:: ../../src/example/c99_simple_example.c
   :language: c
   :start-after: //@![[[ANCHOR:ApplySolver]]]
   :end-before: //@![[[ANCHOR:Exit]]]

Exitting
^^^^^^^^
Since the program is finished, we can exit (we cut some corners here with cleanup).

Running the Example
-------------------

After you perform a CMake-build in a build-directory called ``<build-dir>``, you can execute:

.. code-block:: shell-session

   <build-dir>/examples $ ./c99_simple_example

.. warning::

   The examples make certain assumptions about the location of the input files.
   The examples are only guaranteed to work if both:

      1. you execute the example-binary from the same-directory where the example-binary is found

      2. ``<build-dir>`` is a top-level directory in the grackle repository (e.g. something like ``my-build`` is fine, but choices like ``../my-grackle-build`` and ``my_builds/my-first-build`` are problematic).

.. rubric:: Footnotes

.. [#perfect-ideal] A calorically perfect ideal gas, satisfies :math:`p = (\gamma - 1) \rho e` and the ideal gas law.

.. [#composition-assumptions] To be clear, we're only making these assumptions to simplify the example.
     Grackle is perfectly capable of relaxing these assumptions.

.. [#about-code-units] In this case, :c:var:`!timeU` stores the number of seconds in 1 time unit, :c:var:`!lenU` holds the number of cm in 1 length unit, and :c:var:`densU` holds the number of :math:`{\rm g}\ {\rm cm}^{-3}` in 1 density unit.
     For the uninitiated, it is standard practice for simulation codes to define a custom unit system.
     This practice aims to avoid numerical issue by adopting a unit-system such that the :ref:`the density, length, and time values are always close to 1 <choosing_appropriate_units>`.

.. [#require-reasonable-units] Grackle has special logic to deal with cases where physical quantities are "tiny."
     For historical reasons, Grackle identifies "tiny" quantities using hardcoded values assuming that they specified with an :ref:`appropriate unit system <choosing_appropriate_units>`.
     This is something we would like to change in the future.
