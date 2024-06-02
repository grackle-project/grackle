.. _diagnostics_debug:

Diagnostics and Debugging
=========================


Like with most computational software, one of the first steps of dealing with most (non-obvious) Grackle-related issue is to identify a minimal example that reproduces the issue (whether you want to simply get help OR you try to debug the underlying issue).

Unfortunately, a few factors have historically complicated this process:

  * Grackle has a **LOT** of input parameters

  * Grackle is usually integrated inside an application that is being executed across many processors (e.g. with MPI)

A Debugging Function: Dumping State
-----------------------------------

To help with identifying a minimal example that reproduces your issue, we have introduced the :c:func:`grunstable_h5dump_state` function.
This function is most useful when you can reliably and repeatedly reproduce the underlying issues.

This function essentially dumps information about the current state of Grackle (and the input fields) to an hdf5 file.
Therefore, you should make sure to dump the state right before the error gets triggered.
We also provide some other tools to help read back in that data (so hopefully you can reproduce the issue)

.. warning::

   At this time, the dumps shoud **NOT** be considered a stable-serialization format.
   Our priorities while writing this were to create something useful and "that works" (we didn't focus much on maintainability).
   Currently, you should think of the format as an implementation detail that can/will change at any time in a backwards incompatible manner.
   When the format does change, we will change all of the tooling in a consistent manner (but you may not be able to read a dump created with a different grackle-version).

   **PLEASE LET US KNOW** (via a GitHub Issue, the mailing list, slack, etc.) if having a stable serialization format is something you want (e.g. if you think it would help you with analysis or doing interesting science or something else).
   It might not take much effort to support this. 

How it works
^^^^^^^^^^^^

Fundamentally, issues arise due to some interaction between the input field values and Grackle's internal state.
The dump needs to encode enough information about both of things to let you reproduce the problem outside of your simulation.

Therefore, this tool **directly** dumps the data held by :c:data:`grackle_field_data` object.
There is one minor caveat: right after you first declare an instance of this type (but before you start storing values in it), please invoke be sure that you invoke the (relatively new) :c:data:`gr_initialize_field_data` to properly initialize all of the members.
(If you forget to do this, the dumping function may try to dump data from arbitrary locations in memory for unused fields.


Dumping Grackle's internal state is a little more complex.
Currently, dumping all of the internal state is intractable.
Instead the tool seeks to dump the information needed to reconstruct the current state.
The idea is that we need to dump the inputs that you, the user, specified that controls Grackle's internal state.
These inputs include:

* the chemistry-parameters, stored within the :c:data:`chemistry_data` type.
* the exact state of the :c:data:`code_units` instance initially used to call :c:func:`initialize_chemistry_data` (or to call :c:func:`local_initialize_chemistry_data` if you're using the :ref:`local_functions`)
* the current value of that the :c:data:`code_units` struct (this might be different in a cosmological simulation).

The only tricky thing here is that you need to keep around a copy of your initial :c:data:`code_units` struct that you used during initialization.
(Otherwise, we may not be able to reproduce the precise internal state).

.. note::

   It is worth mentioning that if you made any modifications to the :c:data:`chemistry_data` instance used in the call to :c:func:`initialize_chemistry_data` (or during the call to :c:func:`local_initialize_chemistry_data`), after calling that function then the precise internal state may not be reproducible.
   Likewise, any direct modifications to values stored within the :c:data:`chemistry_data_storage` struct won't be captured.

   With that said, this shouldn't happen under normal usage conditions.
   (The only reason to do this is if you're trying to hack together some desired behavior and you know what you're doing).

Example - Primary Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To help illustrate what this might look like with the :ref:`primary_functions`.

.. code-block:: c++

  #include "grackle.h"

  // it's important that the following header is present in the file where you
  // call grunstable_h5dump_state
  #include "grackle_unstable.h"

  // we might choose to track a global variable that stores the initial values
  // of the code units (that only gets written to once per simulation)
  code_units initial_code_units;


  // here is a hypothetical function where we initialize grackle
  void my_grackle_setup( /* args ... */ ) { 
    // ...

    // theres a call to `set_default_chemistry_parameters`
    // ...
    // then we store some values in the global `grackle_data` variable:
    grackle_data->use_grackle = 1;
    grackle_data->with_radiative_cooling = 1; // cooling on
    // ...

    // now let's imagine that we have a pointer to a code_units instance
    // called cur_units.
    // -> Maybe we allocated it here, passed it in as an arg or
    //    stored it as a global or something else.
    // -> let's assume it's properly initialized...
    // ...

    // now it's time to finish initializing grackle
    if (initialize_chemistry_data(cur_units) == 0) {
      fprintf(stderr, "Error in initialize_chemistry_data.\n");
      abort();
    }

    // before anything else, we should store a copy of cur_units.
    initial_code_units = *cur_units;
    // now we never modify `initial_code_units` again

    // you also shouldn't really modify the object stored within the
    // global grackle_data variable or the global gracke_rates variable
  }


  // ...

  // here is a function where we have grackle do some work where we know the
  // error will occur. For concreteness, let's imagine it happens during
  // solve_chemistry...

  void my_problematic_func(double dt, /* args ... */ ) {

    // ...

    // Assumptions:
    // - like before, cur_units is a pointer to a code_units instance (that
    //   accurately describes current code-units)
    // - my_fields is a correctly configured grackle_field_data struct

    // we need to know the destination where we'll save the file
    // -> we may need to be careful about naming (we don't want to overwrite
    //    an existing file)
    // -> also need to take care to ensure that multiple processes don't try to
    //    write to the same location (maybe you adopt a custom name?
    const char* dump_path = "path/to/dump.h5";

    // here we dump the state to disk
    // -> ideally we can exactly predict when an issue arises. For example, you
    //    know that the problem always crops up in a specific simulation domain
    //    at a specific simulation time. In this case, you would encode this
    //    requirement into the condition (and only write a single file)
    // -> If we can't exactly predict things, you may need to just write out
    //    data before every relevant grackle call and remove the file-dump if
    //    the grackle call succeeds. (You need to cleanup from the dump if you
    //    so that you don't try to write to the same location). We don't show
    //    this example here.
    if (/* condition */) {
      // now we dump to disk
      int dump_code = grunstable_h5dump_state(dump_path, -1, grackle_data,
                                              &initial_code_units,
                                              cur_units, my_fields);
      if (dump_code == 0) {
        fprintf(stderr, "ERROR in grunstable_h5dump_state\n");
        abort();
      }
      // if we expect problems to arise in solve_chemistry, (like in this
      // example), we probably want to print out that value now (to help us 
      // reproduce the problem later on).
      printf("Dump grackle state to file \"%s\" before a call to "
             "solve_chemistry with dt = %.18e\n", dump_path, dt);
    }

    if (solve_chemistry(cur_units, &my_fields, dt) == 0) {
      fprintf(stderr, "Error in solve_chemistry.\n");
      // if you are trying to diagnose some error condition that doesn't need
      // to handle, you may want to explicitly abort your program here (maybe
      // with a call to the abort() function).
      return 0;
    }
  }

Example - Local Functions
^^^^^^^^^^^^^^^^^^^^^^^^^

When using :ref:`local_functions`, the example code would basically look identical, but with 3 changes:
  1. You would be using :c:func:`local_initialize_chemistry_data` in place of :c:func:`initialize_chemistry_data`.
  2. You would be using :c:func:`local_solve_chemistry` in place of :c:func:`solve_chemistry`.
  3. Rather than passing :c:data:`grackle_data` as an argument to :c:func:`grunstable_h5dump_state`, you would pass in the same :c:type:`chemistry_data` pointer that you are using as an argument to :c:func:`local_solve_chemistry`.

Other Considerations
^^^^^^^^^^^^^^^^^^^^
As stated elsewhere, :c:func:`grunstable_h5dump_state` is technically considered unstable.
This means the function may be modified or replaced with something that is more convenient over time.
This shouldn't really affect the usability of this function.
With that said, if you want to leave this function in your simulation code to make it easy to use when when you need it, you may not want to put it in a regular if-else statement (instead, you may want to comment it out when not in use or conditionally compile).

If necessary, you can pass a ``NULL`` pointer in place of the ``initial_code_units`` struct argument.
If you do that, be aware that your issues may not be replicatable in this scenario.

Using the dumps
---------------

Now that you have the dumps in hand, how do you use them?
We have provided some utility functions within pygrackle to help with this.

Specifically, we provide the :py:func:`!load_h5dump_dicts`, :py:func:`!load_chemistry_data_from_h5dump` and :py:func:`!load_FluidContainer_from_h5dump` functions.
To use one of these functions it needs to be imported from :py:mod:`!pygrackle.debug_tools`.
The usage of these functions well documented by the docstrings in the source code.

Of the 3 functions, the easiest thing to do is use the :py:mod:`!load_FluidContainer_from_h5dump` function.
This will construct a fully configured :py:class:`!FluidContainer` from the dump, which you can easily use to try to reproduce the error.
A quick example is shown below:

.. code-block:: python

   from pygrackle.debug_tools import load_FluidContainer_from_h5dump

   path = "path/to/dump.h5"

   fc = load_FluidContainer_from_h5dump(path,
                                        flatten_field_data = True)

   # if the problem relates to solve_chemistry, then calling the following
   # should reproduce your issue
   fc.solve_chemistry(
       dt = # insert some value here
   )

If the problem occurs in some other function other than :c:func:`solve_chemistry`, than you would need to replace the ``fc.solve_chemistry(...)`` part of the above snippet.
For example, if you were diagnosing an issue in :c:func:`calculate_temperature`, you would replace ``fc.solve_chemistry(...)`` with ``fc.calculate_temperature``.

The :py:mod:`!load_FluidContainer_from_h5dump` function accepts a couple useful keyword arguments.
We highlight 2 particular choices:

  1. the ``override_grackle_data_file`` can be used to override the exact path to the grackle data file (this can be useful if you are trying to reproduce on a different machine from the one used to create the dump)

   2. the ``flatten_field_data`` can be used to flatten the shape of the field data to 1D (since pygrackle currently can't handle multi-dimensional inputs).

See the docstrings for more kwargs (and a more complete description).

.. note::

   To maximize the chance of reproducing your error, you want to try to use a version of pygrackle that was built with the exact same version of grackle as the code where the data was dumped from.

   If there have been major changes between the versions of Grackle, the python functions may not work out-of-the-box at all.

.. warning::

   We don't consider these python functions to be strictly stable.
   But, they are more stable than the dumping functions.

   In other words, there is a much lower chance that we will change these in a backwards incompatible manner (and if we do, we'll try to do it in a way that loudly complains when they change).
