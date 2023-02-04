.. include:: ../../CONTRIBUTING.rst

.. _adding-new-params:

Adding a New Parameter
----------------------

  This section provides a short list of tasks to complete when a new field is introduced to the :c:data:`chemistry_data` struct:

  1. Add a new entry to ``src/clib/grackle_chemistry_data_fields.def`` for this field. This file holds a central list of fields in :c:data:`chemistry_data`. Doing this accomplishes the following three things:

     (i) the field is initialized to the appropriate default value.
     (ii) the new field and its value are printed when ``grackle_verbose = 1``
     (iii) the field can be accessed by the functions providing the dynamic access to members of :c:data:`chemistry_data` (see :ref:`dynamic-api`)
     (iv) the new field is accessible from Pygrackle (because that interface uses the functions providing dynamic access)

  2. Update the fortran definition for the ``grackle_chemistry_data`` type (in ``src/clib/grackle_fortran_interface.def``). This type must exactly match the definition of :c:data:`chemistry_data` and is used to integrate Grackle into simulation codes written in Fortran.

  3. Add docummentation in ``doc/source/Parameters.rst`` for your new parameter.
