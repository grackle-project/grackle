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

  3. Add documentation in ``doc/source/Parameters.rst`` for your new parameter.

Instructions for Making a Release
---------------------------------

Before reading this section, please ensure you are familiar with our :ref:`versioning policy <versioning-code>`.

To create a new release:

  1. Draft the changes to the changelog and draft the release notes (GitHub allows you to circulate a draft of the release notes with other developers).

  2. Create a new PR that includes 2 commits (commit A followed by commit B):

     (i) Commit A (the first commit) is the final commit included in the release. It should:

         - update the changelog
	 - update the version number of the c-library. This is currently tracked in with `src/clib/Make.config.assemble` **AND** (for documentation purposes) in `doc/source/conf.py`
         - update the version number of the python module (stored internally in `src/python/setup.py`)

     (ii) Commit B (the second commit) is the first commit for the next version.
          It should just update the version number of the c-library to specify the next development version

	  -  like before, this must be updated in `src/clib/Make.config.assemble` **and** `doc/source/conf.py`
          - do **NOT** touch anything else (even the python version number)

  3. After the PR is merged, perform the release on GitHub.
     Make sure to associate the release (and the new tag) with Commit A (not with Commit B).

  4. Announce the release in an email to the `Grackle mailing list <https://groups.google.com/g/grackle-cooling-users>`__.
