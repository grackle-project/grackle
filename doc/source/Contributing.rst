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

  1. Draft the GitHub Release (using the GitHub website) and draft the changes to the changelog.

     - It's often easiest to start by drafting the GitHub release because GitHub provides the option to automatically generate a list of the titles and links to all pull-requests that have been merged since the previous release.
       That functionality provides a list of all first-time contributors.

       - We recommend looking to older release notes as a guide.
         With that said, our release generally consists of (i) a short summary about the release, (ii) a list of the changes, and (iii) a list of all contributors (that highlights first-time contributors). 

       - The automatically generated list of pull requests is an excellent way to start producing the list of changes.
         But, the list entries may need to be slightly modified.
         For example, PRs that simply update the version number tracked inside the repository should be removed (if present), it may make logical sense to group a series series of closely related PRs into a single entry, or an entry may need to be slightly more descriptive.
         We also like to sort the list-entries into sections (e.g. New Features, Minor Enhancements, Bugfixes, Documentation Updates, etc.).
         Finally, it makes sense to highlight any functions in the public API that have been deprecated or removed.

     - After you draft the GitHub release, you can take advantage of GitHub's feature to save the draft and circulate it to other developers for comments/recommendations.

     - It's fairly straightforward to draft changes to the changelog based on the contents of the GitHub release.

  2. Create a new PR that includes 2 commits will serve as the final PR of the release. It should:

     - update the changelog (the changelog is stored within the ``CHANGELOG`` file at the root of the repository).

     - update the version number of the c-library. This is currently tracked within the ``VERSION`` file (that can be found at the root level of the repository)
     - (if applicable) update the version number of the python module (stored internally in `src/python/setup.py`)

     (ii) Commit B (the second commit) is the first commit for the next version.
          It should just update the version number of the c-library to specify the next development version

          -  like before, this must be updated in  the ``VERSION`` file (that can be found at the root level of the repository)
          - do **NOT** touch anything else (even the python version number)



  3. After the PR is merged, perform the release on GitHub.
     Make sure to associate the release with the proper target (it should include changes from the PR).

  4. Make a new PR composed of a single commit (this is the first commit of the next release).
     It should **only** update the version number of the c-library to specify the next development version

      -  like before, this must be updated in  the ``VERSION`` file (that can be found at the root level of the repository)
      - do **NOT** touch anything else (even the python version number)

  5. Go to the *Read the Docs* `project webpage <https://readthedocs.org/projects/grackle/>`__ and the new release's tag to the list of "Active Versions."
     This ensures that *Read the Docs* will maintain a version of the documentation for this particular release.

  6. Announce the release in an email to the `Grackle mailing list <https://groups.google.com/g/grackle-cooling-users>`__.
     This can be adapted from the GitHub Release notes.
