
.. _versioning-code:

Versioning Policy
=================

This "policy" is primarily intended to codify our historic treatment of Grackle versioning.
**Please reach out** if you have any concerns about our current policy or suggestions about how it can be improved (especially if it affects your willingness to integrate Grackle into your simulation).

We always encourage users to download the latest commit when installing Grackle.
Our general policy to keep Grackle's main branch stable and in working order.

Periodically, we "release" new versions of Grackle.
A release serves 2 purposes:

  1. It's way to designate a development milestone (whether it's new functionality, changes in the public API, bugfixes) at a snapshot in the project's development.
  2. It's a way to announce the latest changes to the library since the last milestone.

Currently, the repository tracks 2 separate version numbers: the primary version number (used by the documentation and the shared/static library) **and** the pygrackle version number.

.. _primary-version-number:

Primary Version Number
----------------------

It is possible to :ref:`query <query-version>` this version number at runtime.

In the commit associated with a release, the repository stores the release version number.
This number follows the rules of `semantic versioning <https://semver.org/>`__ and has the form ``<MAJOR>.<MINOR>(.<MICRO>)``.
In this template, ``<MAJOR>``, ``<MINOR>``, and ``<MICRO>`` correspond to major, minor, and micro version numbers (the micro version number may be omitted if it’s zero).

At any time other than a release, the repository stores a "development version number."
Since version the ``3.3`` release, the development version number also follows the rules of semantic-versioning and has the form ``<MAJOR>.<MINOR>(.<MICRO>)-dev``.
More detail is provided :ref:`below <dev-version-numbers>`.

.. COMMENT BLOCK

   If we were following the rules of semantic versioning as strictly
   as possible, we would NEVER remove deprecated functionality except
   when incrementing the major version number.  For example, in the
   release of Grackle 3.0, we marked some functions as deprecated.
   Technically, if we were adhering to the rules of semantic
   versioning as strictly as possible, we would have waited until
   Grackle 4.0 to actually remove those functions.  The main reason to
   do this is if we cared about maintaining ABI compatability between
   minor versions (it would make the shared library numbering make
   more sense)

   With all of that said, it's still technically correct for us to say
   that we follow the semantic versioning specification. An explicit
   exception appears to be carved out for removing deprecated
   functionality.

A change in the major version number typicially denotes significant changes in functionality.
It must also be incremented in a new release if a backwards incompatible change has been made to the public API since the last release (alternatively, parts of the public API can be deprecated in a major release and be removed later on in a minor release).

A change in the minor version number can indicate changes in functionality.
It must be incremented in a new release if a backwards-compatable change has been made to the public API since the last release.

A change in the micro version number denotes a bugfix-release.
In other words, it must only be incremented if the only changes since the last release are bugfixes that have not modified the public API in any way.

.. _dev-version-numbers:

Development Version Number Conventions
++++++++++++++++++++++++++++++++++++++

Starting in version ``3.3``, we shifted strategies for incrementing development version numbers.
Our strategy has 2 strict rules:

1. The development version immediately after a release must take the version number from that release, increments the micro number and appends the ``-dev`` suffix.

   - In other words, ``3.3.1-dev`` **MUST** be the development version number immediately after the ``3.3`` release

   - Likewise, after a hypothetical ``3.3.1`` bugfix-release, ``3.3.2-dev`` **MUST** be the next development version number.

2. We are free to increment (but not required) development versions between releases, but we must be sure that the leading numbers of a development-version provide the lower bound on the next allowed release version.
   For example:

   - In version ``3.3.1-dev``, if we added a major new feature, and we were confident the next release had to be a minor release, we could introduce a new development-version with an incremented minor version.
     The resulting version sequence is: ``3.3``, ``3.3.1-dev``, ``3.4-dev``, ``3.4`` ...

   - As noted above this strategy does not require us to introduce a new development-version (e.g. ``3.3``, ``3.3.1-dev``, ``3.4``, ... is allowed). 
     With that said, we will make an effort to increment development-version numbers when deprecated functionality is removed.

Under the premise that releases always correspond to a single commit on the main Grackle branch, this new strategy ensures that development version numbers always satisfy version ordering requirements of semantic versioning.

Prior, to the ``3.3`` release, development versions could violate the ordering requirements. [#f1]_ 
That versioning strategy was better suited for the alternative development-approach where bugfix releases corresponded to commits that on a branch separated from the main development branch.

pygrackle Version Number
------------------------

This is the version number associated with the python bindings.
We have only recently begun incrementing this number in recent releases.


Other Notes
-----------

We are strongly commited to maintaining backwards compatability with the public API.
Once you write a downstream application that makes use of Grackle, they can freely rebuild the application using newer versions of Grackle, without modifying any of the source code and Grackle will **always** [#f2]_ provide consistent results.


At the time of writing, the current design of Grackle's API is not conducive to supporting a stable ABI. [#f3]_
A stable ABI is **only** needed to support a special shared library feature:  if you compile & link a downstream application against one version of a shared library, you replace the shared library with a different version (that is ABI compatible), and then run the downstream application without any changes.
This all just means that you need to simply recompile your downstream application any time you upgrade Grackle.
This isn't a big loss since simulation codes are already regularly recompiled.
More importantly you should probably be using Grackle as a static library, rather than as a shared library, to maximize speed (in that case, you would need to recompile your application when you upgrade Grackle, even if Grackle supported a stable ABI).



.. rubric:: Footnotes:

.. [#f1] Prior to version ``3.3``, the developmemt version number was allowed to violate the version-ordering rules outlined in the semantic versioning specification.
         These development numbers had the form ``<MAJOR>.<MINOR>(.<MICRO>).dev<DEV_NUM>`` (the suffix was always ``.dev`` whereas our newer strategy ends with ``-dev``).
         For reference, the **ORDERED** list of version numbers from version ``3.0`` through ``3.3`` is: ``3.0``, ``3.1.dev1``, ``3.1``, ``3.2.dev1``, ``3.1.1``, ``3.2.dev2``, ``3.2.0``, ``3.3.dev0``, ``3.2.1``, ``3.3.dev1``, ``3.3``.

.. [#f2] The main exception is if we deprecate part of the API.
         If we do that, we will be sure to announce the deprecation (and provide migration guidelines) in a release **long** before we remove that functionality.

.. [#f3] At the time of writing, the API requires the user to explicitly allocate and deallocate the memory used to hold instances of the :c:type:`chemistry_data`, :c:type:`code_units`, :c:type:`grackle_field_data` and (when using the *local* api) the :c:type:`chemistry_data_storage` types, which are all implemented as structs.
         Minor modifications to Grackle commonly involve introducing new struct-members to the underlying structs.
         As noted above, great care is taken so that you can freely compile an application against a newer versions of Grackle without making any modifications to the application's source code (i.e. the API is backwards compatabile).
         However, it's important to make sure that the header-files shipped with the newer version of Grackle are installed when doing this so that the compiler can properly infer that how much memory to allocate **(If you don't know what this all means, don't worry about it. You need to go out of your way to do this wrong)**.

.. COMMENT BLOCK

   If Grackle provided functions for allocating and freeing the memory
   for these structs, that would go a long ways towards acheiving ABI
   stability. I THINK that would technically be enough as long as we
   always made sure to add new struct-members to the end of the struct
   and we NEVER reordered or removed existing struct-members of ANY
   struct declared in the public header.

   In practice, a more robust solution would involve removing all
   structs from the public headers and then implementing functions
   like the dynamic API (or getter/setters) for accessing the members
   of all structs. libpng does something like this

   
