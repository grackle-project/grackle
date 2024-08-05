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

.. _cmake_buildsystem_design_rationale:

CMake Build-System Design Rationale
-----------------------------------

The design rationale for the CMake build-system tries to walk a fine line between 2 competing philosophies:

1. We should provide a curated, convenient out-of-box experience for most people who are building Grackle.

2. The internals of a modern CMake build-system should **only** specify the firm build requirements.
   All other choices of compilation options should be specified through one of the many standardized "hooks" or options that CMake provides for customization and overriding default behavior.
   This is important for making a project easily consumable (as either a standalone project or an embedded build).

The problem with the second philopsphy is that it assumes that the person building software from source is already well-versed in CMake (e.g. a developer, someone who will package and distribute your software, someone who wants to directly embed your software within their project), or they are a motivated developer/user who can be expected to learn CMake.
It implicitly assumes that most users of a project won't ever need to directly build your software from source; they will instead install prebuilt and packaged copies of the software that are distributed through other channels (e.g. through package managers like apt, dnf, homebrew OR downloadable precompiled binaries through a website OR some kind of installer).
While this implicit assumption may be accurate for most CMake software, it obviously doesn't apply to scientific software).

It may be tempting to dismiss the second philosophy for scientific software.
In fact, it doesn't provide substantial benefits in build-systems that simply build an application (like a simulation code).
However, it provides substantial benefits for build-systems of libraries (like Grackle) that are used as components in downstream software.
It does a lot to make the libraries easier to consume.

Consequently, our CMake build-system pursues a pragmatic compromise between these philosophies.
We clearly specify all build requirements (this is compatible with both).
We provide the ability to include commonly desired compilation options (some of these may be specified hostfiles).
However, these extra (not strictly necessary) options must all be introduced in a manner that they are easily enabled/disabled, without breaking a CMake build.
This compromise allows both regular users and CMake experts (including people who want to embed Grackle into their simulation code) to easily interact with the build system.


.. _pkgconfig_rationale:

Supporting ``pkg-config``
-------------------------

This section assumes you are already familiar with :ref:`this <pkgconfig_grackle_linking>` part of the documentation.
This section exists for posterity and to explain some design choices (especially since there appears to be a dearth of authoritative documentation on this topic).

Overview
++++++++

Our main focus is ``pkg-config``'s ``--static`` flag, in the case where CMake was used to create a Grackle installation that includes the shared version and the static version of Grackle (while we focus on Grackle, this discussion has implications for other libraries).

In this scenario, the man-pages (and online discussions) for ``pkg-config`` suggest that

  * ``pkg-config --cflags grackle`` and ``pkg-config --libs grackle`` respectively supply the compiler and linker flags for consuming grackle as a shared library

  * ``pkg-config --static --cflags grackle`` and ``pkg-config --static --libs grackle`` respectively supply the compiler and linker flags for consuming grackle as a static library

Everything is perfectly fine in the first pair of invocations.
In the second pair of invocations, the compiler flags are fine (other than some unnecessary information leaking through), but there are fundamentally 2 issues associated with the linker flags:

  1. The primary issue relates to how ``pkg-config`` determines the static linker flags.
     For now, we ignore dependencies like hdf5 (more on that shortly).
     In more detail:

     * **grackle.pc** internally tracks the 2 set of linker flags.
       The first set is used in the absence of the ``--static`` flag; in this scenario, these would be used for linking grackle as a shared library.
       This currently might look something like ``-L<path/to/installed/grackle> -lgrackle``.
       The file also separately tracks *private* linker-flags that are intended to be used for this scenario when you use Grackle as a static library.
       These are somewhat specific to the compiler toolchain and machine/OS and depend on how Grackle was compiled.
       At the time of running, they *usually* specify implicit dependencies.
       For concreteness, this might look something like ``-lgomp -lgfortran -lm`` on a Linux machine (hdf5 is deliberately omitted).

     * In this scenario, ``pkg-config --static --libs grackle`` determines the list of flags by ostensibly concatenating both sets of linker flags.
       Importantly, this means that ``-lgrackle`` is specified in both cases.

     * When a linker is provided ``-lgrackle``, it will search for a file called ``libgrackle.so`` (or ``libgrackle.dynlib`` on macOS) OR ``libgrackle.a``.
       If both files are present (as is the case in this scenario), most linkers will prefer the shared library version.

     * The approach to inform the linker of a preference to use the static-library that is *most* portable and easiest-to-express is to prepend a flag like ``-Bstatic`` to the start of the flags returned by ``pkg-config --static --libs grackle``. [#pc1]_
       **Consequently, the linker will prefer to try to link every other required library statically** (we return to this shortly).
       **MORE IMPORTANTLY,** there isn't any completely portable way to do this; some linkers (namely the one on macOS [#pc2]_ ) don't provide **ANY** convenient way to express this preference.

  2. Issues dealing with dependencies like hdf5:

     * There are actually 2 ways to deal with linker flags of private dependencies inside **grackle.pc**: (i) directly specify the private linker flags (as described in the preceeding bullets) or (ii) specify a "private requirement" on the dependency.
       The latter option is only possible if the dependency provides its own \*.pc file (indeed, hdf5 usually comes with a **hdf5.pc**) file.

     * Recall from the above bullets that when using ``pkg-config``'s ``--static`` flag, we effectively need to tell the linker to prefer linking **ALL** of Grackle's dependencies statically. 
       With that in mind, we really should prefer the "private requirement" approach for specifying hdf5 as a private dependency of Grackle.
       In slightly more detail, ``libhdf5.a`` has its own set of nontrivial private dependencies, and by using the "private requirement" approach, ``pkg-config`` should automatically handle those dependencies for us.

     * Unfortunately, common hdf5 installations don't properly specify these private linker flags on commonly used platforms.
       For example, on a Debian or Ubuntu system, the **hdf5.pc** file installed alongside hdf5 with ``apt`` doesn't actually specify **ANY** private linker flags (this seems to be a common occurence).
       Additionally, the **hdf5.pc** file installed with HDF5 version 1.14.3 by Homebrew on macOS includes invalid private linker flags [#pc3]_ .

It's for the above reasons that we have choosen not to support static-linking (via ``pkg-config``) in installations where both shared and static libraries are provided.
If people want to use static libraries via pkg-config, we instead encourage an installation with just a static-library.
In an installation that just provides a static-library, the contents of **grackle.pc** are customized to dynamically link every dependency other than Grackle (and there won't be any linker-issues differentiating between the shared and static libraries). 

A potential workaround
++++++++++++++++++++++

.. COMMENT-BLOCK

   Maybe I should open an issue and record these sentiments instead of putting them into the documentation?

A potential work-around for most of the above difficulties when you have a shared library and static library installed is that we create a new file called **grackle_static.pc** that gets shipped alongside **grackle.pc**.
To get the flags for using grackle as a static library, you would invoke ``pkg-config --cflags grackle_static`` and ``pkg-config --libs grackle_static``.
Essentially, the contents of **grackle_static.pc** would look a lot like the contents of **grackle.pc** when only a static library is installed.

While this wouldn't solve the issue of telling the linker to use libgrackle.a on macOS, it would solve the other issues.
Even if we could overcome that remaining issue, this workaround is undesirable (unless there is strong user-demand) for a number of reasons:

  * the primary reason is we would have more to maintain, and this workaround is totally unnecessary outside of a somewhat pathological case (Users need to go somewhat out of their way to use the CMake build-system to create a single installation that features both the static and shared versions of Grackle).
    The much easier/more sensible/more idiomatic default behavior is to compile Grackle just as a shared library or just as a static library (in which case there isn't **any** problem).

  * this is far less "composable" than the existing alternative of requiring the end-user to compile Grackle just as a shared library or just as a static library.
    Under the existing alternative, the build-system of a downstream application will invoke the same commands commands in either scenario.
    If we shipped **grackle_static.pc**, the build system of a downstream application would need to embed additional logic to proactively alter the arguments passed to ``pkg-config`` based on the preference for shared and static libraries.

  * by the `principle of least surprise <https://en.wikipedia.org/wiki/Principle_of_least_astonishment>`__ there would be a solid argument to also supply the **grackle_static.pc** file in the case where we just install grackle as a static-library.

  * the solution wouldn't scale well.
    As development progresses, it is conceivable that we may want to easily support having multiple versions of grackle installed where we support different "backends."
    For example, we might hypothetically support ``grackle.pc``, ``grackle_openmp.pc`` and ``grackle_cuda.pc``.
    Alternatively, one could also imagine hypothetically imagine supporting ``grackle_float.pc`` (for a version compiled in single-precision).
    If we ever go down any of these roads, we would also be somewhat obligated to manually manage variants of each file with/without the ``_static`` suffix.

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

  2. Create a new PR that will serve as the final PR of the release. The PR should:

     - update the changelog (the changelog is stored within the ``CHANGELOG`` file at the root of the repository).

     - update the version number of the c-library. This is currently tracked within the ``VERSION`` file (that can be found at the root level of the repository)

     - (if applicable) update the version number of the python module (stored internally in `src/python/pygrackle/_version.py`)

  3. After the PR is merged, perform the release on GitHub.
     Make sure to associate the release with the proper target (it should include changes from the PR).

  4. Make a new PR composed of a single commit (this is the first commit of the next release).
     It should **only** update the version number of the c-library to specify the next development version
 
      - like before, this must be updated in  the ``VERSION`` file (that can be found at the root level of the repository).
        Instructions for choosing development version numbers are provided :ref:`here <dev-version-numbers>` (note that this changed in the version 3.3 release).

      - do **NOT** touch anything else (even the python version number)

  5. Go to the *Read the Docs* `project webpage <https://readthedocs.org/projects/grackle/>`__ and the new release's tag to the list of "Active Versions."
     This ensures that *Read the Docs* will maintain a version of the documentation for this particular release.

  6. Announce the release in an email to the `Grackle mailing list <https://groups.google.com/g/grackle-cooling-users>`__.
     This can be adapted from the GitHub Release notes.


.. rubric:: Footnotes

.. [#pc1] It would be ideal if ``pkg-config`` could produce something like ``-Wl,--push-state,-Bstatic -lgrackle -Wl,--pop-state``, in place of ``-lgrackle``, when the ``--static`` flag is specified.
          By doing this you would specify that you would prefer Grackle to be linked statically.
          At the same time, you would not affect the default linking behavior of any other library (the end user would be free to specify whatever they want).
          Unfortunately, this isn't possible to easily support (without requiring the invoker to supply a bunch of extra/custom information).
          Plus it would invoke unexpected, custom behavior.
          It also isn't consistent with the way that pkg-config handles private requirements (e.g. like the way that we handle hdf5).
          Furthermore, some commonly used linkers (like the ones on macOS) don't have flags that support equivalent/consistent logic.

.. [#pc2] Yes, the macOS linker technically supports a flag called ``-static``.
          But if you read the documentation (accessible with ``man ld``), you will find that this flag has some atypical behavior that doesn't really parallel something like the ``-Bstatic`` (`this article provides more details <https://michaelspangler.io/posts/statically-linking-on-macos.html>`__ )

.. [#pc3] There are a few sets of issues, that are present in the upstream HDF5 repository's CMake logic.
          Even if they do get fixed, they would take a while to propagate to a HDF5 release.
          More importantly, more errors of this kind could emerge in the future because using pkg-config with ``--static`` doesn't appear to be a common execution-path.
