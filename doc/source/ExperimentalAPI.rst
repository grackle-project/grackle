.. _experimental-api:

Experimental API
================

This document describes prototypes of a new experimental API we are working on for future versions of Grackle.
The goal is to eventually release a major version of Grackle that shifts to using this API.

All functions in the dynamic API are prefixed with ``grunstable``.
The prefix will change to ``gr`` upon stabilization.

.. caution::

   The experimental api is not yet stable. The interface may change at any time without warning.

.. important::

   As a general rule, the experimental API will not be mixed with the existing API.
   However, an exception may be made for :c:type:`gr_fields`, since we will be doubling the number of existing fields.

key innovation: opaque pointers
-------------------------------

A key innovation of the experimental API is our use of `opaque pointers <https://en.wikipedia.org/wiki/Opaque_pointer>`__.
It's instructive to consider how the existing API uses transparent datatypes: the representation of datatypes (e.g. :c:type:`grackle_field_data`) is exposed and client-code is responsible for directly mutating datatype's contents.
In contrast, an API using an opaque pointer refers to a forward-declared type/struct whose definition is **EXPLICITLY** hidden from the user.
The only way to create, destroy, modify, or query the referenced type is using API functions.

Opaque pointers are a very common method of "information hiding" in C APIs.
The most famous example is ``FILE*`` in the standard library.
But it also shows up in other well-known libraries (e.g. sqlite3, libgit2, zlib).

The use of opaque pointers enables us to write an experimental API that can:

1. maintain ABI-compatability when adding a new runtime parameter or field (we may choose not to do this, but it isn't even an option with the current API) [#f1]_
2. allow simulation codes to use newly added parameters/fields in some configurations without requiring conditional-compilation to maintain compatibility with older grackle in other configurations (when the newer parameters/fields aren't used). [#f2]_
3. directly initialize internal data structures that are used throughout Grackle, while retaining the flexibility to refactor the data structures.
   In principle, the internal datatypes could even be C++ classes.
4. gracefully deprecate problematic runtime parameters over an extended deprecation period.
   Importantly, such an API can support gracefully renaming or consolidation of parameters (e.g. during the deprecation period, we could implement logic for gracefully remapping deprecated parameters).
   Under the current system, we can't do this without producing cryptic compiler errors.

Tradeoffs
---------

In general, the use of opaque pointers introduce 2 tradeoffs:
1. client-code theoretically has less control over memory allocation strategies
2. the extra indirection introduces some amount of overhead.

Tradeoff #1 shouldn't be a serious concern of existing Grackle-user (historically, all of the internal Fortran subroutines always performed a lot of internal heap allocations).
With that said, we could address it by allowing client-code to specify to user-specified allocation/deallocation callback functions (this is a pretty standard solution).

While Tradeoff #2 true, but it should only come up when you are configuring Grackle.
The amount of overhead should be negligible compared to the time that Grackle spends performing calculations in most client code.
However, if your code commonly makes lots of calls to grackle for very small fields in quick succession (you will already see lots of overhead with the current api in this scenario), some extra care may be required while interfacing with :c:type:`gr_fields` type.

.. note::

   In general, we are not opposed to supporting a secondary API that exposes more information in order to reduce overhead.
   But, it is not currently a priority.
   Such an API will require a lot of care (i.e. to ensure we don't lock into bad design decisions) and it can't be finailized until a bunch of active development (e.g. Fortran transcription, stablizing the experimental API, adding GPU support) has progressed.

Barriers to Stablization
------------------------

The biggest barrier to stablization relates to improved error-handling.
A GitHub issue will be created in the near future discussing this topic.
Our decision may alter the return-type of each of the proposed functions.

The ``gr_fields`` type
----------------------

As a replacement for :c:type:`grackle_field_data`, we are proposing the :c:type:`gr_fields` type.
The type is operated on functions with names following the template :c:func:`!grunstableFields_<operation>` (after stablization, the template will become :c:func:`!grFields_<operation>`).

.. c:type:: gr_fields

   An opaque type that aggregates field data.

   .. admonition:: Stablization-Consideration

      Should we rename this :c:type:`!gr_field` (no trailing ``s``)? (We would also need to change the prefix).


Just like :c:type:`grackle_field_data`, the :c:type:`gr_fields` type tracks common properties shared among each field and "slots" that hold buffers of user-specified data for each field.
The common properties are specified by the equivalent of :ref:`"setter-methods." <gr_fields-Setting-Grid-Props>`.

In contrast to the way that :c:type:`grackle_field_data` always provides "slots" for a fixed set of fields (regardless of whether the fields are actually used), :c:type:`gr_fields` **ONLY** provides slots for the fields required by the associated Grackle configuration.
In other words, :c:type:`gr_fields` won't let you access a H\ :sub:`2` density if you configured Grackle to use tabulated cooling.
:c:type:`gr_fields` has associated functions :ref:`to query the fields <gr_fields-Query-Field-Names>` that a given instance supports.
The actual "slots" are accessed through functions that implement a dictionary-like interface that maps field-names to "slots".

.. _gr_fields-Mode-Description:

The way that the "slots" are actually used depends on how :c:type:`gr_fields` is configured.
An instance can be configured in one of 2 ways:

  1. The standard (default), "unmanaged," mode where a :c:type:`gr_fields` instances **DOES NOT** manage the lifetime of the memory referenced by a given "slot."
     In this mode, a simulation code is responsible register pointers to memory that they externally manage in each slot.
     This is consistent with the way that :c:type:`grackle_field_data` works.
  2. An experimental "managed" mode, where a :c:type:`gr_fields` instances automatically manages lifetime of each memory buffer referenced by a "slot."
     In this mode, a simulation code must copy the field data into the buffers tracked by the "slots" (when Grackle is used to integrate chemistry/cooling, the simulation code also need to copy updated values out of the buffers).
     This is primarily intended as a convenience for simulations codes that store fluid data in arrays of structs (it may also simplify testing).

.. _gr_fields-Create-Destroy-Section:

Creation and Destruction
++++++++++++++++++++++++

In this section, we cover the functions responsible for creation and destruction of :c:type:`gr_fields`.
We break this discussion into a description of the :ref:`creation flags<gr_fields-CreationFlags>` and the :ref:`API reference for the relevant functions <gr_fields-Create-Destroy-FuncRef>`.

Before we get to those points, it is instructive to briefly consider some illustrative code-snippets for initializing an instance.
In these snippets, you should assume that ``my_chemistry`` and ``my_rates`` are pointers to :c:type:`chemistry_data` and :c:type:`chemistry_data_storage`.

.. tabs::

   .. tab:: "Unmanaged" Mode

      This snippet illustrates code that creates an instance in the standard, :ref:`unmanaged mode<gr_fields-Mode-Description>`

      .. code-block:: c

         gr_fields* fields = NULL; 
         int rc = grunstableFields_init_from_local(
           &fields, my_chemistry, my_rates, 0, 0
         );
         if (rc != GR_SUCCESS) { /*ERROR-HANDLING GOES HERE*/ }

   .. tab:: "Managed" Mode

      This snippet illustrates code that creates an instance in the :ref:`managed mode<gr_fields-Mode-Description>` (it allocates space for 64 elements per field).

      .. code-block:: c

         gr_fields* fields = NULL; 
         int rc = grunstableFields_init_from_local(
           &fields, my_chemistry, my_rates, GR_FCREATE_MANAGED_FIELDDATA, 64
         );
         if (rc != GR_SUCCESS) { /*ERROR-HANDLING GOES HERE*/ }

   .. tab:: Multiple Creation options

      This snippet illustrates code that creates an instance with multiple creation options.
      This is simple: just use *Bitwise-OR* on the flags.
      For concreteness, let's imagine we want an instance

        - in :ref:`managed mode<gr_fields-Mode-Description>` (space is allocated for 64 elements per field).
        - with whatever property is represented by the dummy :c:macro:`GR_FCREATE_DUMMY` placeholder flag.

      .. code-block:: c

         gr_fields* fields = NULL; 
         int rc = grunstableFields_init_from_local(
           &fields, my_chemistry, my_rates,
           GR_FCREATE_MANAGED_FIELDDATA | GR_FCREATE_DUMMY,
           64
         );
         if (rc != GR_SUCCESS) { /*ERROR-HANDLING GOES HERE*/ }

.. _gr_fields-CreationFlags:

Creation Options
^^^^^^^^^^^^^^^^
If you don't want to use any flags, you should use a value of zero.

.. c:macro:: GR_FCREATE_DUMMY

   This is a dummy placeholder that we present here for illustrative purposes.
   This will never actually be defined by Grackle.

.. c:macro:: GR_FCREATE_MANAGED_FIELDDATA

   Creation flag specifying that the :c:type:`gr_fields` instance should be initialized in the :ref:`managed-mode<gr_fields-Mode-Description>`.
   In the absence of this flag, the instance is initialized in the standard, "unmanaged" mode.

Chaining together flags with the bitwise-OR operation is straightforward.
This is illustrated in the :ref:`earlier code snippet <gr_fields-Create-Destroy-Section>`.

.. caution ::

   Caution is required if you want to use temporary variables to build up an argument one flag at a time.
   The following is a snippet that demonstrates an acceptable way to do it.

   .. code:: c

      uint64_t my_flag = 0;
      my_flag |= GR_FCREATE_MANAGED_FIELDDATA;
      my_flag |= GR_FCREATE_DUMMY;

   In particular you may get erroneous (or at least non-portable) behavior if you try to use a signed integer.


.. admonition:: Stablization Consideration

   There are a couple of considerations:

     1. Do we want to use bitwise flags?
        - I do think we should build in a mechanism, (whether it is bitwise flags or something else), to allow us to add more options to the initialization function over time.
        - While we currently just have one option right now, I could easily imagine a future where we want add flags to specify:

          * whether an instance holds single vs double precision values :ref:`at runtime <gr_fields-SupportingRuntimePrecision>`.

          * the memory space where data is located (e.g. when we add GPU support)

        - The main alternative to using bitwise flags is creating another opaque pointer type dedicated to configuring :c:type:`gr_fields` (it needs to be opaque to let us extend it in the future).

     2. Should we require users to specify a flag to denote they want Grackle to be initialized in standard :ref:`unmanaged mode<gr_fields-Mode-Description>`, rather than just saying that behavior is the default?
        (Such a flag might be called, :c:macro:`!GR_FCREATE_UNMANAGED_FIELDDATA`).


.. _gr_fields-Create-Destroy-FuncRef:

API Reference
^^^^^^^^^^^^^

.. c:function:: int grunstableFields_init_from_local(gr_fields** ptr, const chemistry_data* my_chemistry, const chemistry_data_storage* my_rates, uint64_t flags, int64_t nelements_per_field);

   Initialize an instance of :c:type:`gr_fields` from the types used to drive the :ref:`"local API." <local_functions>`.

   :param gr_fields** ptr: holds the created pointer
   :param const chemistry_data* my_chemistry: :ref:`fully configured <local_setup_data-storage>` run-time parameters
   :param const chemistry_data_storage* my_rates: chemistry and cooling rate data structure
   :param uint64_t flags: Bitwise ORed creation flags.
                          If no flags are used, this should be ``0``.
   :param int64_t nelements_per_field: When the ``GR_FCREATE_MANAGED_FIELDDATA`` flag is provided, this argument specifies the number of elements allocated per field. 
                                       Otherwise, this should be ``0``.
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

   .. note::

      The interface could be improved.
      While we currently just support a single flag, we could imagine a future where we use flags to specify the memory-space or a specialized implementation.
      But maybe we should just simplify the implementation for now?


.. c:function:: int grunstableFields_free(gr_fields* ptr);

   Deallocate an instance of :c:type:`gr_fields`

   :param gr_fields* ptr: the pointer to the instance that will be deleted
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

.. _gr_fields-Setting-Grid-Props:

Setting Grid Properties
+++++++++++++++++++++++

These functions are used to specify grid properties:

.. c:function:: int grunstableField_set_grid_layout(gr_fields* fields, int grid_rank, const int64_t* extent, const int64_t* start, const int64_t* stop);

   Specifies the data layout of each field.
   :param gr_fields* fields: the pointer to the instance that will be modified
   :param int grid_rank: the number of dimensions (it is also the length of each specified array)
   :param const int64_t* extent: specifies the extent along each axis
   :param const int64_t* start: specifies the first index valid along each each axis
   :param const int64_t* stop: specifies the first positive invalid index along each each axis
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

   There are 2 valid ways to use this function:

     1. The caller passes non-``NULL`` arguments for  ``extent``,  ``start``, and ``stop`` arguments
     2. The caller passes a non-``NULL`` arguments for  ``extent`` while passing ``NULL`` to both ``start`` and ``stop``.

   The second case will perform calculations on the entire grid (It is equivalent to passing an array of zeros for ``start`` and the same array to both ``extent`` and ``stop``)


   This function is analogous to manually setting the :c:var:`grid_rank`, :c:var:`grid_dimension`, :c:var:`grid_start`, and :c:var:`grid_end` members of :c:type:`grackle_field_data`.
   **However, there are 2 notable differences:**

     1. Whereas the :c:var:`grid_end` member inclusively tracks the maximum valid index along each axis, :c:func:`grunstableFields_set_grid_stop` specifies the minimum invalid positive index along each axis.

     2. This function copies values out of specified arrays (i.e. the caller does **NOT** need to keep these arrays around for the duration of the :c:type:`gr_fields` instance’s lifetime).

   .. admonition:: Stablization Consideration

      A reasonable alternative to this function might involve breaking this into separate “setter” functions for each quantity.
      The best fully consistent way to do that involves specifying the grid_rank during initialization of a :c:type:`gr_fields` and locking that quantity for the instance’s lifetime.
      That is a totally reasonable choice (you do sacrifice the ability to perform consistency checking between the values of ``extend``, ``start`` and ``stop``)

      All other ways require directly storing the user’s specified pointers (or storing a bunch of extra data and doing some extra work before every call calculation).

.. c:function:: int grunstableFields_set_grid_dx(gr_fields* fields, double dx);

   Set the cell-width of the field (in code units).
   We don't recommend setting this unless you need it.

   :param gr_fields* fields: the pointer to the instance that will be modified
   :param double dx: the width of each cell
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

.. _gr_fields-Query-Field-Names:

Querying available fields
+++++++++++++++++++++++++

These functions are used to query the fields that are made available by a :c:type:`gr_fields` instance.
For convenience, we associate each field name with an index.

.. warning::

   While the order of field names are fixed for a given a :c:type:`gr_fields` instance, it is allowed to change for other instances of :c:type:`gr_fields` (when the set of supported fields differ).
   It is also allowed to change for different compilation settings and in different Grackle versions.

.. c:function:: int grunstableFields_num_fields(const gr_fields* fields, int64_t* n_fields);

   Queries the number of fields that a given instance supports.

   :param const gr_fields* fields: pointer to the instance being queried
   :param int64_t* n_fields: pointer where the result will be stored
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

.. c:function:: int grunstableFields_name_by_idx(const gr_fields* fields, int64_t idx, const char** name);

   Queries the field name associated with the index (in the specified :c:type:`gr_fields` instance)

   :param const gr_fields* fields: pointer to the instance being queried
   :param int64_t n_fields: Transient index that identifies the field.
   :param const char** name: pointer where the result will be stored
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

.. c:function:: int grunstableFields_idx_by_name(const gr_fields* fields, const char* name, int64_t* idx);

   Queries the index associated with field name (in the specified :c:type:`gr_fields` instance).
   An index of ``-1`` indicates that the field name is not tracked by the specified instance.

   :param const gr_fields* fields: pointer to the instance being queried
   :param const char* name: the field name
   :param int64_t idx: Pointer where the field will be stored.
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful.
             The case where a field called ``name`` isn't tracked is still considered a success.

Accessing Field "Slots"
+++++++++++++++++++++++

In this section we describe the functions used to access/modify the "slots" of a given :c:type:`gr_fields` instance.

.. _gr_fields-SupportingRuntimePrecision:

.. admonition:: Stablization Consideration

   With the ongoing transcription to transcribe all of Grackle's internal Fortran, we are not far from a future where we could a simulation code's floating-point precision could be specified at runtime.
   Is this something we would be interested in supporting?

   If so, then we should probably replace each function in this section with 2 versions that are respectively suffixed by ``_f32`` and ``_f64``.
   If we wanted to get fancy, we could create function-like macro for each function with the ``_grflt`` suffix that dispatches to the appropriate version based on how ``gr_float`` is defined).

   In this scenario, we would also need to add :ref:`creation flags <gr_fields-CreationFlags>` to accomplish this (maybe called :c:macro:`!GR_FCREATE_F32` and :c:macro:`!GR_FCREATE_F64`?)

.. c:function:: int grunstableFields_store_by_name(gr_fields* fields, const char* name, gr_float* ptr);

   Store the specified pointer in the field-slot associated with the specified name.

   .. caution::

      It is an error to call this for a :c:type:`gr_fields` instance that was initialized in the :ref:`managed mode<gr_fields-Mode-Description>`.

   :param const gr_fields* fields: pointer to the instance being queried
   :param int64_t n_fields: Transient index that identifies the field.
   :param const char** name: pointer where the result will be stored
   :rtype: int
   :returns: ``GR_SUCCESS`` if successful

.. c:function:: gr_float** grunstableFields_get_slot_by_name(gr_fields* fields, const char* name)

   Return a pointer to the specified slot of the provided :c:type:`gr_fields` instance.

   .. warning::

      If the provided :c:type:`gr_fields` instance is configured in :ref:`managed mode<gr_fields-Mode-Description>`, then you are only allowed to modify the data in the buffer currently specified by the slot.
      It is illegal to change the referenced buffer that is referenced by the slot in this mode.
      (It is a little expensive to check this internally, so the implementation has to assume that external code "does the right thing")

   .. admonition:: Stablization Consideration

      The function exists for 2 purposes:

        1. It may be useful with a given :c:type:`gr_fields` instance, configured in the standard, :ref:`unmanaged mode<gr_fields-Mode-Description>`, is resued between multiple different calculations and each "slot's" contents must be overwritten before a calculation.
           In this case, it may be slightly faster (compared to repeatedly calling :c:func:`grunstableFields_store_by_name`) to record the pointers to the "slots" (fetched by this function) and use the recorded pointers to modify the "slots."

        2. In the :ref:`managed mode<gr_fields-Mode-Description>`, this function is used to fetch the managed buffers that the simulation code must copy data into or out of.

      Should we split this into 2 distinct functions so that it can't be misuesed in :ref:`managed mode<gr_fields-Mode-Description>`?

   .. note::

      This is the only function in this API that doesn't return a result code.
      The only alternative implementation to would involve accepting ``gr_float***`` as a parameter, which could be confusing. 

   :param const gr_fields* fields: pointer to the instance being queried
   :param const char* name: name of the associated field
   :rtype: gr_float**
   :returns: Pointer to the field-slot (the slot has a type :c:texpr:`gr_float*`.
             If unsuccesfuly, this is `NULL`.


Incomplete Areas of the Experimental API
----------------------------------------

In the near future, there are plans to introduce types called :c:type:`!gr_params` (a counter-part to :c:type:`chemistry_data`) and a type called :c:type:`!gr_solver`


.. rubric:: Footnotes

.. [#f1] A subtle, but key way that an API using Opaque Pointer pointers support a stable ABI is the fact that they perform memory allocations and deallocations internally.
         This allows the size (and alignment requirements) of the internal representation to change between versions of a library (there are other ways to do this, but need to take extreme care to ensure portability across languages).
         In contrast, an API using transparent types generally embeds the size and layout of the types into the client code, which is used for allocating and interacting with the datatype (and so ABI is broken if transparent type's representation changes at all).

.. [#f2] A common scenario is when a simulation code may be modified to allow usage of a new Grackle parameter/field for a particular simulation, all end-users of that simulation code are forced to update Grackle, even if it doesn't affect the user's simulation.
         Thus, simulation code developers may be incentivized to avoid using newer Grackle features to avoid inconveniencing the code's end-users.
         Avoiding this scenario is technically possible with transparent datatype if the simulation code's developer is very careful and uses specialized API functions (e.g. Enzo-E uses the :ref:`dynamic_api_functions` to avoid this situation while supporting new runtime parameters).
         However, a well-written API using opaque pointers can make it trivial to avoid this situation by default.
