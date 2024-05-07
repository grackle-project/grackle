## Overview

Files in this directory are used to define "interoperable functions." These are C functions that may need to be called from Fortran subroutines (usually these C functions define functionality that was previously transcribed from Fortran subroutines).

Currently, the C function-prototypes for these "interoperable functions" are defined in [interop_funcs.h](./interop_funcs.h) while the corresponding Fortran interfaces are declared within [interop_funcs.def](./interop_funcs.def).

## Calling These Functions from Fortran

To be able to call the C functions from Fortran, we make use of Fortran 2003's ```ISO_C_BINDING`` module (which is backported on most Fortran compilers).

When you're writing a Fortran subroutine (or function) that needs to call one or more of these "interoperable functions", it's critical that your function/routine properly include the interface. To accomplish this, you **MUST** make sure your subroutine includes lines with the following directives:

1. It must contain a line that reads ``USE ISO_C_BINDING`` (this usually comes just before the line that reads ``implicit NONE``)

2. After the line that reads ``implicit NONE``, you should insert your include directives. There are 2 relevant include directives (order matters!):

   1. First, you should have an include-directive for ``grackle_fortran_types.def`` file[^1]

   2. Next, you should have an include-directive for the ``interop_funcs.def`` file

> [!CAUTION]
> If you try to call an interoperable function without following the above instructions, the problems will arise. Specifically, the Fortran compiler may not make proper guesses about name-mangling **AND** it has no way of knowing how to pass certain arguments into the function (whether it should be passed by value or by reference...)

For the sake of concreteness, the following codeblock illustrates the general structure that a Fortran subroutine, defined in `src/clib/`, would need to have in order for it to call one of these interoperable functions:

```fortran
      subroutine my_subroutine( arg1, arg2, arg3 )

!  < docstring sumarizing purpose and describing args ... >

      USE ISO_C_BINDING
      implicit NONE
#include "grackle_fortran_types.def"
#include "interop/interop_funcs.def"

!  < declare arg types and then declare local variables (with their types) ... >

!  < do actual work work... >

      return
      end
```


## Fortran Shims

There are some issues with the datatype of the iteration mask and interoperability. At it's core, this issue boils down to the fact that the Fortran doesn't guarantee any compatibility between its default ``LOGICAL`` type and any C type.

As a short-term solution, we introduce shims to convert the logical mask into a temporary integer mask. This is somewhat undesirable given that memory is allocated each time we call a shim.

There are essentially 2 long-term solutions:

1. the simplest solution is to replace all occurences of the ``LOGICAL`` type in the Fortran code with ``logical(C_BOOL)``.

   - In more detail, ``C_BOOL`` is a "KIND type parameters"  provided by the ``ISO_C_BINDING`` module. Using ``logical(C_BOOL)``, ensures that the datatype is compatible with the ``_Bool`` datatype introduced in C99.[^2]

   - However, this may produce problems if we want to compile the code with a C++ compiler (this is essentially required to support GPUs). In more detail, C++'s ``bool`` type is [NOT guaranteed to be compatible](https://stackoverflow.com/q/40020423) with C's ``_Bool`` type (in practice, this probably isn't much of an issue).

2. The alternative, is to define some consistent integer datatype that is used in to hold the mask values in the Fortran and C layers. This approach might look like the following:

   - inserting something like ``#define gr_mask_int int32_t`` into ``grackle_macros.h`` and ``#define GR_MASK_INT integer*4`` into ``grackle_fortran_types.def``

   - replacing all usage of ``LOGICAL`` with ``GR_MASK_INT`` (all boolean operations would need to be replaced with comparisons against 0)

[^1]: (this may not be strictly necessary right now, but it will probably be necessary in the future).

[^2]: This is footnote is provided for the sake of clarity since people often get confused about C's boolean type. The 1999 C standard introduced ``_Bool`` as a standard type for holding booleans (i.e. any compliant C compiler always defines this type, even when no headers are included). The type is only capable of holding ``1`` or ``0``. The 1999 standandard also introduced the [<stdbool.h>](https://en.cppreference.com/w/c/types/boolean) component to the "standard library" to define useful macros. These macros include ``bool``, which expands to ``_Bool``, as well as the ``true`` and ``false`` macros, which expand to the integer constants ``1`` and ``0``. (As an aside, things are a little different when a compiler the 2023 standard, but that's not relevant to us)