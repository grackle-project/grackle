## Overview

Files in this directory hold routines that have been transcribed from Fortran to C, but still are called within other Fortran subroutines.

## Calling Functions from Fortran

To be able to call the C functions from Fortran, we make use of Fortran 2003's ```ISO_C_BINDING`` module (which is backported on most Fortran compilers).

Currently the interfaces for the function are conditionally defined in header files (next to the declaration of the C prototype). To call these functions from a Fortran subroutine/function, you should add

  1. ``USE ISO_C_BINDING`` to the subroutine/function

  2. On a subsequent line of the subroutine/function, you should add ``#define FORTRAN_INTERFACE``

  3. After that you should then write an include statement for the header file. (this should probably come after an include statement for ``grackle_fortran_types.def``).

If you don't include the Fortran interface declaration, then the Fortran compiler has no way of knowing how to pass certain arguments into the function (whether it should be passed by value or by reference...)

This practice of putting the C and Fortran declarations of a function's interface into the same file is a little hacky and is probably not a great idea... We probably want to revisit this concept

## Fortran Shims

There are some issues with the datatype of the iteration mask and interoperability. At it's core, this issue boils down to the fact that the Fortran doesn't guarantee any compatibility between any C type and the LOGICAL type. The sole exception is C99's ``_Bool`` type and ``logical(C_BOOL)``.

While one long-term solution is to use ``logical(C_BOOL)`` throughout the codebase. However, this may be undesirable if we ever want to use C++ in the core codebase (which is probably necessary for supporting GPUs), since C++'s ``bool`` type is not [guaranteed to be compatible](https://stackoverflow.com/q/40020423) with ``_Bool``.

As a short-term solution, we introduce shims to convert the logical mask into a temporary integer mask.