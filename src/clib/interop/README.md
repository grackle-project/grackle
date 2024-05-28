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

[^1]: (this may not be strictly necessary right now, but it will probably be necessary in the future).