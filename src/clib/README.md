# Basic Code Organization

This directory contains all of Grackle's C and Fortran source code files. It also contains all private header files.

The directory located at `../include` (i.e. the `include` directory at `src/include`) contains Grackle's public header files.

We review the differences between these headers further down this page.

# Include-Directive Conventions

<!-- I feel like this needs to be near the top of the page to increase the
     chance that people see this
     -> I'm giving this advice to make it easy for us to move things around
        later.
 -->

When including a public header file inside one of the source files contained by this directory, please directly specify the name of the header file and **NOT** a relative path.

Suppose we had a file called ``./my_source.c``. To include the ``grackle.h`` file:
- you should write ``include "grackle.h"``
- you should NOT write ``include "../include/grackle.h"``

(we use the `-I` flag to tell the compiler where to search for the public headers)

# More info and some conventions

**What is the difference between public and private headers?**

For the uninitiated:

- public header files get installed with Grackle.
  - These are the header files that are made available to downstream applications.
  - Declarations for all of the functions and types that are part of our public API [described here](https://grackle.readthedocs.io/en/latest/Reference.html) must be stored in the so that files.

- private header files are only used internal during the compilation process of Grackle. They include functions and types that are only used inside of Grackle. You should think of these as functions as private implementation details.

*NOTE: just because a function/type is in the public header does not mean it is a part of the public API (the functions/types that are part of the public API are explicitly listed in the documentation at the above link). For example, we reserve the write to alter the contents of the * `chemistry_data_storage` * struct and it's contained structs*.

Going forward, new private functions or types should generally be declared in the private headers.

**How do we name files?**
In general, all public header files should probably include ``grackle`` at the start of their names (to avoid name-conflicts during installation). There's no need for a private header file to do this (in fact, it's probably preferable if they don't do this to make it easier to distinguish whether its public or private. A public header must **NOT** share a name with a private header.
