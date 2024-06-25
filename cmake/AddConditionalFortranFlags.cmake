# this module primarily defines a function that "conditionally" applies fortran
# flags to a given target (that are compiler-specific).
# -> this is "conditional" in the sense that generator-expressions are used to
#    ensure that the flags are only applied to Fortran files
#
# the compiler will also provide some basic warnings about potential
# compiler-specific issues


# first, we define functions to print commonly-formatted compiler-specific
# warning messages

function(fc_warn_untested)
  # compared to a C compiler, it makes a lot of sense to warn about untested
  # fortran compilers because
  # - we adopt a fairly custom fortran dialect
  # - default name-mangling behavior of Fortran symbols can vary a lot between
  #   different compilers
  message(WARNING
    "${CMAKE_Fortran_COMPILER_ID} Fortran compiler is NOT well tested")
endfunction()

function(fc_warn_ninja)
  # this creates a warning about using the current compiler with the ninja build
  # generator (only gets shown if configured with the ninja build generator)

  # optionally specify a description of the fortran compiler
  if (ARGC EQUAL 0)
    set(compiler_descr "the ${CMAKE_Fortran_COMPILER_ID} Fortran compiler")
  elseif (ARGC EQUAL 1)
    set(compiler_descr "${ARGV0}")
  else()
    message(FATAL_ERROR "Too many args were provided")
  endif()

  # gfortran usually appears to work fine... not sure what goes wrong with the
  # other compilers... (is it ninja related? Or does CMake infer more about how
  # to use gfortran?)

  if (CMAKE_GENERATOR MATCHES "Ninja")
    # We use MATCHES to catch extra Ninja generators like Ninja Multi-Config,
    # or (recently deprecated) IDE generators that wrap Ninja

    message(WARNING
      "Compilation issues have been encountered while using ${compiler_descr} "
      "with the Ninja build-generator. If you encounter problems, consider "
      "using the Makefile generator instead.")
  endif()
endfunction()


# down below, we implement the main function provided by this module

function(add_conditional_fortran_flags target)
  # this function should ONLY introduce the minimal set of flags required to
  # succesfully compile the program
  # -> in other words, these flags only address formatting and name-mangling
  # -> it's considered an anti-pattern to hardcode other kinds of flags (it
  #    creates challenges for the end-user to specify desired flags)

  # report to user what we are doing: we follow convention and repeat the
  # description of what we are doing alongside the result
  set(_grackle_fc_flag_msg "Identify compiler-specific Fortran flags")
  message(STATUS "${_grackle_fc_flag_msg}")

  # do some (non-exhaustive) handling of the fortran-specific flags!
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "INTEL")
    fc_warn_ninja() # issues encountered with ninja while using cmake 3.16.3
    set(GRACKLE_FC_FLAG "") # no flags needed!

  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "^PGI|NVHPC$")
    # context: nvhpc is the new branding of pgi compilers
    # - the flags were inherited from old Machine-Makefile, but unclear whether
    #   they will succesfully compile via cmake
    # - stuff does compile with nvfortran (but tests have not been run)
    fc_warn_untested()
    fc_warn_ninja() # issues encountered with nvhpc and ninja while using cmake
                    # version 3.27.7
    set(GRACKLE_FC_FLAG "-Mnosecond_underscore;-Mextend")

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "XL") # IBM compiler
    fc_warn_untested() # things compile, but tests haven't been run
    fc_warn_ninja() # issues encountered with ninja with cmake version 3.27.7

    # it's worrisome to pass -qhalt=s (it suppress certain classes of errors),
    # but it's currently necessary for cool1d_cloudy_old_tables_g.F
    set(GRACKLE_FC_FLAG "-qfixed=132;-qhalt=s")

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(GRACKLE_FC_FLAG "-fno-second-underscore;-ffixed-line-length-132")

  else()
    message(STATUS
      "No explicit logic for ${CMAKE_Fortran_COMPILER_ID} Fortran compiler, testing GNU flags"
    )
    fc_warn_ninja("certain Fortran compilers")

    set(GRACKLE_FC_FLAG "")
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag(-fno-second-underscore no_second_underscore)
    check_fortran_compiler_flag(-ffixed-line-length-132 fixed_line_length_132)
    list(APPEND GRACKLE_FC_FLAG
      "$<$<BOOL:${no_second_underscore}>:-fno-second-underscore>"
      "$<$<BOOL:${fixed_line_length_132}>:-ffixed-line-length-132>"
    )

  endif()

  # report the results of this function
  if (GRACKLE_FC_FLAG STREQUAL "")
    message(STATUS "${_grackle_fc_flag_msg} - none needed")
  else()
    message(STATUS "${_grackle_fc_flag_msg} - ${GRACKLE_FC_FLAG}")
    target_compile_options(${target} PRIVATE
      "$<$<COMPILE_LANGUAGE:Fortran>:${GRACKLE_FC_FLAG}>"
    )
  endif()
endfunction()
