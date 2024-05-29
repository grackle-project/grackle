# this module defines a function that "conditionally" applies fortran flags to
# a given target (that are compiler-specific).
# -> this is "conditional" in the sense that generator-expressions are used to
#    ensure that the flags are only applied to Fortran files

function(add_conditional_fortran_flags target)
  # this function should ONLY introduce the minimal set of flags required to
  # succesfully compile the program (in other words, they should only address
  # formatting and name-mangling)

  # report to user what we are doing: we follow convention and repeat the
  # description of what we are doing alongside the result
  set(_grackle_fc_flag_msg "Identify compiler-specific Fortran flags")
  message(STATUS "${_grackle_fc_flag_msg}")

  # do some (non-exhaustive) handling of the fortran-specific flags!
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "INTEL")
    set(GRACKLE_FC_FLAG "") # no flags needed!
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(GRACKLE_FC_FLAG "-Mnosecond_underscore;-Mextend")
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(GRACKLE_FC_FLAG "-fno-second-underscore;-ffixed-line-length-132")
  else()
    message(STATUS
      "No explicit logic for ${CMAKE_Fortran_COMPILER_ID} Fortran compiler, testing GNU flags"
    )
    set(GRACKLE_FC_FLAG "")
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag(-fno-second-underscore no_second_underscore)
    check_fortran_compiler_flag(-ffixed-line-length-132 fixed_line_length_132)
    list(APPEND GRACKLE_FC_FLAG
      "$<$<BOOL:${no_second_underscore}>:-fno-second-underscore>"
      "$<$<BOOL:${fixed_line_length_132}>:-ffixed-line-length-132>"
    )
  endif()

  if (GRACKLE_FC_FLAG STREQUAL "")
    message(STATUS "${_grackle_fc_flag_msg} - none needed")
  else()
    message(STATUS "${_grackle_fc_flag_msg} - ${GRACKLE_FC_FLAG}")
    target_compile_options(${target} PRIVATE
      "$<$<COMPILE_LANGUAGE:Fortran>:${GRACKLE_FC_FLAG}>"
    )
  endif()

endfunction()
