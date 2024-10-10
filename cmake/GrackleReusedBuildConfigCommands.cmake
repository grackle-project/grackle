# ========
# Overview
# ========
#
# This is a cmake module that defines commands that are used in 2 cases:
# 1. as part of the build logic
# 2. as part of the configuration logic
#
# Since commands have global scope, these commands (functions or macros)
# defined in this file will always be visible to CMake projects that consume
# Grackle (regardless of how Grackle is consumed). Only a subset of commands
# are intended to be used directly by consumers. These commands are explicitly
# documented in the documentations
#
# IMPORTANT
# ---------
# All code in this file must be executable by the minimum supported cmake
# version described in `GrackleConfig.cmake.in` (at this time of writing, that
# is version 3.3)


# Commands Motivated by the grdata tool
# -------------------------------------
# Here, we turn our attention to commands originally created with the intention
# of supporting the grdata cli tool (the commands are named generically enough
# that they can be used to support other things in the future, too).
#
# Some background:
# - the grdata program can be very useful for downstream codes for testing
#   purposes (and for the hypothetical scenario where we support downloading
#   precompiled copies of Grackle). But it's a little weird in a CMake context
#   - we want people to think of it essentially as a command line program (and
#     be indifferent to the fact that it is implemented as an executable python
#     executable python script behind the scenes)
#   - since we want it to be available to downstream CMake applications whether
#     people embed grackle as part of their build or use find_package, we
#     probably shouldn't declare it as an IMPORTED executable. If we just cared
#     about supporting the latter case, then this would be absolutely fine.
#     But, IMPORTED executables aren't really intended for the former case
#     (especially since we are producing the file with configure_file). While
#     it could work out in the former case, I worry about unforseen side
#     effects.
# - Prior to CMake 3.0, we might have made this information available in a
#   package variable (like GRACKLE_GRDATA_TOOL_PATH), but the Professional
#   CMake book, written by Craig Scott (one of the primary CMake developer),
#   makes it clear we should prefer to avoid package-variables in order to make
#   the project as consumable as possible.
#   - This advice is offered in multiple versions of the book, but you can
#     specifically find it in section 40.4 of the 18th edition
#   - essentially, using package variables produces complication when people
#     want to use cmake 3.24 or newer (specifically related to the integration
#     between find_package and FetchContent).
#   - Following the book's suggestion, we instead make this information
#     available through the grackle_get_info command

function(_grackleprivate_get_propertyprefix outVar) # private helper function
  set("${outVar}" "_PrivateGrackleGlobalPropDoNotUse_" PARENT_SCOPE)
endfunction()

function(_grackleprivate_get_info_setup grdata_tool_path)
  # this is a private command to setup for the grackle_get_info command
  # ARGS:
  #   grdata_tool_path: specifies the path to the grdata tool

  # - the whole point here is to make the grackle_get_info without reading any
  #   information related to grdata_tool_path from a global/package variable.
  #   (It might be possible to define grackle_get_info so that it does read
  #   this info from global/package variables and is portable across all
  #   contexts, but that very much isn't clear).
  # - This leaves us with a few alternatives (th
  #   1. Write a file (maybe to ${PROJECT_BINARY_DIR}/__grackle_get_info.cmake)
  #      that defines the grackle_get_info function (while substituting the
  #      grdata_tool_path information) and then include that file in order to
  #      declare the function (we would need to make this into a macro)
  #   2. Store the grdata_tool_path information as an attribute on a target
  #      a) we could do this on an existing target, like Grackle::Grackle (in
  #         this case, we'd name the property so its clear that it's not part
  #         of the public interface)
  #      b) we could create a new interface target (like Grackle::_PRIVATE) for
  #         this explicit purpose
  # - It really doesn't which choice we pick since this is an implementation
  #   detail we could change at any time. For simplicity, we go with 2a.

  set(_target_name "Grackle::Grackle")
  get_target_property(_aliased "${_target_name}" ALIASED_TARGET)
  if(_aliased)
    set(_target_name "${_aliased}")
  endif()

  _grackleprivate_get_propertyprefix( "_prefix" )

  set_property(TARGET "${_target_name}"
    PROPERTY "${_prefix}GRACKLE_GRDATA_TOOL_PATH" "${grdata_tool_path}"
  )
endfunction()


# a public function for accessing extra grackle information
#
# An empty string is returned if the property is not known
#
# ARGUMENTS
# ---------
# key: the name of the grackle information to be accessed
# outVar: the name of the variable where the accessed info is stored
function(grackle_get_info key outVar)
  _grackleprivate_get_propertyprefix( "_prefix" )
  get_target_property(val Grackle::Grackle "${_prefix}${key}")
  set("${outVar}" "${val}" PARENT_SCOPE)
endfunction()
