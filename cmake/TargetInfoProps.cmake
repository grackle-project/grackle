# to make our package easily consumable and provide a uniform interface whether
# a downstream cmake process directly builds our product as a part of their
# own build or they find an existing installation with a "CMake Package
# Config" file.
# -> a recommended practice to achieve that is to convey relevant information
#    through target properties. (The alternative is to use global variables,
#    but that is now discouraged)
# -> but this means that we need to maintain 2 lists of target properties (one
#    where declare/build things AND the other where we compile things)
# -> the idea is to provide a wrapper function around set_target_properties
#    that records this information, cache the list of properties, and provide
#    a second function to get the properties later (to use it to configure the
#    "CMake Package Config" file)
# -> while we're at it, it probably makes sense to support reporting this 
#    information in pkg-config files for non-cmake builds
#    -> in that case, we must provide it as global variables
#    -> if we're doing this, it makes sense to distinguish between BOOL and
#       STRING dtypes (since cmake has an expansive definition of a BOOL dtype)



# this function is just like set_target_properties. It has the signature
#
# set_target_typed_info_properties(<targets> ...
#                                  (STRING_PROPERTIES|BOOL_PROPERTIES)
#                                  <prop1> <value1>
#                                  [<prop2> <value2>] ...)
function(set_target_typed_info_properties)
  set(err_prefix "set_target_typed_info_properties called with")
 
  if(ARGC LESS_EQUAL 2)
    message(FATAL_ERROR "${err_prefix} incorect number of args") 
  endif()

  math(EXPR loop_stop_inclusive "${ARGC} - 1")
  foreach(i RANGE 1 ${loop_stop_inclusive})
    if("${ARGV${i}}" MATCHES "^((BOOL)|(STRING))_PROPERTIES$")
      math(EXPR properties_start "${i} + 1")
      set(dtype "${CMAKE_MATCH_1}")
      break()
    endif()
  endforeach()

  if (NOT DEFINED properties_start)
    message(FATAL_ERROR "{err_prefix} no (STRING)|(BOOL)_PROPERTIES specifier")
  elseif(properties_start GREATER 2)
    message(FATAL_ERROR "set_target_properties only supports 1 target arg")
  endif()
 
  # confirm that remaing args can be organized into name-value pairs
  math(EXPR is_odd "((${loop_stop_inclusive}+1) - ${properties_start}) % 2")
  if (is_odd)
    message(FATAL_ERROR "${arg_num_err_msg} incorrect number of args")
  endif()

  set(target "${ARGV0}")

  # iterate over all name-value pairs
  foreach(i RANGE ${properties_start} ${loop_stop_inclusive} 2)
    math(EXPR ip1 "${i} + 1")

    if (${dtype} STREQUAL "STRING")
      set(val "${ARGV${ip1}}")
    elseif(${ARGV${ip1}})
      set(val 1) # sanitized TRUE boolean
    else()
      set(val 0) # sanitized FALSE boolean
    endif()

    set_target_properties(${target} PROPERTIES "${ARGV${i}}" ${val})
    list(APPEND property_list "${ARGV${i}}")
  endforeach()

  # now we save the list of information in a special (global) cache variable
  set(varName _TargetInfoProps__${target}__info_props)
  # build a combined list if the global variable already defined
  set(property_list ${${varName}} ${property_list})
  list(REMOVE_DUPLICATES property_list)
  set(${varName} ${property_list} CACHE INTERNAL "list of info props" FORCE)
endfunction()

# returns info properties in format that is embedded into output-config files.
# valid export_format choices include:
# - PKG_CONFIG (in this case the string declares module variables)
# - CMAKE_CONFIG (in this case the string is put into set_target_properties)
function(get_info_properties_export_str target export_format out_varname)

  set(varName _TargetInfoProps__${target}__info_props)
  if(NOT DEFINED CACHE{${varName}})
    message(FATAL_ERROR "No registered info properties for ${target}")
  elseif("${export_format}" STREQUAL "PKG_CONFIG")
    set(separator "=")
  elseif("${export_format}" STREQUAL "CMAKE_CONFIG")
    set(separator " ")
  else()
    message(FATAL_ERROR "export_format must be PKG_CONFIG or CMAKE_CONFIG")
  endif()

  set(output "")
  foreach(property_name IN LISTS "${varName}")
    get_target_property(val ${target} ${property_name})
    string(APPEND output "${property_name}${separator}${val}\n")
  endforeach()
  set(${out_varname} "${output}" PARENT_SCOPE)

endfunction()
