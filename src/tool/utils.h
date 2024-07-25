
#ifndef TOOL_UTILS_H
#define TOOL_UTILS_H

#include <string_view>

// general-stuff that gets reused a lot

// this could be faster!
inline bool starts_with(std::string_view base, std::string_view substr) {
  return base.substr(0,substr.size()) == substr;
}

// there isn't an equivalent to python's (?:...) group (where the parentheses
// are non-capturing)
#define FLT_PATTERN "[-+]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][-+]?\\d+)?"

// regex pattern for a double-quote enclosed string that allows escaped quotes
// https://stackoverflow.com/a/5696141
#define STRING_PATTERN R"REGEX("[^"\\]*(?:\\.[^"\\]*)*")REGEX"


// the rest of this file is designed to help us format errors

namespace grackle::internal {

/// helper function that prints an error message & aborts the program
[[noreturn]] void Abort_With_Err_(const char* func_name, const char* file_name, int line_num, const char* msg, ...);

} // namespace grackle::internal

/// @def      __GRCLI_FUNC__
/// @brief    a magic contant like __LINE__ or __FILE__ used to specify the name
///           of the current function
///
/// In more detail:
/// - The C++11 standard ensures __func__ is provided on all platforms, but it
///   only provides limited information (just the name of the function).
/// - note that __func__ is technically not a macro. It's a static constant
///   string implicitly defined by the compiler within each function definition
/// - Where available, we prefer to use compiler-specific features that provide
///   more information about the function (like the scope of the function, the
///   the function signature, any template specialization, etc.).
#ifdef __GNUG__
  #define __GRCLI_PRETTY_FUNC__ __PRETTY_FUNCTION__
#else
  #define __GRCLI_PRETTY_FUNC__ __func__
#endif


/// @def GRCLI_ERROR
/// @brief function-like macro that handles a (lethal) error message
///
/// This macro should be treated as a function with the signature:
///
///   [[noreturn]] void GRCLI_ERROR(const char* fmt, ...);
///
/// The ``msg`` arg is printf-style format argument specifying the error
/// message. The remaining args arguments are used to format error message
///
/// @note
/// the ``msg`` string is part of the variadic args so that there is always
/// at least 1 variadic argument (even in cases when ``msg`` doesn't format
/// any arguments). There is no portable way around this until C++ 20.
///
/// @note
/// the ``msg`` string is part of the variadic args so that there is always
/// at least 1 variadic argument (even in cases when ``msg`` doesn't format
/// any arguments). There is no portable way around this until C++ 20.
#define GRCLI_ERROR(...)                                                      \
  { grackle::internal::Abort_With_Err_                                        \
      (__GRCLI_PRETTY_FUNC__, __FILE__, __LINE__, __VA_ARGS__); }


//----------------------------------------------------------------------
/// @def GRCLI_REQUIRE
/// @brief implements functionality analogous to the assert() macro
///
/// if the condition is false, print an error-message (with printf
/// formatting) & abort the program.
///
/// This macro should be treated as a function with the signature:
///
///   void GRCLI_REQUIRE(bool cond, const char* fmt, ...);
///
/// - The 1st arg is a boolean condition. When true, this does nothing
/// - The 2nd arg is printf-style format argument specifying the error message
/// - The remaining args arguments are used to format error message
///
/// @note
/// The behavior is independent of the ``NDEBUG`` macro
#define GRCLI_REQUIRE(cond, ...)                                              \
  {  if (!(cond))                                                             \
      { grackle::internal::Abort_With_Err_                                    \
         (__GRCLI_PRETTY_FUNC__, __FILE__, __LINE__, __VA_ARGS__); } }


#endif /* TOOL_UTILS_H */
