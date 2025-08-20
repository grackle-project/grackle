// This existence of this header file is considered an implementation detail
// - it's an error for external project to directly include this file. Any
//   external project directly include this file can/will start encountering
//   problems @ an arbitrary point in the future
// - the only reason we don't currently abort with an error is so we can
//   gracefully warn about deprecated headers without introducing lots of
//   complex preprocessor logic


#ifndef GRACKLE_MISC_H
#define GRACKLE_MISC_H

/// @def GRIMPL_COMPTIME_WARNING(MSG)
/// @brief Internal macro used to produces a compile-time warning with the
///     specified message (it should be a string literal). This exists because
///     `#warning` is a non-standard extension (before C23 and C++23). On less
///     popular compilers, this is a no-op
///
/// @note
/// The use of GRIMPL_DO_PRAGMA is based on gcc documentation
/// https://gcc.gnu.org/onlinedocs/cpp/Pragmas.html

#define GRIMPL_DO_PRAGMA(x) _Pragma(#x)
#ifdef __GNUC__
  // Modern versions of mainstream compilers that masquerade as gcc (like
  // clang or intel), understand the pragma used to implement this macro.
  // - Any compiler that don't understand it will treat it as a no-op (as
  //   mandated by the C and C++ standards).
  #define GRIMPL_COMPTIME_WARNING(MSG) GRIMPL_DO_PRAGMA(GCC warning MSG)
#else
  #define GRIMPL_COMPTIME_WARNING(MSG) /* no-op */
#endif

#endif /* GRACKLE_MISC_H */
