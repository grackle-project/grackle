// See LICENSE file for license and copyright information

/// @file status_reporting.h
/// @brief Declares the tools used for result-reporting
///
/// Purpose
/// =======
/// Define machinery to let us standardize our internal logic for
/// status-reporting. This will:
/// - reduce boilerplate code
/// - allow us to provide more context (e.g. filename/lineno/function name)
/// - allow us to more easily change how we communicate this information with
///   the downstream application in the future (if we decide that is something
///   we ultimately want to do -- we motivate this in the next section)
///
/// Current Status
/// ==============
/// In general, status-reporting is an area where Grackle could really improve.
/// At the time of writing, our error/status reporting leaves a lot to be
/// desired. Essentially, it consists of:
/// - returning a binary code GR_SUCCESS or GR_FAIL
/// - printing some extra output information:
///   - we use fprintf(stderr, ...) in the C routines when errors arise
///   - we print some context information in the Fortran
///
/// Ideally, we would adopt a solution that provides the downstream application
/// (i.e. the simulation-code or Gracklepy) control over whether/how this info
/// is consumed/recorded.
///
/// Motivating Perspective
/// ======================
/// Must of us are developers of simulation-codes, but status reporting in the
/// context of a library is different. In particular, considerations for a
/// numerical library are very different from other libraries.
///
/// For this discussion, we ignore warnings.
///
/// The statuses that we communicate primarily fall under 2: (i) unrecoverable
/// internal errors and (ii) function results.
///
/// It is also important to realize that there are certain errors we just can't
/// check without crippling Grackle's performance. Documentation needs to be
/// exceptionally clear in these cases.
///
/// I. Unrecoverable Internal Errors
/// --------------------------------
/// The idea is that if Grackle detects that it is a state where it seems like
/// it may produce bad/unreliable results/behavior, we print an error message
/// and immediately abort the program.
///
/// Generally, this is the result of a sanity check or a broken invariant.
///
/// As a rule of thumb:
/// - if the sanity check is a direct result of the application's action, we
///   should prefer to gracefully communicate this to the application
/// - always consider, should this error crash an interactive Jupyter session
///   using Gracklepy?
/// - loudly failing is ALWAYS better than silently failing. Since these errors
///   are quick and easy to add, these can be a quick temporary solution to
///   address a problematic code path (and an indication that we need to more
///   gracefully handle the error mode in the future)
///
/// Function Results
/// ----------------
/// There are generally 4 kinds of results/statuses a grackle-like library
/// might want to report. 2 are currently relevant. 2 more may become relevant
/// in the future. These result/statuses may include
/// - the API function successfully completed an operation
/// - (may become relevant with GPUs) a function, which can be configured to be
///   either synchronous or asynchronous (@ compile-time or runtime), wants to
///   report that it is being invoked in an "asynchronous mode" and that it
///   succesfully launched a GPU kernel. The premise is that the function would
///   report the successful-completion status if the function had been
///   configured in its "synchronous mode," and the whole operation had been a
///   success. (Honestly, we may choose to simply denote success in both cases)
/// - (may become relevant with GPUs) a temporary error occurred that could be
///   overcome by trying again. For example, a GPU had too much pending work
///   (honestly, we may want to configure Grackle so it knows to try again).
/// - A generic error occurred that requires human-intervention
///
/// It is useful to highlight categories of these generic errors (the last of
/// the above kinds). There is some overlap here. Categories may include:
/// 1. resource errors like out-of-memory. While we can't really deal with OOM
///    errors on CPUs (due to "memory overcommitment") it may be worth
///    supporting on GPUs.
/// 2. invalid parameter values
///    - depending on the application and parameter this may result from
///      end-user or application
///    - there's a lot of value to providing lots of info about these errors
///      (to make debugging much more pleasant)
/// 3. file-system errors when reading the datafile. This could happen if the
///    specified path isn't valid (e.g. it doesn't point to a valid hdf5 file,
///    there is a permissions issue, etc.). Or it could happen do to a system
///    error (a network filesystem is down or too many files are open)
/// 4. An obviously invalid argument is detected (e.g. passing a nullptr). This
///    is clearly an application error.
/// 5. A broken precondition indicating that the application has made an error
///    using the API
/// 6. Generic calculation error (e.g. exceeding the iteration limit or
///    detecting a NaN). This could arise because the application provided
///    nonsensical field values. It could also arise simply because the current
///    logic doesn't handle a particular case very well. An application might
///    want to know about this class of errors so that it could log extra
///    information (e.g. the current simulation cycle and location) so that the
///    error can be reproduced later.
///
/// Hypothetically, it could also be cool to allow more recoverable
/// error-handling when a calculation error only affects a small subset of
/// field values (e.g. the code might overwrite those zones with values from
/// neighboring cells). From a practical perspective, this may not be a good
/// use of time.
///
/// How this can be improved
/// ========================
/// An issue will be created that describes various strategies for improved
/// reporting. A common strategy for dealing with simple errors involves exit
/// codes...

#ifndef STATUS_REPORTING_H
#define STATUS_REPORTING_H

#ifdef __cplusplus
extern "C" {
#endif

// define attributes to annotate functionsA:
// 1. ERRFMT_ATTR_(fmt_pos) tells the compiler argument number `fmt_pos`
//    expects a printf-style format-string. Compilers supporting this will
//    know to check consistency between arg-types and the string @ compile-time
// 2. suppress warnings that may arise about not returning in control-flows
//    where an internal error aborts the program with NORETURN_ATTR_
// 3. compiler raises warnings when the value returned by a function annotated
//    with NODISCARD_ATTR_ isn't used. We use this for a function where this is
//    almost certainly indicative of programming error

#if 0
  // this branch will be used once we transition to compiling all non-fortran
  // files with C++17. (We use the official syntax for specifying attributes)
  // - starting in C++17, implementations know to ignore attributes that they
  //   don't recognize
  #define ERRFMT_ATTR_(fmt_pos) [[gnu::format(__printf__, fmt_pos, fmt_pos+1)]]
  #define NORETURN_ATTR_ [[noreturn]]
  #define NODISCARD_ATTR_ [[nodiscard]]

#elif defined(__GNUG__)
  #define ERRFMT_ATTR_(fmt_pos)                                               \
    __attribute__((format(__printf__, fmt_pos, fmt_pos+1)))
  #define NORETURN_ATTR_ __attribute__((noreturn))
  // unlike [[nodiscard]], warnings associated with the following might not be
  // suppressed by casting the result to (void). (this is ok for this file)
  #define NODISCARD_ATTR_ __attribute__((warn_unused_result))

#else
  #define NORETURN_ATTR_ /* ... */
  #define ERRFMT_ATTR_(fmt_pos) /* ... */
  #define NODISCARD_ATTR_ /* ... */

#endif

// ---------------------------------------

/// @def      __GRIMPL_PRETTY_FUNC__
/// @brief    a magic contant like __LINE__ or __FILE__ used to specify the
///           name of the current function
///
/// In more detail:
/// - The C99 and C++11 standards ensures __func__ is provided on all
///   platforms, but it only provides limited information (just the name of the
///   function).
/// - note that __func__ is technically not a macro. It's a static constant
///   string implicitly defined by the compiler within each function definition
/// - Where available, we prefer to use compiler-specific features that provide
///   more information about the function (like the scope of the function, the
///   the function signature, any template specialization, etc.).
#ifdef __GNUG__
  #define __GRIMPL_PRETTY_FUNC__ __PRETTY_FUNCTION__
#else
  #define __GRIMPL_PRETTY_FUNC__ __func__
#endif

struct grimpl_source_location_{
  const char* file;
  int lineno;
  const char* fn_name;
};

/// This is a helper function used to help implement __GRIMPL_SRCLOC__
///
/// @note
/// static is required to use inline with C
static inline struct grimpl_source_location_ get_src_location_(
  const char* file, int lineno, const char* fn_name
) {
  struct grimpl_source_location_ out;
  out.file = file;
  out.lineno = lineno;
  out.fn_name = fn_name;
  return out;
}

/// @def __GRIMPL_SRCLOC__
/// @brief Roughly equivalent to __FILE__, __LINE__, etc. But, it gathers the
///        info for us in a very concise manner
#define __GRIMPL_SRCLOC__                                                   \
  get_src_location_(__FILE__, __LINE__, __GRIMPL_PRETTY_FUNC__)

/// helper function that helps implement GR_INTERNAL_ERROR and
ERRFMT_ATTR_(2) NORETURN_ATTR_ void grimpl_abort_with_internal_err_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
);

/// @def GR_INTERNAL_ERROR
/// @brief function-like macro that handles a (lethal) error message
///
/// This macro should be treated as a function with the signature:
///
///   [[noreturn]] void GR_INTERNAL_ERROR(const char* fmt, ...);
///
/// The ``fmt`` arg is a printf-style format argument specifying the error
/// message. The remaining args arguments are used to format error message
///
/// @note
/// the ``fmt`` string is part of the variadic args so that there is always
/// at least 1 variadic argument (even in cases when ``msg`` doesn't format
/// any arguments). There is no portable way around this until C++ 20.
#define GR_INTERNAL_ERROR(...)                                            \
  { grimpl_abort_with_internal_err_(__GRIMPL_SRCLOC__, __VA_ARGS__); }
// we define GRIMPL_ERROR to avoid merge conflicts. The plan is to remove it in
// the future (after avoiding merge conflicts)
#define GRIMPL_ERROR(...)                                                 \
  { grimpl_abort_with_internal_err_(__GRIMPL_SRCLOC__, __VA_ARGS__); }

/// @def GR_INTERNAL_REQUIRE
/// @brief implements functionality analogous to the assert() macro
///
/// if the condition is false, print an error-message (with printf
/// formatting) & abort the program.
///
/// This macro should be treated as a function with the signature:
///
///   void GR_INTERNAL_REQUIRE(bool cond, const char* fmt, ...);
///
/// - The 1st arg is a boolean condition. When true, this does nothing
/// - The 2nd arg is printf-style format argument specifying the error message
/// - The remaining args arguments are used to format error message
///
/// @note
/// The behavior is independent of the ``NDEBUG`` macro
#define GR_INTERNAL_REQUIRE(cond, ...)                                     \
  {  if (!(cond))                                                              \
      { grimpl_abort_with_internal_err_(__GRIMPL_SRCLOC__, __VA_ARGS__); } }
// we define GRIMPL_REQUIRE to avoid merge conflicts. The plan is to remove it
// in the future (after avoiding merge conflicts)
#define GRIMPL_REQUIRE(cond, ...)                                             \
  {  if (!(cond))                                                             \
      { grimpl_abort_with_internal_err_(__GRIMPL_SRCLOC__, __VA_ARGS__); } }

// helper function
ERRFMT_ATTR_(2) NODISCARD_ATTR_ int grimpl_print_and_return_err_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
);

/// @def GrPrintAndReturnErr
/// @brief prints the error message & returns the appropriate status-value
///
/// This macro should be treated as a function with the signature:
///
///   [[nodiscard]] int GrPrintAndReturnErr(const char* fmt, ...);
///
/// The ``fmt`` arg is a printf-style format argument specifying the error
/// message. The remaining args arguments are used to format error message
///
/// This is intended to be used in scenarios like
/// ```{.cpp}
///   return GrPrintAndReturnErr("My err. %g is a bad val", val);
/// ```
/// or
/// ```{.cpp}
///   int ec = GrPrintAndReturnErr("My err. %g is a bad val", val);
///   // do some cleanup...
///   return ec;
/// ```
///
/// The compiler will issue a warning if you don't do anything with the exit
/// code. If you don't care about the exit-code, you should use the
/// `GrPrintErrMsg` function-like macro instead
///
/// @note
/// the ``fmt`` string is part of the variadic args so that there is always
/// at least 1 variadic argument (even in cases when ``fmt`` doesn't format
/// any arguments). There is no portable way around this until C++20.
///
/// @note
/// This macro has been designed so that we have a uniform/easily grep-able
/// interface that we can easily replace in the future if/when we improve error
/// reporting
#define GrPrintAndReturnErr(...)                                             \
  grimpl_print_and_return_err_(__GRIMPL_SRCLOC__, __VA_ARGS__);


// helper function
ERRFMT_ATTR_(2) void grimpl_print_err_msg_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
);

/// @def GrPrintErrMsg
/// @brief prints the appropriate error message.
///
/// > [!important]
/// > In any situation where you ultimately return an error code, you should
/// > prefer to use GrPrintAndReturnErr. This macro is for the subset of cases
/// > where you need to do something different.
///
/// This macro should be treated as a function with the signature:
///
///   void GrPrintErrMsg(const char* fmt, ...);
///
/// The ``fmt`` arg is a printf-style format argument specifying the error
/// message. The remaining args arguments are used to format error message
#define GrPrintErrMsg(...)                              \
  grimpl_print_err_msg_(__GRIMPL_SRCLOC__, __VA_ARGS__);


// undefine the attributes so we avoid leaking them
// ------------------------------------------------
#undef ERRFMT_ATTR_
#undef NORETURN_ATTR_
#undef NODISCARD_ATTR_

// I don't think we can undef __GRIMPL_PRETTY_FUNC__ without causing issues

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* STATUS_REPORTING */
