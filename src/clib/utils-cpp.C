/// @file utils.C
/// @brief Implements various C++ utilities used in the routines transcribed
///        from Fortran

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "utils-cpp.hpp"

int eprintf(const char* format, ...) {
  std::va_list args;
  va_start(args, format);
  int out = vfprintf(stderr, format, args);
  va_end(args);
  return out;
}

/// helper function that prints an error message & aborts the program
[[noreturn]] void grackle::impl::Abort_With_Err_
    (const char* func_name, const char* file_name, int line_num, const char* msg, ...)
{
  const char* santized_func_name = (func_name == nullptr)
    ? "{unspecified}" : func_name;

  std::vector<char> msg_buf;
  if (msg == nullptr) {
    msg_buf = std::vector<char>(80);
    std::snprintf(msg_buf.data(), msg_buf.size(), "{nullptr encountered instead of error message}");
  } else {
    std::va_list args, args_copy;
    va_start(args, msg);
    va_copy(args_copy, args);

    std::size_t msg_len = std::vsnprintf(nullptr, 0, msg, args) + 1;
    va_end(args);

    msg_buf = std::vector<char>(msg_len);
    std::vsnprintf(msg_buf.data(), msg_len, msg, args);
    va_end(args_copy);
  }


  // now write the error and exit
  std::fprintf(stderr,
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "Error occurred in %s on line %d\n"
               "Function: %s\n"
               "Message: %s\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n",
               file_name, line_num, santized_func_name, msg_buf.data());
  std::fflush(stderr);  // may be unnecessary for stderr
  std::abort();
}
