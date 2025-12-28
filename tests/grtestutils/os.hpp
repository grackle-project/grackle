// This defines a series of os-related utilities to help test grackle

#ifndef GRTEST_OS_HPP
#define GRTEST_OS_HPP

#include <memory>  // std::unique_ptr
#include <stdio.h>  // FILE*
#include <string>

namespace grtest {

/// Machinery to capture a stream (namely stderr)
///
/// This is similar in concept to pytest's capfd entity (I haven't looked under
/// the hood of that type, but I imagine they work in a similar manner)
///
/// @note
/// I'm not thrilled about the class's interface. But, it's good enough to
/// start.
class CaptureSink {
  FILE* redirected_stream_ = nullptr;
  int redirected_fd_ = -1;
  int backup_fd_ = -1;

  CaptureSink() = default;
  CaptureSink(const CaptureSink&) = delete;
  CaptureSink& operator=(const CaptureSink&) = delete;
  CaptureSink(CaptureSink&&) = delete;
  CaptureSink& operator=(CaptureSink&&) = delete;

public:

  ~CaptureSink() { end_capture(false); }
  /// create an instance that wraps stream
  static std::unique_ptr<CaptureSink> create(FILE* stream);
  /// end the capture and return the object
  std::string end_capture(bool ret_captured = true);
};

} // namespace grtest

#endif /* GRTEST_OS_HPP */
