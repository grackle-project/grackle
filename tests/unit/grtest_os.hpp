// This defines a series of os-related utilities to help test grackle

#ifndef GRTEST_OS_HPP
#define GRTEST_OS_HPP

#include <filesystem>
#include <map>
#include <memory>  // std::unique_ptr
#include <optional>
#include <stdio.h>  // FILE*
#include <string>

namespace grtest {

/// Represents a temporary directory (the directory is closed at construction)
///
/// @note
/// This class isn't only meant for writing test-cases and NOT intended for
/// general production (For general production, we need to come up with an
/// alternative formalism, where we don't invoke any operations that can fail
/// in the destructor).
class TempDir {
  std::filesystem::path path;
  bool is_open = false;

  TempDir() = default;  // force people to use the create factory method

public:
  // normal copy operations are forbidden
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;

  // move constructor and assignment
  TempDir(TempDir&& other) noexcept : TempDir() { swap(other); }
  TempDir& operator=(TempDir&& other) noexcept {
    this->swap(other);
    return *this;
  }

  /// swap contents with other
  void swap(TempDir& other) noexcept;

  /// indicates whether the file actually exists
  explicit operator bool() const noexcept { return is_open; }

  /// Factory function that create a scratch directory. The last 6 characters
  /// of the dirname will be randomly generated
  static TempDir create(const std::string& prefix = "");

  ~TempDir() { close(); }

  // recursively deletes the scratch directory and everything in it
  // (if this was previously called, then this does nothing)
  void close() noexcept;

  const std::filesystem::path& get_path() const noexcept { return path; }
};


/// Machinery to overwrite env variables (and revert the environemnt when the
/// instance leaves scope)
///
/// @note
/// This is intended to be within GTest tests. Only a single instance of this
/// type should exist at any given time (it functions extremely well as a
/// component of a fixture)
class EnvManager {
  std::map<std::string, std::optional<std::string>> orig_vals_;
  EnvManager() = default;
public:
  EnvManager(const EnvManager&) = delete;
  EnvManager& operator=(const EnvManager&) = delete;
  EnvManager(EnvManager&&) = delete;
  EnvManager& operator=(EnvManager&&) = delete;

  /// returns nullptr if we can't make an EnvManager on current platform
  static std::unique_ptr<EnvManager> create();
  ~EnvManager() noexcept;

  void override_envvar(const std::string& envvar, const std::string& val) {
    override_envvar(envvar, std::make_optional(val));
  }

  void override_envvar(const std::string& envvar,
                       const std::optional<std::string>& val);
};

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
