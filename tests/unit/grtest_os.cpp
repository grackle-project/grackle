#include "grtest_os.hpp"
#include "grtest_utils.hpp"
#include "status_reporting.h"

// the following 2 headers are the c versions of the headers (rather than the
// C++ versions) since it seems more likely that posix-specific functions are
// provided in this version
#include <stdio.h> // fflush, tmpfile, fread, fileno
#include <stdlib.h> // mkstemp, getenv

#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <random>
#include <utility>  // std::swap

void grtest::TempDir::swap(TempDir& other) noexcept {
  std::swap(this->path, other.path);
  std::swap(this->is_open, other.is_open);
}

grtest::TempDir grtest::TempDir::create(const std::string& prefix) {
  // NOTE: at some point we should give more thought to whether we want to use
  // exceptions. Grackle-proper doesn't use exceptions. But, it can be tricky
  // to use std::filesystem without exceptions (as I understand it, even
  // though some functions use std::error_code you can still encounter
  // exceptions)
  std::filesystem::path p;
  if (prefix.empty()) {
    p = std::filesystem::temp_directory_path() / "my-tmpdir-";
  } else {
    p = std::filesystem::absolute(prefix);
    char last = prefix[prefix.size() - 1];
    if (last == '/' || last == std::filesystem::path::preferred_separator) {
      p = p / "my-tmpdir-";
    }
    GR_INTERNAL_REQUIRE(std::filesystem::is_directory(p.parent_path()),
                        "prefix doesn't specify a path in an existing dir");
  }

  std::string prefix_str = p.string();
  std::minstd_rand rng(std::random_device{}());
  std::string full_path;

  // I think python does something similar
  for (int i = 0; i < 100; i++) {
    int random_int =
        static_cast<int>(999999 * grtest::random::uniform_dist_transform(rng));
    full_path = prefix_str + std::to_string(random_int);

    std::error_code ec;
    std::filesystem::create_directory(full_path, ec);
    if (!ec) {
      break;
    }
    full_path.clear();
  }

  GR_INTERNAL_REQUIRE(!full_path.empty(), "could not create directory");
  grtest::TempDir out;
  out.path = full_path;
  out.is_open = true;
  return out;
}

void grtest::TempDir::close() noexcept {
  // we are just going to abort if an exception comes up
  if (is_open) {
    // technically the next function can raise an exception about memory
    // allocations, but noexcept will mean we just abort the program
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
    if (ec) {
      std::string path_string = path;
      GR_INTERNAL_ERROR(
          "there was a problem recursively removing contents of %s",
          path_string.c_str());
    }
    is_open = false;
  }
}

#ifndef PLATFORM_GENERIC_UNIX
// provide dummy implementations that will never work
#define unsetenv(var) /* ... */
#define setenv(var, val, overwrite) /* ... */
std::unique_ptr<grtest::EnvManager> grtest::EnvManager::create() {
  return nullptr;
}
#else
std::unique_ptr<grtest::EnvManager> grtest::EnvManager::create() {
  return std::unique_ptr<grtest::EnvManager>(new grtest::EnvManager());
}
#endif



void grtest::EnvManager::override_envvar(const std::string& envvar,
                                         const std::optional<std::string>& val) {
  const char* cur_val = getenv(envvar.c_str());

  // record the original value so we can restore it at the end of the test
  if (orig_vals_.find(envvar) != orig_vals_.end()) {
    // cur_val isn't the original value (the original is already recorded)
  } else {
    orig_vals_[envvar] = (cur_val == nullptr)
                             ? std::nullopt
                             : std::make_optional(std::string(cur_val));
  }

  if (cur_val == nullptr && !val.has_value()) {
    // nothing needs to be done
  } else if (!val.has_value()) {
    unsetenv(envvar.c_str());
  } else {
    setenv(envvar.c_str(), val.value().c_str(), /*overwrite:*/ 1);
  }
}

grtest::EnvManager::~EnvManager() noexcept{
  for (const auto& kv_pair : orig_vals_) {
    const std::string& envvar = kv_pair.first;
    bool currently_has_value = getenv(envvar.c_str()) != nullptr;
    bool originally_has_value = kv_pair.second.has_value();
    if (currently_has_value && !originally_has_value) {
      unsetenv(envvar.c_str());
    } else if (!originally_has_value) {  // && !currently_has_value
      // do nothing!
    } else {
      setenv(envvar.c_str(), kv_pair.second.value().c_str(),
             /*overwrite:*/ 1);
    }
  }
}


#ifndef PLATFORM_GENERIC_UNIX

// provide dummy implementations that will never work
std::unique_ptr<grtest::CaptureSink> grtest::CaptureSink::create(FILE*) {
  return nullptr;
}
std::string grtest::CaptureSink::end_capture(bool) { return {}; }

#else
#include <unistd.h> // dup, dup2, lseek, SEEK_SET, SEEK_CUR, close

// the logic we are invoking here is VERY standard. Lots and lots of examples
// of it exist online

/// closes the current file descriptor associated with fd and then opens a new
/// file descriptor (reusing the same integer) that points to a temporary file.
/// The temporary file will be deleted once all file descriptors to it are
/// closed.
///
/// @returns
/// true and false indicates success and failure
static bool redirect_fd_to_tempfile_(int fd) {
  // create a temporary file and get a descriptor to it
  char* prefix = getenv("TMPDIR");
  std::string path = (prefix != nullptr) ? prefix : "/tmp";
  path += "/grackle-test-XXXXXX";
  int temporary_fd = mkstemp(path.data());
  if (temporary_fd == -1) { return false; }

  // remove the temporary file's name from the filesystem (so it is deleted
  // once there are no more file descriptors refer to it)
  unlink(path.data()); // <-- nothing we can reasonably do if this fails

  // in a single atomic transaction, the next line:
  // 1. closes the file descriptor tracked by fd (i.e. the integer value will
  //    no longer is associated with any file and can be used to create a new
  //    file descriptor)
  // 2. allocates a brand new file descriptor, which both reuses the file
  //    descriptor number from fd AND refers to the file currently
  //    referenced by temporary_fd
  if (dup2(temporary_fd, fd) == -1) {
    close(temporary_fd);
    return false;
  }

  // finally, close the temporary_fd
  close(temporary_fd); // <-- nothing we can reasonably do if this fails

  return true;
}

// if we want to avoid heap-allocations, we could use std::optional
std::unique_ptr<grtest::CaptureSink> grtest::CaptureSink::create(FILE* stream) {
  if (stream == nullptr) { return nullptr; }

  // 1. flush the stream (this is critical for the general case)!
  if ( fflush(stream) != 0 ) { return nullptr; }

  // 2. get the file descriptor that we want to redirect
  int underlying_fd = fileno(stream);
  if (underlying_fd == -1) { return nullptr; }

  // 3. create a new file descriptor that references the same file as
  //    `underlying_fd` (a "backup" we use later for restoring properties)
  int backup_fd = dup(underlying_fd);
  if (backup_fd == -1) { return nullptr; }

  // 4. closes the current file descriptor associated with underlying_fd and
  //    then opens a new file descriptor (reusing the same integer) that points
  //    to a temporary file. Upon success, any writes to stream will write to
  //    the temporary file -- rather than the original destination
  if (!redirect_fd_to_tempfile_(underlying_fd)) {
    close(backup_fd);
    return nullptr;
  }

  // 5. prepare the container
  std::unique_ptr<grtest::CaptureSink> out{new grtest::CaptureSink};
  out->redirected_stream_ = stream;
  out->redirected_fd_ = underlying_fd;
  out->backup_fd_ = backup_fd;
  return out;
}

// since the class is primary intended to capture stderr, we eagerly abort if
// things go wrong (why even go on if testing is crippled?)
std::string grtest::CaptureSink::end_capture(bool ret_captured) {
  std::string out;
  if (this->redirected_fd_ < 0) { return out; }

  // flushing is critical for the general case!
  if (fflush(this->redirected_stream_) != 0) { abort(); };

  // determine the number of captured bytes we wish to read
  off_t nbytes = 0;
  if (ret_captured) {
    nbytes = lseek(this->redirected_fd_, 0, SEEK_END);
    if (nbytes == -1) { abort(); }
    // restore position in file to the very start
    if (lseek(this->redirected_fd_, 0, SEEK_SET) == -1) { abort(); }
  }

  // actually read in bytes
  if (nbytes > 0) {
    out = std::string(nbytes, ' '); // <- allocate the necessary space
    char* buf = out.data();
    if (read(this->redirected_fd_, buf, nbytes) != nbytes) { abort(); }
  }

  // in a single atomic transaction, the next line:
  // 1. closes the file descriptor tracked by redirected_fd_ (i.e. the integer
  //    value no longer is associated with any file and can be when creating
  //    a new file descriptor)
  // 2. allocates a brand new file descriptor, which both reuses the file
  //    descriptor number from redirected_fd_ AND refers to the file currently
  //    referenced by backup_fd_
  // At the time of writing, redirected_fd_ was the last descriptor refering to
  // a temporary file, whose name was already deleted from the file system.
  // Thus, the temporary file no longer exists after this operation.
  if (dup2(this->backup_fd_, this->redirected_fd_) == -1) { abort(); }

  // we no longer need to keep backup_fd_
  close(this->backup_fd_);

  // overwrite variables so we know not to execute this logic again
  redirected_stream_ = NULL;
  this->redirected_fd_ = -1;
  this->backup_fd_ = -1;

  return out;
}

#endif
