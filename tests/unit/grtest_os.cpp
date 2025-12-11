#include "grtest_os.hpp"

// the following 2 headers are the c versions of the headers (rather than the
// C++ versions) since it seems more likely that posix-specific functions are
// provided in this version
#include <stdio.h> // fflush, tmpfile, fread, fileno
#include <stdlib.h> // mkstemp, getenv

#include <memory>
#include <string>

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

  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
  std::unique_ptr<grtest::CaptureSink> out(new grtest::CaptureSink);
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
