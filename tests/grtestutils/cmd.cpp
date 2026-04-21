#include "cmd.hpp"

// internal Grackle routine:
#include "status_reporting.h"  // GR_INTERNAL_ERROR, GR_INTERNAL_REQUIRE

#include <stdio.h>  // FILE, fgets, fprintf, stderr, (popen/pclose on POSIX)

#include <sstream>  // std::stringstream
#include <string>

#ifdef PLATFORM_GENERIC_UNIX

#define TEMP_BUF_SIZE 128

grtest::ProcessStatusAndStdout grtest::capture_status_and_output(
  const std::string& command
) {
  // note: there are negligible portability benefits to using the `system` 
  //   function instead of `popen` since we would need the `WIFEXITED` &
  //   `WEXITSTATUS` macros (provided by POSIX) to interpret the value returned
  //   by `system`
  std::stringstream buf;
  char temp_buf[TEMP_BUF_SIZE];

  // fp represents a buffered pipe to the standard output of the command.
  FILE* fp = popen(command.data(), "r");
  GR_INTERNAL_REQUIRE(fp != nullptr, "there's a problem launching the command");

  // if our reads from the pipe outpace the rate at which the command writes to
  // the pipe, fgets won't return until more data becomes available. If the
  // processess ends, fgets encounters EOF (causing fgets to return)
  while (fgets(temp_buf, TEMP_BUF_SIZE, fp) != nullptr) { buf << temp_buf; }
  int status = pclose(fp);
  GR_INTERNAL_REQUIRE(status != -1, "there was a problem closing the command");
  return {status, buf.str()};
}

#else

grtest::ProcessStatusAndStdout grtest::capture_status_and_output(
  const std::string& command
) {
  GR_INTERNAL_ERROR("not implemented on this platform");
}

#endif /* PLATFORM_GENERIC_UNIX */
