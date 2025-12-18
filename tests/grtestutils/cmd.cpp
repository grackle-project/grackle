#include "cmd.hpp"

#include <stdio.h>  // FILE, fgets, fprintf, stderr, (popen/pclose on POSIX)

#include <cstdio>  // std::exit
#include <sstream>  // std::stringstream
#include <string>

// TODO: replace with Grackle's existing err machinery
[[noreturn]] static void err_(std::string msg="") {
  const char* ptr = (msg.empty()) ? "<unspecified>" : msg.c_str();
  fprintf(stderr, "ERROR: %s\n", ptr);
  std::exit(1);
}

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
  if (fp == nullptr) { err_("there was a problem launching the command"); }

  // if our reads from the pipe outpace the rate at which the command writes to
  // the pipe, fgets won't return until more data becomes available. If the
  // processess ends, fgets encounters EOF (causing fgets to return)
  while (fgets(temp_buf, TEMP_BUF_SIZE, fp) != nullptr) { buf << temp_buf; }
  int status = pclose(fp);
  if (status == -1) { err_("there was a problem closing the command"); }
  return {status, buf.str()};
}

#else

grtest::ProcessStatusAndStdout grtest::capture_status_and_output(
  const std::string& command
) {
  err_("not implemented on this platform");
}

#endif /* PLATFORM_GENERIC_UNIX */
