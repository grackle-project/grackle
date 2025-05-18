#ifndef GRTEST_CMD_H
#define GRTEST_CMD_H

#include <string>

namespace grtest {

/// Encapsulates the result of launching a process
struct ProcessStatusAndStdout {
  /// negative val indicates the process didn't terminate normally. Otherwise,
  /// it holds the process's exit-code
  int status;
  /// data captured from stdout
  std::string stdout_str;
};


/// Executes the specified command (it is passed to the Shell), waits for it to
/// complete, & returns the captured status and stdout.
///
/// @note
/// The child process's stderr is inherited from the parent
ProcessStatusAndStdout capture_status_and_output(const std::string& command);

} // grtest

#endif /* GRTEST_CMD_H */
