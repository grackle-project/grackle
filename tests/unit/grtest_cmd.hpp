#ifndef GRTEST_CMD_H
#define GRTEST_CMD_H

#include <string>
#include <vector>

#include <sys/types.h> // pid_t

namespace grtest {

/// Creates a scratch dir. On destruction/close the dir & contents are removed
class ScratchDir{
  bool closed_ = true;
  std::string path_;

public:

  std::string get_path() const noexcept { return path_; }

  /// recursively deletes the scratch directory and everything in it
  /// (if this was previously called, then this does nothing)
  void close() noexcept;

  ~ScratchDir() { close(); }

  /// when specified the last 6 characters of suffix_template must be XXXXXX
  ScratchDir(std::string dir = "", std::string basename_template = "");

  ScratchDir(ScratchDir&&)=delete;  // could implement in future
  ScratchDir& operator=(ScratchDir&&)=delete;  // could implement in future
  ScratchDir(const ScratchDir&)=delete;  // incompatible with type
  ScratchDir& operator=(const ScratchDir&)=delete;  // incompatible with type

};


/// Encapsulates the result of launching a process
struct ProcessOutput {
  /// negative val indicates the process didn't terminate normally. Otherwise,
  /// it holds the process's exit-code
  int status;
  /// data captured from stdout
  std::string stdout_str;
  /// data captured from stderr
  std::string stderr_str;
};

/// Tracks component of a command
///
/// If we need to do anything more complex (e.g. more complex stream
/// redirection or ENV variable manipulation), I think we should shift to a
/// builder-pattern loosely inspired by Rust's std::process::Command type
struct Command {
  /// program to be executed
  std::string program;
  /// list of arguments
  std::vector<std::string> args;
  /// optionally holds a path to a file where stdout is redirected
  std::string stdout_path;
  /// optionally holds a path to a file where stderr is redirected
  std::string stderr_path;
};

/// Executes the specified command, waits for it to complete, & returns result
ProcessOutput process_output(const Command &cmd);

/// return a string holding summary of a command and its result
std::string summarize_cmd_and_rslt(const Command&, const ProcessOutput&);

} // grtest

#endif /* GRTEST_CMD_H */
