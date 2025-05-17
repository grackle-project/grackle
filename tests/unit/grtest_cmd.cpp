#include "grtest_cmd.hpp"

#include <cstdio>  // std::fprintf, stderr, std::fopen, std::fread, ...

#include <filesystem>
#include <optional>
#include <sstream> // std::stringstream
#include <stdlib.h> // exit, system, (mkdtemp on linux)
#include <string>
#include <vector>


// TODO: replace with Grackle's existing err machinery
[[noreturn]] static void err_(std::string msg="") {
  const char* ptr = (msg.size() == 0) ? "<unspecified>" : msg.c_str();
  std::fprintf(stderr, "ERROR: %s\n", ptr);
  exit(1);
}

#ifdef PLATFORM_GENERIC_UNIX

#include <unistd.h> // mkdtemp on some systems (e.g. MacOS)
#include <sys/wait.h> // WIFEXITED, WEXITSTATUS

static bool my_mkdtemp(std::string& path_template)
{ return (mkdtemp(path_template.data()) != nullptr); }

/// negative val indicates the process didn't terminate normally. Otherwise,
/// it holds the process's exit-code
static int my_system(const std::string& command) {
  int rv = system(command.c_str());
  if (WIFEXITED(rv)) { return WEXITSTATUS(rv); }
  return -1;
}

#else

static bool my_mkdtemp(std::string& path_template)
{ err_("not implemented on this platform"); }

static int my_system(const std::string& command)
{ err_("not implemented on this platform"); }

#endif /* PLATFORM_GENERIC_UNIX */


grtest::ScratchDir::ScratchDir(std::string dir, std::string suffix_template)
  : closed_(false), path_()
{
  if (dir.size() == 0) {
    path_ = std::string(std::filesystem::temp_directory_path());
  } else {
    path_ = dir;
  }
  path_.push_back('/');
  if (suffix_template.size() == 0) {
    path_.append("scratch-XXXXXX");
  } else {
    path_.append(suffix_template);
  }

  // modify path_ & create temporary directory
  if (!my_mkdtemp(path_)) { err_("mkdtemp(\"" + path_ + "\")"); }
}

void grtest::ScratchDir::close() noexcept {
  if (!closed_){
    try {
      std::filesystem::remove_all(path_);
    } catch (...){
      err_(std::string("problem recursively removing contents of") + path_);
    }
    closed_ = true;
  }
}


// load file contents into a string
static std::string load_to_string(std::string path) {
  std::FILE* f = fopen(path.c_str(), "r");
  if (!f) { return ""; };
  std::fseek(f, 0, SEEK_END);  // set file-position to the end of the file
  const std::size_t nchars = std::ftell(f);
  std::fseek(f, 0, SEEK_SET);  // set file-position to the start of the file

  std::string out = std::string(nchars, ' ');
  std::fread(out.data(), sizeof(char), nchars, f);

  std::fclose(f);
  return out;
}


grtest::ProcessOutput grtest::process_output(const grtest::Command &cmd) {
  // This is a very hacky implementation. A more robust implementation would
  // use file descriptors and posix_spawn rather than defering to the shell.
  //
  // But, this is good enough for running tests

  grtest::ProcessOutput out;
  if (cmd.program.size() == 0) { err_("No program was specified"); }

  // figure out how to handle stdout & stderr (e.g. do we capture them?)
  bool cap_out = cmd.stdout_path.size() == 0;
  bool cap_err = cmd.stderr_path.size() == 0;
  std::optional<ScratchDir> tmpdir{};
  if (cap_out || cap_err) { tmpdir.emplace(); }
  std::string o_path, e_path;
  o_path = (cap_out) ? (tmpdir->get_path() + "/out") : cmd.stdout_path;
  e_path = (cap_err) ? (tmpdir->get_path() + "/err") : cmd.stderr_path;

  // build up the command string (the allocation are very inefficient)
  std::string command = cmd.program;
  for (const std::string& arg : cmd.args) {
    command.push_back(' ');
    command.append(arg);
  }
  command.append(" 1>");
  command.append(o_path);
  command.append(" 2>");
  command.append(e_path);

  // execute command (and load any captured output)
  out.status = my_system(command);
  if (cap_out) { out.stdout_str = load_to_string(o_path); }
  if (cap_err) { out.stderr_str = load_to_string(e_path); }

  return out;
}

std::string grtest::summarize_cmd_and_rslt(const grtest::Command& cmd,
                                           const grtest::ProcessOutput& rslt)
{
  // this is semi-redundant with vector-formatting logic in a pending PR
  auto fmt_vec = [](const std::vector<std::string>& vec) -> std::string
  {
    std::stringstream buf;
    for (std::size_t i = 0; i < vec.size(); i++){
      if (i != 0){ buf << ", "; }
      buf << '"' << vec[i] << '"';
    }
    return buf.str();
  };

  std::string stdout_descr = (cmd.stdout_path.size() == 0)
    ? "\n" + rslt.stdout_str : " directed to " + cmd.stdout_path;
  std::string stderr_descr = (cmd.stderr_path.size() == 0)
    ? "\n" + rslt.stderr_str : " directed to " + cmd.stderr_path;

  std::stringstream buf;
  buf << "  program: " << cmd.program << '\n'
      << "  args:\n"
      << "    " << fmt_vec(cmd.args) << '\n'
      << "  status: " << rslt.status << '\n'
      << "  stdout: " << stdout_descr << '\n'
      << "  stderr: " << stdout_descr << '\n';
  return buf.str();

}
