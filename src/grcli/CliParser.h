/***********************************************************************
/
/ Declare class used to represent the CLI (command-line-interface) parser.
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef GRCLI_CLIPARSER_H
#define GRCLI_CLIPARSER_H

#include <cstdlib>
#include <cstdio>
#include <string_view>
#include <string>

#include "utils.h"


// declare a couple common utilities
void print_help(const char* bin_name);

inline void try_parse_help(std::string_view arg,
                           const char* bin_name) {
  if ((arg == "-h") || (arg == "--help")) {
    print_help(bin_name);
    std::exit(0);
  }
}

[[noreturn]] inline void err_unrecognized_arg(const char* arg) {
  std::fprintf(stderr, "unrecognized argument: %s\n",  arg);
  std::exit(1);
}

[[noreturn]] inline void err_unrecognized_arg(std::string_view arg) {
  std::string arg_str(arg);
  err_unrecognized_arg(arg_str.c_str());
}


/// Encapsulates the command line interface
///
/// @note
/// I have a lot of ideas for improvements
class CliParser {
  char * const * next_;
  char * const * end_;
  char* bin_name_;

public:
  CliParser() = delete;

  CliParser(int argc, char** argv) : next_(argv+1), end_(argv+argc)
  {
    GRCLI_REQUIRE(argc > 0,
      "Something is wrong, argc is always expected to be positive");
    bin_name_ = argv[0];
  }

  char * bin_name() const { return bin_name_; }

  bool has_next() const { 
    bool out = next_ < end_;
    if (out && (std::string_view(*next_) == "--")) {
      std::fprintf(stderr,
          "We don't respect the \"--\" cli argument since we don't "
          "currently support positional arguments.");
      std::exit(1);
    }
    return out;
  } 

  const char * peek() const { return has_next() ? *next_ : nullptr; }

  const char * next() { return has_next() ? *(next_++) : nullptr; }

  /// just like next, but instead of returning the next c-string, we return
  /// the pointer to the next c-string stored in argv
  char * const * next_argv_ptr() { return has_next() ? next_++ : nullptr; }

};




#endif /* GRCLI_CLIPARSER_H */
