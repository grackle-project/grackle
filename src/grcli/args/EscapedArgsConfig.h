#ifndef ESCAPED_ARGS_CONFIG_H
#define ESCAPED_ARGS_CONFIG_H

#include <cstdlib>
#include <cstdio>
#include <string_view>
#include <vector>

#include "../args/CliParser.h"
#include "../utils.h"

namespace args {


/// This struct is used to help with parsing a batch of escaped arguments
struct EscapedArgsSpec {
  char short_opt_character;
  std::string_view group_arg_start;
  std::string_view group_arg_sentinel;

  template<typename Fn>
  bool try_parse(const char* leading_arg, CliParser& parser, Fn fn) {

    char short_opt[3] = {'-', short_opt_character, '\0'};

    if (starts_with(leading_arg, short_opt)) {
      const char * escaped = leading_arg + 2;
      if (leading_arg[2] == '\0') {
        escaped = (parser.has_next()) ? parser.next() : nullptr;
      }

      if (escaped == nullptr) {
        std::fprintf(stderr, "the %s arg is missing an option\n", leading_arg);
        std::exit(1);
      }
      fn(escaped);
      return true;

    } else if (leading_arg == this->group_arg_start) {
      int num_encountered = 0;
      const char * elem = parser.has_next() ? parser.next() : nullptr;
      while ((elem != nullptr) && (elem != this->group_arg_sentinel)) {
        fn(elem);
        num_encountered++;
        elem = parser.has_next() ? parser.next() : nullptr;
      }

      if (num_encountered == 0) {
        std::fprintf(stderr,
                     "the \"%s\" flag doesn't start a group of parameters\n",
                     leading_arg);
        std::exit(1);
      }
      return true;

    }
    return false;
  }
};


/// class used to collect the arguments that are forwarded to a backend
class BackendConfig {

  /// this will be forwarded to the backend
  std::vector<char *> arg_vec_;

public:

  /// this performs a copy of arg_vec_
  ///
  /// @note
  /// This is NOT a deep copy. Modifying the contents of the elements will
  /// affect the value in various locations
  std::vector<char *> args_copy() const { return arg_vec_; }

  /// Try to parse the argument
  bool try_parse(const char* leading_arg, CliParser& parser) { 
    EscapedArgsSpec arg_spec{'B', "--backend-start", "--backend-stop"};

    auto fn = [this](const char* arg) {
      this->arg_vec_.push_back((char*)arg); // this cast is bad...
    };
    return arg_spec.try_parse(leading_arg, parser, fn);
  }

};

} // namespace args

#endif /* ESCAPED_ARGS_CONFIG_H */
