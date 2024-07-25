#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <regex>
#include <string>
#include <string_view>

#include "handle_gr_parameters.h"
#include "utils.h"

namespace { // stuff inside anonymous namespace is local to this file

// the white-space on either side of the equal-sign is ONLY meaningful if the
// user properly quotes the key-value pair
const std::regex kv_regex("^([a-zA-Z_0-9]+)[\\s]*=[\\s]*(.*)$");

const std::regex flt_regex("^" FLT_PATTERN "$");
const std::regex int_regex("^\\d+$");
const std::regex str_regex("^" STRING_PATTERN "$");


// the returned string is empty if there aren't any problems
template<typename Fn>
std::string for_each_char_(std::string_view s, Fn fn) {
  // we assume that s[0] == '"' and s[s.size()-1] == '"'
  int i = 1;
  int stop = s.size()-1;
  while (i < stop) {
    if (s[i] == '\\'){
      if ((i+1) == stop) return {"lone backslash character"};

      switch (s[i+1]) {
        case '"':  { fn('"'); break; }
        case '\\': { fn('\\'); break; }
        case 'n':  { fn('\n'); break; }
        case 't':  { fn('\t'); break; }
        default:
          return (std::string("invalid escape-sequence: \\") + s[i+1]);
      }
      i+=2;
    } else {
      fn(s[i]);
      i++;
    }
  }
  return {};
}

/// returns true for success
///
/// when this returns false and rslt is empty, then s simply didn't correspond
/// to a string. When this returns false and rslt isn't empty, then rslt
/// specifies an error.
bool parse_string_val_(std::string s, std::string& rslt){
  if ((s.size() > 2) && !std::regex_match(s, str_regex)) {
    rslt.clear();
    return false;
  }

  // count the number of characters in the string
  int count = 0;
  std::string err_msg = for_each_char_(s, [&count](char){ count++; });

  if (!err_msg.empty()) {
    rslt = err_msg;
    return false;
  } else {
    rslt.clear();
    rslt.reserve(count);
    for_each_char_(s, [&rslt](char ch){ rslt.push_back(ch); });
    return true;
  }
}

} // anonymous namespace


CliParamSpec::KVPair CliParamSpec::parse_param_(std::string_view token){

  // copy the token
  std::string token_str(token);
  //std::printf("parsing %s\n", token_str.c_str());

  CliParamSpec::KVPair out;
  std::smatch kv_match;
  if (!std::regex_match(token_str, kv_match, kv_regex)) {
    std::fprintf(stderr,
        "grackle-parameters should have the form `<key>=<value>`. The "
        "`%s` argument does not meet these expectations.", token_str.c_str());
    std::exit(1);
  }

  // match 0 is full string
  // extract the key-name
  std::string_view key(token.data() + kv_match.position(1), kv_match.length(1));
  
  // extract the substring corresponding to the value and then parse it
  std::string value_str = kv_match[2].str();

  value_variant value;
  if (std::regex_match(value_str, int_regex)) {
    int tmp_val = std::stoi(value_str);
    value = tmp_val;
  } else if (std::regex_match(value_str, flt_regex)) {
    double tmp_val = std::stod(value_str);
    value = tmp_val;
  } else {
    std::string tmp_val;
    bool is_str = parse_string_val_(value_str, tmp_val);

    if (is_str) {
      value = tmp_val;
    } else if (!tmp_val.empty()) {
      // in this case, tmp_val encodes an error message
      std::fprintf(stderr,
          "An error occured while parsing the right-hand-side of the `%s` "
          "grackle-parameter argument as a string: %s\n",
          token_str.c_str(), tmp_val.c_str());
      std::exit(1);
    } else {
      std::string key_str{key};
      std::fprintf(stderr,
          "The expression on the right-hand-side of the `%s`\n"
          "grackle-parameter argument is invalid.\n"
          "-> It must be an integer, floating-point or doubly-quoted str\n"
          "-> If you are using bash and the argument looked like\n"
          "       %s=\"%s\"\n"
          "   then there's a problem with quote-escaping. You may want to "
          "try passing:\n"
          "       '%s=\"%s\"'\n",
          token_str.c_str(), key_str.c_str(), value_str.c_str(),
          key_str.c_str(), value_str.c_str());
      std::exit(1);
    }
  }

  return {key, value};
}
