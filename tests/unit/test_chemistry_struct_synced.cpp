#include "grtest_cmd.hpp"

#include <grackle.h>

#include <gtest/gtest.h>

#include <cstdio>
#include <map>
#include <optional>
#include <regex>
#include <string>
#include <utility>

// maps the names of types to a set of member names with that type
using MemberTypeNameMap = std::map<std::string, std::string>;

#if defined(stringify) || defined(stringify_helper)
  #error "stringify or stringify_helper was already defined"
#elif !defined(READER_PATH) || !defined(XML_PATH)
  #error "either the READER_PATH or XML_PATH macro isn't defined"
#else
  #define stringify(s) stringify_helper(s)
  #define stringify_helper(s) #s
  static const char* GLOBAL_reader_path = stringify(READER_PATH);
  static const char* GLOBAL_xml_path = stringify(XML_PATH);
#endif


struct MemberTypeNameMapResult {
  std::string err;
  MemberTypeNameMap map;
};


/// try to parse a string that holds a file where each line has the format
/// "<member_name>, <member_type>" and create a MemberTypeNameMap
///
/// @note
/// although std::regex is known to be **VERY** slow, it is convenient
static MemberTypeNameMapResult try_parse_(const std::string& s) {
  MemberTypeNameMap mapping;
  std::regex myregex("^([^,]+), *([^,]+) *$");
  std::size_t len = s.size();
  std::size_t pos = 0;
  while (pos < len) {
    std::size_t lineend = s.find('\n', pos);
    if (pos+1 >= lineend || lineend == s.npos) { break; }

    std::smatch match;
    if (std::regex_match(s.begin()+pos, s.begin()+lineend, match, myregex)) {
      std::string member_name = match[1].str();
      std::string member_type = match[2].str();
      if (member_type == "char*") { member_type = "string"; }
      mapping[member_name] = member_type;
    } else {
      MemberTypeNameMapResult out;
      out.err = "problem converting the reader's output to a map!";
      return out;
    }
    pos = lineend+1;
  }

  return MemberTypeNameMapResult{"", mapping};
}

/// Reads the previously generated xml-file (produced by the build-system with
/// the castxml tool) to determine the data-members of a struct. Returns a
/// std::map that maps the names of data-members to their type. The result of
/// this function is used for checking the correctness of functionallity that
/// effectively implements "manual reflection" of the struct's members.
///
/// @returns a non-empty error-message to denote **any** errors
///
/// In more detail, this function has 2 phases:
/// 1. we check if a valid file was written by build-system
/// 2. a python script is invoked to extract relevant information from the xml
///    file (it uses python's builtin xml-parsing library)
/// 3. we read in the extracted information and put it into the output format
///
/// @note
/// In the future, we'll be able to replace this with a simplified alternative:
/// - For example, we could use 1 of the following lightweight librarie to
///   directly infer members of aggregates (this includes C structs,
///   https://en.cppreference.com/w/cpp/language/aggregate_initialization) at
///   compile-time (using compiler-extensions)
///   - https://github.com/boostorg/pfr  (this doesn't need the rest of boost)
///   - https://github.com/stephenberry/glaze
/// - Alteratively, C++26 is on track to make "reflection" of fields a standard
///   C++ feature.
MemberTypeNameMapResult query_struct_members(std::string struct_name)
{
  // step 0: basic setup
  // ===================
  // -> this is the path to file created by the build-system
  std::string xml_path = GLOBAL_xml_path;
  // -> this is the path to the executable script to parse the file
  std::string xml_reader = GLOBAL_reader_path;

  // step 1: check that xml exists
  // =============================
  std::FILE* f = std::fopen(xml_path.c_str(), "r");
  if (f == NULL) {
    return MemberTypeNameMapResult{
      "The build-system failed to create the xml-file. (THIS SHOULDN'T HAPPEN)",
      {}
    };
  } else if (std::fgetc(f) == EOF) {
    std::string msg(
      "The build-system wrote an empty xml-file, which usually means that "
      "castxml wasn't installed. If you install castxml, you need to instruct "
      "CMake to rebuild of the target containing this test (to trigger a rebuild "
      "the XML file). It's usually easiest to clean things up and build again."
    );
    return MemberTypeNameMapResult{msg, {}};
  }
  std::fclose(f);

  // step 2: parse the xml
  // =====================
  grtest::Command xmlread_cmd;
  xmlread_cmd.program = xml_reader;
  xmlread_cmd.args = {"--input", xml_path, "--type", struct_name};
  xmlread_cmd.stdout_path = "";  // capture stdout
  xmlread_cmd.stderr_path = "";  // capture stderr
  grtest::ProcessOutput xmlread_rslt = grtest::process_output(xmlread_cmd);
  if (xmlread_rslt.status != 0) {
    std::string tmp = grtest::summarize_cmd_and_rslt(xmlread_cmd, xmlread_rslt);
    std::fprintf(stderr, "unexpected outcome:\n%s\n", tmp.c_str());
    return MemberTypeNameMapResult{"Something weird happened", {}};
  } else if (xmlread_rslt.stderr_str.size() != 0) {
    std::fprintf(stderr, "%s\n", xmlread_rslt.stderr_str.c_str());
  }

  // step 3: create the output object
  // ================================
  return try_parse_(xmlread_rslt.stdout_str);
}

/// an empty string indicates that the mappings are consistent
testing::AssertionResult HasExpectedMembers (
  std::string struct_name,
  const MemberTypeNameMap& inferred_map,
  const MemberTypeNameMap& ref_map
) {
  // use a range-based for loop with de-structuring
  for (const auto& [member_name, inferred_type]: inferred_map) {
    MemberTypeNameMap::const_iterator search = ref_map.find(member_name);
    if (search == ref_map.end()) {
      return ::testing::AssertionFailure()
        << "the \"manual-reflection machinery\" incorrectly indicates the "
        << struct_name << " type has a member called " << member_name
        << " (of type " << inferred_type << "), but no such type exists";
    } else if (inferred_type != search->second) {
      return ::testing::AssertionFailure()
        << "the \"manual-reflection machinery\" incorrectly indicates the "
        << struct_name << "::" << member_name << " data-member has type "
        << inferred_type << ". Actual type: " << search->second;
    }
  }

  if (inferred_map.size() != ref_map.size()) {
    for (const auto& [member_name, ref_type]: ref_map) {
      MemberTypeNameMap::const_iterator search = inferred_map.find(member_name);
      if (search == inferred_map.end()) {
        return ::testing::AssertionFailure()
          << "the \"manual-reflection machinery\" doesn't denote the existence "
          << "of the " << struct_name << "::" << member_name << "data-member.";
      }
    }
  }

  return ::testing::AssertionSuccess();
}


#define EXPECT_MEMBERS(member_map, struct_name)                               \
  /* use determine the actual members */                                      \
  MemberTypeNameMapResult ref_TEMP_ = query_struct_members( struct_name );    \
  if (ref_TEMP_.err.size() > 0 ) {                                            \
    GTEST_SKIP() << "Issue inferring members of " << struct_name << ". "      \
                 << ref_TEMP_.err;                                            \
  }                                                                           \
                                                                              \
  GTEST_EXPECT_TRUE(HasExpectedMembers(struct_name, member_map, ref_TEMP_.map))

/// this is used down below!
extern "C" { typedef const char* ParamNameFn(unsigned int); }

TEST(ChemistryData, SynchronizedMembers) {
  // infer {"member_name" : "member_type"} mapping from our "manual reflection"
  // functions
  std::vector<std::pair<const char*, ParamNameFn*>> pairs = {
    {"int", &param_name_int},
    {"double", &param_name_double},
    {"string", &param_name_string}
  };
  MemberTypeNameMap member_map;
  for (const auto& [type_name, fn_ptr]: pairs) {
    const unsigned int n_params = grackle_num_params(type_name);
    for (unsigned int i = 0; i < n_params; i++) {
      member_map[fn_ptr(i)] = type_name;
    }
  }

  // here we check
  EXPECT_MEMBERS(member_map, "chemistry_data");
}
