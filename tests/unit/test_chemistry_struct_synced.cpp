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
#elif !defined(GR_SRC_ROOT_PATH)
#error "GR_SRC_ROOT_PATH macro isn't defined"
#elif !defined(IFLAG_DIR_GENERATED)
#else
#define stringify(s) stringify_helper(s)
#define stringify_helper(s) #s

static const char* GLOBAL_grackle_src_root_path = stringify(GR_SRC_ROOT_PATH);
static const char* GLOBAL_iflag_0 = "-I" stringify(IFLAG_DIR_GENERATED);
#endif

bool program_exists(std::string program) {
  grtest::Command cmd;
  cmd.program = "command";
  cmd.args = {"-v", program};
  cmd.stdout_path = "/dev/null";
  cmd.stderr_path = "/dev/null";

  return grtest::process_output(cmd).status == 0;
}

/// try to parse a string that holds a file where each line has the format
/// "<member_name>, <member_type>" and create a MemberTypeNameMap
///
/// @note
/// although std::regex is known to be **VERY** slow, it is convenient
static std::optional<MemberTypeNameMap> try_parse_(const std::string& s) {
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
      std::fprintf(stderr, "encountered parsing problem!\n");
      return {};
    }
    pos = lineend+1;
  }
  return {mapping};
}

/// Query the data-members of a struct. Returns a std::map where the keys are
/// the names of types and the values are the set of all member-names that have
/// that type. The result of this function is used for checking that efforts for
/// "manual reflection" of a struct's members are correct
///
/// @returns an empty option to denote **any** errors
///
/// In more detail, this function has 3 phases:
/// 1. castxml is used to parse a c++ header file (or compatible C header file)
///    & spit out an xml file that summarizes all of the contained declarations
/// 2. a python script is invoked to extract relevant information from the xml
///    file (it uses python's builtin xml-parsing library)
/// 3. we read in the extracted information and put it into the output format
///
/// @note
/// There are some issues with running this test if the source directory is
/// modified after this file is compile. In the future, a compelling case could
/// be made that the build-system should directly invoke castxml. For now, this
/// is a significant improvement over the python tests we are replacing.
///
/// @note
/// In the future, we'll be able to replace this library with a simplified
/// alternative. For example, we could use 1 of the following lightweight
/// librarie to directly infer members of aggregates (this includes C structs,
/// https://en.cppreference.com/w/cpp/language/aggregate_initialization) at
/// compile-time (using compiler-extensions)
/// - https://github.com/boostorg/pfr  (this doesn't need the rest of boost)
/// - https://github.com/stephenberry/glaze
/// Alteratively, C++26 is on track to make "reflection" of fields a standard
/// C++ feature.
std::optional<MemberTypeNameMap> query_struct_members(
  std::string struct_name, std::string path_relative_to_src_root
)
{
  // step 0: basic setup
  // ===================
  std::string src_root = GLOBAL_grackle_src_root_path + std::string("/");
  // -> get path to header file
  std::string hdr_path = src_root + path_relative_to_src_root;
  // -> get path to xml-reader
  std::string xml_reader = src_root + "tests/scripts/castxml_output_reader.py";
  // -> make a temporary directory (cleaned up when leaving function)
  grtest::ScratchDir temp_dir;
  // -> determine path to intermediate xml file
  std::string xml_path = temp_dir.get_path() + "/castxml-output.xml";

  // step 1: run castxml
  // ===================
  grtest::Command castxml_cmd;
  castxml_cmd.program = "castxml";
  castxml_cmd.args = {
    // tell castxml's compiler to preprocess/compile but don't link
    "-c",
    // tell castxml about IFLAGs
    GLOBAL_iflag_0,
    // tell castxml's compiler that the language is c++
    "-x", "c++",
    // the next required option configure castxml's internal Clang compiler.
    // The second part needs to specify an installed compiler (we may need
    // to support configuration of that path)
    "--castxml-cc-gnu", "g++",
    // specify the output xml format version
    "--castxml-output=1",
    // specify the output location
    "-o", xml_path,
    // specify the file to process
    hdr_path
  };
  castxml_cmd.stdout_path = ""; // <- capture stdout
  castxml_cmd.stderr_path = ""; // <- capture stderr
  grtest::ProcessOutput castxml_rslt = grtest::process_output(castxml_cmd);
  if ((castxml_rslt.status != 0) || (castxml_rslt.stdout_str.size() != 0)) {
    std::string tmp = grtest::summarize_cmd_and_rslt(castxml_cmd, castxml_rslt);
    std::fprintf(stderr, "unexpected outcome:\n%s\n", tmp.c_str());
    return {};
  } else if (castxml_rslt.stderr_str.size() != 0) {
    std::fprintf(stderr, "%s\n", castxml_rslt.stderr_str.c_str());
  }

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
    return {};
  } else if (castxml_rslt.stderr_str.size() != 0) {
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


#define EXPECT_MEMBERS(member_map, struct_name, path_relative_to_srcroot)     \
  if (! program_exists("castxml") ) {                                         \
    GTEST_SKIP() << "the castxml tool is required to run this test";          \
  }                                                                           \
                                                                              \
  /* use castxml to determine the actual members */                           \
  std::optional<MemberTypeNameMap> ref_map_TEMP_ =                            \
    query_struct_members( struct_name, path_relative_to_srcroot );            \
  if (! bool(ref_map_TEMP_) ) {                                               \
    GTEST_SKIP() << "issue inferring the members of " << struct_name;         \
  }                                                                           \
                                                                              \
  GTEST_EXPECT_TRUE(HasExpectedMembers(struct_name, member_map, *ref_map_TEMP_))

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
  EXPECT_MEMBERS(member_map, "chemistry_data", "src/include/grackle.h");
}
