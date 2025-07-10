#include <cstdlib>
#include <filesystem>
#include <map>
#include <optional>
#include <string>

#include "gtest/gtest.h"

#include "grtest_os.hpp"
#include "os_utils.h"
#include "status_reporting.h"

// UNCOMMENT in PR 246
// #include "grtest_cmd.hpp"

/// a very standard shorthand for abbreviating std::filesystem
namespace fs = std::filesystem;

using optional_str = std::optional<std::string>;
using optional_path = std::optional<std::filesystem::path>;

enum class DataDirImpl {
  get_data_dir_,  /// corresponds to the C function
  grdata_tool     /// corresponds to the external grdata tool
};

/// calls get_data_dir and wraps the result in a std::filesystem::path
///
/// The purpose is 3-fold:
/// 1. It makes string-comparisons easier
/// 2. We don't have to worry about memory leaks if a test fails
/// 3. We avoid polluting the output with error-messages when we expect the
///    code to fail
static std::optional<fs::path> get_dir_(
    enum platform_kind kind, DataDirImpl impl = DataDirImpl::get_data_dir_) {
  // capture all stderr within the function
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);
  switch (impl) {
    case DataDirImpl::get_data_dir_: {
      char* tmp = get_data_dir_(kind);
      if (tmp == nullptr) {
        return {};
      }
      fs::path out{tmp};
      free(tmp);
      return {out};
    }
    case DataDirImpl::grdata_tool: {
      return {};  // this explicitly will produce an error
      // UNCOMMENT & modify in PR #246 (I already did some basic checks)
      // grtest::ProcessStatusAndStdout rslt =
      //   grtest::capture_status_and_output("@PATH_TO_GRDATA@ show-data-dir");
      // if (rslt.status != 0) { return {}; }
      // // remove any trailing whitespace
      // std::string tmp = rslt.stdout_str.substr(
      //   0, rslt.stdout_str.find_last_not_of(" \t\n")+1);
      // return {std::filesystem::path(tmp)};
    }
  }
  GR_INTERNAL_ERROR("THIS SHOULD BE UNREACHABLE!");
}

/// a test-fixture for testing get_data_dir_
///
/// The fixture provides the override_envvar to let you override the value of
/// an environment variable and at the end of the test, it restores the
/// environment to the original state
class DataDirTest : public testing::TestWithParam<DataDirImpl> {
  std::unique_ptr<grtest::EnvManager> env_manager;

protected:
  void SetUp() override {
    env_manager = grtest::EnvManager::create();
    if (env_manager == nullptr) {
      GTEST_SKIP()
        << "EnvManager doesn't seem to be supported on the current platform";
    }
  }

  /// overrides the value of the environment variable.
  void override_envvar(const std::string& envvar,
                       const std::optional<std::string>& val) {
    env_manager->override_envvar(envvar, val);
  }

  void override_envvar(const std::string& envvar, const std::string& val) {
    env_manager->override_envvar(envvar, val);
  }

  void override_envvar(const std::string& envvar, const char* val) {
    override_envvar(envvar, std::make_optional(val));
  }
};

TEST_P(DataDirTest, NonUnix) {
  override_envvar("GRACKLE_DATA_DIR", "/my/path1");
  override_envvar("XDG_DATA_HOME", "/my/path2");
  override_envvar("HOME", "/my/path3");
  optional_path rslt = get_dir_(platform_kind_unknown, GetParam());
  ASSERT_FALSE(rslt.has_value())
      << "the functionality should currently error (gracefully) on non-Unix "
      << "platforms (regardless of the env variable values";
}

TEST_P(DataDirTest, NoEnvVarsUnix) {
  override_envvar("GRACKLE_DATA_DIR", std::nullopt);
  override_envvar("XDG_DATA_HOME", std::nullopt);
  override_envvar("HOME", std::nullopt);
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_FALSE(rslt.has_value())
      << "you shouldn't be able to determine a when none of the relevant ENV "
      << "variables are defined";
}

TEST_P(DataDirTest, PrimaryEnvVar) {
  override_envvar("GRACKLE_DATA_DIR", "/path/to/grackle-data");
  override_envvar("XDG_DATA_HOME", "/path/to/data-home");
  override_envvar("HOME", "/path/to/home");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/path/to/grackle-data"))
      << "unexpected path when GRACKLE_DATA_DIR is specified";

  // confirm that changing other env-vars has no impact
  override_envvar("XDG_DATA_HOME", std::nullopt);
  override_envvar("HOME", "rel-path/to/home");
  rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/path/to/grackle-data"))
      << "other env-vars shouldn't matter when GRACKLE_DATA_DIR is specified";
}

TEST_P(DataDirTest, PrimaryEnvVarRelPath) {
  override_envvar("GRACKLE_DATA_DIR", "relpath/to/grackle-data");
  override_envvar("XDG_DATA_HOME", "/path/to/data-home");
  override_envvar("HOME", "/path/to/home");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_FALSE(rslt.has_value())
      << "setting GRACKLE_DATA_DIR to a relative path should be a problem";
}

TEST_P(DataDirTest, PrimaryEnvVarTilde) {
  override_envvar("GRACKLE_DATA_DIR", "~/path/to/grackle-data");
  override_envvar("XDG_DATA_HOME", "/path/to/data-home");
  override_envvar("HOME", "/path/to/home");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_FALSE(rslt.has_value())
      << "envvar paths starting with '~' should be erroneous";
}

TEST_P(DataDirTest, XdgDataHome) {
  override_envvar("GRACKLE_DATA_DIR", std::nullopt);
  override_envvar("XDG_DATA_HOME", "/abs/path");
  override_envvar("HOME", "/wrong/abs/path");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/abs/path/grackle"))
      << "unexpected path when XDG_DATA_HOME (but not GRACKLE_DATA_DIR) is set";

  override_envvar("XDG_DATA_HOME", "/abs/path/alt/");
  rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/abs/path/alt/grackle"))
      << "unexpected path when GRACKLE_DATA_DIR is unset and XDG_DATA_HOME "
      << "is assigned a path with a trailing slash";
}

TEST_P(DataDirTest, XdgDataHomeRelPath) {
  override_envvar("GRACKLE_DATA_DIR", std::nullopt);
  override_envvar("XDG_DATA_HOME", "rel/path");
  override_envvar("HOME", "/path/to/home");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/path/to/home/.local/share/grackle"))
      << "unexpected path when GRACKLE_DATA_DIR is unset and XDG_DATA_HOME "
      << "is assigned a relative path (i.e. it should use the value of HOME)";

  override_envvar("HOME", std::nullopt);
  rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_FALSE(rslt.has_value())
      << "when GRACKLE_DATA_DIR and HOME are unset, and XDG_DATA_HOME "
      << "is assigned a relative path, there should be an error";
}

TEST_P(DataDirTest, Home) {
  override_envvar("GRACKLE_DATA_DIR", std::nullopt);
  override_envvar("XDG_DATA_HOME", std::nullopt);
  override_envvar("HOME", "/path/to/home");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/path/to/home/.local/share/grackle"))
      << "unexpected path when HOME is the only relevant env var specified";

  // confirm that adding a trailing '/' isn't a problem
  override_envvar("HOME", "/path/to/home-alt/");
  rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_EQ(rslt.value(), fs::path("/path/to/home-alt/.local/share/grackle"))
      << "unexpected path when HOME is the only relevant env var specified & "
      << "has a trailing slash";
}

TEST_P(DataDirTest, HomeRelPath) {
  override_envvar("GRACKLE_DATA_DIR", std::nullopt);
  override_envvar("XDG_DATA_HOME", std::nullopt);
  override_envvar("HOME", "rel-path/to/home");
  optional_path rslt = get_dir_(platform_kind_generic_unix, GetParam());
  ASSERT_FALSE(rslt.has_value())
      << "unexpected path when HOME is the only relevant env var specified & "
      << "it holds a relative path";
}

// In PR #246, Uncomment DataDirImpl::grdata_tool
INSTANTIATE_TEST_SUITE_P(
    , /* <- leaving Instantiation name empty */
    DataDirTest,
    testing::Values(DataDirImpl::get_data_dir_ /*, DataDirImpl::grdata_tool */),
    [](const testing::TestParamInfo<DataDirTest::ParamType>& info) {
      // this function describes custom formatting
      switch (info.param) {
        case DataDirImpl::get_data_dir_:
          return std::string("CFn_get_data_dir");
        case DataDirImpl::grdata_tool: {
          return std::string("grdata_cli_tool");
        }
      }
      return std::string("SomethingIsWrong");
    });