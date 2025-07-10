
#include "gtest/gtest.h"
#include "grtest_os.hpp"

#include "grackle.h"
#include "data_file_utils.h"
#include "os_utils.h"

#include <cstdio>
#include <filesystem>
#include <memory>  // std::unique_ptr, std::make_unique
#include <string>
#include <vector>

struct SampleFileCase {
  std::string basename;
  std::string contents;
  std::string cksum;
};

// setup some simple cases
static std::vector<SampleFileCase> get_file_cases_() {
  return {
      {"abc.txt", "abc",
       "sha256:"
       "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
      {"longer.txt", "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
       "sha256:"
       "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"}};
}

class FileRegistryMemStorage {
  std::vector<std::string> line_vec;  // this holds the actual storage
  std::vector<const char*> ptr_vec;   // ptr_vec.data() is passed to the C funcs
public:
  // we delete the copy constructors since they invalidate the pointer
  // representation of the registry
  FileRegistryMemStorage(const FileRegistryMemStorage&) = delete;
  FileRegistryMemStorage& operator=(const FileRegistryMemStorage&) = delete;

  explicit FileRegistryMemStorage(
      const std::vector<SampleFileCase>& file_cases) {
    for (const SampleFileCase& file_case : file_cases) {  // fill up line_vec
      line_vec.push_back(file_case.basename + "  " + file_case.cksum);
    }
    // fill up ptr_vec once line_vec's contents are frozen (this is important
    // to avoid pointer invalidation due to std::string's small string
    // optimization)
    for (std::string& line : line_vec) {
      ptr_vec.push_back(line.data());
    }
    ptr_vec.push_back(nullptr);
  }

  [[nodiscard]] const char** get_registry() { return ptr_vec.data(); }
};

/// This is used for testing the automatic file machinery
/// - Basically, this fixture creates a dummy data directory with dummy data
///   files that have known checksums
/// - It provides a method to create the file-registry for the data directory
class AutoFileTest : public testing::Test {
  std::vector<SampleFileCase> file_cases;
  std::unique_ptr<grtest::EnvManager> env_manager;
  std::unique_ptr<grtest::TempDir> tmp_dir;

protected:
  void SetUp() override {
    if (get_platform_() == platform_kind_unknown) {
      GTEST_SKIP() << "the test can't be run on an unknown platform";
    }
    tmp_dir = std::make_unique<grtest::TempDir>(grtest::TempDir::create());
    std::filesystem::path data_dir = tmp_dir->get_path();
    grackle_version version = get_grackle_version();
    std::filesystem::path data_file_dir =
        (data_dir / "data-store-v1" / version.version);
    std::filesystem::create_directories(data_file_dir);
    file_cases = get_file_cases_();
    for (const auto& file_case : file_cases) {
      std::filesystem::path path = data_file_dir / file_case.basename;
      std::string path_str = path.string();
      std::FILE* fp = std::fopen(path.c_str(), "w");
      std::fprintf(fp, "%s", file_case.contents.c_str());
      std::fclose(fp);
    }

    env_manager = grtest::EnvManager::create();
    if (env_manager == nullptr) {
      GTEST_SKIP() << "EnvManager isn't supported on current platform";
    }
    env_manager->override_envvar("GRACKLE_DATA_DIR", data_dir.string());
  }

  /// creates the registry
  [[nodiscard]] FileRegistryMemStorage get_file_registry(
      bool swap_fnames) const {
    if (swap_fnames) {
      std::vector<SampleFileCase> tmp = file_cases;
      std::size_t size = file_cases.size();
      for (std::size_t i = 0; (i + 1) < size; ++i) {
        tmp[i].basename = file_cases[i + 1].basename;
      }
      tmp[size - 1].basename = file_cases[0].basename;
      return FileRegistryMemStorage(tmp);
    } else {
      return FileRegistryMemStorage(file_cases);
    }
  }
};

TEST_F(AutoFileTest, NullArg) {
  // we supress (really capture) the error messages since we're testing that
  // the function properly fails
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);

  FileRegistryMemStorage tmp = get_file_registry(false);

  generic_file_props file_props =
      determine_data_file_(nullptr, GR_DFOPT_MANAGED, tmp.get_registry());
  EXPECT_EQ(file_props.path, nullptr)
      << "determine_data_file_ should set path to NULL when we lookup a file "
      << "not in the registry";
  free_generic_file_props_(&file_props);
}

TEST_F(AutoFileTest, ManagedWithCksum) {
  FileRegistryMemStorage tmp = get_file_registry(false);

  std::vector<std::string> basenames = {"abc.txt", "longer.txt"};
  for (const std::string& name : basenames) {
    generic_file_props file_props = determine_data_file_(
        name.c_str(), GR_DFOPT_MANAGED, tmp.get_registry());
    if (file_props.path == nullptr) {
      free_generic_file_props_(&file_props);
      GTEST_FAIL() << "there was a problem retrieving the path for " << name;
    }
    ASSERT_TRUE(std::filesystem::exists(file_props.path))
        << "the retrieved path for " << name << ", " << file_props.path
        << " does not exist";
    free_generic_file_props_(&file_props);
  }
}

TEST_F(AutoFileTest, ManagedWithCksumNotAFile) {
  // we supress (really capture) the error messages since we're testing that
  // the function properly fails
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);

  FileRegistryMemStorage tmp = get_file_registry(false);

  generic_file_props file_props =
      determine_data_file_("not-a-file", GR_DFOPT_MANAGED, tmp.get_registry());
  EXPECT_EQ(file_props.path, nullptr)
      << "determine_data_file_ should set path to NULL when we lookup a file "
      << "not in the registry";
  free_generic_file_props_(&file_props);
}

TEST_F(AutoFileTest, ManagedWithCksumWrongCksum) {
  // we supress (really capture) the error messages since we're testing that
  // the function properly fails
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);

  // we are intentionally making a registry with incorrect cksums
  FileRegistryMemStorage tmp = get_file_registry(true);

  std::vector<std::string> basenames = {"abc.txt", "longer.txt"};
  for (const std::string& name : basenames) {
    generic_file_props file_props = determine_data_file_(
        name.c_str(), GR_DFOPT_MANAGED, tmp.get_registry());
    EXPECT_EQ(file_props.path, nullptr)
        << "determine_data_file_ should set path to NULL when we intentionally "
        << "create a registry where we mismatch file names and checksums";
    free_generic_file_props_(&file_props);
  }
}

TEST_F(AutoFileTest, ManagedNoCksum) {
  std::unique_ptr<grtest::CaptureSink> capstderr =
      grtest::CaptureSink::create(stderr);

  // we are intentionally making a registry with incorrect cksums
  FileRegistryMemStorage wrong_reg = get_file_registry(true);

  std::vector<std::string> basenames = {"abc.txt", "longer.txt"};
  for (const std::string& name : basenames) {
    generic_file_props file_props = determine_data_file_(
        name.c_str(), GR_DFOPT_MANAGED_NO_CKSUM, wrong_reg.get_registry());
    if (file_props.path == nullptr) {
      free_generic_file_props_(&file_props);
      GTEST_FAIL() << "there was a problem retrieving the path for " << name
                   << "while running with GR_DFOPT_MANAGED_NO_CKSUM";
    }
    ASSERT_TRUE(std::filesystem::exists(file_props.path))
        << "the retrieved path for " << name << ", " << file_props.path
        << " does not exist";
    free_generic_file_props_(&file_props);
  }
}
