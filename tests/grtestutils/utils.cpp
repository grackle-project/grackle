#include "utils.hpp"

#include <cstring>

// the goal here is to make it easy to initialize chemistry_data without
// proliferating the GR_DATADIR macro throughout the test-code
// -> the GR_DATADIR macro is specified as a compiler arg and it specifies the
//    path to the data-directory

#define stringify(s) stringify_helper(s)
#define stringify_helper(s) #s

#define N_STANDARD_DATAFILES 7

static const char* const standard_data_files[N_STANDARD_DATAFILES] = {
  stringify(GR_DATADIR) "/CloudyData_UVB=FG2011.h5",
  stringify(GR_DATADIR) "/CloudyData_UVB=FG2011_shielded.h5",
  stringify(GR_DATADIR) "/CloudyData_UVB=HM2012.h5",
  stringify(GR_DATADIR) "/CloudyData_UVB=HM2012_high_density.h5",
  stringify(GR_DATADIR) "/CloudyData_UVB=HM2012_shielded.h5",
  stringify(GR_DATADIR) "/CloudyData_noUVB.h5",
  stringify(GR_DATADIR) "/cloudy_metals_2008_3D.h5"
};

/// returns the string-literal for the standard data file
static const char* standard_datafile_literal_(const char* datafile) {

  if (datafile==nullptr) {
    return nullptr;
  }

  // we get the number of characters in the prefix-path
  // -> strlen does not include the null character
  // -> we add 1 because we always insert a '/' between the prefix and the
  //    file's basename
  std::size_t prefix_len = 1 + std::strlen(stringify(GR_DATADIR));

  for (int i = 0; i < N_STANDARD_DATAFILES; i++){
    if (std::strcmp(datafile, standard_data_files[i]+prefix_len) == 0) {
      return standard_data_files[i];
    }
  }
  return nullptr;
}

std::optional<std::string> grtest::get_standard_datafile(const char* datafile) {
  const char* tmp = standard_datafile_literal_(datafile);
  if (tmp==nullptr) {
    return std::nullopt;
  }
  return std::optional<std::string>(tmp);
}

std::string grtest::to_pretty_string(double val) {
  char buf[30];
  snprintf(buf, 30, "%g", val);
  return std::string(buf);
}