#include "grtest_utils.hpp"

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

bool grtest::set_standard_datafile(
  chemistry_data& my_chemistry, const char* datafile
) {

  if (datafile==NULL) {
    return false; // we should probably abort the program with an error
  }

  // we get the number of characters in the prefix-path
  // -> strlen does not include the null character
  // -> we add 1 because we always insert a '/' between the prefix and the
  //    file's basename
  std::size_t prefix_len = 1 + std::strlen(stringify(GR_DATADIR));

  for (int i = 0; i < N_STANDARD_DATAFILES; i++){
    if (std::strcmp(datafile, standard_data_files[i]+prefix_len) == 0) {
      my_chemistry.grackle_data_file = standard_data_files[i];
      return true;
    }
  }
  return false;
}
