/***********************************************************************
/
/ Initialize Cloudy cooling data
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include <cmath> // std::log10, std::pow
#include <cstdio> // std::fprintf, std::snprintf, stderr, stdio
#include <cstring> // std::strcmp
#include "hdf5.h"
#include "grackle.h"
#include "internal_units.h"
#include "initialize_cloudy_data.hpp"
#include "support/h5io.hpp"

#define SMALL_LOG_VALUE (-99.0)
#define MAX_PARAMETER_NAME_LENGTH (512)


namespace { // stuff inside of an anonymous namespace is read-only

/// Initializes an empty #cloudy_data struct with zeros and nullptrs.
void initialize_empty_cloudy_data_struct(cloudy_data *my_cloudy)
{
  my_cloudy->grid_rank = 0LL;
  for (long long i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++){
    my_cloudy->grid_dimension[i] = 0LL;
    my_cloudy->grid_parameters[i] = nullptr;
  }
  my_cloudy->heating_data = nullptr;
  my_cloudy->cooling_data = nullptr;
  my_cloudy->mmw_data = nullptr;
  my_cloudy->data_size = 0LL;
}

double* load_heatcool_data(hid_t file_id, const char* dset_name,
                           double CoolUnit,
                           grackle::impl::h5io::ArrayShape expected_shape) {
  long long data_size =
    grackle::impl::h5io::ArrayShape_elem_count(expected_shape);
  double* tmp_data = new double[data_size];

  if (grackle_verbose) {
    std::fprintf(stdout,"Reading from \"%s\" dataset.\n", dset_name);
  }
  if (grackle::impl::h5io::read_dataset(file_id, dset_name, tmp_data,
                                        &expected_shape) != GR_SUCCESS) {
    delete[] tmp_data;
    return nullptr;
  }

  double* out = new double[data_size];
  for (long long q = 0LL; q < data_size; q++) {
    out[q] = tmp_data[q] > 0 ?
      (double) std::log10(tmp_data[q]) : (double) SMALL_LOG_VALUE;

    // Convert to code units.
    out[q] -= std::log10(CoolUnit);
  }
  delete[] tmp_data;

  return out;
}

}  // anonymous namespace


// initialize cloudy cooling data
int grackle::impl::initialize_cloudy_data(
    chemistry_data *my_chemistry,
    chemistry_data_storage *my_rates,
    cloudy_data *my_cloudy, const char *group_name,
    code_units *my_units, int read_data)
{

  char dset_name[MAX_PARAMETER_NAME_LENGTH];
  const std::size_t name_bufsize =
    static_cast<std::size_t>(MAX_PARAMETER_NAME_LENGTH);

  // Initialize things (to the null-state) even if cloudy cooling is not used.
  initialize_empty_cloudy_data_struct(my_cloudy);

  if (read_data == 0) { return GR_SUCCESS; }

  if (grackle_verbose) {
    std::fprintf(stdout,"Initializing Cloudy cooling: %s.\n", group_name);
    std::fprintf(stdout,"cloudy_table_file: %s.\n",my_chemistry->grackle_data_file);
  }

  // get the conversion factor to code units for the cooling rate
  InternalGrUnits internalu = new_internalu_legacy_init_cloudy_data_(my_units);
  double CoolUnit = internalu.coolunit;

  // Read cooling data in from hdf5 file.
  hid_t file_id = H5Fopen(my_chemistry->grackle_data_file,
                          H5F_ACC_RDONLY, H5P_DEFAULT);

  if (H5Aexists(file_id, "old_style")) {
    my_rates->cloudy_data_new = 0;
    if (grackle_verbose)
      std::fprintf(stdout, "Loading old-style Cloudy tables.\n");
  }

  // Determine the name of the dataset
  std::snprintf(dset_name, name_bufsize, "/CoolingRates/%s/Cooling",
                group_name);

  // Parse the grid properties from the dataset's attributes
  h5io::GridTableProps grid_props = h5io::parse_GridTableProps(file_id,
                                                               dset_name);
  if (!h5io::GridTableProps_is_valid(grid_props)) {
    h5io::drop_GridTableProps(&grid_props);
    H5Fclose (file_id);
    // error messages were already printed by h5io::parse_GridTableProps
    return GR_FAIL;
  }

  // use grid_props to initialize parts of my_cloudy
  // TODO: it would be really nice if we explicitly validated that each
  //       grid_props.axes[i] held the expected values
  my_cloudy->grid_rank = static_cast<long long>(grid_props.table_shape.ndim);
  for (int i = 0; i < grid_props.table_shape.ndim; i++) {
    my_cloudy->grid_dimension[i]
      = static_cast<long long>(grid_props.table_shape.shape[i]);

    bool is_temperature =
      std::strcmp("Temperature", grid_props.axes[i].name) == 0;

    double* buf = new double[my_cloudy->grid_dimension[i]];
    for (long long w = 0LL; w < my_cloudy->grid_dimension[i]; w++) {
      if (is_temperature) {
        buf[w] = std::log10(grid_props.axes[i].values[w]);
      } else {
        buf[w] = grid_props.axes[i].values[w];
      }
    }
    my_cloudy->grid_parameters[i] = buf;
  }

  if (grackle_verbose) {
    std::fprintf(stdout, "Cloudy cooling grid = {\n");
    std::fprintf(stdout, "  rank: %lld,\n", my_cloudy->grid_rank);
    for (long long i = 0LL; i < my_cloudy->grid_rank; i++) {
      std::fprintf(stdout,
          "  axis %lld (\"%s\"): %g to %g (%lld steps),\n",
          i,
          grid_props.axes[i].name,
          my_cloudy->grid_parameters[i][0],
          my_cloudy->grid_parameters[i][my_cloudy->grid_dimension[i]-1],
          my_cloudy->grid_dimension[i]);
    }
    fprintf(stdout, "}\n");
  }
  h5io::ArrayShape expected_shape = grid_props.table_shape;

  // Read Cooling data.
  my_cloudy->data_size = h5io::ArrayShape_elem_count(expected_shape);
  my_cloudy->cooling_data = load_heatcool_data(file_id, dset_name, CoolUnit,
                                               expected_shape);

  if (my_cloudy->cooling_data == nullptr) {
    h5io::drop_GridTableProps(&grid_props);
    H5Fclose(file_id);
    return GR_FAIL;
  }

  // Read Heating data.
  if (my_chemistry->UVbackground == 1) {
    std::snprintf(dset_name, name_bufsize, "/CoolingRates/%s/Heating",
                  group_name);

    // validate that Heating table has have identical GridTableProps
    if (h5io::assert_has_consistent_GridTableProps(file_id, dset_name,
                                                   grid_props) != GR_SUCCESS){
      h5io::drop_GridTableProps(&grid_props);
      H5Fclose(file_id);
      return GR_FAIL;
    }
    my_cloudy->heating_data = load_heatcool_data(file_id, dset_name, CoolUnit,
                                                 expected_shape);
    if (my_cloudy->heating_data == nullptr) {
      h5io::drop_GridTableProps(&grid_props);
      H5Fclose(file_id);
      return GR_FAIL;
    }
  }

  // Read MMW data.
  if (my_chemistry->primordial_chemistry == 0 &&
      std::strcmp(group_name, "Primordial") == 0) {
    const char* mmw_dset_name = "/CoolingRates/Primordial/MMW";

    // validate that MMW table has have identical GridTableProps
    if (h5io::assert_has_consistent_GridTableProps(file_id, mmw_dset_name,
                                                   grid_props) != GR_SUCCESS){
      h5io::drop_GridTableProps(&grid_props);
      H5Fclose(file_id);
      return GR_FAIL;
    }

    my_cloudy->mmw_data = new double[my_cloudy->data_size];
    if (h5io::read_dataset(file_id, mmw_dset_name, my_cloudy->mmw_data,
                           &expected_shape) != GR_SUCCESS) {
      h5io::drop_GridTableProps(&grid_props);
      H5Fclose(file_id);
      return GR_FAIL;
    }

  }

  H5Fclose (file_id);

  if (my_cloudy->grid_rank > GRACKLE_CLOUDY_TABLE_MAX_DIMENSION) {
    std::fprintf(
      stderr,
      "Error: rank of Cloudy cooling data must be less than or equal to %d.\n",
	    GRACKLE_CLOUDY_TABLE_MAX_DIMENSION);
    return GR_FAIL;
  }

  return GR_SUCCESS;
}

int grackle::impl::free_cloudy_data(cloudy_data *my_cloudy,
                                    chemistry_data *my_chemistry,
                                    int primordial) {

  for(int i = 0; i < my_cloudy->grid_rank; i++) {
    if (my_cloudy->grid_parameters[i] != nullptr) {
      delete[] my_cloudy->grid_parameters[i];
      my_cloudy->grid_parameters[i] = nullptr;
    }
  }

  if (my_cloudy->cooling_data != nullptr) {
    delete[] my_cloudy->cooling_data;
    my_cloudy->cooling_data = nullptr;
  }

  if (my_cloudy->heating_data != nullptr) {
    delete[] my_cloudy->heating_data;
    my_cloudy->heating_data = nullptr;
  }

  if (my_cloudy->mmw_data != nullptr) {
    delete[] my_cloudy->mmw_data;
    my_cloudy->mmw_data = nullptr;
  }
  return GR_SUCCESS;
}

