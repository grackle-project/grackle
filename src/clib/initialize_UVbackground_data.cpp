/***********************************************************************
/
/ Initialize UV background data
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include <cstdio>
#include <vector>
#include "hdf5.h"
#include "grackle.h"
#include "grackle_macros.h"
#include "utils/h5io.hpp"

#include "initialize_UVbackground_data.hpp"

namespace { // stuff inside an anonymous namespace is only used in this file

/// Initializes an empty #UVBtable struct with zeros and nullptrs.
void initialize_empty_UVBtable_struct(UVBtable *table)
{
  table->Nz     = 0LL;
  table->z      = nullptr;
  table->k24    = nullptr;
  table->k25    = nullptr;
  table->k26    = nullptr;
  table->k27    = nullptr;
  table->k28    = nullptr;
  table->k29    = nullptr;
  table->k30    = nullptr;
  table->k31    = nullptr;
  table->piHI   = nullptr;
  table->piHeI  = nullptr;
  table->piHeII = nullptr;
  table->crsHI  = nullptr;
  table->crsHeII = nullptr;
  table->crsHeI = nullptr;
}

} // anonymous namespace

// Initialize UV Background data
int grackle::impl::initialize_UVbackground_data(chemistry_data *my_chemistry,
                                                chemistry_data_storage *my_rates)
{
  initialize_empty_UVBtable_struct(&my_rates->UVbackground_table);

  // Return if no UV background selected or using fully tabulated cooling.
  if (my_chemistry->UVbackground == 0 ||
      my_chemistry->primordial_chemistry == 0)
    return GR_SUCCESS;


  if (grackle_verbose)
    std::fprintf(stdout, "Initializing UV background.\n");


  // Read in UV background data from hdf5 file.

  hid_t       file_id, dset_id;
  herr_t      status;
  herr_t      h5_error = -1;

  if (grackle_verbose)
    std::fprintf(stdout, "Reading UV background data from %s.\n",
            my_chemistry->grackle_data_file);
  file_id = H5Fopen(my_chemistry->grackle_data_file,
                    H5F_ACC_RDONLY, H5P_DEFAULT);


  // Read Info dataset

  dset_id =  H5Dopen(file_id, "/UVBRates/Info");
  if (dset_id == h5_error) {
    std::fprintf(stderr, "Can't open 'Info' dataset in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  int strlen = (int)(H5Dget_storage_size(dset_id));
  std::vector<char> info_string(strlen+1);

  hid_t memtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(memtype, strlen+1);

  status = H5Dread(dset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   info_string.data());
  if (status == h5_error) {
    std::fprintf(stderr, "Failed to read dataset 'Info'.\n");
    return GR_FAIL;
  }

  H5Tclose(memtype);
  H5Dclose(dset_id);

  // Open redshift dataset and get number of elements

  const h5io::ArrayShape common_shape
    = h5io::read_dataset_shape(file_id, "/UVBRates/z");
  if (!h5io::ArrayShape_is_valid(common_shape)) {
    return GR_FAIL; // error messages are already printed
  } else if (common_shape.ndim != 1 || common_shape.shape[0] < 0) {
    std::fprintf(
        stderr,
        "Redshift dataset (\"/UVBRates/z\") in %s has inappropriate shape\n",
        my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  long long Nz = static_cast<long long>(common_shape.shape[0]);

  // Now allocate memory for UV background table.
  my_rates->UVbackground_table.Nz = Nz;

  my_rates->UVbackground_table.z = new double[Nz];
  my_rates->UVbackground_table.k24 = new double[Nz];
  my_rates->UVbackground_table.k25 = new double[Nz];
  my_rates->UVbackground_table.k26 = new double[Nz];

  if (my_chemistry->primordial_chemistry > 1) {
    my_rates->UVbackground_table.k27 = new double[Nz];
    my_rates->UVbackground_table.k28 = new double[Nz];
    my_rates->UVbackground_table.k29 = new double[Nz];
    my_rates->UVbackground_table.k30 = new double[Nz];
    my_rates->UVbackground_table.k31 = new double[Nz];
  }

  my_rates->UVbackground_table.piHI = new double[Nz];
  my_rates->UVbackground_table.piHeII = new double[Nz];
  my_rates->UVbackground_table.piHeI = new double[Nz];

  if (my_chemistry->self_shielding_method > 0){
    my_rates->UVbackground_table.crsHI   = new double[Nz];
    my_rates->UVbackground_table.crsHeII = new double[Nz];
    my_rates->UVbackground_table.crsHeI  = new double[Nz];
  }


  // Now read everything.

  using ::grackle::impl::h5io::read_dataset;

  // *** Redshift ***
  if(! read_dataset(file_id, "/UVBRates/z",
                    my_rates->UVbackground_table.z,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset 'z' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  // *** k24 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k24",
                    my_rates->UVbackground_table.k24,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k24' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  // *** k25 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k25",
                    my_rates->UVbackground_table.k25,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k25' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  // *** k26 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k26",
                    my_rates->UVbackground_table.k26,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k26' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  if (my_chemistry->primordial_chemistry > 1) {

    // *** k27 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k27",
                      my_rates->UVbackground_table.k27,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k27' in %s.\n",
              my_chemistry->grackle_data_file);
      return GR_FAIL;
    }

    // *** k28 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k28",
                      my_rates->UVbackground_table.k28,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k28' in %s.\n",
              my_chemistry->grackle_data_file);
      return GR_FAIL;
    }

    // *** k29 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k29",
                      my_rates->UVbackground_table.k29,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k29' in %s.\n",
              my_chemistry->grackle_data_file);
      return GR_FAIL;
    }

    // *** k30 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k30",
                      my_rates->UVbackground_table.k30,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k30' in %s.\n",
              my_chemistry->grackle_data_file);
      return GR_FAIL;
    }

    // *** k31 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k31",
                      my_rates->UVbackground_table.k31,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k31' in %s.\n",
              my_chemistry->grackle_data_file);
      return GR_FAIL;
    }

  }

  // *** piHI ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHI",
                    my_rates->UVbackground_table.piHI,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHI' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  // *** piHeII ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHeII",
                    my_rates->UVbackground_table.piHeII,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHeII' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }

  // *** piHeI ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHeI",
                    my_rates->UVbackground_table.piHeI,
                    /* expected_shape = */ &common_shape) ) {
    std::fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHeI' in %s.\n",
            my_chemistry->grackle_data_file);
    return GR_FAIL;
  }


  if (my_chemistry->self_shielding_method > 0) {
    // *** crsHI ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/hi_avg_crs",
                      my_rates->UVbackground_table.crsHI,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/hi_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      std::fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return GR_FAIL;
    }

    // *** crsHeII ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/heii_avg_crs",
                      my_rates->UVbackground_table.crsHeII,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/heii_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      std::fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return GR_FAIL;
    }

    // *** crsHeI ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/hei_avg_crs",
                      my_rates->UVbackground_table.crsHeI,
                      /* expected_shape = */ &common_shape) ) {
      std::fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/hei_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      std::fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return GR_FAIL;
    }
  }


  H5Fclose(file_id);


  // Get min/max of redshift vector
  my_rates->UVbackground_table.zmin = my_rates->UVbackground_table.z[0];
  my_rates->UVbackground_table.zmax = my_rates->UVbackground_table.z[Nz-1];

  // Print out some information about the dataset just read in.
  if (grackle_verbose) {
    std::fprintf(stdout, "UV background information:\n");
    std::fprintf(stdout, "  %s\n",info_string.data());
    std::fprintf(stdout, "  z_min = %6.3f\n  z_max = %6.3f\n",
            my_rates->UVbackground_table.zmin,
            my_rates->UVbackground_table.zmax);
  }

  // Set redshift on/off flags from data if not set.
  if (my_chemistry->UVbackground_redshift_on <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_on =
      my_rates->UVbackground_table.zmax;
    if (grackle_verbose)
      std::fprintf(stdout, "Setting UVbackground_redshift_on to %f.\n",
              my_chemistry->UVbackground_redshift_on);
  }
  if (my_chemistry->UVbackground_redshift_fullon <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_fullon =
      my_rates->UVbackground_table.zmax;
    if (grackle_verbose)
      std::fprintf(stdout, "Setting UVbackground_redshift_fullon to %f.\n",
              my_chemistry->UVbackground_redshift_fullon);
  }
  if (my_chemistry->UVbackground_redshift_drop <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_drop =
      my_rates->UVbackground_table.zmin;
    if (grackle_verbose)
      std::fprintf(stdout, "Setting UVbackground_redshift_drop to %f.\n",
              my_chemistry->UVbackground_redshift_drop);
  }
  if (my_chemistry->UVbackground_redshift_off <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_off =
      my_rates->UVbackground_table.zmin;
    if (grackle_verbose)
      std::fprintf(stdout, "Setting UVbackground_redshift_off to %f.\n",
              my_chemistry->UVbackground_redshift_off);
  }

  return GR_SUCCESS;
}

void grackle::impl::free_UVBtable(UVBtable *table)
{
  // define a function that loosely approximates GRACKLE_FREE
  auto cleanup_fn = [](double*& ptr) {
    if (ptr != nullptr) {
      delete[] ptr;
      ptr = nullptr;
    }
  };

  cleanup_fn(table->z);
  cleanup_fn(table->k24);
  cleanup_fn(table->k25);
  cleanup_fn(table->k26);
  cleanup_fn(table->k27);
  cleanup_fn(table->k28);
  cleanup_fn(table->k29);
  cleanup_fn(table->k30);
  cleanup_fn(table->k31);
  cleanup_fn(table->piHI);
  cleanup_fn(table->piHeII);
  cleanup_fn(table->piHeI);
  cleanup_fn(table->crsHI);
  cleanup_fn(table->crsHeII);
  cleanup_fn(table->crsHeI);
}
