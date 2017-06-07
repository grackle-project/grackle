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

#include <stdlib.h>
#include <math.h>
#include "hdf5.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

extern int grackle_verbose;

// function prototypes
int read_dataset(hid_t file_id, char *dset_name, double *buffer);

// Initialize UV Background data
int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates)
{
  long long Nz;

  // Return if no UV background selected or using fully tabulated cooling.
  if (my_chemistry->UVbackground == 0 ||
      my_chemistry->primordial_chemistry == 0)
    return SUCCESS;


  if (grackle_verbose)
    fprintf(stdout, "Initializing UV background.\n");


  // Read in UV background data from hdf5 file.

  hid_t       file_id, dset_id, dspace_id;
  herr_t      status;
  herr_t      h5_error = -1;

  if (grackle_verbose)
    fprintf(stdout, "Reading UV background data from %s.\n",
            my_chemistry->grackle_data_file);
  file_id = H5Fopen(my_chemistry->grackle_data_file, 
                    H5F_ACC_RDONLY, H5P_DEFAULT);


  // Read Info dataset

  dset_id =  H5Dopen(file_id, "/UVBRates/Info");
  if (dset_id == h5_error) {
    fprintf(stderr, "Can't open 'Info' dataset in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  int strlen = (int)(H5Dget_storage_size(dset_id));
  char info_string[strlen+1];

  hid_t memtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(memtype, strlen+1);

  status = H5Dread(dset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, info_string);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read dataset 'Info'.\n");
    return FAIL;
  }

  H5Tclose(memtype); 
  H5Dclose(dset_id);



  // Open redshift dataset and get number of elements

  dset_id =  H5Dopen(file_id, "/UVBRates/z");
  if (dset_id == h5_error) {
    fprintf(stderr, "Can't open redshift dataset ('z') in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error) {
    fprintf(stderr, "Error opening dataspace for dataset 'z' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  Nz = H5Sget_simple_extent_npoints(dspace_id);
  if(Nz <= 0) {
    fprintf(stderr, "Redshift dataset ('z') has inappropriate size = %lld in %s.\n",
            Nz, my_chemistry->grackle_data_file);
    return FAIL;
  }

  H5Sclose(dspace_id);
  H5Dclose(dset_id);



  // Now allocate memory for UV background table.
  my_rates->UVbackground_table.Nz = Nz;

  my_rates->UVbackground_table.z = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.k24 = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.k25 = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.k26 = malloc(Nz * sizeof(double));

  if (my_chemistry->primordial_chemistry > 1) {
    my_rates->UVbackground_table.k27 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k28 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k29 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k30 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k31 = malloc(Nz * sizeof(double));
  }    

  my_rates->UVbackground_table.piHI = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.piHeII = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.piHeI = malloc(Nz * sizeof(double));

  if (my_chemistry->self_shielding_method > 0){
    my_rates->UVbackground_table.crsHI   = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.crsHeII = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.crsHeI  = malloc(Nz * sizeof(double));
  }


  // Now read everything.


  // *** Redshift ***
  if(! read_dataset(file_id, "/UVBRates/z",
                    my_rates->UVbackground_table.z) ) {
    fprintf(stderr, "Error reading dataset 'z' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** k24 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k24",
                    my_rates->UVbackground_table.k24) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k24' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** k25 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k25",
                    my_rates->UVbackground_table.k25) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k25' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** k26 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k26",
                    my_rates->UVbackground_table.k26) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k26' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  if (my_chemistry->primordial_chemistry > 1) {

    // *** k27 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k27",
                      my_rates->UVbackground_table.k27) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k27' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k28 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k28",
                      my_rates->UVbackground_table.k28) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k28' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k29 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k29",
                      my_rates->UVbackground_table.k29) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k29' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k30 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k30",
                      my_rates->UVbackground_table.k30) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k30' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k31 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k31",
                      my_rates->UVbackground_table.k31) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k31' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }
    
  }

  // *** piHI ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHI",
                    my_rates->UVbackground_table.piHI) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHI' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** piHeII ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHeII",
                    my_rates->UVbackground_table.piHeII) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHeII' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** piHeI ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHeI",
                    my_rates->UVbackground_table.piHeI) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHeI' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }


  if (my_chemistry->self_shielding_method > 0) {
    // *** crsHI ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/hi_avg_crs",
                      my_rates->UVbackground_table.crsHI) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/hi_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return FAIL;
    }

    // *** crsHeII ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/heii_avg_crs",
                    my_rates->UVbackground_table.crsHeII) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/heii_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return FAIL;
    }

    // *** crsHeI ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/hei_avg_crs",
                      my_rates->UVbackground_table.crsHeI) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/hei_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return FAIL;
    }
  }


  H5Fclose(file_id);


  // Get min/max of redshift vector
  my_rates->UVbackground_table.zmin = my_rates->UVbackground_table.z[0];
  my_rates->UVbackground_table.zmax = my_rates->UVbackground_table.z[Nz-1];

  // Print out some information about the dataset just read in.
  if (grackle_verbose) {
    fprintf(stdout, "UV background information:\n");
    fprintf(stdout, "  %s\n",info_string);
    fprintf(stdout, "  z_min = %6.3f\n  z_max = %6.3f\n",
            my_rates->UVbackground_table.zmin,
            my_rates->UVbackground_table.zmax);
  }

  // Set redshift on/off flags from data if not set.
  if (my_chemistry->UVbackground_redshift_on <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_on =
      my_rates->UVbackground_table.zmax;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_on to %f.\n",
              my_chemistry->UVbackground_redshift_on);
  }
  if (my_chemistry->UVbackground_redshift_fullon <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_fullon =
      my_rates->UVbackground_table.zmax;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_fullon to %f.\n",
              my_chemistry->UVbackground_redshift_fullon);
  }
  if (my_chemistry->UVbackground_redshift_drop <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_drop =
      my_rates->UVbackground_table.zmin;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_drop to %f.\n",
              my_chemistry->UVbackground_redshift_drop);
  }
  if (my_chemistry->UVbackground_redshift_off <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_off =
      my_rates->UVbackground_table.zmin;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_off to %f.\n",
              my_chemistry->UVbackground_redshift_off);
  }

  return SUCCESS;
}



int read_dataset(hid_t file_id, char *dset_name, double *buffer) {
  hid_t dset_id;
  herr_t status;
  herr_t h5_error = -1;

  dset_id =  H5Dopen(file_id, dset_name);
  if (dset_id == h5_error) {
    fprintf(stderr, "Failed to open dataset 'z'.\n");
    return FAIL;
  }

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read dataset 'z'.\n");
    return FAIL;
  }
 
  H5Dclose(dset_id);

  return SUCCESS;
}
