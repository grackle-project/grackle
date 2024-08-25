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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf5.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

#define SMALL_LOG_VALUE -99.0

extern int grackle_verbose;


/**
 * Initializes an empty #cloudy_data struct with zeros and NULLs.
 */
void initialize_empty_cloudy_data_struct(cloudy_data *my_cloudy)
{
  my_cloudy->grid_rank = 0LL;
  for (long long i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++){
    my_cloudy->grid_dimension[i] = 0LL;
    my_cloudy->grid_parameters[i] = NULL;
  }
  my_cloudy->heating_data = NULL;
  my_cloudy->cooling_data = NULL;
  my_cloudy->mmw_data = NULL;
  my_cloudy->data_size = 0LL;
}

// Initialize Cloudy cooling data
int initialize_cloudy_data(const char* path, chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units, int read_data)
{

  long long q, w;
  double *temp_data;
  long long temp_int;
  long long *temp_int_arr;
  char parameter_name[MAX_LINE_LENGTH];

  // Initialize things (to the null-state) even if cloudy cooling is not used.
  initialize_empty_cloudy_data_struct(my_cloudy);

  if (read_data == 0) { return SUCCESS; }

  if (grackle_verbose) {
    fprintf(stdout,"Initializing Cloudy cooling: %s.\n", group_name);
    fprintf(stdout,"cloudy_table_file: %s.\n",path);
  }

  /* Get conversion units. */

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  double tbase1 = my_units->time_units;
  double xbase1 = co_length_units /
    (my_units->a_value * my_units->a_units);
  double dbase1 = co_density_units *
    POW(my_units->a_value * my_units->a_units, 3);
  double mh = 1.67e-24;
  double CoolUnit = (POW(my_units->a_units,5) * POW(xbase1,2) * POW(mh,2)) /
                    (POW(tbase1,3) * dbase1);

  // Read cooling data in from hdf5 file.
  hid_t       file_id, dset_id, attr_id;
  herr_t      status;
  herr_t      h5_error = -1;

  file_id = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (H5Aexists(file_id, "old_style")) {
    my_rates->cloudy_data_new = 0;
    if (grackle_verbose)
      fprintf(stdout, "Loading old-style Cloudy tables.\n");
  }

  // Open cooling dataset and get grid dimensions.

  sprintf(parameter_name, "/CoolingRates/%s/Cooling", group_name);
  dset_id =  H5Dopen(file_id, parameter_name);
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open Cooling in %s.\n", path);
    return FAIL;
  }

  // Grid rank.
  attr_id = H5Aopen_name(dset_id, "Rank");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open Rank attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, &temp_int);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Rank attribute in Cooling dataset.\n");
    return FAIL;
  }
  my_cloudy->grid_rank = (long long) temp_int;
  if (grackle_verbose)
    fprintf(stdout,"Cloudy cooling grid rank: %lld.\n", my_cloudy->grid_rank);
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Rank attribute in Cooling dataset.\n");
    return FAIL;
  }

  // Grid dimension.
  temp_int_arr = malloc(my_cloudy->grid_rank * sizeof(long long));
  attr_id = H5Aopen_name(dset_id, "Dimension");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8,temp_int_arr);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  if (grackle_verbose)
    fprintf(stdout,"Cloudy cooling grid dimensions:");
  for (q = 0;q < my_cloudy->grid_rank;q++) {
    my_cloudy->grid_dimension[q] = (long long) temp_int_arr[q];
    if (grackle_verbose)
      fprintf(stdout," %lld", my_cloudy->grid_dimension[q]);
  }
  if (grackle_verbose)
    fprintf(stdout,".\n");
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  free(temp_int_arr);

  // Grid parameters.
  for (q = 0;q < my_cloudy->grid_rank;q++) {

    if (q < my_cloudy->grid_rank - 1) {
      sprintf(parameter_name,"Parameter%lld",(q+1));
    }
    else {
      sprintf(parameter_name,"Temperature");
    }

    temp_data = malloc(my_cloudy->grid_dimension[q] * sizeof(double));

    attr_id = H5Aopen_name(dset_id, parameter_name);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open %s attribute in Cooling dataset.\n",
              parameter_name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_R8, temp_data);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read %s attribute in Cooling dataset.\n",
              parameter_name);
      return FAIL;
    }

    my_cloudy->grid_parameters[q] = malloc(my_cloudy->grid_dimension[q] *
                                           sizeof(double));
    for (w = 0;w < my_cloudy->grid_dimension[q];w++) {
      if (q < my_cloudy->grid_rank - 1) {
	my_cloudy->grid_parameters[q][w] = (double) temp_data[w];
      }
      else {
	// convert temeperature to log
	my_cloudy->grid_parameters[q][w] = (double) log10(temp_data[w]);
      }

    }
    if (grackle_verbose)
      fprintf(stdout,"%s: %"GSYM" to %"GSYM" (%lld steps).\n",parameter_name,
              my_cloudy->grid_parameters[q][0],
              my_cloudy->grid_parameters[q][my_cloudy->grid_dimension[q]-1],
              my_cloudy->grid_dimension[q]);
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close %s attribute in Cooling dataset.\n",
              parameter_name);
      return FAIL;
    }
    free(temp_data);

  }

  // Read Cooling data.
  my_cloudy->data_size = 1;
  for (q = 0;q < my_cloudy->grid_rank;q++) {
    my_cloudy->data_size *= my_cloudy->grid_dimension[q];
  }
  temp_data = malloc(my_cloudy->data_size * sizeof(double));

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
  if (grackle_verbose)
    fprintf(stdout,"Reading Cloudy Cooling dataset.\n");
  if (status == h5_error) {
    fprintf(stderr,"Failed to read Cooling dataset.\n");
    return FAIL;
  }

  my_cloudy->cooling_data = malloc(my_cloudy->data_size * sizeof(double));
  for (q = 0;q < my_cloudy->data_size;q++) {
    my_cloudy->cooling_data[q] = temp_data[q] > 0 ?
      (double) log10(temp_data[q]) : (double) SMALL_LOG_VALUE;

    // Convert to code units.
    my_cloudy->cooling_data[q] -= log10(CoolUnit);
  }
  free(temp_data);

  status = H5Dclose(dset_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close Cooling dataset.\n");
    return FAIL;
  }

  // Read Heating data.
  if (my_chemistry->UVbackground == 1) {

    temp_data = malloc(my_cloudy->data_size * sizeof(double));

    sprintf(parameter_name, "/CoolingRates/%s/Heating", group_name);
    dset_id =  H5Dopen(file_id, parameter_name);
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open Heating in %s.\n", path);
      return FAIL;
    }

    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    if (grackle_verbose)
      fprintf(stdout,"Reading Cloudy Heating dataset.\n");
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Heating dataset.\n");
      return FAIL;
    }

    my_cloudy->heating_data = malloc(my_cloudy->data_size * sizeof(double));
    for (q = 0;q < my_cloudy->data_size;q++) {
      my_cloudy->heating_data[q] = temp_data[q] > 0 ?
        (double) log10(temp_data[q]) : (double) SMALL_LOG_VALUE;

      // Convert to code units.
      my_cloudy->heating_data[q] -= log10(CoolUnit);
    }
    free(temp_data);

    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Heating dataset.\n");
      return FAIL;
    }
  }

  // Read MMW data.
  if (my_chemistry->primordial_chemistry == 0 &&
      strcmp(group_name, "Primordial") == 0) {

    my_cloudy->mmw_data = malloc(my_cloudy->data_size * sizeof(double));

    sprintf(parameter_name, "/CoolingRates/%s/MMW", group_name);
    dset_id =  H5Dopen(file_id, parameter_name);
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open MMW in %s.\n", path);
      return FAIL;
    }

    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     my_cloudy->mmw_data);
    if (grackle_verbose)
      fprintf(stdout,"Reading Cloudy MMW dataset.\n");
    if (status == h5_error) {
      fprintf(stderr,"Failed to read MMW dataset.\n");
      return FAIL;
    }

    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close MMW dataset.\n");
      return FAIL;
    }
  }

  status = H5Fclose (file_id);

  if (my_cloudy->grid_rank > GRACKLE_CLOUDY_TABLE_MAX_DIMENSION) {
    fprintf(stderr,"Error: rank of Cloudy cooling data must be less than or equal to %d.\n",
	    GRACKLE_CLOUDY_TABLE_MAX_DIMENSION);
    return FAIL;
  }

  return SUCCESS;
}

int _free_cloudy_data(cloudy_data *my_cloudy, chemistry_data *my_chemistry, int primordial) {
  int i;

  for(i = 0; i < my_cloudy->grid_rank; i++) {
    GRACKLE_FREE(my_cloudy->grid_parameters[i]);
  }

  GRACKLE_FREE(my_cloudy->cooling_data);
  if (my_chemistry->UVbackground == 1) {
    GRACKLE_FREE(my_cloudy->heating_data);
  }
  if (my_chemistry->primordial_chemistry == 0 && primordial) {
    GRACKLE_FREE(my_cloudy->mmw_data);
  }
}
