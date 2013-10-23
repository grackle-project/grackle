/***********************************************************************
/
/  INITIALIZE CLOUDY DATA
/
/  written by: Britton Smith
/  date:       November, 2005
/  modified1:  May, 2009
/              Converted Cloudy table format from ascii to hdf5.
/
/  PURPOSE:  Read in heating, cooling, and mean molecular weight values 
/            from file.
/  Rank = 1: interpolate over temperature.
/  Rank = 2: interpolate over density and temperature.
/  Rank = 3: interpolate over density, redshift, and temperature.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h"

#define SMALL_LOG_VALUE -99.0
#define CLOUDY_COOLING_MAX_DIMENSION 3

/**************************** Functions Prototypes ******************************/

// Initialize Cloudy cooling data
int initialize_cloudy_data(chemistry_data &my_chemistry,
                           code_units &my_units, gr_float a_value)
{

  gr_int q, w;
  double *temp_data;
  long long temp_int;
  long long *temp_int_arr;
  char parameter_name[MAX_LINE_LENGTH];
  gr_int debug = 0;

  // Initialize things needed even if cloudy cooling is not used.

  my_chemistry.CloudyCoolingGridParameters = new gr_float*[CLOUDY_COOLING_MAX_DIMENSION];
  my_chemistry.CloudyCoolingGridDimension = new gr_int[CLOUDY_COOLING_MAX_DIMENSION];
  for (q = 0;q < CLOUDY_COOLING_MAX_DIMENSION;q++) {
    my_chemistry.CloudyCoolingGridDimension[q] = 0;
  }

  // Zero arrays if cloudy cooling not used.

  if (!my_chemistry.metal_cooling) {
    my_chemistry.CloudyCoolingGridRank = 0;
    return SUCCESS;
  }

  fprintf(stderr,"Initializing Cloudy cooling.\n");
  fprintf(stderr,"cloudy_table_file: %s.\n",my_chemistry.grackle_data_file);
  fprintf(stderr,"include_metal_heating: %"ISYM".\n",my_chemistry.include_metal_heating);
  fprintf(stderr,"cmb_temperature_floor: %"ISYM".\n",my_chemistry.cmb_temperature_floor);

  /* Get conversion units. */

  gr_float co_length_units, co_density_units;
  if (my_units.comoving_coordinates == TRUE) {
    co_length_units = my_units.length_units;
    co_density_units = my_units.density_units;
  }
  else {
    co_length_units = my_units.length_units *
      a_value * my_units.a_units;
    co_density_units = my_units.density_units /
      POW(a_value * my_units.a_units, 3);
  }

  double tbase1 = my_units.time_units;
  double xbase1 = co_length_units/(a_value * my_units.a_units);
  double dbase1 = co_density_units * POW(a_value * my_units.a_units, 3);
  double mh = 1.67e-24;
  double CoolUnit = (POW(my_units.a_units,5) * POW(xbase1,2) * POW(mh,2)) /
                    (POW(tbase1,3) * dbase1);

  // Read cooling data in from hdf5 file.

  hid_t       file_id, dset_id, attr_id; 
  herr_t      status;
  herr_t      h5_error = -1;

  fprintf(stderr,"Reading Cloudy data from %s.\n", 
          my_chemistry.grackle_data_file);
  file_id = H5Fopen(my_chemistry.grackle_data_file, 
                    H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open cooling dataset and get grid dimensions.

  dset_id =  H5Dopen(file_id, "/CloudyRates/Cooling");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open Cooling in %s.\n",my_chemistry.grackle_data_file);
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
  my_chemistry.CloudyCoolingGridRank = (int) temp_int;
  fprintf(stderr,"Cloudy cooling grid rank: %"ISYM".\n",my_chemistry.CloudyCoolingGridRank);
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Rank attribute in Cooling dataset.\n");
    return FAIL;
  }

  // Grid dimension.
  temp_int_arr = new long long[my_chemistry.CloudyCoolingGridRank];
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
  fprintf(stderr,"Cloudy cooling grid dimensions:");
  for (q = 0;q < my_chemistry.CloudyCoolingGridRank;q++) {
    my_chemistry.CloudyCoolingGridDimension[q] = (int) temp_int_arr[q];
    fprintf(stderr," %"ISYM,my_chemistry.CloudyCoolingGridDimension[q]);
  }
  fprintf(stderr,".\n");
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  delete [] temp_int_arr;

  // Grid parameters.
  for (q = 0;q < my_chemistry.CloudyCoolingGridRank;q++) {

    if (q < my_chemistry.CloudyCoolingGridRank - 1) {
      sprintf(parameter_name,"Parameter%"ISYM,(q+1));
    }
    else {
      sprintf(parameter_name,"Temperature");
    }

    temp_data = new double[my_chemistry.CloudyCoolingGridDimension[q]];

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

    my_chemistry.CloudyCoolingGridParameters[q] = 
      new gr_float[my_chemistry.CloudyCoolingGridDimension[q]];
    for (w = 0;w < my_chemistry.CloudyCoolingGridDimension[q];w++) {
      if (q < my_chemistry.CloudyCoolingGridRank - 1) {
	my_chemistry.CloudyCoolingGridParameters[q][w] = (float) temp_data[w];
      }
      else {
	// convert temeperature to log
	my_chemistry.CloudyCoolingGridParameters[q][w] = (float) log10(temp_data[w]);
      }

    }
    fprintf(stderr,"%s: %"GSYM" to %"GSYM" (%"ISYM" steps).\n",parameter_name,
            my_chemistry.CloudyCoolingGridParameters[q][0],
            my_chemistry.CloudyCoolingGridParameters[q][my_chemistry.CloudyCoolingGridDimension[q]-1],
            my_chemistry.CloudyCoolingGridDimension[q]);
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close %s attribute in Cooling dataset.\n",
              parameter_name);
      return FAIL;
    }
    delete [] temp_data;

  }

  // Read Cooling data.
  my_chemistry.CloudyDataSize = 1;
  for (q = 0;q < my_chemistry.CloudyCoolingGridRank;q++) {
    my_chemistry.CloudyDataSize *= my_chemistry.CloudyCoolingGridDimension[q];
  }
  temp_data = new double[my_chemistry.CloudyDataSize];

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
  fprintf(stderr,"Reading Cloudy Cooling dataset.\n");
  if (status == h5_error) {
    fprintf(stderr,"Failed to read Cooling dataset.\n");
    return FAIL;
  }

  my_chemistry.CloudyCooling = new gr_float[my_chemistry.CloudyDataSize];
  for (q = 0;q < my_chemistry.CloudyDataSize;q++) {
    my_chemistry.CloudyCooling[q] = temp_data[q] > 0 ? (float) log10(temp_data[q]) : (float) SMALL_LOG_VALUE;

    // Convert to code units.
    my_chemistry.CloudyCooling[q] -= log10(CoolUnit);
  }
  delete [] temp_data;

  status = H5Dclose(dset_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close Cooling dataset.\n");
    return FAIL;
  }

  // Read Heating data.
  if (my_chemistry.include_metal_heating) {

    temp_data = new double[my_chemistry.CloudyDataSize];

    dset_id =  H5Dopen(file_id, "/CloudyRates/Heating");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open Heating in %s.\n",my_chemistry.grackle_data_file);
      return FAIL;
    }

    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    fprintf(stderr,"Reading Cloudy Heating dataset.\n");
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Heating dataset.\n");
      return FAIL;
    }

    my_chemistry.CloudyHeating = new gr_float[my_chemistry.CloudyDataSize];
    for (q = 0;q < my_chemistry.CloudyDataSize;q++) {
      my_chemistry.CloudyHeating[q] = temp_data[q] > 0 ? (float) log10(temp_data[q]) : (float) SMALL_LOG_VALUE;

      // Convert to code units.
      my_chemistry.CloudyHeating[q] -= log10(CoolUnit);
    }
    delete [] temp_data;

    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Heating dataset.\n");
      return FAIL;
    }
  }

  status = H5Fclose (file_id);

  if (my_chemistry.CloudyCoolingGridRank > CLOUDY_COOLING_MAX_DIMENSION) {
    fprintf(stderr,"Error: rank of Cloudy cooling data must be less than or equal to %"ISYM".\n",
	    CLOUDY_COOLING_MAX_DIMENSION);
    return FAIL;
  }

  return SUCCESS;
}
