/***********************************************************************
/
/  INITIALIZE UV BACKGROUND DATA
/
/  written by: Michael Kuhlen
/  date:       October, 2012
/  modified1:  
/
/  PURPOSE:  Read in photo-ionization and photo-heating rates from
/            a user-specified HDF5 file.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"

// function prototypes
gr_int read_dataset(hid_t file_id, char *dset_name, gr_float *buffer);


// Initialize UV Background data
int initialize_UVbackground_data(chemistry_data &my_chemistry,
				 code_units &my_units, gr_float a_value)
{
  gr_int Nz, i;

  // Return if no UV background selected.
  if (my_chemistry.UVbackground == 0)
    return SUCCESS;


  fprintf(stderr,"Initializing UV background.\n");


  // Read in UV background data from hdf5 file.

  hid_t       file_id, dset_id, dspace_id;
  herr_t      status;
  herr_t      h5_error = -1;

  fprintf(stderr,"Reading UV background data from %s.\n", 
          my_chemistry.UVbackground_file);
  file_id = H5Fopen(my_chemistry.UVbackground_file, 
                    H5F_ACC_RDONLY, H5P_DEFAULT);


  // Read Info dataset

  dset_id =  H5Dopen(file_id, "/Info");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open 'Info' dataset in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  int strlen = (int)(H5Dget_storage_size(dset_id));
  char info_string[strlen+1];

  hid_t memtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(memtype, strlen+1);

  status = H5Dread(dset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, info_string);
  if (status == h5_error) {
    fprintf(stderr,"Failed to read dataset 'Info'.\n");
    return FAIL;
  }

  H5Tclose(memtype); 
  H5Dclose(dset_id);



  // Open redshift dataset and get number of elements

  dset_id =  H5Dopen(file_id, "/z");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open redshift dataset ('z') in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error) {
    fprintf(stderr,"Error opening dataspace for dataset 'z' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  Nz = H5Sget_simple_extent_npoints(dspace_id);
  if(Nz <= 0) {
    fprintf(stderr,"Redshift dataset ('z') has inappropriate size = %lld in %s.\n",Nz,my_chemistry.UVbackground_file);
    return FAIL;
  }

  H5Sclose(dspace_id);
  H5Dclose(dset_id);



  // Now allocate memory for UV background table.
  my_chemistry.UVbackground_table.Nz = Nz;

  my_chemistry.UVbackground_table.z = new gr_float[Nz];
  my_chemistry.UVbackground_table.k24 = new gr_float[Nz];
  my_chemistry.UVbackground_table.k25 = new gr_float[Nz];
  my_chemistry.UVbackground_table.k26 = new gr_float[Nz];

  if (my_chemistry.primordial_chemistry > 1) {
    my_chemistry.UVbackground_table.k27 = new gr_float[Nz];
    my_chemistry.UVbackground_table.k28 = new gr_float[Nz];
    my_chemistry.UVbackground_table.k29 = new gr_float[Nz];
    my_chemistry.UVbackground_table.k30 = new gr_float[Nz];
    my_chemistry.UVbackground_table.k31 = new gr_float[Nz];
  }    

  my_chemistry.UVbackground_table.piHI = new gr_float[Nz];
  my_chemistry.UVbackground_table.piHeII = new gr_float[Nz];
  my_chemistry.UVbackground_table.piHeI = new gr_float[Nz];


  // Now read everything.


  // *** Redshift ***
  if(! read_dataset(file_id, "z", my_chemistry.UVbackground_table.z) ) {
    fprintf(stderr,"Error reading dataset 'z' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  // *** k24 ***
  if(! read_dataset(file_id, "k24", my_chemistry.UVbackground_table.k24) ) {
    fprintf(stderr,"Error reading dataset 'k24' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  // *** k25 ***
  if(! read_dataset(file_id, "k25", my_chemistry.UVbackground_table.k25) ) {
    fprintf(stderr,"Error reading dataset 'k25' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  // *** k26 ***
  if(! read_dataset(file_id, "k26", my_chemistry.UVbackground_table.k26) ) {
    fprintf(stderr,"Error reading dataset 'k26' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  if (my_chemistry.primordial_chemistry > 1) {

    // *** k27 ***
    if(! read_dataset(file_id, "k27", my_chemistry.UVbackground_table.k27) ) {
      fprintf(stderr,"Error reading dataset 'k27' in %s.\n",my_chemistry.UVbackground_file);
      return FAIL;      
    }

    // *** k28 ***
    if(! read_dataset(file_id, "k28", my_chemistry.UVbackground_table.k28) ) {
      fprintf(stderr,"Error reading dataset 'k28' in %s.\n",my_chemistry.UVbackground_file);
      return FAIL;      
    }

    // *** k29 ***
    if(! read_dataset(file_id, "k29", my_chemistry.UVbackground_table.k29) ) {
      fprintf(stderr,"Error reading dataset 'k29' in %s.\n",my_chemistry.UVbackground_file);
      return FAIL;      
    }

    // *** k30 ***
    if(! read_dataset(file_id, "k30", my_chemistry.UVbackground_table.k30) ) {
      fprintf(stderr,"Error reading dataset 'k30' in %s.\n",my_chemistry.UVbackground_file);
      return FAIL;      
    }

    // *** k31 ***
    if(! read_dataset(file_id, "k31", my_chemistry.UVbackground_table.k31) ) {
      fprintf(stderr,"Error reading dataset 'k31' in %s.\n",my_chemistry.UVbackground_file);
      return FAIL;      
    }
    
  }

  // *** piHI ***
  if(! read_dataset(file_id, "piHI", my_chemistry.UVbackground_table.piHI) ) {
    fprintf(stderr,"Error reading dataset 'piHI' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  // *** piHeII ***
  if(! read_dataset(file_id, "piHeII", my_chemistry.UVbackground_table.piHeII) ) {
    fprintf(stderr,"Error reading dataset 'piHeII' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  // *** piHeI ***
  if(! read_dataset(file_id, "piHeI", my_chemistry.UVbackground_table.piHeI) ) {
    fprintf(stderr,"Error reading dataset 'piHeI' in %s.\n",my_chemistry.UVbackground_file);
    return FAIL;
  }

  
  H5Fclose(file_id);


  // Now convert the rates to code units.

  /* Get conversion units. */
  double tbase1 = my_units.time_units;
  double xbase1 = my_units.length_units/(a_value * my_units.a_units);
  double dbase1 = my_units.density_units*POW(a_value * my_units.a_units, 3);
  double mh     = 1.67e-24;
  double CoolingUnits = (POW(my_units.a_units, 5) * xbase1*xbase1 * mh*mh) /
    (POW(tbase1, 3) * dbase1);
  
  for(i=0;i<Nz;i++) {
    my_chemistry.UVbackground_table.k24[i] *= my_units.time_units;
    my_chemistry.UVbackground_table.k25[i] *= my_units.time_units;
    my_chemistry.UVbackground_table.k26[i] *= my_units.time_units;

    if (my_chemistry.primordial_chemistry > 1) {
      my_chemistry.UVbackground_table.k27[i] *= my_units.time_units;
      my_chemistry.UVbackground_table.k28[i] *= my_units.time_units;
      my_chemistry.UVbackground_table.k29[i] *= my_units.time_units;
      my_chemistry.UVbackground_table.k30[i] *= my_units.time_units;
      my_chemistry.UVbackground_table.k31[i] *= my_units.time_units;
    }

    my_chemistry.UVbackground_table.piHI[i] /= CoolingUnits;
    my_chemistry.UVbackground_table.piHeII[i] /= CoolingUnits;
    my_chemistry.UVbackground_table.piHeI[i] /= CoolingUnits;
  }


  // Get min/max of redshift vector
  my_chemistry.UVbackground_table.zmin = my_chemistry.UVbackground_table.z[0];
  my_chemistry.UVbackground_table.zmax = my_chemistry.UVbackground_table.z[Nz-1];

  // Print out some information about the dataset just read in.
  printf("UV background information:\n");
  printf("  %s\n",info_string);
  printf("  z_min = %6.3f\n  z_max = %6.3f\n",my_chemistry.UVbackground_table.zmin,my_chemistry.UVbackground_table.zmax);


  return SUCCESS;
}



gr_int read_dataset(hid_t file_id, char *dset_name, gr_float *buffer) {
  hid_t dset_id;
  herr_t status;
  herr_t h5_error = -1;

  dset_id =  H5Dopen(file_id, dset_name);
  if (dset_id == h5_error) {
    fprintf(stderr,"Failed to open dataset 'z'.\n");
    return FAIL;
  }

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  if (status == h5_error) {
    fprintf(stderr,"Failed to read dataset 'z'.\n");
    return FAIL;
  }
 
  H5Dclose(dset_id);

  return SUCCESS;
}
