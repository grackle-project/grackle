#include <grackle.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "grackle_wrapper.h"

gr_float a_value;
gr_float temperature_units ;
gr_float velocity_units ;
gr_float cool_units ;
#define griddim 1
#define kboltz  1.3806504e-16
#define mass_h  1.67262171e-24
#define mass_e  9.10938215e-28
#define pi_val  3.14159265
#define hplanck 6.6260693e-27
#define ev2erg  1.60217653e-12
#define c_light 2.99792458e10
#define GravConst 6.6726e-8
#define sigma_sb  5.670373e-5
#define SolarMass 1.9891e33
#define Mpc       3.0857e24
#define kpc       3.0857e21
#define pc        3.0857e18
#define yr_to_s   3.15569e7

// arrays passed to grackle as input and to be filled
#define field_size  1
gr_float density[field_size];
gr_float internal_energy[field_size];
gr_float x_velocity[field_size];
gr_float y_velocity[field_size];
gr_float z_velocity[field_size];
 // for primordial_chemistry >= 1
gr_float HI_density[field_size];
gr_float HII_density[field_size];
gr_float HeI_density[field_size];
gr_float HeII_density[field_size];
gr_float HeIII_density[field_size];
gr_float ie_density[field_size];
 // for primordial_chemistry >= 2
gr_float HM_density[field_size];
gr_float H2I_density[field_size];
gr_float H2II_density[field_size];
 // for primordial_chemistry >= 3
gr_float DI_density[field_size];
gr_float DII_density[field_size];
gr_float HDI_density[field_size];
 // for metal_cooling = 1
gr_float metal_density[field_size];
// calculated
gr_float cooling_time[field_size];
gr_float outgamma[field_size];


// Set grid dimension and size.
// grid_start and grid_end are used to ignore ghost zones.
gr_int grid_rank = griddim;
// If grid rank is less than 3, set the other dimensions, 
// start indices, and end indices to 0.
gr_int grid_dimension[griddim], grid_start[griddim], grid_end[griddim];
static chemistry_data my_chemistry;
code_units my_units;

int wrap_init_cooling(double udensity, double ulength, double utime, 
                       int grackle_chemistry){
	/*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/
  int i;
  // First, set up the units system.
  // These are conversions from code units to cgs.
  my_units.a_units = 1. ;
  my_units.density_units = mass_h;
  my_units.length_units = 1.0;
  my_units.time_units = 1e12;
  velocity_units = my_units.a_units * my_units.length_units /
      my_units.time_units;
  temperature_units = mass_h * velocity_units*velocity_units / kboltz;

  // Second, create a chemistry object for parameters and rate data.
 my_chemistry  = set_default_chemistry_parameters();
  // Set parameter values for chemistry.
  my_chemistry.use_chemistry = 1;          
  my_chemistry.with_radiative_cooling = 1;
  my_chemistry.primordial_chemistry = grackle_chemistry;   // molecular network with H, He, D
  my_chemistry.metal_cooling = 1;          // metal cooling on
  my_chemistry.UVbackground = 1;          
  my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"; 


  // Finally, initialize the chemistry object.
  // snl: a_value is not the true initial ae!! This should get set during update_UVbackground
  gr_float initial_redshift = 100.; 
  a_value = 1. / (1. + initial_redshift);
  if (initialize_chemistry_data(my_chemistry, my_units, a_value) == 0) {
      return 0;
  }
  
  for (i = 0;i < griddim;i++) {
	  grid_dimension[i] = 0; 
	  grid_start[i] = 0;
	  grid_end[i] = 0;
  }
  grid_dimension[0] = field_size;
  grid_end[0] = field_size - 1;

  return 1;
}

int wrap_update_UVbackground_rates(double auni) {
	// The UV background rates must be updated before 
	// calling the other functions.
	/* a_value = auni / my_units.a_units; */
	if (update_UVbackground_rates(my_chemistry, 
	                              my_units, a_value) == 0) {
		return 0;
	}
	return 1;
}

int wrap_get_cooling_time( double rho_cgs, double tempK, double Zs, 
                           double vx_cgs, double vy_cgs, double vz_cgs, 
                           double csize, double auni, double dt_cgs, 
                           double *cooling_time_cgs,
                           double *gamma_cgs){
  int i;
  gr_float tiny_number = 1.0e-20;
  gr_float rho;
  gr_float dt;
  
  rho = rho_cgs/my_units.density_units;
  dt = dt_cgs/my_units.time_units;

  if(field_size != 1){fprintf(stderr,"field_size must currently be set to 1.\n");
      return 0;
  } 
  for (i = 0;i < field_size;i++) {
      density[i] =  rho; 
      HI_density[i] = my_chemistry.HydrogenFractionByMass * rho;
      HII_density[i] = tiny_number * rho;
      HM_density[i] = tiny_number * rho;
        
      HeI_density[i] = (1.0 - my_chemistry.HydrogenFractionByMass) * rho;
      HeII_density[i] = tiny_number * rho;
      HeIII_density[i] = tiny_number * rho;
      
      H2I_density[i] = tiny_number * rho;
      H2II_density[i] = tiny_number * rho;
      
      DI_density[i] = 2.0 * 3.4e-5 * rho;
      DII_density[i] = tiny_number * rho;
      HDI_density[i] = tiny_number * rho;
      
      ie_density[i] = tiny_number * rho; 
      metal_density[i] = Zs * rho;
      
      x_velocity[i] = vx_cgs;
      y_velocity[i] = vy_cgs;
      z_velocity[i] = vz_cgs;
        
      internal_energy[i] = tempK/temperature_units;
  }
  if (solve_chemistry(my_chemistry, my_units,
                      a_value, dt,
                      grid_rank, grid_dimension,
                      grid_start, grid_end,
                      density, internal_energy,
                      x_velocity, y_velocity, z_velocity,
                      HI_density, HII_density, HM_density,
                      HeI_density, HeII_density, HeIII_density,
                      H2I_density, H2II_density,
                      DI_density, DII_density, HDI_density,
                      ie_density, metal_density) == 0) {
      fprintf(stderr,"Error in solve_chemistry.\n");
      return 0;
  }
  if (calculate_cooling_time(my_chemistry, my_units,
                             a_value,
                             grid_rank, grid_dimension,
                             grid_start, grid_end,
                             density, internal_energy,
                             x_velocity, y_velocity, z_velocity,
                             HI_density, HII_density, HM_density,
                             HeI_density, HeII_density, HeIII_density,
                             H2I_density, H2II_density,
                             DI_density, DII_density, HDI_density,
                             ie_density, metal_density, 
                             cooling_time) == 0) {
      fprintf(stderr,"Error in calculate_cooling_time.\n");
      return 0;
  }
  if (my_chemistry.primordial_chemistry > 1){
      if (calculate_gamma(my_chemistry, my_units,
                          grid_rank, grid_dimension,
                          density, internal_energy,
                          HI_density, HII_density, HM_density,
                          HeI_density, HeII_density, HeIII_density,
                          H2I_density, H2II_density,
                          DI_density, DII_density, HDI_density,
                          ie_density, metal_density,
                          outgamma) == 0) {
          fprintf(stderr,"Error in calculate_gamma.\n");
          return 0;
      }
      (*gamma_cgs) = (double)outgamma[0];
      if(*gamma_cgs <1){
          fprintf(stderr,"bad gamma %f < 1 \n", *gamma_cgs);
          return 0;
      }
  }else{
	  (*gamma_cgs) = my_chemistry.Gamma;
  }
  (*cooling_time_cgs) = (double)cooling_time[0]*my_units.time_units;
  return 1;
}
