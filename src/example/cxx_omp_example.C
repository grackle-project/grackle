/***********************************************************************
/
/ Example executable using libgrackle
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
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "sys/time.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include <grackle.h>
}

#define mh     1.67262171e-24
#define kboltz 1.3806504e-16

// random fluctuation
#define FLUC( Amp )   ( ( (gr_float)rand()/RAND_MAX )*2.0*Amp - Amp + 1.0 )

// Grackle routines (TAB for tabulated cooling; CHE for non-equilibrium chemistry)
enum Routine_t { TAB_COOL=1, TAB_TIME=2, TAB_TEMP=3, TAB_PRES=4,
                 CHE_COOL=5, CHE_TIME=6, CHE_TEMP=7, CHE_PRES=8, CHE_GAMA=9 };

// function prototypes
float InvokeGrackle( const Routine_t Routine, const int NThread, const int NIter, const int NCell,
                     const int RSeed, const gr_float D0, const gr_float T0, const gr_float FAmp,
                     code_units *my_units, double a_value, double dt_value,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density, gr_float *output );
void CompareResult( const gr_float **Field_t1, const gr_float **Field_tN, const int NCell, const int NField );




//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :  Main function
//
// Note        :  Test different routines in Grackle using single and multiple OpenMP threads
//-------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

  // read the command-line options
# ifdef _OPENMP
  int NThread=omp_get_max_threads();
# else
  int NThread=1;
# endif
  int NIter=10, NCell1D=16, c;
  int Mode=3; // (1/2/3) => (single-thread/multi-thread/both)

  while ( (c = getopt(argc, argv, "ht:a:n:m:")) != -1 )
  {
    switch ( c )
    {
      case 't': NThread = atoi(optarg);   break;
      case 'a': NIter   = atoi(optarg);   break;
      case 'n': NCell1D = atoi(optarg);   break;
      case 'm': Mode    = atoi(optarg);   break;
      case 'h': 
      case '?': fprintf( stderr, "usage: %s\n", argv[0] );
                fprintf( stderr, "       %s\n", "[-h (for help)]" );
                fprintf( stderr, "       %s\n", "[-t number of OpenMP threads [max]]" );
                fprintf( stderr, "       %s\n", "[-a number of iterations [10]]" );
                fprintf( stderr, "       %s\n", "[-n number of cells along each direction [16]]" );
                fprintf( stderr, "       %s\n", "[-m (1/2/3) => (single-thread/multi-thread/both) [3]]" );
                exit( 1 );
    }
  }

  if ( NThread < 1 )
  {
    fprintf( stderr, "ERROR: incorrect input parameter [-t > 1]\n" );
    exit( EXIT_FAILURE );
  }

  if ( NIter < 1 )
  {
    fprintf( stderr, "ERROR: incorrect input parameter [-a > 1]\n" );
    exit( EXIT_FAILURE );
  }

  if ( NCell1D < 1 )
  {
    fprintf( stderr, "ERROR: incorrect input parameter [-n > 1]\n" );
    exit( EXIT_FAILURE );
  }

  if ( Mode < 1  ||  Mode > 3 )
  {
    fprintf( stderr, "ERROR: incorrect input parameter [-m 1/2/3]\n" );
    exit( EXIT_FAILURE );
  }


  char  FileName[100]  = "Table__OMP_Performance";
  const gr_float D0    = 1.0e0;   // background density (in my_units.density_units)
  const gr_float Temp0 = 1.0e3;   // background temperature (in K)
  const gr_float FAmp  = 0.1;     // perturbation amplitude
  const int      RSeed = 123;     // random seed

  timeval tv1, tv2;               // timing variables

  float T_CoolT_t1=0.0, T_CoolT_tN=0.0; // computation time for different routines
  float T_TimeT_t1=0.0, T_TimeT_tN=0.0;
  float T_TempT_t1=0.0, T_TempT_tN=0.0;
  float T_PresT_t1=0.0, T_PresT_tN=0.0;

  float T_CoolC_t1=0.0, T_CoolC_tN=0.0;
  float T_TimeC_t1=0.0, T_TimeC_tN=0.0;
  float T_TempC_t1=0.0, T_TempC_tN=0.0;
  float T_PresC_t1=0.0, T_PresC_tN=0.0;
  float T_GamaC_t1=0.0, T_GamaC_tN=0.0;


  // Disable output
  grackle_verbose = 0;

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0;    // 1 if cosmological sim, 0 if not
  my_units.density_units        = 1.67e-24;
  my_units.length_units         = 1.0;
  my_units.time_units           = 1.0e12;
  my_units.velocity_units       = my_units.length_units / my_units.time_units;
  my_units.a_units              = 1.0;  // units for the expansion factor

  // Second, create a chemistry object for parameters and rate data.
  if (set_default_chemistry_parameters() == 0) {
    fprintf( stderr, "Error in set_default_chemistry_parameters.\n" );
    return 0;
  }

  // Set parameter values for chemistry.
  grackle_data.use_grackle            = 1;          // chemistry on
  grackle_data.with_radiative_cooling = 1;          // cooling on
  grackle_data.primordial_chemistry   = 0;          // fully tabulated cooling
  grackle_data.metal_cooling          = 1;          // metal cooling on
  grackle_data.UVbackground           = 1;          // UV background on
  grackle_data.grackle_data_file      = "../../input/CloudyData_UVB=HM2012.h5"; // data file
# ifdef _OPENMP
  grackle_data.omp_nthreads           = NThread;    // number of OpenMP threads
                                                    // (remove this line if Grackle is not compiled with OpenMP support)
# endif

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  double initial_redshift = 0;
  double a_value = 1. / (1. + initial_redshift);

  // Finally, initialize the chemistry object.
  if ( initialize_chemistry_data(&my_units, a_value) == 0 ) {
    fprintf( stderr, "Error in initialize_chemistry_data.\n" );
    return EXIT_FAILURE;
  }

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  const int N3        = NCell1D*NCell1D*NCell1D;
  const int grid_rank = 3;

  // If grid rank is less than 3, set the other dimensions,
  // start indices, and end indices to 0.
  int grid_dimension[3], grid_start[3], grid_end[3];
  for (int i=0; i<3; i++)
  {
    grid_dimension[i] = NCell1D; // the active dimension not including ghost zones.
    grid_start    [i] = 0;
    grid_end      [i] = NCell1D-1;
  }

  // Allocate field arrays. (t1/tN for single-/multi-thread results)
  gr_float *density, *x_velocity, *y_velocity, *z_velocity, *metal_density;
  gr_float *energy_t1, *cooling_time_t1, *temperature_t1, *pressure_t1, *gamma_t1;
  gr_float *energy_tN, *cooling_time_tN, *temperature_tN, *pressure_tN, *gamma_tN;
  gr_float *HI_density_t1, *HII_density_t1, *HM_density_t1, *HeI_density_t1, *HeII_density_t1,
           *HeIII_density_t1, *H2I_density_t1, *H2II_density_t1, *DI_density_t1, *DII_density_t1,
           *HDI_density_t1, *e_density_t1;
  gr_float *HI_density_tN, *HII_density_tN, *HM_density_tN, *HeI_density_tN, *HeII_density_tN,
           *HeIII_density_tN, *H2I_density_tN, *H2II_density_tN, *DI_density_tN, *DII_density_tN,
           *HDI_density_tN, *e_density_tN;

  density          = new gr_float[N3];
  x_velocity       = new gr_float[N3];
  y_velocity       = new gr_float[N3];
  z_velocity       = new gr_float[N3];
  metal_density    = new gr_float[N3];

  energy_t1        = new gr_float[N3];
  cooling_time_t1  = new gr_float[N3];
  temperature_t1   = new gr_float[N3];
  pressure_t1      = new gr_float[N3];
  gamma_t1         = new gr_float[N3];

  energy_tN        = new gr_float[N3];
  cooling_time_tN  = new gr_float[N3];
  temperature_tN   = new gr_float[N3];
  pressure_tN      = new gr_float[N3];
  gamma_tN         = new gr_float[N3];

  HI_density_t1    = new gr_float[N3];
  HII_density_t1   = new gr_float[N3];
  HeI_density_t1   = new gr_float[N3];
  HeII_density_t1  = new gr_float[N3];
  HeIII_density_t1 = new gr_float[N3];
  e_density_t1     = new gr_float[N3];
  HM_density_t1    = new gr_float[N3];
  H2I_density_t1   = new gr_float[N3];
  H2II_density_t1  = new gr_float[N3];
  DI_density_t1    = new gr_float[N3];
  DII_density_t1   = new gr_float[N3];
  HDI_density_t1   = new gr_float[N3];

  HI_density_tN    = new gr_float[N3];
  HII_density_tN   = new gr_float[N3];
  HeI_density_tN   = new gr_float[N3];
  HeII_density_tN  = new gr_float[N3];
  HeIII_density_tN = new gr_float[N3];
  e_density_tN     = new gr_float[N3];
  HM_density_tN    = new gr_float[N3];
  H2I_density_tN   = new gr_float[N3];
  H2II_density_tN  = new gr_float[N3];
  DI_density_tN    = new gr_float[N3];
  DII_density_tN   = new gr_float[N3];
  HDI_density_tN   = new gr_float[N3];

  // set temperature units and other constants
  const double temperature_units = mh * pow(my_units.a_units *
                                            my_units.length_units /
                                            my_units.time_units, 2) / kboltz;
  const gr_float T0 = Temp0 / temperature_units;
  const double   dt = 3.15e7 * 1e6 / my_units.time_units;

  // initialize random seed
  srand( RSeed );
  
  // initialize input fields
  // fields to be updated are reset by InvokeGrackle to ensure the same input fields
  for (int i=0; i<N3; i++)
  {
    density      [i] = D0*FLUC(FAmp);
    metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i] * FLUC(FAmp);
    x_velocity   [i] = 0.0; // velocity are not useful currently
    y_velocity   [i] = 0.0;
    z_velocity   [i] = 0.0;
  }


  /*********************************************************************
  / 1. Tabulated cooling 
  *********************************************************************/

  // 1-1. Evolve tabulated cooling for dt
  fprintf( stdout, "Calculating cooling (tabulated) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_CoolT_t1
      = InvokeGrackle( TAB_COOL, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, NULL );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_CoolT_tN
      = InvokeGrackle( TAB_COOL, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, NULL );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { energy_t1 };
    const gr_float *Field_tN[NField] = { energy_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 1-2. Calculate cooling time
  fprintf( stdout, "Calculating cooling time (tabulated) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_TimeT_t1
      = InvokeGrackle( TAB_TIME, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, cooling_time_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_TimeT_tN
      = InvokeGrackle( TAB_TIME, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, cooling_time_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { cooling_time_t1 };
    const gr_float *Field_tN[NField] = { cooling_time_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 1-3. Calculate temperature
  fprintf( stdout, "Calculating temperature (tabulated) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_TempT_t1
      = InvokeGrackle( TAB_TEMP, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, temperature_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_TempT_tN
      = InvokeGrackle( TAB_TEMP, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, temperature_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { temperature_t1 };
    const gr_float *Field_tN[NField] = { temperature_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 1-4. Calculate pressure
  fprintf( stdout, "Calculating pressure (tabulated) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_PresT_t1
      = InvokeGrackle( TAB_PRES, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, pressure_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_PresT_tN
      = InvokeGrackle( TAB_PRES, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, pressure_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { pressure_t1 };
    const gr_float *Field_tN[NField] = { pressure_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }



  /*********************************************************************
  / 2. Non-equilibrium chemistry network
  *********************************************************************/

  // Set parameter values for chemistry and initialize the chemistry object again
  grackle_data.primordial_chemistry = 3;  // molecular network with H, He, D

  if ( initialize_chemistry_data(&my_units, a_value) == 0 ) {
    fprintf( stderr, "Error in initialize_chemistry_data.\n" );
    return EXIT_FAILURE;
  }

  // 2-1. Evolve chemistry for dt
  fprintf( stdout, "Calculating cooling (non-equilibrium chemistry) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_CoolC_t1
      = InvokeGrackle( CHE_COOL, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, NULL );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_CoolC_tN
      = InvokeGrackle( CHE_COOL, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, NULL );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 13;
    const gr_float *Field_t1[NField] =
     { HI_density_t1, HeI_density_t1, HII_density_t1, HM_density_t1,
       HeII_density_t1, HeIII_density_t1, H2I_density_t1, H2II_density_t1,
       DI_density_t1, DII_density_t1, HDI_density_t1, e_density_t1, energy_t1 };
    const gr_float *Field_tN[NField] =
     { HI_density_tN, HeI_density_tN, HII_density_tN, HM_density_tN,
       HeII_density_tN, HeIII_density_tN, H2I_density_tN, H2II_density_tN,
       DI_density_tN, DII_density_tN, HDI_density_tN, e_density_tN, energy_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 2-2. Calculate cooling time
  fprintf( stdout, "Calculating cooling time (non-equilibrium chemistry) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_TimeC_t1
      = InvokeGrackle( CHE_TIME, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, cooling_time_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_TimeC_tN
      = InvokeGrackle( CHE_TIME, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, cooling_time_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { cooling_time_t1 };
    const gr_float *Field_tN[NField] = { cooling_time_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 2-3. Calculate temperature
  fprintf( stdout, "Calculating temperature (non-equilibrium chemistry) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_TempC_t1
      = InvokeGrackle( CHE_TEMP, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, temperature_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_TempC_tN
      = InvokeGrackle( CHE_TEMP, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, temperature_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { temperature_t1 };
    const gr_float *Field_tN[NField] = { temperature_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 2-4. Calculate pressure
  fprintf( stdout, "Calculating pressure (non-equilibrium chemistry) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_PresC_t1
      = InvokeGrackle( CHE_PRES, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, pressure_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_PresC_tN
      = InvokeGrackle( CHE_PRES, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, pressure_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { pressure_t1 };
    const gr_float *Field_tN[NField] = { pressure_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }


  // 2-5. Calculate gamma 
  fprintf( stdout, "Calculating gamma (non-equilibrium chemistry) ...\n" ); fflush( stdout );

  // single-thread
  if ( Mode != 2 )
  {
    T_GamaC_t1
      = InvokeGrackle( CHE_GAMA, 1, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_t1, x_velocity, y_velocity, z_velocity,
                       HI_density_t1, HII_density_t1, HM_density_t1,
                       HeI_density_t1, HeII_density_t1, HeIII_density_t1,
                       H2I_density_t1, H2II_density_t1,
                       DI_density_t1, DII_density_t1, HDI_density_t1,
                       e_density_t1, metal_density, gamma_t1 );
  }

  // multi-thread
  if ( Mode != 1 )
  {
    T_GamaC_tN
      = InvokeGrackle( CHE_GAMA, NThread, NIter, N3, RSeed, D0, T0, FAmp,
                       &my_units, a_value, dt,
                       grid_rank, grid_dimension, grid_start, grid_end,
                       density, energy_tN, x_velocity, y_velocity, z_velocity,
                       HI_density_tN, HII_density_tN, HM_density_tN,
                       HeI_density_tN, HeII_density_tN, HeIII_density_tN,
                       H2I_density_tN, H2II_density_tN,
                       DI_density_tN, DII_density_tN, HDI_density_tN,
                       e_density_tN, metal_density, gamma_tN );
  }

  // check results
  if ( Mode == 3 )
  {
    const int NField = 1;
    const gr_float *Field_t1[NField] = { gamma_t1 };
    const gr_float *Field_tN[NField] = { gamma_tN };

    CompareResult( Field_t1, Field_tN, N3, NField );
  }



  // record header
  if ( fopen(FileName,"r") == NULL )
  {
    FILE *File = fopen( FileName, "w" );

    fprintf( File, "#NC          = # of cells in each grid\n" );
    fprintf( File, "#NI          = # of iterations for estimating the average performance\n" );
    fprintf( File, "#NT          = # of OpenMP threads\n" );
    fprintf( File, "\n" );
    fprintf( File, "#Cool        = cooling\n" );
    fprintf( File, "#Time        = cooling time\n" );
    fprintf( File, "#Temp        = temperature\n" );
    fprintf( File, "#Pres        = pressure\n" );
    fprintf( File, "#Gama        = ratio of specific heat (gamma)\n" );
    fprintf( File, "#Ratio       = multi-thread performance / single-thread performance\n" );
    fprintf( File, "\n" );
    fprintf( File, "#T/C         = for tabulated / 12-species chemistry network\n" );
    fprintf( File, "#_t1/_tN     = single-/multi-thread performance in [# of cells / sec]\n" );
    fprintf( File, "\n" );
    fprintf( File, "#Precision   = %s\n", (sizeof(gr_float)==4)?"single":"double" );
    fprintf( File, "#Density     = %13.7e g/cm^3\n", D0*my_units.density_units );
    fprintf( File, "#Temperature = %13.7e K\n",      Temp0 );
    fprintf( File, "#====================================================================================================\n" );
    fprintf( File, "\n" );
    fprintf( File, "#%6s  %3s  %3s", "NC", "NI", "NT" );
    fprintf( File, "  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s",
             "CoolT_t1", "CoolT_tN", "Ratio", "TimeT_t1", "TimeT_tN", "Ratio",
             "TempT_t1", "TempT_tN", "Ratio", "PresT_t1", "PresT_tN", "Ratio" );
    fprintf( File, "  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s",
             "CoolC_t1", "CoolC_tN", "Ratio", "TimeC_t1", "TimeC_tN", "Ratio",
             "TempC_t1", "TempC_tN", "Ratio", "PresC_t1", "PresC_tN", "Ratio" );
    fprintf( File, "  %8s  %8s  %8s",
             "GamaC_t1", "GamaC_tN", "Ratio" );
    fprintf( File, "\n" );

    fclose( File );
  }

  // record performance
  FILE *File = fopen( FileName, "a" );
  fprintf( File, "%7d  %3d  %3d", N3, NIter, NThread );
  fprintf( File, "  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e",
           N3/T_CoolT_t1, N3/T_CoolT_tN, T_CoolT_t1/T_CoolT_tN,
           N3/T_TimeT_t1, N3/T_TimeT_tN, T_TimeT_t1/T_TimeT_tN,
           N3/T_TempT_t1, N3/T_TempT_tN, T_TempT_t1/T_TempT_tN,
           N3/T_PresT_t1, N3/T_PresT_tN, T_PresT_t1/T_PresT_tN );
  fprintf( File, "  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e",
           N3/T_CoolC_t1, N3/T_CoolC_tN, T_CoolC_t1/T_CoolC_tN,
           N3/T_TimeC_t1, N3/T_TimeC_tN, T_TimeC_t1/T_TimeC_tN,
           N3/T_TempC_t1, N3/T_TempC_tN, T_TempC_t1/T_TempC_tN,
           N3/T_PresC_t1, N3/T_PresC_tN, T_PresC_t1/T_PresC_tN );
  fprintf( File, "  %8.2e  %8.2e  %8.2e",
           N3/T_GamaC_t1, N3/T_GamaC_tN, T_GamaC_t1/T_GamaC_tN );
  fprintf( File, "\n" );
  fclose( File );


  delete [] density;
  delete [] x_velocity;
  delete [] y_velocity;
  delete [] z_velocity;
  delete [] metal_density;

  delete [] energy_t1;
  delete [] cooling_time_t1;
  delete [] temperature_t1;
  delete [] pressure_t1;
  delete [] gamma_t1;

  delete [] energy_tN;
  delete [] cooling_time_tN;
  delete [] temperature_tN;
  delete [] pressure_tN;
  delete [] gamma_tN;

  delete [] HI_density_t1;
  delete [] HII_density_t1;
  delete [] HeI_density_t1;
  delete [] HeII_density_t1;
  delete [] HeIII_density_t1;
  delete [] e_density_t1;
  delete [] HM_density_t1;
  delete [] H2I_density_t1;
  delete [] H2II_density_t1;
  delete [] DI_density_t1;
  delete [] DII_density_t1;
  delete [] HDI_density_t1;

  delete [] HI_density_tN;
  delete [] HII_density_tN;
  delete [] HeI_density_tN;
  delete [] HeII_density_tN;
  delete [] HeIII_density_tN;
  delete [] e_density_tN;
  delete [] HM_density_tN;
  delete [] H2I_density_tN;
  delete [] H2II_density_tN;
  delete [] DI_density_tN;
  delete [] DII_density_tN;
  delete [] HDI_density_tN;

  return EXIT_SUCCESS;

} // FUNCTION : main



//-------------------------------------------------------------------------------------------------------
// Function    :  InvokeGrackle
// Description :  Invoke different routines in Grackle
//
// Note        :  1. The pointer "output" is used to store the updated
//                   cooling_time/pressure/temperature/gamma
//                2. Several input fields are reinitialized here to ensure they are the same
//                   when comparing single- and multi-threads results
//
// Return      :  Elapsed time
//-------------------------------------------------------------------------------------------------------
float InvokeGrackle( const Routine_t Routine, const int NThread, const int NIter, const int NCell,
                     const int RSeed, const gr_float D0, const gr_float T0, const gr_float FAmp,
                     code_units *my_units, double a_value, double dt,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density, gr_float *output )
{

  fprintf( stdout, "   Number of threads = %2d ... ", NThread ); fflush( stdout );

  const gr_float tiny_number = 1.e-20;
  timeval tv1, tv2;   // timing variables
  int     ChkErr;

  // set the number of OpenMP thread
  // Note that for real applications one does not need to invoke "omp_set_num_threads"
  // explicitly. Instead one should set "grackle_data.omp_nthreads" and then call
  // "initialize_chemistry_data".
# ifdef _OPENMP
  omp_set_num_threads( NThread );
# endif

  // reset random seed to ensure the same input fields
  srand( RSeed+(int)Routine );

  gettimeofday( &tv1, NULL );

  for (int t=0; t<NIter; t++)
  {
    switch ( Routine )
    {
      case TAB_COOL:
        // initialize density to fix the round-off errors for comoving coord.
        if ( my_units->comoving_coordinates == 1 )
        for (int i=0; i<NCell; i++)
        {
          density      [i] = D0*FLUC(FAmp);
          metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i] * FLUC(FAmp);
        }

        // initialize the fields to be updated
        for (int i=0; i<NCell; i++) energy[i] = T0*FLUC(FAmp);

        ChkErr = solve_chemistry_table(
                   my_units, a_value, dt,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   x_velocity, y_velocity, z_velocity,
                   metal_density );
        break;

      case TAB_TIME:
        // initialize density to fix the round-off errors for comoving coord.
        if ( my_units->comoving_coordinates == 1 )
        for (int i=0; i<NCell; i++)
        {
          density      [i] = D0*FLUC(FAmp);
          metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i] * FLUC(FAmp);
        }

        ChkErr = calculate_cooling_time_table(
                   my_units, a_value,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   x_velocity, y_velocity, z_velocity,
                   metal_density, output );
        break;

      case TAB_TEMP:
        // initialize density to fix the round-off errors for comoving coord.
        if ( my_units->comoving_coordinates == 1 )
        for (int i=0; i<NCell; i++)
        {
          density      [i] = D0*FLUC(FAmp);
          metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i] * FLUC(FAmp);
        }

        ChkErr = calculate_temperature_table(
                   my_units, a_value,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   metal_density, output );
        break;

      case TAB_PRES:
        ChkErr =  calculate_pressure_table(
                    my_units, a_value,
                    grid_rank, grid_dimension,
                    grid_start, grid_end,
                    density, energy,
                    output );
        break;

      case CHE_COOL:
        // initialize density to fix the round-off errors for comoving coord.
        if ( my_units->comoving_coordinates == 1 )
        for (int i=0; i<NCell; i++)
        {
          density      [i] = D0*FLUC(FAmp);
          metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i] * FLUC(FAmp);
        }

        // initialize the fields to be updated
        for (int i=0; i<NCell; i++)
        {
          HI_density   [i] = grackle_data.HydrogenFractionByMass * density[i];
          HeI_density  [i] = density[i] - HI_density[i];
          HII_density  [i] = tiny_number * density[i];
          HM_density   [i] = tiny_number * density[i];
          HeII_density [i] = tiny_number * density[i];
          HeIII_density[i] = tiny_number * density[i];
          H2I_density  [i] = tiny_number * density[i];
          H2II_density [i] = tiny_number * density[i];
          DI_density   [i] = 2.0*3.4e-5  * density[i];
          DII_density  [i] = tiny_number * density[i];
          HDI_density  [i] = tiny_number * density[i];
          e_density    [i] = tiny_number * density[i];
          energy       [i] = T0*FLUC(FAmp);
        }

        ChkErr = solve_chemistry(
                   my_units, a_value, dt,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   x_velocity, y_velocity, z_velocity,
                   HI_density, HII_density, HM_density,
                   HeI_density, HeII_density, HeIII_density,
                   H2I_density, H2II_density,
                   DI_density, DII_density, HDI_density,
                   e_density, metal_density );
        break;

      case CHE_TIME:
        // initialize density to fix the round-off errors for comoving coord.
        if ( my_units->comoving_coordinates == 1 )
        for (int i=0; i<NCell; i++)
        {
          density      [i] = D0*FLUC(FAmp);
          metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i] * FLUC(FAmp);
        }

        ChkErr = calculate_cooling_time(
                   my_units, a_value,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   x_velocity, y_velocity, z_velocity,
                   HI_density, HII_density, HM_density,
                   HeI_density, HeII_density, HeIII_density,
                   H2I_density, H2II_density,
                   DI_density, DII_density, HDI_density,
                   e_density, metal_density, output );
        break;

      case CHE_TEMP:
        ChkErr = calculate_temperature(
                   my_units, a_value,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   HI_density, HII_density, HM_density,
                   HeI_density, HeII_density, HeIII_density,
                   H2I_density, H2II_density,
                   DI_density, DII_density, HDI_density,
                   e_density, metal_density, output );
        break;

      case CHE_PRES:
        ChkErr = calculate_pressure(
                   my_units, a_value,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   HI_density, HII_density, HM_density,
                   HeI_density, HeII_density, HeIII_density,
                   H2I_density, H2II_density,
                   DI_density, DII_density, HDI_density,
                   e_density, metal_density, output );
        break;

      case CHE_GAMA:
        ChkErr = calculate_gamma(
                   my_units, a_value,
                   grid_rank, grid_dimension,
                   grid_start, grid_end,
                   density, energy,
                   HI_density, HII_density, HM_density,
                   HeI_density, HeII_density, HeIII_density,
                   H2I_density, H2II_density,
                   DI_density, DII_density, HDI_density,
                   e_density, metal_density, output );
        break;

      default:
        fprintf( stderr, "ERROR: unknown Grackle routine (%d) !!\n", Routine );
        exit( EXIT_FAILURE );
    } // switch ( Routine )

    if ( ChkErr == 0 )
    {
      fprintf( stderr, "ERROR: routine %u !!\n", Routine );
      return EXIT_FAILURE;
    }
  } // for (int t=0; t<NIter; t++)

  gettimeofday( &tv2, NULL );

  fprintf( stdout, "done\n" ); fflush( stdout );

  // return the average elapsed time in seconds
  return ( ( tv2.tv_sec*1000000 + tv2.tv_usec )
          -( tv1.tv_sec*1000000 + tv1.tv_usec ) )*1.0e-6/(float)NIter;

} // FUNCTION : InvokeGrackle



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareResult
// Description :  Compare the single- and multi-threads results
//
// Note        :  1. We expect even the round-off errors are the same
//                2. Program will be terminated if any inconsistency is found
//-------------------------------------------------------------------------------------------------------
void CompareResult( const gr_float **Field_t1, const gr_float **Field_tN, const int NCell, const int NField )
{

    fprintf( stdout, "   Checking results ... " ); fflush( stdout );

    for (int v=0; v<NField; v++)
    for (int i=0; i<NCell; i++)
    {
      // we expect even the round-off errors are the same
      if ( Field_t1[v][i] != Field_tN[v][i] )
      {
        fprintf( stderr, "\nERROR: inconsistent results for field %d, cell %d\n", v, i );
        fprintf( stderr, "   single-thread %21.14e, multi-thread %21.14e, error %21.14e !!\n\n",
                 Field_t1[v][i], Field_tN[v][i], fabs(Field_tN[v][i]-Field_t1[v][i])/fabs(Field_t1[v][i]) );
        exit( EXIT_FAILURE );
      }
    }

    fprintf( stdout, "passed\n" ); fflush( stdout );

} // FUNCTION : CompareResult
