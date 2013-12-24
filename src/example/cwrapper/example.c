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

#include "grackle_wrapper.h" 

//How much chemistry will be used by Grackle for cooling? 1) HI, HII, HeI, HeII, HeIII, e; 2) H2, H2I+, H-; 3) D, D+, HD");
int grackle_chemistry=3;
double kboltz = 1.3806504e-16;
double mass_h = 1.67262171e-24;
double Myear = 3.14e13;
double Zsol = 0.02041;
                                    
void get_cooling_time( double rho_cgs, double tempK, double Zs, 
                       double vx_cgs, double vy_cgs, double vz_cgs, 
                       double dx_cell_cgs, double auni, double dt_cgs,
                       double *cooling_time,
                       double *gamma){
	if(wrap_get_cooling_time( rho_cgs, tempK, Zs, 
	                          vx_cgs, vy_cgs, vz_cgs, 
	                          dx_cell_cgs, 
	                          auni, dt_cgs, 
	                          cooling_time,
	                          gamma) != 1){
		fprintf(stderr,"get_cooling_time failed\n");
	    exit(-1);
    }
}
int main(){
	double a, rhogl, T_g, Z_met;
	double Tl_min, Tl_max, d_Tl;
	int nt;
	int it;
	FILE *output;
	char filename[128];
	double dlmin, dlmax, dld, Zlmin, Zlmax, dlZ;
	double rhog, cooling_time, gamma, Zl_met; 
	double vx_cgs=0, vy_cgs=0, vz_cgs=0; 
	double dx_cell_cgs = 1;
	double dt_cgs = Myear;

	/*NOTE: junk units are reset in the wrapper but they could be passed from here */
	double units_density=1e-99, units_length=1e-99, units_time=1e-99; 
	if(wrap_init_cooling(units_density, units_length, units_time, 
	                     grackle_chemistry) != 1){
		fprintf(stderr,"Error in initialize_chemistry_data.");
		exit(-1);
	}
	
	/* dlmin = log10(mass_h)-2.0;  */
	/* dlmax = log10(mass_h)+2.0;  */
	/* dld = (dlmax-dlmin)/8.0; */
	dlmin =log10(mass_h);
	dlmax =log10(mass_h);
	dld = 100;

	/* Zlmin = -3; */
	/* Zlmax = 1; */
	/* dlZ = (Zlmax-Zlmin)/8.0; */
	Zlmin = 0;
	Zlmax = 0;
	dlZ = 1;

	Tl_min = 1.0;
	Tl_max = 9.0;
	d_Tl = 0.1;

		
	for ( a = 1.0; a <= 1.0; a += 0.1 ) { /*NOTE: aexpn in UVbackground rates is currently RESET TO INITIAL VALUE in the wrapper */
		if(wrap_update_UVbackground_rates(a) != 1){
			fprintf(stderr,"Error in update_UBbackground_rates.\n");
			exit(-1);
		}
		
		sprintf(filename, "coolingrate_%6.4f.dat", a );
		output = fopen( filename, "w" );
		
		for ( rhogl = dlmin; rhogl <= dlmax; rhogl += dld ) {
			rhog = pow(10.0, rhogl);
			nt = (int)((Tl_max - Tl_min)/d_Tl);
			for ( it = 0; it < nt; it++ ) {
				T_g = pow(10.0, Tl_min + (float)(it)*d_Tl);
				for ( Zl_met = Zlmin; Zl_met <= Zlmax; Zl_met += dlZ ) {
					Z_met = pow(10.0, Zl_met)* Zsol ;
					get_cooling_time( rhog, T_g, Z_met, 
					                  vx_cgs, vy_cgs, vz_cgs, 
					                  dx_cell_cgs, a, dt_cgs,
					                  &cooling_time,
					                  &gamma);
					printf("coolingtime %e %e %e %e \n", rhog, T_g, Z_met, cooling_time);
					fprintf(output, "%e %e %e %e \n", rhog, T_g, Z_met, cooling_time);
					//ie = 1.5*kboltz*T_g/rhog*mass_h;
					//cr = ie/cooling_time; //erg/s*cm^3
				}
			}
		}
		fclose(output);
	}
	return 0;
}

