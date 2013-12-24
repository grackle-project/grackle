/***********************************************************************
/
/ Grackle c wrapper
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __WRAPGRACKLE_H__
#define __WRAPGRACKLE_H__

#ifdef __cplusplus
extern "C" {
#endif 
int wrap_init_cooling(double udensity, double ulength, double utime, 
                       int grackle_chemistry);
int wrap_update_UVbackground_rates(double auni);
int wrap_get_cooling_time(double rho_cgs, double ie_density_cgs, double Zs, 
                           double vx_cgs, double vy_cgs, double vz_cgs, 
                           double size, double auni, double dt_cgs, 
                           double *cooling_time,
                           double *gamma_cgs);
#ifdef __cplusplus
}
#endif
    

#endif /*__WRAPGRACKLE_H__*/
