/***********************************************************************
/
/ Declare functions used internally for formatting data (across routines)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef GRACKLE_SHOW_H
#define GRACKLE_SHOW_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "grackle.h"

void show_parameters(FILE *fp, chemistry_data *my_chemistry);
void show_version(FILE *fp);

static inline void show_units(FILE *fp, code_units* my_units) {
  fprintf(fp, "comoving_coordinates = %d\n", my_units->comoving_coordinates);
  fprintf(fp, "density_units        = %g\n", my_units->density_units);
  fprintf(fp, "length_units         = %g\n", my_units->length_units);
  fprintf(fp, "time_units           = %g\n", my_units->time_units);
  fprintf(fp, "a_units              = %g\n", my_units->a_units);
  fprintf(fp, "a_value              = %g\n", my_units->a_value);
}


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* GRACKLE_SHOW_H */
