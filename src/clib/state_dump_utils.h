/***********************************************************************
/
/ Declare utility functions defined for dumping Grackle's internal state (for
/ debugging purposes)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdio.h>

/// write each member of my_chemistry to fp
void show_parameters_(FILE *fp, const chemistry_data *my_chemistry);

/// write each member of gversion to fp
void show_version_(FILE *fp, const grackle_version* gversion);
