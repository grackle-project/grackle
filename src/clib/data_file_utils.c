/***********************************************************************
/
/ Implement logic used internally by Grackle to determine data files.
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
#include <stdlib.h>

#include "data_file_utils.h"

struct generic_file_props determine_data_file_(const char* grackle_data_file,
                                               int grackle_data_file_options)
{
  // initialize output struct in a format that will denote an error (if is
  // never modified)
  struct generic_file_props out = {NULL, 0, NULL, 0};

  if (grackle_data_file == NULL) {
    fprintf(stderr, "grackle_data_file must not be NULL\n");
    return out;

  } else if (grackle_data_file_options == -1) { // the legacy case!
    out.path = grackle_data_file;
    return out;

  } else {
    fprintf(stderr, "grackle_data_file_options has an unexpected value: %d\n",
            grackle_data_file_options);
    return out;

  }


}

void free_generic_file_props_(struct generic_file_props* ptr) {
  if (ptr != NULL) {
    if (ptr->path_requires_dealloc){
      free((char*)ptr->path);
    }
    ptr->path = NULL;
    if (ptr->checksum_requires_dealloc){
      free((char*)ptr->checksum);
    }
    ptr->checksum = NULL;
  }
}
