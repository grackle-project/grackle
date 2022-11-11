#include <stdio.h>
#include "grackle_types.h"

// the following macros are auto-generated:
#define AUTO_VERSION "3.2.dev2"
#define AUTO_BRANCH "genchiaki_merge"
#define AUTO_REVISION "29fe448d7726edbb31a079426adf2ebfadb61f32"

// test that ensures that all macros were correctly defined:
#if !(defined(AUTO_VERSION) && defined(AUTO_BRANCH) && defined(AUTO_REVISION))
#error "Something went wrong while auto-generating macros"
#endif

grackle_version get_grackle_version(void) {
  grackle_version out;
  out.version = AUTO_VERSION;
  out.branch = AUTO_BRANCH;
  out.revision = AUTO_REVISION;
  return out;
}
