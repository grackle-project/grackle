/***********************************************************************
/
/ Declare function used to run the bench subcommand
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef GRCLI_CMDBENCH_H
#define GRCLI_CMDBENCH_H

#include "CliParser.h"

namespace cmd::bench {

/// execute the bench command
[[noreturn]] void run(CliParser& parser);

} // namespace cmd::bench

#endif /* GRCLI_CMDBENCH_H */
