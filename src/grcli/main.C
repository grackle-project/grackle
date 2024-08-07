#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string_view>

#include "CliParser.h"
#include "cmd_bench.h"
#include "cmd_show.h"
#include "utils.h"

static void print_name_section_(FILE* f, const char* bin_name) {
  std::fprintf(f,
      R"HELP(NAME:

    %s  - A command line interface for benchmarking and debugging Grackle
)HELP",
      bin_name);
}


static void print_synopsis_section_(FILE* f, const char* bin_name) {
  const char* pre_subcommand_opts = "[--version] [--precision] [-h | --help]";

  std::fprintf(f,
      R"HELP(
SYNOPSIS:

    %s %s bench -S... PARAMETER_OPTIONS -O...

    %s %s show-parameters PARAMETER_OPTIONS

    %s %s show-initial-units
)HELP",
      bin_name, pre_subcommand_opts, bin_name, pre_subcommand_opts,
      bin_name, pre_subcommand_opts);
}

static const char* help_message_ = R"HELP(
DESCRIPTION:
    A simple command-line tool intended for debugging and basic benchmarking of
    grackle.

    We summarize the subcommands down below. Be note that the current
    organization is experimental and is subject to change

      -> bench sets up an example problem and benchmarks it

      -> show-parameters is intended to display the configured parameters.

      -> show-initial-units is intended to display the initial units used
         during configuration of the grackle-solver. In the future, we may add
         support for customizing these units and using a different set of units
         during Grackle operations

OPTIONS:
    Some Flags are grouped into categories

    -S...
        Prefix for all "scenario arguments". The "scenario" is the algorithm
        for initializing all grackle fields that are used in a "test-problem."

        At the moment, this might include
          * -Sgrid.ax<I>=<quantity>;<seq>
        Where the quantities in brackets get replaced by various things. For
        example:
          * <I> is replaced with 0, 1, or 2
          * <quantity> is replaced with temperature, density, metallicity
          * <seq> is replaced by arguments that specify a sequence of 1 or more
            values. In practice, this is either a lone floating-point value OR
            geomspace;<start>;<end>;<num>

    -O...
        Prefix for all "operation arguments". Essentially, this denotes the
        operation that is executed in a given "test-problem." Currently, you
        can exclusively specify
          * -Ocalc-cooling-time  Computes the cooling time
          * -Ocalc-dust-temperature  Computes the dust-temperature
          * -Ocalc-pressure  Computes the pressure
          * -Ocalc-temperature  Computes the temperature
          * -Osolve-chemistry-dt=<dt-value>  Solves chemistry for the specified
            timestep (given in code-units).

    -h; --help
        Prints this message

    --version
        Prints the version of Grackle

    --precision
        Prints the floating-point precision selected while compiling Grackle

  PARAMETER OPTIONS

    --par-start
        Every "word" after this option (until the --par-stop option) is parsed
        as a key-value pair that is used as a grackle parameter.

    --par-stop
        The sentinel value that must follow --par-start
)HELP";

/// prints the help message
///
/// @note
/// when somebody passes an -h flag, a help message is the expected result, so
/// we print to stdout (rather than stderr)
void print_help(const char* bin_name) {

  // adjust bin_name so that it only specifies the base name of the program
  std::string_view bin_name_view(bin_name);
  std::size_t pos = bin_name_view.rfind('/');
  if ((pos != std::string_view::npos) && ((pos + 1) < bin_name_view.size())) {
    bin_name = bin_name + (pos + 1);
  }

  print_name_section_(stdout, bin_name);
  print_synopsis_section_(stdout, bin_name);
  std::printf("%s", help_message_);
}

int main(int argc, char* argv[]) {
  CliParser parser(argc, argv);

  // parse the first argument (i.e. the subcommand)
  const char* arg = parser.next();
  if (arg == nullptr) {
    std::fprintf(stderr, "%s is missing a subcommand", parser.bin_name());
    std::exit(1);
  }
  // handle the case where arg is a -h or --help
  try_parse_help(arg, parser.bin_name());

  std::string_view subcommand(arg);

  if (subcommand == "bench") {
    cmd::bench::run(parser);
  } else if (subcommand == "show-parameters") {
    cmd::show::run(parser, cmd::show::Kind::parameters);
  } else if (subcommand == "show-initial-units") {
    cmd::show::run(parser, cmd::show::Kind::initial_units);
  } else if (subcommand == "--precision") {
    cmd::show::run(parser, cmd::show::Kind::precision);
  } else if (subcommand == "--version") {
    cmd::show::run(parser, cmd::show::Kind::version);
  } else {
    std::fprintf(stderr, "unknown subcommand: %s\n", arg);
    std::exit(1);
  }

  GRCLI_ERROR("Something is wrong, this shouldn't be executed!");
}



