#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string_view>


#include "grid_problem.h"
#include "handle_gr_parameters.h"
#include "utils.h"


static void print_usage_(FILE* f, const char* bin_name) {
  std::fprintf(f,
      R"HELP(
USAGE:

    %s --version
    %s --precision
    %s bench -S... PARAMETER_OPTIONS
    %s show-parameters PARAMETER_OPTIONS
)HELP",
      bin_name, bin_name, bin_name, bin_name);
}

static const char* help_message_ = R"HELP(
DESCRIPTION:
    A simple command-line tool intended for debugging and basic benchmarking of
    grackle.


    The different modes are currently experimental:
      -> bench sets up an example problem and tries to run a benchmark
      -> show-parameters shows the parameters... (this may change a little in the future)

OPTIONS:
    Some Flags are grouped into categories

    -S
        Prefix for all "scenario arguments".

        At the moment, this might include
          * -Sgrid.ax<I>=<quantity>;<seq>
        Where the quantities in brackets get replaced by various things. For
        example:
          * <I> is replaced with 0, 1, or 2
          * <quantity> is replaced with temperature, density, metallicity
          * <seq> is replaced by arguments that specify a sequence of 1 or more
            values. In practice, this is either a lone floating-point value OR
            geomspace;<start>;<end>;<num>

        
    -h; --help
        Prints this message

    --precision
        Prints whether the program was compiled with single or double precision

    --version
        Prints the program's version information

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
void print_help_(const char* bin_name) {
  print_usage_(stdout, bin_name);
  std::printf("%s", help_message_);
}


struct Cli {
  std::optional<CliParamSpec> gr_param_spec;
  //std::optional<GridScenarioSpec> grid_spec_scenario;
};

Cli parse(int argc, char** argv) {
  GRCLI_REQUIRE(argc > 0,
                "Something is wrong, argc is always expected to be positive");
  const char* bin_name = argv[0];

  // this is a weird workaround to use pointers of the type `const char**`,
  // since it is illegal to directly cast `char**` to `const char**`
  const char* const * ptr = argv + 1;
  const char* const * stop = argv + argc;

  std::string_view subcommand{(argc < 2) ? "" : argv[1]};
  if (subcommand != "bench") {
    std::fprintf(stderr, "at the moment, we always expect `bench` to be the first arg\n");
    std::exit(1);
  }
  ptr++;

  Cli out;
  while (ptr < stop) {
    const std::string_view cur(*ptr);

    if ((cur == "-h") || (cur == "--help")) {
      print_help_(bin_name);
      std::exit(0);
    } else if (cur == "--") {
      std::fprintf(stderr,
                   "We don't respect the \"--\" cli argument since we don't "
                   "currently support positional arguments.");
      std::exit(1);
    } else if ((cur == "--version") || (cur == "--precision" )) {
      std::fprintf(stderr, "We don't support the \"%s\" flag yet!\n", *ptr);
      std::exit(1);
    } else if (check_start_cli_paramspec(*ptr)){
      GRCLI_REQUIRE(!out.gr_param_spec.has_value(),
                    "Grackle parameters cli group was provided more than once");
      out.gr_param_spec = parse_cli_paramspec(ptr, stop);

      if (out.gr_param_spec->end == stop) {
        break; // in this case, nothing comes after the parameter
      } else {
        ptr = (1 + out.gr_param_spec->end);
      }

    } else {
      std::fprintf(stderr, "unrecognized argument: %s\n", *ptr);
      std::exit(1);
    }
  }

  return out;

}


int main(int argc, char* argv[]) {

  //test();
  Cli cli_args = parse(argc, argv);

  if (cli_args.gr_param_spec.has_value()) {
    ChemistryData my_chem;
    init_gr_params(*cli_args.gr_param_spec, my_chem);

    // in the future, we will support customization of units!
    double initial_redshift = 1.0;
    code_units my_units;
    my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
    my_units.density_units = 1.67e-24;
    my_units.length_units = 3.0857e18; // ~ 1 parsec
    my_units.time_units = 3.15e13; // ~ 1 Myr
    my_units.a_units = 1.0; // units for the expansion factor
    // Set expansion factor to 1 for non-cosmological simulation.
    my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;

    
  }

  /*
  seq_variant seq = get_seq("geomspace;2.43243e-32;+10.231e34;3");

  double lo = std::visit([](auto&& arg) -> double { return arg.lo; }, seq);
  printf("start = %g\n", lo);
  */
  return 0;

}
