/// @file evolve.hpp
/// @brief implement the integrator logic
///
/// This is a quick and dirty implementation. Things will probably be refactored
/// if we start using the machinery from C++ for tests/benchmarks (that's OK
/// since this is a private API)
///
/// The idea here is that users employ a ``IntegratorBuilder`` to specify the
/// integrator's configuration. This is a lightweight mutable object. Invoking
/// the ``IntegratorBuilder::build`` method will provide the user with a
/// ``IntegratorFn`` object that does the heavy lifting.

#ifndef GRTEST_EVOLVE_HPP
#define GRTEST_EVOLVE_HPP

#include <grackle.h>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <utility> // std::pair
#include <vector>

namespace grtest {

/// Characterizes the current state of the integrator
struct IntegrationState {
  unsigned int cycle = 0;
  double t = 0.0;
};

/// Characterizes whether an error occurred
///
/// @note
/// This could obviously be improved to carry more information
struct EvolveRslt {
  const char* problem_phase;
  int task_id;
  bool triggered_stopping;

  std::string to_err_msg() const {
    if (problem_phase == nullptr) { return ""; }
    std::string phase(problem_phase);
    std::string msg("Error during the \"");
    msg.append(phase);
    msg.append("\" integration phase");
    if (phase == "task") {
      msg.append(", taskid = ");
      msg.append(std::to_string(task_id));
    }
    return msg;
  }
};


/// This is the integrator type that external is handed. In practice this is a
/// simple wrapper type.
///
/// The original plan was to use
/// ``typedef std::function<EvolveRslt(IntegrationState&, unsigned int)> ...``.
/// Unfortunately, that doesn't play well with Cython. Thus we introduce this
/// really lightweight wrapper type.
///
/// @par Background
/// The purpose of this type is to perform "type-erasure":
/// - different integrator-configurations can be implemented using distinct
///   types that have specialized implementations to improve performance (this
///   can be done with templates or more manually). To make it possible
///   to write code (without templates) that treats the various configurations
///   in a consistent manner, our factories return a "wrapper-type".
/// - The implementation is extremely trivial since the integrator acts like
///   a function.
struct IntegratorFn {
  typedef std::function<EvolveRslt(IntegrationState&, unsigned int)> InnerFn;
  InnerFn fn;

  /// evaluate the integrator
  EvolveRslt operator()(IntegrationState& state, unsigned int n_cycles) {
    return fn(state, n_cycles);
  }
};


namespace impl {

enum class FieldStopThresh {
  None, MinTemperature, MaxDensity
};

struct FreeFallConfig_{
  double gravitational_constant_cgs;
  bool presure_free;
  std::vector<std::string> density_fields;
  bool strict_density_field_match;
};

} // namespace impl


struct GrackleTypePack {
  chemistry_data* my_chem;
  chemistry_data_storage* my_rates;
  code_units* my_units;
  grackle_field_data* my_fields;
};


/// This type is used to construct an ``IntegratorFn`` instance via a builder
/// pattern.
///
/// The builder pattern is common in languages without named kwargs. The idea
/// is that the builder is mutable, relatively lightweight and very
/// configurable, and the result is immutable.
///
/// @note
/// Many of the builder methods return IntegratorBuilder& to make it possible
/// to chain methods together.
class IntegratorBuilder {
  /// indicates preference for recording at the start or end of cycle.
  bool prefer_trailing_record_ = true;

  /// the courant safety factor
  double safety_factor_ = 0.01;

  /// optionally specifies an output time
  std::optional<double> max_time_ = {};

  /// specifies a limiting field-value
  impl::FieldStopThresh stop_field_kind_ = impl::FieldStopThresh::None;
  double field_thresh_val_cgs_ = -1.0;

  /// when this is empty, we evolve with constant density
  std::optional<impl::FreeFallConfig_> ff_config_;

public:

  /// set preference for placement of the recording-task
  IntegratorBuilder& prefer_trailing_record(bool value)
  {
    prefer_trailing_record_ = value;
    return *this;
  }

  /// set the safety factor
  IntegratorBuilder& safety_factor(double value)
  { 
    safety_factor_ = value;
    return *this;
  }

  /// set the max time (in coude units)
  IntegratorBuilder& max_time(double value)
  { 
    max_time_ = {value};
    return *this;
  }

  IntegratorBuilder& min_stopping_temperature_cgs(double value)
  { 
    stop_field_kind_ = impl::FieldStopThresh::MinTemperature;
    field_thresh_val_cgs_ = value;
    return *this;
  }

  IntegratorBuilder& max_stopping_density_cgs(double value)
  { 
    stop_field_kind_ = impl::FieldStopThresh::MaxDensity;
    field_thresh_val_cgs_ = value;
    return *this;
  }

  IntegratorBuilder& config_freefall(double gravitational_constant_cgs,
                                     bool pressure_free,
                                     std::vector<std::string> density_fields,
                                     bool strict_density_field_match);

  /// Creates an Integrator Function that evolves the specified fields (assuming
  /// constant density)
  ///
  /// @param pack specifies the quantities that will be evolved
  /// @param chem_register_map maps names to pointers to the registers in
  ///   pack.my_fields. This is only to help with copying values into output
  ///   buffers during each cycle.
  /// @param output_bufs specifies output buffers to use for various data
  ///   registers.
  /// @param buffer_len The size of each output buffer
  ///
  /// On success, the std::string is empty. Otherwise the std::string provides
  /// an error message
  ///
  /// @note
  /// If/when PR #271 & #280 are merged we can get rid of the chem_register_map
  /// argument (since Grackle's public API will be able to map between strings
  /// and grackle_field pointers)
  std::pair<std::string, IntegratorFn> build(
    const GrackleTypePack& pack,
    std::map<std::string, gr_float*> chem_register_map,
    std::map<std::string, gr_float*> output_bufs,
    int buffer_len
  ) const;

};


} // namespace grtest

#endif /* GRTEST_EVOLVE_HPP */
