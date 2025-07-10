/// @file evolve.cpp
/// @brief implement the integrator logic

#include "evolve.hpp"
#include "misc.hpp"
#include "task.hpp"

#include <cmath>
#include <functional>
#include <utility> // std::pair, std::move

namespace { // stuff inside of an anonymous namespace is local to the
            // current translation unit


/// Implements the actual integrator
class IntegratorImpl {
  // attributes
  std::function<double()> calc_dt_fn_;
  std::vector<grtest::impl::TaskFn> task_vec_;
  double max_time_;
  grtest::impl::FieldStopTrigger field_stop_trigger_;

public:
  /// constructor
  IntegratorImpl(std::function<double()> calc_dt_fn,
                 std::vector<grtest::impl::TaskFn> task_vec,
                 double max_time,
                 grtest::impl::FieldStopTrigger field_stop_trigger)
    : calc_dt_fn_(calc_dt_fn),
      task_vec_(task_vec),
      max_time_(max_time),
      field_stop_trigger_(field_stop_trigger)
  { }

  /// Execute the integrator for up to stop_cycle
  grtest::EvolveRslt operator()(grtest::IntegrationState& state,
                                unsigned int stop_cycle = -1); // set to max val
};

grtest::EvolveRslt IntegratorImpl::operator()(grtest::IntegrationState& state,
                                              unsigned int stop_cycle) {
  const std::size_t num_tasks = task_vec_.size();

  // this is the "event loop" (the only way to exit is to return)
  for (unsigned int i = 0;; i++) {
    // check if we are ready to exit
    bool cycle_stop = i >= stop_cycle;
    std::optional<bool> exceed_thresh_rslt = field_stop_trigger_();
    if (!exceed_thresh_rslt) { return {"field_stop_trigger", 0, false}; }
    bool triggered_stopping = (state.t >= max_time_) || *exceed_thresh_rslt;
    if (cycle_stop || triggered_stopping) {
      return {nullptr, 0, triggered_stopping};
    }

    // calculate timestep
    double dt = calc_dt_fn_();
    if (dt < 0) { return {"calc_dt", 0, false}; }

    // apply each of the tasks
    for (std::size_t j = 0; j < num_tasks; j++) {
      int rv = task_vec_[j](state, dt);
      if (rv != GR_SUCCESS) {return {"task", static_cast<int>(j), false}; }
    }

    // update the state
    state.cycle++;
    state.t += dt;
  }
}

} // anonymous namespace


std::pair<std::string, grtest::IntegratorFn> grtest::IntegratorBuilder::build(
  const grtest::GrackleTypePack& pack,
  std::map<std::string, gr_float*> chem_registers,
  std::map<std::string, gr_float*> output_bufs,
  int buffer_len
) const
{
  namespace grimpl = ::grtest::impl;

  if (buffer_len < 0) {return {"buffer_len must be positive", {}}; }

  int register_len = grimpl::elements_per_field_ptr(*pack.my_fields);

  std::map<std::string, gr_float*> extra_registers;
  std::function<double()> calc_dt_fn;
  std::vector<grimpl::TaskFn> task_vec;

  if (ff_config_.has_value()) {
    // extra_registers["force"] = ...
    return {"NOT IMPLEMENTED YET", {}};
  } else {
    // get a function to compute the timestep from the minimum cooling time
    std::function<double()> tmp_fn = grimpl::make_tcool_dt_func(
      safety_factor_, pack
    );
    // we just want to work with the value based on
    // the initial timestep
    double initial_dt = tmp_fn();
    if (initial_dt < 0) { return {"something went wrong computing tcool", {}}; }

    // the function for calculating the timestep will always return initial_dt
    calc_dt_fn = [initial_dt]() { return initial_dt; };

    // now set up the tasks
    grimpl::SolveChemistryTask chem_task{pack};
    task_vec.push_back(chem_task);
  }

  // setup RecordTask
  bool any_bufs = output_bufs.size() > 0;
  if (any_bufs && buffer_len == 0) {
    return {"buffer_len can't be zero when buffers are provided", {}};

  } else if ((!any_bufs) && buffer_len == 0) {
    return {"buffers can't be provided when buffer_len is 0", {}};

  } else {
    bool end_of_timestep = prefer_trailing_record_;
    std::string err_msg; // only used if there is an error
    std::optional<grimpl::RecordTask> maybe_task = grimpl::RecordTask::create(
      register_len, buffer_len, end_of_timestep, chem_registers,
      extra_registers, output_bufs, err_msg
    );
    if (!maybe_task) {
      return {err_msg, {}};
    } else if (prefer_trailing_record_) {
      task_vec.push_back(std::move(*maybe_task));
    } else {
      task_vec.insert(task_vec.begin(), std::move(*maybe_task));
    }
  }

  // setup the machinery to detect whether a field-value (such as density or
  // temperature) triggered the stopping condition
  std::pair<std::string, grimpl::FieldStopTrigger> trigger_rslt =
    grimpl::FieldStopTrigger::create(stop_field_kind_,
                                           field_thresh_val_cgs_,
                                           pack);
  if (trigger_rslt.first.size() > 0) { return {trigger_rslt.first, {}}; }

  IntegratorImpl tmp(
    calc_dt_fn,
    std::move(task_vec),
    max_time_.value_or(std::numeric_limits<double>::max()),
    trigger_rslt.second
  );

  grtest::IntegratorFn integrator_fn;
  integrator_fn.fn = grtest::IntegratorFn::InnerFn(tmp);
  return {std::string(), integrator_fn};
}
