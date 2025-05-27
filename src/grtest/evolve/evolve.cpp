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
  if (buffer_len < 0) {return {"buffer_len must be positive", {}}; }

  int register_len = grtest::impl::elements_per_field_ptr(*pack.my_fields);

  std::map<std::string, gr_float*> extra_registers;
  std::function<double()> calc_dt_fn;
  std::vector<grtest::impl::TaskFn> task_vec;

  if (ff_config_.has_value()) {
    // extra_registers["force"] = ...
    return {"NOT IMPLEMENTED YET", {}};
  } else {
    // get a function to compute the timestep from the minimum cooling time
    std::function<double()> tmp_fn = grtest::impl::make_tcool_dt_func(
      safety_factor_, pack
    );
    // we just want to work with the value based on
    // the initial timestep
    double initial_dt = tmp_fn();
    if (initial_dt < 0) { return {"something went wrong computing tcool", {}}; }

    // the function for calculating the timestep will always return initial_dt
    calc_dt_fn = [initial_dt]() { return initial_dt; };

    // now set up the tasks
    grtest::impl::SolveChemistryTask chem_task{pack};
    task_vec.push_back(chem_task);
  }

  bool any_bufs = output_bufs.size() > 0;
  if (any_bufs && buffer_len == 0) {
    return {"buffer_len can't be zero when buffers are provided", {}};
  } else if ((!any_bufs) && buffer_len == 0) {
    return {"buffers can't be provided when buffer_len is 0", {}};
  } else {
    std::pair<std::string, grtest::impl::RecordTask> pair =
      grtest::impl::RecordTask::create(
        register_len, buffer_len, chem_registers, extra_registers, output_bufs);
    if (pair.first.size() > 0) { return {pair.first, {}}; }

    task_vec.push_back(std::move(pair.second));
  }

  std::pair<std::string, grtest::impl::FieldStopTrigger> trigger_rslt =
    grtest::impl::FieldStopTrigger::create(stop_field_kind_,
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
  integrator_fn.fn_ = grtest::IntegratorFn::InnerFn(tmp);
  return {std::string(), integrator_fn};
}
