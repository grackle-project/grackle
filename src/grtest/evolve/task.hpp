/// @file task.hpp
/// @brief implement the tasks used within the integrator

#ifndef GRTEST_EVOLVE_TASK_HPP
#define GRTEST_EVOLVE_TASK_HPP

#include <grackle.h>
#include <string>
#include <utility>  // std::move, std::pair
#include <vector>

#include "evolve.hpp"

namespace grtest::impl {

// all Tasks can be wrapped by:
using TaskFn = std::function<int(const grtest::IntegrationState&, double)>;

/// helper function
bool chain_find_(gr_float*& out,
                 const std::string& key,
                 const std::map<std::string, gr_float*>& primary_map,
                 const std::map<std::string, gr_float*>& secondary_map)
{
  auto search_primary = primary_map.find(key);
  if (search_primary != primary_map.end()) {
    out = search_primary->second;
    return true;
  }
  auto search_secondary = secondary_map.find(key);
  if (search_secondary != secondary_map.end()) {
    out = search_secondary->second;
    return true;
  }
  return false;
}

/// A callable that is used to record data every cycle
class RecordTask {
  using RegisterBufferPair = std::pair<gr_float*, gr_float*>;

  int register_len;
  int batches_before_wrap;
  std::vector<RegisterBufferPair> pairs;
  gr_float* dt_buffer = nullptr;

  RecordTask() = default;  // this makes the default constructor private
public:

  /// record the current values to the registerred buffers
  int operator()(const grtest::IntegrationState& state, double dt) const {
    int batch = state.cycle % batches_before_wrap;
    // dt_buffer is special because it always holds a single element
    if (dt_buffer != nullptr) {
      dt_buffer[batch] = static_cast<gr_float>(state.t + dt);
    }

    int buffer_start = batch * register_len;
    for (const RegisterBufferPair& pair: pairs) {
      const gr_float* register_ptr = pair.first;
      gr_float* buffer_ptr = pair.second;
      for (int i = 0; i < register_len; i++) {
        buffer_ptr[buffer_start + i] = register_ptr[i];
      }
    }
    return GR_SUCCESS;
  }

  /// factory method that creates a new RecordTask
  static std::pair<std::string, RecordTask> create(
    int register_len, int buffer_len,
    std::map<std::string, gr_float*> chem_registers,
    std::map<std::string, gr_float*> extra_registers,
    std::map<std::string, gr_float*> buffer_map
  ) {
    RecordTask task;
    if ((register_len <= 0) || (buffer_len < register_len)) {
      return {"register_len must be positive & can't exceed buffer_len", task};
    }
    task.register_len = register_len;
    task.batches_before_wrap = buffer_len / register_len;

    // extract buffers from buffer_map
    for (const auto& [key, buffer_ptr]: buffer_map) {
      gr_float* tmp_reg = nullptr;

      if (buffer_ptr == nullptr) {
        return {"nullptr buffer for: " + key, task};

      } else if (chain_find_(tmp_reg, key, chem_registers, extra_registers)) {
        if (tmp_reg == nullptr) { return {"register is nullptr: " + key, task}; }
        task.pairs.push_back({tmp_reg, buffer_ptr});

      } else if (key == "time") {
        task.dt_buffer = buffer_ptr;

      } else {
        return {"unknown register: " + key, task};
      }
    }
    return {"", task};
  }

};

/// this is the "task" for solving chemistry
struct SolveChemistryTask {
  grtest::GrackleTypePack pack;

  int operator()(const grtest::IntegrationState& state, double dt) const {
    return local_solve_chemistry(
      pack.my_chem, pack.my_rates, pack.my_units, pack.my_fields, dt
    );
  }
};

} // namespace grtest::impl

#endif /* GRTEST_EVOLVE_TASK_HPP */
