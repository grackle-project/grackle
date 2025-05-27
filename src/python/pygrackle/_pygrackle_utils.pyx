# cimports:
from libcpp cimport bool as CppBool
from libcpp.map cimport map as CppMap
from libcpp.string cimport string as CppString
from libcpp.utility cimport pair
from libcpp.vector cimport vector

from .grackle_defs cimport (
    gr_float,
    c_chemistry_data,
    c_chemistry_data_storage,
    c_code_units,
    c_field_data
)
from .grackle_wrapper cimport GrackleTypePack_ExtType

# runtime imports:
from collections import defaultdict
import numpy as np
from .fluid_container import FluidContainer, _calculated_fields, _fc_calculated_fields
from .grackle_wrapper import _get_grackle_type_pack
from .utilities.physical_constants import sec_per_year

cdef extern from "grtest/evolve/evolve.hpp" namespace "grtest":
    cdef cppclass IntegrationState:
        IntegrationState()
        unsigned int cycle
        double t

    cdef cppclass GrackleTypePack:
        GrackleTypePack()
        c_chemistry_data* my_chem
        c_chemistry_data_storage* my_rates
        c_code_units* my_units
        c_field_data* my_fields

    cdef cppclass IntegratorBuilder:
        IntegratorBuilder()
        IntegratorBuilder& safety_factor(double)
        IntegratorBuilder& max_time(double)
        IntegratorBuilder& min_stopping_temperature_cgs(double)
        IntegratorBuilder& max_stopping_density_cgs(double)
        IntegratorBuilder& config_freefall(double gravitational_constant_cgs,
                                           CppBool pressure_free,
                                           vector[CppString] density_fields,
                                           CppBool strict_density_field_match)

cdef extern from *:
    """
    #include <map>
    #include <string>
    #include <utility>

    #include "grtest/evolve/evolve.hpp"
    typedef int (*callbackFn)(void* callback_ctx, unsigned int cur_cycle);

    std::string create_and_exec_integrator(
      grtest::IntegrationState& state,
      const grtest::IntegratorBuilder& builder,
      const grtest::GrackleTypePack& pack,
      std::map<std::string, gr_float*> chem_register_map,
      std::map<std::string, gr_float*> output_bufs,
      int buffer_len,
      unsigned int callback_every_n,
      callbackFn callback,
      void* callback_ctx
    ) {
      if ((callback == nullptr) || (callback_ctx == nullptr)) {
        return "callback and callback_ctx must be non-NULL";
      }

      std::pair<std::string, grtest::IntegratorFn> rslt_pair = builder.build(
        pack, chem_register_map, output_bufs, buffer_len
      );

      if (rslt_pair.first.size() != 0) {
        return "Error while building integrator:" + rslt_pair.first;
      }
      grtest::IntegratorFn& integrator_fn = rslt_pair.second;

      while (true) {
        grtest::EvolveRslt rslt = integrator_fn(state, callback_every_n);
        if (rslt.problem_phase != nullptr) {
          return rslt.to_err_msg();
        } else if ((callback(callback_ctx, state.cycle) != 1) ||
                   rslt.triggered_stopping) {
          return "";
        }
      }
    }
    """
    ctypedef int (*callbackFn)(void*, unsigned int)

    cdef CppString create_and_exec_integrator(
        const IntegrationState& state,
        const IntegratorBuilder& builder,
        const GrackleTypePack& pack,
        CppMap[CppString, gr_float*] chem_register_map,
        CppMap[CppString, gr_float*] output_bufs,
        int buffer_len,
        unsigned int callback_every_n,
        callbackFn callback,
        void* callback_ctx
    )

cdef int invoke_callable(void* python_callable,
                         unsigned int current_cycle) noexcept:
    cdef object rslt = (<object>python_callable)(current_cycle)
    return rslt


_DERIVED_FIELDS = _calculated_fields + _fc_calculated_fields

class DataFlushCallback:
    """
    This is a callback for data recording with an integrator.

    Every cycle, the integrator will directly record values to each of the
    buffers, held in the `tmp_buffer` attribute. Every `n_cylces_per_flush`,
    this callback is invoked to copy the data out of `tmp_buffer` into the
    `dict` attribute (i.e. we flush the data).

    Attributes
    ----------
    n_cycles_per_flush: int
        Specifies the cadence for flushing the buffers
        The integrator fills the buffers, held by the tmp_buffers
        attribute ). Then every
        `n_cycles_per_flush`, this callback is invoked and moves the data out
        of self.tmp_buffers into self.data.
    compute_derived_quan: bool
        When True, we computed derived chemistry quantities during each flush.
    """

    def __init__(
        self,
        fc,
        n_cycles_per_flush=1,
        compute_derived_quan=False,
        extra_buffers = None
    ):
        """
        Initializes a new instance

        Parameters
        ----------
        fc: FluidContainer
            This holds the information used during integration
        n_cycles_per_flush: int
            The cadence for flushing data
        compute_derived_quan: bool
            When True, we computed derived chemistry quantities during flushes.
        extra_buffers
            A sequence of non-chemistry quantity buffers that should be copied.
            This should not contain "time" or any derived quantites
        """
        self.data = defaultdict(list)
        if n_cycles_per_flush <= 0:
            raise ValueError("n_cycles_per_flush must be positive")
        self.n_cycles_per_flush = n_cycles_per_flush

        shape = (fc["density"].size * n_cycles_per_flush,)
        dtype = fc["density"].dtype
        if compute_derived_quan:
            self._secondary_fc = FluidContainer(
                fc.chemistry_data, shape[0], dtype=dtype
            )
            self.tmp_buffers = {k: self._secondary_fc[k] for k in fc.input_fields}
        else:
            self._secondary_fc = None
            self.tmp_buffers = {
                k: np.empty(shape=shape, dtype=dtype) for k in fc.input_fields
            }

        self.tmp_buffers["time"] = np.empty(shape=shape, dtype=dtype)
        if extra_buffers is not None:
            for key in extra_buffers:
                if key in self.tmp_buffers or key in _DERIVED_FIELDS:
                    raise ValueError(f"extra_buffers can't hold {key}")
                self.tmp_buffers[key] = np.empty(shape=shape, dtype=dtype)

    @property
    def compute_derived_quan(self):
        """When True, we compute derived quantities during each flush."""
        return self._secondary_fc is not None

    @property
    def buffer_len(self):
        return self.tmp_buffers["density"].size

    def __call__(self, cycle):
        """
        This is the "flush". We copy the values out of self.tmp_buffers into
        self.data
        """
        copy_n_cycles = self.n_cycles_per_flush - (cycle % self.n_cycles_per_flush)
        stop = copy_n_cycles * self.buffer_len
        for key, buf in self.tmp_buffers.items():
            self.data[key].extend(buf[:stop])

        if self.compute_derived_quan:
            is_first_flush = len(self.data["temperature"]) == 0
            if stop != self.n_cycles_per_flush and is_first_flush:
                for key in self._secondary_fc.input_fields:
                    self.tmp_buffers[key][stop:] = self.tmp_buffers[key][0]
            for key in _DERIVED_FIELDS:
                func = getattr(self._secondary_fc, f"calculate_{key}")
                if func is None:
                    raise RuntimeError(f"No function for calculating {key}.")
                func()
                buf = self._secondary_fc[key]
                self.data[key].extend(self._secondary_fc[key][:stop])

class CallbackLogAndFlush:
    def __init__(self, log_prefix, flusher, chemistry_data):
        self.log_prefix = log_prefix
        self.flusher = flusher
        self.time_units = chemistry_data.time_units
        self.density_units = chemistry_data.density_units
        self.exception = None

    def __call__(self, cycle):
        try:
            self.flusher(cycle)
            t = self.flusher.data["time"][-1] * self.time_units / sec_per_year
            rho = self.flusher.data["density"][-1] * self.density_units
            T = self.flusher.data["temperature"][-1]
            print(f"{self.log_prefix}t: {t:e} yr, rho: {rho:e} g/cm^3, T: {T:e} K.")
            return 1
        except BaseException as e:
            exception = e
            return 0


cdef CppMap[CppString, gr_float*] construct_ptr_map(object d, object keys=None):
    cdef CppMap[CppString, gr_float*] out
    cdef CppString tmp_string
    cdef gr_float[::1] tmp_view
    cdef object it = d.keys() if keys is None else keys
    for key in it:
        tmp_string = key.encode("ascii")
        tmp_view = d[key]
        out[tmp_string] = <gr_float *> &tmp_view[0]
    return out


cdef object invoke_integrator(
    object fc,
    const IntegratorBuilder& builder,
    int n_cycles_per_flush=1,
    object compute_derived_quan=False,
    object extra_buffers=None,
    object log_prefix=""
) except+:
    """
    Parameters
    ----------
    fc: FluidContainer
        The fluid container
    """
    cdef IntegrationState state
    state.cycle = 0
    state.t = 0.0

    # now create and initialize GrackleTypePack
    cdef GrackleTypePack_ExtType pack_ext_type = _get_grackle_type_pack(fc)
    cdef GrackleTypePack pack
    pack.my_chem = pack_ext_type.my_chem
    pack.my_rates = pack_ext_type.my_rates
    pack.my_units = pack_ext_type.my_units
    pack.my_fields = &pack_ext_type.my_fields

    # in the future, chem_register_map won't be necessary
    cdef CppMap[CppString, gr_float*] chem_register_map = construct_ptr_map(
        fc, keys=fc.input_fields
    )

    cdef object flusher = DataFlushCallback(
        fc,
        n_cycles_per_flush=n_cycles_per_flush,
        compute_derived_quan=compute_derived_quan,
        extra_buffers=extra_buffers
    )

    cdef object cb = CallbackLogAndFlush(log_prefix, flusher, fc.chemistry_data)

    cdef CppString result_str = create_and_exec_integrator(
        state=state,
        builder=builder,
        pack=pack,
        chem_register_map=chem_register_map,
        output_bufs=construct_ptr_map(flusher.tmp_buffers),
        buffer_len=flusher.buffer_len,
        callback_every_n=n_cycles_per_flush,
        callback=invoke_callable,
        callback_ctx=<void*>cb
    )

    if (result_str.size() != 0):
        raise RuntimeError(result_str.decode("ascii"))
    elif cb.exception is not None:
        raise RuntimeError(
            "The callback used during integration encountered an exception"
        ) from cb.exception
    else:
        return flusher.data

def evolve_constant_density(
    fc, final_temperature = None, final_time = None, safety_factor = 0.01,
):
    """
    Evolves the fluid_container while holding density constant
    """
    cdef IntegratorBuilder builder

    if final_temperature is not None:
        if final_temperature <= 0:
            raise ValueError("final_temperature must exceed 0")
        builder.min_stopping_temperature_cgs(final_temperature)

    if final_time is not None:
        if final_time <= 0:
            raise ValueError("final_time must exceed 0")
        builder.max_time(final_time)

    builder.safety_factor(safety_factor)

    prefix = "Evolve constant density - "
    data = invoke_integrator(
        fc, builder, 1, compute_derived_quan=True, extra_buffers=[], log_prefix=prefix
    )
    for field in data:
        data[field] = np.squeeze(np.array(data[field]))
    return fc.finalize_data(data=data)
