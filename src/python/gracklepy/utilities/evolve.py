########################################################################
#
# Functions for evolving a fluid container using Grackle
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from collections import defaultdict
from collections.abc import Callable, Mapping
import math

import numpy as np

from gracklepy.fluid_container import FluidContainer, _indirectly_calculated_fields
from gracklepy import chemistry_data

from .physical_constants import \
    gravitational_constant_cgs, \
    sec_per_year

_MAX_64BIT_SIGNED_INT = 2**63 - 1



class StateRecorder:
    """
    Aggregate many batches of field-states while evolving a FluidContainer

    The premise is that we can subsequently compute the derived properties
    all at once (for sufficiently complex systems, computing all derived
    properties after a call to fc.solve_chemistry can be more expensive than
    the call to fc.solve_chemistry)
    """

    # number of recordings before object must be reset
    n_recordings_per_reset_: int
    # current number of recordings
    n_recordings_: int

    # number of values copied out of each field when a `FluidContainer` is recorded
    # (if we plan to record a FluidContainer, `fc`, this is `fc.n_vals`)
    field_register_len_: int
    # buffers for recorded field data
    field_data_: dict[str, np.ndarray]
    # tracks recorded times
    time_data_: np.ndarray

    def __init__(
        self,
        n_recordings_per_reset: int,
        field_register_len: int,
        field_names: list[str]
    ):
        # field_names specifies the field names that will be recorded. This should
        # generally just be the input field names (the recorder will NOT compute
        # derived fields)
        assert n_recordings_per_reset > 0
        self.n_recordings_per_reset_ = n_recordings_per_reset
        assert field_register_len > 0
        self.field_register_len_ = field_register_len
        self.n_recordings_ = 0

        field_shape = (n_recordings_per_reset * field_register_len,)
        self.field_data_ = {name: np.empty(shape=field_shape, dtype=np.float64)
                            for name in field_names}
        self.time_data_ = np.empty(shape=field_shape, dtype=np.float64)

    def is_full(self) -> bool:
        return self.n_recordings_ == self.n_recordings_per_reset_

    def is_empty(self) -> bool:
        return self.n_recordings_ == 0

    def record(self, time: float, fc: FluidContainer) -> bool:
        """Record the current state

        This does NOT compute derived fields.
        """
        assert fc.n_vals == self.field_register_len_
        if self.is_full():
            return False

        dst_start = self.n_recordings_ * self.field_register_len_
        dst_stop = dst_start + self.field_register_len_
        for key, dst_buf in self.field_data_.items():
            dst_buf[dst_start:dst_stop] = fc[key]

        self.time_data_[self.n_recordings_] = time

        self.n_recordings_+=1
        return True

    def reset(self):
        self.n_recordings_=0

    def recorded_times_and_fields(self) -> tuple[np.ndarray, dict[str,np.ndarray]]:
        times = self.time_data_[:self.n_recordings_]
        slc = slice(0, self.n_recordings_*self.field_register_len_)
        fields = {k:buf[slc] for k,buf in self.field_data_.items()}
        return (times, fields)

class _Evolver:
    """
    Solves evolution of isolated zones

    The system can be integrated until time reaches `t_max` or an
    exit-condition optionally encoded in a callback passed to `calc_dt`
    has been met.

    This has only been tested with evolving a single zone, but there's no
    reason it can't support multiple zones

    All callback functions take the same set of function arguments:
    `(time, dt, fc)`.

    Parameters
    ----------
    t0
        Initial time.
    fc0
        Encodes the initial state of the system. The state will be copied out
        of the provided object (i.e. updates in the state will **NOT** be
        reflected in the object passed as this argument)
    t_max
        The max integration time. Passing ``math.inf`` is means that the max
        integration time is unbounded.
    calc_dt
        Specifies a fixed timestep or a callable for computing that computes
        a timestep from ``fc`` (it should *NOT* mutate fc). If this is a
        callable, the ``dt`` argument should be ignored. When a callable
        returns a negative value that's considered an indication that
        evolution should stop.
    prechem_update
        Optional callback function called before `solve_chemistry` during each
        timestep. It is expected to mutate `fc`
    recorder
        A StateRecorder object that will be filled

    Notes
    -----
    A lot of implementation choices for this class were made by a desire
    to easily convert this core logic into C++.
    """

    def __init__(
        self,
        t0: float,
        fc0: FluidContainer,
        t_max: float,
        calc_dt: Callable[[float, float, FluidContainer], float] | float,
        *,
        prechem_update: Callable[[float, float, FluidContainer], None] | None = None,
        recorder: StateRecorder | None = None,
    ):
        self._t = float(t0)
        self._t_max = float(t_max)
        if (self._t > self._t_max):
            raise ValueError("t0 cannot exceed t_max")
        self._fc = FluidContainer(fc0.chemistry_data, n_vals=fc0.n_vals)
        for field in self._fc.input_fields:
            self._fc[field][...] = fc0[field][...]

        if callable(calc_dt):
            self._calc_dt = calc_dt
        elif float(calc_dt) <= 0:
            raise ValueError("when calc_dt isn't callable, it must exceed 0")
        else:
            self._calc_dt = _DummyCallback(ret_val=float(calc_dt))

        self._prechem_update = (
            _DummyCallback() if prechem_update is None else prechem_update
        )
        self._recorder = recorder

    def t(self) -> float:
        """Returns the current time"""
        return self._t

    def record_state(self):
        """Invoke observer callback function on the current state
        """
        # This is primarily intended to support observing the first timestep
        if self._recorder is not None:
            self._recorder.record(self.t(), self._fc)

    def log_current_prop(self, msg_prefix):
        # frankly, I don't love that this is attached as a method, but we can revisit
        # this choice in the future

        self._fc.calculate_temperature()
        my_chemistry = self._fc.chemistry_data
        t = self._t * my_chemistry.time_units / sec_per_year
        rho = float(self._fc["density"][0] * my_chemistry.density_units)
        T = float(self._fc["temperature"][0])
        print(f"{msg_prefix} - t: {t:e} yr, rho: {rho:e} g/cm^3, T: {T:e} K")


    def copy_state(
        self, dst: Mapping[str, np.ndarray] | None = None
    ) -> Mapping[str, np.ndarray]:
        """Return a copy of the current state
        """
        f_list = self._fc.input_fields # <- fields to copy
        if dst is None:
            return {f: np.copy(self._fc[f]) for f in f_list}
        else:
            shape = (self._fc.n_vals,)
            for f in f_list:
                if f not in dst:
                    raise ValueError(f"dst doesn't have a key {f!r}")
                elif getattr(dst[f], "shape") != shape:
                    raise ValueError(f"dst[{f!r}] isn't an array with shape: {shape}")
            for f in f_list:
                dst[f][...] = self._fc[f][...]
        return dst

    def step(self, n_steps: int | None = None) -> int:
        """Advance evolution for n_steps or until stopping condition is met

        Returns the number of elapsed steps
        """
        _up_bound = _MAX_64BIT_SIGNED_INT  # <- for future c++ compat
        n_steps = _up_bound if n_steps is None else n_steps
        assert int(n_steps) == n_steps
        assert 0 < n_steps <= _up_bound

        has_recorder = self._recorder is not None

        for i in range(n_steps):
            t_at_tmax = self._t >= self._t_max
            dt = -1 if t_at_tmax else self._calc_dt(self._t, math.nan, self._fc)
            if t_at_tmax or (dt < 0):
                return i # <- intentionally don't add 1 to this

            self._prechem_update(self._t, dt, self._fc)  # <-pre-chemistry update
            self._fc.solve_chemistry(dt)  # <- perform the chemistry update

            # does it make more sense to call observer after updating timestep?
            if has_recorder:
                self._recorder.record(self._t, self._fc)
            self._t += dt
        return i + 1  # add one since the last step was complete

class _DummyCallback:
    def __init__(self, ret_val = None):
        self.ret_val = ret_val
    def __call__(self, *args, **kwargs):
        return self.ret_val

def evolve_freefall(fc, final_density, safety_factor=0.01,
                    include_pressure=True):
    my_chemistry = fc.chemistry_data

    n_recordings_per_reset = 100
    state_recorder = StateRecorder(
        n_recordings_per_reset=n_recordings_per_reset,
        field_register_len=fc.n_vals,
        field_names=fc.input_fields
    )

    data = defaultdict(list)


    ff_updater = FFDensityAndEnergyUpdate(
        n_vals=fc.n_vals, include_pressure=include_pressure, my_chemistry=fc.chemistry_data
    )

    force_factor_list = []
    dt_list = []

    def prechem_update(time, dt, fc):
        dt_list.append(dt)
        # perform analytic update from our free-fall model. Then, record force_factor
        force_factor = ff_updater(time, dt, fc)
        force_factor_list.append(force_factor)


    gravitational_constant = ff_updater.gravitational_constant

    def calc_dt(time, dt_IGNORE, fc):
        if fc["density"][0] * my_chemistry.density_units >= final_density:
            return -1  # <- indicates that evolution should stop now
        return safety_factor * \
          np.power(((3. * np.pi) /
                    (32. * gravitational_constant *
                     fc["density"][0])), 0.5)

    evolver = _Evolver(
        t0=0.0,
        fc0=fc,
        t_max=math.inf,
        calc_dt=calc_dt,
        prechem_update=prechem_update,
        recorder=state_recorder
    )
    evolver.record_state() # <- observe t=0 state
    add_to_data(fc.chemistry_data, state_recorder, data)
    state_recorder.reset() # <- clear out contents

    while evolver.step(n_steps=n_recordings_per_reset) > 0:
        add_to_data(fc.chemistry_data, state_recorder, data)
        state_recorder.reset() # <- clear out contents
        evolver.log_current_prop(msg_prefix="Evolve Freefall")

    evolver.copy_state(dst=fc)

    data["time"][1:] += np.array(dt_list)
    for field in data:
        data[field] = np.squeeze(np.array(data[field]))[:-1]
    data["force_factor"] = np.squeeze(np.array(force_factor_list))
    return fc.finalize_data(data=data)

class FFDensityAndEnergyUpdate:

    def __init__(self, n_vals: int, include_pressure: bool, my_chemistry: chemistry_data):
        self.last_n_recorded = 0
        self.last3_density = np.full(shape=(3, n_vals), fill_value=np.nan, dtype=np.float64)
        self.last3_pressure = np.full(shape=(3, n_vals), fill_value=np.nan, dtype=np.float64)
        self.include_pressure = include_pressure

        # Set units of gravitational constant
        self.gravitational_constant = (
            4.0 * np.pi * gravitational_constant_cgs *
            my_chemistry.density_units * my_chemistry.time_units**2)

        # some constants for the analytical free-fall solution
        self.freefall_time_constant = np.power(((32. * self.gravitational_constant) /
                                                (3. * np.pi)), 0.5)

    def __call__(
        self, current_time: float, dt: float, fc: FluidContainer
    ) -> np.ndarray:
        # returns the force-factor
        my_chemistry = fc.chemistry_data

        if self.include_pressure:
            self.last3_density[-3:-1, :] = self.last3_density[-2:, :]
            self.last3_pressure[-3:-1, :] = self.last3_pressure[-2:, :]
            fc.calculate_pressure()
            self.last3_density[-1, :] = fc["density"]
            self.last3_pressure[-1, :] = fc["pressure"]

            if self.last_n_recorded < 3:
                self.last_n_recorded += 1

            force_factor = calculate_collapse_factor(
                self.last3_pressure[-self.last_n_recorded:],
                self.last3_density[-self.last_n_recorded:])
        else:
            force_factor = np.zeros(shape=(fc.n_vals,), dtype=np.float64)
        
        # compute the new density using the modified
        # free-fall collapse as per Omukai et al. (2005)

        # calculate new density from altered free-fall solution
        new_density = np.power((np.power(fc["density"][0], -0.5) -
                                (0.5 * self.freefall_time_constant * dt *
                                 np.power((1 - force_factor), 0.5))), -2.)

        # use this to multiply by elemental densities if you are tracking those
        density_ratio = new_density / fc["density"][0]

        # update densities
        for field in fc.density_fields:
            fc[field] *= density_ratio

        # now update energy for adiabatic heating from collapse
        fc["internal_energy"][0] += (my_chemistry.Gamma - 1.) * \
          fc["internal_energy"][0] * self.freefall_time_constant * \
          np.power(fc["density"][0], 0.5) * dt

        return force_factor


def calculate_collapse_factor(pressure, density):
    # Calculate the effective adiabatic index, dlog(p)/dlog(rho).

    if len(pressure) < 3:
        return np.array([0.])

    # compute dlog(p) / dlog(rho) using last two timesteps
    gamma_eff = np.log10(pressure[-1] / pressure[-2]) / \
        np.log10(density[-1] / density[-2])

    # compute a higher order derivative if more than two points available
    if len(pressure) > 2:
        gamma_eff += 0.5 * ((np.log10(pressure[-2] / pressure[-3]) /
                             np.log10(density[-2] / density[-3])) - gamma_eff)

    gamma_eff = np.clip(gamma_eff, a_min=0, a_max=4/3)

    # Equation 9 of Omukai et al. (2005)
    if gamma_eff < 0.83:
        force_factor = gamma_eff * 0
    elif gamma_eff < 1.0:
        force_factor = 0.6 + 2.5 * (gamma_eff - 1) - \
            6.0 * np.power((gamma_eff - 1.0), 2.)
    else:
        force_factor = 1.0 + 0.2 * (gamma_eff - (4./3.)) - \
            2.9 * np.power((gamma_eff - (4./3.)), 2.)
    force_factor = np.clip(force_factor, a_min=0, a_max=0.95)
    return force_factor

def evolve_constant_density(fc, final_temperature=None,
                            final_time=None, safety_factor=0.01):
    my_chemistry = fc.chemistry_data

    if final_temperature is None and final_time is None:
        raise RuntimeError("Must specify either final_temperature " +
                           "or final_time.")

    n_recordings_per_reset = 100
    data = defaultdict(list)

    state_recorder = StateRecorder(
        n_recordings_per_reset=n_recordings_per_reset,
        field_register_len=fc.n_vals,
        field_names=fc.input_fields
    )

    # compute the fixed timestep
    fc.calculate_cooling_time()
    fixed_dt = safety_factor * np.abs(fc["cooling_time"][0])

    def calc_dt(time, dt_IGNORE, fc):
        if final_temperature is not None:
            fc.calculate_temperature()
            if fc["temperature"][0] <= final_temperature:
                return -1  # <- indicates that evolution should stop now
        return fixed_dt


    evolver = _Evolver(
        t0=0.0,
        fc0=fc,
        t_max=math.inf if final_time is None else final_time,
        calc_dt=calc_dt,
        recorder=state_recorder,
    )
    while evolver.step(n_steps=n_recordings_per_reset) > 0:
        add_to_data(fc.chemistry_data, state_recorder, data)
        state_recorder.reset() # <- clear out contents
        evolver.log_current_prop(msg_prefix="Evolve constant density")

    evolver.copy_state(dst=fc)

    for field in data:
        data[field] = np.squeeze(np.array(data[field]))
    return fc.finalize_data(data=data)


def add_to_data(chemistry_data, state_recorder, data):
    """
    Add current state_recorder values to the data structure.
    """

    if state_recorder.is_empty():
        return

    times, field_data = state_recorder.recorded_times_and_fields()
    n_vals = field_data["density"].size

    fc = FluidContainer(chemistry_data=chemistry_data, n_vals=n_vals)
    for field_name, buf in field_data.items():
        fc[field_name][:] = buf

    # for indirectly calculated fields, we get more than one field from a call.
    # keep a list of what we've called so we don't have to call them again.
    called = []

    for field in fc.all_fields:
        if field not in fc.input_fields:
            if field in _indirectly_calculated_fields:
                fname = _indirectly_calculated_fields[field]
            else:
                fname = f"calculate_{field}"

            func = getattr(fc, fname, None)
            if func is None:
                raise RuntimeError(f"No function for calculating {field}.")

            if fname not in called:
                func()
                called.append(fname)

        data[field].extend(fc[field])
    data["time"].extend(times)
