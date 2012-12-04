import numpy as np

_fluid_names = {}
_fluid_names[0] = ["HI", "HII", "HeI", "HeII", "HeIII", "de",
    "x-velocity", "y-velocity", "z-velocity",
    "temperature", "gamma", "energy", "density", "metal_density"]
_fluid_names[1] = _fluid_names[0] + \
    ["H2I", "H2II", "HM"]
_fluid_names[2] = _fluid_names[1] + \
    ["DI", "DII", "HDI"]

class FluidContainer(dict):
    def __init__(self, chemistry_data, n_vals, dtype="float64"):
        super(FluidContainer, self).__init__()
        self.dtype = dtype
        self.chemistry_data = chemistry_data
        self.n_vals = n_vals
        for fluid in _fluid_names[self.chemistry_data.primordial_chemistry]:
            self._setup_fluid(fluid)

    def _setup_fluid(self, fluid_name):
        self[fluid_name] = np.zeros(self.n_vals, self.dtype)
