import numpy as np

from yt.fields.field_detector import \
    FieldDetector
from yt.frontends.enzo.data_structures import \
    EnzoDataset
from yt.utilities.physical_constants import \
    me, mp
from pygrackle import \
    FluidContainer, \
    chemistry_data

_parameter_map = {}
_parameter_map[EnzoDataset] = {
    "use_grackle": "use_grackle",
    "Gamma": "Gamma",
    "primordial_chemistry": "MultiSpecies",
    "metal_cooling": "MetalCooling",
    "h2_on_dust": "H2FormationOnDust",
    "cmb_temperature_floor": "CMBTemperatureFloor",
    "three_body_rate": "ThreeBodyRate",
    "cie_cooling": "CIECooling",
    "h2_optical_depth_approximation": "H2OpticalDepthApproximation",
    "photoelectric_heating": "PhotoelectricHeating",
    "photoelectric_heating_rate": "PhotoelectricHeatingRate",
    "NumberOfTemperatureBins": "NumberOfTemperatureBins",
    "CaseBRecombination": "CaseBRecombination",
    "TemperatureStart": "TemperatureStart",
    "TemperatureEnd": "TemperatureEnd",
    "NumberOfDustTemperatureBins": "NumberOfDustTemperatureBins",
    "DustTemperatureStart": "DustTemperatureStart",
    "DustTemperatureEnd": "DustTemperatureEnd",
    "HydrogenFractionByMass": "HydrogenFractionByMass",
    "DeuteriumToHydrogenRatio": "DeuteriumToHydrogenRatio",
    "SolarMetalFractionByMass": "SolarMetalFractionByMass",
    "UVbackground_redshift_on": "RadiationRedshiftOn",
    "UVbackground_redshift_off": "RadiationRedshiftOff",
    "UVbackground_redshift_fullon": "RadiationRedshiftFullOn",
    "UVbackground_redshift_drop": "RadiationRedshiftDropOff",
    "use_radiative_transfer": "RadiativeTransfer",
    "radiative_transfer_coupled_rate_solver": "RadiativeTransferCoupledRateSolver",
    "radiative_transfer_hydrogen_only": "RadiativeTransferHydrogenOnly",
    "with_radiative_cooling": "with_radiative_cooling",
    "use_volumetric_heating_rate": "use_volumetric_heating_rate",
    "use_specific_heating_rate": "use_specific_heating_rate",
    "self_shielding_method": "self_shielding_method",
    "H2_self_shielding": "H2_self_shielding",
    "grackle_data_file": "grackle_data_file",
    "UVbackground": "UVbackground",
    "Compton_xray_heating": "Compton_xray_heating",
    "LWbackground_intensity": "LWbackground_intensity",
    "LWbackground_sawtooth_suppression": "LWbackground_sawtooth_suppression"
}

_field_map = {
    'density': (('gas', 'density'), 'code_mass / code_length**3'),
    'HI_density': (('gas', 'H_p0_density'), 'code_mass / code_length**3'),
    'HII_density': (('gas', 'H_p1_density'), 'code_mass / code_length**3'),
    'HM_density': (('gas', 'H_m1_density'), 'code_mass / code_length**3'),
    'HeI_density': (('gas', 'He_p0_density'), 'code_mass / code_length**3'),
    'HeII_density': (('gas', 'He_p1_density'), 'code_mass / code_length**3'),
    'HeIII_density': (('gas', 'He_p2_density'), 'code_mass / code_length**3'),
    'H2I_density': (('gas', 'H2_p0_density'), 'code_mass / code_length**3'),
    'H2II_density': (('gas', 'H2_p1_density'), 'code_mass / code_length**3'),
    'DI_density': (('gas', 'D_p0_density'), 'code_mass / code_length**3'),
    'DII_density': (('gas', 'D_p1_density'), 'code_mass / code_length**3'),
    'HDI_density': (('gas', 'HD_p0_density'), 'code_mass / code_length**3'),
    'e_density': (('gas', 'El_density'), 'code_mass / code_length**3'),
    'metal_density': (('gas', 'total_metal_density'), 'code_mass / code_length**3'),
    'dust_density': (('gas', 'dust_density'), 'code_mass / code_length**3'),
    'x_velocity': (('gas', 'velocity_x'), 'code_velocity'),
    'y_velocity': (('gas', 'velocity_y'), 'code_velocity'),
    'z_velocity': (('gas', 'velocity_z'), 'code_velocity'),
    'internal_energy': (('gas', 'specific_thermal_energy'), 'code_velocity**2'),
    'RT_heating_rate': (('gas', 'photo_gamma'), 'erg/s')
}

def _data_to_fc(data, size=None, fc=None):
    if size is None:
        size = data['gas', 'density'].size
    if fc is None:
        fc = FluidContainer(data.ds.grackle_data, size)

    flatten = len(data['gas', 'density'].shape) > 1

    for gfield in fc.input_fields:
        yfield, units = _field_map[gfield]
        fdata = data[yfield].to(units)

        if flatten:
            fdata = fdata.flatten()
        fc[gfield][:] = fdata

    if 'e_density' in fc:
        fc['e_density'] *= (mp/me)

    return fc

def prepare_grackle_data(ds, parameters=None):
    sim_type = type(ds)
    par_map = _parameter_map.get(sim_type)
    if par_map is None:
        print (f"Warning: cannot get Grackle parameters from {sim_type}.\n"
               "All parameters need to be supplied through the \"parameters\" keyword.")
        par_map = {}

    all_parameters = \
      dict([(gpar, ds.parameters[dpar])
            for gpar, dpar in par_map.items()
            if dpar in ds.parameters])
    all_parameters['use_grackle'] = 1

    if parameters is None:
        parameters = {}
    all_parameters.update(parameters)

    my_chemistry = chemistry_data()
    for gpar, val in all_parameters.items():
        if val is None:
            continue
        if isinstance(val, str):
            sval = bytes(val, 'utf-8')
            setattr(my_chemistry, gpar, sval)
        else:
            setattr(my_chemistry, gpar, val)

    my_chemistry.comoving_coordinates = ds.cosmological_simulation
    my_chemistry.density_units = (ds.mass_unit / ds.length_unit**3).in_cgs().d
    my_chemistry.length_units = ds.length_unit.in_cgs().d
    my_chemistry.time_units = ds.time_unit.in_cgs().d
    my_chemistry.a_units = 1 / (1 + ds.parameters.get('CosmologyInitialRedshift', 0))
    my_chemistry.a_value = 1 / (1 + ds.current_redshift) / my_chemistry.a_units
    my_chemistry.velocity_units = ds.velocity_unit.in_cgs().d
    my_chemistry.initialize()
    ds.grackle_data = my_chemistry


_grackle_fields = {
    'cooling_time': 'code_time',
    'dust_temperature': 'K',
    'gamma': '',
    'mean_molecular_weight': '',
    'pressure': 'code_mass * code_velocity**2 / code_length**3',
    'temperature': 'K',
    }

def _grackle_field(field, data):
    gfield = field.name[1][len("grackle_"):]
    units = _grackle_fields[gfield]

    if not hasattr(data.ds, "grackle_data"):
        raise RuntimeError("Grackle has not been initialized.")

    fc = _data_to_fc(data)
    if isinstance(data, FieldDetector):
        if gfield not in fc:
            fc._setup_fluid(gfield)
    else:
        getattr(fc, f"calculate_{gfield}")()

    fdata = fc[gfield]
    if hasattr(data, 'ActiveDimensions'):
        fdata = fdata.reshape(data.ActiveDimensions)

    return fdata * data.ds.quan(1, units).in_cgs()

def _total_metal_density(field, data):
    field_data = data.ds.arr(np.zeros(data['index', 'ones'].shape),
                             'code_mass / code_length**3')
    fields = [
        ("enzo", "Metal_Density"),
        ("enzo", "SN_Colour")]

    for field in fields:
        if field not in data.ds.field_list:
            continue
        field_data += data[field]
    return field_data

def add_grackle_fields(ds, parameters=None):
    ds.add_field(('gas', 'total_metal_density'),
                 function=_total_metal_density,
                 units='g/cm**3',
                 sampling_type='cell')

    prepare_grackle_data(ds, parameters=parameters)
    for field, units in _grackle_fields.items():
        fname = f"grackle_{field}"
        funits = str(ds.quan(1, units).in_cgs().units)
        ds.add_field(('gas', fname), function=_grackle_field,
                     sampling_type="cell", units=funits)
