########################################################################
#
# Cooling cell example script
#
#  This will initialize a single cell at a given temperature,
#  iterate the cooling solver for a fixed time, and output the
#  temperature vs. time.
#
#
# Copyright (c) 2015-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import os
import sys
import yt

from gracklepy import \
    chemistry_data, \
    evolve_constant_density, \
    setup_fluid_container
from gracklepy.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

from gracklepy.utilities.data_path import grackle_data_dir
from gracklepy.utilities.model_tests import \
    get_test_variables

_MODEL_NAME = os.path.basename(__file__[:-3]) # strip off ".py"

def gen_plot(fc,  data, fname):
    p1, = pyplot.loglog(data["time"].to("Myr"),
                        data["temperature"],
                        color="black", label="T")
    pyplot.xlabel("Time [Myr]")
    pyplot.ylabel("T [K]")
    pyplot.twinx()
    p2, = pyplot.semilogx(data["time"].to("Myr"),
                          data["mean_molecular_weight"],
                          color="red", label="$\\mu$")
    pyplot.ylabel("$\\mu$")
    pyplot.legend([p1,p2],["T","$\\mu$"], fancybox=True,
                  loc="center left")
    pyplot.tight_layout()
    pyplot.savefig(fname)


def main(args=None):
    args = sys.argv[1:] if args is None else args
    if len(args) != 0:  # we are using the testing framework
        my_vars = get_test_variables(_MODEL_NAME, args)

        metallicity = my_vars["metallicity"]
        redshift = my_vars["redshift"]
        extra_attrs = my_vars["extra_attrs"]
        my_chemistry = my_vars["my_chemistry"]
        output_name = my_vars["output_name"]

        in_testing_framework = True

    else:  # Just run the script as is.
        metallicity = 0.01 # Solar
        redshift = 0.
        # dictionary to store extra information in output dataset
        extra_attrs = {}

        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.grackle_data_file = \
          os.path.join(grackle_data_dir, "cloudy_metals_2008_3D.h5")
        
        def set_opts(obj, **kwargs):
            for k, v in kwargs.items():
                if hasattr(obj, k):
                    setattr(obj, k, v)
                else:
                    print(f"# warning: chemistry_data has no field '{k}' (skipped)")

        set_opts(
            my_chemistry,

            with_radiative_cooling=1,
            primordial_chemistry=4,
            dust_chemistry = 1,
            metal_cooling=1,
            UVbackground=0,
            cmb_temperature_floor=1,
            Gamma=1.66667,
            h2_on_dust=1,
            use_dust_density_field=1,
            metal_chemistry=1,
            multi_metals=0,
            metal_abundances=0,
            dust_species=3,
            # dust_temperature_multi=0,
            dust_sublimation=1,
            grain_growth=0,
            photoelectric_heating=0,
            photoelectric_heating_rate=0,
            use_isrf_field=0,
            interstellar_radiation_field=0,
            use_volumetric_heating_rate=0,
            use_specific_heating_rate=0,
            three_body_rate=1,
            cie_cooling=1,
            h2_optical_depth_approximation=1,
            ih2co=1,
            ipiht=1,

            HydrogenFractionByMass=0.76,
            DeuteriumToHydrogenRatio=6.8e-05,
            SolarMetalFractionByMass=0.01295,
            local_dust_to_gas_ratio=0.009387,
            NumberOfTemperatureBins=600,
            CaseBRecombination=1,
            TemperatureStart=1.0,
            TemperatureEnd=1.0e9,
            NumberOfDustTemperatureBins=250,
            DustTemperatureStart=1.0,
            DustTemperatureEnd=1500.0,
            Compton_xray_heating=0,
            LWbackground_sawtooth_suppression=0,
            LWbackground_intensity=0,
            UVbackground_redshift_on=-99999,
            UVbackground_redshift_off=-99999,
            UVbackground_redshift_fullon=-99999,
            UVbackground_redshift_drop=-99999,
            cloudy_electron_fraction_factor=0.00915396,
            use_radiative_transfer=1,
            radiative_transfer_coupled_rate_solver=0,
            radiative_transfer_intermediate_step=0,
            radiative_transfer_hydrogen_only=0,
            self_shielding_method=0,
            H2_custom_shielding=0,
            H2_self_shielding=0,
            radiative_transfer_H2II_diss=1,
            radiative_transfer_HDI_dissociation=1,
            radiative_transfer_metal_ionization=1,
            radiative_transfer_metal_dissociation=1,

            use_multiple_dust_temperatures=1
        )

        output_name = _MODEL_NAME
        in_testing_framework = False

    density = 0.1 * mass_hydrogen_cgs # g /cm^3
    temperature = 50000 # K
    final_time = 100. # Myr

    # Set units
    my_chemistry.comoving_coordinates = 0
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs
    my_chemistry.length_units = cm_per_mpc
    my_chemistry.time_units = sec_per_Myr
    my_chemistry.set_velocity_units()

    metal_mass_fraction = metallicity * my_chemistry.SolarMetalFractionByMass
    fc = setup_fluid_container(
        my_chemistry,
        density=density,
        temperature=temperature,
        metal_mass_fraction=metal_mass_fraction,
        state="ionized",
        converge=True)

    # evolve gas at constant density
    data = evolve_constant_density(
        fc, final_time=final_time,
        safety_factor=0.01)

    if not in_testing_framework:
        gen_plot(fc, data, fname=f"{output_name}.png")

    # save data arrays as a yt dataset
    yt.save_as_dataset({}, f"{output_name}.h5",
                       data=data, extra_attrs=extra_attrs)
    return 0


if __name__ == '__main__':
    sys.exit(main())