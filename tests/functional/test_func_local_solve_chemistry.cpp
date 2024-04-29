#include "gtest/gtest.h"

extern "C" {
    #include "grackle.h"
}

TEST(LocalSolveChemistryTest, LocalSolveChemistry) {

    code_units units;
    chemistry_data *chem_data;
    chemistry_data_storage rates;
    grackle_field_data fields;

    int field_size = 1;
    double initial_redshift = 0.;
    double tiny_number = 1.e-20;
    double dt; // Computed below

    grackle_verbose = 0;

    // Set units
    units.comoving_coordinates = 0;
    units.density_units = 1.67e-24;
    units.length_units = 1.0;
    units.time_units = 1.0e12;
    units.a_units = 1.0;
    units.a_value = 1.0;
    set_velocity_units(&units);

    // Create chemistry data variable and set default vaules
    chem_data = new chemistry_data;
    if (set_default_chemistry_parameters(chem_data) == 0) {
        throw std::runtime_error("Error in set_default_chemistry_parameters");
    }

    // Set cheistrydata
    chem_data->use_grackle = 1;
    chem_data->with_radiative_cooling = 1;
    chem_data->primordial_chemistry = 3;
    chem_data->dust_chemistry = 1;
    chem_data->metal_cooling = 1;
    chem_data->UVbackground = 1;
    chem_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5";

    // Run chemistry inizialization
    if (local_initialize_chemistry_data(chem_data, &rates, &units) == 0) {
        throw std::runtime_error("Error in local_initialize_chemistry_data");
    }

    // Set fields data
    fields.grid_rank = 3;
    fields.grid_dimension = new int[fields.grid_rank];
    fields.grid_start = new int[fields.grid_rank];
    fields.grid_end = new int[fields.grid_rank];
    fields.grid_dx = 0.0;
    for(int i=0;i<3;i++) {
        fields.grid_dimension[i] = 1;
        fields.grid_start[i] = 0;
        fields.grid_end[i] = 0;
    }
    fields.grid_dimension[0] = field_size;
    fields.grid_end[0] = field_size - 1;
    fields.density         = new gr_float[field_size];
    fields.internal_energy = new gr_float[field_size];
    fields.x_velocity      = new gr_float[field_size];
    fields.y_velocity      = new gr_float[field_size];
    fields.z_velocity      = new gr_float[field_size];
    fields.HI_density      = new gr_float[field_size];
    fields.HII_density     = new gr_float[field_size];
    fields.HeI_density     = new gr_float[field_size];
    fields.HeII_density    = new gr_float[field_size];
    fields.HeIII_density   = new gr_float[field_size];
    fields.e_density       = new gr_float[field_size];
    fields.HM_density      = new gr_float[field_size];
    fields.H2I_density     = new gr_float[field_size];
    fields.H2II_density    = new gr_float[field_size];
    fields.DI_density      = new gr_float[field_size];
    fields.DII_density     = new gr_float[field_size];
    fields.HDI_density     = new gr_float[field_size];
    fields.metal_density   = new gr_float[field_size];
    fields.volumetric_heating_rate = new gr_float[field_size];
    fields.specific_heating_rate = new gr_float[field_size];
    fields.RT_HI_ionization_rate = new gr_float[field_size];
    fields.RT_HeI_ionization_rate = new gr_float[field_size];
    fields.RT_HeII_ionization_rate = new gr_float[field_size];
    fields.RT_H2_dissociation_rate = new gr_float[field_size];
    fields.RT_heating_rate = new gr_float[field_size];

    double temperature_units = get_temperature_units(&units);

    for(int i=0;i<field_size;i++) {
        fields.density[i] = 1.0;
        fields.HI_density[i] = chem_data->HydrogenFractionByMass*fields.density[i];
        fields.HII_density[i] = tiny_number*fields.density[i];
        fields.HM_density[i] = tiny_number*fields.density[i];
        fields.HeI_density[i] = (1.0 - chem_data->HydrogenFractionByMass)*fields.density[i];
        fields.HeII_density[i] = tiny_number*fields.density[i];
        fields.HeIII_density[i] = tiny_number*fields.density[i];
        fields.H2I_density[i] = tiny_number*fields.density[i];
        fields.H2II_density[i] = tiny_number*fields.density[i];
        fields.DI_density[i] = 2.0*3.4e-5*fields.density[i];
        fields.DII_density[i] = tiny_number*fields.density[i];
        fields.HDI_density[i] = tiny_number*fields.density[i];
        fields.e_density[i] = tiny_number*fields.density[i];
        fields.metal_density[i] = chem_data->SolarMetalFractionByMass*fields.density[i];
        fields.x_velocity[i] = 0.0;
        fields.y_velocity[i] = 0.0;
        fields.z_velocity[i] = 0.0;

        fields.internal_energy[i] = 1000./temperature_units;
        fields.volumetric_heating_rate[i] = 0.0;
        fields.specific_heating_rate[i] = 0.0;

        fields.RT_HI_ionization_rate[i] = 0.0;
        fields.RT_HeI_ionization_rate[i] = 0.0;
        fields.RT_HeII_ionization_rate[i] = 0.0;
        fields.RT_H2_dissociation_rate[i] = 0.0;
        fields.RT_heating_rate[i] = 0.0;
    }

    // Run chemical solver
    dt = 3.15e7*1e6/units.time_units;
    if (local_solve_chemistry(chem_data, &rates, &units, &fields, dt) == 0) {
        throw std::runtime_error("local_solve_chemistry");
    }

    EXPECT_NEAR(fields.HI_density[0]/fields.density[0], 0.630705, 1e-6);
    EXPECT_NEAR(fields.internal_energy[0], 2.95159e+35, 1e30);

    delete [] fields.grid_dimension;
    delete [] fields.grid_start;
    delete [] fields.grid_end;

    delete chem_data;

    delete [] fields.density;
    delete [] fields.internal_energy;
    delete [] fields.x_velocity;
    delete [] fields.y_velocity;
    delete [] fields.z_velocity;
    delete [] fields.HI_density;
    delete [] fields.HII_density;
    delete [] fields.HeI_density;
    delete [] fields.HeII_density;
    delete [] fields.HeIII_density;
    delete [] fields.e_density;
    delete [] fields.HM_density;
    delete [] fields.H2I_density;
    delete [] fields.H2II_density;
    delete [] fields.DI_density;
    delete [] fields.DII_density;
    delete [] fields.HDI_density;
    delete [] fields.metal_density;
    delete [] fields.volumetric_heating_rate;
    delete [] fields.specific_heating_rate;
    delete [] fields.RT_HI_ionization_rate;
    delete [] fields.RT_HeI_ionization_rate;
    delete [] fields.RT_HeII_ionization_rate;
    delete [] fields.RT_H2_dissociation_rate;
    delete [] fields.RT_heating_rate;

}
