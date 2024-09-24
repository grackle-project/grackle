########################################################################
#
# Tests rates calculated when chemistry is initialized.
#
# Written by Ewan Jones, 19/03/2021
#
########################################################################

#Standard modules
import h5py
import os

#Chemistry_data struct from grackle
from pygrackle import chemistry_data
#Necessary constants from grackle
from pygrackle.utilities.physical_constants import mass_hydrogen_cgs
from pygrackle.utilities.testing import assert_allclose

#* Function which returns chemistry_data instance with default initialisation settings.
def get_defChem():
    #* Initialise chemistry_data instance
    my_chemistry = chemistry_data()

    #* Set parameters
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 0
    my_chemistry.primordial_chemistry = 1
    my_chemistry.metal_cooling = 0
    my_chemistry.UVbackground = 0
    my_chemistry.comoving_coordinates = 0
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0
    my_chemistry.density_units = mass_hydrogen_cgs
    my_chemistry.length_units = 1.0
    my_chemistry.time_units = 1.0
    my_chemistry.velocity_units = my_chemistry.length_units / my_chemistry.time_units

    return my_chemistry


#* Function which sets parameters for a given parameter set.
def set_parameters(parSet, my_chemistry):
    #Default parameter set.
    my_chemistry.CaseBRecombination = 0
    my_chemistry.h2_charge_exchange_rate = 1
    my_chemistry.h2_dust_rate = 1
    my_chemistry.collisional_excitation_rates = 1
    my_chemistry.collisional_ionisation_rates = 1
    my_chemistry.recombination_cooling_rates = 1
    my_chemistry.bremsstrahlung_cooling_rates = 1
    my_chemistry.h2_h_cooling_rate = 1
    my_chemistry.photoelectric_heating = -1
    my_chemistry.dust_chemistry = 0
    my_chemistry.three_body_rate = 0
    if parSet == 1:
        return True
    #Alternate parameter set.
    elif parSet == 2:
        my_chemistry.CaseBRecombination = 1
        my_chemistry.h2_charge_exchange_rate = 2
        my_chemistry.three_body_rate = 1
        my_chemistry.h2_dust_rate = 0
        my_chemistry.collisional_excitation_rates = 0
        my_chemistry.collisional_ionisation_rates = 0
        my_chemistry.recombination_cooling_rates = 0
        my_chemistry.bremsstrahlung_cooling_rates = 0
        my_chemistry.h2_h_cooling_rate = 2
        my_chemistry.photoelectric_heating = 2
        return True
    #Default parameter sets with other three_body_rate values.
    #Set 3 checks caseBrecombination rates.
    elif parSet == 3:
        my_chemistry.three_body_rate = 2
        my_chemistry.CaseBRecombination = 1
        my_chemistry.recombination_cooling_rates = 1
        return True
    elif parSet == 4:
        my_chemistry.three_body_rate = 3
        return True
    elif parSet == 5:
        my_chemistry.three_body_rate = 4
        return True
    elif parSet == 6:
        my_chemistry.three_body_rate = 5
        return True
    # Default parameter set with dust chemistry enabled.
    elif parSet == 7:
        #Load required data file.
        current_path = os.path.abspath(__file__)
        dirs = "..,..,..,grackle_data_files,input,cloudy_metals_2008_3D.h5".split(",")
        data_file_path = os.path.join(os.path.dirname(current_path), *dirs)
        my_chemistry.grackle_data_file = data_file_path
        #Set parameters.
        my_chemistry.dust_chemistry = 1
        my_chemistry.metal_cooling = 1
        return True
    #Invalid parameter set.
    else:
        return False


#* Function which tests that the rates have been initialised correctly for each parameter set.
def test_rate_initialisation(answertestspec,
                             printParameters=False,
                             printOOMdiscrepanices=False,
                             testCustomFile=False,
                             parSets=[1,2,3,4,5,6,7],
                             fileName="rate_coefficients.h5"):
    """
    Test that the rate tables are initialized correctly.

    If writeHDF5 is set to true, all rate coefficients calculated in the initialisation routine will
    be saved into a hdf5 file. This should be used for testing and is not meant for frequent use.

    printParameters (bool) --> If set to True will print the parameter settings for each parameter set.
    """

    #* Navigate to the directory where the file is located.
    filePath = os.path.abspath(__file__)
    dirPath = os.path.dirname(filePath)
    os.chdir(dirPath)

    #* List of all rate variable names which will be checked.
    testRates = "k1,k3,k4,k2,k5,k6,k7,k8,k9,k10,k11,k12,k14,k15,k16,k17,k18,k19,k20,k23,"\
                "k13dd,k13,k21,k22,k50,k51,k52,k53,k54,k55,k56,k57,k58,h2dust,n_cr_n,n_cr_d1,"\
                "n_cr_d2,ceHI,ceHeI,ceHeII,ciHeIS,ciHI,ciHeI,ciHeII,reHII,reHeII1,reHeII2,"\
                "reHeIII,brem,vibh,hyd01k,h2k01,rotl,roth,GP99LowDensityLimit,GP99HighDensityLimit,"\
                "GAHI,GAH2,GAHe,GAHp,GAel,H2LTE,HDlte,HDlow,cieco,comp,gammah,"\
                "gamma_isrf,gas_grain,regr".split(',')

    #* Calculate rates for each parameter set and write to hdf5 file
    #Create and open file. If the file already exists this will overwrite it.
    if answertestspec.generate_answers:
        fileName = os.path.join(answertestspec.answer_dir, fileName)
    f = h5py.File(fileName, "w")

    #Iterate over parameter sets.
    for parSet in parSets:
        my_chemistry = get_defChem()
        #Set chemistry parameters.
        if not set_parameters(parSet, my_chemistry):
            raise RuntimeError("Invalid parameter set encountered.")
        
        #Initialise the rate coefficients.
        if not my_chemistry.initialize():
            raise RuntimeError("Failed to initialize chemistry_data")

        #If the parameter set is 7 (dust enabled) then only write those that depend on the dust setting.
        if parSet == 7:
            f.create_dataset("gas_grain", data=getattr(my_chemistry, "gas_grain"))
            f.create_dataset("regr", data=getattr(my_chemistry, "regr"))
        else:
            #Write rates to file.
            for rate_key in testRates:
                if rate_key not in "regr,gas_grain".split(","):
                    f.create_dataset(rate_key + f"_{parSet}", data=getattr(my_chemistry, rate_key))

    #Close the file.
    f.close()

    # Just generate results and leave.
    if answertestspec.generate_answers:
        return

    #* Compare rates with the expected (correct) ones which are stored and check they are in agreement
    expectedRates = h5py.File(os.path.join(answertestspec.answer_dir, fileName), "r")
    initialisedRates = h5py.File(fileName, "r")


    #*Check all rates for each parameter set.
    #All rates except dust parameters
    for rate_key in testRates:
        for parSet in parSets:
            #Check dust-enabled parameters individually and only once.
            if rate_key in "regr,gas_grain".split(",") or parSet == 7:
                None
            else:
                rate_name = rate_key + f"_{parSet}"
                #Check rates agree to what we deem is an acceptable relative tolerance.
                assert_allclose(expectedRates[rate_name][()], initialisedRates[rate_name][()], rtol=1e-7, atol=0,
                                err_msg=f"Rate Coefficients for {rate_name} do not agree.")
                                    
    #Dust parameters.
    for rate_key in "regr,gas_grain".split(","):
        assert_allclose(expectedRates[rate_key][()], initialisedRates[rate_key][()], rtol=1e-7, atol=0,
                                err_msg=f"Rate Coefficients for {rate_key} do not agree.")
        

    #Close files.
    expectedRates.close()
    initialisedRates.close()
