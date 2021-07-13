import sys
sys.path.insert(0, "./../")
sys.path.insert(0, "./")

import os
import subprocess
import json
import warnings

from netCDF4 import Dataset

import pytest
import pprint as pp


import main as TurbulenceConvection
import common as cmn

# list of possible test cases
case_list = ["Bomex", "life_cycle_Tan2018", "Soares", "Rico", "TRMM_LBA",
             "ARM_SGP", "GATE_III", "DYCOMS_RF01", "GABLS", "SP"]

@pytest.fixture(scope="module")
def data(request):

    # dictionary where simulation results will be stored
    data = {}

    # loop over all test cases
    for case in case_list:

        # generate namelist and paramlist
        setup = cmn.simulation_setup(case)

        # run for 2 hours, output only at the end
        setup["namelist"]["time_stepping"]["t_max"] = 2*60*60
        setup["namelist"]["stats_io"]["frequency"] = setup["namelist"]["time_stepping"]["t_max"]

        # run TurbulenceConvection
        TurbulenceConvection.main1d(setup["namelist"], setup["paramlist"])

        # simulation results
        data[case] = Dataset(setup["outfile"], "r")

    request.addfinalizer(cmn.removing_files)

    return data

@pytest.mark.parametrize("case", case_list)
def test_mean_qt_after_2hr(data, case):
    """
    Check if the mean qt is equal to updraft_area * updraft_qt + (1 - updraft_area) * env_q
    """

    eps_dict = {"SP":                 7e-3,\
                "DYCOMS_RF01":        5e-3,\
                "GABLS":              5e-3,\
                "TRMM_LBA":           4e-3,\
                "life_cycle_Tan2018": 2e-3,\
                "Soares":             2e-3,\
                "ARM_SGP":            7e-4,\
                "GATE_III":           6e-4,\
                "Bomex":              4e-4,\
                "Rico":               4e-4}

    # read in the data
    qt_mean  = np.array(data[str(case)]["profiles/qt_mean"][:,:])
    udr_area = np.array(data[str(case)]["profiles/updraft_area"][:,:])
    udr_qt   = np.array(data[str(case)]["profiles/updraft_qt"][:,:])
    env_qt   = np.array(data[str(case)]["profiles/env_qt"][:,:])

    # calculate the mean qt from updraft and environmet means
    tmp_mean = udr_area * udr_qt + (1 - udr_area) * env_qt

    # absolute error
    abs_error = np.abs(qt_mean - tmp_mean)
    # relative error
    rel_error = np.zeros_like(abs_error)
    idx_not_zero = tmp_mean != 0
    np.place(rel_error, idx_not_zero, abs_error[idx_not_zero] / tmp_mean[idx_not_zero])

    # output for humans
    print "max_abs_error = ", np.max(abs_error), "max_rel_error = ", np.max(rel_error)

    # test
    assert np.allclose(qt_mean[end,:], tmp_mean[end,:], rtol = eps_dict[case], atol=0),\
           "qt_mean != updraft_area * updraft_qt_mean + (1-updraft_area) * env_qt_mean"
