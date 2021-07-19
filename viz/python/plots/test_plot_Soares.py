import sys
sys.path.insert(0, "./")
sys.path.insert(0, "../")

import os
import subprocess
from pathlib import Path

from netCDF4 import Dataset

import pytest

import main as TurbulenceConvection
import common as cmn
import plot_scripts as pls

@pytest.fixture(scope="module")
def sim_data(request):

    # remove netcdf file from previous failed test
    request.addfinalizer(cmn.removing_files)
    # generate namelists and paramlists
    setup = cmn.simulation_setup("Soares")

    # run TurbulenceConvection
    subprocess.call("python setup.py build_ext --inplace", shell=True, cwd="../")
    TurbulenceConvection.main1d(setup["namelist"], setup["paramlist"])

    # simulation results
    sim_data = Dataset(setup["outfile"], "r")

    # remove netcdf file after tests
    request.addfinalizer(cmn.removing_files)

    return sim_data


def test_plot_Soares(sim_data):
    """
    plot Soares timeseries
    """
    # make directory
    localpath = Path.cwd()
    (localpath / "plots/output/Soares/all_variables/").mkdir(parents=True, exist_ok=True)
    les_data_path = localpath / "les_data/Soares.nc"
    if not les_data_path.is_file():
        url_ = r"https://drive.google.com/uc?export=download&id=1D6DKEGnIFeH9G4y77_BhhZNq1edC9gat"
        os.system(f"curl -sLo {les_data_path} "{url_}"")
    les_data = Dataset(les_data_path, "r")

    f1 = "plots/output/Soares/"
    f2 = f1 + "all_variables/"
    cn = "Soares_"
    t0 = 6
    t1 = 8
    zmin = 0.0
    zmax = 3.0
    cb_min = [0., 0.]
    cb_max = [0.01, 2.5]
    fixed_cbar = True
    cb_min_t = [300, 300, 300, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 4.5,\
                -0.4, 0, -0.7, -0.45,\
                0, -0.03, 0,\
                -0.045, -0.06, 0,\
                0, 0, -0.04]
    cb_max_t = [305, 305, 302, 1, 1, 1, 1, 1, 100, 1, 100, 5, 5, 5,\
                0, 2.4, 0.7, 0.45,\
                0.28, 0.01, 1.6,\
                0.06, 0.045, 0.1,\
                0.01, 0.1, 0.08]

    scm_dict = cmn.read_scm_data(sim_data)
    les_dict = cmn.read_les_data(les_data)

    scm_dict_t = cmn.read_scm_data_timeseries(sim_data)
    les_dict_t = cmn.read_les_data_timeseries(les_data)

    pls.plot_closures(scm_dict, les_dict, t0, t1, zmin, zmax, cn+"closures.pdf", folder=f1)
    pls.plot_spec_hum(scm_dict, les_dict, t0, t1, zmin, zmax,  cn+"humidities.pdf", folder=f1)
    pls.plot_upd_prop(scm_dict, les_dict, t0, t1, zmin, zmax,  cn+"updraft_properties.pdf", folder=f1)
    pls.plot_fluxes(scm_dict, les_dict, t0, t1, zmin, zmax,  cn+"mean_fluxes.pdf", folder=f1)
    pls.plot_tke_comp(scm_dict, les_dict, t0, t1, zmin, zmax,  cn+"tke_components.pdf", folder=f1)

    pls.plot_cvar_mean(scm_dict, les_dict, t0, t1, zmin, zmax,  cn+"var_covar_mean.pdf", folder=f2)
    pls.plot_cvar_comp(scm_dict, t0, t1, zmin, zmax,  cn+"var_covar_components.pdf", folder=f2)
    pls.plot_tke_break(scm_dict, les_dict, t0, t1, zmin, zmax, cn+"tke_breakdown.pdf",folder=f2)

    pls.plot_contour_t(scm_dict, les_dict, fixed_cbar, cb_min_t, cb_max_t, zmin, zmax, folder=f2)
    pls.plot_mean_prof(scm_dict, les_dict, t0, t1,  zmin, zmax, folder=f2)

    pls.plot_main(scm_dict_t, les_dict_t, scm_dict, les_dict,
                  cn+"main_timeseries.pdf", cb_min, cb_max, zmin, zmax, folder=f1)

    pls.plot_1D(scm_dict_t, les_dict_t, cn, folder=f2)
