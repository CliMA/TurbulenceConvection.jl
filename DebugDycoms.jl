# Compare TurbulenceConvection.jl with SCAMPy
using TurbulenceConvection
tc_dir = dirname(dirname(pathof(TurbulenceConvection)));
include(joinpath("integration_tests", "utils", "compute_mse_scampy.jl"))

ds_filename = joinpath(tc_dir, "Output.DYCOMS_RF01.01", "stats", "Stats.DYCOMS_RF01.nc")

scampy_dir = joinpath(tc_dir, "..", "SCAMPyRepos", "SCAMPyMaster", "SCAMPy")
using Glob
try
    cd(scampy_dir) do
        for f in glob("*.c")
             rm(f)
        end
        for f in glob("*.so")
             rm(f)
        end
        Base.run(`python generate_namelist.py DYCOMS_RF01`)
        Base.run(`python generate_paramlist.py DYCOMS_RF01`)
        Base.run(`python setup.py build_ext --inplace`)
        Base.run(`python main.py DYCOMS_RF01.in paramlist_DYCOMS_RF01.in`)
    end
catch
    println("catch")
end
best_mse = OrderedDict()
# best_mse["rho0"] = 3.3444957379390586e-02 looks good
# best_mse["p0"] = 3.3444957379390586e-02 looks good
# best_mse["alpha0"] = 3.3444957379390586e-02 looks good
# best_mse["rho0_half"] = 3.3444957379390586e-02 looks good
# best_mse["p0_half"] = 3.3444957379390586e-02 looks good
# best_mse["alpha0_half"] = 3.3444957379390586e-02 looks good
best_mse["qt_mean"] = 3.3444957379390586e-02
best_mse["ql_mean"] = 3.3444957379390586e-02
best_mse["updraft_area"] = 2.5801657716690602e+01
best_mse["updraft_w"] = 5.3015266139161090e+00
best_mse["updraft_qt"] = 1.9044631066492126e+00
best_mse["updraft_thetal"] = 4.6206573044439402e+01
best_mse["v_mean"] = 1.9191349011646031e+04
best_mse["u_mean"] = 7.3347017004521811e+04
best_mse["tke_mean"] = 2.2114155090507843e+01

ds_py_filename = joinpath(scampy_dir, "Output.DYCOMS_RF01.01", "stats", "Stats.DYCOMS_RF01.nc")
Dataset(ds_py_filename) do ds_py
    Dataset(ds_filename) do ds_tc

        compute_mse_scampy(
            ds_tc,
            ds_py,
            "DYCOMS_RF01",
            best_mse,
            "DebuggingDycoms";
            plot_comparison=true,
        )
    end
end


