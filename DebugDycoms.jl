# Compare TurbulenceConvection.jl with SCAMPy
using TurbulenceConvection
tc_dir = dirname(dirname(pathof(TurbulenceConvection)));
include(joinpath("integration_tests", "utils", "compute_mse.jl"))

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
best_mse["rho0"] = 3.3444957379390586e-02
best_mse["p0"] = 3.3444957379390586e-02
best_mse["alpha0"] = 3.3444957379390586e-02
best_mse["rho0_half"] = 3.3444957379390586e-02
best_mse["p0_half"] = 3.3444957379390586e-02
best_mse["alpha0_half"] = 3.3444957379390586e-02
best_mse["qt_mean"] = 3.4203549443317979e-02
best_mse["ql_mean"] = 1.8556088804370651e+02
best_mse["updraft_area"] = 2.1254899281263803e+02
best_mse["updraft_w"] = 3.3720682408389879e+00
best_mse["updraft_qt"] = 7.3758208297447347e-01
best_mse["updraft_thetal"] = 1.2774041797120240e+01
best_mse["v_mean"] = 3.9947573428069575e+01
best_mse["u_mean"] = 3.5718662854892038e+01
best_mse["tke_mean"] = 2.9756050060113374e+01

best_mse["updraft_buoyancy"] = 1.0
best_mse["eddy_viscosity"] = 1.0
best_mse["mixing_length"] = 1.0
best_mse["entrainment_sc"] = 1.0
best_mse["detrainment_sc"] = 1.0
best_mse["rad_dTdt"] = 1.0
best_mse["subsidence"] = 1.0
best_mse["QT_tendencies"] = 1.0


best_mse["RH_mean"] = 1.0
best_mse["thetal_mean"] = 1.0
best_mse["temperature_mean"] = 1.0
best_mse["buoyancy_mean"] = 1.0
best_mse["Hvar_mean"] = 1.0
best_mse["QTvar_mean"] = 1.0
best_mse["HQTcov_mean"] = 1.0
best_mse["W_third_m"] = 1.0
best_mse["H_third_m"] = 1.0
best_mse["QT_third_m"] = 1.0
best_mse["cloud_fraction_mean"] = 1.0
best_mse["QT_tendencies"] = 1.0
best_mse["rad_dTdt"] = 1.0
best_mse["subsidence"] = 1.0
best_mse["rad_flux"] = 1.0
best_mse["updraft_area"] = 1.0
best_mse["updraft_w"] = 1.0
best_mse["updraft_qt"] = 1.0
best_mse["updraft_ql"] = 1.0
best_mse["updraft_RH"] = 1.0
best_mse["updraft_thetal"] = 1.0
best_mse["updraft_temperature"] = 1.0
best_mse["updraft_buoyancy"] = 1.0
best_mse["updraft_cloud_fraction"] = 1.0
best_mse["env_w"] = 1.0
best_mse["env_qt"] = 1.0
best_mse["env_ql"] = 1.0
best_mse["env_area"] = 1.0
best_mse["env_temperature"] = 1.0
best_mse["env_RH"] = 1.0
best_mse["env_thetal"] = 1.0
best_mse["env_tke"] = 1.0
best_mse["env_Hvar"] = 1.0
best_mse["env_QTvar"] = 1.0
best_mse["env_HQTcov"] = 1.0
best_mse["env_cloud_fraction"] = 1.0
best_mse["qr_mean"] = 1.0
best_mse["updraft_qr"] = 1.0
best_mse["env_qr"] = 1.0
best_mse["rain_area"] = 1.0
best_mse["updraft_rain_area"] = 1.0
best_mse["env_rain_area"] = 1.0
best_mse["eddy_viscosity"] = 1.0
best_mse["eddy_diffusivity"] = 1.0
best_mse["entrainment_sc"] = 1.0
best_mse["detrainment_sc"] = 1.0
best_mse["nh_pressure"] = 1.0
best_mse["nh_pressure_adv"] = 1.0
best_mse["nh_pressure_drag"] = 1.0
best_mse["nh_pressure_b"] = 1.0
best_mse["asp_ratio"] = 1.0
best_mse["b_coeff"] = 1.0
best_mse["horizontal_KM"] = 1.0
best_mse["horizontal_KH"] = 1.0
best_mse["sorting_function"] = 1.0
best_mse["b_mix"] = 1.0
best_mse["turbulent_entrainment"] = 1.0
best_mse["turbulent_entrainment_full"] = 1.0
best_mse["turbulent_entrainment_W"] = 1.0
best_mse["turbulent_entrainment_H"] = 1.0
best_mse["turbulent_entrainment_QT"] = 1.0
best_mse["massflux"] = 1.0
best_mse["massflux_h"] = 1.0
best_mse["massflux_qt"] = 1.0
best_mse["massflux_tendency_h"] = 1.0
best_mse["massflux_tendency_qt"] = 1.0
best_mse["diffusive_flux_h"] = 1.0
best_mse["diffusive_flux_u"] = 1.0
best_mse["diffusive_flux_v"] = 1.0
best_mse["diffusive_flux_qt"] = 1.0
best_mse["diffusive_tendency_h"] = 1.0
best_mse["diffusive_tendency_qt"] = 1.0
best_mse["total_flux_h"] = 1.0
best_mse["total_flux_qt"] = 1.0
best_mse["mixing_length"] = 1.0
best_mse["updraft_qt_precip"] = 1.0
best_mse["updraft_thetal_precip"] = 1.0
best_mse["ed_length_scheme"] = 1.0
best_mse["mixing_length_ratio"] = 1.0
best_mse["entdet_balance_length"] = 1.0
best_mse["interdomain_tke_t"] = 1.0
best_mse["tke_buoy"] = 1.0
best_mse["tke_dissipation"] = 1.0
best_mse["tke_entr_gain"] = 1.0
best_mse["tke_detr_loss"] = 1.0
best_mse["tke_shear"] = 1.0
best_mse["tke_pressure"] = 1.0
best_mse["tke_interdomain"] = 1.0
best_mse["tke_transport"] = 1.0
best_mse["tke_advection"] = 1.0
best_mse["Hvar_dissipation"] = 1.0
best_mse["QTvar_dissipation"] = 1.0
best_mse["HQTcov_dissipation"] = 1.0
best_mse["Hvar_entr_gain"] = 1.0
best_mse["QTvar_entr_gain"] = 1.0
best_mse["Hvar_detr_loss"] = 1.0
best_mse["QTvar_detr_loss"] = 1.0
best_mse["HQTcov_detr_loss"] = 1.0
best_mse["HQTcov_entr_gain"] = 1.0
best_mse["Hvar_shear"] = 1.0
best_mse["QTvar_shear"] = 1.0
best_mse["HQTcov_shear"] = 1.0
best_mse["Hvar_rain"] = 1.0
best_mse["QTvar_rain"] = 1.0
best_mse["HQTcov_rain"] = 1.0
best_mse["Hvar_interdomain"] = 1.0
best_mse["QTvar_interdomain"] = 1.0
best_mse["HQTcov_interdomain"] = 1.0
best_mse["prec_source_qt"] = 1.0
best_mse["rain_evap_source_qt"] = 1.0
best_mse["QT_mf_update"] = 1.0
best_mse["QT_new"] = 1.0

ds_scm_filename = joinpath(scampy_dir, "Output.DYCOMS_RF01.01", "stats", "Stats.DYCOMS_RF01.nc")

computed_mse = Dataset(ds_filename, "r") do ds
    Dataset(ds_scm_filename, "r") do ds_scampy
        compute_mse(
            "DYCOMS_RF01",
            best_mse,
            joinpath(dirname(ds_filename), "DebuggingDycoms_comparison");
            ds_turb_conv=ds,
            ds_scampy=ds_scampy,
            plot_comparison=false
        )
    end
end

nothing
