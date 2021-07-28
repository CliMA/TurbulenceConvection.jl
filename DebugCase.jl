# Compare TurbulenceConvection.jl with SCAMPy
using TurbulenceConvection
# Need to define case_name, first. e.g., `case_name = "Nieuwstadt"`
tc_dir = dirname(dirname(pathof(TurbulenceConvection)));
include(joinpath("integration_tests", "utils", "compute_mse.jl"))

ds_filename = joinpath(tc_dir, "Output.$(case_name).01", "stats", "Stats.$(case_name).nc")

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
        Base.run(`python generate_namelist.py $(case_name)`)
        Base.run(`python generate_paramlist.py $(case_name)`)
        Base.run(`python setup.py build_ext --inplace`)
        Base.run(`python main.py $(case_name).in paramlist_$(case_name).in`)
    end
catch
    println("catch")
end
best_mse = OrderedDict()
best_mse["rho0"] = 1.0
best_mse["p0"] = 1.0
best_mse["alpha0"] = 1.0
best_mse["rho0_half"] = 1.0
best_mse["p0_half"] = 1.0
best_mse["alpha0_half"] = 1.0
best_mse["qt_mean"] = 1.0
best_mse["ql_mean"] = 1.0
best_mse["updraft_area"] = 1.0
best_mse["updraft_w"] = 1.0
best_mse["updraft_qt"] = 1.0
best_mse["updraft_thetal"] = 1.0
best_mse["v_mean"] = 1.0
best_mse["u_mean"] = 1.0
best_mse["tke_mean"] = 1.0
best_mse["tke_values"] = 1.0

best_mse["updraft_buoyancy_values"] = 1.0
best_mse["updraft_buoyancy"] = 1.0
best_mse["eddy_viscosity"] = 1.0
best_mse["mixing_length"] = 1.0
best_mse["entrainment_sc"] = 1.0
best_mse["detrainment_sc"] = 1.0
best_mse["rad_dTdt"] = 1.0
best_mse["subsidence"] = 1.0
best_mse["QT_tendencies"] = 1.0


best_mse["T_values"] = 1.0
best_mse["QT_values"] = 1.0
best_mse["QL_values"] = 1.0
best_mse["RH_values"] = 1.0
best_mse["H_values"] = 1.0
best_mse["Area_values"] = 1.0
best_mse["W_values"] = 1.0
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

ds_scm_filename = joinpath(scampy_dir, "Output.$(case_name).01", "stats", "Stats.$(case_name).nc")
# ds_scm_filename = joinpath(SCAMPy_output_dataset_path, "$(case_name).nc")

computed_mse = Dataset(ds_filename, "r") do ds
    Dataset(ds_scm_filename, "r") do ds_scampy
        compute_mse(
            case_name,
            best_mse,
            joinpath(dirname(ds_filename), "Debugging_$(case_name)_comparison");
            ds_turb_conv=ds,
            ds_scampy=ds_scampy,
            plot_comparison=true
        )
    end
end

nothing
