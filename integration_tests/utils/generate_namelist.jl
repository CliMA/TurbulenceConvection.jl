module NameList
# See Table ? of Cohen et al, 2020
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.13
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 0.51
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.0
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.015
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.25
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_sigma"] = 10.0
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.004
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] = 2.0
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["aspect_ratio"] = 0.2

# See Table ? of Lopez Gomez et al, 2020
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.14
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.22
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.4
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["lambda_stab"] = 0.9

# See Table ? of He et al, 2021
# namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12 ==> alpha_b (scaling constant for virtual mass term)
# namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] = 0.0 ==> alpha_b (scaling constant for virtual mass term)
# namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.1 alpha_a (scaling constant for advection term)
# namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 10.0 ==> alpha_d (scaling constant for drag term)
# namelist["turbulence"]["EDMF_PrognosticTKE"]["pressure_plume_spacing"] ==> r_d (horizontal length scale of plume spacing)

#NB: except for Bomex and life_cycle_Tan2018 cases, the parameters listed have not been thoroughly tuned/tested
# and should be regarded as placeholders only. Optimal parameters may also depend on namelist options, such as
# entrainment/detrainment rate formulation, diagnostic vs. prognostic updrafts, and vertical resolution
export default_namelist

using ArgParse
import JSON

function parse_commandline()
    s = ArgParseSettings(; description = "namelist Generator")

    @add_arg_table! s begin
        "case_name"
        help = "The case name"
        arg_type = String
        required = true
    end

    return parse_args(s)
end

function default_namelist(::Nothing)

    args = parse_commandline()
    case_name = args["case_name"]
    return default_namelist(case_name)
end

function default_namelist(case_name::String)

    namelist_defaults = Dict()
    namelist_defaults["meta"] = Dict()

    namelist_defaults["turbulence"] = Dict()
    namelist_defaults["turbulence"]["Ri_bulk_crit"] = 0.2
    namelist_defaults["turbulence"]["prandtl_number_0"] = 0.74

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"] = Dict()
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = 0.9
    # mixing_length
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.14
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.22
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.4
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["lambda_stab"] = 0.9
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_scale"] = 53.0 / 13.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_ub"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_rm"] = 1.5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["l_max"] = 1.0e6
    # entrainment
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.13
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 0.51
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.015
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_sigma"] = 10.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.004
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] = 2.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_upd_velocity"] = 0.001
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["aspect_ratio"] = 0.2
    # pressure 
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_updraft_top"] = 500.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 10.0

    # From namelist
    namelist_defaults["grid"] = Dict()
    namelist_defaults["grid"]["dims"] = 1
    namelist_defaults["grid"]["gw"] = 2

    namelist_defaults["thermodynamics"] = Dict()
    namelist_defaults["thermodynamics"]["thermal_variable"] = "thetal"
    namelist_defaults["thermodynamics"]["sgs"] = "quadrature"
    namelist_defaults["thermodynamics"]["quadrature_order"] = 3
    namelist_defaults["thermodynamics"]["quadrature_type"] = "log-normal" #"gaussian" or "log-normal"

    namelist_defaults["time_stepping"] = Dict()
    namelist_defaults["time_stepping"]["dt"] = 3.0

    namelist_defaults["microphysics"] = Dict()
    namelist_defaults["microphysics"]["rain_model"] = "None"

    namelist_defaults["turbulence"]["scheme"] = "EDMF_PrognosticTKE"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = 1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "moisture_deficit"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["extrapolate_buoyancy"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_local_micro"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["constant_area"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["mixing_length"] = "sbtd_eq"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_asp_label"] = "const"

    namelist_defaults["output"] = Dict()
    namelist_defaults["output"]["output_root"] = "./"

    namelist_defaults["stats_io"] = Dict()
    namelist_defaults["stats_io"]["stats_dir"] = "stats"
    namelist_defaults["stats_io"]["frequency"] = 60.0

    if case_name == "Soares"
        namelist = Soares(namelist_defaults)
    elseif case_name == "Nieuwstadt"
        namelist = Nieuwstadt(namelist_defaults)
    elseif case_name == "Bomex"
        namelist = Bomex(namelist_defaults)
    elseif case_name == "life_cycle_Tan2018"
        namelist = life_cycle_Tan2018(namelist_defaults)
    elseif case_name == "Rico"
        namelist = Rico(namelist_defaults)
    elseif case_name == "TRMM_LBA"
        namelist = TRMM_LBA(namelist_defaults)
    elseif case_name == "ARM_SGP"
        namelist = ARM_SGP(namelist_defaults)
    elseif case_name == "GATE_III"
        namelist = GATE_III(namelist_defaults)
    elseif case_name == "DYCOMS_RF01"
        namelist = DYCOMS_RF01(namelist_defaults)
    elseif case_name == "GABLS"
        namelist = GABLS(namelist_defaults)
    elseif case_name == "SP"
        namelist = SP(namelist_defaults)
    elseif case_name == "DryBubble"
        namelist = DryBubble(namelist_defaults)
    else
        error("Not a valid case name")
    end

    write_file(namelist)
    return namelist
end
function Soares(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Soares"

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 8 * 3600.0

    namelist["meta"]["simname"] = "Soares"
    namelist["meta"]["casename"] = "Soares"

    return namelist
end
function Nieuwstadt(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Nieuwstadt"

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 8 * 3600.0

    namelist["meta"]["simname"] = "Nieuwstadt"
    namelist["meta"]["casename"] = "Nieuwstadt"

    return namelist
end
function Bomex(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Bomex"

    namelist["grid"]["nz"] = 60
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 21600.0

    namelist["meta"]["simname"] = "Bomex"
    namelist["meta"]["casename"] = "Bomex"

    return namelist
end
function life_cycle_Tan2018(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "life_cycle_Tan2018"

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 40.0

    namelist["time_stepping"]["t_max"] = 6 * 3600.0
    namelist["meta"]["simname"] = "life_cycle_Tan2018"
    namelist["meta"]["casename"] = "life_cycle_Tan2018"

    return namelist
end
function Rico(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "Rico"

    namelist["grid"]["nz"] = 120
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 86400.0

    # namelist["microphysics"]["rain_model"] = "cutoff"
    namelist["microphysics"]["rain_model"] = "clima_1m"
    namelist["microphysics"]["tau_acnv"] = 1e4

    namelist["meta"]["simname"] = "Rico"
    namelist["meta"]["casename"] = "Rico"
    return namelist
end
function TRMM_LBA(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "TRMM_LBA"

    namelist["grid"]["nz"] = 320
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 21600.0

    namelist["microphysics"]["rain_model"] = "cutoff"
    # namelist["microphysics"]["rain_model"] = "clima_1m"

    namelist["meta"]["simname"] = "TRMM_LBA"
    namelist["meta"]["casename"] = "TRMM_LBA"

    return namelist
end
function ARM_SGP(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "ARM_SGP"

    namelist["grid"]["nz"] = 88
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 3600.0 * 14.5
    namelist["meta"]["simname"] = "ARM_SGP"
    namelist["meta"]["casename"] = "ARM_SGP"

    return namelist
end
function GATE_III(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "GATE_III"

    namelist["grid"]["nz"] = 1700
    namelist["grid"]["dz"] = 10

    namelist["time_stepping"]["t_max"] = 3600.0 * 24.0
    namelist["meta"]["simname"] = "GATE_III"
    namelist["meta"]["casename"] = "GATE_III"

    return namelist
end
function DYCOMS_RF01(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "DYCOMS_RF01"

    namelist["grid"]["nz"] = 30
    namelist["grid"]["dz"] = 50

    namelist["time_stepping"]["t_max"] = 60 * 60 * 16.0
    namelist["meta"]["simname"] = "DYCOMS_RF01"
    namelist["meta"]["casename"] = "DYCOMS_RF01"

    return namelist
end
function GABLS(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "GABLS"

    namelist["grid"]["nz"] = 8
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 9 * 3600.0
    namelist["meta"]["simname"] = "GABLS"
    namelist["meta"]["casename"] = "GABLS"

    return namelist
end
# Not fully implemented yet - Ignacio
function SP(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "SP"

    namelist["grid"]["nz"] = 256
    namelist["grid"]["dz"] = 8

    namelist["time_stepping"]["t_max"] = 7200.0
    namelist["meta"]["simname"] = "SP"
    namelist["meta"]["casename"] = "SP"

    return namelist
end
function DryBubble(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "DryBubble"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.0

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 0.1

    namelist["grid"]["nz"] = 200
    namelist["grid"]["dz"] = 50.0

    namelist["stats_io"]["frequency"] = 10.0
    namelist["time_stepping"]["t_max"] = 1000.0
    namelist["meta"]["simname"] = "DryBubble"
    namelist["meta"]["casename"] = "DryBubble"

    namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.4

    return namelist
end

function write_file(namelist)

    @assert haskey(namelist, "meta")
    @assert haskey(namelist["meta"], "simname")

    namelist["meta"]["uuid"] = basename(tempname())

    open("namelist_" * namelist["meta"]["casename"] * ".in", "w") do io
        JSON.print(io, namelist, 4)
    end

    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_namelist(nothing)
end

end
