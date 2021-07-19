module ParamList
# See Table 1 of Tan et al, 2018
#paramlist["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] ==> c_k (scaling constant for eddy diffusivity/viscosity
#paramlist["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] == > c_e (scaling constant for tke dissipation)
#paramlist["turbulence"]["EDMF_PrognosticTKE"]["pressure_buoy_coeff"] ==> alpha_b (scaling constant for virtual mass term)
#paramlist["turbulence"]["EDMF_PrognosticTKE"]["pressure_drag_coeff"] ==> alpha_d (scaling constant for drag term)
# paramlist["turbulence"]["EDMF_PrognosticTKE"]["pressure_plume_spacing"] ==> r_d (horizontal length scale of plume spacing)

# Parameters below can be used to multiply any entrainment rate for quick tuning/experimentation
# (NOTE: these are not c_epsilon, c_delta,0 defined in Tan et al 2018)
# paramlist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.1
# paramlist["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 1.0

#NB: except for Bomex and life_cycle_Tan2018 cases, the parameters listed have not been thoroughly tuned/tested
# and should be regarded as placeholders only. Optimal parameters may also depend on namelist options, such as
# entrainment/detrainment rate formulation, diagnostic vs. prognostic updrafts, and vertical resolution
export default_paramlist

using ArgParse
import JSON

function parse_commandline()
    s = ArgParseSettings(; description = "Paramlist Generator")

    @add_arg_table! s begin
        "case_name"
        help = "The case name"
        arg_type = String
        required = true
    end

    return parse_args(s)
end

function default_paramlist(::Nothing)

    args = parse_commandline()
    case_name = args["case_name"]
    return default_paramlist(case_name)
end

function default_paramlist(case_name::String)

    paramlist_defaults = Dict()
    paramlist_defaults["meta"] = Dict()

    paramlist_defaults["turbulence"] = Dict()
    paramlist_defaults["turbulence"]["Ri_bulk_crit"] = 0.2
    paramlist_defaults["turbulence"]["prandtl_number_0"] = 0.74

    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"] = Dict()
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.1
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.14
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.22
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.4
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["lambda_stab"] = 0.9
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = 0.9
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.13
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 0.51
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.4
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.015
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_ed_mf_sigma"] = 50.0
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.25
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_sigma"] = 10.0
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.004
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] = 2.0
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["aspect_ratio"] = 0.2
    # This constant_plume_spacing corresponds to plume_spacing/alpha_d in the Tan et al paper,
    #with values plume_spacing=500.0, alpha_d = 0.375
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["constant_plume_spacing"] = 1333.0
    # TODO: merge the tan18 buoyancy forluma into normalmode formula -> simply set buoy_coeff1 as 1./3. and buoy_coeff2 as 0.
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_buoy_coeff"] = 1.0 / 3.0

    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] = 0.0
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.1
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 10.0

    if case_name == "Soares"
        paramlist = Soares(paramlist_defaults)
    elseif case_name == "Nieuwstadt"
        paramlist = Nieuwstadt(paramlist_defaults)
    elseif case_name == "Bomex"
        paramlist = Bomex(paramlist_defaults)
    elseif case_name == "life_cycle_Tan2018"
        paramlist = life_cycle_Tan2018(paramlist_defaults)
    elseif case_name == "Rico"
        paramlist = Rico(paramlist_defaults)
    elseif case_name == "TRMM_LBA"
        paramlist = TRMM_LBA(paramlist_defaults)
    elseif case_name == "ARM_SGP"
        paramlist = ARM_SGP(paramlist_defaults)
    elseif case_name == "GATE_III"
        paramlist = GATE_III(paramlist_defaults)
    elseif case_name == "DYCOMS_RF01"
        paramlist = DYCOMS_RF01(paramlist_defaults)
    elseif case_name == "GABLS"
        paramlist = GABLS(paramlist_defaults)
    elseif case_name == "SP"
        paramlist = SP(paramlist_defaults)
    elseif case_name == "DryBubble"
        paramlist = DryBubble(paramlist_defaults)
    else
        error("Not a valid case name")
    end

    write_file(paramlist)
    return paramlist
end
function Soares(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "Soares"

    return paramlist
end
function Nieuwstadt(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "Nieuwstadt"

    return paramlist
end
function Bomex(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "Bomex"

    return paramlist
end
function life_cycle_Tan2018(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "life_cycle_Tan2018"

    return paramlist
end
function Rico(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)

    paramlist["meta"]["casename"] = "Rico"

    return paramlist
end
function TRMM_LBA(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)

    paramlist["meta"]["casename"] = "TRMM_LBA"

    return paramlist
end
function ARM_SGP(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "ARM_SGP"

    return paramlist
end
function GATE_III(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "GATE_III"

    return paramlist
end
function DYCOMS_RF01(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)

    paramlist["meta"]["casename"] = "DYCOMS_RF01"

    return paramlist
end
function GABLS(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)

    paramlist["meta"]["casename"] = "GABLS"
    return paramlist
end
# Not fully implemented yet - Ignacio
function SP(paramlist_defaults)

    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "SP"

    return paramlist
end
function DryBubble(paramlist_defaults)
    paramlist = deepcopy(paramlist_defaults)
    paramlist["meta"]["casename"] = "DryBubble"
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.0

    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.25
    paramlist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 0.1

    return paramlist
end

function write_file(paramlist)

    open("paramlist_" * paramlist["meta"]["casename"] * ".in", "w") do io
        JSON.print(io, paramlist, 4)
    end

    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_paramlist(nothing)
end

end
