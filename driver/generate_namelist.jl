module NameList
# See Table 2 of Cohen et al, 2020
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.13
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 0.51
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.0
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.075
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.25
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.0004
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] = 2.0

# See Table 1 of Lopez Gomez et al, 2020
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.14
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.22
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.4 # Square of value in the paper
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] = 3.75 # Square of value in the paper
# namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] = 0.74

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

import StaticArrays
const SA = StaticArrays

include(joinpath(@__DIR__, "..", "integration_tests", "artifact_funcs.jl"))

import Random

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

function default_namelist(
    case_name::String;
    root::String = ".",
    write::Bool = true,
    set_seed::Bool = true,
    seed::Int = 2022,
)

    if set_seed
        Random.seed!(seed)
    end

    namelist_defaults = Dict()
    namelist_defaults["meta"] = Dict()
    namelist_defaults["meta"]["uuid"] = basename(tempname())

    namelist_defaults["turbulence"] = Dict()

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"] = Dict()
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = 0.9
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = 1e-5

    # mixing_length
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.14
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.22
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.4
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] = 3.75
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_scale"] = 53.0 / 13.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] = 0.74
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_ub"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_rm"] = 1.5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["l_max"] = 1.0e6
    # entrainment
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.13
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 0.51

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.075
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 10.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 3.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.0004
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] = 2.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_upd_velocity"] = 0.001
    # pressure
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_updraft_top"] = 500.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 10.0

    # From namelist
    namelist_defaults["grid"] = Dict()
    namelist_defaults["grid"]["dims"] = 1
    namelist_defaults["grid"]["stretch"] = Dict()
    namelist_defaults["grid"]["stretch"]["flag"] = false
    namelist_defaults["grid"]["stretch"]["nz"] = 55
    namelist_defaults["grid"]["stretch"]["dz_surf"] = 30.0
    namelist_defaults["grid"]["stretch"]["dz_toa"] = 8000.0
    namelist_defaults["grid"]["stretch"]["z_toa"] = 45000.0

    namelist_defaults["thermodynamics"] = Dict()
    namelist_defaults["thermodynamics"]["thermal_variable"] = "thetal"
    namelist_defaults["thermodynamics"]["moisture_model"] = "equilibrium" #"nonequilibrium"
    namelist_defaults["thermodynamics"]["sgs"] = "quadrature"
    namelist_defaults["thermodynamics"]["quadrature_order"] = 3
    namelist_defaults["thermodynamics"]["quadrature_type"] = "log-normal" #"gaussian" or "log-normal"

    namelist_defaults["time_stepping"] = Dict()
    namelist_defaults["time_stepping"]["dt_max"] = 2.0
    namelist_defaults["time_stepping"]["dt_min"] = 1.0
    namelist_defaults["time_stepping"]["adapt_dt"] = true
    namelist_defaults["time_stepping"]["cfl_limit"] = 0.5

    namelist_defaults["microphysics"] = Dict()
    namelist_defaults["microphysics"]["precipitation_model"] = "None"

    namelist_defaults["turbulence"]["scheme"] = "EDMF_PrognosticTKE"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = 1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "RF"  # {"moisture_deficit", "NN", "NN_nonlocal", "Linear", "FNO", "RF"}
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] = "buoy_vel" # {"buoy_vel", "inv_z", "none"}
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_pi_subset"] = ntuple(i -> i, 4) # or, e.g., (1, 3, 6)
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pi_norm_consts"] = [478.298, 1.0, 1.0, 1.0, 1.0, 1.0] # normalization constants for Pi groups
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] = "deterministic"  # {"deterministic", "noisy_relaxation_process", "lognormal_scaling", "prognostic_noisy_relaxation_process"}
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_local_micro"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["constant_area"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["mixing_length"] = "sbtd_eq"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"] = "normalmode"

    namelist_defaults["output"] = Dict()
    namelist_defaults["output"]["output_root"] = "./"

    namelist_defaults["stats_io"] = Dict()
    namelist_defaults["stats_io"]["stats_dir"] = "stats"
    namelist_defaults["stats_io"]["frequency"] = 60.0
    namelist_defaults["stats_io"]["skip"] = false
    namelist_defaults["stats_io"]["calibrate_io"] = false # limit io for calibration when `true`

    # nn parameters
    #! format: off
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["nn_arc"] = (6, 5, 4, 2) # [#inputs, #neurons in L1, #neurons in L2, ...., #outputs]

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["nn_ent_params"] =
        SA.SVector{58}(rand(58))

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["nn_ent_biases"] = false

    # m=100 random features, d=6 input Pi groups
    # RF: parameters to optimize, 2 x (m + 1 + d)
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["rf_opt_ent_params"] =
        vec([-0.37202122293049483, -0.6918575291754739, 0.35265093911995316, -1.2676542400276773, 0.3634954695178507, -1.2152100958722636, -0.7318865530445821, -0.20490609825171713, 0.4362212821398628, 1.9229172469951967, -1.8462922138249722, 0.5284051268096007, 0.8834774018182467, 1.3981008690984473,
             0.987369919895047, -0.8185539483252741, 0.5213957436573822, -0.9383219360570113, -0.16980975651279465, -0.7893995981968782, 0.7593451662499965, -1.941397120422341, 0.32401189333909924, 1.182500970485504, 0.6162448000852512, -0.12412174568480978, -0.7484836738586393, -1.2318976813652487,
             -0.43927100837036703, -2.6124394102311834, 0.378458816357903, 0.17722284611138117, 0.8352945611841436, 1.0594900911903236, -0.007101450661217078, -0.23449146053875652, -0.04848665714763242, -0.1906590175353914, -1.6413139801382026, -0.5862430780400537, -0.09651442030984692, -1.4334941737354499, -0.3648667483936556, -0.6480311432494812,
             1.6618153118350076, 2.8264791354231806, 1.1846562064154362, -0.5712779601534338, -0.7952460470546551, 2.0877871984477046
            ])
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["rf_fix_ent_params"] = vec([3.671073078569706,
                1.8372209890762072, 3.9228808428598323, 5.182737854497501, 5.895379472008546, 0.9263736906562064, 2.4098464460423887, 3.388135105376523, 2.0131065823396868, 3.1345674541927253, 1.9690673131135776, 6.0866997244685965,
                1.8505751646645416, 5.29869907287348, 3.075020935116759, 4.763084961207747, 0.4168428725881591, 4.563444347669573, 5.8938410276605255, 3.475737135794858, 2.682838532871888, 4.528858505426556, 2.502376089909426, 3.25309267377295, 0.3981977799918744, 1.5490358824309094, 1.605392207533615, 2.9276174450011774, 5.176023132824261, 4.667496500928188,
                1.4236895674455206, 0.4699829750444357, 0.18966103260989964, 6.116881845235636, 2.5456697105062878, 0.8195055516902616, 0.5839032935357767, 1.3320666393453522, 4.00921385818987, 1.2885309542722414, 0.2669323490457898, 1.4846125920567927, 1.3321375588014992, -1.0385205544988825, -1.255718489134564, -0.762581861398507, 0.23349697833423816, -1.4214953587602666, 0.16600711040551805, 0.38079977954047267, 0.08272575553575008, -0.5615305504116471, 1.3659766000268307, -0.7206954265023001,
                -0.34111848138559286, 0.9492802510411192, 0.9696926710212368, 0.9514530768928352, 3.20925079740019, -1.3033896440740826, 0.2767827299976193, -1.1612348919954187, 0.396448190226119, 0.10291336869644606, 0.595677796971947, -0.20902535182069515, 1.0078499189274748, -0.6652825044185814, 0.8377053255180636, -0.3752654927981753, -0.07700315438659211, 1.0125454638264197, 0.14378557962449287, 0.31726441097031294, -0.4329975746974675, -1.0398260290039212, -1.4501778593806574,
                0.8940147876416912, -0.6271218996797788, 0.6388019832355858, 0.34767221087212374, 0.10881231870129565, 0.9238424631683654, 0.028906967986920593, -0.999931925575978, 0.3280366268086224, -0.27572776646487696, -0.9281902110856616, -0.2797903644082264, -0.44592105764345563, -0.16884911353490914, 1.768020028428105, 0.06285533024920076, 0.8978704094466262, 1.1675974166241327, -0.5956955403537542, -1.9491948074566552, -0.0351527475883048, 0.6011609607494409, -0.06847762506277805,
                -0.807702188710335, 0.2751655770483313, 0.2870696415156405, -1.8278998853374089, 0.5708441577491565, -0.7502040617590306, 1.5481226493043265, -0.7532689382609125, -0.8439846032667974, 0.31958743110203486, -0.4972383360844897, 0.06563466828128733, -0.5926193445951025, -0.3425295593954885, 0.10481679368604506, -0.7663427835996953, -0.45117880260394183, -1.0195407087739505, -1.1649080581517623, 0.09855272107031215, -0.6370574781928262,
                -0.30441228195058734,-1.4785611365428615,2.6653809232155594,-0.6906175890897475,0.515481284421803,0.6771994837291646,0.3622428652096128,1.346461113110066,-0.3233109713629789,0.672464115364142,0.2443557547497113,1.42134607333706,-0.34187838136807186,1.537560494740988,0.2864192485202193,-0.9103022225613995,-0.45768058682677387,0.3766506042943305,-1.634708422936025,
                0.011770117214530158,-1.5851861120391288,0.8981172646323772,0.19522845710532188,1.7674908605938122,0.27359799530949114,-0.9430199565078057,-1.4153179394284294,0.6701106931374562,-1.214084221734669,-0.29069920571689345,0.40513186349856933,0.15539511274069823,-0.3577958094790661,-0.11667209188284637,0.36153198836215833,0.3436519297296048,-1.5359275491680546,0.6620070392468738,0.4066790657918409,0.16942397768731038,-0.43046331207511795,
                0.13015712229530746,0.775237635625121,-1.414132285437142,-1.0072264707274916,-0.3858910021039693,1.0954014935056957,-0.7567040019714246,1.6326401183754162,-0.0921318398276079,0.46353691147040454,1.106628164834488,-0.9120287698005056,-0.22062067140142333,-1.4903166907308258,-0.9468628826313928,1.4592520680833163,0.779895415667417,-0.8256911317233132,
                -0.3543526941863613, -0.14591846381835635, 0.864482614222563, 0.9234241594754865, 0.03187074254078362, -0.09793151831367036, -0.20015033788922665, -0.9790777117759897, 0.1051114571305638, 0.06438865461463174, 0.9364045653898226, 0.39765164743902587, -0.2558983497157908, -1.4702702926299065, 0.17316645733119101, -0.8642893069811461, 0.31972250350656883, -0.3509046560560277, -1.0483464458554193, 0.04517806532859413])

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["linear_ent_params"] =
        SA.SVector{14}(rand(14))

    # General stochastic entrainment/detrainment parameters
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["general_stochastic_ent_params"] =
    SA.SVector(
        0.1, 0.1,   # ε_σ², δ_σ²
        0.05, 0.05  # ε_λ, δ_λ
    )

    # For FNO add here
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_width"] = 2
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_n_modes"] = 2
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_params"] =
        SA.SVector{50}(rand(50))

    #calibrated RF
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] =  0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["stochastic_entrainment"] =  "deterministic"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"] =  "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["l_max"] =  1.0e6
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] =  0.033664101288398834
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_n_modes"] =  2
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] =  10.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] =  "buoy_vel"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_width"] =  2
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] =  25.86429119909126
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] =  0.0829262449532631
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["max_area"] =  0.9
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["nn_ent_biases"] =  false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] =  0.05941752121234454
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] =  0.3
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["constant_area"] =  false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_area"] =  1.0e-5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"] =  "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_upd_velocity"] =  0.001
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] =  "RF"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_rm"] =  1.5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_local_micro"] =  true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] =  0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] =  0.0004
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] =  0.242459512532276
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_ub"] =  0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] =  4.623714442670591
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] =  0.34698667231319186
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] =  0.0007648469415345847
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] =  1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"] =  0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"] =  true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] =  0.51
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] =  2.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] =  0.13
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["mixing_length"] =  "sbtd_eq"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_scale"] =  4.076923076923077
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] =  0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_updraft_top"] =  500.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] =  0.8750439088617863
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] =  0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.07341823210211777

    #! format: on

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
    elseif case_name == "DYCOMS_RF02"
        namelist = DYCOMS_RF02(namelist_defaults)
    elseif case_name == "GABLS"
        namelist = GABLS(namelist_defaults)
    elseif case_name == "SP"
        namelist = SP(namelist_defaults)
    elseif case_name == "DryBubble"
        namelist = DryBubble(namelist_defaults)
    elseif case_name == "LES_driven_SCM"
        namelist = LES_driven_SCM(namelist_defaults)
    else
        error("Not a valid case name")
    end

    if write
        write_file(namelist, root)
    end
    return namelist
end
function Soares(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Soares"

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 8 * 3600.0
    namelist["time_stepping"]["dt_min"] = 1.0

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
    namelist["time_stepping"]["dt_min"] = 1.2

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
    namelist["time_stepping"]["dt_min"] = 6.0

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
    namelist["time_stepping"]["dt_min"] = 10.0
    namelist["meta"]["simname"] = "life_cycle_Tan2018"
    namelist["meta"]["casename"] = "life_cycle_Tan2018"

    return namelist
end
function Rico(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "Rico"

    namelist["grid"]["nz"] = 80
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["adapt_dt"] = false
    namelist["time_stepping"]["t_max"] = 86400.0
    #namelist["time_stepping"]["dt_max"] = 5.0
    namelist["time_stepping"]["dt_min"] = 1.5

    namelist["microphysics"]["precipitation_model"] = "clima_1m"
    namelist["microphysics"]["τ_acnv_rai"] = 2500.0
    namelist["microphysics"]["τ_acnv_sno"] = 100.0
    namelist["microphysics"]["q_liq_threshold"] = 0.5e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-6
    namelist["microphysics"]["microph_scaling"] = 1.0
    namelist["microphysics"]["microph_scaling_dep_sub"] = 1.0
    namelist["microphysics"]["microph_scaling_melt"] = 1.0
    namelist["microphysics"]["E_liq_rai"] = 0.8
    namelist["microphysics"]["E_liq_sno"] = 0.1
    namelist["microphysics"]["E_ice_rai"] = 1.0
    namelist["microphysics"]["E_ice_sno"] = 0.1
    namelist["microphysics"]["E_rai_sno"] = 1.0

    namelist["meta"]["simname"] = "Rico"
    namelist["meta"]["casename"] = "Rico"
    return namelist
end
function TRMM_LBA(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "TRMM_LBA"

    namelist["grid"]["nz"] = 80
    namelist["grid"]["dz"] = 200

    namelist["time_stepping"]["adapt_dt"] = true
    namelist["time_stepping"]["t_max"] = 60 * 60 * 6.0
    namelist["time_stepping"]["dt_max"] = 5.0
    namelist["time_stepping"]["dt_min"] = 1.0

    namelist["microphysics"]["precipitation_model"] = "clima_1m" #"cutoff"
    namelist["microphysics"]["τ_acnv_rai"] = 2500.0
    namelist["microphysics"]["τ_acnv_sno"] = 100.0
    namelist["microphysics"]["q_liq_threshold"] = 0.5e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-6
    namelist["microphysics"]["microph_scaling"] = 1.0
    namelist["microphysics"]["microph_scaling_dep_sub"] = 1.0
    namelist["microphysics"]["microph_scaling_melt"] = 1.0
    namelist["microphysics"]["E_liq_rai"] = 0.8
    namelist["microphysics"]["E_liq_sno"] = 0.1
    namelist["microphysics"]["E_ice_rai"] = 1.0
    namelist["microphysics"]["E_ice_sno"] = 0.1
    namelist["microphysics"]["E_rai_sno"] = 1.0

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
    namelist["time_stepping"]["dt_min"] = 2.0

    namelist["meta"]["simname"] = "ARM_SGP"
    namelist["meta"]["casename"] = "ARM_SGP"

    return namelist
end
function GATE_III(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "GATE_III"

    namelist["grid"]["nz"] = 200 # 1700
    namelist["grid"]["dz"] = 85  # 10

    namelist["time_stepping"]["adapt_dt"] = false
    namelist["time_stepping"]["t_max"] = 3600.0 * 24.0
    namelist["time_stepping"]["dt_max"] = 5.0
    namelist["time_stepping"]["dt_min"] = 2.0

    namelist["microphysics"]["precipitation_model"] = "clima_1m" #"cutoff"

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
    namelist["time_stepping"]["dt_min"] = 6.0

    namelist["meta"]["simname"] = "DYCOMS_RF01"
    namelist["meta"]["casename"] = "DYCOMS_RF01"

    return namelist
end
function DYCOMS_RF02(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "DYCOMS_RF02"

    namelist["grid"]["nz"] = 30
    namelist["grid"]["dz"] = 50

    namelist["time_stepping"]["adapt_dt"] = true
    namelist["time_stepping"]["t_max"] = 60 * 60 * 6.0
    namelist["time_stepping"]["dt_max"] = 4.0
    namelist["time_stepping"]["dt_min"] = 1.0

    namelist["microphysics"]["precipitation_model"] = "clima_1m" #"cutoff"
    namelist["microphysics"]["τ_acnv_rai"] = 2500.0
    namelist["microphysics"]["τ_acnv_sno"] = 100.0
    namelist["microphysics"]["q_liq_threshold"] = 0.5e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-6
    namelist["microphysics"]["microph_scaling"] = 1.0
    namelist["microphysics"]["microph_scaling_dep_sub"] = 1.0
    namelist["microphysics"]["microph_scaling_melt"] = 1.0
    namelist["microphysics"]["E_liq_rai"] = 0.8
    namelist["microphysics"]["E_liq_sno"] = 0.1
    namelist["microphysics"]["E_ice_rai"] = 1.0
    namelist["microphysics"]["E_ice_sno"] = 0.1
    namelist["microphysics"]["E_rai_sno"] = 1.0

    namelist["meta"]["simname"] = "DYCOMS_RF02"
    namelist["meta"]["casename"] = "DYCOMS_RF02"

    return namelist
end
function GABLS(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["meta"]["casename"] = "GABLS"

    namelist["grid"]["nz"] = 8
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["t_max"] = 9 * 3600.0
    namelist["time_stepping"]["dt_min"] = 4.0
    namelist["time_stepping"]["dt_max"] = 8.0
    namelist["meta"]["simname"] = "GABLS"
    namelist["meta"]["casename"] = "GABLS"

    return namelist
end
# Not fully implemented yet - Ignacio
function SP(namelist_defaults)

    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "SP"

    # this case is resolution dependent, we should check why
    namelist["grid"]["nz"] = 64
    namelist["grid"]["dz"] = 32

    namelist["time_stepping"]["t_max"] = 7200.0
    namelist["time_stepping"]["dt_min"] = 1.0
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
    namelist["time_stepping"]["dt_min"] = 0.5

    namelist["meta"]["simname"] = "DryBubble"
    namelist["meta"]["casename"] = "DryBubble"

    namelist["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.4

    return namelist
end

function LES_driven_SCM(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    # Only one can be defined by user
    # namelist["grid"]["dz"] = 50.0
    namelist["grid"]["nz"] = 80

    namelist["stats_io"]["frequency"] = 10.0
    namelist["time_stepping"]["t_max"] = 3600.0 * 6
    namelist["time_stepping"]["dt_min"] = 1.0

    # use last 6 hours of LES simulation to drive LES
    namelist["t_interval_from_end_s"] = 3600.0 * 6
    # average in 1 hour interval around `t_interval_from_end_s`
    namelist["initial_condition_averaging_window_s"] = 3600.0

    namelist["meta"]["lesfile"] =
        joinpath(les_driven_scm_data_folder(), "Stats.cfsite23_HadGEM2-A_amip_2004-2008.07.nc")
    namelist["meta"]["simname"] = "LES_driven_SCM"
    namelist["meta"]["casename"] = "LES_driven_SCM"
    namelist["forcing"] = Dict()
    namelist["forcing"]["nudging_timescale"] = 6.0 * 3600.0

    return namelist
end

function write_file(namelist, root::String = ".")
    mkpath(root)

    @assert haskey(namelist, "meta")
    @assert haskey(namelist["meta"], "simname")

    casename = namelist["meta"]["casename"]
    open(joinpath(root, "namelist_$casename.in"), "w") do io
        JSON.print(io, namelist, 4)
    end

    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_namelist(nothing)
end

end
