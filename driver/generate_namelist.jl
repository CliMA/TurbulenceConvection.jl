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
    truncate_stack_trace::Bool = false,
)

    if set_seed
        Random.seed!(seed)
    end

    namelist_defaults = Dict()
    namelist_defaults["meta"] = Dict()
    namelist_defaults["meta"]["uuid"] = basename(tempname())

    namelist_defaults["logging"] = Dict()
    namelist_defaults["logging"]["truncate_stack_trace"] = truncate_stack_trace

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
    namelist_defaults["thermodynamics"]["thermo_covariance_model"] = "diagnostic" #"prognostic" or "diagnostic"
    namelist_defaults["thermodynamics"]["diagnostic_covar_limiter"] = 1e-3 # this controls the magnitude of the spike in covariance
    namelist_defaults["thermodynamics"]["sgs"] = "quadrature"
    namelist_defaults["thermodynamics"]["quadrature_order"] = 3
    namelist_defaults["thermodynamics"]["quadrature_type"] = "log-normal" #"gaussian" or "log-normal"

    namelist_defaults["time_stepping"] = Dict()
    namelist_defaults["time_stepping"]["dt_max"] = 12.0
    namelist_defaults["time_stepping"]["dt_min"] = 1.0
    namelist_defaults["time_stepping"]["adapt_dt"] = true
    namelist_defaults["time_stepping"]["cfl_limit"] = 0.5

    namelist_defaults["microphysics"] = Dict()
    namelist_defaults["microphysics"]["precipitation_model"] = "None"

    namelist_defaults["turbulence"]["scheme"] = "EDMF_PrognosticTKE"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = 1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "FNO"  # {"moisture_deficit", "NN", "NN_nonlocal", "Linear", "FNO", "RF"}
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] = "buoy_vel" # {"buoy_vel", "inv_z", "none"}
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_pi_subset"] = ntuple(i -> i, 6) # or, e.g., (1, 3, 6)
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
        vec(cat(randn(2,100), # vec(cat(randn(2, m),
                    ones(2,7), dims=2)) # ones(2, d + 1), dims=2))

    # RF: fixed realizations of random variables, 2 x m x (1 + d)
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["rf_fix_ent_params"] =
        vec(cat(2*pi*rand(2,100,1), # vec(cat(2*pi*rand(2, m, 1),
                    randn(2,100,6), dims=3)) # randn(2, m, d), dims=3))

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
        SA.SVector{210}([-0.3092221779242847, -0.04155275036902362, -0.00828794655237279, -0.5338731510547093, 0.4337737876532, -1.25890074544991, 0.606077977802638, 0.18995753697331583, 1.4144247823696445, 1.4972682930298096, 0.6239691813143565, -0.16072912349133647, -1.1364649910759521, 0.663839008107084, 0.3846410604617908, -0.13936331886633976, 0.8777236531444704, -0.9006563241598223, -0.346441743076754, 0.0977301179945955, -1.1566731682161442, 0.46743213746280216, -1.7746504317932341, -0.24108917864654167, 0.8585866414049047,
                0.19026342379447853, -0.2920422536072628, -1.2852235434151387, -0.18689101058901092, -0.5523339498285424, -0.16771886329624072, 0.12480272166747983, 0.011825836446342962, -0.5222928868949412, 1.4783121088075504, -0.3320952074745621, 1.2137883997703947, -0.6127232124238309, -0.22638845702364926, -0.07148033572755433, 1.1180071317347073, 0.38432796906286343, -0.05367941157422737, -0.10640518340683039, 0.33548181588546316, 0.27410330969155816, -0.9380214659775208, 1.6621125150711187, -0.5793384406149311, -0.06283234681502381,
                -0.8633092482481071, -0.06369229975170931, 0.22051167499026436, 1.5928941290515206, -1.5127270447243122, 1.2634125542652817, 1.003726075762215, 0.3014207849001762, 0.49668476449241905, -0.04269048192496615, 1.9305036251199659, -0.27591920711806667, 0.5663875774894052, -0.2751961428620821, -1.7417602324043233, 0.18765319931893748, -1.1204494312597422, 0.4327657366916974, -1.826807978407353, -0.9703377185346372, 0.8791549293364562, -0.6818009900664804, -0.5323801335637169, -1.0628235854219685, 1.056069167120661, 0.3199868930514552, 0.9264098071252714, -0.1170885893354422, -0.35555398345411227,
                0.7428221341430343, 0.47380646463392967, 0.8737401427720888, 0.8642616339776521, 0.12484111900406361, -0.7440970976612922, -0.6993859013579365, 0.6638363238849514, -0.22478469350986266, -0.30710915867836525, -0.5892339510368617, 1.3146488830780119, -0.7236407501729057, 0.21211606924164236, -1.3098999919751972, -1.2149636728508648, 0.08503692928233067, -1.7475133994619434, -0.252856314445811, -0.14588454851153937, -1.2328220521402051, -0.40567417293016633, -0.4764081570622447, 1.3304532799153597, 1.26072133851544, -1.3038353746581601, -0.16018479659455254,
                0.26247783123906226, 1.426935464572116, -1.6389089786160433, 0.12908329549810757, -0.33636081712366817, 0.08206464846297183, 0.9691192028096582, -1.603939392440497, -1.024227878616787, 0.5491589735136971, -0.01099755575516874, 0.38939349022747133, 0.055591934750767916, 1.0393956906316726, 0.08097920965866261, 0.5598944983496505, -1.2524345624370565, -1.1112281928128391, -0.0027949341851488344, -0.5086639442793711, 1.718189374561086, 0.28358787066853963, -1.5573511294539306, 0.08517045779650749, -0.5503375453034913, -1.5253436739305422, 2.198804116518643, -1.402499956550675,
                -0.43064223762230797, 0.15605989490809563, 1.1942548418998111, -0.9418318983986724, -2.5527036231748204, 1.6405329305312868, 0.7441625854053707, -0.520195383437397, -0.6711712819250466, -1.4444998569935863, 0.08166490684799627, 1.549149526734603, -0.2863489349364617, 0.1765185325093699, -0.2233995112783648, -0.9867625865834102, -1.2011668705699439, 0.6823992293448087, 1.8879937393684079, 0.46791872532279194, -0.9608780017604908, 0.05891384238684695, 0.25152529786426925, -0.244629393411706, -0.748963128734426, 1.3698824534234966, 1.1001135565215783, 0.6328809527806487,
                -0.8970421978234012, 0.8191172437372284, 0.301915694921396, 0.35690486203724225, 0.7428814693382113, 0.6874090096071277, -0.733828144741926, -0.1147827810192663, 0.6051601084466574, -1.1202310925370933, -0.6958621452285669, 0.012271891058833626, 0.5309492044847048, 0.035929355631467, 0.07627754231380779, -1.043657412103859, 1.3915059041863, 0.5866310563780519, 0.3209227641070763, -0.6088395158519854, 0.11372833072050852, 0.03635598368335872, 0.291510961029302, -3.062681168290446, 1.0017041578102255, -1.2274463144423111, 0.3629740012826612, 0.7945701039625718,
                0.45047497376696866, -0.3022569521413312, 1.5389061881483015, 0.48668556769276644, -0.3097869203771747, -0.26491273121242376, 0.506330978715088, -0.07841516108517377, -1.1938141746947968, 1.0229268913212626, -0.6903850181394218, -0.5159400548859849, 1.740158466654292, 1.781087094344552, 0.8738212415894402, -0.1970975451710268, 0.38250155525597274, 1.1257496154703959, 1.0975376626020097, 0.6563731316253598
            ])

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["l_max"] = 1.0e6
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.10220377250699608
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_n_modes"] = 4
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 10.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_dim_scale"] = "inv_z"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["fno_ent_width"] = 4
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] = 15.781660336189047
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.13800888724279345
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = 0.9
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["nn_ent_biases"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.1062537310824627
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = 1.0e-5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_upd_velocity"] = 0.001
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "FNO"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_rm"] = 1.5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_local_micro"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.0004
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.14115523161055626
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_ub"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] = 1.982975035816953
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.5380391199229009
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.0011509841897409344
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entr_pi_subset"] = [ 1, 2, 3, 4, 5, 6]
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["nn_arc"] = [ 6, 5, 4, 2]
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_scale"] = 4.076923076923077
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pi_norm_consts"] = [ 478.298, 1.0, 1.0, 1.0, 1.0, 1.0]
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] = 0.8616622353153358
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.0989308955867611
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

    namelist["microphysics"]["precipitation_model"] = "clima_1m" # "cutoff"
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

    # LES filename should follow pattern:
    # Stats.cfsite<SITE-NUMBER>_<FORCING-MODEL>_<EXPERIMENT>_2004-2008.<MONTH>.nc
    namelist["meta"]["lesfile"] =
        joinpath(les_driven_scm_data_folder(), "Stats.cfsite23_HadGEM2-A_amip_2004-2008.07.nc")
    namelist["meta"]["simname"] = "LES_driven_SCM"
    namelist["meta"]["casename"] = "LES_driven_SCM"

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
