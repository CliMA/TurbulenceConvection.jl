"""
Get timescales for relaxation of the supersaturation


Throughout, mean_r_factor exists so you can choose dfiferent estimates of <r> based on the q and N you're using.
For example in gamma distribution, <r> = 1/2λ, q = N_0 Γ(5) / λ^4 , and N = N_0/λ so <r> = (32q/N)^1/3. Then to go from r = (q/(4/3πρN))^(1/3) to  (32q/N)^1/3 we need to multiply by (24/(πρ))^(1/3)
"""

# N can be arbitrarily small so D * N_l could yield 0, times r = Inf could yield NaN... always group (N * r) to avoid this, can't have Inf and 0 here bc of condition in r_from_qN()

# ======================================================================================================================================== #
# D_func(T::FT, p::FT) where {FT} = FT((2.11 * 1e-5) * (T / FT(273.15))^1.94 * (p / 101325))::FT  # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
D_func(T::FT, p::FT) where {FT} = FT((2.11 * 1e-5) * (T / FT(273.15))^1.94 * (FT(101325)/p))::FT  # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226) [e.g. 12.1 in https://apps.dtic.mil/sti/tr/pdf/ADA440352.pdf]
# ======================================================================================================================================== #

function get_τ(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::AbstractRelaxationTimescaleType, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    τ_liq, τ_ice =  get_τ_helper(param_set, microphys_params, relaxation_timescale_type, q, T, p, ρ, w, z)

    # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster... change to be dt related) [used to have FT(1)]
    τ_liq = clamp(τ_liq, relaxation_timescale_type.args.min_τ_liq, relaxation_timescale_type.args.max_τ_liq)
    τ_ice = clamp(τ_ice, relaxation_timescale_type.args.min_τ_ice, relaxation_timescale_type.args.max_τ_ice)

    if isnan(τ_liq) || isnan(τ_ice)  # chasing down NaNs
        error("Timescale calculation failed. Got τ_liq: $τ_liq, τ_ice: $τ_ice for inputs q = $q; T = $T; p = $p; ρ = $ρ; w = $w; z = $z; relaxation_timescale_type = $relaxation_timescale_type")
    end

    return τ_liq, τ_ice
end

function get_τ(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::Symbol, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    return get_τ(param_set, microphys_params, get_relaxation_timescale_type(relaxation_timescale_type, param_set), q, T, p, ρ, w, z)
end

# ---------------------------------------------------------------------------------------------------------------------------------------- #

# :Base
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    return relaxation_timescale_type.τ_liq, relaxation_timescale_type.τ_ice
end

# :exponential_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::ExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    # T_fr = TCP.T_freeze(param_set)
    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius

    r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
    τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))
    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
    return τ_liq, τ_ice
end

# :exponential_T_scaling_ice_raw
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::ExponentialTScalingIceRawRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)
    q_0 = FT(1e-6)
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    τ_ice = FT(1 / (c_1 * exp(c_2 * (T - T_fr))))
    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
    return τ_liq, τ_ice
end

# :powerlaw_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::PowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)
    if T >= T_fr
        τ_ice = FT(Inf)
    else
        ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
        r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
        N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
        τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))
    end
    τ_liq = CMNe.τ_relax(microphys_params, liq_type)
    return τ_liq, τ_ice
end

# :exponential_times_powerlaw_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::ExponentialTimesPowerlawScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    error("NotImplmentedError: This relaxation_timescale_type functionality has not been implemented yet")
end

# :geometric_liq__geometric_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::GeometricLiqGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius

    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_l = r_from_qN(q.liq, N_l, ρ_l, r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_liq)
    τ_liq = FT(1 / (4 * π * D * (N_l * r_l)))

    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
    τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))
    return τ_liq, τ_ice
end

# :geometric_liq__exponential_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::GeometricLiqExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius

    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_liq)
    τ_liq = FT(1 / (4 * π * D * (N_l * r_l)))

    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
    τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))
    return τ_liq, τ_ice
end

# :geometric_liq__powerlaw_T_scaling_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::GeometricLiqPowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius

    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_liq)
    τ_liq = FT(1 / (4 * π * D * (N_l * r_l)))

    if T >= T_fr
        τ_ice = FT(Inf)
    else
        N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
        r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
        τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))
    end
    return τ_liq, τ_ice
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)

    ρ_l = CMP.ρ_cloud_liq(microphys_params) # CLIMAParameters default for cloud_liquid   
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice

    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius

    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_liq)
    τ_liq = FT(1 / (4 * π * D * (N_l * r_l)))

    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
    τ_ice = 1 / (4 * π * D * (N_i * r_i))

    return τ_liq, τ_ice
end

# :linear_combination
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::LinearCombinationRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)

    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius

    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_liq)
    τ_liq = FT(1 / (4 * π * D * (N_l * r_l)))

    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
    τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))

    return τ_liq, τ_ice
end

# :linear_combination_with_w
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::LinearCombinationWithWRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    T_fr = TCP.T_freeze(param_set)
    ρ_l = CMP.ρ_cloud_liq(microphys_params) # kg m^-3, CLIMAParameters default for cloud_liquid
    ρ_i = CMP.ρ_cloud_ice(microphys_params) # CLIMAParameters default for cloud_ice
    r_0 = param_set.user_params.particle_min_radius # CLIMAParameters default for particle_min_radius
    #
    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_l = r_from_qN(q.liq, N_l, ρ_l; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_liq)
    τ_liq = FT(1 / (4 * π * D * (N_l * r_l))) 
    #

    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
    r_i = r_from_qN(q.ice, N_i, ρ_i; r_min = r_0, mean_r_factor = param_set.user_params.mean_r_factor_ice)
    τ_ice = FT(1 / (4 * π * D * (N_i * r_i)))
    return τ_liq, τ_ice
end

# :neural_network
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    # neural_network_params = param_set.user_params.neural_microphysics_relaxation_network # training weights for tau (made so could pass in things in TrainTau.jl but now I use CalibrateEDMF)
    # neural_network_params = Float32.(collect(neural_network_params)) # collect from ntuple, then convert to Float32 for NN since that's what it's supposed to be (save eltype in the jld2?)
    # model_x_0_characteristic = TCP.get_isbits_nt(param_set.user_params, :model_x_0_characteristic, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)
    model_x_0_characteristic = relaxation_timescale_type.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale_type
    # TC = TurbulenceConvection
    # if isdefined(TC, :neural_network) # do this so we don't have to read from disk and recreate the network evry time bc that's super slow
        # TC.neural_network = TC.re(neural_network_params) # set the parameters to the ones we just read in                       
    # else
        # model_re_location = TCP.get_isbits_nt(param_set.user_params, :model_re_location, nothing) # store the model structure so you can create the model from the parameters in the namelist (the parameters must be in namelist for CalibrateEDMF.jl)
        # we could also red this direct from relaxation_timescale_type right?
        # model_re_location = string(TCP.unwrap_val(relaxation_timescale_type.model_re_location)) # get the model_re_location from the relaxation_timescale_type
        # if !isnothing(model_re_location)
            # model_re_location = model_re_location
        # end

        # re = model_destructure_re_from_file(model_re_location) # get the reconstruction function from the file ( this is probably slow, we could pass the string repr)
        # neural_network = vec_to_NN(neural_network_params, re) # construct the NN from the parameters
        # TC.re = re # store it in the TC so we don't have to reconstruct it every time
        # TC.neural_network = neural_network # store it in the TC so we don't have to reconstruct it every time
    # end

    τ_liq, τ_ice, _, _ = predict_τ(ρ, T, q, w, neural_network, model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?

    return τ_liq, τ_ice
end

# :raymond_ice_test
function get_τ_helper(param_set::APS, microphys_params::ACMP, relaxation_timescale_type::RaymondIceTestRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, ρ::FT, w::FT, z::FT,) where {FT}
    T_fr = TCP.T_freeze(param_set)

    if (N_0 = relaxation_timescale_type.N_0) isa NamedTuple
        N_0 = pyinterp(z, N_0.z, N_0.values)
    end

    N_l = get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
    N_i = get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)

    if relaxation_timescale_type.N_r_closure isa InhomogeneousNRClosure
        thermo_params = TCP.thermodynamics_params(param_set)
        ts_LCL = error("ts_LCL not implemented yet...")
        N_l, R_liq = NR_inhomogeneous_mixing_liquid(thermo_params, N_l, p, q.liq, ts_LCL) # testing inhomogeneous mixing would have r fixed and then let N vary... set r based on adiabatic extrapolation from cloud base 
        _, R_ice = NR_monodisperse(N_i, q.ice)
        # maybe look into using a mixed model where N/R are partly towards the inhomogenous value depending on the true entrainment/mixing params... see literature on this
        # NOTE, ON DYCOMS ADIABATIC R W/ ORIGINAL N_0 WORKED BETTER... HMMM (though that's not using DYCOMS N)

    elseif relaxation_timescale_type.N_r_closure isa MonodisperseNRClosure # uniform size for all droplets, liquid and Ice I guess
        _, R_liq = NR_monodisperse(N_l, q.liq)
        _, R_ice = NR_monodisperse(N_i, q.ice)
    else
        error("Unsupported size distribution closure (N_r_closure): $(relaxation_timescale_type.N_r_closure)")
    end

    D = D_func(T, p) # m2 s**-1 for T in Kelvin, p in Pa (from Pruppacher and Klett 1997) # D_ref = FT(0.0000226)
    base = FT(1 / (4 * π * D)) # as q goes up, R goes 
    τ_liq = base / (N_l * R_liq)
    τ_ice = base / (N_i * R_ice)
    # limit effective tau to at least one second for stability (stable with dt == limit timescale is unstable, gotta go faster...
    return τ_liq, τ_ice
end

# ======================================================================================================================================== #




# ======================================================================================================================================== #

function get_N_i(param_set::APS, relaxation_timescale_type::AbstractRelaxationTimescaleType, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    N_i = get_N_i_helper(param_set, relaxation_timescale_type, q, T, ρ, w)
    if get_adjust_ice_N(relaxation_timescale_type)
        microphys_params::ACMP = TCP.microphysics_params(param_set)
        N_i = adjust_ice_N(param_set, microphys_params, N_i, q; r_min =  param_set.user_params.particle_min_radius)
    end
    N_i = clamp(N_i, relaxation_timescale_type.args.min_N_ice, relaxation_timescale_type.args.max_N_ice)
    return N_i
end

function get_N_l(param_set::APS, relaxation_timescale_type::AbstractRelaxationTimescaleType, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    N_l =  get_N_l_helper(param_set, relaxation_timescale_type, q, T, ρ, w)
    N_l = clamp(N_l, relaxation_timescale_type.args.min_N_liq, relaxation_timescale_type.args.max_N_liq)
    return N_l
end

# --------------------------------------- #

function get_N_i(param_set::APS, relaxation_timescale_type::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT} # should not use this if possible bc it constructs another relaxation tiemscale object...
    return get_N_i(param_set, get_relaxation_timescale_type(relaxation_timescale_type, param_set), ts, w)
end
function get_N_i(param_set::APS, relaxation_timescale_type::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_i(param_set, relaxation_timescale_type, q, T, ρ, w)
end

function get_N_l(param_set::APS, relaxation_timescale_type::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT}
    return get_N_l(param_set, get_relaxation_timescale_type(relaxation_timescale_type, param_set), ts, w)
end
function get_N_l(param_set::APS, relaxation_timescale_type::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT) where {FT}
    thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts)
    T::FT = TD.air_temperature(thermo_params, ts)
    # p::FT = TD.air_pressure(thermo_params, ts) 
    ρ::FT = TD.air_density(thermo_params, ts)
    return get_N_l(param_set, relaxation_timescale_type, q, T, ρ, w)
end

get_N_i(param_set::APS, Ni::FT, ts::TD.ThermodynamicState, w::FT) where {FT} = Ni
get_N_l(param_set::APS, Nl::FT, ts::TD.ThermodynamicState, w::FT) where {FT} = Nl
get_N_i(param_set::APS, Ni::FT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = Ni
get_N_l(param_set::APS, Nl::FT, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = Nl

# ---------------------------------------------------------------------------------------------------------------------------------------- #

# :Base
get_N_i_helper(param_set::APS, relaxation_timescale_type::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)    
get_N_l_helper(param_set::APS, relaxation_timescale_type::BaseRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# :exponential_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale_type::ExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)

    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    return FT(c_1 * exp(c_2 * (T - T_fr)))
end
get_N_l_helper(param_set::APS, relaxation_timescale_type::ExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# :exponential_T_scaling_ice_raw
get_N_i_helper(param_set::APS, relaxation_timescale_type::ExponentialTScalingIceRawRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)
get_N_l_helper(param_set::APS, relaxation_timescale_type::ExponentialTScalingIceRawRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# :powerlaw_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale_type::PowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)

    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i

    if T >= T_fr
        N_i = FT(0) # cant do powerlaw otherwise...
    else
        N_i = (10^c_1) * (-(T - T_fr))^c_2
    end

    return N_i
end
get_N_l_helper(param_set::APS, relaxation_timescale_type::PowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT} =  FT(NaN)

# :exponential_times_powerlaw_scaling_ice
get_N_i_helper(param_set::APS, relaxation_timescale_type::ExponentialTimesPowerlawScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = error("NotImplmentedError: This relaxation_timescale_type functionality has not been implemented yet")
get_N_l_helper(param_set::APS, relaxation_timescale_type::ExponentialTimesPowerlawScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = error("NotImplmentedError: This relaxation_timescale_type functionality has not been implemented yet")


# :geometric_liq__geometric_ice
function get_N_i_helper(param_set::APS, relaxation_timescale_type::GeometricLiqGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    c_3 = relaxation_timescale_type.c_3i
    c_4 = relaxation_timescale_type.c_4i
    return min( FT(c_1 * (q.ice/1e-7)^(c_2) + 10^c_3), FT(10^c_3 + 10^c_4))
end
function get_N_l_helper(param_set::APS, relaxation_timescale_type::GeometricLiqGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale_type.c_1l
    c_2 = relaxation_timescale_type.c_2l
    c_3 = relaxation_timescale_type.c_3l
    c_4 = relaxation_timescale_type.c_4l
    return min(FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3), FT(10^c_3 + 10^c_4))
end

# :geometric_liq__exponential_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale_type::GeometricLiqExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    return FT(c_1 * exp(c_2 * (T - T_fr)))
end
function get_N_l_helper(param_set::APS, relaxation_timescale_type::GeometricLiqExponentialTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale_type.c_1l
    c_2 = relaxation_timescale_type.c_2l
    c_3 = relaxation_timescale_type.c_3l
    c_4 = relaxation_timescale_type.c_4l
    return min( FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3) , FT(10^c_3 + 10^c_4))
end

# :geometric_liq__powerlaw_T_scaling_ice
function get_N_i_helper(param_set::APS, relaxation_timescale_type::GeometricLiqPowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    if T >= T_fr
        N_i = FT(0) # cant do powerlaw otherwise...
    else
        N_i = (10^c_1) * (-(T - T_fr))^c_2
    end
    return N_i
end
function get_N_l_helper(param_set::APS, relaxation_timescale_type::GeometricLiqPowerlawTScalingIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale_type.c_1l
    c_2 = relaxation_timescale_type.c_2l
    c_3 = relaxation_timescale_type.c_3l
    c_4 = relaxation_timescale_type.c_4l
    return min(FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3), FT(10^c_3 + 10^c_4))
end

# :geometric_liq__exponential_T_scaling_and_geometric_ice
# function get_N_i_helper(param_set::APS, relaxation_timescale_type::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, p::FT, w::FT) where {FT}
#     T_fr = TCP.T_freeze(param_set)
#     c_1 = relaxation_timescale_type.c_1i
#     c_2 = relaxation_timescale_type.c_2i
#     c_3 = relaxation_timescale_type.c_3i
#     c_4 = relaxation_timescale_type.c_4i
#     return (c_1 * (q.ice/1e-7)^(c_2) + c_3) * exp(c_4 * (T - T_fr))
# end
function get_N_i_helper(param_set::APS, relaxation_timescale_type::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    c_3 = relaxation_timescale_type.c_3i
    c_4 = relaxation_timescale_type.c_4i
    c_5 = relaxation_timescale_type.c_5i
    N_i_T = c_4 * exp(c_5 * (T - T_fr))
    N_i_q = c_1 * (q.ice/1e-7)^(c_2) + (10^c_3) * N_i_T
    return min(N_i_q, N_i_T)
end

function get_N_l_helper(param_set::APS, relaxation_timescale_type::GeometricLiqExponentialTScalingAndGeometricIceRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    c_1 = relaxation_timescale_type.c_1l
    c_2 = relaxation_timescale_type.c_2l
    c_3 = relaxation_timescale_type.c_3l
    c_4 = relaxation_timescale_type.c_4l
    return min( FT(c_1 * (q.liq/1e-5)^(c_2) + 10^c_3) , FT(10^c_3 + 10^c_4))
end

# :linear_combination
function get_N_i_helper(param_set::APS, relaxation_timescale_type::LinearCombinationRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    c_3 = relaxation_timescale_type.c_3i
    return exp(c_1 + c_2 * (T - T_fr) + c_3 * q.ice / 1e-7)
end
function get_N_l_helper(param_set::APS, relaxation_timescale_type::LinearCombinationRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    c_1 = relaxation_timescale_type.c_1l
    c_2 = relaxation_timescale_type.c_2l
    c_3 = relaxation_timescale_type.c_3l
    N_l = exp(c_1 + c_2 * (T - T_fr) + c_3 * q.liq )
end

# :linear_combination_with_w
function get_N_i_helper(param_set::APS, relaxation_timescale_type::LinearCombinationWithWRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    w_0 = FT(1e-3) # 1 mm/s
    c_1 = relaxation_timescale_type.c_1i
    c_2 = relaxation_timescale_type.c_2i
    c_3 = relaxation_timescale_type.c_3i
    c_4 = relaxation_timescale_type.c_4i
    N_i = exp(c_1 + c_2 * (T - T_fr) + c_3 * (q.ice / 1e-7) + c_4 * w / w_0)
end
function get_N_l_helper(param_set::APS, relaxation_timescale_type::LinearCombinationWithWRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    T_fr = TCP.T_freeze(param_set)
    w_0 = FT(1e-3) # 1 mm/s
    c_1 = relaxation_timescale_type.c_1l
    c_2 = relaxation_timescale_type.c_2l
    c_3 = relaxation_timescale_type.c_3l
    c_4 = relaxation_timescale_type.c_4l
    return exp(c_1 + c_2 * (T - T_fr) + c_3 * q.liq + c_4 * w / w_0)
end

# :neural_network
function get_N_i_helper(param_set::APS, relaxation_timescale_type::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    model_x_0_characteristic = relaxation_timescale_type.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale_type
    _, _, _, N_i = predict_τ(ρ, T, q, w, neural_network, model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?
    return N_i
end

function get_N_l_helper(param_set::APS, relaxation_timescale_type::NeuralNetworkRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT}
    model_x_0_characteristic = relaxation_timescale_type.model_x_0_characteristic # get the model_x_0_characteristic from the relaxation_timescale_type
    _, _, N_l, _ = predict_τ(ρ, T, q, w, neural_network, model_x_0_characteristic) # pass in the NN and get the τs out, neural_network should be a global variable?
    return N_l
end

# :raymond_ice_test
get_N_i_helper(param_set::APS, relaxation_timescale_type::RaymondIceTestRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)
get_N_l_helper(param_set::APS, relaxation_timescale_type::RaymondIceTestRelaxationTimescale, q::TD.PhasePartition, T::FT, ρ::FT, w::FT) where {FT} = FT(NaN)

# ======================================================================================================================================== #

function get_Ns(param_set::APS, relaxation_timescale_type::Symbol, ts::TD.ThermodynamicState, w::FT) where {FT}
    return get_Ns(param_set, get_relaxation_timescale_type(relaxation_timescale_type, param_set), ts, w)
end

function get_Ns(param_set::APS, relaxation_timescale_type::AbstractRelaxationTimescaleType, ts::TD.ThermodynamicState, w::FT) where {FT}
    N_l::FT = get_N_l(param_set, relaxation_timescale_type, ts, w)
    N_i::FT = get_N_i(param_set, relaxation_timescale_type, ts, w)
    return (; liq = N_l, ice = N_i)
end

# ======================================================================================================================================== #
