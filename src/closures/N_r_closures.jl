
# maybe add an oversider dispatch fcn here?

function NR_monodisperse(N::FT, q::FT; r_0::FT = 0.2 * 1e-6, ρ_w = 1000) where {FT}
    """
    constant N
    """
    r = (q / ((4 / 3) * π * ρ_w * N))^(1 / 3)
    r = max(r, 0.2 * 1e-6) # bound to be at least ~micron size...something like kohler crit radius (e.g. if the cloud is new, diagnosed r is 0 for selecting tau in principle)
    return (; N, r)
end


function NR_fixed_radius(r::FT, q::FT; ρ_w = 1000) where {FT}
    """
    constant r
    """
    N = q / ((4 / 3) * π * r^3 * ρ_w)
    return (; N, r)
end




function NR_inhomogeneous_mixing_liquid(
    param_set::TD.APS,
    N::FT,
    p::FT,
    q_liq::FT,
    p_LCL::FT,
    T_LCL::FT,
    q_LCL::FT,
) where {FT}
    """
    # add some DOI references

    # Diagnose q_l from adiabatic parcel ascent from LCL to this altitude
    # Calculate r from this q_l and our assumed N
    # Calculate real n given our true q_l and this r -- implication is we're subadiabatic and n is depleted by entrainment.

    In a mixed-phase cloud then, WBF would deplete q_l even further leading to reduced n -- an interesting idea considering the Tan et al style pocketing idea of mixed phase clouds, maybe this is better than reducing r homogeneously.

    We have no similar method for ice because we have no expectation of an N -- INP are continually activated as a function of temperature.
    The ice N/R relationship maybe is best left as being uniform  since earlier particles would have grown more but there's less of them (don't know relative scaling)... or some other data informed parameterization...
    A potential avenue for thinking about equilibrium adiabat q_l, q_i is Korolev, Field 2006 https://doi.org/10.1175/2007JAS2355.1, but then it would be horrible to reconvert that back to some representation of anything for ice...

    What we'd really want is a lagrangian parcel model that evolves the ice size distribution that we can compare to -- couldn't find obs of ice size dist as fcn of height...



    what happens if the previous LCL is above where we are now? is the calculation order correct to resolve this? (then the cloud base might be at the toa... not ideal...)

    """
    # theta_li is conserved so we just need to back out states based on theta_li and our two pressures
    if p_LCL < p # LCL is above current condensation, so new cloud... use default monodisperse relation
        return NR_monodisperse(N, q_liq)
    else
        q_LCL = TD.PhasePartition(q_LCL) # assume all vapor at cloud base (or should I read in the real value?) would have to make q_pt an argument rather than q_LCL or only accept ts for example...
        θ_liq_ice = TD.liquid_ice_pottemp_given_pressure(param_set, T_LCL, p_LCL, q_LCL)
        # @show(p, θ_liq_ice)
        q_adiabatic = PhaseEquil_pθq_q(param_set, p, θ_liq_ice, q_LCL.tot)# can we ascend the adiabat here without having to do a full solve? -- also how do I partition liquid and ice in this framework? assume one phase?e
        ρ_w = FT(1)
        r_adiabatic = (q_adiabatic.liq / (4 / 3 * π * ρ_w * N))^(1 / 3)
        r_adiabatic = max(r_adiabatic, 0.2 * 1e-6) # bound to be at least ~micron size...something like kohler crit radius (e.g. if the cloud is new, diagnosed r is 0 for selecting tau in principle) (also if q_liq is 0 so is r_liq so rly need some bounds)
        N = q_liq / (4 / 3 * π * r_adiabatic^3 * ρ_w)
        return (; N, r = r_adiabatic)
    end
end

function NR_inhomogeneous_mixing_liquid(
    param_set::TD.APS,
    N::FT,
    p::FT,
    q_liq::FT,
    ts_LCL::TD.ThermodynamicState,
) where {FT}
    p_LCL = TD.air_pressure(param_set, ts_LCL)
    T_LCL = TD.air_temperature(param_set, ts_LCL)
    q_LCL = TD.total_specific_humidity(param_set, ts_LCL)
    return NR_inhomogeneous_mixing_liquid(param_set, N, p, q_liq, p_LCL, T_LCL, q_LCL)
end

function NR_inhomogeneous_mixing_liquid(
    param_set::TD.APS,
    N::FT,
    ts::TD.ThermodynamicState,
    q_liq::FT,
    ts_LCL::TD.ThermodynamicState,
) where {FT}
    p = TD.air_pressure(param_set, ts)
    p_LCL = TD.air_pressure(param_set, ts_LCL)
    T_LCL = TD.air_temperature(param_set, ts_LCL)
    q_LCL = TD.total_specific_humidity(param_set, ts_LCL)
    return NR_inhomogeneous_mixing_liquid(param_set, N, p, q_liq, p_LCL, T_LCL, q_LCL)
end


"""
Phase equil but the output is always liquid, pretend ice doesnt exist even below freezing...
"""
function PhaseEquil_pθq_q(
    # param_set::APS,
    param_set::TD.APS, # a thermo_params to pass in
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    maxiter::IT = nothing,
    relative_temperature_tol::FTT = nothing,
    ::Type{sat_adjust_method} = TD.RS.SecantMethod,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, IT <: TD.ITERTYPE, FTT <: TD.TOLTYPE(FT), sat_adjust_method}
    maxiter === nothing && (maxiter = 50)
    relative_temperature_tol === nothing && (relative_temperature_tol = FT(1e-4))
    phase_type = TD.PhaseEquil{FT} # equil or non equil
    q_tot_safe = clamp(q_tot, FT(0), FT(1))
    T = saturation_adjustment_given_pθq( # use my local version... force liq_frac to 1
        sat_adjust_method,
        param_set,
        p,
        θ_liq_ice,
        q_tot_safe,
        phase_type,
        maxiter,
        relative_temperature_tol,
        T_guess,
    ) # finds T based on liq/ice potential temperature and q_tot ... afterwards it's up to you to back out the phase partitioning I guess
    q_pt = PhasePartition_equil_given_p(param_set, T, p, q_tot_safe, phase_type) # our version of this uses liq_frac always 1, and in principle I guess T from saturation adjustment didn't care about what the phase partition looked like? (that doesnt feel right but...) edit -- replaced with my own version whew
    return q_pt
    # ρ = TD.air_density(param_set, T, p, q_pt)
    # e_int = TD.internal_energy(param_set, T, q_pt)
    # return TD.PhaseEquil{FT}(ρ, p, e_int, q_tot_safe, T) # for our output we want the 
end


function PhasePartition_equil_given_p(
    param_set::TD.APS,
    T::FT,
    p::FT,
    q_tot::FT,
    ::Type{phase_type},
) where {FT <: Real, phase_type <: TD.ThermodynamicState}

    q_v_sat = TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
    # _liquid_frac = liquid_fraction(param_set, T, phase_type;ramp=false)
    _liquid_frac = 1 #force liquid for our adiabtic parcel test...
    q_c = q_tot - q_v_sat
    q_liq = _liquid_frac * q_c
    q_ice = (1 - _liquid_frac) * q_c
    return TD.PhasePartition(q_tot, q_liq, q_ice)
end


# adjust this one to make sure it uses our local PhasePartition_equil_given_p which always has liq_frac = 1
function saturation_adjustment_given_pθq(
    ::Type{sat_adjust_method},
    param_set::TD.APS,
    p::FT,
    θ_liq_ice::FT,
    q_tot::FT,
    ::Type{phase_type},
    maxiter::Int,
    relative_temperature_tol::FT,
    T_guess::Union{FT, Nothing} = nothing,
) where {FT <: Real, sat_adjust_method, phase_type <: TD.PhaseEquil}
    tol = TD.RS.RelativeSolutionTolerance(relative_temperature_tol)
    T_min::FT = TD.TP.T_min(param_set)
    T_freeze::FT = TD.TP.T_freeze(param_set)
    cp_d::FT = TD.TP.cp_d(param_set)
    cp_v::FT = TD.TP.cp_v(param_set)
    air_temp(q) = TD.air_temperature_given_pθq(param_set, p, θ_liq_ice, q)
    function θ_liq_ice_closure(T)
        # use local version of PhasePartition_equil_given_p which always has liq_frac = 1
        q_pt = PhasePartition_equil_given_p(param_set, T, oftype(T, p), oftype(T, q_tot), phase_type)
        return TD.liquid_ice_pottemp_given_pressure(param_set, T, oftype(T, p), q_pt)
    end
    q_vap_sat(T) = TD.q_vap_saturation_from_pressure(param_set, q_tot, p, T, phase_type)
    T_1 = max(T_min, air_temp(TD.PhasePartition(q_tot))) # Assume all vapor
    q_v_sat_1 = q_vap_sat(T_1)
    unsaturated = q_tot <= q_v_sat_1
    if unsaturated && T_1 > T_min
        return T_1
    end
    temperature_tol = T_freeze * relative_temperature_tol
    θ_liq_ice_upper = θ_liq_ice_closure(T_freeze + temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    θ_liq_ice_lower = θ_liq_ice_closure(T_freeze - temperature_tol / 2) # /2 => resulting interval is `temperature_tol` wide
    if θ_liq_ice_lower < θ_liq_ice < θ_liq_ice_upper
        return T_freeze
    end
    roots(T) = oftype(T, θ_liq_ice) - θ_liq_ice_closure(T)
    sol = TD.RS.find_zero(
        roots,
        TD.sa_numerical_method_pθq(sat_adjust_method, param_set, p, θ_liq_ice, q_tot, phase_type, T_guess),
        TD.RS.CompactSolution(),
        tol,
        maxiter,
    )
    if !sol.converged
        if TD.print_warning()
            TD.KA.@print("-----------------------------------------\n")
            TD.KA.@print("maxiter reached in saturation_adjustment_given_pθq:\n")
            TD.print_numerical_method(sat_adjust_method)
            TD.print_T_guess(sat_adjust_method, T_guess)
            TD.KA.@print(", p=", p)
            TD.KA.@print(", θ_liq_ice=", θ_liq_ice)
            TD.KA.@print(", q_tot=", q_tot)
            TD.KA.@print(", T=", sol.root)
            TD.KA.@print(", maxiter=", maxiter)
            TD.KA.@print(", tol=", tol.tol, "\n")
        end
        if TD.error_on_non_convergence()
            error("Failed to converge with printed set of inputs.")
        end
    end
    return sol.root
end


function cloud_base(aux, grid, ts, mode)
    """
    adopted from driver/compute_diagnostics.jl

    really just need to map from from cloud base p,T,q to local p,T,q

    either runs in each updraft or in the environment so doesn't need to worry about that structure
    """


    k_end = collect(real_center_indices(grid))[end]

    if mode === :env
        cloud_base = zc_toa(grid).z
        cloud_base_ts = ts[k_end] # i think it goes bottom to top but need to check...
        @inbounds for k in real_center_indices(grid)
            if aux.area[k] > 1e-3
                if TD.has_condensate(aux.q_liq[k] + aux.q_ice[k])
                    if cloud_base > grid.zc[k].z
                        cloud_base = grid.zc[k].z
                        cloud_base_ts = ts[k]
                    end
                end
            end
        end
        # Note definition of cloud cover : each updraft is associated with
        # a cloud cover equal to the maximum area fraction of the updraft
        # where ql > 0. Each updraft is assumed to have maximum overlap with
        # respect to itup (i.e. no consideration of tilting due to shear)
        # while the updraft classes are assumed to have no overlap at all.
        # Thus total updraft cover is the sum of each updraft's cover

    elseif mode === :up
        cloud_base = zc_toa(grid).z
        cloud_base_ts = ts[k_end]
        @inbounds for k in real_center_indices(grid)
            if TD.has_condensate(aux.q_liq[k] + aux.q_ice[k]) && aux.area[k] > 1e-6
                if cloud_base > grid.zc[k].z
                    cloud_base = grid.zc[k].z
                    cloud_base_ts = ts[k]
                end
            end
        end
    end

    return (; cloud_base_ts, cloud_base)

end

function q_is(
    microphys_params::ACMP,
    ::CMT.IceType,
    r_min::FT = FT(0),
    )  where {FT}
    """
    q_is = 4/3 * π * (r_is^3 - r_min^3) * ρ_i
    """
    return (4 / 3) * π * (CMP.r_ice_snow(microphys_params)^3 - r_min^3) * CMP.ρ_cloud_ice(microphys_params)
end
# no equiv for liquid

function r_ice_acnv(
    microphys_params::ACMP,
    r_ice_acnv_scaling_factor::FT = FT(1.0),
    ) where {FT}
    """
    r_ice_acnv = r_is * r_ice_acnv_scaling_factor
    """
    return CMP.r_ice_snow(microphys_params) * r_ice_acnv_scaling_factor
end
r_ice_acnv(param_set::APS) =
    r_ice_acnv(TCP.microphysics_params(param_set), param_set.user_params.r_ice_acnv_scaling_factor)

function q_acnv_0(
    microphys_params::ACMP,
    ::CMT.IceType,
    r_min::FT = FT(0),
    r_ice_acnv_scaling_factor::FT = FT(1.0),
    ) where {FT}
    """
    q_acnv_0 = 4/3 * π * (r_acnv^3 - r_min^3) * ρ_i
    """
    return (4 / 3) * π * (r_ice_acnv(microphys_params, r_ice_acnv_scaling_factor)^3 - r_min^3) * CMP.ρ_cloud_ice(microphys_params)
end

function q_acnv_0(param_set::APS, type::CMT.IceType, r_min = 0)
    FT = eltype(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    return q_acnv_0(microphys_params, type, FT(r_min), FT(param_set.user_params.r_ice_acnv_scaling_factor))
end

"""

"""
function get_N_threshold(
    param_set::APS, # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.IceType,
    q::TD.PhasePartition,
    T::FT,
    p::FT
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N::FT = FT(NaN),
) where {FT}
    if isnan(N)
        # Let N be the number of droplets of radius r_is that means q = q_threshold
        microphys_params::ACMP = TCP.microphysics_params(param_set)
        q_is_here::FT = q_is(microphys_params, type)
        q_threshold::FT = CMP.q_ice_threshold(microphys_params)
        N = q_threshold / q_is_here # the implied number concentration at threshold for the largest possible droplets... so N and threshold are linked... ( i guess that makes sense? idk...)
    end

    # N = get_N_i(param_set, microphys_params, relaxation_timescale_type, q, T, p, w) # I think this is bad... bc high N w/ high r should lead to more autoconversion..., but this raises the threshold in an uncontrolled way...

    return N

end

function get_N_threshold(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.IceType,
    ts::TD.ThermodynamicState,
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N = NaN,
)
    FT = eltype(param_set)
    thermo_params = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts) 
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts)
    return get_N_threshold(param_set, type, q, T, p; N=N)
end


function get_q_threshold(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    ::CMT.LiquidType,
    q::TD.PhasePartition,
    T::FT,
    p::FT
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N::FT = FT(NaN),
) where {FT}
    return CMP.q_liq_threshold(TCP.microphysics_params(param_set))::FT
end

function get_q_threshold(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.LiquidType,
    ts::TD.ThermodynamicState,
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N = NaN,
)
    return CMP.q_liq_threshold(TCP.microphysics_params(param_set))::eltype(param_set)
end

function get_q_threshold(param_set::APS, type::CMT.LiquidType; N = NaN) # shorthand bc we don't have a fancy q_liq_threshold so we can save some arguments
    # FT = eltype(param_set)
    return CMP.q_liq_threshold(TCP.microphysics_params(param_set))::eltype(param_set)
end

"""
Get the implied q threshold based on r_is, and q 
"""
function get_q_threshold(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.IceType,
    q::TD.PhasePartition,
    T::FT,
    p::FT
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N::FT = FT(NaN),
    scaling_factor::FT = FT(1.0),
) where {FT}

    microphys_params::ACMP = TCP.microphysics_params(param_set)

    if isnan(N)
        return CMP.q_ice_threshold(microphys_params)::FT
    else
        N_thresh::FT = get_N_threshold(param_set, type, q, T, p; N = N) # if N was nothing, get_N_threshold() makes N_thresh * q_is = q_threshold (we could add an if block check for that to save fcn calls...)
        q_is_here::FT = q_is(microphys_params, type)
        return N_thresh * q_is_here * scaling_factor # q_is is the volume of a single ice crystal, so N * q_is = q_threshold
    end
end


function get_q_threshold(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.IceType,
    ts::TD.ThermodynamicState,
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N::FT = eltype(param_set)(NaN),
) where {FT}
    # FT = eltype(param_set)
    thermo_params = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts) 
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts)
    return get_q_threshold(param_set, type, q, T, p; N=FT(N))
end


"""
Get the implied q threshold based on r_is*r_ice_acnv_scaling_factor, and q 
"""
function get_q_threshold_acnv(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.IceType,
    q::TD.PhasePartition,
    T::FT,
    p::FT
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N::FT = FT(NaN),
    assume_N_is::Bool = true,
) where {FT}

    microphys_params::ACMP = TCP.microphysics_params(param_set)

    # get N_is * q_acnv_0, if it's greater than q_ice_threshold, then use that, otherwise use q_ice_threshold

    q_ice_threshold = CMP.q_ice_threshold(microphys_params)

    if isnan(N) && !assume_N_is # if N is NaN, but we don't assume we're at N_is, then we need to use the q_ice_threshold
        # since we assume we're at N_is, the thresh is just q * r_is_acnv_scaling_factor^3
        
        return q_ice_threshold # just take this as given

    else
        # it's not a straight ice_dep_acnv_scaling_factor scaling bc we allow for an r_min that impacts N
        N_thresh::FT = get_N_threshold(param_set, type, q, T, p; N = N) # either take the given N or use the threshold N that assumes we're already at r_is
        q_acnv_0_here::FT = q_acnv_0(param_set, type)
        q_thresh_acnv::FT = N_thresh * q_acnv_0_here # q_is is the volume of a single ice crystal, so N * q_is = q_threshold
    end

    return max(q_thresh_acnv, q_ice_threshold) # make sure we don't get too low of a threshold... (this is a bit arbitrary but I think it makes sense to not let the threshold go below the q_ice threshold)
end

function get_q_threshold_acnv(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    type::CMT.IceType,
    ts::TD.ThermodynamicState,
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    ;
    N = NaN,
    assume_N_is::Bool = true,
)
    FT = eltype(param_set)
    thermo_params = TCP.thermodynamics_params(param_set)
    q::TD.PhasePartition = TD.PhasePartition(thermo_params, ts) 
    T::FT = TD.air_temperature(thermo_params, ts)
    p::FT = TD.air_pressure(thermo_params, ts)
    return get_q_threshold_acnv(param_set, type, q, T, p; N=N, assume_N_is=assume_N_is)
end





# copy of CM1.conv_q_liq_to_q_rai that uses our local get_q_threshold
function my_conv_q_liq_to_q_rai(
    param_set::APS,   # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition,
    T::FT, 
    p::FT;
    N::FT = FT(NaN),
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    return max(0, q.liq - get_q_threshold(param_set, CMT.LiquidType(), q, T, p; N=N)) /
           CMP.τ_acnv_rai(microphys_params)
end

# copy of CM1.conv_q_ice_to_q_sno_no_supersat that uses our local get_q_threshold if we supply N otherwise just uses the default q_ice_threshold
function my_conv_q_ice_to_q_sno_no_supersat(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition,
    T::FT,
    p::FT;
    N::FT = FT(NaN),
    τ::FT = FT(NaN),
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    assume_N_is::Bool = true,
) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    if isnan(τ)
        τ = CMP.τ_acnv_sno(microphys_params)
    end
    return max(0, q.ice - get_q_threshold_acnv(param_set, CMT.IceType(), q, T, p; N=N, assume_N_is = assume_N_is)) / τ
end


function my_conv_q_ice_to_q_sno(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition,
    T::FT,
    p::FT;
    N::FT = FT(NaN),
    τ::FT = FT(NaN),
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
) where {FT}

    acnv_rate = FT(0) 
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)

    ρ::FT = TD.air_density(thermo_params, T, p, q)
    _S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())

    # if any(isnan, (T, p, q.ice, N, Dmin, Dmax))
    #     @warn "NaN in my_conv_q_ice_to_q_sno; T = $T; p = $p; q.ice = $(q.ice); N = $N; τ = $τ; Dmin = $Dmin; Dmax = $Dmax;"
    # end

    if (q.ice > FT(0) && _S > FT(0))

        _G::FT = CM.Common.G_func(microphys_params, T, TD.Ice())


        _r_ice_snow::FT = CMP.r_ice_snow(microphys_params)
        # _n0::FT = n0(microphys_params, FT(0), ρ, ice_type)
        _n0::FT = isnan(N) ? n0(microphys_params, q.ice, ρ, ice_type) : n0(microphys_params, q.ice, ρ, N, ice_type) # use wanted N If given
        _me::FT = me(microphys_params, ice_type)
        _Δm::FT = Δm(microphys_params, ice_type)
        # _λ::FT = lambda(microphys_params, CT.IceType(), q.ice, ρ)
        _λ::FT = lambda(microphys_params, ice_type, q.ice, ρ, N, Dmin, Dmax; _n0 = _n0)

        acnv_rate =
            4 * FT(π) * _S * _G * _n0 / ρ *
            exp(-_λ * _r_ice_snow) *
            (
                _r_ice_snow^FT(2) / (_me + _Δm) +
                (_r_ice_snow * _λ + FT(1)) / _λ^FT(2)
            )

        # if isnan(acnv_rate)
        #     @warn "NaN in my_conv_q_ice_to_q_sno; acnv_rate = $acnv_rate; T = $T; p = $p; q.ice = $(q.ice); N = $N; τ = $τ; Dmin = $Dmin; Dmax = $Dmax; _λ = $_λ; _r_ice_snow = $_r_ice_snow; _n0 = $_n0; _me = $_me; _Δm = $_Δm;"
        # end

    end



    # return acnv_rate
    return resolve_nan(acnv_rate, FT(0)) # if we get a NaN, just return 0... (n0 and lambda can blow up for very small q which map to no conv)
end


"""
Separate fcn for type stability.
    As long as we employ a separate acnv threshold, knowing what that thresh is is important for the limiter.
    Otherwise you can get into a scenario where if τ is very fast, you build up q and then suddenly demolish it all at once because of a fast τ and long enough timestep... this looks like oscillatory `pistoning` and is bad.
   
    You could try to separate out the part of the threshold that is q based even if we have the fallback q_ice_threshold...
    I guess you could argue we shouldn't limit a q based threshold bc it will just decay to 0, but really it's an exponential decay to 0 so it feels like we should stick to our `you can't do everything in one timestep` rule.
    and if the N is from T we really shouldn't go below it in 1 timestep ever.


    To make it simpler and reduce code duplication, although it's slower, we'll make this the main fcn and have my_conv_q_ice_to_q_sno_thresh() just only return the source
"""
function get_thresh_and_my_conv_q_ice_to_q_sno_thresh(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition,
    T::FT,
    p::FT;
    N::FT = FT(NaN),
    τ::FT = FT(NaN),  # consider making this just one timestep?
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    # v_ice::FT = FT(0),
    # v_snow::FT = FT(0),
) where {FT}
    
    thermo_params = TCP.thermodynamics_params(param_set)
    microphys_params::ACMP = TCP.microphysics_params(param_set)

    if isnan(τ)
        τ = CMP.τ_acnv_sno(microphys_params)
    end

    S_qs::FT = FT(0)
    q_thresh::FT = FT(0) # We use this in limiters. if q = 0, S_qi is 0, so it doesn't matter too much but let's just also say 0.

    if q.ice > FT(0)
        
        #=
        Integrating from r_ice_snow to Dmax w/ q_int() doesnt work bc we cant actually create the distribution with q and N below r_ice_snow, only going to infinity (bc incomplete gamma fcns) so we end up getting wild guesses for the distributions that still mostly seem to scale w/ q...
        
        i.e. this doesn't work out:
        μ_ice::FT = FT(0)
        S_qs = q_int(microphys_params, ice_type, q.ice, ρ, μ_ice, Nt=N, Dmin=_r_ice_snow, Dmax=Dmax) # from src/closures/terminal_velocity
        
        Instead we do something like the threshold (q-q_thresh) / τ, but to assess what is crossing the size thresh we scale it by (r/r_is)**3 (maybe the power could be calibrateable but i think it's 2 for area for coll and one for itself being of size to cross r_is)
        
        Note that if we don't/can't predict N, we just assume we can either:
            - assume we're always at r_is and just get 0 agg driven acnv (or use r_acnv scaled down from r_is))
            - assume a fixed thresh.
            The fixed thresh seems a bit subpar but beats adding another parameter I suppose.
        =#
        if isnan(N) 
            # we could also use r_acnv here w/
            # [ assume_N_is means we're always at r_is which means q_thresh is just q but we use q_acnv here... (maybe look at switching that idk, might be out of date now) ]
            # [ not assuming N_is means q_thresh is just q_ice_threshold ]
            S_qs = my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p; N=N, assume_N_is=false) # this is the rate of change of q_ice with respect to N, so we need to divide by τ to get the rate of change of N
            q_thresh = get_q_threshold_acnv(param_set, CMT.IceType(), q, T, p; N=N, assume_N_is = false)
            # r = CMP.r_ice_snow(microphys_params)
        else
            S_qs = my_conv_q_ice_to_q_sno_no_supersat(param_set, q, T, p; N=N, assume_N_is=false) # this is the rate of change of q_ice with respect to N, so we need to divide by τ to get the rate of change of N
            q_thresh = get_q_threshold_acnv(param_set, CMT.IceType(), q, T, p; N=N, assume_N_is = false)
            # r = r_from_qN(q.ice, N, CMP.ρ_cloud_ice(microphys_params); r_min = param_set.user_params.particle_min_radius, mean_r_factor = param_set.user_params.mean_r_factor_ice)
        end

        q_0 = FT(1e-7)
        # S_qs *= S_qs * τ / q_0 # do (q-q_thresh)^2 / τ by squaring the original (q-q_thresh)/τ and multiplying by τ the r/r_is thing also saturates at 1 so no good, only really helps deplete when not in geometric regime. We're looking to instead boost depletion, note squaring could require a huge boost in timescale, since q is small.. we multiply by 1e7 to combat this... it's like 1/q_0

        test = S_qs * (S_qs * τ / q_0)^(param_set.user_params.ice_acnv_power-1) 


        if S_qs > FT(0) # avoid 0^neg power = Inf... (S_qs is 0 or positive)
            S_qs *= (S_qs * τ / q_0)^(param_set.user_params.ice_acnv_power-1) # we want q_0 ((q-q_thresh)/q_0)^p_acnv / τ and we have S_qs = (q-q_thresh)/τ . (S_qs * τ / q_0) ^ (p_acnv-1) = ((q - q_thresh)/q_0)^(p_acnv-1) so we'll multiply and get  (q - q_thresh)^(p_acnv) / (q_0^(p_acnv-1) τ) = q_0 * ((q - q_thresh)/q_0)^(p_acnv) / τ
        end

        if !isfinite(S_qs)
            @warn "test is not finite, S_qs = $S_qs; τ = $τ; q_0 = $q_0; q.ice = $(q.ice); N = $N; T = $T; param_set.user_params.ice_acnv_power = $(param_set.user_params.ice_acnv_power); Dmin = $Dmin; Dmax = $Dmax;"
        end 
        
        # S_qs *= (r / _r_ice_snow)^3 # scale by the size of the particles crossing the threshold (this is a bit arbitrary but I think it makes sense to not let the threshold go below the q_ice threshold)
    end

    return S_qs, q_thresh # return the source term and the threshold
end

"""
    To make it simpler and reduce code duplication, although it's slower, we'll make this the secondary fcn and only return the source
    
    Having duplicate code in my_conv_q_ice_to_q_sno_thresh() and get_thresh_and_my_conv_q_ice_to_q_sno_thresh() would increase the risk of making a mistake trying to keep them in sync, especially if we make changes in my_conv_q_ice_to_q_sno_no_supersat() or get_q_threshold_acnv()
"""
function my_conv_q_ice_to_q_sno_thresh(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    q::TD.PhasePartition,
    T::FT,
    p::FT;
    N::FT = FT(NaN),
    τ::FT = FT(NaN),  # consider making this just one timestep?
    # relaxation_timescale_type::Union{Symbol, AbstractNonEquillibriumSourcesType};
    Dmin::FT = FT(0),
    Dmax::FT = FT(Inf),
    # v_ice::FT = FT(0),
    # v_snow::FT = FT(0),
) where {FT}

    S_qs, _ = get_thresh_and_my_conv_q_ice_to_q_sno_thresh(param_set, q, T, p; N=N, τ=τ, Dmin=Dmin, Dmax=Dmax)
    return S_qs
end



function adjust_ice_N(
    param_set::APS,  # consider changing to just microphys_params since we're not using user_params anymore
    microphys_params::ACMP,
    N_i::FT, 
    q::TD.PhasePartition; 
    r_min::FT=FT(0)
    ) where {FT}
    r_is = CMP.r_ice_snow(microphys_params)
    ρ_i = CMP.ρ_cloud_ice(microphys_params)

    q_is = FT(ρ_i * 4 / 3 * π * (r_is^3 - r_min^3)) # r_is radius ice sphere --
    N_i = max(N_i, q.ice / q_is) # make sure existing ice contributes to the size, particularly w/ sedimentation bringing in new ice

    #= limit N based on q autoconversion threshold =#
    # -- unclear if to do this part
    # We are optin NOT do this rn bc we'll assume if we haven't autoconverted yet, then it really is still in the cloud category, and should affect the timescale, otherwise it's nowhere. powerlaw_T_scaling_ice seemed to work just fine like that... can revisit

    # The argument to do it is bc otherwise, N can grow arbitrarily large... and r can't respnod as fast... lack of autoconversion is more due to r than N surely? but the N could have come from above... so no way around that, and autoconv should eventually take care of it...


    # N_thresh = get_N_threshold(param_set, CMT.IceType(), q, T, p, relaxation_timescale_type; N=nothing) # testing a threshold on how large N_i can be -- basically saying that we can't have more droplets than q_thresh / q(r_is) bc then they should be snow... w/o this limit we just assume autoconv will catch up and take care of it 
    # N_i = min(N_i, N_thresh) # make sure we don't get too many ice crystals

    #= scale N_i by any scaling factor after limiting =#
    #= I don't think this is necessary bc either we're scaling N_i directly which we could just use a constant for or we're scaling the threshold value which is also calibrated =#
    # N_i *= TCP.get_isbits_nt(param_set.user_params, :adjusted_ice_N_scaling_factor, FT(1.0)) # scale N_i by any scaling factor after limiting

    return N_i
end


# r_from_qN(q::FT, N::FT, ρ::FT; r_min::FT=0) = FT(max((q / ((4/3) * π * N * ρ))^(1/3), r_min))
"""
Calculate r from q, N, and ρ.
mean_r_factor exists so you can choose dfiferent estimates of <r> based on the q and N you're using.
For example in gamma distribution, <r> = 1/2λ, q = π ρ N_0 / λ^4 and N = N_0 / λ so r = ( q / ( 8 π N ρ))^(1/3)

We are currently using q = 4/3 π r^3 N ρ, so r = (q / (4/3 π N ρ))^(1/3)

The difference in factor is then  FT(((4/3) / 8)^(1/3)) ≈ 0.55, that is, our r estimates will be double the size they should be, and our timescales will be about half of what they should be.
This sounds small but it doubles any sources... if you rely on having a good prediction of N... this could result in say deposition sources being 2x too large, which isn't insignificant.

"""
function r_from_qN(q::FT, N::FT, ρ::FT; r_min::FT = FT(0), mean_r_factor::FT = FT(1) ) where {FT}
    if iszero(N) || isinf(N)
        return FT(r_min) # N/N will give NaN...
        # even if N = 0 and q is not 0 (in which case N shouldn't be 0), we need some fix...
        # if N = Inf, also bad but don't let r go below r_min, 1/Nr = 1/Inf will yield timescale = 0
        # If N = 0, we'll just return r_min, then 1/Nr = 1/0 will yield timescale = Inf
    end

    return FT((q + 4 / 3 * π * FT(r_min)^3 * ρ * N) / (4 / 3 * π * N * ρ))^(1 / 3) * mean_r_factor
end



"""
We have no way of passing in a minimum aerosol size into λ, N etc in terminal_velocity...
Instead, we calculate what q would be if there was a r_min radius sphere in the center of every drop.
q0 = liquid in one droplet

but then q and N won't match no? you'll have one giant drop... which will have a very high terminal velocity... so adusting N really is the right call...

You would then pass this q into sedimentation to get a more realistic speed for extreme N
"""
function q_effective(
    q::FT,
    N::FT,
    ρ::FT,
    r_min::FT # 
) where {FT}
    """
    effective q for a given N, r, ρ that basically shoves the r_min sphere into the center of every drop
    """
    return q + 4 / 3 * π * ρ * N * (r_min^3) 
end

q_effective(param_set::APS, q::FT, N::FT, ρ::FT) where {FT} = q_effective(q, N, ρ, FT(param_set.user_params.particle_min_radius))


#= 
q if there were no r_min sized aerosols in the middle of every drop, but we rinstead replaced that volume with liquid...
Useful for things that only care about droplet size but estimate it from a q and N only (like my sedimentation terminal velocity implementation)
I suppose you could try to figure how how they're doing that q <--> N mapping, but when we supply our own N, things can get wonky for very large/small N so hopefully this helps.
=#
q_effective_nan_N_safe(q::FT, N::FT, ρ::FT, r_min::FT) where{FT} = isnan(N) ? q : q_effective(q, N, ρ, r_min)
q_effective_nan_N_safe(param_set::APS, q::FT, N::FT, ρ::FT) where{FT} = isnan(N) ? q :  q_effective(q, N, ρ, FT(param_set.user_params.particle_min_radius))
