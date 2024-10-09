

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
    N0::FT,
    p::FT,
    q_liq::FT,
    p_LCL::FT,
    T_LCL::FT,
    q_LCL::FT,
) where {FT}
    """
    # add some DOI references

    # Diagnose q_l from adiabatic parcel ascent from LCL to this altitude
    # Calculate r from this q_l and our assumed N0
    # Calculate real n given our true q_l and this r -- implication is we're subadiabatic and n is depleted by entrainment.

    In a mixed-phase cloud then, WBF would deplete q_l even further leading to reduced n -- an interesting idea considering the Tan et al style pocketing idea of mixed phase clouds, maybe this is better than reducing r homogeneously.

    We have no similar method for ice because we have no expectation of an N0 -- INP are continually activated as a function of temperature.
    The ice N/R relationship maybe is best left as being uniform  since earlier particles would have grown more but there's less of them (don't know relative scaling)... or some other data informed parameterization...
    A potential avenue for thinking about equilibrium adiabat q_l, q_i is Korolev, Field 2006 https://doi.org/10.1175/2007JAS2355.1, but then it would be horrible to reconvert that back to some representation of anything for ice...

    What we'd really want is a lagrangian parcel model that evolves the ice size distribution that we can compare to -- couldn't find obs of ice size dist as fcn of height...



    what happens if the previous LCL is above where we are now? is the calculation order correct to resolve this? (then the cloud base might be at the toa... not ideal...)

    """
    # theta_li is conserved so we just need to back out states based on theta_li and our two pressures
    if p_LCL < p # LCL is above current condensation, so new cloud... use default monodisperse relation
        return NR_monodisperse(N0, q_liq)
    else
        q_LCL = TD.PhasePartition(q_LCL) # assume all vapor at cloud base (or should I read in the real value?) would have to make q_pt an argument rather than q_LCL or only accept ts for example...
        θ_liq_ice = TD.liquid_ice_pottemp_given_pressure(param_set, T_LCL, p_LCL, q_LCL)
        # @show(p, θ_liq_ice)
        q_adiabatic = PhaseEquil_pθq_q(param_set, p, θ_liq_ice, q_LCL.tot)# can we ascend the adiabat here without having to do a full solve? -- also how do I partition liquid and ice in this framework? assume one phase?e
        ρ_w = FT(1)
        r_adiabatic = (q_adiabatic.liq / (4 / 3 * π * ρ_w * N0))^(1 / 3)
        r_adiabatic = max(r_adiabatic, 0.2 * 1e-6) # bound to be at least ~micron size...something like kohler crit radius (e.g. if the cloud is new, diagnosed r is 0 for selecting tau in principle) (also if q_liq is 0 so is r_liq so rly need some bounds)
        N = q_liq / (4 / 3 * π * r_adiabatic^3 * ρ_w)
        return (; N, r = r_adiabatic)
    end
end

function NR_inhomogeneous_mixing_liquid(
    param_set::TD.APS,
    N0::FT,
    p::FT,
    q_liq::FT,
    ts_LCL::TD.ThermodynamicState,
) where {FT}
    p_LCL = TD.air_pressure(param_set, ts_LCL)
    T_LCL = TD.air_temperature(param_set, ts_LCL)
    q_LCL = TD.total_specific_humidity(param_set, ts_LCL)
    return NR_inhomogeneous_mixing_liquid(param_set, N0, p, q_liq, p_LCL, T_LCL, q_LCL)
end

function NR_inhomogeneous_mixing_liquid(
    param_set::TD.APS,
    N0::FT,
    ts::TD.ThermodynamicState,
    q_liq::FT,
    ts_LCL::TD.ThermodynamicState,
) where {FT}
    p = TD.air_pressure(param_set, ts)
    p_LCL = TD.air_pressure(param_set, ts_LCL)
    T_LCL = TD.air_temperature(param_set, ts_LCL)
    q_LCL = TD.total_specific_humidity(param_set, ts_LCL)
    return NR_inhomogeneous_mixing_liquid(param_set, N0, p, q_liq, p_LCL, T_LCL, q_LCL)
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

    # @show(real_center_indices(grid))
    # @show(collect(real_center_indices(grid)))
    # @show(collect(real_center_indices(grid))[end])

    k_end = collect(real_center_indices(grid))[end]

    if mode == :env
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

    elseif mode == :up
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

function q_is(microphys_params::ACMP, ::CMT.IceType)
    return (4 / 3) * π * CMP.r_ice_snow(microphys_params)^3 * CMP.ρ_cloud_ice(microphys_params)
end
# no equiv for liquid

"""

"""
function get_N_threshold(
    param_set::APS,
    type::CMT.IceType,
    q::TD.PhasePartition,
    T::FT,
    p::FT,
    supersat_type::Symbol;
    N0::FT = NaN,
) where {FT}
    if isnothing(N0)
        # Let N0 be the number of droplets of radius r_is that means q = q_threshold
        microphys_params::ACMP = TCP.microphysics_params(param_set)
        q_is_here::FT = q_is(microphys_params, type)
        q_threshold::FT = CMP.q_ice_threshold(microphys_params)
        N0 = q_threshold / q_is_here # the implied number concentration at threshold for the largest possible droplets... so N and threshold are linked... ( i guess that makes sense? idk...)
    end

    # N = get_N_i(param_set, microphys_params, supersat_type, q, T, p, w) # I think this is bad... bc high N w/ high r should lead to more autoconversion..., but this raises the threshold in an uncontrolled way...
    # return  isnothing(N) ? N0 : max(N, N0)

    return N0

end

function get_q_threshold(
    param_set::APS,
    ::CMT.LiquidType,
    q::TD.PhasePartition,
    T::FT,
    p::FT,
    supersat_type::Symbol;
    N0::FT = NaN,
) where {FT}
    return CMP.q_liq_threshold(TCP.microphysics_params(param_set))::FT
end

"""
Get the implied q threshold based on r_is, and q 
"""
function get_q_threshold(
    param_set::APS,
    type::CMT.IceType,
    q::TD.PhasePartition,
    T::FT,
    p::FT,
    supersat_type::Symbol;
    N0::FT = NaN,
) where {FT}

    microphys_params::ACMP = TCP.microphysics_params(param_set)

    if isnan(N0)
        return CMP.q_ice_threshold(microphys_params)::FT
    else
        N_thresh::FT = get_N_threshold(param_set, type, q, T, p, supersat_type; N0 = N0) # if N0 was nothing, get_N_threshold() makes N_thresh * q_is = q_threshold (we could add an if block check for that to save fcn calls...)
        q_is_here::FT = q_is(microphys_params, type)
        return N_thresh * q_is_here
    end
end

# copy of CM1.conv_q_liq_to_q_rai that uses our local get_q_threshold
function my_conv_q_liq_to_q_rai(param_set::APS, q::TD.PhasePartition, T::FT, p::FT, supersat_type::Symbol) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    return max(0, q.liq - get_q_threshold(param_set, CMT.LiquidType(), q, T, p, supersat_type)) /
           CMP.τ_acnv_sno(microphys_params)
end

# copy of CM1.conv_q_ice_to_q_sno_no_supersat that uses our local get_q_threshold
function my_conv_q_ice_to_q_sno_no_supersat(
    param_set::APS,
    q::TD.PhasePartition,
    T::FT,
    p::FT,
    supersat_type::Symbol,
) where {FT}
    microphys_params::ACMP = TCP.microphysics_params(param_set)
    return max(0, q.ice - get_q_threshold(param_set, CMT.IceType(), q, T, p, supersat_type)) /
           CMP.τ_acnv_sno(microphys_params)
end
