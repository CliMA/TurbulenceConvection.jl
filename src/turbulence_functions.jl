# convective velocity scale
get_wstar(bflux, zi) = cbrt(max(bflux * zi, 0))

function buoyancy_c(param_set::APS, ρ::FT, ρ_i::FT) where {FT}
    g::FT = TCP.grav(param_set)
    return g * (ρ - ρ_i) / ρ
end

# BL height
function get_inversion(state::State, param_set::APS, Ri_bulk_crit)
    grid = Grid(state)
    FT = float_type(state)
    g::FT = TCP.grav(param_set)
    kc_surf = kc_surface(grid)
    θ_virt = center_aux_grid_mean(state).θ_virt
    u = physical_grid_mean_u(state)
    v = physical_grid_mean_v(state)
    Ri_bulk = center_aux_grid_mean(state).Ri
    θ_virt_b = θ_virt[kc_surf]
    z_c = grid.zc
    ∇c = CCO.DivergenceF2C()
    wvec = CC.Geometry.WVector

    # test if we need to look at the free convective limit
    if (u[kc_surf]^2 + v[kc_surf]^2) <= 0.01
        ∇θ_virt = center_aux_turbconv(state).temporary_1
        k_star = findlast_center(k -> θ_virt[k] > θ_virt_b, grid)
        LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(θ_virt[kc_surf]))
        @. ∇θ_virt = ∇c(wvec(LB(θ_virt)))
        h = (θ_virt_b - θ_virt[k_star - 1]) / ∇θ_virt[k_star] + z_c[k_star - 1].z
    else
        ∇Ri_bulk = center_aux_turbconv(state).temporary_1
        Ri_bulk_fn(k) = g * (θ_virt[k] - θ_virt_b) * z_c[k].z / θ_virt_b / (u[k] * u[k] + v[k] * v[k])

        @inbounds for k in real_center_indices(grid)
            Ri_bulk[k] = Ri_bulk_fn(k)
        end
        k_star = findlast_center(k -> Ri_bulk_fn(k) > Ri_bulk_crit, grid)
        LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(Ri_bulk[kc_surf]))
        @. ∇Ri_bulk = ∇c(wvec(LB(Ri_bulk)))
        h = (Ri_bulk_crit - Ri_bulk[k_star - 1]) / ∇Ri_bulk[k_star] + z_c[k_star - 1].z
    end

    return h
end
# Teixiera convective tau
function get_mixing_tau(zi::FT, wstar::FT) where {FT}
    # return 0.5 * zi / wstar
    #return zi / (max(wstar, FT(1e-5)))
    return zi / (wstar + FT(0.001))
end

# MO scaling of near surface tke and scalar variance
function get_surface_tke(mixing_length_params, ustar::FT, zLL::FT, oblength::FT) where {FT}
    κ_star² = mixing_length_params.κ_star²
    if oblength < 0
        return ((κ_star² + cbrt(zLL / oblength * zLL / oblength)) * ustar * ustar)
    else
        return κ_star² * ustar * ustar
    end
end

function get_surface_tke_turb_flux(mixing_length_params, ustar::FT, ρa_e_surf::FT, U_surf_norm::FT) where {FT}
    κ_star² = mixing_length_params.κ_star²
    c_d = mixing_length_params.c_d
    c_m = mixing_length_params.c_m
    return ρa_e_surf * (1 - c_d * c_m * κ_star²^2) * ustar^2 * U_surf_norm
end

function get_surface_variance(flux1::FT, flux2::FT, ustar::FT, zLL::FT, oblength::FT) where {FT}
    c_star1 = -flux1 / ustar
    c_star2 = -flux2 / ustar
    if oblength < 0
        return 4 * c_star1 * c_star2 * (1 - FT(8.3) * zLL / oblength)^(-FT(2 / 3))
    else
        return 4 * c_star1 * c_star2
    end
end

function gradient_Richardson_number(mixing_length_params, ∂b∂z::FT, Shear²::FT, ϵ::FT) where {FT}
    Ri_c::FT = mixing_length_params.Ri_c
    return min(∂b∂z / max(Shear², ϵ), Ri_c)
end

# Turbulent Prandtl number given the obukhov length sign and the gradient Richardson number
function turbulent_Prandtl_number(mixing_length_params, obukhov_length::FT, ∇Ri::FT) where {FT}
    ω_pr = mixing_length_params.ω_pr
    Pr_n = mixing_length_params.Pr_n
    if obukhov_length > 0 && ∇Ri > 0 #stable
        # CSB (Dan Li, 2019, eq. 75), where ω_pr = ω_1 + 1 = 53.0/13.0
        prandtl_nvec = Pr_n * (2 * ∇Ri / (1 + ω_pr * ∇Ri - sqrt((1 + ω_pr * ∇Ri)^2 - 4 * ∇Ri)))
    else
        prandtl_nvec = Pr_n
    end
    return prandtl_nvec
end

"""
    This allows for MSE gradients to drive instability...
        In this way it is `theoretical` and not immediately apparent...
        but can help drive TKE formation in conditionally unstable moist environments that are globally stable prior to microphysical processes.

    It does assume moisture is present, either for latent heating (condensation/deposition/fusion) or cooling (evaporation/sublimation/melting)
    so you don't need to be at saturation or anything to use it, but you do need to ensure the relevant conditions are met

    ... short derivation here ... 
"""
function moist_gradient_Richardson_number_helper(mixing_length_params, param_set::APS, ∂b∂z::FT, ∂MSE∂z::FT, Shear²::FT, ϵ::FT, T::FT) where {FT}
    g::FT = TCP.grav(param_set)
    c_p::FT = TCP.cp_d(param_set)
    Ri_c::FT = mixing_length_params.Ri_c

    ∂b∂z_moist =  g/(c_p * T) * ∂MSE∂z
    return min(min(∂b∂z, ∂b∂z_moist) / max(Shear², ϵ), Ri_c) # the same as usual but the less stable of the two gradients
end

"""
This one checks which state we're in... for example using moist is ok if subsat given we have condensate but not if we dont.

Perhaps what you'd want to pass in to avoid some complexity is the latent heating rate... but that's not strictly available in Eq ....

"""
function moist_gradient_Richardson_number(mixing_length_params, param_set::APS, ts::TD.ThermodynamicState, qr::FT, qs::FT, ∂b∂z::FT, ∂MSE∂z::FT, Shear²::FT, ϵ::FT; is_noneq::Bool=true, RH_liq::FT=FT(NaN), RH_ice::FT=FT(NaN), T::FT = FT(NaN)) where {FT}
    thermo_params = TCP.thermodynamics_params(param_set)
    T_freeze::FT = TCP.T_freeze(param_set)
    T::FT = isnan(T) ? TD.air_temperature(thermo_params, ts) : T
    q = TD.PhasePartition(thermo_params, ts)

    below_freezing = (T < T_freeze)
    q_min = FT(1e-10) # a reasonable minimum for thermodynamic activity

    if is_noneq
        # This only works in NonEq... I guess in 
        RH_ice = isnan(RH_ice) ? relative_humidity_over_ice(thermo_params, ts) : RH_ice
        RH_liq = isnan(RH_liq) ? relative_humidity_over_liquid(thermo_params, ts) : RH_liq
        if (
            (RH_ice < FT(1)) ||
            ((RH_liq > FT(1)) && below_freezing) || 
            (RH_ice < FT(1) && ((q.ice + qs) < q_min)) ||
            (RH_liq < FT(1) && (q.liq + qr) < q_min)
        )
            return moist_gradient_Richardson_number_helper(mixing_length_params, param_set, ∂b∂z, ∂MSE∂z, Shear², ϵ, T)
        else
            return gradient_Richardson_number(mixing_length_params, ∂b∂z, Shear², ϵ)
        end
    else
        # For the Eq case we can't check for supersat...

        if ((q.ice + qs) < q_min) && ((q.liq + qr) < q_min) # Having condensate implies not subsat in Eq. Hard to know what direction we're going then.
            return gradient_Richardson_number(mixing_length_params, ∂b∂z, Shear², ϵ)
        # saturated and MSE gradient I guess is all we have, there's no supersaturation or lack of condensate to check for
        else # we have condensate. so we're at sat... presumably...
            return moist_gradient_Richardson_number_helper(mixing_length_params, param_set, ∂b∂z, ∂MSE∂z, Shear², ϵ, T)
        end
    end
end