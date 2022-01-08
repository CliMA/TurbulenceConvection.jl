# convective velocity scale
get_wstar(bflux, zi) = cbrt(max(bflux * zi, 0))

function buoyancy_c(param_set::APS, ρ0::FT, ρ::FT) where {FT}
    g::FT = CPP.grav(param_set)
    return g * (ρ0 - ρ) / ρ0
end

# BL height
function get_inversion(grid::Grid{FT}, state::State, param_set::APS, Ri_bulk_crit::FT) where {FT}
    g::FT = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    θ_virt = center_aux_turbconv(state).θ_virt
    u = center_prog_grid_mean(state).u
    v = center_prog_grid_mean(state).v
    Ri_bulk = center_aux_bulk(state).Ri
    θ_virt_b = θ_virt[kc_surf]
    z_c = grid.zc
    ∇c = CCO.DivergenceF2C()
    wvec = CC.Geometry.WVector

    # test if we need to look at the free convective limit
    if (u[kc_surf]^2 + v[kc_surf]^2) <= 0.01
        ∇θ_virt = center_aux_turbconv(state).ϕ_temporary
        k_star = findlast_center(k -> θ_virt[k] > θ_virt_b, grid)
        LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(θ_virt[kc_surf]))
        @. ∇θ_virt = ∇c(wvec(LB(θ_virt)))
        h = (θ_virt_b - θ_virt[k_star - 1]) / ∇θ_virt[k_star] + z_c[k_star - 1]
    else
        ∇Ri_bulk = center_aux_turbconv(state).ϕ_temporary
        Ri_bulk_fn(k) = g * (θ_virt[k] - θ_virt_b) * z_c[k] / θ_virt_b / (u[k] * u[k] + v[k] * v[k])

        @inbounds for k in real_center_indices(grid)
            Ri_bulk[k] = Ri_bulk_fn(k)
        end
        k_star = findlast_center(k -> Ri_bulk_fn(k) > Ri_bulk_crit, grid)
        LB = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(Ri_bulk[kc_surf]))
        @. ∇Ri_bulk = ∇c(wvec(LB(Ri_bulk)))
        h = (Ri_bulk_crit - Ri_bulk[k_star - 1]) / ∇Ri_bulk[k_star] + z_c[k_star - 1]
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

function get_surface_tke(ustar::FT, zLL::FT, oblength::FT) where {FT}
    if oblength < 0
        return ((FT(3.75) + cbrt(zLL / oblength * zLL / oblength)) * ustar * ustar)
    else
        return (FT(3.75) * ustar * ustar)
    end
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

function gradient_Richardson_number(∂b∂z::FT, Shear²::FT, ϵ::FT) where {FT}
    return min(∂b∂z / max(Shear², ϵ), FT(0.25))
end

function turbulent_Prandtl_number(param_set::APS, obukhov_length::FT, ∇Ri::FT) where {FT}
    ω_pr::FT = CPEDMF.ω_pr(param_set)
    Pr_n::FT = CPEDMF.Pr_n(param_set)
    if obukhov_length > 0 && ∇Ri > 0 #stable
        # CSB (Dan Li, 2019), with Pr_neutral=0.74 and w1=40.0/13.0
        prandtl_nvec = Pr_n * (2 * ∇Ri / (1 + ω_pr * ∇Ri - sqrt((1 + ω_pr * ∇Ri)^2 - 4 * ∇Ri)))
    else
        prandtl_nvec = Pr_n
    end
    return prandtl_nvec
end
