# convective velocity scale
function get_wstar(bflux, zi)
    return cbrt(max(bflux * zi, 0.0))
end

# BL height
function get_inversion(param_set, θ_ρ, u, v, grid::Grid, Ri_bulk_crit)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    θ_ρ_b = θ_ρ[kc_surf]
    Ri_bulk = center_field(grid)
    z_c = grid.zc

    # test if we need to look at the free convective limit
    if (u[kc_surf]^2 + v[kc_surf]^2) <= 0.01
        kmask = map(k -> (k, θ_ρ[k] > θ_ρ_b), real_center_indices(grid))
        k_star = first(kmask[findlast(km -> km[2], kmask)])
        ∇θ_ρ = c∇_upwind(θ_ρ, grid, k_star; bottom = SetGradient(0), top = FreeBoundary())
        h = (θ_ρ_b - θ_ρ[k_star - 1]) / ∇θ_ρ + z_c[k_star - 1]
    else
        Ri_bulk_fn(k) = g * (θ_ρ[k] - θ_ρ_b) * z_c[k] / θ_ρ_b / (u[k] * u[k] + v[k] * v[k])

        Ri_bulk[real_center_indices(grid)] .= map(k -> Ri_bulk_fn(k), real_center_indices(grid))
        kmask = map(k -> (k, Ri_bulk_fn(k) > Ri_bulk_crit), real_center_indices(grid))
        k_star = first(kmask[findlast(km -> km[2], kmask)])
        ∇Ri_bulk = c∇_upwind(Ri_bulk, grid, k_star; bottom = SetGradient(0), top = FreeBoundary())
        h = (Ri_bulk_crit - Ri_bulk[k_star - 1]) / ∇Ri_bulk + z_c[k_star - 1]
    end

    return h
end
# Teixiera convective tau
function get_mixing_tau(zi, wstar)
    # return 0.5 * zi / wstar
    #return zi / (max(wstar, 1e-5))
    return zi / (wstar + 0.001)
end

# MO scaling of near surface tke and scalar variance

function get_surface_tke(ustar, wstar, zLL, oblength)
    if oblength < 0.0
        return ((3.75 + cbrt(zLL / oblength * zLL / oblength)) * ustar * ustar)
    else
        return (3.75 * ustar * ustar)
    end
end
function get_surface_variance(flux1, flux2, ustar, zLL, oblength)
    c_star1 = -flux1 / ustar
    c_star2 = -flux2 / ustar
    if oblength < 0.0
        return 4.0 * c_star1 * c_star2 * (1.0 - 8.3 * zLL / oblength)^(-2.0 / 3.0)
    else
        return 4.0 * c_star1 * c_star2
    end
end

function construct_tridiag_diffusion(grid, dt, ρ_ae_K_m, ρ_0, ae)
    a = center_field(grid) # for tridiag solver
    b = center_field(grid) # for tridiag solver
    c = center_field(grid) # for tridiag solver
    Δzi = grid.Δzi
    @inbounds for k in real_center_indices(grid)
        X = ρ_0[k] * ae[k] / dt
        Y = ρ_ae_K_m[k + 1] * Δzi * Δzi
        Z = ρ_ae_K_m[k] * Δzi * Δzi
        if is_surface_center(grid, k)
            Z = 0.0
        elseif is_toa_center(grid, k)
            Y = 0.0
        end
        a[k] = -Z / X
        b[k] = 1.0 + Y / X + Z / X
        c[k] = -Y / X
    end
    A = LinearAlgebra.Tridiagonal(a[2:end], vec(b), c[1:(end - 1)])
    return A
end

tridiag_solve(b_rhs, A) = A \ b_rhs

# Still need this temporarily
function tridiag_solve(b_rhs, a, b, c)
    A = LinearAlgebra.Tridiagonal(a[2:end], parent(b), c[1:(end - 1)])
    return A \ parent(b_rhs)
end

# Dustbin

function set_cloudbase_flag(ql, current_flag)
    if ql > 1.0e-8
        new_flag = true
    else
        new_flag = current_flag
    end
    return new_flag
end

function gradient_Richardson_number(∂b∂θ_l, ∂b∂q_tot, Shear², ϵ)
    return min(∂b∂θ_l / max(Shear², ϵ) + ∂b∂q_tot / max(Shear², ϵ), 0.25)
end

function turbulent_Prandtl_number(obukhov_length, ∇Ri, Pr, ω_pr)
    if obukhov_length > 0.0 && ∇Ri > 0.0 #stable
        # CSB (Dan Li, 2019), with Pr_neutral=0.74 and w1=40.0/13.0
        prandtl_nvec = Pr * (2 * ∇Ri / (1 + ω_pr * ∇Ri - sqrt((1 + ω_pr * ∇Ri)^2 - 4 * ∇Ri)))
    else
        prandtl_nvec = Pr
    end
    return prandtl_nvec
end
