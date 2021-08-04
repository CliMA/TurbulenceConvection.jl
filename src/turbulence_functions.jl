# Entrainment Rates

function entr_detr_env_moisture_deficit_b_ED_MF(entr_in::entr_struct)
    _ret = entr_struct()

    moisture_deficit_d =
        (fmax(
            (entr_in.RH_upd / 100.0)^entr_in.sort_pow - (entr_in.RH_env / 100.0)^entr_in.sort_pow,
            0.0,
        ))^(1.0 / entr_in.sort_pow)
    moisture_deficit_e =
        (fmax(
            (entr_in.RH_env / 100.0)^entr_in.sort_pow - (entr_in.RH_upd / 100.0)^entr_in.sort_pow,
            0.0,
        ))^(1.0 / entr_in.sort_pow)
    _ret.sorting_function = moisture_deficit_e
    c_det = entr_in.c_det
    if (entr_in.ql_up + entr_in.ql_env) ≈ 0.0
        c_det = 0.0
    end

    dw = entr_in.w_upd - entr_in.w_env
    if dw < 0.0
        dw -= 0.001
    else
        dw += 0.001
    end

    db = (entr_in.b_upd - entr_in.b_env)
    mu = entr_in.c_mu / entr_in.c_mu0

    inv_timescale = fabs(db / dw)
    logistic_e = 1.0 / (1.0 + exp(-mu * db / dw * (entr_in.chi_upd - entr_in.a_upd / (entr_in.a_upd + entr_in.a_env))))
    logistic_d = 1.0 / (1.0 + exp(mu * db / dw * (entr_in.chi_upd - entr_in.a_upd / (entr_in.a_upd + entr_in.a_env))))

    #Logistic of buoyancy fluxes
    inv_timescale = fabs(db / dw)
    ed_mf_ratio =
        fabs(entr_in.buoy_ed_flux) /
        (fabs(entr_in.a_upd * entr_in.a_env * (entr_in.w_upd - entr_in.w_env) * (entr_in.b_upd - entr_in.b_env)) + 1e-8)
    logistic_e *= (1.0 / (1.0 + exp(entr_in.c_ed_mf * (ed_mf_ratio - 1.0))))
    _ret.entr_sc = inv_timescale / dw * (entr_in.c_ent * logistic_e + c_det * moisture_deficit_e)
    _ret.detr_sc = inv_timescale / dw * (entr_in.c_ent * logistic_d + c_det * moisture_deficit_d)

    return _ret
end

function entr_detr_env_moisture_deficit(entr_in::entr_in_struct)
    _ret = entr_struct()
    l = pyzeros(2)

    moisture_deficit_d =
        (fmax(
            (entr_in.RH_upd / 100.0)^entr_in.sort_pow - (entr_in.RH_env / 100.0)^entr_in.sort_pow,
            0.0,
        ))^(1.0 / entr_in.sort_pow)
    moisture_deficit_e =
        (fmax(
            (entr_in.RH_env / 100.0)^entr_in.sort_pow - (entr_in.RH_upd / 100.0)^entr_in.sort_pow,
            0.0,
        ))^(1.0 / entr_in.sort_pow)
    _ret.sorting_function = moisture_deficit_e
    c_det = entr_in.c_det
    if (entr_in.ql_up + entr_in.ql_env) == 0.0
        c_det = 0.0
    end

    dw = entr_in.w_upd - entr_in.w_env
    if dw < 0.0
        dw -= 0.001
    else
        dw += 0.001
    end

    db = (entr_in.b_upd - entr_in.b_env)
    mu = entr_in.c_mu / entr_in.c_mu0

    inv_timescale = fabs(db / dw)
    logistic_e = 1.0 / (1.0 + exp(-mu * db / dw * (entr_in.chi_upd - entr_in.a_upd / (entr_in.a_upd + entr_in.a_env))))
    logistic_d = 1.0 / (1.0 + exp(mu * db / dw * (entr_in.chi_upd - entr_in.a_upd / (entr_in.a_upd + entr_in.a_env))))

    #smooth min
    l[0] = entr_in.tke_coef * fabs(db / sqrt(entr_in.tke + 1e-8))
    l[1] = fabs(db / dw)
    inv_timescale = lamb_smooth_minimum(l, 0.1, 0.0005)
    _ret.entr_sc = inv_timescale / dw * (entr_in.c_ent * logistic_e + c_det * moisture_deficit_e)
    _ret.detr_sc = inv_timescale / dw * (entr_in.c_ent * logistic_d + c_det * moisture_deficit_d)

    return _ret
end

function entr_detr_env_moisture_deficit_div(entr_in::entr_struct)
    _ret = entr_struct()
    l = pyzeros(2)

    moisture_deficit_d =
        (fmax(
            (entr_in.RH_upd / 100.0)^entr_in.sort_pow - (entr_in.RH_env / 100.0)^entr_in.sort_pow,
            0.0,
        ))^(1.0 / entr_in.sort_pow)
    moisture_deficit_e =
        (fmax(
            (entr_in.RH_env / 100.0)^entr_in.sort_pow - (entr_in.RH_upd / 100.0)^entr_in.sort_pow,
            0.0,
        ))^(1.0 / entr_in.sort_pow)
    _ret.sorting_function = moisture_deficit_e
    c_det = entr_in.c_det
    if (entr_in.ql_up + entr_in.ql_env) == 0.0
        c_det = 0.0
    end

    dw = entr_in.w_upd - entr_in.w_env
    if dw < 0.0
        dw -= 0.001
    else
        dw += 0.001
    end

    db = (entr_in.b_upd - entr_in.b_env)
    mu = entr_in.c_mu / entr_in.c_mu0

    inv_timescale = fabs(db / dw)
    logistic_e = 1.0 / (1.0 + exp(-mu * db / dw * (entr_in.chi_upd - entr_in.a_upd / (entr_in.a_upd + entr_in.a_env))))
    logistic_d = 1.0 / (1.0 + exp(mu * db / dw * (entr_in.chi_upd - entr_in.a_upd / (entr_in.a_upd + entr_in.a_env))))

    entr_MdMdz = fmax(entr_in.dMdz / fmax(entr_in.M, 1e-12), 0.0)
    detr_MdMdz = fmax(-entr_in.dMdz / fmax(entr_in.M, 1e-12), 0.0)


    #smooth min
    l[0] = entr_in.tke_coef * fabs(db / sqrt(entr_in.tke + 1e-8))
    l[1] = fabs(db / dw)
    inv_timescale = lamb_smooth_minimum(l, 0.1, 0.0005)

    _ret.entr_sc =
        inv_timescale / dw * (entr_in.c_ent * logistic_e + c_det * moisture_deficit_e) + entr_MdMdz * entr_in.c_div
    _ret.detr_sc =
        inv_timescale / dw * (entr_in.c_ent * logistic_d + c_det * moisture_deficit_d) + detr_MdMdz * entr_in.c_div


    return _ret
end


function pressure_normalmode_buoy(press_in::pressure_in_struct)
    _ret = pressure_buoy_struct()

    _ret.b_coeff = press_in.alpha1 / (1 + press_in.alpha2 * press_in.asp_ratio^2)
    _ret.nh_pressure_b = -1.0 * press_in.rho0_kfull * press_in.a_kfull * press_in.b_kfull * _ret.b_coeff

    return _ret
end

function pressure_normalmode_drag(press_in::pressure_in_struct)
    _ret = pressure_drag_struct()

    _ret.nh_pressure_adv =
        press_in.rho0_kfull *
        press_in.a_kfull *
        press_in.beta1 *
        press_in.w_kfull *
        (press_in.w_kfull - press_in.w_kmfull) *
        press_in.dzi

    # drag as w_dif and account for downdrafts
    _ret.nh_pressure_drag =
        -1.0 *
        press_in.rho0_kfull *
        press_in.a_kfull *
        press_in.beta2 *
        (press_in.w_kfull - press_in.w_kenv) *
        fabs(press_in.w_kfull - press_in.w_kenv) / fmax(press_in.updraft_top, 500.0)

    return _ret
end

# convective velocity scale
function get_wstar(bflux, zi)
    return cbrt(fmax(bflux * zi, 0.0))
end

# BL height
function get_inversion(theta_rho, u, v, grid::Grid, Ri_bulk_crit)
    gw = grid.gw
    kmin = gw
    theta_rho_b = theta_rho[kmin]
    k = kmin
    Ri_bulk = 0.0
    Ri_bulk_low = 0.0
    z_half = grid.z_half

    k_last = last(real_center_indicies(grid))
    # test if we need to look at the free convective limit
    if (u[kmin] * u[kmin] + v[kmin] * v[kmin]) <= 0.01
        @inbounds for k in real_center_indicies(grid)
            if theta_rho[k] > theta_rho_b
                k_last = k
                break
            end
        end
        k = k_last
        h =
            (z_half[k] - z_half[k - 1]) / (theta_rho[k] - theta_rho[k - 1]) * (theta_rho_b - theta_rho[k - 1]) +
            z_half[k - 1]
    else
        @inbounds for k in real_center_indicies(grid)
            Ri_bulk_low = Ri_bulk
            Ri_bulk = g * (theta_rho[k] - theta_rho_b) * z_half[k] / theta_rho_b / (u[k] * u[k] + v[k] * v[k])
            if Ri_bulk > Ri_bulk_crit
                k_last = k
                break
            end
        end
        k = k_last
        h = (z_half[k] - z_half[k - 1]) / (Ri_bulk - Ri_bulk_low) * (Ri_bulk_crit - Ri_bulk_low) + z_half[k - 1]
    end

    return h
end
# Teixiera convective tau
function get_mixing_tau(zi, wstar)
    # return 0.5 * zi / wstar
    #return zi / (fmax(wstar, 1e-5))
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
        return 4.0 * c_star1 * c_star2 * pow(1.0 - 8.3 * zLL / oblength, -2.0 / 3.0)
    else
        return 4.0 * c_star1 * c_star2
    end
end

function construct_tridiag_diffusion(grid, dt, rho_ae_K_m, rho, ae, a, b, c)
    dzi = grid.dzi
    @inbounds for k in real_center_indicies(grid)
        X = rho[k] * ae[k] / dt
        Y = rho_ae_K_m[k] * dzi * dzi
        Z = rho_ae_K_m[k - 1] * dzi * dzi
        if is_surface_bc_centers(grid, k)
            Z = 0.0
        elseif is_toa_bc_centers(grid, k)
            Y = 0.0
        end
        a[k] = -Z / X
        b[k] = 1.0 + Y / X + Z / X
        c[k] = -Y / X
    end
    return
end

function tridiag_solve(b_rhs, a, b, c)
    # Note that `1:end` is zero-based indexing.
    A = Tridiagonal(a[1:end], parent(b), c[0:(end - 1)])
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
