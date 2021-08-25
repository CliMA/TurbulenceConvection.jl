function mixing_length(param_set, ml_model::MinDisspLen)
    l = zeros(3)
    c_m = CPEDMF.c_m(param_set)
    c_d = CPEDMF.c_d(param_set)
    smin_ub = CPEDMF.smin_ub(param_set)
    smin_rm = CPEDMF.smin_rm(param_set)
    l_max = ICP.l_max(param_set)
    c_b = ICP.static_stab_coeff(param_set)
    g = CPP.grav(param_set)

    # kz scale (surface layer)
    if ml_model.obukhov_length < 0.0 #unstable
        l_W =
            ml_model.κ_vk * ml_model.z / (sqrt(ml_model.tke_surf / ml_model.ustar / ml_model.ustar) * c_m) *
            min((1.0 - 100.0 * ml_model.z / ml_model.obukhov_length)^0.2, 1.0 / ml_model.κ_vk)
    else # neutral or stable
        l_W = ml_model.κ_vk * ml_model.z / (sqrt(ml_model.tke_surf / ml_model.ustar / ml_model.ustar) * c_m)
    end

    # compute l_TKE - the production/destruction term
    a_pd =
        c_m * (ml_model.Shear² - ml_model.∂b∂z_θl / ml_model.Pr - ml_model.∂b∂z_qt / ml_model.Pr) * sqrt(ml_model.tke)
    # Dissipation term
    c_neg = c_d * ml_model.tke * sqrt(ml_model.tke)
    # Subdomain exchange term
    b_exch = 0.0
    for i in xrange(ml_model.N_up)
        wc_upd_nn = ml_model.wc_up[i]
        wc_env = ml_model.wc_en
        b_exch +=
            ml_model.a_up[i] * wc_upd_nn * ml_model.δ_dyn[i] / ml_model.a_en *
            ((wc_upd_nn - wc_env) * (wc_upd_nn - wc_env) / 2.0 - ml_model.tke) -
            ml_model.a_up[i] * wc_upd_nn * (wc_upd_nn - wc_env) * ml_model.ε_turb[i] * wc_env / ml_model.a_en
    end

    if abs(a_pd) > eps(0.0) && 4.0 * a_pd * c_neg > -b_exch * b_exch
        l_TKE = max(-b_exch / 2.0 / a_pd + sqrt(b_exch * b_exch + 4.0 * a_pd * c_neg) / 2.0 / a_pd, 0.0)
    elseif abs(a_pd) < eps(0.0) && abs(b_exch) > eps(0.0)
        l_TKE = c_neg / b_exch
    else
        l_TKE = 0.0
    end

    # compute l_N² - the effective static stability using environmental mean.
    # Set lambda for now to environmental cloud_fraction (TBD: Rain)
    N²_eff =
        (1.0 - ml_model.en_cld_frac) * ml_model.∂θv∂z +
        ml_model.en_cld_frac * (
            1.0 / exp(-latent_heat(ml_model.T_en) * ml_model.ql_en / cpm_c(ml_model.qt_en) / ml_model.T_en) * (
                (1.0 + (eps_vi - 1.0) * ml_model.qt_en) * ml_model.∂θl∂z +
                (eps_vi - 1.0) * ml_model.θ_li_en * ml_model.∂qt∂z
            )
        )

    N² = sqrt(max(g / ml_model.θv * N²_eff, 0.0))
    if N² > 0.0
        l_N² = min(sqrt(max(c_b * ml_model.tke, 0.0)) / N², l_max)
    else
        l_N² = l_max
    end

    l[1] = l_N²
    l[2] = l_TKE
    l[3] = l_W

    # add limiters 
    l[(l .< eps(0.0)) .| (l .> l_max)] .= l_max

    # get soft minimum
    min_len, min_len_ind = findmin(l)
    mixing_length = lamb_smooth_minimum(l, smin_ub, smin_rm)
    ml_ratio = mixing_length / min_len
    return MixLen(min_len_ind, mixing_length, ml_ratio)
end
