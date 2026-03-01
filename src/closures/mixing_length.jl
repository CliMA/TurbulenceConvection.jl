function mixing_length(mix_len_params, param_set, ml_model::MinDisspLen{FT}, convective_tke_production::FT) where {FT} # technically if we add a separe convective production we should allow for a separate convective dissipation....
    c_m = mix_len_params.c_m
    c_d = mix_len_params.c_d
    smin_ub = mix_len_params.smin_ub
    l_max = mix_len_params.l_max
    unstable_ref_factor = mix_len_params.unstable_ref_factor
    c_b = mix_len_params.c_b
    g::FT = TCP.grav(param_set)
    molmass_ratio::FT = TCP.molmass_ratio(param_set)
    vkc::FT = TCP.von_karman_const(param_set)
    ustar = ml_model.ustar
    z = ml_model.z
    tke_surf = ml_model.tke_surf
    ∂b∂z = ml_model.∇b.∂b∂z
    tke = ml_model.tke

    # kz scale (surface layer)
    if ml_model.obukhov_length < 0.0 #unstable
        l_W =
            vkc * z / (sqrt(tke_surf / ustar / ustar) * c_m) *
            min((1 - 100 * z / ml_model.obukhov_length)^FT(0.2), 1 / vkc)
    else # neutral or stable
        l_W = vkc * z / (sqrt(tke_surf / ustar / ustar) * c_m)
    end

    # # compute l_TKE - the production-dissipation balanced length scale [[ we use buoy = KH * ∂b∂z, KH = KM / Pr, KM = c_m * l_TKE * sqrt(tke) ]]
    # a_pd = c_m * (ml_model.Shear² - ∂b∂z / ml_model.Pr) * sqrt(tke)
    # # Dissipation term
    # c_neg = c_d * tke * sqrt(tke)
    # # Subdomain exchange term
    # b_exch = ml_model.b_exch

    # if abs(a_pd) > eps(FT) && 4 * a_pd * c_neg > -b_exch * b_exch
    #     l_TKE = max(-b_exch / 2 / a_pd + sqrt(b_exch * b_exch + 4 * a_pd * c_neg) / 2 / a_pd, 0)
    # elseif abs(a_pd) < eps(FT) && abs(b_exch) > eps(FT)
    #     l_TKE = c_neg / b_exch
    # else
    #     l_TKE = FT(0)
    # end

    # now we really have tke production = buoy + tke_production, so we need to add that in
    a_pd = c_m * (ml_model.Shear² - ∂b∂z / ml_model.Pr) * sqrt(tke)
    c_neg = c_d * tke * sqrt(tke)
    b_lin = ml_model.b_exch # + convective_tke_production # this seems to break things since the closure assumes shorter mixing lengths to counter larger production for a given tke amount...
    if abs(a_pd) > eps(FT) && 4 * a_pd * c_neg > -b_lin * b_lin
        l_TKE = max(-b_lin / (2 * a_pd) + sqrt(b_lin * b_lin + 4 * a_pd * c_neg) / (2 * a_pd), 0)
    elseif abs(a_pd) < eps(FT) && abs(b_lin) > eps(FT)
        l_TKE = c_neg / b_lin
    else
        l_TKE = FT(0)
    end




    # # compute l_N - the effective static stability length scale.
    # N_eff = sqrt(max(∂b∂z, 0))
    # if N_eff > 0.0
    #     l_N = min(sqrt(max(c_b * tke, 0)) / N_eff, l_max)
    # else
    #     l_N = l_max
    # end

    # compute l_N - the effective static stability length scale.

    # 1. Configuration (Tunable)
    #    ∂b∂z_ref: The stable gradient where we trust the standard mixing length.
    #    unstable_ref_factor:   How much larger mixing is at Neutral (0) vs Stable (∂b∂z_ref).
    ∂b∂z_ref = FT(1e-5)
    # unstable_ref_factor = FT(1.1)

    N_eff = if ∂b∂z > ∂b∂z_ref
        sqrt(∂b∂z)
    elseif FT(0) <= ∂b∂z <= ∂b∂z_ref # At ∂b∂z = 0, N_eff = sqrt(∂b∂z_ref) / unstable_ref_factor
        # interpolate
        sqrt(∂b∂z_ref) * (1 / unstable_ref_factor + ∂b∂z / ∂b∂z_ref * (1 - 1 / unstable_ref_factor)) # Linearly scales from (1/factor) at 0 up to (1.0) at ∂b∂z_ref
    else # ∂b∂z < 0
        # extrapolate
        # sqrt(∂b∂z_ref) * (1 / factor + ∂b∂z / ∂b∂z_ref * (1 - 1 / unstable_ref_factor)) # Linearly scales from (1/factor) at 0 up to (1.0) at ∂b∂z_ref
        max(sqrt(∂b∂z_ref) * (1 / unstable_ref_factor + ∂b∂z / ∂b∂z_ref * (1 - 1 / unstable_ref_factor)), eps(FT))
    end


    # 3. Calculate Length
    #    This formula now applies to ALL regimes.
    #    - If TKE ~ 0: Numerator is 0 -> l_N = 0. (SOLVES THE JUMP)
    #    - If Unstable: Denominator is small (eps) -> l_N_val is huge -> caps at l_max.
    if tke > eps(FT)
        l_N = min(sqrt(max(c_b * tke, 0)) / N_eff, l_max)
    else
        l_N = FT(0)
    end

    # add limiters
    l = SA.SVector(
        (l_N < eps(FT) || l_N > l_max) ? l_max : l_N,
        (l_TKE < eps(FT) || l_TKE > l_max) ? l_max : l_TKE,
        (l_W < eps(FT) || l_W > l_max) ? l_max : l_W,
    )

    # get soft minimum
    min_len, min_len_ind = findmin(l)
    mix_len = lamb_smooth_minimum(l, smin_ub)
    ml_ratio = mix_len / min_len
    return MixLen(min_len_ind, mix_len, ml_ratio)
end
