function buoyancy_gradients(param_set, bg_model::Tan2018)

    g = CPP.grav(param_set)
    molmass_ratio = CPP.molmass_ratio(param_set)
    R_d = CPP.R_d(param_set)
    R_v = CPP.R_v(param_set)

    # buoyancy_gradients
    phase_part = TD.PhasePartition(0.0, 0.0, 0.0) # assuming R_d = R_m
    Π = TD.exner_given_pressure(param_set, bg_model.p0, phase_part)

    ∂b∂θv = g * (R_d / bg_model.alpha0 / bg_model.p0) * Π
    ∂b∂θl_dry = ∂b∂θv * (1 + (molmass_ratio - 1) * bg_model.qt_dry)
    ∂b∂qt_dry = ∂b∂θv * bg_model.th_dry * (molmass_ratio - 1)

    if bg_model.en_cld_frac > 0.0
        ts_cloudy = TD.PhaseEquil_pθq(param_set, bg_model.p0, bg_model.th_cloudy, bg_model.qt_cloudy)
        lh = TD.latent_heat_vapor(param_set, bg_model.t_cloudy)
        cp_m = TD.cp_m(ts_cloudy)
        ∂b∂θl_cld = (
            ∂b∂θv * (1 + molmass_ratio * (1 + lh / R_v / bg_model.t_cloudy) * bg_model.qv_cloudy - bg_model.qt_cloudy) /
            (1 + lh * lh / cp_m / R_v / bg_model.t_cloudy / bg_model.t_cloudy * bg_model.qv_cloudy)
        )
        ∂b∂qt_cld = (lh / cp_m / bg_model.t_cloudy * ∂b∂θl_cld - ∂b∂θv) * bg_model.th_cloudy
    else
        ∂b∂θl_cld = 0
        ∂b∂qt_cld = 0
    end
    ∂b∂θl = (bg_model.en_cld_frac * ∂b∂θl_cld + (1 - bg_model.en_cld_frac) * ∂b∂θl_dry)
    ∂b∂qt = (bg_model.en_cld_frac * ∂b∂qt_cld + (1 - bg_model.en_cld_frac) * ∂b∂qt_dry)

    ∂b∂z_θl = ∂b∂θl * bg_model.∂θl∂z
    ∂b∂z_qt = ∂b∂qt * bg_model.∂qt∂z

    return GradBuoy(∂b∂θl, ∂b∂qt, ∂b∂z_θl, ∂b∂z_qt)
end
