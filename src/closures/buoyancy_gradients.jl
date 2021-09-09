function buoyancy_gradients(param_set, bg_model::Tan2018)

    g = CPP.grav(param_set)

    # buoyancy_gradients
    phase_part = TD.PhasePartition(0.0, 0.0, 0.0) # assuming Rd=Rm
    Π = TD.exner_given_pressure(param_set, bg_model.p0, phase_part)

    prefactor = g * (Rd / bg_model.alpha0 / bg_model.p0) * Π
    ∂b∂θl_dry = prefactor * (1 + (eps_vi - 1) * bg_model.qt_dry)
    ∂b∂qt_dry = prefactor * bg_model.th_dry * (eps_vi - 1)

    if bg_model.en_cld_frac > 0.0
        ts_cloudy = TD.PhaseEquil_pθq_anelastic(param_set, bg_model.p0, bg_model.th_cloudy, bg_model.qt_cloudy)
        lh = TD.latent_heat_vapor(param_set, bg_model.t_cloudy)
        cpm = TD.cp_m(ts_cloudy)
        ∂b∂θl_cld = (
            prefactor * (1 + eps_vi * (1 + lh / Rv / bg_model.t_cloudy) * bg_model.qv_cloudy - bg_model.qt_cloudy) /
            (1 + lh * lh / cpm / Rv / bg_model.t_cloudy / bg_model.t_cloudy * bg_model.qv_cloudy)
        )
        ∂b∂qt_cld = (lh / cpm / bg_model.t_cloudy * ∂b∂θl_cld - prefactor) * bg_model.th_cloudy
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
