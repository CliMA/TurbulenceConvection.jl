"""
    buoyancy_gradients(param_set, bg_model::EnvBuoyGrad{FT, EBG}) where {FT <: Real, EBG <: EnvBuoyGradClosure}

Returns the vertical buoyancy gradients in the environment, as well as in its dry and cloudy volume fractions.
The dispatch on EnvBuoyGrad type is performed at the EnvBuoyGrad construction time, and the analytical solutions
used here are consistent for both mean fields and conditional fields obtained from assumed distributions
over the conserved thermodynamic variables.
"""
function buoyancy_gradients(param_set, bg_model::EnvBuoyGrad{FT, EBG}) where {FT <: Real, EBG <: EnvBuoyGradClosure}

    g = CPP.grav(param_set)
    molmass_ratio = CPP.molmass_ratio(param_set)
    R_d = CPP.R_d(param_set)
    R_v = CPP.R_v(param_set)

    phase_part = TD.PhasePartition(0.0, 0.0, 0.0) # assuming R_d = R_m
    Π = TD.exner_given_pressure(param_set, bg_model.p0, phase_part)

    ∂b∂θv = g * (R_d / bg_model.alpha0 / bg_model.p0) * Π

    if bg_model.en_cld_frac > 0.0
        ts_cloudy = thermo_state_pθq(param_set, bg_model.p0, bg_model.θ_liq_ice_cloudy, bg_model.qt_cloudy)
        phase_part = TD.PhasePartition(ts_cloudy)
        lh = TD.latent_heat_liq_ice(param_set, phase_part)
        cp_m = TD.cp_m(ts_cloudy)
        ∂b∂θl_cld = (
            ∂b∂θv * (1 + molmass_ratio * (1 + lh / R_v / bg_model.t_cloudy) * bg_model.qv_cloudy - bg_model.qt_cloudy) /
            (1 + lh * lh / cp_m / R_v / bg_model.t_cloudy / bg_model.t_cloudy * bg_model.qv_cloudy)
        )
        ∂b∂qt_cld = (lh / cp_m / bg_model.t_cloudy * ∂b∂θl_cld - ∂b∂θv) * bg_model.θ_cloudy
    else
        ∂b∂θl_cld = FT(0)
        ∂b∂qt_cld = FT(0)
    end

    ∂b∂z, ∂b∂z_dry, ∂b∂z_cloudy = buoyancy_gradient_chain_rule(bg_model, ∂b∂θv, ∂b∂θl_cld, ∂b∂qt_cld)
    return GradBuoy(∂b∂z, ∂b∂z_dry, ∂b∂z_cloudy)
end

"""
    buoyancy_gradient_chain_rule(
        bg_model::EnvBuoyGrad{FT, EBG},
        ∂b∂θv::FT,
        ∂b∂θl_cld::FT,
        ∂b∂qt_cld::FT,
    ) where {FT <: Real, EBG <: EnvBuoyGradClosure}

Returns the vertical buoyancy gradients in the environment, as well as in its dry and cloudy volume fractions,
from the partial derivatives with respect to thermodynamic variables in dry and cloudy volumes.
"""
function buoyancy_gradient_chain_rule(
    bg_model::EnvBuoyGrad{FT, EBG},
    ∂b∂θv::FT,
    ∂b∂θl_cld::FT,
    ∂b∂qt_cld::FT,
) where {FT <: Real, EBG <: EnvBuoyGradClosure}
    if bg_model.en_cld_frac > FT(0)
        ∂b∂z_θl_cld = ∂b∂θl_cld * bg_model.∂θl∂z_cloudy
        ∂b∂z_qt_cld = ∂b∂qt_cld * bg_model.∂qt∂z_cloudy
    else
        ∂b∂z_θl_cld = FT(0)
        ∂b∂z_qt_cld = FT(0)
    end

    ∂b∂z_dry = bg_model.en_cld_frac < FT(1) ? ∂b∂θv * bg_model.∂θv∂z_dry : FT(0)

    ∂b∂z_cloudy = ∂b∂z_θl_cld + ∂b∂z_qt_cld
    ∂b∂z = (1 - bg_model.en_cld_frac) * ∂b∂z_dry + bg_model.en_cld_frac * ∂b∂z_cloudy

    return ∂b∂z, ∂b∂z_dry, ∂b∂z_cloudy
end
