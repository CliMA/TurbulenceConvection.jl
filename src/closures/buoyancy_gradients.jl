"""
    buoyancy_gradients(
        param_set,
        bg_model::EnvBuoyGrad{FT, EBG}
    ) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}

Returns the vertical buoyancy gradients in the environment, as well as in its dry and cloudy volume fractions.
The dispatch on EnvBuoyGrad type is performed at the EnvBuoyGrad construction time, and the analytical solutions
used here are consistent for both mean fields and conditional fields obtained from assumed distributions
over the conserved thermodynamic variables.
"""
function buoyancy_gradients(
    param_set::APS,
    bg_model::EnvBuoyGrad{FT, EBG},
) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}

    g = ICP.grav(param_set)
    molmass_ratio = ICP.molmass_ratio(param_set)
    R_d = ICP.R_d(param_set)
    R_v = ICP.R_v(param_set)

    phase_part = TD.PhasePartition(0.0, 0.0, 0.0) # assuming R_d = R_m
    Π = TD.exner_given_pressure(param_set, bg_model.p0, phase_part)

    ∂b∂θv = g * (R_d / bg_model.alpha0 / bg_model.p0) * Π

    if bg_model.en_cld_frac > 0.0
        ts_sat = thermo_state_pθq(param_set, bg_model.p0, bg_model.θ_liq_ice_sat, bg_model.qt_sat)
        phase_part = TD.PhasePartition(param_set, ts_sat)
        lh = TD.latent_heat_liq_ice(param_set, phase_part)
        cp_m = TD.cp_m(param_set, ts_sat)
        ∂b∂θl_sat = (
            ∂b∂θv * (1 + molmass_ratio * (1 + lh / R_v / bg_model.t_sat) * bg_model.qv_sat - bg_model.qt_sat) /
            (1 + lh * lh / cp_m / R_v / bg_model.t_sat / bg_model.t_sat * bg_model.qv_sat)
        )
        ∂b∂qt_sat = (lh / cp_m / bg_model.t_sat * ∂b∂θl_sat - ∂b∂θv) * bg_model.θ_sat
    else
        ∂b∂θl_sat = FT(0)
        ∂b∂qt_sat = FT(0)
    end

    ∂b∂z, ∂b∂z_unsat, ∂b∂z_sat = buoyancy_gradient_chain_rule(bg_model, ∂b∂θv, ∂b∂θl_sat, ∂b∂qt_sat)
    return GradBuoy{FT}(∂b∂z, ∂b∂z_unsat, ∂b∂z_sat)
end

"""
    buoyancy_gradient_chain_rule(
        bg_model::EnvBuoyGrad{FT, EBG},
        ∂b∂θv::FT,
        ∂b∂θl_sat::FT,
        ∂b∂qt_sat::FT,
    ) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}

Returns the vertical buoyancy gradients in the environment, as well as in its dry and cloudy volume fractions,
from the partial derivatives with respect to thermodynamic variables in dry and cloudy volumes.
"""
function buoyancy_gradient_chain_rule(
    bg_model::EnvBuoyGrad{FT, EBG},
    ∂b∂θv::FT,
    ∂b∂θl_sat::FT,
    ∂b∂qt_sat::FT,
) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}
    if bg_model.en_cld_frac > FT(0)
        ∂b∂z_θl_sat = ∂b∂θl_sat * bg_model.∂θl∂z_sat
        ∂b∂z_qt_sat = ∂b∂qt_sat * bg_model.∂qt∂z_sat
    else
        ∂b∂z_θl_sat = FT(0)
        ∂b∂z_qt_sat = FT(0)
    end

    ∂b∂z_unsat = bg_model.en_cld_frac < FT(1) ? ∂b∂θv * bg_model.∂θv∂z_unsat : FT(0)

    ∂b∂z_sat = ∂b∂z_θl_sat + ∂b∂z_qt_sat
    ∂b∂z = (1 - bg_model.en_cld_frac) * ∂b∂z_unsat + bg_model.en_cld_frac * ∂b∂z_sat

    return ∂b∂z, ∂b∂z_unsat, ∂b∂z_sat
end
