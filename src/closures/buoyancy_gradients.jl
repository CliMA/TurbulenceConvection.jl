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
    g = TCP.grav(param_set)
    R_d = TCP.R_d(param_set)
    R_v = TCP.R_v(param_set)
    T_0 = TCP.T_0(param_set)
    L_v0 = TCP.LH_v0(param_set)
    L_s0 = TCP.LH_s0(param_set)

    if bg_model.en_cld_frac < FT(1)
        T = bg_model.T_unsat
        Tv_gm = bg_model.Tv_gm
        ts_unsat = TD.PhaseEquil_pTq(param_set, bg_model.p, T, bg_model.qt_unsat)
        R_m = TD.gas_constant_air(param_set, ts_unsat)
        cp_m = TD.cp_m(param_set, ts_unsat)

        ∂b∂M_unsat = (g / Tv_gm) * (R_m / R_d / cp_m)
        ∂b∂qt_unsat = (g / Tv_gm / R_d) * (R_m / cp_m * (R_d * T_0 - L_v0) + T * (R_v - R_d))
    else
        ∂b∂M_unsat = FT(0)
        ∂b∂qt_unsat = FT(0)

    end

    if bg_model.en_cld_frac > FT(0)
        T = bg_model.T_sat
        Tv_gm = bg_model.Tv_gm
        ts_sat = TD.PhaseEquil_pTq(param_set, bg_model.p, T, bg_model.qt_sat)
        L_v = TD.latent_heat_vapor(param_set, ts_sat)
        L_s = TD.latent_heat_sublim(param_set, ts_sat)
        L_f = TD.latent_heat_fusion(param_set, ts_sat)
        R_m = TD.gas_constant_air(param_set, ts_sat)
        cp_m = TD.cp_m(param_set, ts_sat)
        qv_star = TD.vapor_specific_humidity(param_set, ts_sat)
        λ = TD.liquid_fraction(param_set, ts_sat)
        L = λ * L_v + (1 - λ) * L_s

        ∂b∂M_sat = (g / Tv_gm) / R_d * (R_m / cp_m + L * qv_star / (cp_m * T)) / (1 + L^2 * qv_star / (R_v * T^2))
        ∂b∂qt_sat = (g / Tv_gm) * ((R_m / cp_m - L * qv_star / (cp_m * T)) * (R_d * T_0 - L_v0) + R_d * T)
    else
        ∂b∂M_sat = FT(0)
        ∂b∂qt_sat = FT(0)
    end

    ∂b∂z, ∂b∂z_unsat, ∂b∂z_sat = buoyancy_gradient_chain_rule(bg_model, ∂b∂M_sat, ∂b∂qt_sat, ∂b∂M_unsat, ∂b∂qt_unsat)
    return GradBuoy{FT}(∂b∂z, ∂b∂z_unsat, ∂b∂z_sat)
end

"""
    buoyancy_gradient_chain_rule(
        bg_model::EnvBuoyGrad{FT, EBG},
        ∂b∂M_sat::FT,
        ∂b∂qt_sat::FT,
        ∂b∂M_unsat::FT,
        ∂b∂qt_unsat::FT,
    ) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}

Returns the vertical buoyancy gradients in the environment, as well as in its dry and cloudy volume fractions,
from the partial derivatives with respect to thermodynamic variables in dry and cloudy volumes.
"""
function buoyancy_gradient_chain_rule(
    bg_model::EnvBuoyGrad{FT, EBG},
    ∂b∂M_sat::FT,
    ∂b∂qt_sat::FT,
    ∂b∂M_unsat::FT,
    ∂b∂qt_unsat::FT,
) where {FT <: Real, EBG <: AbstractEnvBuoyGradClosure}
    if bg_model.en_cld_frac > FT(0)
        ∂b∂z_M_sat = ∂b∂M_sat * bg_model.∂M∂z_sat
        ∂b∂z_qt_sat = ∂b∂qt_sat * bg_model.∂qt∂z_sat
    else
        ∂b∂z_M_sat = FT(0)
        ∂b∂z_qt_sat = FT(0)
    end

    if bg_model.en_cld_frac < FT(1)
        ∂b∂z_M_unsat = ∂b∂M_unsat * bg_model.∂M∂z_unsat
        ∂b∂z_qt_unsat = ∂b∂qt_unsat * bg_model.∂qt∂z_unsat
    else
        ∂b∂z_M_unsat = FT(0)
        ∂b∂z_qt_unsat = FT(0)
    end

    ∂b∂z_sat = ∂b∂z_M_sat + ∂b∂z_qt_sat
    ∂b∂z_unsat = ∂b∂z_M_unsat + ∂b∂z_qt_unsat
    ∂b∂z = (1 - bg_model.en_cld_frac) * ∂b∂z_unsat + bg_model.en_cld_frac * ∂b∂z_sat

    return ∂b∂z, ∂b∂z_unsat, ∂b∂z_sat
end
