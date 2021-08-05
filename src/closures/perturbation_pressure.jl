"""
    perturbation_pressure(
                updraft_top,
                a_up,
                b_up,
                ρ0_k,
                α₁,
                α₂,
                β₁,
                β₂,
                w_up,
                ∇w_up,
                w_en,
                asp_ratio
                )
Returns the value of perturbation pressure gradient
for updraft i following He et al. (JAMES, 2020), given:
 - `updraft_top,`, the heigh of the updraft in the previous timestep
 - `a_up,`, updraft area
 - `b_up,`, updraft buoyancy
 - `ρ0_k,`, reference denbsity
 - `α₁,`, pressure closure free parameter
 - `α₂,`, pressure closure free parameter
 - `β₁,`, pressure closure free parameter
 - `β₂,`, pressure closure free parameter
 - `w_up,`, updraft vertical velocity
 - `∇w_up,`, updraft divergence of vertical velocity
 - `w_en,`, environment vertical velocity
 - `asp_ratio`, the aspect radio of the updraft
"""
function perturbation_pressure(
                updraft_top,
                a_up,
                b_up,
                ρ0_k,
                α₁,
                α₂,
                β₁,
                β₂,
                w_up,
                ∇w_up,
                w_en,
                asp_ratio
                )

    b_coeff = α₁ / (1 + α₂ * asp_ratio^2)
    nh_press_buoy = -b_coeff * b_up * ρ0_k * a_up * b_up

    nh_pressure_adv = ρ0_k * a_up * β₁ * w_up * ∇w_up

    # drag as w_dif and account for downdrafts
    nh_pressure_drag = -1.0 * ρ0_k * a_up * β₂ * (w_up - w_en) * fabs(w_up - w_en) / fmax(updraft_top, 500.0)

    return nh_press_buoy, nh_pressure_adv, nh_pressure_drag, b_coeff
end;
