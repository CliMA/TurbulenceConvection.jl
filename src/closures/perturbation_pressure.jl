"""
    perturbation_pressure(
                updraft_top,
                min_updraft_top,
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
 - `updraft_top`: the height of the updraft in the previous timestep
 - `min_updraft_top`: the minimal height of the updraft to avoid zero devision
 - `a_up`: updraft area
 - `b_up`: updraft buoyancy
 - `ρ0_k`: reference density
 - `α₁`: pressure closure free parameter
 - `α₂`: pressure closure free parameter
 - `β₁`: pressure closure free parameter
 - `β₂`: pressure closure free parameter
 - `w_up`: updraft vertical velocity
 - `∇w_up`: updraft divergence of vertical velocity
 - `w_en`: environment vertical velocity
 -  esp_ratio`: the specific aspect ratio of the updraft
"""
function perturbation_pressure(
    updraft_top,
    min_updraft_top,
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
    asp_ratio,
)

    nh_press_buoy = -α₁ / (1 + α₂ * asp_ratio^2) * ρ0_k * a_up * b_up

    nh_pressure_adv = ρ0_k * a_up * β₁ * w_up * ∇w_up
    # drag as w_dif and account for downdrafts
    nh_pressure_drag = -1.0 * ρ0_k * a_up * β₂ * (w_up - w_en) * fabs(w_up - w_en) / fmax(updraft_top, min_updraft_top)

    return nh_press_buoy, nh_pressure_adv, nh_pressure_drag
end;
