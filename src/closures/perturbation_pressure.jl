"""
    perturbation_pressure(
        param_set,
        updraft_top,
        a_up,
        b_up,
        ρ0_k,
        w_up,
        ∇w_up,
        w_en,
        asp_ratio
    )
Returns the value of perturbation pressure gradient
for updraft i following He et al. (JAMES, 2020), given:
 - `updraft_top`: the height of the updraft in the previous timestep
 - `a_up`: updraft area
 - `b_up`: updraft buoyancy
 - `ρ0_k`: reference density
 - `w_up`: updraft vertical velocity
 - `∇w_up`: updraft divergence of vertical velocity
 - `w_en`: environment vertical velocity
 -  esp_ratio`: the specific aspect ratio of the updraft
"""
function perturbation_pressure(
    param_set::APS,
    updraft_top::FT,
    a_up::FT,
    b_up::FT,
    ρ0_k::FT,
    w_up::FT,
    ∇w_up::FT,
    w_en::FT,
    asp_ratio::FT,
) where {FT}

    α_b::FT = CPEDMF.α_b(param_set)
    α_a::FT = CPEDMF.α_a(param_set)
    α_d::FT = CPEDMF.α_d(param_set)
    H_up_min::FT = CPEDMF.H_up_min(param_set)
    α₂ = FT(0)

    nh_press_buoy = -α_b / (1 + α₂ * asp_ratio^2) * ρ0_k * a_up * b_up

    nh_pressure_adv = ρ0_k * a_up * α_a * w_up * ∇w_up
    # drag as w_dif and account for downdrafts
    nh_pressure_drag = -1 * ρ0_k * a_up * α_d * (w_up - w_en) * abs(w_up - w_en) / max(updraft_top, H_up_min)

    return nh_press_buoy, nh_pressure_adv, nh_pressure_drag
end;
