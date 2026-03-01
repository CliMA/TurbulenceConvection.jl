


"""
    compute_CAPE_old!(edmf::EDMFModel, state::State, param_set::APS, grid::Grid; terminate_on_neg_buoyancy::Bool=false)

Legacy CAPE calculation.
Updates:
- **Saturation Preference**: Adds a safety check preventing the swap of a Saturated Parcel to an Unsaturated Environment, regardless of MSE.
- **Debug Block**: Preserved.
"""
function compute_CAPE_old!(
    edmf::EDMFModel,
    state::State,
    param_set::APS,
    grid::Grid;
    terminate_on_neg_buoyancy::Bool = false,
)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)

    # Aliases
    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    ts_env = aux_en.ts

    # Initialize Source with Surface Parcel
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)

    # Track the current parcel being lifted
    current_parcel_θ_liq = TD.liquid_ice_pottemp(thermo_params, ts_env[kc_surf])
    current_parcel_q_tot = TD.total_specific_humidity(thermo_params, ts_env[kc_surf])

    # Track the running energy integral (Kinetic Energy density)
    # We allow this to go negative transiently (CIN), which triggers a restart.
    running_cape = zero(FT)

    print_debug_token = rand() < 1e-3

    aux_gm = center_aux_grid_mean(state)

    # Track previous MSE for gradient calculation
    prev_mse_env = TD.moist_static_energy(thermo_params, ts_env[kc_surf], grid.zc.z[kc_surf])

    # =========================================================================
    # 1. ENVIRONMENTAL CAPE (Restart Logic) - O(N)
    # =========================================================================
    @inbounds for k in real_center_indices(grid)

        if k == kc_toa # we cannot look up to the next layer at the top
            aux_en.CAPE[k] = running_cape
            break
        end

        # A. LIFT PARCEL to current level
        # ---------------------------------------------------------------------
        # p_env = TD.air_pressure(thermo_params, ts_env[k])
        p_env = aux_gm.p[k]
        # if p_env != p_c
        # @warn "In compute_CAPE!: p_env ($p_env) != p_c ($p_c) at k=$(k.i). This may indicate inconsistency between thermodynamic state and pressure field."
        # end

        if edmf.moisture_model isa EquilibriumMoisture
            ts_parcel = TD.PhaseEquil_pθq(thermo_params, p_env, current_parcel_θ_liq, current_parcel_q_tot)
        else
            # CAPE is defined by equilibrium potential :: Note here, `pow_icenuc` matters. `pow_icenuc = FT(Inf) leads to all ice below freezing (more unstable), pow_icenuc = FT(0) leads to no ice (less unstable)`
            # Given slow ice timescales and for simplicity, in this calcluation we want to err to the zero side.
            # ts_parcel = TD.PhaseEquil_pθq(thermo_params, p_env, current_parcel_θ_liq, current_parcel_q_tot)
            # ts_parcel = PhaseEquil_pθq_given_liquid_fraction(thermo_params, p_env, current_parcel_θ_liq, current_parcel_q_tot, one(FT)) # assume all liquid for CAPE calculation, since we do not use pow_icenuc anywhere else in NonEq, and the ice timescale is slower.

            # Non-Equilibrium Logic (Restored)
            # liq_frac = zero(FT) # Assume all ice
            liq_frac = one(FT) # Assume all liquid
            # liq_frac = TD.liquid_fraction(thermo_params, ts_env[k]) # Match environment

            ts_parcel = PhaseEquil_pθq_given_liquid_fraction(
                thermo_params,
                p_env,
                current_parcel_θ_liq,
                current_parcel_q_tot,
                liq_frac,
            ) # assume all liquid for CAPE calculation, since we do not use pow_icenuc anywhere else in NonEq, and the ice timescale is slower.
            q_parcel = PhasePartition_given_liquid_fraction(thermo_params, ts_parcel, liq_frac)
            ts_parcel = TD.PhaseNonEquil_pθq(thermo_params, p_env, current_parcel_θ_liq, q_parcel) # stop pow_icenuc from coming back into the calculation by enforcing nonequil again
        end

        # B. COMPUTE BUOYANCY
        # ---------------------------------------------------------------------
        ρ_env = TD.air_density(thermo_params, ts_env[k])
        ρ_parcel = TD.air_density(thermo_params, ts_parcel)

        # buoyancy_c > 0 if parcel is lighter (ρ_parcel < ρ_env)
        b = buoyancy_c(param_set, ρ_env, ρ_parcel)

        # C. INTEGRATE & RESTART
        # ---------------------------------------------------------------------
        kf = CCO.PlusHalf(k.i)
        dz = grid.zf.z[kf + 1] - grid.zf.z[kf]  # thickness of current layer

        # Add contribution (Positive or Negative)
        # Note: We integrate BEFORE swapping identity to ensure the parcel pays the buoyancy tax
        Δ_CAPE = b * dz

        # Optional: terminating at any stable layer.
        if terminate_on_neg_buoyancy && Δ_CAPE < zero(FT)
            Δ_CAPE = -running_cape
        end

        running_cape += Δ_CAPE

        # RESTART CHECK:
        # If the updraft has been killed by negative buoyancy (running_cape < 0),
        # we reset the process to the local environment. 
        # This allows us to capture elevated convection layers if the surface parcel dies.
        if running_cape ≤ zero(FT)
            running_cape = zero(FT)

            # Reset parcel properties to the LOCAL environment
            ts_local = ts_env[k] # Restart here
            current_parcel_θ_liq = TD.liquid_ice_pottemp(thermo_params, ts_local)
            current_parcel_q_tot = TD.total_specific_humidity(thermo_params, ts_local)
        else
            # D. MUCAPE SWAP LOGIC (Engineered)
            # ---------------------------------------------------------------------
            z_k = grid.zc.z[k]
            mse_parcel = TD.moist_static_energy(thermo_params, ts_parcel, z_k)
            mse_env = TD.moist_static_energy(thermo_params, ts_env[k], z_k)

            # Check Environmental Stability
            # Instability implies MSE decreases with height (dMSE/dz < 0)
            dMSEdz = (mse_env - prev_mse_env) / dz

            # Check Supersaturation (Latent potential)
            q_sat = TD.q_vap_saturation(thermo_params, ts_env[k])
            q_env = TD.total_specific_humidity(thermo_params, ts_env[k])
            is_env_saturated = q_env > q_sat

            # FIX: Check if current parcel is saturated
            parcel_is_saturated = TD.has_condensate(thermo_params, ts_parcel)

            # SAFETY CHECK: 
            # If Parcel is Saturated and Env is NOT, we BLOCK the swap.
            # We only proceed to the MSE check if the Env is "worthy" (also saturated) OR if the Parcel is already dry.
            safe_to_swap = !parcel_is_saturated || is_env_saturated

            # Swap only if Environment is energetically superior AND (Unstable OR Supersaturated)
            # AND if it passes the Saturation Safety Check
            if safe_to_swap && mse_env > mse_parcel && (dMSEdz < 0 || is_env_saturated)
                current_parcel_θ_liq = TD.liquid_ice_pottemp(thermo_params, ts_env[k])
                current_parcel_q_tot = q_env
                # Inherit running_cape, do not reset
            end
        end

        prev_mse_env = TD.moist_static_energy(thermo_params, ts_env[k], grid.zc.z[k])

        # @info "CAPE: k=$(k.i); dz=$dz; p_env=$p_env; ρ_env=$ρ_env; ρ_parcel=$ρ_parcel; b=$b; ΔCAPE=$(Δ_CAPE); current_parcel_θ_liq=$current_parcel_θ_liq; current_parcel_q_tot=$current_parcel_q_tot; ts_parcel=$(ts_parcel); ts_env=$(ts_env[k]); running_cape=$running_cape"
        # println("-----")

        # # --- DEBUG BLOCK START ---
        # print_debug = true
        # if print_debug && print_debug_token
        #      let
        #          # 1. Strict Width Definition (W=12)
        #          W = 12
        #          fmt(x) = rpad(string(round(x, sigdigits=4)), W)

        #          # 2. Header Generators (Guarantees alignment with data)
        #          #    Triad width = 12 + 1 + 12 + 1 + 12 = 38 chars
        #          head_triad() = "$(rpad("P", W)) $(rpad("E", W)) $(rpad("Diff", W))"

        #          # 3. Data Triad Generator
        #          val_triad(p, e) = "$(fmt(p)) $(fmt(e)) $(fmt(p-e))"

        #          println("\n== CAPE DEBUG: k=$(k.i) | z=$(round(grid.zc.z[k], digits=1))m | dz=$(round(dz, digits=1))m | p=$(round(p_env, digits=0))Pa ==")

        #          # --- ROW 1: MOISTURE ---
        #          # Header Row 1 (Category Labels) - padded to 38 chars
        #          println("            | $(rpad("TOTAL Q (kg/kg)", 38))| $(rpad("LIQUID Q (kg/kg)", 38))| $(rpad("ICE Q (kg/kg)", 38))")
        #          # Header Row 2 (P E Diff) - generated to match data exactly
        #          println("            | $(head_triad())| $(head_triad())| $(head_triad())")
        #          println("------------|---------------------------------------|---------------------------------------|---------------------------------------")
        #          # Data Row
        #          println("MOISTURE    | $(val_triad(current_parcel_q_tot, TD.total_specific_humidity(thermo_params, ts_env[k])))| $(val_triad(q_parcel.liq, TD.PhasePartition(thermo_params, ts_env[k]).liq))| $(val_triad(q_parcel.ice, TD.PhasePartition(thermo_params, ts_env[k]).ice))")

        #          # --- ROW 2: THERMO ---
        #          println("            | $(rpad("THETA_V (K)", 38))| $(rpad("THETA_LIQ (K)", 38))| $(rpad("DENSITY (kg/m3)", 38))")
        #          println("THERMO      | $(val_triad(TD.virtual_pottemp(thermo_params, ts_parcel), TD.virtual_pottemp(thermo_params, ts_env[k])))| $(val_triad(current_parcel_θ_liq, TD.liquid_ice_pottemp(thermo_params, ts_env[k])))| $(val_triad(ρ_parcel, ρ_env))")

        #          # --- ROW 3: ENERGETICS ---
        #          println("------------|---------------------------------------|---------------------------------------|---------------------------------------")
        #          println("ENERGETICS  | Buoy (b): $(fmt(b))        | dCAPE: $(fmt(b * dz))        | RunCAPE: $(fmt(running_cape))")
        #          println("======================================================================================================================================\n")
        #      end
        # end
        # # --- DEBUG BLOCK END ---

        # Store result (CAPE is Potential Energy, so we store the positive component)
        aux_en.CAPE[k] = running_cape
    end

    # =========================================================================
    # 2. UPDRAFT & BULK CAPE - O(N)
    # =========================================================================
    N_up = n_updrafts(edmf)

    # Reset bulk CAPE
    @inbounds for k in real_center_indices(grid)
        aux_bulk.CAPE[k] = zero(FT)
    end

    @inbounds for i in 1:N_up
        up_cape = zero(FT)

        @inbounds for k in real_center_indices(grid)
            a_up = aux_up[i].area[k]

            if a_up > FT(1e-5)
                ts_up = aux_up[i].ts[k]
                ρ_up = TD.air_density(thermo_params, ts_up)
                ρ_env = TD.air_density(thermo_params, ts_env[k])

                b_up = buoyancy_c(param_set, ρ_env, ρ_up)
                kf = CCO.PlusHalf(k.i)
                dz = grid.zf.z[kf + 1] - grid.zf.z[kf]

                # Updrafts also accumulate and die
                up_cape += b_up * dz

                if up_cape < 0
                    up_cape = zero(FT)
                end

                aux_up[i].CAPE[k] = up_cape
                aux_bulk.CAPE[k] += a_up * up_cape
            else
                aux_up[i].CAPE[k] = zero(FT)
                up_cape = zero(FT)
            end
        end
    end

    # Normalize Bulk CAPE
    @inbounds for k in real_center_indices(grid)
        total_up_area = sum(aux_up[i].area[k] for i in 1:N_up)
        if total_up_area > FT(1e-5)
            aux_bulk.CAPE[k] /= total_up_area
        else
            aux_bulk.CAPE[k] = 0.0
        end
    end

    return nothing
end


"""
    compute_CAPE!(edmf::EDMFModel, state::State, param_set::APS, grid::Grid; do_quadrature::Bool=false, terminate_on_neg_buoyancy::Bool=false)

Computes CAPE with Dynamic Ensemble Selection (Tournament).

Tournament Logic Hierarchy (Selection Priority):
1. **Survivability**: Alive Parcels (CAPE >= 0) beat Dead Parcels (CAPE < 0).
2. **Fuel (Saturation)**: Saturated Parcels beat Subsaturated Parcels (preventing "Dry Trap").
3. **Potential (MSE)**: Higher MSE wins ties.
4. **Reset (Zombie Fix)**: Dead Rivals beat Dead Survivors (forcing restart).

Additional Mechanics:
- **Conservation**: Energy from valid "Displaced" survivors is redistributed to winners.
- **Gatekeeper**: Matches MUCAPE logic (allows neutral b>=0 swaps if layer is unstable).

terminate_on_neg_buoyancy :: Bool - If true, any layer with negative buoyancy will reset CAPE to zero, effectively killing the updraft immediately. This is a more aggressive criterion for updraft survivability, as it does not allow for any transient negative buoyancy (CIN). Use with caution, but it may help because neutral buoyancy may imply rapid and strong detrainment rather than slow deceleration, which could be more realistic in some cases and prevent excessive CAPE spread and penetration. One day we can add a relaxation rate rather than just an abrupt cutoff to mimic this detrainment.
"""
function compute_CAPE!(
    edmf::EDMFModel,
    state::State,
    param_set::APS,
    grid::Grid;
    do_quadrature::Bool = false,
    terminate_on_neg_buoyancy::Bool = false,
)
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)

    aux_en = center_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    ts_env = aux_en.ts

    kc_surf = kc_surface(grid)
    aux_en.CAPE[kc_surf] = zero(FT)

    # =========================================================================
    # SGS QUADRATURE KERNEL
    # =========================================================================
    function run_cape_core!(::Val{N_QUAD}) where {N_QUAD}
        NUM_SLOTS = N_QUAD * N_QUAD
        POOL_SIZE = 2 * NUM_SLOTS

        # 1. Static Quadrature Maps (Weights/Nodes)
        # ---------------------------------------------------------------------
        inv_pi = FT(1 / π)
        if N_QUAD == 1
            slot_ξ_q = SA.SVector{1, FT}(0)
            slot_ξ_h = SA.SVector{1, FT}(0)
            slot_w = SA.SVector{1, FT}(1)
        else
            n_raw, w_raw = FastGaussQuadrature.gausshermite(N_QUAD)
            nodes = SA.SVector{N_QUAD, FT}(n_raw)
            weights = SA.SVector{N_QUAD, FT}(w_raw)
            tmp_ξ_q = SA.MVector{NUM_SLOTS, FT}(undef)
            tmp_ξ_h = SA.MVector{NUM_SLOTS, FT}(undef)
            tmp_w = SA.MVector{NUM_SLOTS, FT}(undef)
            idx = 1
            for i in 1:N_QUAD, j in 1:N_QUAD
                tmp_ξ_q[idx] = nodes[i]
                tmp_ξ_h[idx] = nodes[j]
                tmp_w[idx] = (weights[i] * weights[j]) * inv_pi
                idx += 1
            end
            slot_ξ_q = SA.SVector(tmp_ξ_q)
            slot_ξ_h = SA.SVector(tmp_ξ_h)
            slot_w = SA.SVector(tmp_w)
        end

        # 2. Allocations
        s_qt = SA.MVector{NUM_SLOTS, FT}(undef)
        s_θl = SA.MVector{NUM_SLOTS, FT}(undef)
        s_cape = SA.MVector{NUM_SLOTS, FT}(undef)
        s_mse = SA.MVector{NUM_SLOTS, FT}(undef)
        s_w = SA.MVector{NUM_SLOTS, FT}(undef)
        p_qt = SA.MVector{POOL_SIZE, FT}(undef)
        p_θl = SA.MVector{POOL_SIZE, FT}(undef)
        p_cape = SA.MVector{POOL_SIZE, FT}(undef)
        p_mse = SA.MVector{POOL_SIZE, FT}(undef)
        p_w = SA.MVector{POOL_SIZE, FT}(undef)
        indices = SA.MVector{POOL_SIZE, Int}(undef)
        p_is_rival = SA.MVector{POOL_SIZE, Int}(undef)
        p_sat = SA.MVector{POOL_SIZE, Bool}(undef)

        # Helper
        function get_parcel_thermo(p_local, θ_liq, q_tot, ts_env_level)
            if edmf.moisture_model isa EquilibriumMoisture
                return TD.PhaseEquil_pθq(thermo_params, p_local, θ_liq, q_tot)
            else
                liq_frac = one(FT)
                ts_p = PhaseEquil_pθq_given_liquid_fraction(thermo_params, p_local, θ_liq, q_tot, liq_frac)
                q_p = PhasePartition_given_liquid_fraction(thermo_params, ts_p, liq_frac)
                return TD.PhaseNonEquil_pθq(thermo_params, p_local, θ_liq, q_p)
            end
        end

        # 3. Surface Init
        ts_src = ts_env[kc_surf]
        qt_mean = TD.total_specific_humidity(thermo_params, ts_src)
        θl_mean = TD.liquid_ice_pottemp(thermo_params, ts_src)
        p_src = TD.air_pressure(thermo_params, ts_src)
        z_src = grid.zc.z[kc_surf]
        sqrt2 = FT(sqrt(2))
        var_q = aux_en.QTvar[kc_surf]
        var_h = aux_en.Hvar[kc_surf]
        σ_q = sqrt(max(var_q, 0))
        σ_h = sqrt(max(var_h, 0))
        if NUM_SLOTS > 1 && σ_q > eps(FT)
            σ_q = min(σ_q, -qt_mean / (sqrt2 * slot_ξ_q[1]))
        end
        r = zero(FT)
        if σ_q > eps(FT) && σ_h > eps(FT)
            r = clamp(aux_en.HQTcov[kc_surf] / (σ_q * σ_h), -1, 1)
        end
        σ_cond = sqrt(max(1 - r^2, 0)) * σ_h

        for i in 1:NUM_SLOTS
            if NUM_SLOTS == 1
                s_qt[i] = qt_mean
                s_θl[i] = θl_mean
            else
                s_qt[i] = max(qt_mean + sqrt2 * σ_q * slot_ξ_q[i], zero(FT))
                s_θl[i] = θl_mean + sqrt2 * (σ_h * r * slot_ξ_q[i] + σ_cond * slot_ξ_h[i])
            end
            s_cape[i] = zero(FT)
            s_w[i] = slot_w[i]
            ts = get_parcel_thermo(p_src, s_θl[i], s_qt[i], ts_src)
            s_mse[i] = TD.moist_static_energy(thermo_params, ts, z_src)
        end

        # 4. Vertical Integration Loop
        @inbounds for k in real_center_indices(grid)
            if k == kc_surf
                continue
            end
            ts_e = ts_env[k]
            p_env = TD.air_pressure(thermo_params, ts_e)
            ρ_env = TD.air_density(thermo_params, ts_e)
            z_k = grid.zc.z[k]
            kf = CCO.PlusHalf(k.i)
            dz = grid.zf.z[kf + 1] - grid.zf.z[kf]

            ts_prev = ts_env[k - 1]
            z_prev = grid.zc.z[k - 1]
            mse_prev = TD.moist_static_energy(thermo_params, ts_prev, z_prev)
            mse_curr = TD.moist_static_energy(thermo_params, ts_e, z_k)
            dMSEdz = (mse_curr - mse_prev) / (z_k - z_prev)
            q_sat = TD.q_vap_saturation(thermo_params, ts_e)
            q_env = TD.total_specific_humidity(thermo_params, ts_e)
            is_unstable_layer = (dMSEdz < 0 || q_env > q_sat)

            # A. Survivors
            for i in 1:NUM_SLOTS
                ts_p = get_parcel_thermo(p_env, s_θl[i], s_qt[i], ts_e)
                ρ_p = TD.air_density(thermo_params, ts_p)
                b = buoyancy_c(param_set, ρ_env, ρ_p)

                # Check Saturation
                p_sat[i] = TD.has_condensate(thermo_params, ts_p)

                # Update CAPE
                if terminate_on_neg_buoyancy && b < 0
                    p_cape[i] = zero(FT)
                    p_mse[i] = -FT(Inf)
                else
                    p_cape[i] = s_cape[i] + b * dz
                    if p_cape[i] < 0
                        p_cape[i] = zero(FT)
                        p_mse[i] = -FT(Inf)
                    else
                        p_mse[i] = TD.moist_static_energy(thermo_params, ts_p, z_k)
                    end
                end
                p_qt[i] = s_qt[i]
                p_θl[i] = s_θl[i]
                p_w[i] = s_w[i]
                p_is_rival[i] = 0
            end

            # B. Rivals
            qt_env_l = TD.total_specific_humidity(thermo_params, ts_e)
            θl_env_l = TD.liquid_ice_pottemp(thermo_params, ts_e)
            l_σ_q = sqrt(max(aux_en.QTvar[k], 0))
            l_σ_h = sqrt(max(aux_en.Hvar[k], 0))
            if NUM_SLOTS > 1 && l_σ_q > eps(FT)
                l_σ_q = min(l_σ_q, -qt_env_l / (sqrt2 * slot_ξ_q[1]))
            end
            l_r = zero(FT)
            if l_σ_q > eps(FT) && l_σ_h > eps(FT)
                l_r = clamp(aux_en.HQTcov[k] / (l_σ_q * l_σ_h), -1, 1)
            end
            l_σ_cond = sqrt(max(1 - l_r^2, 0)) * l_σ_h

            for i in 1:NUM_SLOTS
                idx = NUM_SLOTS + i
                if NUM_SLOTS == 1
                    r_qt = qt_env_l
                    r_θl = θl_env_l
                else
                    r_qt = max(qt_env_l + sqrt2 * l_σ_q * slot_ξ_q[i], zero(FT))
                    r_θl = θl_env_l + sqrt2 * (l_σ_h * l_r * slot_ξ_q[i] + l_σ_cond * slot_ξ_h[i])
                end

                ts_rival = get_parcel_thermo(p_env, r_θl, r_qt, ts_e)
                ρ_rival = TD.air_density(thermo_params, ts_rival)
                b_rival = buoyancy_c(param_set, ρ_env, ρ_rival)

                p_sat[idx] = TD.has_condensate(thermo_params, ts_rival)

                valid_swap = (b_rival > 0) || (b_rival >= 0 && is_unstable_layer)

                if valid_swap
                    p_mse[idx] = TD.moist_static_energy(thermo_params, ts_rival, z_k)
                else
                    p_mse[idx] = -FT(Inf)
                end
                p_qt[idx] = r_qt
                p_θl[idx] = r_θl
                p_w[idx] = slot_w[i]
                p_cape[idx] = zero(FT)
                p_is_rival[idx] = 1
            end

            # C. Tournament (Selection with Saturation Preference)
            for i in 1:POOL_SIZE
                indices[i] = i
            end
            for i in 1:NUM_SLOTS
                best_val = -FT(Inf)
                best_k = i
                for j in i:POOL_SIZE
                    idx = indices[j]
                    is_better = false

                    candidate_dead = (p_mse[idx] == -FT(Inf))
                    best_dead = (best_val == -FT(Inf))

                    if candidate_dead
                        # Zombie Fix: Prefer Rival if both dead (Reset)
                        if best_dead && p_is_rival[idx] == 1 && p_is_rival[best_k] == 0
                            is_better = true
                        end
                    else
                        if best_dead
                            is_better = true
                        else
                            # Alive vs Alive: Preference Saturation
                            cand_sat = p_sat[idx]
                            best_sat = p_sat[best_k]

                            if cand_sat && !best_sat
                                is_better = true
                            elseif !cand_sat && best_sat
                                is_better = false
                            else
                                # Tie on Saturation: Use MSE
                                if p_mse[idx] > best_val
                                    is_better = true
                                end
                            end
                        end
                    end
                    if is_better
                        best_val = p_mse[idx]
                        best_k = j
                    end
                end
                indices[i], indices[best_k] = indices[best_k], indices[i]
            end

            # D. Conservation
            displaced_energy = zero(FT)
            winning_weight = zero(FT)
            for i in (NUM_SLOTS + 1):POOL_SIZE
                loser = indices[i]
                if p_cape[loser] > 0
                    displaced_energy += p_cape[loser] * p_w[loser]
                end
            end
            for i in 1:NUM_SLOTS
                winner = indices[i]
                if p_mse[winner] > -FT(Inf)
                    winning_weight += p_w[winner]
                end
            end

            total_cape = zero(FT)
            total_weight = zero(FT)
            for i in 1:NUM_SLOTS
                winner = indices[i]
                s_qt[i] = p_qt[winner]
                s_θl[i] = p_θl[winner]
                s_w[i] = p_w[winner]
                s_mse[i] = p_mse[winner]
                if p_mse[winner] > -FT(Inf)
                    s_cape[i] = p_cape[winner]
                    if winning_weight > eps(FT)
                        s_cape[i] += displaced_energy / winning_weight
                    end
                    total_cape += s_cape[i] * s_w[i]
                    total_weight += s_w[i]
                else
                    s_cape[i] = zero(FT)
                end
            end
            aux_en.CAPE[k] = total_weight > eps(FT) ? total_cape / total_weight : zero(FT)
        end
    end

    if do_quadrature
        run_cape_core!(Val(2))
    else
        run_cape_core!(Val(1))
    end

    # Legacy Updraft Block [Unchanged]
    N_up = n_updrafts(edmf)
    @inbounds for k in real_center_indices(grid)
        aux_bulk.CAPE[k] = zero(FT)
    end
    @inbounds for i in 1:N_up
        up_cape = zero(FT)
        @inbounds for k in real_center_indices(grid)
            a_up = aux_up[i].area[k]
            if a_up > FT(1e-5)
                ts_up = aux_up[i].ts[k]
                ρ_up = TD.air_density(thermo_params, ts_up)
                ρ_env = TD.air_density(thermo_params, ts_env[k])
                b_up = buoyancy_c(param_set, ρ_env, ρ_up)
                kf = CCO.PlusHalf(k.i)
                dz = grid.zf.z[kf + 1] - grid.zf.z[kf]
                up_cape += b_up * dz
                if up_cape < 0
                    up_cape = zero(FT)
                end
                aux_up[i].CAPE[k] = up_cape
                aux_bulk.CAPE[k] += a_up * up_cape
            else
                aux_up[i].CAPE[k] = zero(FT)
                up_cape = zero(FT)
            end
        end
    end
    @inbounds for k in real_center_indices(grid)
        total_up_area = sum(aux_up[i].area[k] for i in 1:N_up)
        if total_up_area > FT(1e-5)
            aux_bulk.CAPE[k] /= total_up_area
        else
            aux_bulk.CAPE[k] = 0.0
        end
    end
    return nothing
end
