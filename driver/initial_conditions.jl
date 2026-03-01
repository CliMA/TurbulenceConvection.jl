import TurbulenceConvection as TC
import TurbulenceConvection.Parameters as TCP
const APS = TCP.AbstractTurbulenceConvectionParameters

import Thermodynamics as TD

function initialize_edmf(edmf::TC.EDMFModel, state::TC.State, surf_params, param_set::APS, t::Real, case)
    grid = TC.Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    initialize_covariance(edmf, state)
    aux_gm = TC.center_aux_grid_mean(state)
    ts_gm = aux_gm.ts
    @. aux_gm.θ_virt = TD.virtual_pottemp(thermo_params, ts_gm)
    surf = get_surface(surf_params, state, t, param_set)
    if case isa Cases.DryBubble
        initialize_updrafts_DryBubble(edmf, state, surf)
    elseif case isa Cases.SOCRATES
        initialize_updrafts_SOCRATES(edmf, state, surf, param_set) # testing so i can start updrafts w/o 0 area
    else
        initialize_updrafts(edmf, state, surf)
    end
    TC.set_edmf_surface_bc(edmf, state, surf, param_set)
    return
end

function initialize_covariance(edmf::TC.EDMFModel, state::TC.State)
    grid = TC.Grid(state)
    kc_surf = TC.kc_surface(grid)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_en = TC.center_prog_environment(state)
    aux_en = TC.center_aux_environment(state)
    prog_gm = TC.center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ
    aux_bulk = TC.center_aux_bulk(state)
    ae = 1 .- aux_bulk.area # area of environment

    aux_en.tke .= aux_gm.tke
    aux_en.Hvar .= aux_gm.Hvar
    aux_en.QTvar .= aux_gm.QTvar
    aux_en.HQTcov .= aux_gm.HQTcov

    prog_en.ρatke .= aux_en.tke .* ρ_c .* ae
    if edmf.thermo_covariance_model isa TC.PrognosticThermoCovariances
        prog_en.ρaHvar .= aux_gm.Hvar .* ρ_c .* ae
        prog_en.ρaQTvar .= aux_gm.QTvar .* ρ_c .* ae
        prog_en.ρaHQTcov .= aux_gm.HQTcov .* ρ_c .* ae
    end
    return
end

function initialize_updrafts(edmf, state, surf)
    grid = TC.Grid(state)
    N_up = TC.n_updrafts(edmf)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    aux_up = TC.center_aux_updrafts(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_up = TC.center_prog_updrafts(state)
    prog_up_f = TC.face_prog_updrafts(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    FT = TC.float_type(state)
    @inbounds for i in 1:N_up
        @inbounds for k in TC.real_face_indices(grid)
            aux_up_f[i].w[k] = 0
            prog_up_f[i].ρaw[k] = 0
        end

        @inbounds for k in TC.real_center_indices(grid)
            aux_up[i].buoy[k] = 0
            # Simple treatment for now, revise when multiple updraft closures
            # become more well defined
            aux_up[i].area[k] = 0
            aux_up[i].q_tot[k] = aux_gm.q_tot[k]
            aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
            aux_up[i].q_liq[k] = aux_gm.q_liq[k]
            aux_up[i].q_ice[k] = aux_gm.q_ice[k]
            aux_up[i].T[k] = aux_gm.T[k]
            prog_up[i].ρarea[k] = 0
            prog_up[i].ρaq_tot[k] = 0
            prog_up[i].ρaθ_liq_ice[k] = 0
        end
        if edmf.entr_closure isa TC.PrognosticNoisyRelaxationProcess
            @. prog_up[i].ε_nondim = 0
            @. prog_up[i].δ_nondim = 0
        end

        if edmf.surface_area_bc isa TC.FixedSurfaceAreaBC || edmf.surface_area_bc isa TC.ClosureSurfaceAreaBC
            a_surf = TC.area_surface_bc(surf, edmf, i, edmf.surface_area_bc)
            aux_up[i].area[kc_surf] = a_surf
            prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * a_surf

            # To avoid degeneracy from 0 BC on surface updraft area fraction, add small perturbation
            # to relevant prognostic variables in bottom cell -> This avoids trivial solution with no updrafts.
        elseif edmf.surface_area_bc isa TC.PrognosticSurfaceAreaBC
            a_up_initial = FT(0.1)
            aux_up[i].area[kc_surf] = a_up_initial
            prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * a_up_initial

            aux_up_f[i].w[kf_surf + 1] = FT(0.01) # small velocity perturbation
            prog_up_f[i].ρaw[kf_surf + 1] = ρ_f[kf_surf + 1] * a_up_initial * FT(0.01) # small ρaw perturbation

            prog_up[i].ρaq_tot[kc_surf] = ρ_c[kc_surf] * a_up_initial * aux_gm.q_tot[kc_surf]
            prog_up[i].ρaθ_liq_ice[kc_surf] = ρ_c[kc_surf] * a_up_initial * aux_gm.θ_liq_ice[kc_surf]
        end
    end
    return
end

# my own testing fcn for changing initial area
function initialize_updrafts_SOCRATES(edmf, state, surf, param_set)
    grid = TC.Grid(state)
    N_up = TC.n_updrafts(edmf)
    kc_surf = TC.kc_surface(grid)
    kf_surf = TC.kf_surface(grid)
    aux_up = TC.center_aux_updrafts(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_up = TC.center_prog_updrafts(state)
    prog_up_f = TC.face_prog_updrafts(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    FT = TC.float_type(state)

    C2F = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0))) # face to center interpolation

    kc_toa = TC.kc_top_of_atmos(grid)
    kf_toa = TC.kf_top_of_atmos(grid)


    k_mf = CCO.PlusHalf(kc_surf.i)
    k_pf = CCO.PlusHalf((kc_surf + 1).i)
    @info "kc_surf = $kc_surf; kf_surf = $kf_surf; kc_toa = $kc_toa; kf_toa = $kf_toa; k_mf = $k_mf; k_pf = $k_pf"

    initial_profile_updraft_area::FT = TCP.get_isbits_nt(param_set.user_params, :initial_profile_updraft_area, FT(0))
    stable_updraft_area_reduction_factor::FT =
        TCP.get_isbits_nt(param_set.user_params, :stable_updraft_area_reduction_factor, FT(100))
    @inbounds for i in 1:N_up
        @inbounds for k in TC.real_face_indices(grid)
            # aux_up_f[i].w[k] = 0
            # prog_up_f[i].ρaw[k] = 0

            # aux_up_f[i].w[k] = FT(0.01) # a small velocity perturbation [needed so area doesn't immediately go to 0]
            aux_up_f[i].w[k] = FT(0)
            # prog_up_f[i].ρaw[k] = ρ_f[k] * initial_profile_updraft_area * aux_up_f[i].w[k] # small ρaw perturbation
        end

        Φ = TC.geopotential.(param_set, getfield.(grid.zc, :z))
        MSE_gm = TD.moist_static_energy.(TCP.thermodynamics_params(param_set), aux_gm.ts, Φ)
        # ∇0_bcs = (; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
        # If0 = CCO.InterpolateC2F(; ∇0_bcs...)
        # ∂MSE_gm_∂z = ∇c(wvec(Ifx(MSE_gm)))

        @inbounds for k in TC.real_center_indices(grid)
            aux_up[i].buoy[k] = 0
            # Simple treatment for now, revise when multiple updraft closures become more well defined


            # check for instability by calculating dθ_liq_ice/dz
            if k ∉ (kc_toa, kc_surf) # don't do this for surface or top of atmosphere
                # unstable = aux_gm.θ_liq_ice[k] ≥ aux_gm.θ_liq_ice[k+1] # unstable if θ_liq_ice decreases with height
                unstable = aux_gm.θ_virt[k] ≥ aux_gm.θ_virt[k + 1] # unstable if θ_virt decreases with height [[ would start w/ MSE but idk if it's been calculated yet for env? it's not stored for GM either.. ]]
                # unstable = unstable || (∂MSE_gm_∂z[k] < FT(0)) # also unstable if ∂MSE/∂z < 0
                unstable = unstable || (MSE_gm[k] < MSE_gm[k + 1]) # also unstable if MSE increases with height
                if unstable
                    # aux_up[i].area[k+2] = TC.resolve_nan(aux_up[i].area[k+2], FT(0)) + (((k+2) != kc_toa) ? initial_profile_updraft_area : FT(0)) # socrates init area
                    aux_up[i].area[k + 1] =
                        TC.resolve_nan(aux_up[i].area[k + 1], FT(0)) +
                        (((k + 1) != kc_toa) ? initial_profile_updraft_area : FT(0)) # socrates init area
                    aux_up[i].area[k] = TC.resolve_nan(aux_up[i].area[k], FT(0)) + initial_profile_updraft_area # socrates init area
                    aux_up[i].area[k - 1] =
                        TC.resolve_nan(aux_up[i].area[k - 1], FT(0)) +
                        (((k - 1) != kc_surf) ? initial_profile_updraft_area : FT(0)) # socrates init area
                    # aux_up[i].area[k-2] = TC.resolve_nan(aux_up[i].area[k-2], FT(0)) + (((k-2) != kc_surf) ? initial_profile_updraft_area : FT(0)) # socrates init area

                    k_mf = CCO.PlusHalf(k.i)
                    k_pf = CCO.PlusHalf((k + 1).i)
                    aux_up_f[i].w[k_pf] = TC.resolve_nan(aux_up_f[i].w[k_pf], FT(0)) + FT(0.01) # small velocity perturbation
                    aux_up_f[i].w[k_mf] = TC.resolve_nan(aux_up_f[i].w[k_mf], FT(0)) + FT(0.01) # small velocity perturbation

                    k_m2f = CCO.PlusHalf((k - 1).i)
                    k_p2f = CCO.PlusHalf((k + 2).i)
                    aux_up_f[i].w[k_m2f] =
                        (k_m2f != kf_surf) ? TC.resolve_nan(aux_up_f[i].w[k_m2f], FT(0)) + FT(0.01) :
                        aux_up_f[i].w[k_m2f] # small velocity perturbation
                    aux_up_f[i].w[k_p2f] =
                        (k_p2f != kf_toa) ? TC.resolve_nan(aux_up_f[i].w[k_p2f], FT(0)) + FT(0.01) :
                        aux_up_f[i].w[k_p2f] # small velocity perturbation

                    k_m3f = CCO.PlusHalf((k - 2).i)
                    k_p3f = CCO.PlusHalf((k + 3).i)
                    aux_up_f[i].w[k_p3f] =
                        (k_p3f != kf_toa) ? TC.resolve_nan(aux_up_f[i].w[k_p3f], FT(0)) + FT(0.01) :
                        aux_up_f[i].w[k_p3f] # small velocity perturbation
                    # aux_up_f[i].w[k_m3f] = (k_m3f != kf_surf) ? TC.resolve_nan(aux_up_f[i].w[k_m3f], FT(0)) + FT(0.01) : aux_up_f[i].w[k_m3f] # small velocity perturbation

                    # Update state variables slightly to reflect buoyant updraft
                    aux_up[i].q_tot[k] = aux_gm.q_tot[k] * 1.01
                    aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k] * 1.001 # smaller bc absolute temp
                    aux_up[i].q_liq[k] = aux_gm.q_liq[k] * 1.01
                    aux_up[i].q_ice[k] = aux_gm.q_ice[k] * 1.01
                    aux_up[i].T[k] = aux_gm.T[k] * 1.002 # would need to solve for T given θ_liq_ice and q_tot increase, no?

                else
                    aux_up[i].area[k] = TC.resolve_nan(aux_up[i].area[k], FT(0)) # socrates init area
                    aux_up[i].area[k] += initial_profile_updraft_area / stable_updraft_area_reduction_factor

                    k_mf = CCO.PlusHalf(k.i)
                    k_pf = CCO.PlusHalf((k + 1).i)
                    aux_up_f[i].w[k_mf] += FT(0.01) / stable_updraft_area_reduction_factor
                    aux_up_f[i].w[k_pf] += FT(0.01) / stable_updraft_area_reduction_factor

                    # update state variables, default to grid mean [[ we could divide by stable_updraft_area_reduction_factor, we don't bc we assume the area reduction is enough ]]
                    aux_up[i].q_tot[k] = aux_gm.q_tot[k] * (1 + 0.01 / 1.0)
                    aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k] * (1 + 0.001 / 1.0)
                    aux_up[i].q_liq[k] = aux_gm.q_liq[k] * (1 + 0.01 / 1.0)
                    aux_up[i].q_ice[k] = aux_gm.q_ice[k] * (1 + 0.01 / 1.0)
                    aux_up[i].T[k] = aux_gm.T[k] * (1 + 0.002 / 1.0)
                end
            else
                aux_up[i].area[k] = FT(0)

                # update state variables, default to grid mean
                aux_up[i].q_tot[k] = aux_gm.q_tot[k]
                aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
                aux_up[i].q_liq[k] = aux_gm.q_liq[k]
                aux_up[i].q_ice[k] = aux_gm.q_ice[k]
                aux_up[i].T[k] = aux_gm.T[k]
            end

            # Do not allow initializing above max_area because strong detrainment can occur immediately, which compounded with things like MF_grad detrainment can do strange things...
            aux_up[i].area[k] = min(aux_up[i].area[k], edmf.max_area)


            # aux_up[i].area[k] = initial_profile_updraft_area # socrates init area
            # aux_up[i].q_tot[k] = aux_gm.q_tot[k]
            # aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
            # aux_up[i].q_liq[k] = aux_gm.q_liq[k]
            # aux_up[i].q_ice[k] = aux_gm.q_ice[k]
            # aux_up[i].T[k] = aux_gm.T[k]
            prog_up[i].ρarea[k] = ρ_c[k] * aux_up[i].area[k] # copied from drybubble below
            prog_up[i].ρaq_tot[k] = prog_up[i].ρarea[k] * aux_up[i].q_tot[k] # copied from drybubble below
            prog_up[i].ρaθ_liq_ice[k] = prog_up[i].ρarea[k] * aux_up[i].θ_liq_ice[k] # copied from drybubble below
            if edmf.moisture_model isa TC.NonEquilibriumMoisture
                prog_up[i].ρaq_liq[k] = prog_up[i].ρarea[k] * aux_up[i].q_liq[k] # copied from drybubble below
                prog_up[i].ρaq_ice[k] = prog_up[i].ρarea[k] * aux_up[i].q_ice[k] # copied from drybubble below
            end
        end

        @. prog_up_f[i].ρaw = ρ_f * C2F(aux_up[i].area) * aux_up_f[i].w # small ρaw perturbation


        if edmf.entr_closure isa TC.PrognosticNoisyRelaxationProcess
            @. prog_up[i].ε_nondim = 0
            @. prog_up[i].δ_nondim = 0
        end


        if edmf.surface_area_bc isa TC.FixedSurfaceAreaBC || edmf.surface_area_bc isa TC.ClosureSurfaceAreaBC
            a_surf = TC.area_surface_bc(surf, edmf, i, edmf.surface_area_bc)
            aux_up[i].area[kc_surf] = a_surf
            prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * a_surf

            # To avoid degeneracy from 0 BC on surface updraft area fraction, add small perturbation
            # to relevant prognostic variables in bottom cell -> This avoids trivial solution with no updrafts.
        elseif edmf.surface_area_bc isa TC.PrognosticSurfaceAreaBC
            # a_up_initial = max(FT(0.1), initial_profile_updraft_area)
            a_up_initial = TC.area_surface_bc(surf, edmf, i, edmf.surface_area_bc) # we jerry-rigged it to use surface_area from input parameters...
            aux_up[i].area[kc_surf] = a_up_initial
            prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * a_up_initial

            aux_up_f[i].w[kf_surf + 1] = FT(0.01) # small velocity perturbation
            prog_up_f[i].ρaw[kf_surf + 1] = ρ_f[kf_surf + 1] * a_up_initial * FT(0.01) # small ρaw perturbation

            prog_up[i].ρaq_tot[kc_surf] = ρ_c[kc_surf] * a_up_initial * aux_gm.q_tot[kc_surf]
            prog_up[i].ρaθ_liq_ice[kc_surf] = ρ_c[kc_surf] * a_up_initial * aux_gm.θ_liq_ice[kc_surf]
        end

        # set toa back to 0
        aux_up[i].area[kc_toa] = 0.0 # [ might make an unstable area gradient]
        prog_up[i].ρarea[kc_toa] = 0.0 # [ might make an unstable area gradient]

        aux_up_f[i].w[kf_toa] = 0.0 # [ no penetration boundary condition]
        prog_up_f[i].ρaw[kf_toa] = 0.0 # [ no penetration boundary condition]

    end
    return
end

import AtmosphericProfilesLibrary
const APL = AtmosphericProfilesLibrary
function initialize_updrafts_DryBubble(edmf, state, surf)
    grid = TC.Grid(state)
    # criterion 2: b>1e-4
    aux_up = TC.center_aux_updrafts(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    prog_up = TC.center_prog_updrafts(state)
    prog_up_f = TC.face_prog_updrafts(state)
    kc_surf = TC.kc_surface(grid)
    ρ_0_c = prog_gm.ρ
    ρ_0_f = aux_gm_f.ρ
    N_up = TC.n_updrafts(edmf)
    FT = TC.float_type(state)
    z_in = APL.DryBubble_updrafts_z(FT)
    z_min, z_max = first(z_in), last(z_in)
    prof_θ_liq_ice = APL.DryBubble_updrafts_θ_liq_ice(FT)
    prof_area = APL.DryBubble_updrafts_area(FT)
    prof_w = APL.DryBubble_updrafts_w(FT)
    prof_T = APL.DryBubble_updrafts_T(FT)
    face_bcs = (; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    If = CCO.InterpolateC2F(; face_bcs...)
    @inbounds for i in 1:N_up
        @inbounds for k in TC.real_face_indices(grid)
            if z_min <= grid.zf[k].z <= z_max
                aux_up_f[i].w[k] = eps(FT) # non zero value to aviod zeroing out fields in the filter
            end
        end

        @inbounds for k in TC.real_center_indices(grid)
            z = grid.zc[k].z
            if z_min <= z <= z_max
                aux_up[i].area[k] = prof_area(z)
                aux_up[i].θ_liq_ice[k] = prof_θ_liq_ice(z)
                aux_up[i].q_tot[k] = 0.0
                aux_up[i].q_liq[k] = 0.0
                aux_up[i].q_ice[k] = 0.0

                # for now temperature is provided as diagnostics from LES
                aux_up[i].T[k] = prof_T(z)
                prog_up[i].ρarea[k] = ρ_0_c[k] * aux_up[i].area[k]
                prog_up[i].ρaθ_liq_ice[k] = prog_up[i].ρarea[k] * aux_up[i].θ_liq_ice[k]
                prog_up[i].ρaq_tot[k] = prog_up[i].ρarea[k] * aux_up[i].q_tot[k]
            else
                aux_up[i].area[k] = 0.0
                aux_up[i].θ_liq_ice[k] = aux_gm.θ_liq_ice[k]
                aux_up[i].T[k] = aux_gm.T[k]
                prog_up[i].ρarea[k] = 0.0
                prog_up[i].ρaθ_liq_ice[k] = 0.0
                prog_up[i].ρaq_tot[k] = 0.0
            end
        end
        @. prog_up_f[i].ρaw = If(prog_up[i].ρarea) * aux_up_f[i].w
        a_surf = TC.area_surface_bc(surf, edmf, i, edmf.surface_area_bc)
        aux_up[i].area[kc_surf] = a_surf
        prog_up[i].ρarea[kc_surf] = ρ_0_c[kc_surf] * a_surf
    end
    return nothing
end
