
function update_cloud_frac(edmf::EDMFModel, state::State)
    # update grid-mean cloud fraction and cloud cover
    aux_bulk = center_aux_bulk(state)
    aux_gm = center_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    a_up_bulk = aux_bulk.area
    grid = Grid(state)
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        aux_gm.cloud_fraction[k] = aux_en.area[k] * aux_en.cloud_fraction[k] + a_up_bulk[k] * aux_bulk.cloud_fraction[k]
    end
end

function compute_turbconv_tendencies!(
    edmf::EDMFModel,
    state::State,
    param_set::APS,
    surf::SurfaceBase,
    Δt::Real,
    cfl_limit::Real,
    use_fallback_tendency_limiters::Bool,
)
    compute_up_tendencies!(edmf, state, param_set, Δt, use_fallback_tendency_limiters)
    compute_en_tendencies!(edmf, state, param_set, surf, Val(:tke), Val(:ρatke), Δt, cfl_limit, use_fallback_tendency_limiters)
    if edmf.thermo_covariance_model isa PrognosticThermoCovariances
        compute_en_tendencies!(edmf, state, param_set, surf, Val(:Hvar), Val(:ρaHvar), Δt, cfl_limit, use_fallback_tendency_limiters)
        compute_en_tendencies!(edmf, state, param_set, surf, Val(:QTvar), Val(:ρaQTvar), Δt, cfl_limit, use_fallback_tendency_limiters)
        compute_en_tendencies!(edmf, state, param_set, surf, Val(:HQTcov), Val(:ρaHQTcov), Δt, cfl_limit, use_fallback_tendency_limiters)
    end

    return nothing
end


"""
Apply correction to the fluxes so that the time tendencies of the convergence of the fluxes are second order accurate in time.
This is based on the second order approximation
    dq/dt ≈ -(∇·(wq)) instantaneously but
    dq/dt ≈ -(∇·(wq)) + (Δt/2) * ∇(w·∇(wq))  over a time step Δt

    Over a finite time step Δt, a second-order Taylor expansion in time gives:
        q(t + Δt) ≈ q(t) + Δt * dq/dt + (Δt^2 / 2) * d²q/dt² + O(Δt^3)

    Using dq/dt = -∇·(w q), the second derivative is:
        d²q/dt² ≈ d/dt(-∇·(w q)) ≈ -∇·(w dq/dt) ≈ -∇·(w (-∇·(w q))) = ∇·(w ∇·(w q))

    Substituting this into the Taylor expansion gives the second-order approximation for dq/dt over a time step Δt:
        dq/dt ≈ -∇·(w q) + (Δt / 2) * ∇·(w ∇·(w q)) = -∇(wq + (Δt/2) * w∇(wq)), so the flux goes from [ wq ] to  [ wq - (Δt/2) * w∇(wq) ] = [flux - (Δt/2) * w∇(flux) ]
"""
function apply_second_order_correction(flux::CC.Fields.Field, w::CC.Fields.Field, Δt::FT) where{FT}
    wvec = CC.Geometry.WVector
    # ∇c = CCO.DivergenceF2C()
    # Iff = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    Iff = CCO.InterpolateC2F(; bottom = CCO.SetGradient(wvec(zero(FT))), top = CCO.SetGradient(wvec(zero(FT)))) # we know the bc on flux is 0, but not what the gradient of the flux is. we can extrapolate the gradient or assume it's constant
    return @. flux - (Δt/FT(2)) * w * Iff(∇c(wvec(flux))) # -Δt/2 * wΔ(wq)
end

function apply_second_order_flux_correction(flux::CC.Fields.Field, w::CC.Fields.Field, Δt::FT; apply_diffusive_damping::Bool = false, correction_limit_factor::FT = FT(Inf)) where{FT}
    wvec = CC.Geometry.WVector
    # ∇c = CCO.DivergenceF2C()
    Iff = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # we know the bc on flux is 0, but not what the gradient of the flux is. we can extrapolate the gradient or assume it's constant
    # Iff = CCO.InterpolateC2F(; bottom = CCO.SetGradient(wvec(zero(FT))), top = CCO.SetGradient(wvec(zero(FT)))) # we know the bc on flux is 0, but not what the gradient of the flux is. we can extrapolate the gradient or assume it's constant
    # Iff = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    # @. flux -= (Δt/FT(2)) * w * Iff(∇c(wvec(flux))) # -Δt/2 * wΔ(wq)

    UB = CCO.UpwindBiasedProductC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # upwinding [[ I think this should be stable because it has no zero crossings within each domain...]]
    # C2Fwf = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0))) # not sure but maybe we should use this for stability as it's not biased for your upwind algorithm... no penetration. #
    w∇flux = @. UB(wvec(w), ∇c(wvec(flux))) # works, but output type is different than tendencies type, outer wvec important for getting physical units
    w∇flux_conv = (wvec.(w∇flux)).components.data.:1 # convert to field

    corrected = @. flux - clamp((Δt/FT(2)) * w∇flux_conv, -correction_limit_factor * abs(flux), correction_limit_factor * abs(flux)) # -Δt/2 * wΔ(wq)

    if apply_diffusive_damping
        # we want ν ≈ 0.1 (Δz^2/Δt). However, we don't know Δz here, so instad we use 0.1 w_max^2 . let w_max = 1 here
        ν = FT(0.1) * FT(1.0)^2 / Δt # ν ≈ 0.1 (Δz^2/Δt) with Δz = 1
        q = @. ifelse(iszero(w), FT(0), flux / w) # avoid division by zero [[ this isn't quite correct bc q could be anything... ]]
        ∇q = q # alias
        @. ∇q = ∇c(wvec(q)) # maybe just pass in q
        ∇q .= clamp.(∇q, -FT(1e3), FT(1e3)) # prevent extreme values
        @. corrected -= Iff(ν * ∇q)

        error("not implemented yet. need to pass in q...")
    end

    return corrected
end

function LBC2F_field(field::CC.Fields.Field, kc_surf::Cent{Int64};)
    LBC2F = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(field[kc_surf]))
    return LBC2F #.(field)
end
function RBC2F_field(field::CC.Fields.Field, kc_toa::Cent{Int64})
    RBC2F = CCO.RightBiasedC2F(; top = CCO.SetValue(field[kc_toa]))
    return RBC2F #.(field)
end


function compute_sgs_flux!(edmf::EDMFModel, state::State, surf::SurfaceBase, param_set::APS, Δt::Real)
    grid = Grid(state)
    N_up = n_updrafts(edmf)
    tendencies_gm = center_tendencies_grid_mean(state)
    FT = float_type(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_en = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_up = center_aux_updrafts(state)
    aux_tc = center_aux_turbconv(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_up_f = face_aux_updrafts(state)
    ρ_f = aux_gm_f.ρ
    ρ_c = prog_gm.ρ
    p_c = aux_gm.p
    kf_surf = kf_surface(grid)
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_toa = kf_top_of_atmos(grid)
    massflux = aux_tc_f.massflux
    massflux_en = aux_tc_f.massflux_en
    massflux_h = aux_tc_f.massflux_h
    massflux_qt = aux_tc_f.massflux_qt
    if edmf.moisture_model isa NonEquilibriumMoisture
        massflux_ql = aux_tc_f.massflux_ql
        massflux_qi = aux_tc_f.massflux_qi
        ql_flux_vert_adv = aux_gm_f.ql_flux_vert_adv # load for storage
        qi_flux_vert_adv = aux_gm_f.qi_flux_vert_adv # load for storage
    end


    wvec = CC.Geometry.WVector
    # ∇c = CCO.DivergenceF2C()

    # TODO: we shouldn't need to call parent here
    a_en = aux_en.area
    w_en = aux_en_f.w
    w_gm = prog_gm_f.w
    θ_liq_ice_en = aux_en.θ_liq_ice
    θ_liq_ice_gm = aux_gm.θ_liq_ice
    q_tot_gm = aux_gm.q_tot
    q_tot_en = aux_en.q_tot
    if edmf.moisture_model isa NonEquilibriumMoisture
        q_liq_en = aux_en.q_liq
        q_ice_en = aux_en.q_ice
        q_liq_gm = aux_gm.q_liq
        q_ice_gm = aux_gm.q_ice
    end
    # If = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0))) # I dont think this is right for most variables.... if we were interpolating w sure but w is already on faces. scalar vars are not 0 bounded.
    If = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # I feel like we should be using left and right biased...
    # Ic = CCO.InterpolateF2C()

    # ====================================== #
    # Using these leads to runaway updrafts...
    # ᶠinterp_a_LBF(field::CC.Fields.Field) = LBC2F_field(field, kc_surf)
    # ᶠinterp_a_RBF(field::CC.Fields.Field) = RBC2F_field(field, kc_toa)
    # IfLBF(field::CC.Fields.Field) = LBC2F_field(field, kc_surf)
    # IfRBF(field::CC.Fields.Field) = RBC2F_field(field, kc_toa)

    # above leads to runaway updrafts? [[ Note: just turning off LBF and RBF for area does stabilize the model from runaway updrafts and allow for better mixing of qt downwards]]
    # ᶠinterp_a_LBF(field::CC.Fields.Field) = ᶠinterp_a
    # ᶠinterp_a_RBF(field::CC.Fields.Field) = ᶠinterp_a
    # IfLBF(field::CC.Fields.Field) = If
    # IfRBF(field::CC.Fields.Field) = If

    # if !edmf.area_partition_model.apply_second_order_flux_correction # Allow for LBF and RBF on tracers, which better mixes downwards, but not for area, as that would leads to runaway updrafts
    #     ᶠinterp_a_LBF(field::CC.Fields.Field) = ᶠinterp_a
    #     ᶠinterp_a_RBF(field::CC.Fields.Field) = ᶠinterp_a
    #     IfLBF(field::CC.Fields.Field) = LBC2F_field(field, kc_surf)
    #     IfRBF(field::CC.Fields.Field) = RBC2F_field(field, kc_toa)

    # else # second order flux correction :: Does not need LB or RBF, since the whole point of second order is to break the first order balance
    #     # ᶠinterp_a_LBF(field::CC.Fields.Field) = ᶠinterp_a
    #     # ᶠinterp_a_RBF(field::CC.Fields.Field) = ᶠinterp_a
    #     # IfLBF(field::CC.Fields.Field) = If
    #     # IfRBF(field::CC.Fields.Field) = If

    #     # just testing
    #     ᶠinterp_a_LBF(field::CC.Fields.Field) = ᶠinterp_a
    #     ᶠinterp_a_RBF(field::CC.Fields.Field) = ᶠinterp_a
    #     IfLBF(field::CC.Fields.Field) = LBC2F_field(field, kc_surf)
    #     IfRBF(field::CC.Fields.Field) = RBC2F_field(field, kc_toa)
    # end

    ᶠinterp_a_LBF(::CC.Fields.Field) = ᶠinterp_a
    ᶠinterp_a_RBF(::CC.Fields.Field) = ᶠinterp_a
    # second order flux correction doesnt need any biasing
    # IfLBF(field::CC.Fields.Field) = edmf.area_partition_model.apply_second_order_flux_correction ? If : LBC2F_field(field, kc_surf) 
    # IfRBF(field::CC.Fields.Field) = edmf.area_partition_model.apply_second_order_flux_correction ? If : RBC2F_field(field, kc_toa)

    # # just testing
    # IfLBF(field::CC.Fields.Field) = LBC2F_field(field, kc_surf)
    # IfRBF(field::CC.Fields.Field) = RBC2F_field(field, kc_toa)

    IfLBF(::CC.Fields.Field) = If
    IfRBF(::CC.Fields.Field) = If
 


    # Compute the mass flux and associated scalar fluxes [[ note the mass flux here goes into the grid mean, but the updraft part should be pretty close to the same as compute_up_tendencies!() So the `updraft` and `env` and thus the massflux are `sgs` but not really in the classical way we'd think about turbulence ]]

    if edmf.area_partition_model.apply_second_order_flux_correction
        seclf = edmf.area_partition_model.second_order_correction_limit_factor
    end

    # ====================================== # # en is downdraft only so RBF
    ᶠinterp_a_RBF_a_en = ᶠinterp_a_RBF(a_en)
    #
    IfRBF_θ_liq_ice_en = IfRBF(θ_liq_ice_en)
    IfRBF_q_tot_en = IfRBF(q_tot_en)
    if edmf.moisture_model isa NonEquilibriumMoisture
        IfRBF_q_liq_en = IfRBF(q_liq_en)
        IfRBF_q_ice_en = IfRBF(q_ice_en)
    end
    # ====================================== #


    if edmf.area_partition_model isa StandardAreaPartitionModel

        if edmf.area_partition_model.apply_second_order_flux_correction
            massflux_en .= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en.(a_en) .* w_en, w_en, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux
            massflux_h .= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en.(a_en) .* w_en .* IfRBF_θ_liq_ice_en.(θ_liq_ice_en), w_en, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux
            massflux_qt .= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en.(a_en) .* w_en .* IfRBF_q_tot_en.(q_tot_en), w_en, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux
            if edmf.moisture_model isa NonEquilibriumMoisture
                massflux_ql .= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en.(a_en) .* (w_en) .* IfRBF_q_liq_en.(q_liq_en), w_en, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux
                massflux_qi .= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en.(a_en) .* (w_en) .* IfRBF_q_ice_en.(q_ice_en), w_en, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux

                @. ql_flux_vert_adv = massflux_ql # rn they're the same since we omit gm in the second order correction case
                @. qi_flux_vert_adv = massflux_qi # rn they're the same since we omit gm in the second order correction case
            end
        else
            # @. massflux_en = ρ_f * ᶠinterp_a(a_en) * (w_en - toscalar(w_gm))
            # @. massflux_h = ρ_f * ᶠinterp_a(a_en) * (w_en - toscalar(w_gm)) * (If(θ_liq_ice_en) - If(θ_liq_ice_gm))
            # @. massflux_qt = ρ_f * ᶠinterp_a(a_en) * (w_en - toscalar(w_gm)) * (If(q_tot_en) - If(q_tot_gm))
            massflux_en .= ρ_f .* ᶠinterp_a_RBF_a_en.(a_en) .* w_en
            @. massflux_h = massflux_en .* IfRBF_θ_liq_ice_en(θ_liq_ice_en)
            @. massflux_qt = massflux_en .* IfRBF_q_tot_en(q_tot_en)

            if edmf.moisture_model isa NonEquilibriumMoisture
                # @. massflux_ql = massflux_en * (If(q_liq_en) - If(q_liq_gm))
                # @. massflux_qi = massflux_en * (If(q_ice_en) - If(q_ice_gm))
                @. massflux_ql = massflux_en .* IfRBF_q_liq_en(q_liq_en)
                @. massflux_qi = massflux_en .* IfRBF_q_ice_en(q_ice_en)

                if param_set.user_params.use_convective_tke
                    # ql and qi should be correlated w/ tke. We handle prt of that in diffusivity but they shouldn't feel the netire advective force
                    w_conv = aux_tc_f.temporary_f1
                    @. w_conv = sqrt(max(Ifx(aux_en.tke_convective), FT(0))) # convective velocity scale
                    massflux_vert_adv_en = aux_tc_f.temporary_f2
                    @. massflux_vert_adv_en = ρ_f * ᶠinterp_a_RBF_a_en(a_en) * w_conv # feels mostly the tke scale, subsidence happens elsewhere dry
                    @. ql_flux_vert_adv = massflux_vert_adv_en .* IfRBF_q_liq_en(q_liq_en)
                    @. qi_flux_vert_adv = massflux_vert_adv_en .* IfRBF_q_ice_en(q_ice_en)
                else
                    massflux_vert_adv_en = aux_tc_f.temporary_f1
                    @. massflux_vert_adv_en = ρ_f * ᶠinterp_a_RBF_a_en(a_en) * w_en # (w_en - toscalar(w_gm)) # not sure what this should be, it's just the raw flux right?
                    @. ql_flux_vert_adv = massflux_vert_adv_en .* IfRBF_q_liq_en(q_liq_en)
                    @. qi_flux_vert_adv = massflux_vert_adv_en .* IfRBF_q_ice_en(q_ice_en)
                end
            end
        end

    elseif edmf.area_partition_model isa CoreCloakAreaPartitionModel

        aux_bulk = center_aux_bulk(state)
        aux_bulk_f = face_aux_bulk(state)

        a_cloak_up = aux_en.a_cloak_up
        a_cloak_dn = aux_en.a_cloak_dn
        a_en_remaining = aux_en.a_en_remaining
        w_cloak_up = aux_en_f.w_cloak_up
        w_cloak_dn = aux_en_f.w_cloak_dn
        #
        q_tot_cloak_up = aux_en.q_tot_cloak_up
        q_tot_cloak_dn = aux_en.q_tot_cloak_dn
        h_tot_cloak_up = aux_en.θ_liq_ice_cloak_up
        h_tot_cloak_dn = aux_en.θ_liq_ice_cloak_dn
        # The transport of these is only relevant in noneq since qt transport is also calculated
        if edmf.moisture_model isa NonEquilibriumMoisture
            q_liq_cloak_up = aux_en.q_liq_cloak_up
            q_ice_cloak_up = aux_en.q_ice_cloak_up
            q_liq_cloak_dn = aux_en.q_liq_cloak_dn
            q_ice_cloak_dn = aux_en.q_ice_cloak_dn
        end

        # q_liq_en_remaining = q_liq_en
        # q_ice_en_remaining = q_ice_en
        q_liq_en_remaining = aux_en.q_liq_en_remaining
        q_ice_en_remaining = aux_en.q_ice_en_remaining

        # zero out
        @. massflux_h = FT(0)
        @. massflux_qt = FT(0)
        if edmf.moisture_model isa NonEquilibriumMoisture
            @. ql_flux_vert_adv = FT(0) # Storing updraft part for diagnostic, prolly should rename
            @. qi_flux_vert_adv = FT(0) # Storing updraft part for diagnostic, prolly should rename
            @. massflux_ql = FT(0)
            @. massflux_qi = FT(0)
        end

         # The overall env massflux should remain unchanged [ since it balances the core updraft massflux ]
        # @. massflux_en = ρ_f * ᶠinterp_a_(a_en) * (w_en - toscalar(w_gm))
        @. massflux_en = ρ_f * ᶠinterp_a_RBF_a_en(a_en) * w_en


        # ======================================================================================================================== # # cloak up is updraft so LBF, cloak dn is downdraft so RBF, leftover env is downdraft so RBF
        ᶠinterp_a_LBF_a_cloak_up = ᶠinterp_a_LBF(a_cloak_up)
        ᶠinterp_a_RBF_a_cloak_dn = ᶠinterp_a_RBF(a_cloak_dn)
        ᶠinterp_a_RBF_a_en_remaining = ᶠinterp_a_RBF(a_en_remaining)
        #
        IfLBF_θ_liq_ice_cloak_up = IfLBF(h_tot_cloak_up)
        IfRBF_θ_liq_ice_cloak_dn = IfRBF(h_tot_cloak_dn)
        #
        IfLBF_q_tot_cloak_up = IfLBF(q_tot_cloak_up)
        IfRBF_q_tot_cloak_dn = IfRBF(q_tot_cloak_dn)
        if edmf.moisture_model isa NonEquilibriumMoisture
            IfLBF_q_liq_cloak_up = ᶠinterp_a_LBF(q_liq_cloak_up)
            IfRBF_q_liq_cloak_dn = IfRBF(q_liq_cloak_dn)
            #
            IfLBF_q_ice_cloak_up = IfLBF(q_ice_cloak_up)
            IfRBF_q_ice_cloak_dn = IfRBF(q_ice_cloak_dn)
        end
        # ======================================================================================================================== #

        # calculate fluxes
        if edmf.area_partition_model.apply_second_order_flux_correction
            massflux_h .+=
                apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_cloak_up.(a_cloak_up) .* w_cloak_up .* IfLBF_θ_liq_ice_cloak_up.(h_tot_cloak_up), w_cloak_up, Δt; correction_limit_factor = seclf) .+
                apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_cloak_dn.(a_cloak_dn) .* w_cloak_dn .* IfRBF_θ_liq_ice_cloak_dn.(h_tot_cloak_dn), w_cloak_dn, Δt; correction_limit_factor = seclf)
                if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                    massflux_h .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en_remaining.(a_en_remaining) .* w_en .* IfRBF_θ_liq_ice_en.(θ_liq_ice_en), w_en, Δt; correction_limit_factor = seclf) # add in the remaining env part if we didn't move all the downdraft to cloak
                end
            massflux_qt .+=
                apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_cloak_up.(a_cloak_up) .* w_cloak_up .* IfLBF_q_tot_cloak_up.(q_tot_cloak_up), w_cloak_up, Δt; correction_limit_factor = seclf) .+
                apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_cloak_dn.(a_cloak_dn) .* w_cloak_dn .* IfRBF_q_tot_cloak_dn.(q_tot_cloak_dn), w_cloak_dn, Δt; correction_limit_factor = seclf)
            if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                massflux_qt .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en_remaining.(a_en_remaining) .* w_en .* IfRBF_q_tot_en.(q_tot_en), w_en, Δt; correction_limit_factor = seclf) # add in the remaining env part if we didn't move all the downdraft to cloak
            end
        else
            # @. massflux_h += 
            #     ρ_f * ᶠinterp_a(a_cloak_up) * (w_cloak_up - toscalar(w_gm)) * (If(h_tot_cloak_up) - If(θ_liq_ice_gm)) + 
            #     ρ_f * ᶠinterp_a(a_cloak_dn) * (w_cloak_dn - toscalar(w_gm)) * (If(h_tot_cloak_dn) - If(θ_liq_ice_gm))
            @. massflux_h += ρ_f * ᶠinterp_a_LBF_a_cloak_up(a_cloak_up) * (w_cloak_up) * IfLBF_θ_liq_ice_cloak_up(h_tot_cloak_up) + ρ_f * ᶠinterp_a_RBF_a_cloak_dn(a_cloak_dn) * (w_cloak_dn) * IfRBF_θ_liq_ice_cloak_dn(h_tot_cloak_dn) # RB for downdraft
                if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                    # @. massflux_h += ρ_f * ᶠinterp_a(a_en_remaining) * (w_en - toscalar(w_gm)) * (If(θ_liq_ice_en) - If(θ_liq_ice_gm)) # add in the remaining env part if we didn't move all the downdraft to cloak
                    @. massflux_h += ρ_f * ᶠinterp_a_RBF_a_en_remaining(a_en_remaining) * (w_en) * IfRBF_θ_liq_ice_en(θ_liq_ice_en) # add in the remaining env part if we didn't move all the downdraft to cloak
                end
            
            # @. massflux_qt +=
                # ρ_f * ᶠinterp_a_LB_a_cloak_up(a_cloak_up) * (w_cloak_up - toscalar(w_gm)) * (IfLB_q_tot_cloak_up(q_tot_cloak_up) - IfLB_q_tot_gm(q_tot_gm)) + # LB for updraft
                # ρ_f * ᶠinterp_a_RB_a_cloak_dn(a_cloak_dn) * (w_cloak_dn - toscalar(w_gm)) * (IfRB_q_tot_cloak_dn(q_tot_cloak_dn) - IfRB_q_tot_gm(q_tot_gm)) # RB for downdraft
            @. massflux_qt += ρ_f * ᶠinterp_a_LBF_a_cloak_up(a_cloak_up) * (w_cloak_up) * IfLBF_q_tot_cloak_up(q_tot_cloak_up) + ρ_f * ᶠinterp_a_RBF_a_cloak_dn(a_cloak_dn) * (w_cloak_dn) * IfRBF_q_tot_cloak_dn(q_tot_cloak_dn) # RB for downdraft
            if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                # @. massflux_qt += ρ_f * ᶠinterp_a_RB_a_en_remaining(a_en_remaining) * (w_en - toscalar(w_gm)) * (IfRB_q_tot_en(q_tot_en) - IfRB_q_tot_gm(q_tot_gm)) # add in the remaining env part if we didn't move all the downdraft to cloak
                @. massflux_qt += ρ_f * ᶠinterp_a_RBF_a_en_remaining(a_en_remaining) * (w_en) * IfRBF_q_tot_en(q_tot_en) # add in the remaining env part if we didn't move all the downdraft to cloak
            end
        end

        if edmf.moisture_model isa NonEquilibriumMoisture
            # calculate fluxes
            if edmf.area_partition_model.apply_second_order_flux_correction
                massflux_ql .+=
                    apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_cloak_up.(a_cloak_up) .* w_cloak_up .* IfLBF_q_liq_cloak_up.(q_liq_cloak_up), w_cloak_up, Δt; correction_limit_factor = seclf) .+
                    apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_cloak_dn.(a_cloak_dn) .* w_cloak_dn .* IfRBF_q_liq_cloak_dn.(q_liq_cloak_dn), w_cloak_dn, Δt; correction_limit_factor = seclf)
                if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                    massflux_ql .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en_remaining.(a_en_remaining) .* w_en .* IfRBF_q_liq_en.(q_liq_en_remaining), w_en, Δt; correction_limit_factor = seclf) # add in the remaining env part if we didn't move
                end
                massflux_qi .+=
                    apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_cloak_up.(a_cloak_up) .* w_cloak_up .* IfLBF_q_ice_cloak_up.(q_ice_cloak_up), w_cloak_up, Δt; correction_limit_factor = seclf) .+
                    apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_cloak_dn.(a_cloak_dn) .* w_cloak_dn .* IfRBF_q_ice_cloak_dn.(q_ice_cloak_dn), w_cloak_dn, Δt; correction_limit_factor = seclf)
                if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                    massflux_qi .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_RBF_a_en_remaining.(a_en_remaining) .* w_en .* IfRBF_q_ice_en.(q_ice_en_remaining), w_en, Δt; correction_limit_factor = seclf)
                end

                @. ql_flux_vert_adv = massflux_ql # rn they're the same since we omit gm in the second order correction case
                @. qi_flux_vert_adv = massflux_qi # rn they're the same since we omit gm in the second order correction case
            else
                # @. massflux_ql += ρ_f * ᶠinterp_a(a_cloak_up) * (w_cloak_up - toscalar(w_gm)) * (If(q_liq_cloak_up) - If(q_liq_gm)) + ρ_f * ᶠinterp_a(a_cloak_dn) * (w_cloak_dn - toscalar(w_gm)) * (If(q_liq_cloak_dn) - If(q_liq_gm))
                # @. massflux_qi += ρ_f * ᶠinterp_a(a_cloak_up) * (w_cloak_up - toscalar(w_gm)) * (If(q_ice_cloak_up) - If(q_ice_gm)) + ρ_f * ᶠinterp_a(a_cloak_dn) * (w_cloak_dn - toscalar(w_gm)) * (If(q_ice_cloak_dn) - If(q_ice_gm))
                @. massflux_ql += ρ_f * ᶠinterp_a_LBF_a_cloak_up(a_cloak_up) * w_cloak_up * IfLBF_q_liq_cloak_up(q_liq_cloak_up) + ρ_f * ᶠinterp_a_RBF_a_cloak_dn(a_cloak_dn) * w_cloak_dn * IfRBF_q_liq_cloak_dn(q_liq_cloak_dn) # leave w_en outside the cloaks as 0
                @. massflux_qi += ρ_f * ᶠinterp_a_LBF_a_cloak_up(a_cloak_up) * w_cloak_up * IfLBF_q_ice_cloak_up(q_ice_cloak_up) + ρ_f * ᶠinterp_a_RBF_a_cloak_dn(a_cloak_dn) * w_cloak_dn * IfRBF_q_ice_cloak_dn(q_ice_cloak_dn)
                if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                    # @. massflux_ql += ρ_f * ᶠinterp_a(a_en_remaining) * (w_en - toscalar(w_gm)) * (If(q_liq_en) - If(q_liq_gm)) # add in the remaining env part if we didn't move all the downdraft to cloak
                    # @. massflux_qi += ρ_f * ᶠinterp_a(a_en_remaining) * (w_en - toscalar(w_gm)) * (If(q_ice_en) - If(q_ice_gm))
                    @. massflux_ql += ρ_f * ᶠinterp_a_RBF_a_en_remaining(a_en_remaining) * w_en * IfRBF_q_liq_en(q_liq_en_remaining) # add in the remaining env part if we didn't move all the downdraft to cloak
                    @. massflux_qi += ρ_f * ᶠinterp_a_RBF_a_en_remaining(a_en_remaining) * w_en * IfRBF_q_ice_en(q_ice_en_remaining)
                end

                @. ql_flux_vert_adv += ρ_f * ᶠinterp_a_LBF_a_cloak_up(a_cloak_up) * w_cloak_up * IfLBF_q_liq_cloak_up(q_liq_cloak_up) + ρ_f * ᶠinterp_a_RBF_a_cloak_dn(a_cloak_dn) * w_cloak_dn * IfRBF_q_liq_cloak_dn(q_liq_cloak_dn) # add up all the updraft cloaks, leave w_en outside the cloaks as 0
                @. qi_flux_vert_adv += ρ_f * ᶠinterp_a_LBF_a_cloak_up(a_cloak_up) * w_cloak_up * IfLBF_q_ice_cloak_up(q_ice_cloak_up) + ρ_f * ᶠinterp_a_RBF_a_cloak_dn(a_cloak_dn) * w_cloak_dn * IfRBF_q_ice_cloak_dn(q_ice_cloak_dn) # add up all the updraft cloaks, leave w_en outside the cloaks as 0
                if !edmf.area_partition_model.confine_all_downdraft_to_cloak
                    @. ql_flux_vert_adv += ρ_f * ᶠinterp_a_RBF_a_en_remaining(a_en_remaining) * w_en * IfRBF_q_liq_en(q_liq_en_remaining)
                    @. qi_flux_vert_adv += ρ_f * ᶠinterp_a_RBF_a_en_remaining(a_en_remaining) * w_en * IfRBF_q_ice_en(q_ice_en_remaining)
                end
            end

        end

    end

    # # save en massflux tendencies before adding updraft contributions
    # @. aux_tc.env_qt_tendency_vert_adv = -∇c(wvec(massflux_qt)) / ρ_c
    # @. aux_tc.env_h_tendency_vert_adv = -∇c(wvec(massflux_h)) / ρ_c
    # if edmf.moisture_model isa NonEquilibriumMoisture
    #     @. aux_tc.env_ql_tendency_vert_adv = -∇c(wvec(massflux_ql)) / ρ_c
    #     @. aux_tc.env_qi_tendency_vert_adv = -∇c(wvec(massflux_qi)) / ρ_c
    # end

    # ============================================================================================================================== #
    @inbounds for i in 1:N_up
        massflux_face_i = aux_up_f[i].massflux # is this not redundant
        parent(massflux_face_i) .= 0 # is this not redundant with the line below?
        aux_up_i = aux_up[i]
        a_up = aux_up[i].area
        w_up_i = aux_up_f[i].w
        q_tot_up = aux_up_i.q_tot
        θ_liq_ice_up = aux_up_i.θ_liq_ice
        if edmf.moisture_model isa NonEquilibriumMoisture
            q_liq_up = aux_up_i.q_liq
            q_ice_up = aux_up_i.q_ice
        end
        @. aux_up_f[i].massflux = ρ_f * ᶠinterp_a(a_up) * (w_up_i - toscalar(w_gm))
        massflux_up_i = aux_up_f[i].massflux


        # =========================================== # # updrafts so LBF only
        ᶠinterp_a_LBF_a_up = ᶠinterp_a_LBF(a_up)
        #
        IfLBF_θ_liq_ice_up = IfLBF(θ_liq_ice_up)
        IfLBF_q_tot_up = IfLBF(q_tot_up)
        #
        if edmf.moisture_model isa NonEquilibriumMoisture
            IfLBF_q_liq_up = IfLBF(q_liq_up)
            IfLBF_q_ice_up = IfLBF(q_ice_up)
        end

        # # I think doing Left bias here leads to runaway growth since nothing is there to stop you at the top
        # ᶠinterp_a_LBF_a_up = ᶠinterp_a
        # #
        # IfLBF_θ_liq_ice_up = If
        # IfLBF_q_tot_up = If
        # if edmf.moisture_model isa NonEquilibriumMoisture
        #     IfLBF_q_liq_up = If
        #     IfLBF_q_ice_up = If
        # end
        # =========================================== #

        if edmf.area_partition_model.apply_second_order_flux_correction
            massflux_h .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_up.(a_up) .* w_up_i .* IfLBF_θ_liq_ice_up.(θ_liq_ice_up), w_up_i, Δt; correction_limit_factor = seclf)
            massflux_qt .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_up.(a_up) .* w_up_i .* IfLBF_q_tot_up.(q_tot_up), w_up_i, Δt; correction_limit_factor = seclf)
        else
            # @. massflux_h += massflux_up_i * (If(θ_liq_ice_up) - If(θ_liq_ice_gm))
            # @. massflux_qt += massflux_up_i * (If(q_tot_up) - If(q_tot_gm))

            @. massflux_h += ρ_f * ᶠinterp_a_LBF_a_up(a_up) * (w_up_i) * IfLBF_θ_liq_ice_up(θ_liq_ice_up) # leave w_en outside the updraft as 0
            @. massflux_qt += ρ_f * ᶠinterp_a_LBF_a_up(a_up) * (w_up_i) * IfLBF_q_tot_up(q_tot_up) # using LB interpolation for qt to match the rest of the model, leave w_en outside the updraft as 0
        end
    
        if edmf.moisture_model isa NonEquilibriumMoisture
            q_liq_up = aux_up_i.q_liq
            q_ice_up = aux_up_i.q_ice

            if edmf.area_partition_model.apply_second_order_flux_correction
                massflux_ql .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_up.(a_up) .* w_up_i .* IfLBF_q_liq_up.(q_liq_up), w_up_i, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux
                massflux_qi .+= apply_second_order_flux_correction(ρ_f .* ᶠinterp_a_LBF_a_up.(a_up) .* w_up_i .* IfLBF_q_ice_up.(q_ice_up), w_up_i, Δt; correction_limit_factor = seclf) # note we use w_en not (w_en - w_gm) here since this is the full flux

                @. ql_flux_vert_adv = massflux_ql # rn theyre the same since we took gm out the second order correction
                @. qi_flux_vert_adv = massflux_qi # rn theyre the same since we took gm out the second order correction
            else
                # @. massflux_ql += massflux_up_i * (If(q_liq_up) - If(q_liq_gm)) # including the gm in both en and up means it should cancel out since massflux_en + massflux_bulk = 0. In compute_up_tendencies() we add 
                # @. massflux_qi += massflux_up_i * (If(q_ice_up) - If(q_ice_gm))
                @. massflux_ql += ρ_f * ᶠinterp_a_LBF_a_up(a_up) * (w_up_i) * IfLBF_q_liq_up(q_liq_up) # leave w_en outside the updraft as 0
                @. massflux_qi += ρ_f * ᶠinterp_a_LBF_a_up(a_up) * (w_up_i) * IfLBF_q_ice_up(q_ice_up) # leave w_en outside the updraft as 0

                @. ql_flux_vert_adv += massflux_up_i * IfLBF_q_liq_up(q_liq_up) # storage [[ not diff from gm -- though i feel like it should cancel out?  mean flux x gm = 0 x gm = 0]]
                @. qi_flux_vert_adv += massflux_up_i * IfLBF_q_ice_up(q_ice_up) # storage [[ not diff from gm -- though i feel like it should cancel out?  mean flux x gm = 0 x gm = 0]]
            end

        end
    end
    # ------------- #

    massflux_h[kf_surf] = 0
    massflux_qt[kf_surf] = 0
    if edmf.moisture_model isa NonEquilibriumMoisture
        massflux_ql[kf_surf] = 0
        massflux_qi[kf_surf] = 0
    end

    massflux_tendency_h = aux_tc.massflux_tendency_h
    massflux_tendency_qt = aux_tc.massflux_tendency_qt
    # Compute the  mass flux tendencies
    # Adjust the values of the grid mean variables
    # Prepare the output
    @. massflux_tendency_h = -∇c(wvec(massflux_h)) / ρ_c
    @. massflux_tendency_qt = -∇c(wvec(massflux_qt)) / ρ_c

    diffusive_flux_h = aux_tc_f.diffusive_flux_h
    diffusive_flux_qt = aux_tc_f.diffusive_flux_qt
    diffusive_flux_uₕ = aux_tc_f.diffusive_flux_uₕ

    sgs_flux_θ_liq_ice = aux_gm_f.sgs_flux_θ_liq_ice
    sgs_flux_q_tot = aux_gm_f.sgs_flux_q_tot
    sgs_flux_uₕ = aux_gm_f.sgs_flux_uₕ

    @. sgs_flux_θ_liq_ice = diffusive_flux_h + massflux_h
    @. sgs_flux_q_tot = diffusive_flux_qt + massflux_qt
    @. sgs_flux_uₕ = diffusive_flux_uₕ # + massflux_u

    # apply surface BC as SGS flux at lowest level
    lg_surf = CC.Fields.local_geometry_field(axes(ρ_f))[kf_surf]
    sgs_flux_θ_liq_ice[kf_surf] = surf.ρθ_liq_ice_flux
    sgs_flux_q_tot[kf_surf] = surf.ρq_tot_flux
    sgs_flux_uₕ[kf_surf] =
        CCG.Covariant3Vector(wvec(FT(1)), lg_surf) ⊗
        CCG.Covariant12Vector(CCG.UVVector(surf.ρu_flux, surf.ρv_flux), lg_surf)

    if edmf.moisture_model isa NonEquilibriumMoisture
        massflux_tendency_ql = aux_tc.massflux_tendency_ql
        massflux_tendency_qi = aux_tc.massflux_tendency_qi

        @. massflux_tendency_ql = -∇c(wvec(massflux_ql)) / ρ_c
        @. massflux_tendency_qi = -∇c(wvec(massflux_qi)) / ρ_c

        # my addition
        ql_tendency_vert_adv = aux_gm.ql_tendency_vert_adv # store tendency for vert adv
        qi_tendency_vert_adv = aux_gm.qi_tendency_vert_adv # store tendency for vert adv
        @. ql_tendency_vert_adv = -∇c(wvec(ql_flux_vert_adv)) / ρ_c # store tendency for vert adv
        @. qi_tendency_vert_adv = -∇c(wvec(qi_flux_vert_adv)) / ρ_c # store tendency for vert adv

        diffusive_flux_ql = aux_tc_f.diffusive_flux_ql
        diffusive_flux_qi = aux_tc_f.diffusive_flux_qi

        sgs_flux_q_liq = aux_gm_f.sgs_flux_q_liq
        sgs_flux_q_ice = aux_gm_f.sgs_flux_q_ice

        @. sgs_flux_q_liq = diffusive_flux_ql + massflux_ql
        @. sgs_flux_q_ice = diffusive_flux_qi + massflux_qi

        diffusive_flux_qr = aux_tc_f.diffusive_flux_qr # my addition
        diffusive_flux_qs = aux_tc_f.diffusive_flux_qs # my addition
        sgs_flux_q_rai = aux_gm_f.sgs_flux_q_rai # my addition
        sgs_flux_q_sno = aux_gm_f.sgs_flux_q_sno # my addition
        @. sgs_flux_q_rai = diffusive_flux_qr # we didn't store the massflux portion...
        @. sgs_flux_q_sno = diffusive_flux_qs # we didn't store the massflux portion...

        if !state.calibrate_io
            @. aux_tc.diffusive_tendency_ql = -∇c(wvec(diffusive_flux_ql)) / ρ_c
            @. aux_tc.diffusive_tendency_qi = -∇c(wvec(diffusive_flux_qi)) / ρ_c
            @. aux_tc.diffusive_tendency_qr = -∇c(wvec(diffusive_flux_qr)) / ρ_c # my addition
            @. aux_tc.diffusive_tendency_qs = -∇c(wvec(diffusive_flux_qs)) / ρ_c # my addition
            @. aux_tc.diffusive_tendency_qt = -∇c(wvec(diffusive_flux_qt)) / ρ_c
            @. aux_tc.diffusive_tendency_h = -∇c(wvec(diffusive_flux_h)) / ρ_c
        end

        sgs_flux_q_liq[kf_surf] = surf.ρq_liq_flux
        sgs_flux_q_ice[kf_surf] = surf.ρq_ice_flux
    end
    

    return nothing
end

# function smooth!(x::CC.Fields.Field, grid::Grid; window::Int=5)
#     FT = float_type(x)
#     kf_surf = kf_surface(grid)
#     kf_toa  = kf_top_of_atmos(grid)

#     half = div(window, 2)
#     buf = fill(zero(FT), window)  # circular buffer
#     buf_len = window
#     buf_start_idx = 1
#     sum_win = zero(FT)

#     # --- initialize buffer (repeat bottom value if needed) ---
#     k0 = kf_surf
#     for i in 1:window
#         k_i = k0 + i - half - 1
#         k_i = clamp(k_i, kf_surf, kf_toa)
#         val = x[k_i]
#         val = isnan(val) ? zero(FT) : val
#         buf[i] = val
#         sum_win += val
#     end

#     # --- main smoothing loop ---
#     for k in real_face_indices(grid)
#         # smoothed value
#         @info "replacing x[$(k.i)] = $(x[k]) with smoothed value $(sum_win / buf_len)"
#         x[k] = sum_win / buf_len

#         # stop at top
#         if k == kf_toa
#             break
#         end

#         # determine leaving and entering elements
#         leaving_val = buf[buf_start_idx]
#         entering_idx = k + half + 1
#         entering_idx = entering_idx > kf_toa ? kf_toa : entering_idx
#         entering_val = x[entering_idx]
#         entering_val = isnan(entering_val) ? zero(FT) : entering_val

#         # update rolling sum and buffer
#         sum_win += entering_val - leaving_val
#         buf[buf_start_idx] = entering_val
#         buf_start_idx = buf_start_idx % buf_len + 1
#     end

#     return x
# end






function compute_diffusive_fluxes(edmf::EDMFModel, state::State, surf::SurfaceBase, param_set::APS)
    grid = Grid(state)
    FT = float_type(state)
    aux_bulk = center_aux_bulk(state)
    aux_tc_f = face_aux_turbconv(state)
    aux_en_f = face_aux_environment(state)
    aux_en = center_aux_environment(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    KQ = center_aux_turbconv(state).KQ
    aeKM = center_aux_turbconv(state).temporary_1
    aeKH = center_aux_turbconv(state).temporary_2
    aeKQ = center_aux_turbconv(state).temporary_3
    prog_gm_uₕ = grid_mean_uₕ(state)
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    a_en = aux_en.area
    @. aeKM = a_en * KM
    @. aeKH = a_en * KH
    @. aeKQ = a_en * KQ
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_surf = kf_surface(grid)
    prog_gm = center_prog_grid_mean(state)
    IfKM = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKM[kc_surf]), top = CCO.SetValue(aeKM[kc_toa]))
    IfKH = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKH[kc_surf]), top = CCO.SetValue(aeKH[kc_toa]))
    IfKQ = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKQ[kc_surf]), top = CCO.SetValue(aeKQ[kc_toa]))

    @. aux_tc_f.ρ_ae_KM = IfKM(aeKM) * ρ_f
    @. aux_tc_f.ρ_ae_KH = IfKH(aeKH) * ρ_f
    @. aux_tc_f.ρ_ae_KQ = IfKQ(aeKQ) * ρ_f

    aeKQq_tot_bc = -surf.ρq_tot_flux / aux_tc_f.ρ_ae_KQ[kf_surf]
    aeKHθ_liq_ice_bc = -surf.ρθ_liq_ice_flux / aux_tc_f.ρ_ae_KH[kf_surf]
    aeKMu_bc = -surf.ρu_flux / aux_tc_f.ρ_ae_KM[kf_surf]
    aeKMv_bc = -surf.ρv_flux / aux_tc_f.ρ_ae_KM[kf_surf]

    aeKMuₕ_bc = CCG.UVVector(aeKMu_bc, aeKMv_bc)

    ∇q_tot_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_tot_bc), top = CCO.SetDivergence(FT(0)))
    ∇θ_liq_ice_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKHθ_liq_ice_bc), top = CCO.SetDivergence(FT(0)))
    # CCG.Covariant3Vector(FT(1)) ⊗ CCG.Covariant12Vector(FT(aeKMu_bc),FT(aeKMv_bc))
    local_geometry_surf = CC.Fields.local_geometry_field(axes(ρ_f))[kf_surf]
    wvec = CC.Geometry.WVector
    ∇uₕ_gm = CCO.GradientC2F(;
        bottom = CCO.SetGradient(
            CCG.Covariant3Vector(wvec(FT(1)), local_geometry_surf) ⊗
            CCG.Covariant12Vector(aeKMuₕ_bc, local_geometry_surf),
        ),
        top = CCO.SetGradient(
            CCG.Covariant3Vector(wvec(FT(0)), local_geometry_surf) ⊗ CCG.Covariant12Vector(FT(0), FT(0)),
        ),
    )


    use_separate_convective_tke_ed_coeff = (param_set.user_params.use_convective_tke && !isone(param_set.user_params.convective_tke_ed_scaling_factor))  # tke_ed_coeff    
    if use_separate_convective_tke_ed_coeff 
        #= 
            See in update_aux():
                KM[k] = c_m * ml.mixing_length * sqrt(max(aux_en.tke[k], 0))
                KH[k] = KM[k] / aux_tc.prandtl_nvec[k]
                KQ[k] = KH[k] / Le
        =#

        # IfKM_tke = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKM[kc_surf] * sqrt(max(aux_en.tke[kc_surf] - aux_en.tke_convective[kc_surf], FT(0)))), top = CCO.SetValue(aeKM[kc_toa] * sqrt(max(aux_en.tke[kc_toa] - aux_en.tke_convective[kc_toa], FT(0)))))
        # IfKH_tke = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKH[kc_surf] * sqrt(max(aux_en.tke[kc_surf] - aux_en.tke_convective[kc_surf], FT(0)))), top = CCO.SetValue(aeKH[kc_toa] * sqrt(max(aux_en.tke[kc_toa] - aux_en.tke_convective[kc_toa], FT(0)))))
        # IfKQ_tke = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKQ[kc_surf] * sqrt(max(aux_en.tke[kc_surf] - aux_en.tke_convective[kc_surf], FT(0)))), top = CCO.SetValue(aeKQ[kc_toa] * sqrt(max(aux_en.tke[kc_toa] - aux_en.tke_convective[kc_toa], FT(0)))))
        # #
        # IfKH_tke_conv = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKH[kc_surf] * sqrt(aux_en.tke_convective[kc_surf])), top = CCO.SetValue(aeKH[kc_toa] * sqrt(aux_en.tke_convective[kc_toa])))
        # IfKH_tke_conv = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKH[kc_surf] * sqrt(aux_en.tke_convective[kc_surf])), top = CCO.SetValue(aeKH[kc_toa] * sqrt(aux_en.tke_convective[kc_toa])))
        # IfKQ_tke_conv = CCO.InterpolateC2F(; bottom = CCO.SetValue(aeKQ[kc_surf] * sqrt(aux_en.tke_convective[kc_surf])), top = CCO.SetValue(aeKQ[kc_toa] * sqrt(aux_en.tke_convective[kc_toa])))

        Iftke = CCO.InterpolateC2F(; bottom = CCO.SetValue(max(aux_en.tke[kc_surf] - aux_en.tke_convective[kc_surf], FT(0))), top = CCO.SetValue(max(aux_en.tke[kc_toa] - aux_en.tke_convective[kc_toa], FT(0)))) # total tke. doesnt go to 0, but also don't extrapolate and risk neg values
        Iftke_conv = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.SetValue(FT(0))) # no convective tke at sfc since w = 0



        # ∇q_tot_en_tke_conv = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(FT(0)), top = CCO.SetDivergence(FT(0))) # no convective tke at sfc since w = 0
        # ∇θ_liq_ice_en_tke_conv = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(FT(0)), top = CCO.SetDivergence(FT(0))) # no convective tke at sfc since w = 0
        # ∇uₕ_gm_tke_conv = CCO.GradientC2F(;
        #     bottom = CCO.SetGradient(CCG.Covariant3Vector(wvec(FT(0)), local_geometry_surf) ⊗ CCG.Covariant12Vector(FT(0), FT(0))), # no convective tke at sfc since w = 0
        #     top = CCO.SetGradient(CCG.Covariant3Vector(wvec(FT(0)), local_geometry_surf) ⊗ CCG.Covariant12Vector(FT(0), FT(0))),
        # )

        ∇q_tot_en_tke_conv = ∇q_tot_en
        ∇θ_liq_ice_en_tke_conv = ∇θ_liq_ice_en
        ∇uₕ_gm_tke_conv = ∇uₕ_gm

        
        f_c_m = param_set.user_params.convective_tke_ed_scaling_factor # adjustment to c_m, tke_ed_coeff
        f_K_tke = aux_tc_f.temporary_f1
        @. f_K_tke = ifelse(Iftke(aux_en.tke) > 0, sqrt(max(Iftke(aux_en.tke) - Iftke_conv(aux_en.tke_convective), FT(0))) / sqrt(Iftke(aux_en.tke)), FT(0)) # adjustment factor for K based on convective tke removal (we leave tke in elsewhere to generate covariances... probably should be separated out also lol)
        # @. f_k_tke_conv = ifelse(Iftke(aux_en.tke) > 0, sqrt(Iftke_conv(aux_en.tke_convective)) / sqrt(Iftke(aux_en.tke)), FT(0)) # also probably should be 1 - f_K_tke? should double check..

        aux_tc = center_aux_turbconv(state)
        c_m = mixing_length_params(edmf).c_m
        Le = mixing_length_params(edmf).Le
        N = aux_tc.temporary_4
        @. N = sqrt(max(aux_tc.∂b∂z, FT(1e-6))) # buoyancy frequency, w/N is our convective overshoot mixing length, /2 for drag
        KMconv = aux_tc.temporary_5
        @. KMconv = f_c_m * c_m * ifelse(aux_en.∂MSE∂z < 0,  FT(200), FT(sqrt(2 * aux_en.tke_convective) /(2*max(N, FT(1e-6))))) * sqrt(max(aux_en.tke_convective, FT(0))) # rough convective tke eddy diffusivity for convective tke portion only, sqrt(w_tke/N)
        # KH ::KM / aux_tc.prandtl_nvec # In unstable regimes Prandtl should just be Prandtl_0 in this region... (See `turbulent_Prandtl_number()`)
        # Pr_n < 1 is unstable, we want to default to regular if ∂MSE/∂z > 0, and some enhanced value if ∂MSE/∂z < 0 :: See `/src/turbulence_functions.jl`
        Prconv = aux_tc.temporary_6
        @. Prconv = ifelse(aux_en.∂MSE∂z < 0, mixing_length_params(edmf).Pr_n * FT(1), aux_tc.prandtl_nvec) # reduce prandtl in unstable regions to enhance KH [maybe using nvec isn't right idk... coherent updrfts are still more penetrative but I guss that what's f_c_m was for]
        # KHconv ::  KMconv / Prconv
        # KQconv ::  KHconv / Le
        IfKMconv = CCO.InterpolateC2F(; bottom = CCO.SetValue(a_en[kc_surf] * KMconv[kc_surf]), top = CCO.SetValue(a_en[kc_toa] * KMconv[kc_toa]))
        IfKHconv = CCO.InterpolateC2F(; bottom = CCO.SetValue(a_en[kc_surf] * KMconv[kc_surf]/Prconv[kc_surf]), top = CCO.SetValue(a_en[kc_toa] * KMconv[kc_toa]/Prconv[kc_toa]))
        IfKQconv = CCO.InterpolateC2F(; bottom = CCO.SetValue(a_en[kc_surf] * KMconv[kc_surf]/(Prconv[kc_surf]*Le)), top = CCO.SetValue(a_en[kc_toa] * KMconv[kc_toa]/(Prconv[kc_toa]*Le)))

        f_K_tke .= 1

        #= 
            Sat efficiency for enhanced transport when we cross saturation boundaries for condensate....
            I guess this should be predicated on supersaturation and the condensate having the same gradient but oh well lol.
        =#
        peak_enhancement = FT(2.0) # max enhancement factor for saturated conditions
        sat_efficiency = aux_tc_f.temporary_f2
        @. sat_efficiency = max( Ifx(one(FT) + (peak_enhancement - one(FT)) * (FT(4) * aux_en.frac_supersat * (one(FT) - aux_en.frac_supersat))), one(FT))
    end


    if !use_separate_convective_tke_ed_coeff # regular
        @. aux_tc_f.diffusive_flux_qt = -aux_tc_f.ρ_ae_KQ * ∇q_tot_en(wvec(aux_en.q_tot))
        @. aux_tc_f.diffusive_flux_h = -aux_tc_f.ρ_ae_KH * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice))
        @. aux_tc_f.diffusive_flux_uₕ = -aux_tc_f.ρ_ae_KM * ∇uₕ_gm(prog_gm_uₕ)
    else # with convective tke adjustment [[ we do 2 way advection for stability]]

        qt_var_boost = aux_tc_f.temporary_f3
        @. qt_var_boost = clamp(Ifx(( 1 + sqrt(aux_en.QTvar)/aux_en.q_tot)/(1 - sqrt(aux_en.QTvar)/aux_en.q_tot)), one(FT), FT(3)) # bound by (1 + 0.5) / (1 - 0.5) = 3 
        @. aux_tc_f.diffusive_flux_qt = -((aux_tc_f.ρ_ae_KQ * f_K_tke)) * ∇q_tot_en(wvec(aux_en.q_tot)) -
        #  (f_c_m * qt_var_boost * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_tot_en_tke_conv(wvec(aux_en.q_tot))
            (qt_var_boost * sat_efficiency * ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_tot_en_tke_conv(wvec(aux_en.q_tot))

        h_var_boost = aux_tc_f.temporary_f3
        @. h_var_boost = clamp(Ifx((1 + sqrt(aux_en.Hvar)/aux_en.θ_liq_ice)/(1 - sqrt(aux_en.Hvar)/aux_en.θ_liq_ice)), one(FT), FT(3)) # bound by (1 + 0.5) / (1 - 0.5) = 3
        @. aux_tc_f.diffusive_flux_h = -((aux_tc_f.ρ_ae_KH * f_K_tke)) * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice)) -
        #  (f_c_m * h_var_boost * aux_tc_f.ρ_ae_KH * f_k_tke_conv) * ∇θ_liq_ice_en_tke_conv(wvec(aux_en.θ_liq_ice))
            (h_var_boost * sat_efficiency * ρ_f * IfKHconv(a_en * KMconv/Prconv)) * ∇θ_liq_ice_en_tke_conv(wvec(aux_en.θ_liq_ice))

        @. aux_tc_f.diffusive_flux_uₕ = -((aux_tc_f.ρ_ae_KM * f_K_tke)) * ∇uₕ_gm(prog_gm_uₕ) - 
            # (f_c_m * aux_tc_f.ρ_ae_KM * f_k_tke_conv) * ∇uₕ_gm_tke_conv(prog_gm_uₕ)
            (ρ_f * IfKMconv(a_en * KMconv)) * ∇uₕ_gm_tke_conv(prog_gm_uₕ)
    end

    # we could consider adding a convective tke block here with a different K... c_m by default is very conservative and probably very wrong for convective tke.... we cant fully separate the TKE's because convective tke still does generate correlations... [[ would need a KQ no convective tke too then.. i think that's just * sqrt(tke - convective tke) / sqrt(tke) ? if tke > 0 else 0 ]]

    if edmf.moisture_model isa NonEquilibriumMoisture
        aeKQq_liq_bc = FT(0)
        aeKQq_ice_bc = FT(0)

        ∇q_liq_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_liq_bc), top = CCO.SetDivergence(FT(0)))
        ∇q_ice_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_ice_bc), top = CCO.SetDivergence(FT(0)))

        if !use_separate_convective_tke_ed_coeff # regular
            @. aux_tc_f.diffusive_flux_ql = -aux_tc_f.ρ_ae_KQ * ∇q_liq_en(wvec(aux_en.q_liq))
            @. aux_tc_f.diffusive_flux_qi = -aux_tc_f.ρ_ae_KQ * ∇q_ice_en(wvec(aux_en.q_ice))
        else # with convective tke adjustment
            ∇q_liq_en_tke_conv = ∇q_liq_en # no change bc they already are using aeKQq_liq_bc = 0
            ∇q_ice_en_tke_conv = ∇q_ice_en # no change bc they already are using aeKQq_ice_bc = 0

            @. aux_tc_f.diffusive_flux_ql = -((aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_liq_en(wvec(aux_en.q_liq))) -
            #  (f_c_m * sat_efficiency * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_liq_en_tke_conv(wvec(aux_en.q_liq))
                (sat_efficiency * ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_liq_en_tke_conv(wvec(aux_en.q_liq))

            @. aux_tc_f.diffusive_flux_qi = -((aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_ice_en(wvec(aux_en.q_ice))) - 
            # (f_c_m * sat_efficiency * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_ice_en_tke_conv(wvec(aux_en.q_ice))
                (sat_efficiency * ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_ice_en_tke_conv(wvec(aux_en.q_ice))
        end
    end

    prog_pr = center_prog_precipitation(state)

    aeKQq_rai_bc = FT(0) # rain
    aeKQq_sno_bc = FT(0) # snow
    ∇q_rai_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_rai_bc), top = CCO.SetDivergence(FT(0)))
    ∇q_sno_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_sno_bc), top = CCO.SetDivergence(FT(0)))

    if !use_separate_convective_tke_ed_coeff
        @. aux_tc_f.diffusive_flux_qr = -aux_tc_f.ρ_ae_KQ * ∇q_rai_en(wvec(prog_pr.q_rai)) # add diffusive (SGS) flux for rain
        @. aux_tc_f.diffusive_flux_qs = -aux_tc_f.ρ_ae_KQ * ∇q_sno_en(wvec(prog_pr.q_sno)) # add diffusive (SGS) flux for snow
    else
        ∇q_rai_en_tke_conv = ∇q_rai_en # no change bc they already are using aeKQq_rai_bc = 0
        ∇q_sno_en_tke_conv = ∇q_sno_en # no change bc they already are using aeKQq_sno_bc = 0
        # we forgo any efficiency boosts here since downdrafts should have plenty of precip
        @. aux_tc_f.diffusive_flux_qr = -((aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_rai_en(wvec(prog_pr.q_rai))) -
        #  (f_c_m * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_rai_en_tke_conv(wvec(prog_pr.q_rai))
            (ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_rai_en_tke_conv(wvec(prog_pr.q_rai))

        @. aux_tc_f.diffusive_flux_qs = -((aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_sno_en(wvec(prog_pr.q_sno))) - 
        #  (f_c_m * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_sno_en_tke_conv(wvec(prog_pr.q_sno))
            (ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_sno_en_tke_conv(wvec(prog_pr.q_sno))
    end

    return nothing
end

function filter_small_moisture_vars(edmf::EDMFModel, state::State, param_set::APS)
    # convert small ql, qi to 0
    # this can induce small leaks in the mass/energy budget
    # if cond/evap sub/dep are a function of existing ql, qi, this could unnecessarily stunt what otherwise could be exponential growth depending on the base case
    # but the values we get sometimes of like 1e-100 are diabolical for explicit time stepping
    
    q_min = param_set.user_params.q_min

    if !iszero(q_min)
        aux_en = center_aux_environment(state)
        # prog_en = center_prog_environment(state)
        # prog_en = center_prog_environment_up_gm_version(state)


        aux_gm = center_aux_grid_mean(state)
        prog_gm = center_prog_grid_mean(state)

        prog_up = center_prog_updrafts(state)
        aux_up = center_aux_updrafts(state)
        aux_bulk = center_aux_bulk(state)
        ρ_c = prog_gm.ρ

        a_min = edmf.minimum_area

        # env
        @. aux_en.q_liq = cutoff_small_values_positive(aux_en.q_liq, q_min)
        @. aux_en.q_ice = cutoff_small_values_positive(aux_en.q_ice, q_min)

        # if edmf.moisture_model isa NonEquilibriumMoisture
        #     @. prog_en.ρaq_liq = cutoff_small_values_positive(prog_en.ρaq_liq, q_min * prog_en.ρarea)
        #     @. prog_en.ρaq_ice = cutoff_small_values_positive(prog_en.ρaq_ice, q_min * prog_en.ρarea)
        # end

        # up
        for i in 1:n_updrafts(edmf)
            if edmf.moisture_model isa NonEquilibriumMoisture
                @. prog_up[i].ρaq_liq = cutoff_small_values_positive(prog_up[i].ρaq_liq, q_min * prog_up[i].ρarea)
                @. prog_up[i].ρaq_ice = cutoff_small_values_positive(prog_up[i].ρaq_ice, q_min * prog_up[i].ρarea)
            end

            @. aux_up[i].q_liq = cutoff_small_values_positive(aux_up[i].q_liq, q_min)
            @. aux_up[i].q_ice = cutoff_small_values_positive(aux_up[i].q_ice, q_min)
        end

        @. aux_bulk.q_liq = cutoff_small_values_positive(aux_bulk.q_liq, q_min)
        @. aux_bulk.q_ice = cutoff_small_values_positive(aux_bulk.q_ice, q_min)


        # gm
        @. aux_gm.q_liq = cutoff_small_values_positive(aux_gm.q_liq, q_min)
        @. aux_gm.q_ice = cutoff_small_values_positive(aux_gm.q_ice, q_min)

        if edmf.moisture_model isa NonEquilibriumMoisture
            @. prog_gm.q_liq = cutoff_small_values_positive(prog_gm.q_liq, q_min)
            @. prog_gm.q_ice = cutoff_small_values_positive(prog_gm.q_ice, q_min)
        end

    end
    return nothing

end

function filter_small_vars(edmf::EDMFModel, state::State, param_set::APS)

    if !iszero(param_set.user_params.q_min)
        filter_small_moisture_vars(edmf, state, param_set)
    end
    return nothing
end

function filter_precipitation_vars(state::State)
    #=
        ensure positivity of q_sno, q_rai in center_prog_precipitation(state)
        I'm not sure there are any downsides lol.
    =#
    prog_pr = center_prog_precipitation(state)

    @. prog_pr.q_rai = max(prog_pr.q_rai, 0)
    @. prog_pr.q_sno = max(prog_pr.q_sno, 0)

    return nothing
end

function affect_filter!(edmf::EDMFModel, state::State, param_set::APS, surf::SurfaceBase, cfl_limit::FTT, Δt::FTT) where {FTT <: Real}
    grid = Grid(state)
    prog_en = center_prog_environment(state)
    aux_en = center_aux_environment(state)
    ###
    ### Filters
    ###
    # filter_gm_vars(edmf, state) # my addition, theyre not filtered anywhere... [[ moved to before aux_gm and prog are set in dycore.jl ]]
    filter_updraft_vars(edmf, state, surf)
    set_edmf_surface_bc(edmf, state, surf, param_set)

    filter_precipitation_vars(state)

    # filter env vars
    aux_en = center_aux_environment(state)
    prog_gm = center_prog_grid_mean(state)

    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_f = face_aux_environment(state)

    ρ_c = prog_gm.ρ
    FT = float_type(state)
    (FTT === FT) || throw(ArgumentError("Δt type $FTT does not match state float type $FT")) # hoping this helps inference...
    vmin = FT(0.01) * FT(1e-3) # 0.01 mm/s minimum velocity scale for tke floor
    @inbounds for k in real_center_indices(grid)
        # Lower bound tke. no upper bound bc env area is never 0...
        prog_en.ρatke[k] = max(prog_en.ρatke[k],   FT(0.5) * ρ_c[k] * aux_en.area[k]  * vmin^2)  # min tke based on min velocity scale [[ don't let go to 0 bc it never re-emerges]]

        prog_en.ρatke_convective[k] = max(prog_en.ρatke_convective[k], FT(0)) # min tke based on min velocity scale [[ this one is specifically convective so could go to 0, could also argue there's always some micro-convection going on lol]]

        # CFL Filter
        Δz = get_Δz(prog_gm.ρ)
        w_max = cfl_limit * min(Δz[k], Δz[k+1]) / Δt # conservative...
        ρatke_max = FT(0.5) * ρ_c[k] * aux_en.area[k] * w_max^2

        if param_set.user_params.use_convective_tke # convective tke sums with other tke in aux, so need to limit sum
            total_tke = prog_en.ρatke[k] + prog_en.ρatke_convective[k]
            if total_tke > ρatke_max
                scaling_factor = ρatke_max / total_tke
                prog_en.ρatke[k] *= scaling_factor
                prog_en.ρatke_convective[k] *= scaling_factor
            end
        else # covective tke is lumped in with regular tke so just limit individually
            prog_en.ρatke[k] = min(prog_en.ρatke[k], ρatke_max)
            prog_en.ρatke_convective[k] = min(prog_en.ρatke_convective[k], ρatke_max)
            # aux_en tke is a combination of these two so we need to really limit their sum...
        end



        # if param_set.user_params.use_convective_tke # rn this one keeps them separate but we need to add tke's so it acts on the real qt etc [[ dont do this!!!]]
        #     prog_en.ρatke[k] += prog_en.ρatke_convective[k] # this is bad because we have no way to remove so it builds up forever... i think only apply to aux...
        # end

        if edmf.thermo_covariance_model isa PrognosticThermoCovariances
            # prog_en.ρaHvar[k] = max(prog_en.ρaHvar[k], 0.0)
            # prog_en.ρaQTvar[k] = max(prog_en.ρaQTvar[k], 0.0)

            # Dont let STD exceed half of mean value
            prog_en.ρaHvar[k] = clamp(prog_en.ρaHvar[k], FT(0), (prog_en.ρaH[k] * 0.5)^2)
            prog_en.ρaQTvar[k] = clamp(prog_en.ρaQTvar[k], FT(0), (prog_en.ρaQT[k] * 0.5)^2)

            prog_en.ρaHQTcov[k] = max(prog_en.ρaHQTcov[k], -sqrt(prog_en.ρaHvar[k] * prog_en.ρaQTvar[k]))
            prog_en.ρaHQTcov[k] = min(prog_en.ρaHQTcov[k], sqrt(prog_en.ρaHvar[k] * prog_en.ρaQTvar[k]))
        end
    end

    filter_small_vars(edmf, state, param_set)

    return nothing
end


function set_edmf_surface_bc(edmf::EDMFModel, state::State, surf::SurfaceBase, param_set::APS)
    grid = Grid(state)
    thermo_params = TCP.thermodynamics_params(param_set)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    prog_en = center_prog_environment(state)
    prog_up_f = face_prog_updrafts(state)
    aux_bulk = center_aux_bulk(state)
    ts_gm = aux_gm.ts
    cp = TD.cp_m(thermo_params, ts_gm[kc_surf])
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    ρa_env_surf::FT = ρ_c[kc_surf]
    @inbounds for i in 1:N_up
        θ_surf = θ_surface_bc(surf, state, N_up, i)
        q_surf = q_surface_bc(surf, state, N_up, i)

        ρa_surf = prog_up[i].ρarea[kc_surf]
        prog_up[i].ρaθ_liq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * θ_surf
        prog_up[i].ρaq_tot[kc_surf] = prog_up[i].ρarea[kc_surf] * q_surf
        if edmf.moisture_model isa NonEquilibriumMoisture
            q_liq_surf = FT(0)
            q_ice_surf = FT(0)
            prog_up[i].ρaq_liq[kc_surf] = prog_up[i].ρarea[kc_surf] * q_liq_surf
            prog_up[i].ρaq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * q_ice_surf
        end
        prog_up_f[i].ρaw[kf_surf] = ρ_f[kf_surf] * w_surface_bc(surf)
        ρa_env_surf -= ρa_surf
    end

    flux1 = surf.ρθ_liq_ice_flux
    flux2 = surf.ρq_tot_flux
    zLL::FT = grid.zc[kc_surf].z
    ustar = surf.ustar
    oblength = surf.obukhov_length
    ρLL = prog_gm.ρ[kc_surf]
    mix_len_params = mixing_length_params(edmf)
    if edmf.thermo_covariance_model isa PrognosticThermoCovariances
        prog_en.ρaHvar[kc_surf] = ρa_env_surf * get_surface_variance(flux1 / ρLL, flux1 / ρLL, ustar, zLL, oblength)
        prog_en.ρaQTvar[kc_surf] = ρa_env_surf * get_surface_variance(flux2 / ρLL, flux2 / ρLL, ustar, zLL, oblength)
        prog_en.ρaHQTcov[kc_surf] = ρa_env_surf * get_surface_variance(flux1 / ρLL, flux2 / ρLL, ustar, zLL, oblength)
    end
    return nothing
end

function surface_helper(surf::SurfaceBase, state::State)
    FT = float_type(state)
    grid = Grid(state)
    kc_surf = kc_surface(grid)
    prog_gm = center_prog_grid_mean(state)
    zLL::FT = grid.zc[kc_surf].z
    ustar = surf.ustar
    oblength = surf.obukhov_length
    ρLL = prog_gm.ρ[kc_surf]
    return (; ustar, zLL, oblength, ρLL)
end

const ᶠinterp_a = CCO.InterpolateC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())

function area_surface_bc(surf::SurfaceBase{FT}, edmf::EDMFModel, i::Int, bc::FixedSurfaceAreaBC)::FT where {FT}
    N_up = n_updrafts(edmf)
    surface_area = length(edmf.surface_area) == N_up ? edmf.surface_area[i] : edmf.surface_area[1] / N_up # length and [1] will work on scalars or arrays, so either take the i'th value or assume it's a sum like it did before
    if (surf.bflux > 0) && iszero(surface_area)
        @warn "area_surface_bc: surface area is zero, but bflux > 0; this is likely a bug in the model; surface_area = $(surface_area); bflux = $(surf.bflux);"
    end
    return surf.bflux > 0 ? surface_area : FT(0)
end

function area_surface_bc(surf::SurfaceBase{FT}, edmf::EDMFModel, i::Int, bc::ClosureSurfaceAreaBC)::FT where {FT}
    params = bc.params
    surf_area_bc_pred = params[1] + params[2] * surf.lhf + params[3] * surf.shf
    return min(max(0, FT(surf_area_bc_pred)), edmf.max_area)
end

function area_surface_bc(surf::SurfaceBase{FT}, edmf::EDMFModel, i::Int, bc::PrognosticSurfaceAreaBC)::FT where {FT}
    return area_surface_bc(surf, edmf, i, FixedSurfaceAreaBC())  # helper for init condition to allow a custom starting area
end


function w_surface_bc(::SurfaceBase{FT})::FT where {FT}
    return FT(0)
end

function uₕ_bcs()
    return CCO.InterpolateC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
end

function θ_surface_bc(surf::SurfaceBase{FT}, state::State, N_up::Int, i::Int)::FT where {FT}
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    ρ_c = prog_gm.ρ
    grid = Grid(state)
    kc_surf = kc_surface(grid)
    ts_gm = aux_gm.ts
    UnPack.@unpack ustar, zLL, oblength, ρLL = surface_helper(surf, state)
    # N_up = n_updrafts(edmf)

    # surf.bflux > 0 || return FT(0) # this is bad I believe
    surf.bflux > 0 || return aux_gm.θ_liq_ice[kc_surf]

    if iszero(prog_up[i].ρarea[kc_surf])
        # this can happen while bflux > 0, but only with Prognostic surface area Boundary Condition [though maybe we should use a limiter on that?]
        # if it does, we need to short circuit here bc surface_scalar_coeff will be NaN from going percentile 0 to 0. Fall back to gm
        return aux_gm.θ_liq_ice[kc_surf] # fall back -- surface_scalar_coeff would be integral from percentiles 1-0 to 1 which goes to inf but idk...
    end

    a_total = sum(i -> prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf], 1:N_up)
    a_ = prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf]

    ρθ_liq_ice_flux = surf.ρθ_liq_ice_flux # assuming no ql,qi flux
    h_var = get_surface_variance(ρθ_liq_ice_flux / ρLL, ρθ_liq_ice_flux / ρLL, ustar, zLL, oblength)
    surface_scalar_coeff = percentile_bounds_mean_norm(1 - a_total + (i - 1) * a_, 1 - a_total + i * a_) # with one updraft this is from 1-a to 1
    return aux_gm.θ_liq_ice[kc_surf] + surface_scalar_coeff * sqrt(h_var)
end
function q_surface_bc(surf::SurfaceBase{FT}, state::State, N_up::Int, i::Int)::FT where {FT}
    grid = Grid(state)
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    ρ_c = prog_gm.ρ
    kc_surf = kc_surface(grid)
    # N_up = n_updrafts(edmf)

    surf.bflux > 0 || return aux_gm.q_tot[kc_surf]

    if iszero(prog_up[i].ρarea[kc_surf])
        # this can happen while bflux > 0, but only with Prognostic surface area Boundary Condition [though maybe we should use a limiter on that?]
        # if it does, we need to short circuit here bc surface_scalar_coeff will be NaN from going percentile 0 to 0. Fall back to gm
        return aux_gm.q_tot[kc_surf] # fall back -- surface_scalar_coeff would be integral from percentiles 1-0 to 1 which goes to inf but idk...
    end

    a_total = sum(i -> prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf], 1:N_up)
    a_ = prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf]

    UnPack.@unpack ustar, zLL, oblength, ρLL = surface_helper(surf, state)
    ρq_tot_flux = surf.ρq_tot_flux
    qt_var = get_surface_variance(ρq_tot_flux / ρLL, ρq_tot_flux / ρLL, ustar, zLL, oblength)
    surface_scalar_coeff = percentile_bounds_mean_norm(1 - a_total + (i - 1) * a_, 1 - a_total + i * a_)
    return aux_gm.q_tot[kc_surf] + surface_scalar_coeff * sqrt(qt_var)
end
function ql_surface_bc(surf::SurfaceBase{FT})::FT where {FT}
    return FT(0)
end
function qi_surface_bc(surf::SurfaceBase{FT})::FT where {FT}
    return FT(0)
end

function get_GMV_CoVar(
    edmf::EDMFModel,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_sym},
    ::Val{ψ_sym},
) where {covar_sym, ϕ_sym, ψ_sym}
    grid = Grid(state)
    N_up = n_updrafts(edmf)
    is_tke = covar_sym === :tke
    FT = float_type(state)
    tke_factor = is_tke ? FT(0.5) : 1
    aux_gm_c = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    prog_gm_c = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    aux_up_c = center_aux_updrafts(state)
    aux_up = is_tke ? aux_up_f : aux_up_c
    gmv_covar = getproperty(center_aux_grid_mean(state), covar_sym)
    covar_e = getproperty(center_aux_environment(state), covar_sym)
    gm = is_tke ? prog_gm_f : aux_gm_c
    to_scalar = is_tke ? toscalar : x -> x
    ϕ_gm = getproperty(gm, ϕ_sym)
    ψ_gm = getproperty(gm, ψ_sym)
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    area_en = aux_en_c.area

    Icd = is_tke ? CCO.InterpolateF2C() : x -> x
    @. gmv_covar = tke_factor * area_en * Icd(ϕ_en - to_scalar(ϕ_gm)) * Icd(ψ_en - to_scalar(ψ_gm)) + area_en * covar_e
    @inbounds for i in 1:N_up
        ϕ_up = getproperty(aux_up[i], ϕ_sym)
        ψ_up = getproperty(aux_up[i], ψ_sym)
        @. gmv_covar += tke_factor * aux_up_c[i].area * Icd(ϕ_up - to_scalar(ϕ_gm)) * Icd(ψ_up - to_scalar(ψ_gm))
    end
    return nothing
end

function compute_updraft_top(state::State, i::Int)::eltype(State)
    grid = Grid(state)
    aux_up = center_aux_updrafts(state)
    a_up = aux_up[i].area
    return z_findlast_center(k -> a_up[k] > 1e-3, grid)
end

function compute_plume_scale_height(state::State, H_up_min::FT, i::Int)::FT where {FT}
    grid = Grid(state)
    updraft_top::FT = compute_updraft_top(state, i)
    return max(updraft_top, H_up_min)
end

function compute_up_stoch_tendencies!(edmf::EDMFModel, state::State)
    N_up = n_updrafts(edmf)

    aux_up = center_aux_updrafts(state)
    tendencies_up = center_tendencies_updrafts(state)

    @inbounds for i in 1:N_up
        # prognostic entr/detr
        tends_ε_nondim = tendencies_up[i].ε_nondim
        tends_δ_nondim = tendencies_up[i].δ_nondim

        c_gen_stoch = edmf.entr_closure.c_gen_stoch
        mean_entr = aux_up[i].ε_nondim
        mean_detr = aux_up[i].δ_nondim
        ε_σ² = c_gen_stoch[1]
        δ_σ² = c_gen_stoch[2]
        ε_λ = c_gen_stoch[3]
        δ_λ = c_gen_stoch[4]
        @. tends_ε_nondim = √(2ε_λ * mean_entr * ε_σ²)
        @. tends_δ_nondim = √(2δ_λ * mean_detr * δ_σ²)
    end
end

"""
For a system
    dx1/dt = -c1 * x1 + c2 * x2 (e.g.  updraft <-- detr, entr)
    dx2/dt = c1 * x1 - c2 * x2 (e.g.   env <-- detr, entr)
and given x0 and x1 at time 0,
we calculate the analytical solution for x0 and x1 at time Δt

This actually is wrong... bc it's not area_en we're using to calculate entr... it's the area of the updraft so the ode is truly unbounded...
But useful for considering reformulations
"""
function exchange_limiter_return_x(x1::FT, x2::FT, c1::FT, c2::FT, Δt::FT) where {FT}
    k1 = (x1 + x2) / (c1 + c2)
    k2 = (x1*c1 - x2*c2) / (c1 + c2)

    x1_new = k1 * c2 + k2 * exp(-(c1 + c2) * Δt)
    x2_new = k1 * c1 - k2 * exp(-(c1 + c2) * Δt)

    if isnan(x1_new)
        x1_new = x1
    end
    if isnan(x2_new)
        x2_new = x2
    end

    return x1_new, x2_new
end

"""
Same as exchange_limiter() but returns Δx1 == -Δx2
"""
function exchange_return_dxdt(x1::FT, x2::FT, c1::FT, c2::FT, Δt::FT) where {FT}
    k1 = (x1 + x2) / (c1 + c2)
    k2 = (x1*c1 - x2*c2) / (c1 + c2)

    x1_new = k1 * c2 + k2 * exp(-(c1 + c2) * Δt)
    return isnan(x1_new) ? FT(0) : (x1_new - x1) / Δt
end

# """
# Should this take into account the advective part? I say no for rn...
# """
# function limit_entr_detr(entr_plus_detr::FT, x_up::FT, x_en::FT, Δt::FT ) where {FT}
#     return safe_clamp(entr_plus_detr, -x_up / Δt, x_en / Δt) # this could really just be clamp here I guess...
# end


"""
Stop these ridiculous upgradient fluxes from using separate equations for area and the prognostic variable
    - makes sure the resulting pure prognostic variable (without ρa) comes out between itself and the original grid mean value...

    Based on the idea that things don't gotta be balanced, but you gotta be approaching the mean...

    TODO: Deprecate this !!!

        It can't really work well because other tendencies will de-rail it. So if you boost your tendency to say, send you to grid mean, but you have other tendencies, you're going to be in an undefined spot.
        In a true model, the entr/detr balance might actually counteract these other tendencies (e.g. something gets advected in and mixed out) so turning that off could do arbitrary things.

        E.g. say whenever w = 0, you mix θ all the way to grid mean value. Say that's an increase in θ. And say below you there's still an updraft w>0: then you have also updraft flux from below. Well now you're above the grid mean value..., positively buoyant, etc you can see how this can go wrong. This could even chain react and continue up to the top of the domain!!
        This will only work out if this is applied post-summing up all the other tendencies. So you'd have to do it all the way at the end....

        Also the grid mean values can move anyway via our large-scale forcings, etc...

        At the very least, we'd need to change it to wrap all the other tendencies in ρaprogvar_tendency... and then you wouldn't want to apply it without still adding relevant tendencies to the gridmean_tendency so something like advection can still impact the grid mean.

        So it would be something like 
            - calculate ρaprogvar_tendency tendency w/o entr/detr, use it to calculate a new ρaprogvar
            - calculate ρgridmean_progvar_tendency w/o entr/detr, use it to get new gridmean_progvar value at end of timestep

            - calculate the new area ρarea. Note this could still be decoupled from ρaprogvar, we don't have a fix for that...
            - with the new area, calculate the new progbar

            Now, we assume entr/detr will mix to somewhere between progvar and gridmean_progvar
            - proceed w/ the sol'n as before

        This is very challenging to do in practice, and you'd need to pass in something like `ρaprogvar_tendencies` and `ρgridmean_progvar_tendencies` to the function to do this. So it's not really worth it for now.
"""

function progvar_area_tendency_resolver(
    limiter::AbstractTendencyLimiter,
    ρarea::FT,
    progvar::FT,
    gridmean_progvar::FT,
    ρarea_tendency::FT,
    ρaprogvar_tendency::FT, # e.g. the entr/detr tendency.
    ρ::FT,
    Δt::FT,
    use_tendency_resolver::Bool = true,
    do_now::Bool = true,
    other_ρaprogvar_tendency::FT = FT(0), # these are assumed to affect both the progvar and the gridmean_progvar [think like external advection]
    # gridmean_progvar_tendency::FT, # we can't really know this in TC, it has too many other things in it... [large-scale forcing etc, env microphysics etc]
) where {FT}

    # if iszero(ρarea) # if no area, use the grid mean value bc the progvar value came from bad information potentially. [deprecated in favor of forcing the user to fix -- you couuld hav put in aux for example...]
    #     progvar = gridmean_progvar
    # end


    if !do_now # we're not calculating the tendency now, just return 0. In the code below we use this to choose when to calculate the tendency... sort of a hack I know.
        if use_tendency_resolver
            return FT(0) # we are using the resolver, but we aren't doing it now, so just return 0.
        else
            ; # if we're not using the resolver, there's no point in deferring to a later by returning FT(0), that later call will never come bc use_tendency_resolver will still be false.
        end
    end

    if !use_tendency_resolver # to make it easier to toggle this on/off w/o recompiling
        return ρaprogvar_tendency # test just short circuiting
    end
        
    new_ρa = clamp(ρarea + ρarea_tendency * Δt,  FT(0), FT(ρ)) # only w/ explicit timesteps...
    ρarea_tendency = (new_ρa - ρarea) / Δt # update after clamping.

    # calculate how other tendencies for the prognostic variable should affect it and the grid mean..
    if !iszero(other_ρaprogvar_tendency)
        new_ρaprogvar = max(ρarea*progvar + other_ρaprogvar_tendency * Δt, FT(0))
        new_ρgridmean_progvar = max(ρ*gridmean_progvar + other_ρaprogvar_tendency * Δt, FT(0)) # assume positive definite

        # now entr/detr can work on these new values...
        progvar = new_ρaprogvar / new_ρa
        gridmean_progvar = new_ρgridmean_progvar / ρ
        ρaprogvar = new_ρaprogvar

        # could we instead use exchange_return_dxdt() in this case to summarily change calculate completely new tendency values?
        val_1 = (new_ρa*gridmean_progvar - ρaprogvar) / Δt # going to grid mean
        val_2 = FT(0) # the mixing does no work too move towards grid mean.
    else
        ρaprogvar = progvar * ρarea
        val_1 = (new_ρa*gridmean_progvar - ρaprogvar) / Δt # whatever direction this is, we must not pass it, it could be up or down
        val_2 = ρarea_tendency * progvar # this is no change. We also must not pass this and go the wrong direction. [ same as (new_ρa*progvar - ρaprogvar) / Δt ]
    end
    # ρgridmean_progvar = gridmean_progvar * ρ

    #=
    Updraft raw var is greater than grid mean, so after mixing assuming no upgradient fluxes, the raw value should be greater than the current grid mean
        Thus: ρaprogvar > (new_ρa * gridmean_progvar)
    Note ρarea_tendency could be negative so we don't actually know the sign of the tendency, so we use safe_clamp

        - really, no matter the direction of the area tendency, the raw value should approach the grid mean
        - so increasing area
    =#

    # in case we're using a more conservative limiter
    if val_1 > FT(0) # positive tendency, (new_ρa*gridmean_progvar - ρaprogvar) > 0, don't deplete it
        val_1 = -limit_tendency(limiter, -val_1, new_ρa*gridmean_progvar - ρaprogvar, Δt) # in case we're using a more conservative limiter
    else # negative tendency, so (ρaprogvar - new_ρa*gridmean_progvar) > 0, don't deplete it
        val_1 = limit_tendency(limiter, val_1, ρaprogvar - new_ρa*gridmean_progvar, Δt) # in case we're using a more conservative limiter
    end

    ρaprogvar_tendency = safe_clamp(ρaprogvar_tendency, val_1, val_2) # this is the tendency we want to limit

    return ρaprogvar_tendency
end




# FT version to allow dispatching w/ a CC Field
function progvar_area_tendency_resolver(
    limiter::AbstractTendencyLimiter,
    ρarea::FT,
    progvar::FT,
    gridmean_progvar::FT,
    ρarea_tendency::FT,
    ρaprogvar_tendency::FT, # e.g. the entr/detr tendency.
    ρ::FT,
    Δt::FT,
    use_tendency_resolver::FT, # can't have a default bc the other version already has that method defined.
    do_now::FT, # whehter to perform the operation now or later
    other_ρaprogvar_tendency::FT = FT(0), # these are assumed to affect both the progvar and the gridmean_progvar [think like external advection]
) where {FT}

    use_tendency_resolver = if iszero(use_tendency_resolver)
        false
    elseif isone(use_tendency_resolver)
        true
    else
        error("progvar_area_tendency_resolver: integer use_tendency_resolver must be 1 or 0 to match boolean true or false, not $use_tendency_resolver")
    end

    do_now = if iszero(do_now)
        false
    elseif isone(do_now)
        true
    else
        error("progvar_area_tendency_resolver: integer do_now must be 1 or 0 to match boolean true or false, not $do_now")
    end

    return progvar_area_tendency_resolver(limiter, ρarea, progvar, gridmean_progvar, ρarea_tendency, ρaprogvar_tendency, ρ, Δt, use_tendency_resolver, do_now, other_ρaprogvar_tendency)
end


function compute_up_tendencies!(edmf::EDMFModel, state::State, param_set::APS, Δt::Real, use_fallback_tendency_limiters::Bool)
    N_up = n_updrafts(edmf)
    grid = Grid(state)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = float_type(state)

    aux_up = center_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    aux_bulk = center_aux_bulk(state)
    aux_en_f = face_aux_environment(state)
    aux_up_f = face_aux_updrafts(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    tendencies_en = center_tendencies_environment(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    aux_tc_f = face_aux_turbconv(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    a_max = edmf.max_area
    # edtl = edmf.tendency_limiters.entr_detr_tendency_limiter
    edtl = get_tendency_limiter(edmf.tendency_limiters, Val(:entr_detr), use_fallback_tendency_limiters)

    # use resolver
    use_tendency_resolver = edmf.tendency_limiters.use_tendency_resolver
    use_tendency_resolver_on_full_tendencies = edmf.tendency_limiters.use_tendency_resolver_on_full_tendencies


    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        @. aux_up_i.entr_turb_dyn = aux_up_i.entr_sc + aux_up_i.entr_ml + aux_up_i.frac_turb_entr
        @. aux_up_i.detr_turb_dyn = aux_up_i.detr_sc + aux_up_i.detr_ml + aux_up_i.frac_turb_entr
    end

    UB = CCO.UpwindBiasedProductC2F(bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    # Ic = CCO.InterpolateF2C()

    wvec = CC.Geometry.WVector
    # ∇c = CCO.DivergenceF2C()
    w_bcs = (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    LBF = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))

    # We know that, since W = 0 at z = 0, BCs for entr, detr,
    # and buoyancy should not matter in the end
    zero_bcs = (; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))
    I0f = CCO.InterpolateC2F(; zero_bcs...)
    adv_bcs = (; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0))))
    LBC = CCO.LeftBiasedF2C(; bottom = CCO.SetValue(FT(0)))
    ∇f = CCO.DivergenceC2F(; adv_bcs...)
    a_en = aux_en.area

    # Solve for updraft area fraction
    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        w_up = aux_up_f[i].w
        a_up = aux_up_i.area
        q_tot_up = aux_up_i.q_tot
        q_tot_en = aux_en.q_tot
        θ_liq_ice_en = aux_en.θ_liq_ice
        θ_liq_ice_up = aux_up_i.θ_liq_ice
        entr_turb_dyn = aux_up_i.entr_turb_dyn
        detr_turb_dyn = aux_up_i.detr_turb_dyn
        entr_rate_inv_s = aux_up_i.entr_rate_inv_s
        detr_rate_inv_s = aux_up_i.detr_rate_inv_s
        θ_liq_ice_tendency_precip_formation = aux_up_i.θ_liq_ice_tendency_precip_formation
        qt_tendency_precip_formation = aux_up_i.qt_tendency_precip_formation

        ρarea = prog_up[i].ρarea
        ρaθ_liq_ice = prog_up[i].ρaθ_liq_ice
        ρaq_tot = prog_up[i].ρaq_tot

        tends_ρarea = tendencies_up[i].ρarea
        tends_ρaθ_liq_ice = tendencies_up[i].ρaθ_liq_ice
        tends_ρaq_tot = tendencies_up[i].ρaq_tot
        area_en = aux_en.area


        #=
            We tried some convective-tke based adjustments to entr/detr rate but they didn't work at all.
            Instead we now just at the end directly adjoin to the existing updraft.
        =#

        


        # advection, precipitation, and sedimentation tendencies are somewaht constrained at construction, but the overall system can still be unstable.
        # one method would be to apply clamp and to normalize the tendencies by their individual contributions to the total tendency... would take a lot more tracking though...



        # store area tendency due to entrainment and detrainment bc that is the only part that really affects the limiters later on that push us towards grid mean. the advection part really isn't relevant to that -- advection is handled outside.
        # For the same reason, we also don't know the true new area so it's hard to say the limiter will work but only that in the area part from entr/detr we are pushed towards grid mean.
        
        #= 
            The logic should be detrainment can detrain positive advection have priority over negative advection, and entrainment can entrain positive advection have priority over negative advection.
            so in the net, we should allow net detrainment to remove all incoming advection, but w/ net entrainment, the bound is harder to define. We can entrain to cancel all negative advection, and then up to the limit..

            We of course can't account for the other tendency parts, like ql_tendency_noneq... really you'd want to limit sinks together... since entr/detr is the last tendency, maybe you'd just limit it only? idk...
        =#

        ql_tendency_sedimentation = aux_up_i.ql_tendency_sedimentation
        qi_tendency_sedimentation = aux_up_i.qi_tendency_sedimentation
        qt_tendency_sedimentation = aux_up_i.qt_tendency_sedimentation # this is just the sum of the liq and ice tendencies
        θ_liq_ice_tendency_sedimentation = aux_up_i.θ_liq_ice_tendency_sedimentation

        tends_advec = aux_tc.temporary_1 # Reuse temp [[ no nested calls so ok ]]
        tends_other = aux_tc.temporary_2 # Reuse temp [[ no nested calls so ok ]]

        # thermo_params = TCP.thermodynamics_params(param_set)
        # Π = TD.exner.(thermo_params, aux_gm.ts)

        # f_q = similar(ρarea)
        # if edmf.area_partition_model isa CoreCloakAreaPartitionModel # Ignore w for now
        #     # As with advection, we will again partition everything outside of en_remaining into the updrafts and cloak up, leaving nothing in cloak_dn, and en_remainig unchanged
        #     # Boost factors for ice formation segregation to updraft.
        #     @inbounds for k in real_center_indices(grid)
        #         combined_up_area = aux_bulk.area[k] + aux_en.a_cloak_up[k]
        #         f_q[k] = (combined_up_area > sqrt(eps(FT))) ? ((one(FT) - aux_en.a_en_remaining[k]) / combined_up_area) : one(FT)
        #     end
        # else
        #     @. f_q = one(FT)
        # end

        scalar_nudge_τᵣ = 20*60 # testing hardcoded to see if its' good before going all in [ would neet to pass forcing object to this fcn]

        if edmf.entrainment_type isa FractionalEntrModel
            @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρarea)))
            @. tends_ρarea = 
                tends_advec + 
                limit_tendency(edtl, ρarea * Ic(w_up) * (entr_turb_dyn - detr_turb_dyn), ρarea + max(tends_advec,0), ρ_c * a_max - ρarea + max(-tends_advec,0), Δt)
        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρarea)))
            @. tends_ρarea = 
                tends_advec + 
                limit_tendency(edtl, ρarea * (entr_rate_inv_s - detr_rate_inv_s), ρarea + max(tends_advec,0), ρ_c * a_max - ρarea + max(-tends_advec,0), Δt)
        end

        if edmf.entrainment_type isa FractionalEntrModel
            @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaθ_liq_ice)))
            @. tends_other = (ρ_c * θ_liq_ice_tendency_precip_formation) + (ρ_c * θ_liq_ice_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
            @. tends_ρaθ_liq_ice =
                # -∇c(wvec(LBF(Ic(w_up) * ρaθ_liq_ice))) + 
                tends_advec + 
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaθ_liq_ice / ρarea, prog_gm.ρθ_liq_ice / ρ_c), prog_gm.ρθ_liq_ice / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * Ic(w_up) * entr_turb_dyn * θ_liq_ice_en) - (ρaθ_liq_ice * Ic(w_up) * detr_turb_dyn) , ρaθ_liq_ice + max(tends_advec+tends_other,0), ρ_c * area_en * θ_liq_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                    (ρ_c * θ_liq_ice_tendency_precip_formation) # + 
                    # f_q*max(ρarea * aux_tc.θ_liq_ice_tendency_precip_sinks, 0) # + # qi_sub_dep contribution [[ i think this is needed so significant snow formation doesn't kill updrafts?]]
                    #+ ρarea*aux_up[i].dTdt/Π # test adding nudging in so that nudging doesnt prevent env and updraft from mixing and converging...

            @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaq_tot)))
            @. tends_other = (ρ_c * qt_tendency_precip_formation) + (ρ_c * qt_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
            @. tends_ρaq_tot =
                # -∇c(wvec(LBF(Ic(w_up) * ρaq_tot))) +
                tends_advec +
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_tot / ρarea, prog_gm.ρq_tot / ρ_c), prog_gm.ρq_tot / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * Ic(w_up) * entr_turb_dyn * q_tot_en) - (ρaq_tot * Ic(w_up) * detr_turb_dyn), ρaq_tot + max(tends_advec+tends_other,0), ρ_c * area_en * q_tot_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                    (ρ_c * qt_tendency_precip_formation) 

        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaθ_liq_ice)))
            @. tends_other = (ρ_c * θ_liq_ice_tendency_precip_formation) + (ρ_c * θ_liq_ice_tendency_sedimentation)# external tendencies... we will limit entr/detr as the final tendency to
            @. tends_ρaθ_liq_ice =
                # -∇c(wvec(LBF(Ic(w_up) * ρaθ_liq_ice))) +
                tends_advec +
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaθ_liq_ice / ρarea, prog_gm.ρθ_liq_ice / ρ_c), prog_gm.ρθ_liq_ice / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * entr_rate_inv_s * θ_liq_ice_en) - (ρarea * detr_rate_inv_s * θ_liq_ice_up), ρaθ_liq_ice+max(tends_advec+tends_other,0), ρ_c * area_en * θ_liq_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                (ρ_c * θ_liq_ice_tendency_precip_formation) +
                # f_q*max(ρarea * aux_tc.θ_liq_ice_tendency_precip_sinks, 0) # + # qi_sub_dep contribution [[ i think this is needed so significant snow formation doesn't kill updrafts?]]
                # ifelse(param_set.user_params.apply_nudging_to_updraft && (ρarea/ρ_c > 1.), ρarea * aux_up[i].dTdt/TD.exner(thermo_params, aux_gm.ts), FT(0)) # This version breaks because dTdt contains extraneous stuff and because we're not restoring directly to gm, we can have unresolved drift that simply gets compensated in env.
                ifelse(param_set.user_params.apply_nudging_to_updraft && (ρarea/ρ_c > FT(0.01)), ρarea * (aux_gm.H_nudge - aux_up[i].θ_liq_ice)/scalar_nudge_τᵣ, FT(0)) # dTdt includes hadv and stuff... ... nudging can kill the updraft though...

            @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaq_tot)))
            @. tends_other = (ρ_c * qt_tendency_precip_formation) + (ρ_c * qt_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
            @. tends_ρaq_tot =
                # -∇c(wvec(LBF(Ic(w_up) * ρaq_tot))) + 
                tends_advec +
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_tot / ρarea, prog_gm.ρq_tot / ρ_c), prog_gm.ρq_tot / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * entr_rate_inv_s * q_tot_en) - (ρarea * detr_rate_inv_s * q_tot_up), ρaq_tot + max(tends_advec+tends_other,0), ρ_c * area_en * q_tot_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                    (ρ_c * qt_tendency_precip_formation) +
                    # ifelse(param_set.user_params.apply_nudging_to_updraft, ρarea * aux_gm.dqtdt_nudge, FT(0)) # dqvdt includes hadv and stuff... ... nudging can kill the updraft though... [[ Tghis version has has no restoring force... so it can cause drift... as the updraft area changes... we'd have to nudge updraft and env separate
                    ifelse(param_set.user_params.apply_nudging_to_updraft && (ρarea/ρ_c > FT(0.01)), ρarea * (aux_gm.qt_nudge - aux_up[i].q_tot)/scalar_nudge_τᵣ, FT(0)) # dqvdt includes hadv and stuff... ... nudging can kill the updraft though... [[ This version calculates nudging separately so no drift... ]]
        end

        if edmf.moisture_model isa NonEquilibriumMoisture

            q_liq_up = aux_up_i.q_liq
            q_ice_up = aux_up_i.q_ice
            q_liq_en = aux_en.q_liq
            q_ice_en = aux_en.q_ice

            ql_tendency_noneq = aux_up_i.ql_tendency_noneq
            qi_tendency_noneq = aux_up_i.qi_tendency_noneq
            ql_tendency_precip_formation = aux_up_i.ql_tendency_precip_formation
            qi_tendency_precip_formation = aux_up_i.qi_tendency_precip_formation

            ρaq_liq = prog_up[i].ρaq_liq
            ρaq_ice = prog_up[i].ρaq_ice

            tends_ρaq_liq = tendencies_up[i].ρaq_liq
            tends_ρaq_ice = tendencies_up[i].ρaq_ice




            # TODO: Double check if it's better to use:: limit_tendency(edtl, (ρarea * Ic(w_up) * entr_turb_dyn * q_liq_en) - (ρaq_liq * Ic(w_up) * detr_turb_dyn) , max(ρaq_liq + tends_advec+tends_other,0), max(ρ_c * area_en * q_liq_en + -(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +

    

            if edmf.entrainment_type isa FractionalEntrModel
                
                @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaq_liq)))
                @. tends_other = (ρ_c * ql_tendency_precip_formation + ql_tendency_noneq) + (ρ_c * ql_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_liq =
                    # -∇c(wvec(LBF(Ic(w_up) * ρaq_liq))) + 
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        limit_tendency(edtl, (ρarea * Ic(w_up) * entr_turb_dyn * q_liq_en) - (ρaq_liq * Ic(w_up) * detr_turb_dyn) , ρaq_liq + max(tends_advec+tends_other,0), ρ_c * area_en * q_liq_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (ql_tendency_precip_formation + ql_tendency_noneq))
                        
                @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaq_ice)))
                @. tends_other = (ρ_c * qi_tendency_precip_formation + qi_tendency_noneq) + (ρ_c * qi_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_ice =
                    # -∇c(wvec(LBF(Ic(w_up) * ρaq_ice))) +
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_ice / ρarea, prog_gm.q_ice), prog_gm.q_ice, tends_ρarea,
                        limit_tendency(edtl, (ρarea * Ic(w_up) * entr_turb_dyn * q_ice_en) - (ρaq_ice * Ic(w_up) * detr_turb_dyn) , ρaq_ice + max(tends_advec+tends_other,0), ρ_c * area_en * q_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (qi_tendency_precip_formation + qi_tendency_noneq))
                    
            elseif edmf.entrainment_type isa TotalRateEntrModel
                @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaq_liq)))
                @. tends_other = (ρ_c * ql_tendency_precip_formation + ql_tendency_noneq) + (ρ_c * ql_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_liq =
                    # -∇c(wvec(LBF(Ic(w_up) * ρaq_liq))) + 
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        limit_tendency(edtl, (ρarea * entr_rate_inv_s * q_liq_en) - (ρarea * detr_rate_inv_s * q_liq_up), ρaq_liq + max(tends_advec+tends_other, 0), ρ_c * area_en * q_liq_en + max(-(tends_advec+tends_other), 0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (ql_tendency_precip_formation + ql_tendency_noneq))

                @. tends_advec = -∇c(wvec(LBF(Ic(w_up) * ρaq_ice)))
                @. tends_other = (ρ_c * qi_tendency_precip_formation + qi_tendency_noneq) + (ρ_c * qi_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_ice =
                    # -∇c(wvec(LBF(Ic(w_up) * ρaq_ice))) +
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_ice / ρarea, prog_gm.q_ice), prog_gm.q_ice, tends_ρarea,
                        limit_tendency(edtl, (ρarea * entr_rate_inv_s * q_ice_en) - (ρarea * detr_rate_inv_s * q_ice_up), ρaq_ice + max(tends_advec+tends_other,0), ρ_c * area_en * q_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (qi_tendency_precip_formation + qi_tendency_noneq))
            end


            if edmf.cloud_sedimentation_model isa CloudSedimentationModel # &&  !edmf.cloud_sedimentation_model.grid_mean        
                ql_tendency_sedimentation = aux_up_i.ql_tendency_sedimentation
                qi_tendency_sedimentation = aux_up_i.qi_tendency_sedimentation
                qt_tendency_sedimentation = aux_up_i.qt_tendency_sedimentation # this is just the sum of the liq and ice tendencies
                θ_liq_ice_tendency_sedimentation = aux_up_i.θ_liq_ice_tendency_sedimentation
                # ==  these get backed out after gm tendencies applied in dycore so only add if using in calculation of entr/detr == #
                # ql_tendency_sedimentation_en = aux_en.ql_tendency_sedimentation 
                # qi_tendency_sedimentation_en = aux_en.qi_tendency_sedimentation
                # qt_tendency_sedimentation_en = aux_en.qt_tendency_sedimentation # this is just the sum of the liq and ice tendencies
                # θ_liq_ice_tendency_sedimentation_en = aux_en.θ_liq_ice_tendency_sedimentation
                # ================================================================================================ #

                @. tends_ρaq_liq += (ρ_c * ql_tendency_sedimentation) # + # rolled vertical velocity into the sedimentation tendency (undid)in 
                @. tends_ρaq_ice += (ρ_c * qi_tendency_sedimentation) # + # rolled vertical velocity into the sedimentation tendency (undid)
                
                @. tends_ρaθ_liq_ice += (ρ_c * θ_liq_ice_tendency_sedimentation) # + # rolled vertical velocity into the sedimentation tendency (undid)
                @. tends_ρaq_tot += (ρ_c * qt_tendency_sedimentation) # + # rolled vertical velocity into the sedimentation tendency (undid)

            end

            tends_ρaq_liq[kc_surf] = 0
            tends_ρaq_ice[kc_surf] = 0




        else

            if edmf.cloud_sedimentation_model isa CloudSedimentationModel # &&  !edmf.cloud_sedimentation_model.grid_mean        
                qt_tendency_sedimentation = aux_up_i.qt_tendency_sedimentation # this is just the sum of the liq and ice tendencies
                θ_liq_ice_tendency_sedimentation = aux_up_i.θ_liq_ice_tendency_sedimentation
                # ==  these get backed out after gm tendencies applied in dycore so only add if using in calculation of entr/detr == #
                # ql_tendency_sedimentation_en = aux_en.ql_tendency_sedimentation 
                # qi_tendency_sedimentation_en = aux_en.qi_tendency_sedimentation
                # qt_tendency_sedimentation_en = aux_en.qt_tendency_sedimentation # this is just the sum of the liq and ice tendencies
                # θ_liq_ice_tendency_sedimentation_en = aux_en.θ_liq_ice_tendency_sedimentation
                # ================================================================================================ #                
                @. tends_ρaθ_liq_ice += (ρ_c * θ_liq_ice_tendency_sedimentation) # + # rolled vertical velocity into the sedimentation tendency (undid)
                @. tends_ρaq_tot += (ρ_c * qt_tendency_sedimentation) # + # rolled vertical velocity into the sedimentation tendency (undid)
            end
        end

        # we cannot bound the overall tendency , we are not bound to go towards the mean overall.... just only from entr/detr !!!!

        # prognostic entr/detr
        if edmf.entr_closure isa PrognosticNoisyRelaxationProcess
            c_gen_stoch = edmf.entr_closure.c_gen_stoch
            mean_entr = aux_up[i].ε_nondim
            mean_detr = aux_up[i].δ_nondim
            ε_λ = c_gen_stoch[3]
            δ_λ = c_gen_stoch[4]
            tends_ε_nondim = tendencies_up[i].ε_nondim
            tends_δ_nondim = tendencies_up[i].δ_nondim
            ε_nondim = prog_up[i].ε_nondim
            δ_nondim = prog_up[i].δ_nondim
            @. tends_ε_nondim = ε_λ * (mean_entr - ε_nondim)
            @. tends_δ_nondim = δ_λ * (mean_detr - δ_nondim)
        end

        tends_ρaθ_liq_ice[kc_surf] = 0
        tends_ρaq_tot[kc_surf] = 0


        # How do we account for tendencies_gm = center_tendencies_grid_mean(state) here?
        #=
        area - fine
        θ_liq_ice, q_tot, q_liq, q_ice -- you can take up the entire environment + grid mean tendency * Δt, or you can lose yourself...
        ... However, the gm_tendencies in ∑_tendencies are calculated AFTER the updraft tendencies are calculated...
        Thus, maybe it's best left to filtering after everything is calculated?

        filters are as implemented in affect_filter!() which is in this file as well which calls filter_updraft_vars() which is also in this file
        However, filtering is done at the beginning of ∑tendencies!(), not the end... so any limiter here won't know the gm_tendencies
        =#

        # apply the entr/detr after having calculated all the other tendencies. Not we don't need any limit_tendency() here because the resolver is self limiting.
        if use_tendency_resolver && use_tendency_resolver_on_full_tendencies

            if edmf.entrainment_type isa FractionalEntrModel

                @. tends_ρaθ_liq_ice += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaθ_liq_ice / ρarea, prog_gm.ρθ_liq_ice / ρ_c), prog_gm.ρθ_liq_ice / ρ_c, tends_ρarea,
                    (ρarea * Ic(w_up) * entr_turb_dyn * θ_liq_ice_en) - (ρaθ_liq_ice * Ic(w_up) * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaθ_liq_ice)

                @. tends_ρaq_tot += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_tot / ρarea, prog_gm.ρq_tot / ρ_c), prog_gm.ρq_tot / ρ_c, tends_ρarea,
                    (ρarea * Ic(w_up) * entr_turb_dyn * q_tot_en) - (ρaq_tot * Ic(w_up) * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_tot)

                if edmf.moisture_model isa NonEquilibriumMoisture
                    @. tends_ρaq_liq += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        (ρarea * Ic(w_up) * entr_turb_dyn * q_liq_en) - (ρaq_liq * Ic(w_up) * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_liq)

                    @. tends_ρaq_ice += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_ice / ρarea, prog_gm.q_ice), prog_gm.q_ice, tends_ρarea,
                        (ρarea * Ic(w_up) * entr_turb_dyn * q_ice_en) - (ρaq_ice * Ic(w_up) * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_ice)
                end

            elseif edmf.entrainment_type isa TotalRateEntrModel

                @. tends_ρaθ_liq_ice += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaθ_liq_ice / ρarea, prog_gm.ρθ_liq_ice / ρ_c), prog_gm.ρθ_liq_ice / ρ_c, tends_ρarea,
                    (ρarea * entr_rate_inv_s * θ_liq_ice_en) - (ρarea * detr_rate_inv_s * θ_liq_ice_up), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaθ_liq_ice)
                
                @. tends_ρaq_tot += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_tot / ρarea, prog_gm.ρq_tot / ρ_c), prog_gm.ρq_tot / ρ_c, tends_ρarea,
                    (ρarea * entr_rate_inv_s * q_tot_en) - (ρarea * detr_rate_inv_s * q_tot_up), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_tot)

                if edmf.moisture_model isa NonEquilibriumMoisture
                    @. tends_ρaq_liq += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        (ρarea * entr_rate_inv_s * q_liq_en) - (ρarea * detr_rate_inv_s * q_liq_up), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_liq)

                    @. tends_ρaq_ice += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_ice / ρarea, prog_gm.q_ice), prog_gm.q_ice, tends_ρarea,
                        (ρarea * entr_rate_inv_s * q_ice_en) - (ρarea * detr_rate_inv_s * q_ice_up), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_ice)
                end
                
            end
        end


        # # Solve for updraft velocity
        ρaw = prog_up_f[i].ρaw
        tends_ρaw = tendencies_up_f[i].ρaw
        nh_pressure = aux_up_f[i].nh_pressure
        a_up = aux_up[i].area
        w_up = aux_up_f[i].w
        w_en = aux_en_f.w
        entr_w = aux_up[i].entr_turb_dyn
        detr_w = aux_up[i].detr_turb_dyn
        entr_rate_inv_s = aux_up[i].entr_rate_inv_s
        detr_rate_inv_s = aux_up[i].detr_rate_inv_s
        buoy = aux_up[i].buoy

        @. tends_ρaw = -(∇f(wvec(LBC(ρaw * w_up))))

        # w_en can go negative so these limiters don't make sense..., if it's a problem we can set the limiter based on CFL using new w_en and w_up estimates but we already limit CFL in the updraft so maybe it's best to just do nothing... not sure.
        if edmf.entrainment_type isa FractionalEntrModel
            @. tends_ρaw +=
                (ρaw * (I0f(entr_w) * w_en - I0f(detr_w) * w_up)) + (ρ_f * ᶠinterp_a(a_up) * I0f(buoy)) + nh_pressure
        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. tends_ρaw +=
                ρ_f * ᶠinterp_a(a_up) * (ᶠinterp_a(entr_rate_inv_s) * w_en - ᶠinterp_a(detr_rate_inv_s) * w_up) +
                (ρ_f * ᶠinterp_a(a_up) * I0f(buoy)) +
                nh_pressure
        end

        buoyf = I0f.(buoy)
        a_upf = I0f.(a_up)
        @inbounds for k in real_face_indices(grid) 
            # z = grid.zf[k].z
            if iszero(w_up[k]) && (buoyf[k] > 0) && (a_upf[k] > 0)
                if tends_ρaw[k] ≤ 0

                    # if `pressure_normalmode_buoy_coeff1` > 1, then `nh_pressure` via `nh_pressure_b` will prevent an updraft forming! 
                    # This is because it's defined as  -α_b / (1 + α₂_asp_ratio²) * ρ_f * ᶠinterp_a(a_up) * Ifb(b_up), and α₂_asp_ratio² is set to 0. so you're just left with -α_b * ρ_f * ᶠinterp_a(a_up) * Ifb(b_up) which matches the above
                    #
                    # tends_ρaw[k] -= nh_pressure[k] # remove the nh_pressure contribution # i think this is bad and leads to huge spikes in w...?
                    # e.g. the area might be going down but we can't speculate what happens to w...

                    tends_ρaw[k] = eps(FT) # just some small positive number

                    # @error "We have buoyancy and no updraft velocity yet failed to generate an updraft velocity... inputs were z = $z; w_up = $(w_up[k]); w_en = $(w_en[k]); buoy = $(buoyf[k]); ρaw = $(ρaw[k]); tends_ρaw = $(tends_ρaw[k]); entr_w = $(entr_wf); detr_w = $(detr_wf); entr_rate_inv_s = $(entr_rate_inv_sf); detr_rate_inv_s = $(detr_rate_inv_sf); advection = $(advection[k]); nh_pressure = $(nh_pressure[k]); ρ_f = $(ρ_f[k]); a_upf = $(a_upf[k]);"
               
                
                end
            end
        end

        # TESTING CONVECTIVE TKE [[ NOTE: this will dilute the updrafts, so we should probably have some heuristic around that... ]]
        if param_set.user_params.use_convective_tke # GRAFT onto UPDRAFT
        # if false

            #= 
                If we're generating downdrafts, they operate like updrafts probably (smaller area and all that) so we don't need to broaden updraft.
                However, since they operate in the environment, we don't need to spark any updraft area,
                    the compensating upwards return flow should be broad and slow, and primarily in the env, so we'll put it in normal TKE.
                In a future model with downdrafts, you might consider sparking coherent downdraft areas.

                Since in principle condensate forms in the updrafts, each updraft should in principle risk termination by downdraft formation
                rather than getting to gently fade out as it loses buoyancy and maybe getting a boost from the subsiding env.
                Then again you could have say surface convection rising into air that is dry with falling snow from a layer above or something.
            =#

            entr_rate_inv_s_tke_conv = param_set.user_params.tke_conv_entr_inv_s # FT(1/600) or so seems ok
            a_en_tke_conv = aux_tc.temporary_1
            @. a_en_tke_conv = ifelse(aux_en.latent_heating > 0,  (aux_en.area / 2) * (1 - exp(-aux_en.tke_convective / FT(0.1^2))), FT(0))
            # @. a_en_tke_conv = (aux_en.area / 2) * (1 - exp(-aux_en.tke_convective / FT(0.1^2)))

            # ramp onto max_area as a function of tke_conv [[ fully on by tke = 0.1² ]]
            # @. Δa_en_tke_conv = min(max(a_en_tke_conv - aux_up[i].area, FT(0)), max(edmf.max_area - aux_up_i.area, FT(0))) # ensure we don't entrain more than available area

            Δa_en_tke_conv = a_en_tke_conv # alias
            @inbounds for k in real_center_indices(grid)
                if aux_en.latent_heating[k] > 0
                    Δa_en_tke_conv[k] = min(max(a_en_tke_conv[k] - aux_up[i].area[k], FT(0)), max(edmf.max_area - aux_up_i.area[k], FT(0)))
                else
                    Δa_en_tke_conv[k] = -aux_up_i.area[k]
                end
            end
            tends_ρarea_tke_conv = aux_tc.temporary_2
            @. tends_ρarea_tke_conv = limit_tendency(edtl, ρ_c * Δa_en_tke_conv * entr_rate_inv_s_tke_conv, FT(0), Δa_en_tke_conv * ρ_c, Δt)


            q_tot_en_tke_conv = aux_tc.temporary_1
            θ_liq_ice_en_tke_conv = aux_tc.temporary_3
            @. q_tot_en_tke_conv = aux_en.q_tot + sqrt(aux_en.QTvar)*sqrt(2/π) # entrain the moist unstable air [[ mean of upper half of gaussian ]]
            @. θ_liq_ice_en_tke_conv = aux_en.θ_liq_ice - sqrt(aux_en.Hvar)*sqrt(2/π) # entrain the warm unstable air [[ mean of upper half of gaussian ]]

            @inbounds for k in real_center_indices(grid)
                if aux_en.latent_heating[k] > 0
                    tends_ρarea[k] += tends_ρarea_tke_conv[k]
                    tends_ρaq_tot[k] += limit_tendency(edtl, tends_ρarea_tke_conv[k] * q_tot_en_tke_conv[k], FT(0), FT(Inf), Δt)
                    tends_ρaθ_liq_ice[k] += limit_tendency(edtl, tends_ρarea_tke_conv[k] * θ_liq_ice_en_tke_conv[k], FT(0), FT(Inf), Δt)
                    if edmf.moisture_model isa NonEquilibriumMoisture
                        tends_ρaq_liq[k] += limit_tendency(edtl, tends_ρarea_tke_conv[k] * aux_en.q_liq[k], FT(0), FT(Inf), Δt)
                        tends_ρaq_ice[k] += limit_tendency(edtl, tends_ρarea_tke_conv[k] * aux_en.q_ice[k], FT(0), FT(Inf), Δt)
                    end
                    w_convective = sqrt(2 * aux_en.tke_convective[k])
                    tendencies_en.ρatke_convective[k] -= tends_ρarea_tke_conv[k] * w_convective^2 # remove tke we are sending to updraft
                else
                    tends_ρarea[k] += limit_tendency(edtl, -tends_ρarea_tke_conv[k], prog_up[i].ρarea[k], FT(Inf), Δt) # detraining the tke convective area
                    tends_ρaq_tot[k] += limit_tendency(edtl, -tends_ρarea_tke_conv[k] * q_tot_en_tke_conv[k], prog_up[i].ρaq_tot[k], FT(Inf), Δt)
                    tends_ρaθ_liq_ice[k] += limit_tendency(edtl, -tends_ρarea_tke_conv[k] * θ_liq_ice_en_tke_conv[k], prog_up[i].ρaθ_liq_ice[k], FT(Inf), Δt)
                    if edmf.moisture_model isa NonEquilibriumMoisture
                        tends_ρaq_liq[k] += limit_tendency(edtl, -tends_ρarea_tke_conv[k] * aux_en.q_liq[k], prog_up[i].ρaq_liq[k], FT(Inf), Δt)
                        tends_ρaq_ice[k] += limit_tendency(edtl, -tends_ρarea_tke_conv[k] * aux_en.q_ice[k], prog_up[i].ρaq_ice[k], FT(Inf), Δt)
                    end
                    w_convective = sqrt(2 * aux_en.tke_convective[k])
                    tendencies_en.ρatke_convective[k] += tends_ρarea_tke_conv[k] * w_convective^2 # add back tke we are detraining from updraft
                end
            end

            tends_ρarea_tke_conv_f = @. I0f(tends_ρarea_tke_conv) # face version
            tke_convective_f = aux_tc_f.temporary_f1
            @. tke_convective_f = I0f(aux_en.tke_convective)
            latent_heating_f = aux_tc_f.temporary_f2
            @. latent_heating_f = I0f(aux_en.latent_heating)
            @inbounds for k in real_face_indices(grid)
                w_convective = sqrt(2 * tke_convective_f[k])
                if latent_heating_f[k] > 0
                    tends_ρaw[k] += limit_tendency(edtl, tends_ρarea_tke_conv_f[k] * w_convective, FT(0), FT(Inf), Δt)
                else
                    tends_ρaw[k] += limit_tendency(edtl, -tends_ρarea_tke_conv_f[k] * w_convective, prog_up_f[i].ρaw[k], FT(Inf), Δt)
                end
            end

            
            # @. tends_ρarea += tends_ρarea_tke_conv
            # @. tends_ρaq_tot += limit_tendency(edtl, tends_ρarea_tke_conv * q_tot_en_tke_conv, FT(0), FT(Inf), Δt)
            # @. tends_ρaθ_liq_ice += limit_tendency(edtl, tends_ρarea_tke_conv * θ_liq_ice_en_tke_conv, FT(0), FT(Inf), Δt)
            # if edmf.moisture_model isa NonEquilibriumMoisture
            #     @. tends_ρaq_liq += limit_tendency(edtl, tends_ρarea_tke_conv * aux_en.q_liq, FT(0), FT(Inf), Δt)
            #     @. tends_ρaq_ice += limit_tendency(edtl, tends_ρarea_tke_conv * aux_en.q_ice, FT(0), FT(Inf), Δt)
            # end

            # @. w_convective = sqrt(2 * I0f(aux_en.tke_convective))
            # @. tends_ρaw += limit_tendency(edtl, I0f(tends_ρarea_tke_conv) * w_convective, FT(0), FT(Inf), Δt)
            # @. tendencies_en.ρatke_convective -= tends_ρarea_tke_conv * Ic(w_convective)^2 # remove tke we are sending to updraft
        end

        tends_ρaw[kf_surf] = 0
    end



    return nothing
end




function filter_gm_vars(edmf::EDMFModel, state::State)
    FT = float_type(state)
    
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ

    # force individual gm tendencies to not deplete themselves
    # prog_gm.ρq_tot .= max.(prog_gm.ρq_tot, 0)
    @. prog_gm.ρq_tot = clamp(prog_gm.ρq_tot, FT(0), max(FT(0), prog_gm.ρ - eps(FT))) # specific humidity cannot exceed 1 (at 1 you'll get inf)
    # @. prog_gm.ρθ_liq_ice = max(prog_gm.ρθ_liq_ice, 0) # If you get here, the model is crashing anyway... technically this should be applied ealier
    @. prog_gm.ρθ_liq_ice = max(prog_gm.ρθ_liq_ice, prog_gm.ρ * eps(FT)) # If you get here, the model is crashing anyway... technically this should be applied ealier



    if edmf.moisture_model isa NonEquilibriumMoisture
        # I think if we get here, we're already cooked...
        # @. prog_gm.q_liq = safe_clamp(prog_gm.q_liq, FT(0), FT(0.5) - eps(FT)) # limit to 1/2 so sum doesn't exceed 1
        # @. prog_gm.q_ice = safe_clamp(prog_gm.q_ice, FT(0), FT(0.5) - eps(FT)) # limit to 1/2 so sum doesn't exceed 1

        @. prog_gm.q_liq = safe_clamp(prog_gm.q_liq, FT(0), max(FT(0), prog_gm.ρq_tot / (2ρ_c) - eps(FT))) # limit to 1/2  of avialble qt so we still have q_vap
        @. prog_gm.q_ice = safe_clamp(prog_gm.q_ice, FT(0), max(FT(0), prog_gm.ρq_tot / (2ρ_c) - eps(FT))) # limit to 1/2 of available qt so we still have q_vap
         
        # @. prog_gm.ρq_tot = max.(prog_gm.ρq_tot, ρ_c * (prog_gm.q_liq + prog_gm.q_ice) + eps(FT)) # ensure that the total specific humidity is at least the sum of the liquid and ice specific humidities
    end
    
    return nothing

end

function filter_updraft_vars(edmf::EDMFModel, state::State, surf::SurfaceBase)
    grid = Grid(state)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = float_type(state)
    N_up = n_updrafts(edmf)

    prog_up = center_prog_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    a_min = edmf.minimum_area
    a_max = edmf.max_area


    f_lim_q = FT(3) # 100% error on q_tot is already quite large. but a large updraft could bring large amounts of water up... That would ideally add to the grid mean though...
    f_lim_θ = FT(1.1) # 20% error on θ_li is huge, that's like 50+K Maybe a super deep updraft could pull it off, but updrafts by definition imply thetali isn't wildly varying


    @inbounds for i in 1:N_up
        @. prog_up[i].ρarea = max(prog_up[i].ρarea, 0)
        @. prog_up[i].ρaθ_liq_ice = max(prog_up[i].ρaθ_liq_ice, 0)
        # prog_up[i].ρaq_tot .= max.(prog_up[i].ρaq_tot, 0)
        @. prog_up[i].ρaq_tot = safe_clamp(prog_up[i].ρaq_tot, FT(0),  max(FT(0), prog_up[i].ρarea - eps(FT))) # specific humidity cannot exceed 1 (at 1 you'll get inf)

        # maybe we could filter ρaq_tot and ρaθ_liq_ice so that qt and θ_liq_ice are within 1/f and f times the grid mean? [ we re-enforce this again at the end, but doing it here is good for intermediate calculations]
        @. prog_up[i].ρaq_tot = safe_clamp(prog_up[i].ρaq_tot, prog_up[i].ρarea  * prog_gm.ρq_tot / ρ_c / f_lim_q, prog_up[i].ρarea * prog_gm.ρq_tot / ρ_c * f_lim_q)
        @. prog_up[i].ρaθ_liq_ice = safe_clamp(prog_up[i].ρaθ_liq_ice, prog_up[i].ρarea * prog_gm.ρθ_liq_ice / ρ_c / f_lim_θ, prog_up[i].ρarea * prog_gm.ρθ_liq_ice / ρ_c * f_lim_θ)

        if edmf.entr_closure isa PrognosticNoisyRelaxationProcess
            @. prog_up[i].ε_nondim = max(prog_up[i].ε_nondim, 0)
            @. prog_up[i].δ_nondim = max(prog_up[i].δ_nondim, 0)
        end
        @inbounds for k in real_center_indices(grid)
            # if a_max < FT(.4)
            #     a_max_for_limiting = a_max * 2 # use a slightly looser limit to give the detrainment a chance to work just in case....
            # else
            #     a_max_for_limiting = a_max # use the exact limit
            # end
            a_max_for_limiting = a_max + FT(0.5) * (1-a_max) # use a slightly looser limit to give the detrainment a chance to work just in case....
                
            if !iszero(prog_up[i].ρarea[k]) && (prog_up[i].ρarea[k] > (ρ_c[k] * a_max_for_limiting)) # if the area is less than the minimum, we need to scale it up
                scaling_factor = ρ_c[k] * a_max_for_limiting / prog_up[i].ρarea[k] # i think this is loose, if we scale down the area, we need to be careful to do tracers as well or else we are inducing a tracer change.
                # prog_up[i].ρarea[k] = min(prog_up[i].ρarea[k], ρ_c[k] * a_max) # should we not instead be scaling down the sum of all updrafts to at most a_max?

                prog_up[i].ρarea[k] *= scaling_factor
                # scale other tracers by this factor [[ really i think we oughtta just turn off that limiter completely but we'll see. eithr way this should stop the huge temperature excursions perhaps? ]]
                prog_up[i].ρaθ_liq_ice[k] *= scaling_factor
                prog_up[i].ρaq_tot[k] *= scaling_factor
                if edmf.moisture_model isa NonEquilibriumMoisture
                    prog_up[i].ρaq_liq[k] *= scaling_factor
                    prog_up[i].ρaq_ice[k] *= scaling_factor
                end
            end

            #= If the area does not go to 0 but the tracer does, bad things may be ahead.
            We cannot physically have θ_liq_ice = 0 if ρarea > 0, for example. qt is more debatable but probably should also not be 0.
            So, we will enforce that if the area is not 0, and the tracer is *exactly* 0, we will set it to the grid mean value.
            The *exactly* is important, the tracer value can be less than the grid mean, it's just that if it reaches 0 while the area hasn't something is wrong (and the model will crash w/ T=0 for example).
            =#
            if (prog_up[i].ρarea[k] > 0) && iszero(prog_up[i].ρaθ_liq_ice[k])
                prog_up[i].ρaθ_liq_ice[k] = prog_gm.ρθ_liq_ice[k] / ρ_c[k] * prog_up[i].ρarea[k]
            end
            if (prog_up[i].ρarea[k] > 0) && iszero(prog_up[i].ρaq_tot[k])# this one is more debatable but ok. really qt is never practically 0 but you could imagine some initialization where it is
                prog_up[i].ρaq_tot[k] = prog_gm.ρq_tot[k] / ρ_c[k] * prog_up[i].ρarea[k]
            end

        end
        if edmf.moisture_model isa NonEquilibriumMoisture
            # @. prog_up[i].ρaq_liq = max(prog_up[i].ρaq_liq, 0)
            # @. prog_up[i].ρaq_ice = max(prog_up[di].ρaq_ice, 0)


            @. prog_up[i].ρaq_liq = safe_clamp(prog_up[i].ρaq_liq, FT(0), max(FT(0), prog_up[i].ρaq_tot / 2 - eps(FT))) # try to keep some vapor around but don't pass 1
            @. prog_up[i].ρaq_ice = safe_clamp(prog_up[i].ρaq_ice, FT(0), max(FT(0), prog_up[i].ρaq_tot / 2 - eps(FT))) # try to keep some vapor around but don't pass 1

            # ensure that q_tot is at least more than the sum of q_liq and q_ice
            # @. prog_up[i].ρaq_tot = max(prog_up[i].ρaq_tot, prog_up[i].ρaq_liq + prog_up[i].ρaq_ice + eps(FT)) # specific humidity cannot exceed 1 (at 1 you'll get inf)

        end

        # prog_up_f[i].ρaw .= max.(prog_up_f[i].ρaw, 0)


        # -- my addition to stabilize model ----------------------------------------- #
        # At some point, a large enough external forcing must be thought of as impacting the updraft if necessary...
        # if the value are too large, these filters don't stop that... so the environment has unlimited reduction... 
        # It's hard to do at this point because the environment should be allowed to decrease, you could clamp at this point using gm to negative environment area (does this scale equally between area and tracers?)
        # we know e.g. how entrainment/detrainment behave but don't for example know differences in microphysics, sedimentation, advection etc...
        # the best limit we can have is not exceeding putting everything in the updraft I suppose.

        # Filtering happens in ∑_tendencies before aux is updated, then the next set of tendencies are calculated...
        # Thus, the previous gm tendencies should already be applied...? I think...
        
        @. prog_up[i].ρaθ_liq_ice = safe_clamp(prog_up[i].ρaθ_liq_ice, FT(0), prog_gm.ρθ_liq_ice)
        @. prog_up[i].ρaq_tot = safe_clamp(prog_up[i].ρaq_tot, FT(0), prog_gm.ρq_tot)
        if edmf.moisture_model isa NonEquilibriumMoisture
            @. prog_up[i].ρaq_liq = safe_clamp(prog_up[i].ρaq_liq, FT(0), ρ_c .* prog_gm.q_liq)
            @. prog_up[i].ρaq_ice = safe_clamp(prog_up[i].ρaq_ice, FT(0), ρ_c .* prog_gm.q_ice)
        end
    end
    # apply clipping at 0 and minimum area to ρaw
    @inbounds for i in 1:N_up
        @. prog_up_f[i].ρaw = max.(prog_up_f[i].ρaw, 0)
        @. prog_up_f[i].ρaw = Int(ᶠinterp_a(prog_up[i].ρarea) >= ρ_f * a_min) * prog_up_f[i].ρaw
    end

    # no penetration at the top of the atmosphere
    kf_toa = kf_top_of_atmos(grid)
    kc_toa = kc_top_of_atmos(grid)
    @inbounds for i in 1:N_up
        prog_up_f[i].ρaw[kf_toa] = 0
        prog_up[i].ρarea[kc_toa] = 0
        prog_up[i].ρaθ_liq_ice[kc_toa] = 0
        prog_up[i].ρaq_tot[kc_toa] = 0
        if edmf.moisture_model isa NonEquilibriumMoisture
            prog_up[i].ρaq_liq[kc_toa] = 0
            prog_up[i].ρaq_ice[kc_toa] = 0
        end
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:N_up
            is_surface_center(grid, k) && continue
            prog_up[i].ρaq_tot[k] = max(prog_up[i].ρaq_tot[k], 0)
            # this is needed to make sure Rico is unchanged.
            # TODO : look into it further to see why
            # a similar filtering of ρaθ_liq_ice breaks the simulation
            
            if prog_up[i].ρarea[k] / ρ_c[k] < a_min
                # prog_up[i].ρaq_tot[k] = 0 # is this even a good idea? why should qt be 0? at least maybe we should take the grid mean value?
                # prog_up[i].ρaθ_liq_ice[k] = 0 ######################## NOTE THIS BROKE THE SIMS, SO IF THEY'RE STLL BREAKING CHECK HERE!!!!! ( actually i think this is dumb, bc if you don't set area to 0 it still needs a temperature. we do interpret as 0 in aux, so we don't need to here) also our `factor` limiter won't allow such an excursion anyway
                
                # testing just going to grid mean value.. it's what we do w/ aux anyway... either you do that or you set area here to 0 idk...
                prog_up[i].ρaq_tot[k] = prog_gm.ρq_tot[k] / ρ_c[k] * prog_up[i].ρarea[k] # this is the grid mean value, so we don't have a problem with the updraft being too small
                prog_up[i].ρaθ_liq_ice[k] = prog_gm.ρθ_liq_ice[k] / ρ_c[k] * prog_up[i].ρarea[k]

                if edmf.moisture_model isa NonEquilibriumMoisture
                    prog_up[i].ρaq_liq[k] = 0
                    prog_up[i].ρaq_ice[k] = 0
                end
            end
            #
            # if ρaq_tot or ρaθ_liq_ice or goes to 0 and gm is not 0, we should prolly just detrain (hopefullly this fixes the tiny small area problem)
            if iszero(prog_up[i].ρaq_tot[k]) || iszero(prog_up[i].ρaθ_liq_ice[k])
                if !iszero(prog_gm.ρq_tot[k])
                    prog_up[i].ρarea[k] = 0
                    prog_up[i].ρaq_tot[k] = 0
                    prog_up[i].ρaθ_liq_ice[k] = 0
                    if edmf.moisture_model isa NonEquilibriumMoisture
                        prog_up[i].ρaq_liq[k] = 0
                        prog_up[i].ρaq_ice[k] = 0
                    end
                end
            end
            #
        end
    end


    # Ic = CCO.InterpolateF2C()
    @inbounds for i in 1:N_up
        # test not doing this instant detrainment..... at the least we could maybe to to grid mean? or force detrainment to actually do it... idk... also it breaks things for later mixing only really makes sense at updraft top perhaps but now we have elevated convection.
        # the real problem is that it then spikes the gradients of everything else just bc w went to 0... so then ρarea or ρaq gradients get all jacked up 

        if edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftDoNothing
            ; # do nothing, we will just let the updraft linger
        elseif edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftKill
            @. prog_up[i].ρarea = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρarea)  # trying not doing this [ I wonder if this is still bad now that we've added other fixes... adding it back would fix the drifting updraft area problem... you get if you turn off mix_stalled_updraft_to_grid_mean. keeping the area but setting the values to grid mean is catastropic to the microphysics N/supersat etc]
            @. prog_up[i].ρaθ_liq_ice = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρaθ_liq_ice) # trying not doing this
            @. prog_up[i].ρaq_tot = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρaq_tot) # trying not doing this...
        elseif edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftMixToGridMean
            @. prog_up[i].ρaθ_liq_ice = ifelse(Ic(prog_up_f[i].ρaw) <= 0, prog_gm.ρθ_liq_ice * prog_up[i].ρarea / ρ_c, prog_up[i].ρaθ_liq_ice) # test using gm | not ideal but we don't know detrainment timescale in absence of w_up so to avoid languishing w/o adding another calibrated parameter maybe it's better?
            @. prog_up[i].ρaq_tot = ifelse(Ic(prog_up_f[i].ρaw) <= 0, prog_gm.ρq_tot * prog_up[i].ρarea / ρ_c, prog_up[i].ρaq_tot) # test using gm 
        elseif edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftDetrainDowndrafts # does nothing but ensure that the updraft θ_liq_ice does not decrease, by reducing the area
            # if prog_up[i].ρaθ_liq_ice has gone down, we reduce prog_up[i].ρarea so that prog_up[i].ρaθ_liq_ice / prog_up[i].ρarea = prog_gm.ρθ_liq_ice / ρ_c
            ρaw_c = Ic.(prog_up_f[i].ρaw) 
            @inbounds for k in real_center_indices(grid)
                if ρaw_c[k] <= 0
                    if !iszero(prog_up[i].ρarea[k])
                        if (prog_up[i].ρaθ_liq_ice[k] / prog_up[i].ρarea[k]) < (prog_gm.ρθ_liq_ice[k] / ρ_c[k])
                            # reduce the area so that the ratio is correct, and θ_liq_ice_up = θ_liq_ice_gm
                            scaling_factor = (prog_up[i].ρaθ_liq_ice[k] / prog_up[i].ρarea[k]) / (prog_gm.ρθ_liq_ice[k] / ρ_c[k]) # f = (ρaθ)_up / ((ρa)_up * θ_gm)...
                            # if scaling_factor > 1
                            #     error("Stalled updraft handler StalledUpdraftDetrainDowndrafts tried to increase area, which is not allowed. From inputs ρaθ_liq_ice=$(prog_up[i].ρaθ_liq_ice[k]), ρarea=$(prog_up[i].ρarea[k]), θ_gm=$(prog_gm.ρθ_liq_ice[k] / ρ_c[k]), factor=$factor")
                            # end
                            prog_up[i].ρarea[k] *= scaling_factor
                            # do nothing to prog_up[i].ρaθ_liq_ice, now that we've decreased the area...
                            prog_up[i].ρaq_tot[k] *= scaling_factor # not sure what to do with qt... was the requisite part already detrained and we do nothing? or do we need to decrease qt as we decrease area? I think we do or else if it gets large enough, the environment can be forced negative. w/ a_max = 0.33 this is not entirely impossible.
                            if edmf.moisture_model isa NonEquilibriumMoisture # moved this up here because we won't be able to calculate the factor afterwards
                                prog_up[i].ρaq_liq[k] *= scaling_factor # not sure what to do with qt... was the requisite part already detrained and we do nothing? or do we need to decrease qt as we decrease area?
                                prog_up[i].ρaq_ice[k] *= scaling_factor # not sure what to do with qt... was the requisite part already detrained and we do nothing? or do we need to decrease qt as we decrease area?
                            end
                        else
                            # we can still see gains in buoyancy that don't seem to trigger convection? idk why...
                            # TODO: if theta is going up but w is 0, idk...
                        end
                    end
                end
            end
        else
            error("Unknown stalled updraft handler type: $(edmf.entrainment_type.stalled_updraft_handler)")
        end


        if edmf.moisture_model isa NonEquilibriumMoisture
            if edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftDoNothing
                ;
            elseif edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftKill
                @. prog_up[i].ρaq_liq = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρaq_liq) # trying not doing this
                @. prog_up[i].ρaq_ice = ifelse(Ic(prog_up_f[i].ρaw) <= 0, FT(0), prog_up[i].ρaq_ice) # trying not doing this
            elseif edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftMixToGridMean
                @. prog_up[i].ρaq_liq = ifelse(Ic(prog_up_f[i].ρaw) <= 0, prog_gm.q_liq * prog_up[i].ρarea / ρ_c, prog_up[i].ρaq_liq) # test using gm
                @. prog_up[i].ρaq_ice = ifelse(Ic(prog_up_f[i].ρaw) <= 0, prog_gm.q_ice * prog_up[i].ρarea / ρ_c, prog_up[i].ρaq_ice) # test using gm
            elseif edmf.entrainment_type.stalled_updraft_handler isa StalledUpdraftDetrainDowndrafts 
                ; # moved above bc relies on factors
            end # should be StalledUpdraftDoNothing [we have background detrainment now so it might work? idk]
        end


        # re-enforce this in case any of our other limiting broke it...
        @. prog_up[i].ρaq_tot = safe_clamp(prog_up[i].ρaq_tot, prog_up[i].ρarea * prog_gm.ρq_tot / ρ_c / f_lim_q, prog_up[i].ρarea * prog_gm.ρq_tot / ρ_c * f_lim_q)
        @. prog_up[i].ρaθ_liq_ice = safe_clamp(prog_up[i].ρaθ_liq_ice, prog_up[i].ρarea * prog_gm.ρθ_liq_ice / ρ_c / f_lim_θ, prog_up[i].ρarea * prog_gm.ρθ_liq_ice / ρ_c * f_lim_θ)
        # ---------------------------------------------------------------------------- #


        if edmf.surface_area_bc isa FixedSurfaceAreaBC || edmf.surface_area_bc isa ClosureSurfaceAreaBC
            a_surf = area_surface_bc(surf, edmf, i, edmf.surface_area_bc)
            prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * a_surf
        elseif edmf.surface_area_bc isa PrognosticSurfaceAreaBC
            # don't let the area go to 0 if the bflux is positive [updraft seems to never come back if it ever goes away]
            if iszero(prog_up[i].ρarea[kc_surf]) && (surf.bflux > 0) # not sure if this will work, maybe when it goes to 0 it's because bflux is negative.... [ update, we apparently did trigger this...]
                prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * max(a_min, FT(0.01))
                # set variables to grid mean
                prog_up[i].ρaθ_liq_ice[kc_surf] = prog_gm.ρθ_liq_ice[kc_surf] * (prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf])
                prog_up[i].ρaq_tot[kc_surf] = prog_gm.ρq_tot[kc_surf] * (prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf])
                if edmf.moisture_model isa NonEquilibriumMoisture
                    prog_up[i].ρaq_liq[kc_surf] = prog_gm.q_liq[kc_surf] * (prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf])
                    prog_up[i].ρaq_ice[kc_surf] = prog_gm.q_ice[kc_surf] * (prog_up[i].ρarea[kc_surf] / ρ_c[kc_surf])
                end
            end
        end

        θ_surf = θ_surface_bc(surf, state, N_up, i)
        q_surf = q_surface_bc(surf, state, N_up, i)
        prog_up[i].ρaθ_liq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * θ_surf
        prog_up[i].ρaq_tot[kc_surf] = prog_up[i].ρarea[kc_surf] * q_surf

        if edmf.moisture_model isa NonEquilibriumMoisture
            ql_surf = ql_surface_bc(surf)
            qi_surf = qi_surface_bc(surf)
            prog_up[i].ρaq_liq[kc_surf] = prog_up[i].ρarea[kc_surf] * ql_surf
            prog_up[i].ρaq_ice[kc_surf] = prog_up[i].ρarea[kc_surf] * qi_surf
        end

        # Enforce CFL on updraft velocities. [[ rn we enforce this in update_aux only ... ]]


        
        #=
            When 2 updrafts collide, the middle cell when it goes from 0 to on usually has horrible things happen in w, θli, etc...
            It would be good to investigate why this happens and concoct a cure.
        =#

        # TEST NO UPDRAFT...
        # @. prog_up[i].ρarea = FT(0)
        # @. prog_up[i].ρaθ_liq_ice = FT(0)
        # @. prog_up[i].ρaq_tot = FT(0)
        # if edmf.moisture_model isa NonEquilibriumMoisture
        #     @. prog_up[i].ρaq_liq = FT(0)
        #     @. prog_up[i].ρaq_ice = FT(0)
        # end 
        # @. prog_up_f[i].ρaw = FT(0)

    end
    return nothing
end

function compute_covariance_shear(
    edmf::EDMFModel,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_en_sym},
    ::Val{ψ_en_sym},
) where {covar_sym, ϕ_en_sym, ψ_en_sym}

    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ
    is_tke = covar_sym === :tke
    FT = float_type(state)
    k_eddy = is_tke ? aux_tc.KM : aux_tc.KH
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    uₕ_gm = grid_mean_uₕ(state)

    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    wvec = CC.Geometry.WVector
    ϕ_en = getproperty(aux_en, ϕ_en_sym)
    ψ_en = getproperty(aux_en, ψ_en_sym)

    bcs = (; bottom = CCO.Extrapolate(), top = CCO.SetGradient(wvec(zero(FT))))
    If = CCO.InterpolateC2F(; bcs...)
    area_en = aux_en_c.area
    shear = aux_covar.shear

    C123 = CCG.Covariant123Vector
    local_geometry = CC.Fields.local_geometry_field(axes(ρ_c))
    k̂ = center_aux_turbconv(state).k̂
    @. k̂ = CCG.Contravariant3Vector(CCG.WVector(FT(1)), local_geometry)
    Ifuₕ = uₕ_bcs()
    ∇uvw = CCO.GradientF2C()
    # TODO: k_eddy and Shear² should probably be tensors (Contravariant3 tensor),
    #       so that the result (a contraction) is a scalar.
    if is_tke
        uvw = face_aux_turbconv(state).uvw
        Shear² = center_aux_turbconv(state).Shear²
        @. uvw = C123(Ifuₕ(uₕ_gm)) + C123(wvec(ϕ_en)) # ϕ_en === ψ_en
        @. Shear² = LA.norm_sqr(adjoint(∇uvw(uvw)) * k̂)
        @. shear = ρ_c * area_en * k_eddy * Shear²
    else
        ∇c = CCO.GradientF2C() # This one is Gradient not Divergence so it's not the global one
        @. shear = 2 * ρ_c * area_en * k_eddy * LA.dot(∇c(If(ϕ_en)), k̂) * LA.dot(∇c(If(ψ_en)), k̂)
    end
    return nothing
end

function compute_covariance_interdomain_src(
    # edmf::EDMFModel,
    N_up::Int,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_sym},
    ::Val{ψ_sym},
) where {covar_sym, ϕ_sym, ψ_sym}
    # N_up = n_updrafts(edmf)
    is_tke = covar_sym === :tke
    FT = float_type(state)
    tke_factor = is_tke ? FT(0.5) : 1
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en_2m = center_aux_environment_2m(state)
    interdomain = getproperty(aux_en_2m, covar_sym).interdomain
    prog_up = is_tke ? aux_up_f : aux_up
    aux_en_c = center_aux_environment(state)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    Ic = is_tke ? CCO.InterpolateF2C() : x -> x

    parent(interdomain) .= 0
    @inbounds for i in 1:N_up
        ϕ_up = getproperty(prog_up[i], ϕ_sym)
        ψ_up = getproperty(prog_up[i], ψ_sym)
        a_up = aux_up[i].area
        @. interdomain += tke_factor * a_up * (1 - a_up) * (Ic(ϕ_up) - Ic(ϕ_en)) * (Ic(ψ_up) - Ic(ψ_en))
    end
    return nothing
end

function compute_covariance_entr(
    edmf::EDMFModel,
    state::State,
    ::Val{covar_sym},
    ::Val{ϕ_sym},
    ::Val{ψ_sym},
) where {covar_sym, ϕ_sym, ψ_sym}

    N_up = n_updrafts(edmf)
    FT = float_type(state)
    is_tke = covar_sym === :tke
    tke_factor = is_tke ? FT(0.5) : 1
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_gm_c = center_aux_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ
    gm = is_tke ? prog_gm_f : aux_gm_c
    prog_up = is_tke ? aux_up_f : aux_up
    to_scalar = is_tke ? toscalar : x -> x
    ϕ_gm = getproperty(gm, ϕ_sym)
    ψ_gm = getproperty(gm, ψ_sym)
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_en_c = center_aux_environment(state)
    covar = getproperty(aux_en_c, covar_sym)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    entr_gain = aux_covar.entr_gain
    detr_loss = aux_covar.detr_loss
    # Ic = CCO.InterpolateF2C()
    Idc = is_tke ? Ic : x -> x
    # TODO: we shouldn't need `parent` call here:
    parent(entr_gain) .= 0
    parent(detr_loss) .= 0
    min_area = edmf.minimum_area

    @inbounds for i in 1:N_up
        aux_up_i = aux_up[i]
        frac_turb_entr = aux_up_i.frac_turb_entr
        eps_turb = frac_turb_entr
        detr_sc = aux_up_i.detr_sc
        entr_sc = aux_up_i.entr_sc
        detr_ml = aux_up_i.detr_ml
        entr_ml = aux_up_i.entr_ml
        entr_rate_inv_s = aux_up_i.entr_rate_inv_s
        detr_rate_inv_s = aux_up_i.detr_rate_inv_s
        w_up = aux_up_f[i].w
        prog_up_i = prog_up[i]
        ϕ_up = getproperty(prog_up_i, ϕ_sym)
        ψ_up = getproperty(prog_up_i, ψ_sym)

        a_up = aux_up_i.area


        if edmf.entrainment_type isa FractionalEntrModel
            @. entr_gain +=
                Int(a_up > min_area) * (
                    tke_factor *
                    ρ_c *
                    a_up *
                    abs(Ic(w_up)) *
                    (detr_sc + detr_ml) *
                    (Idc(ϕ_up) - Idc(ϕ_en)) *
                    (Idc(ψ_up) - Idc(ψ_en))
                ) + (
                    tke_factor *
                    ρ_c *
                    a_up *
                    abs(Ic(w_up)) *
                    eps_turb *
                    (
                        (Idc(ϕ_en) - Idc(to_scalar(ϕ_gm))) * (Idc(ψ_up) - Idc(ψ_en)) +
                        (Idc(ψ_en) - Idc(to_scalar(ψ_gm))) * (Idc(ϕ_up) - Idc(ϕ_en))
                    )
                )

            @. detr_loss +=
                Int(a_up > min_area) * tke_factor * ρ_c * a_up * abs(Ic(w_up)) * (entr_sc + entr_ml + eps_turb) * covar

        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. entr_gain +=
                Int(a_up > min_area) *
                (tke_factor * ρ_c * a_up * detr_rate_inv_s * (Idc(ϕ_up) - Idc(ϕ_en)) * (Idc(ψ_up) - Idc(ψ_en)))

            @. detr_loss += Int(a_up > min_area) * tke_factor * ρ_c * a_up * entr_rate_inv_s * covar

        end

    end

    return nothing
end

function compute_covariance_dissipation(
    edmf::EDMFModel,
    state::State,
    ::Val{covar_sym},
    param_set::APS,
) where {covar_sym}
    FT = float_type(state)
    c_d = mixing_length_params(edmf).c_d
    aux_tc = center_aux_turbconv(state)
    prog_en = center_prog_environment(state)
    prog_gm = center_prog_grid_mean(state)
    aux_en = center_aux_environment(state)
    ρ_c = prog_gm.ρ
    aux_en_2m = center_aux_environment_2m(state)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    covar = getproperty(aux_en, covar_sym)
    dissipation = aux_covar.dissipation
    area_en = aux_en.area
    tke_en = aux_en.tke
    mixing_length = aux_tc.mixing_length

    if param_set.user_params.use_convective_tke && (covar_sym === :tke)
        @. dissipation = ρ_c * area_en * max(covar - aux_en.tke_convective, 0) * max(max(tke_en - aux_en.tke_convective, 0), 0)^FT(0.5) / max(mixing_length, FT(1.0e-3)) * c_d # we handle convective tke dissipation separately since it has different mixing length, coherence, etc.
    else
        @. dissipation = ρ_c * area_en * covar * max(tke_en, 0)^FT(0.5) / max(mixing_length, FT(1.0e-3)) * c_d
    end
    return nothing
end

function compute_en_tendencies!(
    edmf::EDMFModel,
    state::State,
    param_set::APS,
    surf::SurfaceBase,
    ::Val{covar_sym},
    ::Val{prog_sym},
    Δt::Real,
    cfl_limit::Real,
    use_fallback_tendency_limiters::Bool,
) where {covar_sym, prog_sym}
    grid = Grid(state)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    aux_gm_f = face_aux_grid_mean(state)
    prog_en = center_prog_environment(state)
    prog_gm = center_prog_grid_mean(state)
    aux_en_2m = center_aux_environment_2m(state)
    tendencies_en = center_tendencies_environment(state)
    tend_covar = getproperty(tendencies_en, prog_sym)
    prog_covar = getproperty(prog_en, prog_sym)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    prog_pr = center_prog_precipitation(state)
    covar = getproperty(aux_en, covar_sym)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_up = center_aux_updrafts(state)
    w_en_f = face_aux_environment(state).w
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    c_d = mixing_length_params(edmf).c_d
    mix_len_params = mixing_length_params(edmf)
    is_tke = covar_sym === :tke
    FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)

    ρ_ae_K = face_aux_turbconv(state).ρ_ae_K
    KM = center_aux_turbconv(state).KM
    KH = center_aux_turbconv(state).KH
    aux_tc = center_aux_turbconv(state)
    aux_bulk = center_aux_bulk(state)
    D_env = aux_tc.temporary_1
    aeK = aux_tc.temporary_2
    a_bulk = aux_bulk.area
    if is_tke
        @. aeK = (1 - a_bulk) * KM
    else
        @. aeK = (1 - a_bulk) * KH
    end

    press = aux_covar.press
    buoy = aux_covar.buoy
    shear = aux_covar.shear
    entr_gain = aux_covar.entr_gain
    rain_src = aux_covar.rain_src

    wvec = CC.Geometry.WVector
    aeK_bcs = (; bottom = CCO.Extrapolate(), top = CCO.Extrapolate())

    if edmf.thermo_covariance_model isa PrognosticThermoCovariances
        prog_bcs = (; bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))
    else
        prog_bcs = (; bottom = CCO.Extrapolate(), top = CCO.SetGradient(wvec(FT(0))))
    end

    If = CCO.InterpolateC2F(; aeK_bcs...)
    ∇f = CCO.GradientC2F(; prog_bcs...)
    # ∇c = CCO.DivergenceF2C()

    # compute bottom BC for TKE
    if is_tke
        ρa_e_surf = ρ_c[kc_surf] * aux_en.area[kc_surf]
        u_surf = physical_grid_mean_u(state)[kc_surf]
        v_surf = physical_grid_mean_v(state)[kc_surf]
        U_surf_norm = sqrt(u_surf^2 + v_surf^2)
        ustar = surf.ustar
        surface_tke_turb_flux = get_surface_tke_turb_flux(mix_len_params, ustar, ρa_e_surf, U_surf_norm)
        ∇c_turb = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(surface_tke_turb_flux)))
    else
        ∇c_turb = CCO.DivergenceF2C()
    end

    mixing_length = aux_tc.mixing_length
    min_area = edmf.minimum_area

    area_en = aux_en.area
    tke_en = aux_en.tke

    parent(D_env) .= 0

    @inbounds for i in 1:N_up
        turb_entr = aux_up[i].frac_turb_entr
        entr_sc = aux_up[i].entr_sc
        entr_ml = aux_up[i].entr_ml
        entr_rate_inv_s = aux_up[i].entr_rate_inv_s
        w_up = aux_up_f[i].w
        a_up = aux_up[i].area
        # TODO: using `Int(bool) *` means that NaNs can propagate
        # into the solution. Could we somehow call `ifelse` instead?
        if edmf.entrainment_type isa FractionalEntrModel
            @. D_env += Int(a_up > min_area) * ρ_c * a_up * Ic(w_up) * (entr_sc + entr_ml + turb_entr)
        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. D_env += Int(a_up > min_area) * ρ_c * a_up * entr_rate_inv_s
        end
    end

    RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0)))
    @. tend_covar =
        press + buoy + shear + entr_gain + rain_src - D_env * covar -
        (c_d * sqrt(max(tke_en, 0)) / max(mixing_length, 1)) * prog_covar - ∇c(wvec(RB(prog_covar * Ic(w_en_f)))) -
        ∇c_turb(-1 * ρ_f * If(aeK) * ∇f(covar))


    if is_tke && param_set.user_params.use_convective_tke_production_only # let it advect and diffuse on its own.... [[ this one is not tke dependent ]]
        # production due to instability [[should we use static stability ? ]]
        ρatke_convective_production = aux_en.tke_convective_production # for some reason these have ρa in ρatke but are named this way so we follow
        K_buoy = param_set.user_params.convective_tke_buoyancy_coeff
        a_en = aux_en.area
        ∂MSE∂z = aux_en.∂MSE∂z
        @. ρatke_convective_production = K_buoy * ρ_c * a_en * max(-∂MSE∂z, 0) # only positive values of convective TKE [[ not sure if existing tke should matter...]]
        @. tend_covar += ρatke_convective_production # maybe don't combine?

        # Here we let normal dissipation take care of it... we find it is very slow to diffuse though... it mainly stays where it is created...
    end

    if is_tke && param_set.user_params.use_convective_tke

        (!param_set.user_params.use_convective_tke_production_only) || error("Cannot use convective TKE while also using convective TKE production only")
        Ifw = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))

        a_en = aux_en.area
        # We want here to keep our `convective TKE` and move it around etc, but dont dissipate it unless there's strong stability.
        ρatke_convective = prog_en.ρatke_convective
        tend_ρatke_convective = tendencies_en.ρatke_convective
        ∂MSE∂z = aux_en.∂MSE∂z



        #=
        Compute theoretical max TKE (J/kg) from local MSE stratification.
        Derived from potential energy available over one mixing length:
            Emax = (g / (2(h_0-Φ))) * (-∂h/∂z) * ℓ_mix^2
        where:
            h0 ≈ cp*T + L*qv   (reference moist static energy)
            ∂h/∂z < 0 indicates unstable stratification
            ℓ_mix   is the mixing length (m)
        Used to cap TKE growth at physically realizable values.
        =#
        MSE = aux_en.MSE
        # Φ = TC.geopotential(param_set, grid.zc[k].z)
        Φ = geopotential.(param_set, getfield.(grid.zc, :z))
        g = TCP.grav(param_set)
        # mixing length is set by dry static stability so we prolly don't wanna use that...
        tke_max = aux_tc.temporary_1 # 1 and 2 used above but avaiable again
        @. tke_max = (g/(2*max(MSE - Φ, eps(FT)))) * max(-∂MSE∂z, 0) * ifelse((∂MSE∂z < FT(0)) , max(aux_tc.mixing_length, FT(200)) , aux_tc.mixing_length)^2 * param_set.user_params.tke_convective_max_scaling_factor # max tke we could get from total conversion of MSE gradient to tke


        # confine tke_max to CFL
        Δz = get_Δz(prog_gm.ρ) # This is just a regular array not a CC field so we 
        @inbounds for k in real_center_indices(grid) # because get_Δz returns a regular array not a CC Field, we do this for loop
            if tke_max[k] > 0
                w_max = cfl_limit * Δz[k] / Δt
                tke_max[k] = min(tke_max[k], FT(0.5) * ρ_c[k] * w_max^2) # kinetic energy from max w
            end
        end


        # -- instability calculation -- #
        #= unstable if:
            - subsaturated and dry static instability [[ dry static energy gradient looks like ]]
            - saturated and moist static instability
              - saturated wrt ice if below freezing, liquid if above
        =#

        instability = aux_tc.temporary_2
        @. instability = zero(FT)

        @inbounds for k in real_center_indices(grid)
            ts = aux_en.ts[k]
            T = aux_en.T[k]
            T_freeze = TCP.T_freeze(param_set)
            below_freezing = (T < T_freeze)


            std_qt = sqrt(aux_en.QTvar[k])
            if std_qt > sqrt(eps(FT))
                q_tot_hi = aux_en.q_tot[k] + sqrt(aux_en.QTvar[k])
                # q_tot_lo = aux_en.q_tot[k] - sqrt(aux_en.QTvar[k])
                ts_hi = (edmf.moisture_model isa EquilibriumMoisture) ? thermo_state_pθq(param_set, aux_en.p[k], aux_en.θ_liq_ice[k], q_tot_hi) : thermo_state_pθq(param_set, aux_en.p[k], aux_en.θ_liq_ice[k], q_tot_hi, aux_en.q_liq[k], aux_en.q_ice[k])
                # ts_lo = (edmf.moisture_model isa EquilibriumMoisture) ? thermo_state_pθq(param_set, aux_en.p[k], aux_en.θ_liq_ice[k], q_tot_lo) : thermo_state_pθq(param_set, aux_en.p[k], aux_en.θ_liq_ice[k], q_tot_lo, aux_en.q_liq[k], aux_en.q_ice[k])
                RH_ice = relative_humidity_over_ice(thermo_params, ts_hi)
                RH_liq = relative_humidity_over_liquid(thermo_params, ts_hi)
            else
                if !state.calibrate_io
                    RH_ice = aux_en.RH_ice[k]
                    RH_liq = aux_en.RH_liq[k]
                else
                    RH_ice = relative_humidity_over_ice(thermo_params, ts)
                    RH_liq = relative_humidity_over_liquid(thermo_params, ts)
                end
            end


            # frac_supersat
            # check fraction of qt distribution above saturation..
            q_vap_sat = below_freezing ? TD.q_vap_saturation_generic(thermo_params, T, TD.air_density(thermo_params, ts), TD.Ice()) : TD.q_vap_saturation_generic(thermo_params, T, TD.air_density(thermo_params, ts), TD.Liquid())
            qt_std = sqrt(aux_en.QTvar[k])
            qt_mean = aux_en.q_tot[k]
            cov_qtθ = aux_en.HQTcov[k]
            # θli_mean = aux_en.θ_liq_ice[k]
            θli_std = aux_en.Hvar[k]
            p = aux_en.p[k]
            p0 = TD.TP.p_ref_theta(thermo_params)
            R_d = TCP.R_d(param_set)
            R_v = TCP.R_v(param_set)
            c_p = TCP.cp_d(param_set)
            L_v = TD.latent_heat_vapor(thermo_params, T)
            L_f = TD.latent_heat_fusion(thermo_params, T)
            L_s = TD.latent_heat_sublim(thermo_params, T)

            ql_mean = aux_en.q_liq[k]
            qi_mean = aux_en.q_ice[k]
            
            # frac_supersat_old = (qt_std > 0) ? (one(FT) - Distributions.cdf(Distributions.Normal(qt_mean, qt_std), q_vap_sat)) : ((qt_mean > q_vap_sat) ? one(FT) : zero(FT)) # doesn't handle covar
                            
            β = (p / p0)^(R_d/c_p) # θl -> T conversion factor
            dq_sat_dT = below_freezing ? L_f * q_vap_sat / (R_v * T^2) : L_v * q_vap_sat / (R_v * T^2) # Clausius–Clapeyron slope for linearization
            α = dq_sat_dT * β
            μ_X = qt_mean - q_vap_sat # linear combination mean and variance
            σ2_X = qt_std^2 + (α*θli_std)^2 - 2*α*cov_qtθ
            σ_X = sqrt(max(σ2_X, 0.0))
            frac_supersat = σ_X > 0 ? 1 - Distributions.cdf(Distributions.Normal(0, σ_X), -μ_X) : (μ_X > 0 ? one(FT) : zero(FT)) # fraction supersaturated
            frac_subsat = one(FT) - frac_supersat

            aux_en.frac_supersat[k] = frac_supersat # for debugging

            # the engine efficiency peaks when we are half saturated and half subsaturated... then we can accelerate both updrafts and downdrafts... We aim for 0.25 at frac_supersat = 0 and 1, 1 at frac_supersat = 0.5, efficiency = 1.
            peak_enhancement = FT(8)
            efficiency = one(FT) + (peak_enhancement - one(FT)) * (FT(4) * frac_supersat * (one(FT) - frac_supersat))





            latent_heating = aux_en.latent_heating[k] # to be precise you'd want a split into something in the subsat and something in the supersat part, but this is close enough for now...

            latent_heating = sign(latent_heating) * min(abs(latent_heating), FT(3)) # avoid excessive latent heating (i.e. oscillations)

            if (below_freezing && (RH_ice < one(FT)) || (!below_freezing && (RH_liq < one(FT)))) # subsaturated
                # I'm not sure if you should evaporation to contribue and use MSE or just stick to SE

                if frac_supersat > 0
                    instability[k] = -max(∂MSE∂z[k], -1) * frac_supersat * abs(latent_heating) # partially saturated
                else
                    instability[k] = FT(0) # no instability
                end

                if (aux_en.q_liq[k] + aux_en.q_ice[k] + prog_pr.q_rai[k] + prog_pr.q_sno[k]) > FT(1e-10) # Might not need this since we have 
                    # instability[k] = -∂MSE∂z[k]
                    # instability[k] = -∂DSE∂z[k]
                    # instability[k] += -FT(Inf) # don't allow, DSE is handled fine in I think by normal tke... not sure...

                    # We are unstable to evaporation fueled downdrafts... [[ technically liquid should release more energy but... ]]
                    instability[k] = min(instability[k], 0)
                    instability[k] += (-max(∂MSE∂z[k], -1) * abs(latent_heating)) * frac_subsat
                else # use DSE
                    # Do nothing so as to not overwrite any frac_supersat contribution
                    # instability[k] = -∂DSE∂z[k]
                    # c_p = TCP.cp_d(param_set)
                    # instability[k] = c_p * -dθ_virt_dz[k] # dry static stability
                    # instability[k] += -FT(Inf) # disallow
                end
            else # saturated
                instability[k] = -max(∂MSE∂z[k], -1) * abs(latent_heating) * frac_supersat

                if (frac_subsat > zero(FT)) && ((aux_en.q_liq[k] + aux_en.q_ice[k] + prog_pr.q_rai[k] + prog_pr.q_sno[k]) > FT(1e-10)) # partially subsaturated with condensate
                    instability[k] += (-max(∂MSE∂z[k], -1) * abs(latent_heating)) * frac_subsat
                end
            end
            instability[k] *= efficiency # when the downdrafts/updrafts don't have to fight their natural buoyancy you unlock a world of extra instability... its essentially undamped. you could think of it as having the circulation tucked into 1/2 the domain at 1/2 (KE down by 4, area down by 2) vs having the full domain. You get to cycle and reuse the same latent heat indefinitely
        end

        g = TCP.grav(param_set)
        c_p = TCP.cp_d(param_set)
        T = aux_en.T
        @. instability *= (g/(c_p * T)) # dw/dt = g/(c_p T)


        # -- production -- #

        #=
            Derivation of dw/dt and dKE/dt from a single MSE profile

            1. Start from buoyancy definition for a parcel in an environment:
                b = g/(c_p * T) * (MSE_parcel - MSE_environment)
            Since we only have one MSE profile, we treat the local environmental MSE as reference,
            and the parcel is assumed to be a small perturbation moving through the profile.

            2. As the parcel moves upward/downward, it gains/loses buoyancy from two sources:
            a) Vertical gradient of MSE: moving into lower MSE (dMSE/dz < 0) increases b.
                Using finite differences, the local contribution is:
                    b_gradient ≈ - g/(c_p * T) * dMSE/dz
            b) Latent heating (LH) adds or removes energy at the level, directly changing buoyancy:
                    b_heating ≈ g/(c_p * T) * LH

            3. The total instantaneous acceleration of the parcel is the sum of these effects:
                dw/dt = b = g/(c_p * T) * (-dMSE/dz + LH)

            4. Kinetic energy per unit mass is KE = 0.5 * w^2, so its rate of change follows:
                dKE/dt = w * dw/dt

            5. Summary: at each vertical level, compute dw/dt from the MSE gradient and LH,
            then compute dKE/dt from the local w. Step these forward to evolve the profile.

            At KE = 0, we can use the mean over one timestep, ΔKE = 0.5 * (dw/dt * Δt)^2 to initialize KE.  0.5 * (dw/dt * Δt)^2 = 0.5( (g/(c_p * T) * (-dMSE/dz + LH)) * Δt )^2
        =#

        ρatke_convective_production = aux_en.tke_convective_production # for some reason these have ρa in ρatke but are named this way so we follow
 
        K_buoy = param_set.user_params.convective_tke_buoyancy_coeff
        @. ρatke_convective_production = K_buoy * ρ_c * a_en * max(instability, 0) # only positive values of convective TKE [[ not sure if existing tke should matter...]]

        # @. ρatke_convective_production = ifelse(iszero(ρatke_convective), FT(0.5) * ρ_c * a_en * K_buoy * max(instability, 0)^2 * Δt^2, sqrt(2 * ρatke_convective/(ρ_c * a_en)) * K_buoy * ρ_c * a_en * max(instability, 0)) # if no existing tke use 0.5 * (a*Δt)^2 else use w* a


        #= 
            Production should be upper bounded by buoyancy. buoyant acceleration follows B = g * (ρ_0 - ρ)/ρ_0. For us this is essentially always below 0.1

            Δw(t+Δt) = a * Δt < 0.1 * Δt ==> ΔKE = 1/2 * Δw^2 < 0.5 * (0.1 * Δt)^2 = 0.005 * Δt^2 J/kg
            For Δt = 10s, that gives  0.5
        =#
        b_max = FT(0.1) # m/s^2 # maximum buoyancy
        @. ρatke_convective_production = min(ρatke_convective_production, FT(0.5) * ρ_c * a_en * (b_max * Δt)^2)


        if any(!isfinite, ρatke_convective_production)
            @warn "Non-finite values in convective TKE production computation."
            error("K_buoy = $(K_buoy); a_en = $(a_en); instability = $(parent(instability)); tke = $(parent(aux_en.tke))")
        end

        
        # -- advection/diffusion -- #
        # ∇c_tkec = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0)))) # we have no sure surface value so we could just extrapolate. However since our sfc generates no tke and we leave that to updraft and shear (for now, could revisit if it's problematic), we'll just go w/ 0 on both (technically tke could go out the top but we enforce w=0 so...)
        # advection (is diffusive, like pressure diffusion I guess...) # [[ they use ::         @. nh_press_adv = Int(ᶠinterp_a(a_up) > 0) * ρ_f * ᶠinterp_a(a_up) * α_a * w_up * ∇(wvec(Ifc(w_up)))      ]]
        # We'll have a fast advection part and a slow diffusion part. diffusion is up gradient, advection is with w_convective [[ could in principle go both directions but should be biased up for convection so go up ]]
        ρatke_convective_advection = aux_en.tke_convective_advection # for some reason these have ρa in ρatke but are named this way so we follow
        w_convective = aux_tc.temporary_3 # 1 and 2 already taken but avaiable again
        @. w_convective = sqrt(2 * ρatke_convective/(ρ_c * a_en)) # tke = 1/2 ρ w^2  => w = sqrt(2*tke) , tke is per unit mass here
        k_adv = param_set.user_params.convective_tke_advection_coeff
        # @. tke_convective_advection = ρ_c * a_en * k_adv * w_convective * ∇c(wvec(Ifw(tke_convective)))
        LBw = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))
        RBw = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0)))
        T_freeze = TCP.T_freeze(param_set)
        RH = aux_tc.temporary_4 # 1 and 2 already taken

        @. RH = ifelse(aux_en.T < T_freeze, state.calibrate_io ?  relative_humidity_over_ice(thermo_params, aux_en.ts) : aux_en.RH_ice, state.calibrate_io ?  relative_humidity_over_liquid(thermo_params, aux_en.ts) : aux_en.RH_liq) # for enhancement of updrafts/downdrafts
        @. ρatke_convective_advection = -∇c(wvec(LBw(k_adv * (1+(RH-1)) * w_convective * ρatke_convective/2 ))) / (ρ_c * a_en/2) # negative divergence of flux [[ up -- faster if supersat]
        @. ρatke_convective_advection += -∇c(wvec(RBw(k_adv * (1-(RH-1)) * -w_convective * ρatke_convective/2 ))) / (ρ_c * a_en/2) # some going the other way for diffusion [[ down -- faster if subsat]]

        # k_diff = param_set.user_params.convective_tke_advection_coeff # rename to diffusiong if we go with this option...
        # @. ρatke_convective_advection = -∇c_tkec(-1 * ρ_f * ᶠinterp_a(a_en) * k_diss * (Ifw(w_convective))^2 * ∇f(ρatke_convective))
        # diffusive_flux_tke = -ρ_f * ᶠinterp_a(a_en) * k_diss * Ifw(w_convective)^2 * ∇f(ρatke_convective) # :: Recall, @. aux_tc_f.diffusive_flux_qt = -aux_tc_f.ρ_ae_KQ * ∇q_tot_en(wvec(aux_en.q_tot))
        # @. ρatke_convective_advection = -∇c_tkec(wvec( -1 * ρ_f * ᶠinterp_a(a_en) * k_diff * (Ifw(w_convective))^2 * ∇f(ρatke_convective) )) # (w_convective)^2 comes from KE = 1/2 ρ w^2 → w = sqrt(2*KE/ρ); flux ∼ K_turb * ∇(KE) with K_turb ∼ ℓ_mix * w ∼ w^2 scaling for convective transport [[ I think because of the ∇f we need to specify a BC? idk... the output is on F so it's not super clear to me why...]]

        
        # -- dissipation -- #
        # dissipation due to stability [[ should we use MSE for static stability? part of me also feels like we should be using moist static energy here since it's actually constant in the cloud but idk.... ]]
        ρatke_convective_dissipation = aux_en.tke_convective_dissipation # for some reason these have ρa in ρatke but are named this way so we follow
        ∂MSE∂z = aux_en.∂MSE∂z
        k_diss = param_set.user_params.convective_tke_dissipation_coeff
        K_diss = aux_tc.temporary_5
        # @. K_diss = -k_diss * ∇c(ρ_f * wvec(Ifw(min(∂MSE∂z,0)))) # if MSE decreases with height, we dissipate convective TKE
        @. K_diss = k_diss * max(∂MSE∂z, 0)^2 # if MSE increases with height, we dissipate convective TKE squared goes like brunt vaisala
        @. ρatke_convective_dissipation = K_diss * (ρatke_convective/(ρ_c * area_en))^(3/2) * ρ_c * area_en / FT(1000) # dissipation  goes as tke^(3/2)/l

        # self dissipation [[ calculated in  `compute_covariance_dissipation()`.. we apply here or else tke can build up for the entire run if it's not advected away ]]
        # @. ρatke_convective_dissipation += (aux_en_2m.tke.dissipation * 1) * ifelse.( aux_en.tke > eps(FT), aux_en.tke_convective / aux_en.tke , FT(0)) # scale by fraction of convective tke in total tke to get dissipation of convective part
        k_self_diss = param_set.user_params.convective_tke_self_dissipation_coeff
        @. ρatke_convective_dissipation += k_self_diss * ρ_c * area_en * max(aux_en.tke_convective, 0)^FT(3/2) / FT(1000) # should have mixing_length maybe 1000+ . For them, c_d was < 1 so this made dissipation smaller... unclear where our c_d lands... but it's bigger but we want less disipation I think..

        # -- limit production to enforce tke max -- #
        # limit production and net if we're exceeding tke max (smooth limiter that ramps on as we exceed tke max)

        # current KE per unit mass
        # KE_current :: ρatke_convective / (ρ_c * a_en)
        # KE increment from this timestep (production + dissipation + advection)
        # ΔKE_timestep :: (ρatke_convective_production - ρatke_convective_dissipation + ρatke_convective_advection) * Δt / (ρ_c * a_en)
        # overshoot due to this timestep’s contribution only
        # @. overage_prod = min(KE_current + ΔKE_timestep - tke_max, ρatke_convective_production * Δt/(ρ_c * a_en)) # only consider production contribution to overshoot
        # @. overage_prod = max(min(KE_current + ΔKE_timestep - tke_max, ρatke_convective_production * Δt / (ρ_c * a_en)), 0) # GPT
        overage_prod = aux_tc.temporary_4 # used for RH but avail again
        @. overage_prod = min(
            max((ρatke_convective / (ρ_c * a_en)) +  # KE_current
            (ρatke_convective_production - ρatke_convective_dissipation + ρatke_convective_advection) * Δt / (ρ_c * a_en) - # ΔKE_timestep
            tke_max, 0), max(ρatke_convective_production, 0) * Δt/(ρ_c * a_en)) # my fix, i think ρatke_convective_production is already positive

        # smooth exponential decay applied only to the overshoot fraction
        # ϵ = FT(0.1)
        ϵ = FT(4) # get rid of like 90% of overage
        @. ρatke_convective_production -= (ρ_c * a_en)/Δt * overage_prod * (1 - exp(-ϵ))
        @. ρatke_convective_production = max(ρatke_convective_production, 0) # no negative production [[just a double check safety net, should still be positive]]

        if any(!isfinite, ρatke_convective_production)
            @warn "Non-finite values in convective TKE production limiter computation."
            error("KE_current = $(parent(KE_current)); ΔKE_timestep = $(parent(ΔKE_timestep)); overage_prod = $(parent(overage_prod)); tke_max = $(parent(tke_max)); tke = $(parent(aux_en.tke))")
        end

        # turn off production at surface, surface +1, toa, toa -1 
        ρatke_convective_production[kc_surf] = 0
        ρatke_convective_production[kc_surf + 1] = 0
        ρatke_convective_production[kc_toa] = 0
        ρatke_convective_production[kc_toa - 1] = 0



        # -- combine tendencies -- #
        @. tend_ρatke_convective = ρatke_convective_production - ρatke_convective_dissipation + ρatke_convective_advection  # prolly should have some kinda limiter but...
        @. tend_ρatke_convective = max(tend_ρatke_convective, -ρatke_convective / Δt) # can't dissipate more than we have
        # we do not combine tendencies [  tend_covar += tend_ρatke_convective  ] here because we track ρatke_convective separately and we cant separate them later after adding

        # if any(!isfinite, tend_ρatke_convective)
        #     full_print("ρatke_convective = $(parent(ρatke_convective))")
        #     full_print("ρatke_convective_production = $(parent(ρatke_convective_production))")
        #     full_print("ρatke_convective_advection = $(parent(ρatke_convective_advection))")
        #     full_print("ρatke_convective_dissipation = $(parent(ρatke_convective_dissipation))")
        #     full_print("tend_ρatke_convective = $(parent(tend_ρatke_convective))")
        #     full_print("KE_current = $(parent(KE_current))")
        #     full_print("ΔKE_timestep = $(parent(ΔKE_timestep))")
        #     full_print("overage_prod = $(parent(overage_prod))")
        #     full_print("tke_max = $(parent(tke_max))")
        #     full_print("tke = $(parent(aux_en.tke))")

        #     error("Non-finite values in convective TKE tendency computation. Details: 
        #         ρatke_convective = $((ρatke_convective))
        #         ρatke_convective_production = $((ρatke_convective_production))
        #         ρatke_convective_advection = $((ρatke_convective_advection))
        #         ρatke_convective_dissipation = $((ρatke_convective_dissipation))
        #         tend_ρatke_convective = $((tend_ρatke_convective))
        #         KE_current = $((KE_current))
        #         ΔKE_timestep = $((ΔKE_timestep))
        #         overage_prod = $((overage_prod))
        #         tke_max = $((tke_max))
        #         tke = $((aux_en.tke))
        #     ")
        # end


        # if rand() < 1e-5
        #     @warn "1) Maybe switch to advection (strong) down MSE gradient, and diffusion (weaker). 2) Figure out how we want to combine this with regular ρatke. rn we combine tendencies, but idk if that's ideal for bookkeeping and ρatke and ρatke_convective could get out of sync. esp since ρatke already faces its own transport, dissipation, etc... Maybe we should just combine prog tke later, and limit prog tke to respect convective tke and vice versa..."
        #     @warn "1) at the same time this tke should be self-dissipatting as it removes the MSE gradient which is different from being damped, so maybe we need two dampings? Or just let the ρatke damping occur as normal and then our damping is only the stability damping"
        #     @warn "If we keep convective TKE separate I think we're ok as long as diffusion/advection come to balance production... it still might be good to add a self-loss mechanism so in areas with persistent MSE gradients we dont simply generate more and more tke forever but it is still generated fast enough for export upwards"
        #     @warn "The true ρatke decouples pretty quickly out of the main production region from dissipation... also it's getting double advected and stuff... so idk... the advection and dissipation from a separate convection_tke cancel out... while the real tke is already 0. i think either we just put MSE generation into real tke or we keep convective tke completely separate...."
        #     @warn "We still have a Δt hydraulic jump somewhere"
        #     @warn "Add an option to just increase K and/or mixing length in statically unstable regions (using θ_v or θ in unsat, and θ_e == MSE in saturated regions)"
        #     @warn "If we're not relying on MF alone, the massflux boost should also allow for boosts from convective TKE transport... It could probably be tied to some estimate of mixing length or w... idk. dT/dz * mixing_length = ΔT or something though that's prone to variability... idk."
        # end

    end

    return nothing
end



function update_diagnostic_covariances!(
    edmf::EDMFModel,
    state::State,
    ::Val{covar_sym},
) where {covar_sym}
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    grid = Grid(state)
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ
    aux_en_2m = center_aux_environment_2m(state)
    aux_up_f = face_aux_updrafts(state)
    aux_en = center_aux_environment(state)
    covar = getproperty(aux_en, covar_sym)
    aux_covar = getproperty(aux_en_2m, covar_sym)
    aux_up = center_aux_updrafts(state)
    w_en_f = face_aux_environment(state).w
    c_d = mixing_length_params(edmf).c_d
    covar_lim = edmf.thermo_covariance_model.covar_lim

    ρ_ae_K = face_aux_turbconv(state).ρ_ae_K
    KH = center_aux_turbconv(state).KH
    aux_tc = center_aux_turbconv(state)
    aux_bulk = center_aux_bulk(state)
    D_env = aux_tc.temporary_1
    a_bulk = aux_bulk.area
    tke_en = aux_en.tke

    shear = aux_covar.shear
    entr_gain = aux_covar.entr_gain
    rain_src = aux_covar.rain_src
    mixing_length = aux_tc.mixing_length
    min_area = edmf.minimum_area

    area_en = aux_en.area

    parent(D_env) .= 0
    @inbounds for i in 1:N_up
        turb_entr = aux_up[i].frac_turb_entr
        entr_sc = aux_up[i].entr_sc
        entr_ml = aux_up[i].entr_ml
        entr_rate_inv_s = aux_up[i].entr_rate_inv_s
        w_up = aux_up_f[i].w
        a_up = aux_up[i].area
        # TODO: using `Int(bool) *` means that NaNs can propagate
        # into the solution. Could we somehow call `ifelse` instead?
        if edmf.entrainment_type isa FractionalEntrModel
            @. D_env += Int(a_up > min_area) * ρ_c * a_up * Ic(w_up) * (entr_sc + entr_ml + turb_entr)
        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. D_env += Int(a_up > min_area) * ρ_c * a_up * entr_rate_inv_s
        end
    end

    @. covar =
        (shear + entr_gain + rain_src) /
        max(D_env + ρ_c * area_en * c_d * sqrt(max(tke_en, 0)) / max(mixing_length, 1), covar_lim)
    return nothing
end


function GMV_third_m(
    # edmf::EDMFModel,
    N_up::Int,
    state::State,
    ::Val{covar_en_sym},
    ::Val{var},
    ::Val{gm_third_m_sym},
) where {covar_en_sym, var, gm_third_m_sym}
    grid = Grid(state)
    # N_up = n_updrafts(edmf)
    gm_third_m = getproperty(center_aux_grid_mean(state), gm_third_m_sym)
    kc_surf = kc_surface(grid)
    FT = float_type(state)

    aux_bulk = center_aux_bulk(state)
    aux_up_f = face_aux_updrafts(state)
    is_tke = covar_en_sym === :tke
    aux_en_c = center_aux_environment(state)
    covar_en = getproperty(aux_en_c, covar_en_sym)
    aux_en_f = face_aux_environment(state)
    aux_en = is_tke ? aux_en_f : aux_en_c
    aux_up_c = center_aux_updrafts(state)
    aux_tc = center_aux_turbconv(state)
    ϕ_gm = aux_tc.ϕ_gm
    ϕ_gm_cov = aux_tc.ϕ_gm_cov
    ϕ_en_cov = aux_tc.ϕ_en_cov
    ϕ_up_cubed = aux_tc.ϕ_up_cubed
    aux_up = is_tke ? aux_up_f : aux_up_c
    var_en = getproperty(aux_en, var)
    area_en = aux_en_c.area
    Ic = is_tke ? CCO.InterpolateF2C() : x -> x
    wvec = CC.Geometry.WVector
    # ∇c = CCO.DivergenceF2C()
    w_en = aux_en_f.w

    @. ϕ_gm = area_en * Ic(var_en)
    @inbounds for i in 1:N_up
        a_up = aux_up_c[i].area
        var_up = getproperty(aux_up[i], var)
        @. ϕ_gm += a_up * Ic(var_up)
    end

    # w'w' ≈ 2/3 TKE (isotropic turbulence assumption)
    if is_tke
        @. ϕ_en_cov = FT(2 / 3) * covar_en
    else
        @. ϕ_en_cov = covar_en
    end

    parent(ϕ_up_cubed) .= 0
    @. ϕ_gm_cov = area_en * (ϕ_en_cov + (Ic(var_en) - ϕ_gm)^2)
    @inbounds for i in 1:N_up
        a_up = aux_up_c[i].area
        var_up = getproperty(aux_up[i], var)
        @. ϕ_gm_cov += a_up * (Ic(var_up) - ϕ_gm)^2
        @. ϕ_up_cubed += a_up * Ic(var_up)^3
    end

    @. gm_third_m = ϕ_up_cubed + area_en * (Ic(var_en)^3 + 3 * Ic(var_en) * ϕ_en_cov) - ϕ_gm^3 - 3 * ϕ_gm_cov * ϕ_gm

    gm_third_m[kc_surf] = 0
    return nothing
end



