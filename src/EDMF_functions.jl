
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
    Δt::FT,
    cfl_limit::FT,
    use_fallback_tendency_limiters::Bool,
) where {FT}
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
    error("Second order flux correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")
end

function apply_second_order_flux_correction(flux::CC.Fields.Field, w::CC.Fields.Field, Δt::FT; apply_diffusive_damping::Bool = false, correction_limit_factor::FT = FT(Inf)) where{FT}
    error("Second order flux correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")

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
    # If = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # I feel like we should be using left and right biased...
    If = Ifx
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

    # if edmf.area_partition_model.apply_second_order_flux_correction
        # seclf = edmf.area_partition_model.second_order_correction_limit_factor
    # end

    # ====================================== # # en is downdraft only so RBF
    ᶠinterp_a_RBF_a_en = ᶠinterp_a_RBF(a_en)
    #
    IfRBF_θ_liq_ice_en = IfRBF(θ_liq_ice_en)
    IfRBF_q_tot_en = IfRBF(q_tot_en)
    if edmf.moisture_model isa NonEquilibriumMoisture
        IfRBF_q_liq_en = IfRBF(q_liq_en)
        IfRBF_q_ice_en = IfRBF(q_ice_en)
        IfLBF_q_liq_en = IfLBF(q_liq_en)
        IfLBF_q_ice_en = IfLBF(q_ice_en)
    end
    # ====================================== #

    #=
        Ideally we'd use LBF for up and RBF for down to match the updraft fluxes but that doesnt work for a few reasons:
            1. That offset induces strong vertical motion and turbulence that isn't real. It becomes very hard to maintain quiescent air.
            It's not `unstable` per se, but you become susceptible to the tiniest of perturbations growing into periodic updrafts.
            This only gets worse at coarser resolutions when the difference between LBF and RBF is larger.
        
        The downside becomes that very strong downdrafts are untenable -- since the environment is not right biased, it is ripe for odd-even decoupling (checkerboarding).
        I'm not sure if there's an easy way around this that doesnt involve having at least 1 coherent downdraft.

        Thus, everything here uses co-located fluxes (no LBF or RBF), meaning the induced tendencies in env will not match this exactly but will be controlled somewhat by what LBF in compute_up_tendencies!() does for the updrafts.
    =#

    if (edmf.convective_tke_handler isa ConvectiveTKE)
        w_convective = aux_tc.temporary_3
        @. w_convective = sqrt(2 * max(aux_en.tke_convective, FT(0))) # convective velocity scale

        # 1. Define Operator (Same as before)
        adv_bcs = (; bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))
        ∇f = CCO.GradientC2F(; adv_bcs...)

        ℓ_mix = FT(500) # mixing length scale for convective fluxes -- could be edmf dependent
    end

    if edmf.area_partition_model isa StandardAreaPartitionModel

        if edmf.area_partition_model.apply_second_order_flux_correction
            error("second order correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")
        else
            @. massflux_en = ρ_f * ᶠinterp_a_RBF_a_en(a_en) * w_en

            if (edmf.convective_tke_handler isa ConvectiveTKE) && (edmf.convective_tke_handler.transport_conserved_by_advection || edmf.convective_tke_handler.transport_condensed_by_advection)
                massflux_en_conv = aux_tc_f.temporary_f1
                @. massflux_en_conv = ρ_f * ᶠinterp_a_RBF_a_en(a_en)/2 * Ifx(w_convective) # feels mostly the tke scale, subsidence happens elsewhere dry
            end


            if !((edmf.convective_tke_handler isa ConvectiveTKE) && (edmf.convective_tke_handler.transport_conserved_by_advection))
                @. massflux_h = massflux_en * IfRBF_θ_liq_ice_en(θ_liq_ice_en)
                @. massflux_qt = massflux_en * IfRBF_q_tot_en(q_tot_en)
            else # :: Use advection for conserved variables ::


                # -=- Simple DownGradient Perturbation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

                # # 2. Calculate Flux with 2nd Order Correction
                # efficiency = FT(0.25)
                # perturbation = aux_tc_f.temporary_f2
                # @. perturbation = clamp(-1 * ℓ_mix * toscalar(wvec(∇f(θ_liq_ice_en))), -sqrt(max(If(aux_en.Hvar), 0)), sqrt(max(If(aux_en.Hvar), 0))) # Value: -1 * L * Gradient_at_Face | Min/Max: Interp Center -> Face

                # @. massflux_h =  # Logic: Flux = (Stream 1: M * p) - (Stream 2: M * -p) = 2 * M * p
                #     massflux_en * If(θ_liq_ice_en) +                 # Mean Transport
                #     (2 * massflux_en_conv * efficiency) * perturbation # Turbulent Mixing (Downgradient)

                # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=- #

                #=
                Derivation (θ_liq_ice ≡ h):

                We model the unresolved convective turbulent flux of h as an exchange between two
                subgrid streams with equal-and-opposite perturbations ±h′ carried by a convective
                exchange mass flux massflux_en_conv:

                F_h,turb ≈ (massflux_en_conv * (+h′)) + (massflux_en_conv * (−h′))
                        = 2 * massflux_en_conv * h′

                We choose h′ from a local mixing-length estimate, with a variance-based bound:

                h′_ML  = − ℓ_mix * ∂z h   (evaluated at faces)
                |h′| ≤ √Hvar             (variance gives a physically admissible amplitude)

                At a stable inversion, vertical motions cannot efficiently exchange properties across
                faces. We represent that by multiplying the turbulent exchange by a transmission factor:

                transmission_u = KE / (KE + PE_barrier)
                KE ~ 0.5 * w_base^2, with w_base^2 ~ 2 * tke_convective
                PE_barrier ~ (g/θ_virt) * Δθ_virt * ℓ_mix   (stable jump only)

                So the final flux is:

                massflux_h = massflux_en * If(h_en) + 2 * massflux_en_conv * efficiency
                            * transmission_u * clamp(h′_ML, ±√Hvar)

                This is (i) variance-bounded, (ii) downgradient in the absence of inversion blocking,
                and (iii) shuts down across strong inversions in a physically controlled way.
                =#

                efficiency = FT(0.25)
                g = TCP.grav(param_set)

                LBF = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0)))
                RBF = CCO.RightBiasedC2F(; top    = CCO.SetValue(FT(0)))

                w_base         = aux_tc_f.temporary_f1   # face
                transmission_u = aux_tc_f.temporary_f2   # face
                barrier_energy = aux_tc_f.temporary_f3   # face
                perturbation   = aux_tc_f.temporary_f4   # face

                # velocity scale from convective TKE (face-consistent)
                @. w_base = If(w_convective)

                # inversion barrier proxy (stable jump only)
                @. barrier_energy = max(
                    (g / If(aux_en.θ_virt)) * (RBF(aux_en.θ_virt) - LBF(aux_en.θ_virt)) * ℓ_mix,
                    FT(0)
                )

                # transmission KE/(KE+PE)
                @. transmission_u =
                    (FT(0.5) * w_base^2) / (FT(0.5) * w_base^2 + barrier_energy + FT(1e-10))

                # mixing-length perturbation, variance-clipped
                @. perturbation = clamp(
                    -ℓ_mix * toscalar(wvec(∇f(θ_liq_ice_en))),
                    -sqrt(max(If(aux_en.Hvar), FT(0))),
                    sqrt(max(If(aux_en.Hvar), FT(0)))
                )

                # mean + turbulent exchange (blocked by inversion)
                @. massflux_h =
                    massflux_en * If(θ_liq_ice_en) +
                    (FT(2) * massflux_en_conv * efficiency) * (transmission_u * perturbation)

                # -=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

                if !(edmf.moisture_model isa NonEquilibriumMoisture) || ((edmf.moisture_model isa NonEquilibriumMoisture) && !edmf.convective_tke_handler.transport_conserved_by_advection) # if Eq or, despite having convective tke, not using advection for conserved quantities
                    # This is biased towards stability -- moister air goes up, drier air goes down
                    # @. massflux_qt = (massflux_en_conv) * If(q_tot_en + (2*frac_supersat - 1)*sqrt(aux_en.QTvar)) + (massflux_en - massflux_en_conv) * If(q_tot_en - (2*frac_supersat - 1)*sqrt(aux_en.QTvar))
                    @. massflux_qt = (massflux_en_conv * efficiency) * If(q_tot_en + perturbation) + (massflux_en - massflux_en_conv * efficiency) * If(q_tot_en - perturbation)
                end
            
            end

            if edmf.moisture_model isa NonEquilibriumMoisture

                if !((edmf.convective_tke_handler isa ConvectiveTKE) && edmf.convective_tke_handler.transport_condensed_by_advection)
                    @. massflux_ql = massflux_en * IfRBF_q_liq_en(q_liq_en)
                    @. massflux_qi = massflux_en * IfRBF_q_ice_en(q_ice_en)

                    @. ql_flux_vert_adv = massflux_ql # same here
                    @. qi_flux_vert_adv = massflux_qi # same here
                else
                    # ql and qi should be correlated w/ tke. We handle prt of that in diffusivity but they shouldn't feel the netire advective force
                    @. massflux_ql = (massflux_en_conv * efficiency) * IfLBF_q_liq_en(q_liq_en) + (massflux_en) * IfRBF_q_liq_en(q_liq_en) # assume liq/ice are in the upper part of the distribution, and the down side has little flux. (massflux_en is small anyway, and significant conv downdrafts should dry quickly)
                    @. massflux_qi = (massflux_en_conv * efficiency) * IfLBF_q_ice_en(q_ice_en) + (massflux_en) * IfRBF_q_ice_en(q_ice_en) # assume liq/ice are in the upper part of the distribution, and the down side has little flux. (massflux_en is small anyway, and significant conv downdrafts should dry quickly)
                    

                    @. ql_flux_vert_adv = massflux_ql # same here
                    @. qi_flux_vert_adv = massflux_qi # same here

                    # since we introduce bias in ql, qi transport, we need to adjust qt transport accordingly and calculate each term separately

                    if edmf.convective_tke_handler.transport_conserved_by_advection  # if qt was also transported by advection
                        # @. massflux_qt = massflux_ql + massflux_qi + 
                        #     (massflux_en/2 + massflux_en_conv) * If((q_tot_en - q_liq_en - q_ice_en) + (2*frac_supersat - 1)*sqrt(aux_en.QTvar)) +
                        #     (massflux_en/2 - massflux_en_conv) * If((q_tot_en - q_liq_en - q_ice_en) - (2*frac_supersat - 1)*sqrt(aux_en.QTvar))

                        # 2. Calculate Scalar Perturbation (phi')
                        # This is the ANOMALY (e.g. +2 degrees), not the flux.


                        # perturbation_old = aux_tc.temporary_6
                        # ℓ_mix = FT(500)
                        # @. perturbation_old = clamp(-1 * ℓ_mix * (∇c(wvec(If(q_tot_en - q_liq_en - q_ice_en)))), -sqrt(max(aux_en.QTvar, 0)), sqrt(max(aux_en.QTvar, 0)))

                        perturbation = aux_tc_f.temporary_f2
                        @. perturbation = clamp(-1 * ℓ_mix * toscalar(wvec(∇f(q_tot_en - q_liq_en - q_ice_en))), -sqrt(max(If(aux_en.QTvar), 0)), sqrt(max(If(aux_en.QTvar), 0))) # Value: -1 * L * Gradient_at_Face | Min/Max: Interp Center -> Face

                        efficiency = FT(0.25)
                        # @. massflux_qt = massflux_ql + massflux_qi + 
                        #     # --- Stream 1, 2 Advection ---
                        #     (massflux_en_conv * efficiency) * If((q_tot_en - q_liq_en - q_ice_en) + perturbation_old) +
                        #     (massflux_en - massflux_en_conv * efficiency) * If((q_tot_en - q_liq_en - q_ice_en) - perturbation_old) #-
                        #     # # --- Second Order Correction (Fused) --- [ not needed if we can generate variances correctly ]
                        #     # (
                        #     #     0.5 * Δt * (1 / ρ_f) * (
                        #     #         (massflux_en_conv)^2 + 
                        #     #         (massflux_en - massflux_en_conv)^2
                        #     #     ) * toscalar(wvec(∇f(q_tot_en - q_liq_en - q_ice_en))) 
                            # )


                        # 3. Calculate Flux
                        # Now we use the clamped perturbation directly. 
                        # Logic: Flux = (Stream 1: M * p) - (Stream 2: M * -p) = 2 * M * p
                        @. massflux_qt = massflux_ql + massflux_qi + 
                            massflux_en * If(q_tot_en - q_liq_en - q_ice_en) +                 # Mean Transport
                            (2 * massflux_en_conv * efficiency) * perturbation # Turbulent Mixing

                        # @. perturbation = clamp(-1 * ℓ_mix * toscalar(wvec(∇f(q_tot_en))), -sqrt(max(If(aux_en.QTvar), 0)), sqrt(max(If(aux_en.QTvar), 0))) # Value: -1 * L * Gradient_at_Face | Min/Max: Interp Center -> Face
                        # @. massflux_qt = massflux_en * If(q_tot_en) +                 # Mean Transport
                        #     (2 * massflux_en_conv * efficiency) * perturbation # Turbulent Mixing

                        # add some very small diffusion
                        # K_floor = FT(1.0)
                        # @. massflux_qt -= (  
                        #     ρ_f * FT(1e-1) * (If(w_convective .* ℓ_mix) + K_floor) * toscalar(wvec(∇f(q_tot_en)))         
                        # )
                
                    end


                    #=
                    Derivation (q_t):

                    Write the total-water flux as mean transport plus an unresolved convective turbulent
                    exchange term:

                    massflux_qt = mean transport of qt  +  F_qt,turb

                    Mean transport in your code is already assembled as:
                    massflux_ql + massflux_qi + massflux_en * If(q_v)
                    with q_v = q_tot_en - q_liq_en - q_ice_en,
                    so that the sum represents transport of total qt.

                    For the turbulent part we use the same exchange-mass-flux form as for h:
                    F_qt,turb ≈ 2 * massflux_en_conv * efficiency * q′

                    The key difference from h is that at cloud top the LES signature is:
                    QTvar spikes, but QTOFLXR → 0 or slightly negative.
                    That requires that *as the inversion strengthens*, the (upward) exchange shuts down
                    while a (downward) entrainment return can remain.

                    We implement this with a *single bounded perturbation q′*:

                    q′ = transmission_u * q′_ML  −  entrainment_ratio * (1 − transmission_u) * σ_q

                    where
                    σ_q = √QTvar  (amplitude of moist–dry anomalies),
                    q′_ML = −ℓ_mix * ∂z q_t (face mixing-length estimate, signed),
                    and the same transmission_u as in the h block is used.

                    Properties:
                    - When inversion is weak: transmission_u ≈ 1 ⇒ q′ ≈ q′_ML (stable, downgradient core)
                    - When inversion is strong: transmission_u ≈ 0 ⇒ q′ ≈ −entrainment_ratio * σ_q (can drive
                    small negative flux even if ∂z qt < 0), matching LES “flux collapse / slight negative”
                    - Always variance-bounded: |q′_ML| ≤ σ_q
                    =#

                    # ======================================================================== #
                    # :: QT MASS FLUX (mean + derived 2-stream SGS exchange, barrier-limited) ::
                    #    Barrier is grid-independent: E_b = 0.5 * N2 * ℓ_mix^2
                    # ======================================================================== #

                    # g     = TCP.grav(param_set)
                    # ℓ_mix = FT(500)

                    # aux_tc_f = face_aux_turbconv(state)

                    # # Face temporaries (reuse; no allocations)
                    # w_base         = aux_tc_f.temporary_f1   # face: base velocity scale
                    # barrier_energy = aux_tc_f.temporary_f2   # face: E_b
                    # w_trans        = aux_tc_f.temporary_f3   # face: transmitted velocity (then factor)
                    # tmp_f          = aux_tc_f.temporary_f4   # face: σ_q then q′

                    # # --- A) Base convective velocity scale at faces ---
                    # @. w_base = max(If(w_convective), FT(0))  # face-interpolated w_convective

                    # # --- B) Grid-independent barrier: E_b = 0.5 * N2 * ℓ_mix^2, with N2 = (g/θv) * ∂z θv ---
                    # @. barrier_energy = max(
                    #     FT(0.5) * ℓ_mix^2 *
                    #     (g / If(aux_en.θ_virt)) *
                    #     toscalar(wvec(∇f(aux_en.θ_virt))),
                    #     FT(0)
                    # )

                    # # --- C) Transmitted velocity from KE balance: 0.5 w_base^2 -> 0.5 w_trans^2 + E_b ---
                    # @. w_trans = sqrt(max(w_base^2 - FT(2) * barrier_energy, FT(0)))

                    # # --- D) q′ from mixing-length displacement, bounded by σ_q = √QTvar ---
                    # @. tmp_f = sqrt(max(If(aux_en.QTvar), FT(0)))  # tmp_f := σ_q

                    # @. tmp_f = clamp(
                    #     -ℓ_mix * toscalar(wvec(∇f(q_tot_en))),   # q′_ML
                    #     -tmp_f,                                  # -σ_q
                    #     tmp_f                                   # +σ_q
                    # )

                    # # --- E) Scale massflux_en_conv by transmitted/base velocity (since your M_u ∝ w_base) ---
                    # @. w_trans = w_trans / (w_base + FT(1e-10))    # now w_trans := factor in [0,1]

                    # # --- F) Total qt mass flux ---
                    # @. massflux_qt =
                    #     massflux_ql + massflux_qi +
                    #     massflux_en * If(q_tot_en - q_liq_en - q_ice_en) +
                    #     (FT(2) * massflux_en_conv * w_trans) * tmp_f




                    
                
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
            error("second order correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")
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
                error("second order correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")
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
            error("second order correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")
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
                error("second order correction was deprecated for now since it's not trivial to implement performantly and doesn't help that much")
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
    # Compute the mass flux tendencies
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
    prog_pr = center_prog_precipitation(state)
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


    if edmf.moisture_model isa NonEquilibriumMoisture
        aeKQq_liq_bc = FT(0)
        aeKQq_ice_bc = FT(0)
        ∇q_liq_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_liq_bc), top = CCO.SetDivergence(FT(0)))
        ∇q_ice_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_ice_bc), top = CCO.SetDivergence(FT(0)))
    end

    aeKQq_rai_bc = FT(0) # rain
    aeKQq_sno_bc = FT(0) # snow
    ∇q_rai_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_rai_bc), top = CCO.SetDivergence(FT(0)))
    ∇q_sno_en = CCO.DivergenceC2F(; bottom = CCO.SetDivergence(aeKQq_sno_bc), top = CCO.SetDivergence(FT(0)))


    #=
        If we have convective TKE
            - If we are doing TKE transport by advection, we just need to remove the convective TKE from the total TKE in the eddy diffusivity calculation and then proceed normally.
            - If we are not (so we use diffusion)
                - check if we have a different scaling factor for convective TKE eddy diffusivity
                    - if we do, then we need to handle the convective and regular TKE eddy diffusivities separately
                    - if we don't, we should just proceed with the full TKE as normal
        If we don't have convective TKE
            - proceed as normal
    =#


    use_separate_tke_conserved::Bool = (edmf.convective_tke_handler isa ConvectiveTKE) && ((edmf.convective_tke_handler.transport_conserved_by_advection) || !isone(edmf.convective_tke_handler.ed_scaling_factor))
    use_separate_tke_condensed::Bool = (edmf.convective_tke_handler isa ConvectiveTKE) && ((edmf.convective_tke_handler.transport_condensed_by_advection) || !isone(edmf.convective_tke_handler.ed_scaling_factor))

    # Prep if condensed
    if (use_separate_tke_conserved || use_separate_tke_condensed)
        aux_tc = center_aux_turbconv(state)
        
        # Remove the convective TKE

        f_c_m::FT = FT(edmf.convective_tke_handler.ed_scaling_factor) # adjustment to c_m, tke_ed_coeff
        f_K_tke = aux_tc_f.temporary_f1
        @. f_K_tke = FT(1)      # looking for bugs

        #=
            Now this is where all our biasing, saturation efficiencies, etc come into play

            See in update_aux():
                KM[k] = c_m * ml.mixing_length * sqrt(max(aux_en.tke[k], 0))
                KH[k] = KM[k] / aux_tc.prandtl_nvec[k]
                KQ[k] = KH[k] / Le

            Thus, we want to make some adjustments to our convective part.
        =#

        c_m = mixing_length_params(edmf).c_m
        Le = mixing_length_params(edmf).Le
        N = aux_tc.temporary_4
        @. N = sqrt(max(aux_tc.∂b∂z, FT(1e-6))) # buoyancy frequency, w/N is our convective overshoot mixing length, /2 for drag
        KMconv = aux_tc.temporary_5
        @. KMconv = f_c_m * c_m * ifelse(aux_en.∂MSE∂z < 0,  FT(200), FT(sqrt(2 * aux_en.tke_convective) /(2*max(N, FT(1e-6))))) * sqrt(max(aux_en.tke_convective, FT(0))) # rough convective tke eddy diffusivity for convective tke portion only, sqrt(w_tke/N)
        # KH ::KM / aux_tc.prandtl_nvec # In unstable regimes Prandtl should just be Prandtl_0 in this region... (See `turbulent_Prandtl_number()`)
        # Pr_n < 1 is unstable, we want to default to regular if ∂MSE/∂z > 0, and some enhanced value if ∂MSE/∂z < 0 :: See `/src/turbulence_functions.jl`
        Prconv = aux_tc.temporary_6
        # @. Prconv = ifelse(aux_en.∂MSE∂z < 0, mixing_length_params(edmf).Pr_n * FT(1), aux_tc.prandtl_nvec) # reduce prandtl in unstable regions to enhance KH [maybe using nvec isn't right idk... coherent updrfts are still more penetrative but I guss that what's f_c_m was for]
        @. Prconv = ifelse(aux_en.∂MSE∂z < 0, mixing_length_params(edmf).Pr_n * FT(1), FT(Inf)) # I think going all the way to absolute stability is a better plan... TKE at that point stops transporting (up and down fluxes match)? It's the closest we can get to turning into an up-gradient flux

        IfKMconv = CCO.InterpolateC2F(; bottom = CCO.SetValue(a_en[kc_surf] * KMconv[kc_surf]), top = CCO.SetValue(a_en[kc_toa] * KMconv[kc_toa]))
        IfKHconv = CCO.InterpolateC2F(; bottom = CCO.SetValue(a_en[kc_surf] * KMconv[kc_surf]/Prconv[kc_surf]), top = CCO.SetValue(a_en[kc_toa] * KMconv[kc_toa]/Prconv[kc_toa]))
        IfKQconv = CCO.InterpolateC2F(; bottom = CCO.SetValue(a_en[kc_surf] * KMconv[kc_surf]/(Prconv[kc_surf]*Le)), top = CCO.SetValue(a_en[kc_toa] * KMconv[kc_toa]/(Prconv[kc_toa]*Le)))

        #= 
            Sat efficiency for enhanced transport when we cross saturation boundaries for condensate....
            I guess this should be predicated on supersaturation and the condensate having the same gradient but oh well lol.
        =#
        peak_enhancement = FT(2.0) # max enhancement factor for saturated conditions
        sat_efficiency = aux_tc_f.temporary_f2
        @. sat_efficiency = max( Ifx(one(FT) + (peak_enhancement - one(FT)) * (FT(4) * aux_en.frac_supersat * (one(FT) - aux_en.frac_supersat))), one(FT))

        # This one-line version implements the new ramp logic inside the Ifx interpolator
        # @. sat_efficiency = max(Ifx(one(FT) + (peak_enhancement - one(FT)) * clamp((aux_en.frac_supersat - f_crit) / (one(FT) - f_crit), zero(FT), one(FT))), one(FT)) # This one has `1+` ramp up`

    end

    # == Conserved species == #
    if use_separate_tke_conserved

        if edmf.convective_tke_handler.transport_conserved_by_advection # we just need to remove convective tke
            @. aux_tc_f.diffusive_flux_qt = -aux_tc_f.ρ_ae_KQ * f_K_tke * ∇q_tot_en(wvec(aux_en.q_tot))
            @. aux_tc_f.diffusive_flux_h = -aux_tc_f.ρ_ae_KH * f_K_tke * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice))
            @. aux_tc_f.diffusive_flux_uₕ = -aux_tc_f.ρ_ae_KM * f_K_tke * ∇uₕ_gm(prog_gm_uₕ)
        else # diffusion [ # we also need to add in convective tke contributions separately to this diffusion ]
            ∇q_tot_en_tke_conv = ∇q_tot_en
            ∇θ_liq_ice_en_tke_conv = ∇θ_liq_ice_en
            ∇uₕ_gm_tke_conv = ∇uₕ_gm        
            
            # === Diffusive flux qt, h, uₕ === # [ We do 2 way biased advection for stability attm ]
            ∇qc = CCO.DivergenceF2C(; bottom = CCO.SetDivergence(FT(0)), top = CCO.SetDivergence(FT(0))) # placeholder


            qt_var_boost = aux_tc_f.temporary_f3
            @. qt_var_boost = clamp(Ifx(( 1 + sqrt(aux_en.QTvar)/aux_en.q_tot)/(1 - sqrt(aux_en.QTvar)/aux_en.q_tot)), one(FT), FT(3)) # bound by (1 + 0.5) / (1 - 0.5) = 3 
            @. aux_tc_f.diffusive_flux_qt = -((aux_tc_f.ρ_ae_KQ * f_K_tke)) * ∇q_tot_en(wvec(aux_en.q_tot)) -
            #  (f_c_m * qt_var_boost * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_tot_en_tke_conv(wvec(aux_en.q_tot))
                # (qt_var_boost * sat_efficiency * ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_tot_en_tke_conv(wvec(aux_en.q_tot)) # unstable
                (qt_var_boost * sat_efficiency * Ifx(ρ_c * a_en * KMconv/(Prconv*Le) * ∇qc(wvec(Ifx(aux_en.q_tot))))) # stable

            h_var_boost = aux_tc_f.temporary_f3
            @. h_var_boost = clamp(Ifx((1 + sqrt(aux_en.Hvar)/aux_en.θ_liq_ice)/(1 - sqrt(aux_en.Hvar)/aux_en.θ_liq_ice)), one(FT), FT(3)) # bound by (1 + 0.5) / (1 - 0.5) = 3
            @. aux_tc_f.diffusive_flux_h = -((aux_tc_f.ρ_ae_KH * f_K_tke)) * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice)) -
            #  (f_c_m * h_var_boost * aux_tc_f.ρ_ae_KH * f_k_tke_conv) * ∇θ_liq_ice_en_tke_conv(wvec(aux_en.θ_liq_ice))
                (h_var_boost * sat_efficiency * ρ_f * IfKHconv(a_en * KMconv/Prconv)) * ∇θ_liq_ice_en_tke_conv(wvec(aux_en.θ_liq_ice))

            @. aux_tc_f.diffusive_flux_uₕ = -((aux_tc_f.ρ_ae_KM * f_K_tke)) * ∇uₕ_gm(prog_gm_uₕ) - 
                # (f_c_m * aux_tc_f.ρ_ae_KM * f_k_tke_conv) * ∇uₕ_gm_tke_conv(prog_gm_uₕ)
                (ρ_f * IfKMconv(a_en * KMconv)) * ∇uₕ_gm_tke_conv(prog_gm_uₕ)

        end
    else # ::: Normal :::
        # @. aux_tc_f.diffusive_flux_qt = -aux_tc_f.ρ_ae_KQ * ∇q_tot_en(wvec(aux_en.q_tot))
        # @. aux_tc_f.diffusive_flux_h = -aux_tc_f.ρ_ae_KH * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice))
        corr_w_qt = FT(+1.0) # correlation between w and qt fluctuations
        corr_w_h  = FT(-1.0) # correlation between w and h fluctuations
        # We assume QTvar includes the induced correlations by edies. if QTvar is smaller than that would imply, then we assume there's some upgradient process.
        # Flux = w'q' which we assume to be K -∂q/∂z where eddies are inducing transport in accordance with displacements over mixing length l_mix with background gradient ∂q/∂z.
        # At the same time we prognose/diagnose some existing QTvar. if that exceeds what eddies would induce, then our flux is not high enough, if it does not exceed, then our flux is too high and there's some upgradient process going on..
        # σ_explained should be propotional to K / σ_w * ∂q/∂z = K/sqrt(2 * tke)/(∂q/∂z), where tke ≈ 1/2 w'w' = 1/2 σ_w^2, but K / σ_w should be roughly l_mix.

        @. aux_tc_f.diffusive_flux_qt = -aux_tc_f.ρ_ae_KQ * ∇q_tot_en(wvec(aux_en.q_tot))
        # correction but dont allow sign change from down gradient
        @. aux_tc_f.diffusive_flux_qt += ifelse(aux_tc_f.diffusive_flux_qt > 0,
            max(corr_w_qt * Ifx(ρ_c * a_en * sqrt(aux_en.tke) * sqrt(max(aux_en.QTvar, zero(FT)))), -aux_tc_f.diffusive_flux_qt),
            min(corr_w_qt * Ifx(ρ_c * a_en * sqrt(aux_en.tke) * sqrt(max(aux_en.QTvar, zero(FT)))), -aux_tc_f.diffusive_flux_qt)
        )

        @. aux_tc_f.diffusive_flux_h = -aux_tc_f.ρ_ae_KH * ∇θ_liq_ice_en(wvec(aux_en.θ_liq_ice))
        # correction but dont allow sign change from down gradient
        @. aux_tc_f.diffusive_flux_h += ifelse(aux_tc_f.diffusive_flux_h > 0,
            max(corr_w_h * Ifx(ρ_c * a_en * sqrt(aux_en.tke) * sqrt(max(aux_en.Hvar, zero(FT)))), -aux_tc_f.diffusive_flux_h),
            min(corr_w_h * Ifx(ρ_c * a_en * sqrt(aux_en.tke) * sqrt(max(aux_en.Hvar, zero(FT)))), -aux_tc_f.diffusive_flux_h)
        )   

        water_advection_factor = mixing_length_params(edmf).c_KTKEqt
        @. aux_tc_f.diffusive_flux_qt += ρ_f * Ifx(a_en * sqrt(aux_en.tke) * water_advection_factor * (2 * sqrt(max(aux_en.QTvar, zero(FT))))) # (mean + σ) - (mean - σ) = 2σ

        h_advection_factor = mixing_length_params(edmf).c_KTKEh
        @. aux_tc_f.diffusive_flux_h += ρ_f * Ifx(a_en * sqrt(aux_en.tke) * h_advection_factor * (2 * sqrt(max(aux_en.Hvar, zero(FT))))) # (mean + σ) - (mean - σ) = 2σ


        @. aux_tc_f.diffusive_flux_uₕ = -aux_tc_f.ρ_ae_KM * ∇uₕ_gm(prog_gm_uₕ)
    end

    # == Condensed species == #
    if use_separate_tke_condensed

        if  edmf.convective_tke_handler.transport_condensed_by_advection # we just need to remove convective tke
            # ∇qc = CCO.DivergenceF2C(; bottom = CCO.SetDivergence(FT(0)), top = CCO.SetDivergence(FT(0))) # placeholder
            # @. aux_tc_f.diffusive_flux_qt = Ifx(-Ic(aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇qc(wvec(Ifx(aux_en.q_tot)))) # convective tke removed from eddy diffusivity

            if edmf.moisture_model isa NonEquilibriumMoisture
                c_KQl = mixing_length_params(edmf).c_KQl
                c_KQi = mixing_length_params(edmf).c_KQi
                @. aux_tc_f.diffusive_flux_ql = -(c_KQl * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_liq_en(wvec(aux_en.q_liq))
                @. aux_tc_f.diffusive_flux_qi = -(c_KQi * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_ice_en(wvec(aux_en.q_ice))
            end
            c_KQr = mixing_length_params(edmf).c_KQr
            c_KQs = mixing_length_params(edmf).c_KQs
            @. aux_tc_f.diffusive_flux_qr = -(c_KQr * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_rai_en(wvec(prog_pr.q_rai)) # add diffusive (SGS) flux for rain
            @. aux_tc_f.diffusive_flux_qs = -(c_KQs * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_sno_en(wvec(prog_pr.q_sno)) # add diffusive (SGS) flux for snow
        else # diffusion [ # we also need to add in convective tke contributions separately to this diffusion ]

            if edmf.moisture_model isa NonEquilibriumMoisture
                ∇q_liq_en_tke_conv = ∇q_liq_en # no change bc they already are using aeKQq_liq_bc = 0
                ∇q_ice_en_tke_conv = ∇q_ice_en # no change bc they already are using aeKQq_ice_bc = 0
            end
            ∇q_rai_en_tke_conv = ∇q_rai_en # no change bc they already are using aeKQq_rai_bc = 0
            ∇q_sno_en_tke_conv = ∇q_sno_en # no change bc they already are using aeKQq_sno_bc = 0

            # === Diffusive flux ql, qi === #
            if edmf.moisture_model isa NonEquilibriumMoisture
                c_KQl = mixing_length_params(edmf).c_KQl
                c_KQi = mixing_length_params(edmf).c_KQi

                @. aux_tc_f.diffusive_flux_ql = -((c_KQl * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_liq_en(wvec(aux_en.q_liq))) -
                #  (f_c_m * sat_efficiency * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_liq_en_tke_conv(wvec(aux_en.q_liq))
                    (sat_efficiency * ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_liq_en_tke_conv(wvec(aux_en.q_liq))

                @. aux_tc_f.diffusive_flux_qi = -((aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_ice_en(wvec(aux_en.q_ice))) - 
                # (f_c_m * sat_efficiency * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_ice_en_tke_conv(wvec(aux_en.q_ice))
                    (sat_efficiency * ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_ice_en_tke_conv(wvec(aux_en.q_ice))
            end


            # === Diffusive flux qr, qs === #
            c_KQr = mixing_length_params(edmf).c_KQr
            c_KQs = mixing_length_params(edmf).c_KQs
            @. aux_tc_f.diffusive_flux_qr = -((c_KQr * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_rai_en(wvec(prog_pr.q_rai))) -
            #  (f_c_m * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_rai_en_tke_conv(wvec(prog_pr.q_rai))
                (ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_rai_en_tke_conv(wvec(prog_pr.q_rai))

            ∇qc = CCO.DivergenceF2C(; bottom = CCO.SetDivergence(FT(0)), top = CCO.SetDivergence(FT(0))) # placeholder
            @. aux_tc_f.diffusive_flux_qs = -((c_KQs * aux_tc_f.ρ_ae_KQ * f_K_tke) * ∇q_sno_en(wvec(prog_pr.q_sno))) - 
            #  (f_c_m * aux_tc_f.ρ_ae_KQ * f_k_tke_conv) * ∇q_sno_en_tke_conv(wvec(prog_pr.q_sno))
                # (ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_rai_en_tke_conv(wvec(prog_pr.q_rai)) # unstable [ needs higher K ]
                (Ifx(ρ_c * a_en * KMconv/(Prconv*Le) * ∇qc(wvec(Ifx(prog_pr.q_sno))))) # stable [ needs higher K ]
                # (ρ_f * Ifx(a_en * aux_en.tke_convective) * 1000 ) * ∇q_sno_en_tke_conv(wvec(prog_pr.q_sno))  # unstable
                # (Ifx(ρ_c * a_en * aux_en.tke_convective * 700  * ∇qc(wvec(Ifx(prog_pr.q_sno))))) # stable


                # (ρ_f * IfKQconv(a_en * KMconv/(Prconv*Le))) * ∇q_sno_en_tke_conv(wvec(prog_pr.q_sno)) 
               
        end

    else # ::: Normal :::

        if edmf.moisture_model isa NonEquilibriumMoisture
            c_KQl = mixing_length_params(edmf).c_KQl
            c_KQi = mixing_length_params(edmf).c_KQi
            @. aux_tc_f.diffusive_flux_ql = -(c_KQl * aux_tc_f.ρ_ae_KQ) * ∇q_liq_en(wvec(aux_en.q_liq))
            @. aux_tc_f.diffusive_flux_qi = -(c_KQi * aux_tc_f.ρ_ae_KQ) * ∇q_ice_en(wvec(aux_en.q_ice))

            # advective adjustment for condensible species
            liquid_advection_factor = mixing_length_params(edmf).c_KTKEql
            ice_advection_factor = mixing_length_params(edmf).c_KTKEqi
            @. aux_tc_f.diffusive_flux_ql += ρ_f * Ifx(a_en * sqrt(aux_en.tke) * liquid_advection_factor * aux_en.q_liq)
            @. aux_tc_f.diffusive_flux_qi += ρ_f * Ifx(a_en * sqrt(aux_en.tke) * ice_advection_factor * aux_en.q_ice)

        end
        c_KQr = mixing_length_params(edmf).c_KQr # maybe these would be shared but I think sedimentation might change things idk
        c_KQs = mixing_length_params(edmf).c_KQs # maybe these would be shared but I think sedimentation might change things idk
        @. aux_tc_f.diffusive_flux_qr = -(c_KQr * aux_tc_f.ρ_ae_KQ) * ∇q_rai_en(wvec(prog_pr.q_rai)) # add diffusive (SGS) flux for rain
        @. aux_tc_f.diffusive_flux_qs = -(c_KQs * aux_tc_f.ρ_ae_KQ) * ∇q_sno_en(wvec(prog_pr.q_sno)) # add diffusive (SGS) flux for snow

    end

    return nothing
end

function filter_small_moisture_vars(edmf::EDMFModel, state::State, param_set::APS)
    # convert small ql, qi to 0
    # this can induce small leaks in the mass/energy budget
    # if cond/evap sub/dep are a function of existing ql, qi, this could unnecessarily stunt what otherwise could be exponential growth depending on the base case
    # but the values we get sometimes of like 1e-100 are diabolical for explicit time stepping
    
    FT = float_type(state)
    q_min::FT = param_set.user_params.q_min

    if !iszero(q_min)
        aux_en = center_aux_environment(state)

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

        # Precip 
        prog_pr = center_prog_precipitation(state)
        @. prog_pr.q_rai = cutoff_small_values_positive(prog_pr.q_rai, q_min)
        @. prog_pr.q_sno = cutoff_small_values_positive(prog_pr.q_sno, q_min)

    end
    return nothing

end

function filter_small_vars(edmf::EDMFModel, state::State, param_set::APS)

    if !iszero(param_set.user_params.q_min)
        filter_small_moisture_vars(edmf, state, param_set)
    end

    # if !iszero(param_set.user_params.var_min)
        # terminal velocities [[ these might need to be don at point of creation, check order of operations... a reasonable q limit should help tho]]
    # end

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
    filter_gm_vars(edmf, state) # my addition, theyre not filtered anywhere... [[ moved to before aux_gm and prog are set in dycore.jl ]]
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

        if edmf.convective_tke_handler isa ConvectiveTKE  # convective tke sums with other tke in aux, so need to limit sum
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



        # if edmf.convective_tke_handler isa ConvectiveTKE # rn this one keeps them separate but we need to add tke's so it acts on the real qt etc [[ dont do this!!!]]
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

    filter_small_vars(edmf, state, param_set) # This is a little bit redundant w/ filter_precipitation_vars(), filter_updraft_vars()... should prolly combine them.

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
    to_scalar = is_tke ? toscalar : Base.identity
    ϕ_gm = getproperty(gm, ϕ_sym)
    ψ_gm = getproperty(gm, ψ_sym)
    ϕ_en = getproperty(aux_en, ϕ_sym)
    ψ_en = getproperty(aux_en, ψ_sym)
    area_en = aux_en_c.area

    Icd = is_tke ? CCO.InterpolateF2C() : Base.identity
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

    w_up_c = aux_tc.w_up_c

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

        @. w_up_c = Ic(aux_up_f[i].w) # this one gets updated every time it's used

        # thermo_params = TCP.thermodynamics_params(param_set)
        # Π = TD.exner.(thermo_params, aux_gm.ts)


        scalar_nudge_τᵣ = 20*60 # testing hardcoded to see if its' good before going all in [ would neet to pass forcing object to this fcn]

        if edmf.entrainment_type isa FractionalEntrModel
            @. tends_advec = -∇c(wvec(LBF(w_up_c * ρarea)))
            @. tends_ρarea = 
                tends_advec + 
                limit_tendency(edtl, ρarea * w_up_c * (entr_turb_dyn - detr_turb_dyn), ρarea + max(tends_advec,0), ρ_c * a_max - ρarea + max(-tends_advec,0), Δt)
        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. tends_advec = -∇c(wvec(LBF(w_up_c * ρarea)))
            @. tends_ρarea = 
                tends_advec + 
                limit_tendency(edtl, ρarea * (entr_rate_inv_s - detr_rate_inv_s), ρarea + max(tends_advec,0), ρ_c * a_max - ρarea + max(-tends_advec,0), Δt)
        end

        if edmf.entrainment_type isa FractionalEntrModel
            @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaθ_liq_ice)))
            @. tends_other = (ρ_c * θ_liq_ice_tendency_precip_formation) + (ρ_c * θ_liq_ice_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
            @. tends_ρaθ_liq_ice =
                tends_advec + 
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaθ_liq_ice / ρarea, prog_gm.ρθ_liq_ice / ρ_c), prog_gm.ρθ_liq_ice / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * w_up_c * entr_turb_dyn * θ_liq_ice_en) - (ρaθ_liq_ice * w_up_c * detr_turb_dyn) , ρaθ_liq_ice + max(tends_advec+tends_other,0), ρ_c * area_en * θ_liq_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                    (ρ_c * θ_liq_ice_tendency_precip_formation) # + 
                    # f_q*max(ρarea * aux_tc.θ_liq_ice_tendency_precip_sinks, 0) # + # qi_sub_dep contribution [[ i think this is needed so significant snow formation doesn't kill updrafts?]]
                    #+ ρarea*aux_up[i].dTdt/Π # test adding nudging in so that nudging doesnt prevent env and updraft from mixing and converging...

            @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaq_tot)))
            @. tends_other = (ρ_c * qt_tendency_precip_formation) + (ρ_c * qt_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
            @. tends_ρaq_tot =
                tends_advec +
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_tot / ρarea, prog_gm.ρq_tot / ρ_c), prog_gm.ρq_tot / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * w_up_c * entr_turb_dyn * q_tot_en) - (ρaq_tot * w_up_c * detr_turb_dyn), ρaq_tot + max(tends_advec+tends_other,0), ρ_c * area_en * q_tot_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                    (ρ_c * qt_tendency_precip_formation) 

        elseif edmf.entrainment_type isa TotalRateEntrModel
            @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaθ_liq_ice)))
            @. tends_other = (ρ_c * θ_liq_ice_tendency_precip_formation) + (ρ_c * θ_liq_ice_tendency_sedimentation)# external tendencies... we will limit entr/detr as the final tendency to
            @. tends_ρaθ_liq_ice =
                tends_advec +
                progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaθ_liq_ice / ρarea, prog_gm.ρθ_liq_ice / ρ_c), prog_gm.ρθ_liq_ice / ρ_c, tends_ρarea,
                    limit_tendency(edtl, (ρarea * entr_rate_inv_s * θ_liq_ice_en) - (ρarea * detr_rate_inv_s * θ_liq_ice_up), ρaθ_liq_ice+max(tends_advec+tends_other,0), ρ_c * area_en * θ_liq_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                (ρ_c * θ_liq_ice_tendency_precip_formation) +
                # f_q*max(ρarea * aux_tc.θ_liq_ice_tendency_precip_sinks, 0) # + # qi_sub_dep contribution [[ i think this is needed so significant snow formation doesn't kill updrafts?]]
                # ifelse(param_set.user_params.apply_nudging_to_updraft && (ρarea/ρ_c > 1.), ρarea * aux_up[i].dTdt/TD.exner(thermo_params, aux_gm.ts), FT(0)) # This version breaks because dTdt contains extraneous stuff and because we're not restoring directly to gm, we can have unresolved drift that simply gets compensated in env.
                ifelse(param_set.user_params.apply_nudging_to_updraft && (ρarea/ρ_c > FT(0.01)), ρarea * (aux_gm.H_nudge - aux_up[i].θ_liq_ice)/scalar_nudge_τᵣ, FT(0)) # dTdt includes hadv and stuff... ... nudging can kill the updraft though...

            @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaq_tot))) # this is ok ... we would use this in the grid mean/env as well but the offset tendencies lead to extreme transport, gradient sharpening etc.
            @. tends_other = (ρ_c * qt_tendency_precip_formation) + (ρ_c * qt_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
            @. tends_ρaq_tot =
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




            # TODO: Double check if it's better to use:: limit_tendency(edtl, (ρarea * w_up_c * entr_turb_dyn * q_liq_en) - (ρaq_liq * w_up_c * detr_turb_dyn) , max(ρaq_liq + tends_advec+tends_other,0), max(ρ_c * area_en * q_liq_en + -(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +

    

            if edmf.entrainment_type isa FractionalEntrModel
                
                @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaq_liq)))
                @. tends_other = (ρ_c * ql_tendency_precip_formation + ql_tendency_noneq) + (ρ_c * ql_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_liq =
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        limit_tendency(edtl, (ρarea * w_up_c * entr_turb_dyn * q_liq_en) - (ρaq_liq * w_up_c * detr_turb_dyn) , ρaq_liq + max(tends_advec+tends_other,0), ρ_c * area_en * q_liq_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (ql_tendency_precip_formation + ql_tendency_noneq))
                        
                @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaq_ice)))
                @. tends_other = (ρ_c * qi_tendency_precip_formation + qi_tendency_noneq) + (ρ_c * qi_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_ice =
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_ice / ρarea, prog_gm.q_ice), prog_gm.q_ice, tends_ρarea,
                        limit_tendency(edtl, (ρarea * w_up_c * entr_turb_dyn * q_ice_en) - (ρaq_ice * w_up_c * detr_turb_dyn) , ρaq_ice + max(tends_advec+tends_other,0), ρ_c * area_en * q_ice_en + max(-(tends_advec+tends_other),0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (qi_tendency_precip_formation + qi_tendency_noneq))
                    
            elseif edmf.entrainment_type isa TotalRateEntrModel
                @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaq_liq)))
                @. tends_other = (ρ_c * ql_tendency_precip_formation + ql_tendency_noneq) + (ρ_c * ql_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_liq =
                    tends_advec +
                    progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        limit_tendency(edtl, (ρarea * entr_rate_inv_s * q_liq_en) - (ρarea * detr_rate_inv_s * q_liq_up), ρaq_liq + max(tends_advec+tends_other, 0), ρ_c * area_en * q_liq_en + max(-(tends_advec+tends_other), 0), Δt), ρ_c, Δt, use_tendency_resolver, !use_tendency_resolver_on_full_tendencies) +
                        (ρ_c * (ql_tendency_precip_formation + ql_tendency_noneq))

                @. tends_advec = -∇c(wvec(LBF(w_up_c * ρaq_ice)))
                @. tends_other = (ρ_c * qi_tendency_precip_formation + qi_tendency_noneq) + (ρ_c * qi_tendency_sedimentation) # external tendencies... we will limit entr/detr as the final tendency to not violate this contract...
                @. tends_ρaq_ice =
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
                    (ρarea * w_up_c * entr_turb_dyn * θ_liq_ice_en) - (ρaθ_liq_ice * w_up_c * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaθ_liq_ice)

                @. tends_ρaq_tot += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_tot / ρarea, prog_gm.ρq_tot / ρ_c), prog_gm.ρq_tot / ρ_c, tends_ρarea,
                    (ρarea * w_up_c * entr_turb_dyn * q_tot_en) - (ρaq_tot * w_up_c * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_tot)

                if edmf.moisture_model isa NonEquilibriumMoisture
                    @. tends_ρaq_liq += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_liq / ρarea, prog_gm.q_liq), prog_gm.q_liq, tends_ρarea,
                        (ρarea * w_up_c * entr_turb_dyn * q_liq_en) - (ρaq_liq * w_up_c * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_liq)

                    @. tends_ρaq_ice += progvar_area_tendency_resolver(edtl, ρarea, ifelse(ρarea > 0, prog_up[i].ρaq_ice / ρarea, prog_gm.q_ice), prog_gm.q_ice, tends_ρarea,
                        (ρarea * w_up_c * entr_turb_dyn * q_ice_en) - (ρaq_ice * w_up_c * detr_turb_dyn), ρ_c, Δt, use_tendency_resolver, use_tendency_resolver_on_full_tendencies, tends_ρaq_ice)
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

        buoyf = aux_tc_f.temporary_f1
        @. buoyf = I0f(buoy)
        a_upf = aux_tc_f.temporary_f2
        @. a_upf = I0f(a_up)
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

        # === GRAFT onto updraft === #
        if (edmf.convective_tke_handler isa AbstractYesConvectiveTKEHandler) && !iszero(local tke_conv_entr_detr_rate_inv_s::FT = edmf.convective_tke_handler.entr_detr_rate_inv_s) # local helps make the assignment somehow
            prog_en = center_prog_environment(state)    

            #= 
                If we're generating downdrafts, they operate like updrafts probably (smaller area and all that) so we don't need to broaden updraft.
                However, since they operate in the environment, we don't need to spark any updraft area,
                    the compensating upwards return flow should be broad and slow, and primarily in the env, so we'll put it in normal TKE.
                In a future model with downdrafts, you might consider sparking coherent downdraft areas.

                Since in principle condensate forms in the updrafts, each updraft should in principle risk termination by downdraft formation
                rather than getting to gently fade out as it loses buoyancy and maybe getting a boost from the subsiding env.
                Then again you could have say surface convection rising into air that is dry with falling snow from a layer above or something.
            =#

            a_en_tke_conv = aux_tc.temporary_1
            @. a_en_tke_conv = FT(0.75) - aux_en.frac_supersat/2 # high supersat, smaller, stronger updrafts. bound to 0.25 to 0.75 of max area
            Δa_en_tke_conv = a_en_tke_conv # alias
            tends_ρarea_tke_conv = aux_tc.temporary_2
            zero_field!(tends_ρarea_tke_conv)
            w_up_c = aux_tc.w_up_c
            @inbounds for k in real_center_indices(grid)
                if (aux_en.tke_convective_production[k] > eps(FT)) #&& (aux_en.frac_supersat[k] > 0.5)
                    Δa_en_tke_conv[k] = min(max(a_en_tke_conv[k] - aux_up[i].area[k], FT(0)), max(edmf.max_area - aux_up_i.area[k], FT(0)))
                    area_growth_rate = aux_en.tke_convective_production[k] / max(sqrt(w_up_c[k]^2), FT(0.1)^2)  # m²/s growth rate to match TKE production. If w is small, assume small w of 0.1 m/s to avoid huge growth rates
                    tends_ρarea_tke_conv[k] = -limit_tendency(edtl, -ρ_c[k] * Δa_en_tke_conv[k] * tke_conv_entr_detr_rate_inv_s * area_growth_rate, Δa_en_tke_conv[k] * ρ_c[k], Δt) # TKE changing 0.5 in 500m is reasonable so x 1000
                end
            end


            @inbounds for k in real_center_indices(grid)

                ρa_up  = prog_up[i].ρarea[k]
                inv_ρa = inv(max(ρa_up, eps(FT)))
                q_up   = prog_up[i].ρaq_tot[k] * inv_ρa
                w2 = max(w_up_c[k]^2, FT(0.1)^2)
                λ  = tke_conv_entr_detr_rate_inv_s * (aux_en.tke_convective_production[k] / w2)
                ṁ_mix = ρ_c[k] * aux_up[i].area[k] * λ

                # Entrain moist air
                q_up   = prog_up[i].ρaq_tot[k] * inv_ρa
                q_eff = aux_en.q_tot[k] + sqrt(eps(FT)) + sqrt(aux_en.QTvar[k]) * sqrt(FT(2) / FT(π))
                # tends_ρaq_tot[k] += tends_ρarea_tke_conv[k] * q_eff + ṁ_mix * (q_eff - q_up)
                tends_ρaq_tot[k] += -limit_tendency(edtl, -((tends_ρarea_tke_conv[k] * q_eff) + (ṁ_mix * (q_eff - q_up))), prog_up[i].ρaq_tot[k], Δt)

                # Old existing air is likely to have -HVar due to being displaced vertically and θli increasing with height. However, new air being entrained is likely to be +HVar buoyant air.
                θ_up   = prog_up[i].ρaθ_liq_ice[k] * inv_ρa 
                θ_eff = aux_en.θ_liq_ice[k] - sqrt(eps(FT)) + sqrt(aux_en.Hvar[k]) * sqrt(FT(2) / FT(π)); # old air is more likely to be displaced, so it is negative like the typical thetali gradient would indicate from perturbations
                # tends_ρaθ_liq_ice[k] += tends_ρarea_tke_conv[k] * θ_eff + ṁ_mix * (θ_eff - θ_up)
                tends_ρaθ_liq_ice[k] += -limit_tendency(edtl, -((tends_ρarea_tke_conv[k] * θ_eff) + (ṁ_mix * (θ_eff - θ_up))), prog_up[i].ρaθ_liq_ice[k], Δt)


                if tends_ρarea_tke_conv[k] > zero(FT) # Generate updraft area
                    # Give to updraft as we take from environment 
                    tends_ρarea[k] += tends_ρarea_tke_conv[k]
                    if edmf.moisture_model isa NonEquilibriumMoisture
                        tends_ρaq_liq[k] += -limit_tendency(edtl, -tends_ρarea_tke_conv[k] * aux_en.q_liq[k], ρ_c[k] * area_en[k] * aux_en.q_liq[k], Δt)
                        tends_ρaq_ice[k] += -limit_tendency(edtl, -tends_ρarea_tke_conv[k] * aux_en.q_ice[k], ρ_c[k] * area_en[k] * aux_en.q_ice[k], Δt)
                    end
                    w_convective = sqrt(2 * aux_en.tke_convective[k])
                    tendencies_en.ρatke_convective[k] -= tends_ρarea_tke_conv[k] * w_convective^2 # remove tke we are sending to updraft
                end
                # if tends_ρarea_tke_conv[k] < zero(FT)  # Destroy updraft area
                #     # do nothing right now, we only graft for growth
                # end
            
            end

            #=
                We tried before to boost w directly but found it very difficult.
                The main problem relates to how to distribute w from the environment to the updraft since it's the same parcels. Doing it gradually can lead to stunted transport, incorrect anvil placements etc.

                You tend to end up with half the energy in the env and half in the updraft, both too weak to efficiently transport mass
                We've switched to just sparking area instead, which seems to work better overall, so for now we can leave this out
            =#

            # tends_ρarea_tke_conv_f = aux_tc_f.temporary_f1
            # @. tends_ρarea_tke_conv_f = I0f(tends_ρarea_tke_conv) # face version
            # tke_convective_f = aux_tc_f.temporary_f2
            # @. tke_convective_f = I0f(aux_en.tke_convective)

            # tke_convective_production_f = aux_tc_f.temporary_f3
            # @. tke_convective_production_f = I0f(aux_en.tke_convective_production)
            # a_en_f = aux_tc_f.temporary_f4
            # @. a_en_f = I0f(area_en)
            # a_upf = I0f.(aux_up[i].area) # need to figure out a sol'n here
            # @inbounds for k in real_face_indices(grid)
            #     if tends_ρarea_tke_conv_f[k] > FT(0)  # Generate updraft area
            #         # tends_ρaw[k] += max(-limit_tendency(edtl, -tends_ρarea_tke_conv_f[k] * w_convective, ρ_f[k] * a_en_f[k]/2 * w_convective, Δt), FT(0)) # max because w/ the upppper bound on area, we can have area tendency 0 inside the latent heating areas, and then the negatives outside leak in.
            #         # tends_ρaw[k] += sqrt(tke_convective_production_f[k] / (ρ_f[k] * a_en_f[k] + eps(FT)) ) * ρ_f[k] * (a_upf[k] + tke_convective_production_f[k] * Δt)  # small boost to get it going

            #         # since we have the area limiter, we want add w to existing area as well [[ probably should have also done to the other tracers.. idk ]]
            #         # δ_ρaw = max(ρ_f[k] * a_upf[k] * w_convective - prog_up_f[i].ρaw[k], FT(0))
            #         # tends_ρaw[k] += -limit_tendency(edtl, -δ_ρaw * tke_conv_entr_detr_rate_inv_s, δ_ρaw, Δt)
            #         # tends_ρaw[k] += δ_ρaw/(Δt) # add the rest directly... should be small

            #     end
            #     if tends_ρarea_tke_conv_f[k] < FT(0)  # Destroy updraft area
            #         # we want to destroy w at a proportional rate to the rate of loss in env... that is tke_dissipation / tke_conv should match our loss of w / w_up
            #         # δ_ρaw_equiv = prog_up_f[i].ρaw[k] * -( (I0f.(aux_en.tke_convective_dissipation))[k] / max((I0f.(prog_en.ρatke_convective))[k], FT(1e-10)) ) # how much ρaw equivalent to lose based on tke dissipation rate
            #         # tends_ρaw[k] += min(limit_tendency(edtl, tends_ρarea_tke_conv_f[k] * w_convective + δ_ρaw_equiv, prog_up_f[i].ρaw[k] + max(tends_ρaw[k], FT(0))*Δt, Δt), FT(0)) # min because leakage
            #         # tends_ρaw[k] += sqrt(tke_convective_production_f[k] / (ρ_f[k] * a_en_f[k] + eps(FT)) ) * ρ_f[k] * a_upf[k]  # small boost to get it going
            #     end
            # end

            # aux_en.tke_convective_production[k] = FT(0)
            # aux_en.tke_convective_dissipation[k] = FT(0)
            # aux_en.tke_convective_advection[k] = FT(0)
            # tendencies_en.ρatke_convective[k] = FT(0)
            # prog_en.ρatke_convective[k] = FT(0)
            # aux_en.tke_convective[k] = FT(0)
            zero_field!(aux_en.tke_convective)
            zero_field!(tendencies_en.ρatke_convective)
            zero_field!(prog_en.ρatke_convective)
            zero_field!(aux_en.tke_convective_production)
            zero_field!(aux_en.tke_convective_dissipation)
            zero_field!(aux_en.tke_convective_advection)

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
    is_tke = (covar_sym === :tke)
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
        ∇c_shear = CCO.GradientF2C() # This one is Gradient not Divergence so it's not the global one
        @. shear = 2 * ρ_c * area_en * k_eddy * LA.dot(∇c_shear(If(ϕ_en)), k̂) * LA.dot(∇c_shear(If(ψ_en)), k̂)

        # == Get convective TKE contribution to covariance == #
        # if (edmf.convective_tke_handler isa ConvectiveTKE) && (covar_sym ∈ (:QTvar, :Hvar, :HQTcov))
        if (edmf.convective_tke_handler isa ConvectiveTKE) && !(is_tke) # [ shear production not relevant for convective tke ]
            ℓ_mix_conv = FT(500) 
            tke_conv = aux_en.tke_convective
            # Removed local stability damping (Ri) to prevent grid-scale stripes. 
            # Convective mixing scale is driven by TKE_conv, not local db/dz.
            @. shear += 2 * ρ_c * area_en * (ℓ_mix_conv * sqrt(2 * max(tke_conv, 0))) * LA.dot(∇c_shear(If(ϕ_en)), k̂) * LA.dot(∇c_shear(If(ψ_en)), k̂) * FT(0.1)
        end
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
    Ic = is_tke ? CCO.InterpolateF2C() : Base.identity

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
    to_scalar = is_tke ? toscalar : Base.identity
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
    Idc = is_tke ? Ic : Base.identity
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

"""
This seems to be diagnostic only so it doesnt really do anything... (there are no uses of .dissipation anywhere)
"""
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

    if (edmf.convective_tke_handler isa ConvectiveTKE) && (covar_sym === :tke)
        @. dissipation = ρ_c * area_en * covar * max(tke_en, 0)^FT(0.5) / max(mixing_length, FT(1.0e-3)) * c_d
    else
        @. dissipation = ρ_c * area_en * covar * max(tke_en, 0)^FT(0.5) / max(mixing_length, FT(1.0e-3)) * c_d
    end

    is_tke = (covar_sym === :tke)

    # == Add convective TKE contribution to dissipation ==#
    if (edmf.convective_tke_handler isa ConvectiveTKE) && !is_tke 
        FT = float_type(state)
        ℓ_mix_base = FT(500) # Matches base length in shear
        
        # dbdz = aux_tc.∂b∂z # should we be using this or the moist one?
        # dbdz = instability
        
        # A. Calculate Velocity Scale
        # w = sqrt(2 * TKE)
        w_conv = aux_tc.temporary_1 # Reusing temp array
        @. w_conv = sqrt(2 * max(aux_en.tke_convective, 0))

        # B. Calculate Effective Length (ℓ_eff)
        # MUST MATCH SHEAR LOGIC EXACTLY:
        # If stability (N^2) is high relative to w^2, the eddies shrink.
        # This ensures the dissipation speeds up (rate increases) in inversions.
        ℓ_eff = aux_tc.temporary_2 # Reusing temp array
        # @. ℓ_eff = ℓ_mix_base / (1 + ( (sqrt(max(dbdz, 0)) * ℓ_mix_base) / (w_conv + eps(FT)) )^2)
        @. ℓ_eff = FT(ℓ_mix_base)

        # C. Add to Dissipation Rate
        # Rate = w / ℓ_eff
        @. dissipation += ρ_c * area_en * covar * w_conv / max(ℓ_eff, FT(1.0e-3)) * c_d
    end

    return nothing
end


# """
# We have no way for buoyancy to generate TKE outside of the updraft-env split... ideally we'd have some way
# The buoy variable should already exist so we just have to use it...
# """
# function compute_covariance_buoyancy(
#     edmf::EDMFModel,
#     state::State,
#     ::Val{covar_sym},
# ) where {covar_sym}
#     FT = float_type(state)
#     aux_en_2m = center_aux_environment_2m(state)
#     aux_covar = getproperty(aux_en_2m, covar_sym)
#     covar = getproperty(center_aux_environment(state), covar_sym)
#     buoy = aux_covar.buoy

#     is_tke = (covar_sym === :tke)

#     if (edmf.convective_tke_handler isa ConvectiveTKE) && !is_tke
#         error("Buoyancy contribution to covariance is not implemented yet. Right now the diagnostic just uses shear production")

#     end
       
#     return nothing
# end

function compute_en_tendencies!(
    edmf::EDMFModel,
    state::State,
    # param_set::APS,
    param_set::TCP.TurbulenceConvectionParameters{FT},
    surf::SurfaceBase,
    ::Val{covar_sym},
    ::Val{prog_sym},
    Δt::FT,
    cfl_limit::FT,
    use_fallback_tendency_limiters::Bool,
) where {covar_sym, prog_sym, FT <: Real}
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
    # FT = float_type(state)
    thermo_params = TCP.thermodynamics_params(param_set)

    Ifw = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.SetValue(FT(0)))

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
        # prog_bcs = (; bottom = CCO.Extrapolate(), top = CCO.SetGradient(wvec(FT(0)))) # Apparently Extrapolate is an error
        prog_bcs = (; bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))
    end

    If = CCO.InterpolateC2F(; aeK_bcs...)
    ∇f = CCO.GradientC2F(; prog_bcs...)

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
        
        if (edmf.convective_tke_handler isa ConvectiveTKE) 
            # pass (we don't have anything here yet... we dont use prognostic tke but we'd neeed to figure out the rate of tke generation by convective tke using the tracer gradient and assuming maybe the same half up half down split)
            error("Prognostic covariance tendencies for convective TKE not implemented yet.")
        end

    end

    mixing_length = aux_tc.mixing_length
    min_area = edmf.minimum_area

    area_en = aux_en.area
    tke_en = aux_en.tke

    # parent(D_env) .= 0
    zero_field!(D_env)

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


    if is_tke && (edmf.convective_tke_handler isa AbstractYesConvectiveTKEHandler)

        a_en = aux_en.area
        # We want here to keep our `convective TKE` and move it around etc, but dont dissipate it unless there's strong stability.
        ρatke_convective = prog_en.ρatke_convective
        tend_ρatke_convective = tendencies_en.ρatke_convective
        ∂MSE∂z = aux_en.∂MSE∂z

        # -=- W convective =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        w_convective = aux_tc.temporary_3 # 1 and 2 already taken but avaiable again
        @. w_convective = sqrt(2 * ρatke_convective/(ρ_c * a_en)) # tke = 1/2 ρ w^2  => w = sqrt(2*tke) , tke is per unit mass here
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #


        # -- instability calculation -- #
        #= unstable if:
            - subsaturated and dry static instability [[ dry static energy gradient looks like ]]
            - saturated and moist static instability
            - saturated wrt ice if below freezing, liquid if above

            We track stability and latent heating, and use them to calculate instability.
            Stability is tracked separately however for the purposes of convective TKE dissipation.
        =#
        instability = aux_en.instability
        @. instability = zero(FT)
        stability = aux_en.stability
        @. stability = zero(FT)

        g = TCP.grav(param_set)
        c_p = TCP.cp_d(param_set)
        R_d = TCP.R_d(param_set)
        f_crit = FT(0.50)
        ∂MSE∂z = aux_en.∂MSE∂z
        @inbounds for k in real_center_indices(grid)
            T = aux_en.T[k]
            ∂MSE∂z_here = ∂MSE∂z[k] # biased downwards to avoid spurious instability at sharp gradients

            frac_supersat = aux_en.frac_supersat[k]
            frac_subsat = one(FT) - frac_supersat
            ts_here = aux_en.ts[k]

            # This ramps up efficiency from 0 at 50% subsat to 1 at 100% sat
            peak_enhancement = FT(1); #
            sat_efficiency = peak_enhancement * clamp((frac_supersat - f_crit) / (one(FT) - f_crit), zero(FT),  one(FT)) # I think this double counts the effect of saturation and overly kills TKE generation..
            # sat_efficiency = one(FT) # for testing [[ i think this one isn't overly onerous for reducing tke production]]

            Π_here = TD.exner_given_pressure(thermo_params, aux_en.p[k], TD.PhasePartition(thermo_params, ts_here))
            ρ_here = TD.air_density(thermo_params, ts_here)
            ∂b∂θv = g * (R_d * ρ_here / aux_en.p[k]) * Π_here
            ∂b∂z_raw = ∂b∂θv * aux_tc.∂θv∂z[k] # we could use bg_model but it has its own sat/unsat breakdown [[ b_grad as saved contains the moist part already ]]

            latent_heating = aux_en.latent_heating_pos[k] + -aux_en.latent_heating_neg[k] # to be precise you'd want a split into something in the subsat and something in the supersat part, but this is close enough for now...
            latent_heating = sign(latent_heating) * min(abs(latent_heating), FT(3)) # avoid excessive latent heating (i.e. oscillations)

            # supersat part
            if !iszero(frac_supersat)
                # Moist static energy gradient
                instability[k] += (-max(∂MSE∂z_here, -1) * (g/(c_p * aux_en.T[k]))) * frac_supersat * abs(latent_heating) # * sat_efficiency # partially saturated
                stability[k] += FT(0)
            end

            # subsat part
            if !iszero(frac_subsat)

                # Dry evap contribution...
                if (aux_en.q_liq[k] + aux_en.q_ice[k] + prog_pr.q_rai[k] + prog_pr.q_sno[k]) > FT(1e-10) # Have condensate that can evaporate. We are not unstable to these since they can't spark vigorous updrafts and grow the cycle.
                    # We are unstable to evaporation fueled downdrafts... We don't let these generate TKE but they can lessen dissipation[[ technically liquid should release more energy but... ]]
                    # instability[k] += frac_subsat * -∂b∂z_raw *  (((latent_heating < 0) && (∂b∂z_raw < 0)) ? latent_heating : zero(FT)) * sat_efficiency
                    # instability[k] += frac_subsat * min(-∂b∂z_raw, zero(FT)) # error("not implemented yet") I cant figure how to properly reduce the stability without necessarily going to ofar.., esp since latent heating can exceed 1.
                    # instability[k] += zero(FT)
                    instability[k] += (-max(∂MSE∂z_here, -1) * (g/(c_p * aux_en.T[k]))) * frac_subsat * -abs(latent_heating) * sat_efficiency  # partially saturated

                    stability[k] += frac_subsat * max(∂b∂z_raw, zero(FT))
                else # We have no condensate and are subsat. So it's just dry stability. We don't generate TKE based on dry static instability alone since we lack the latent heat enhancement engine and regular tke should handle this
                    # instability[k] += frac_subsat * min(-∂b∂z_raw, zero(FT)) # dry static instability
                    instability[k] += zero(FT)
                    stability[k] += frac_subsat * max(∂b∂z_raw, zero(FT))
                end
            end
        end

       
        # =-=-=- Production calculation -=-=-= #
        ρatke_convective_production = aux_en.tke_convective_production # for some reason these have ρa in ρatke but are named this way so we follow
        K_buoy::FT = edmf.convective_tke_handler.buoyancy_coeff

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

        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        # @. ρatke_convective_production = (ρ_c * a_en) * K_buoy * max(instability, 0) * sqrt(2 * max(ρatke_convective / (ρ_c * a_en + eps(FT)), eps(FT)))
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #



        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        #   UNIFIED BUOYANCY PRODUCTION: The "Spark" (Mean) + The "Engine" (SGS)
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

        # 1. Constants
        g       = TCP.grav(param_set)
        c_p     = TCP.cp_d(param_set)
        R_d     = TCP.R_d(param_set)
        ℓ_mix   = FT(500) 


        @inbounds for k in real_center_indices(grid)
            
            # --- A. Thermodynamics ---
            T_here  = aux_en.T[k]
            ts_here = aux_en.ts[k]
            p_here  = aux_en.p[k]
            ρ_here  = TD.air_density(thermo_params, ts_here)
            
            # --- B. Velocity (With Floor) ---
            ρ_c_here  = ρ_here 
            a_en_here = aux_en.area[k]
            tke_val   = ρatke_convective[k]
            
            if (ρ_c_here * a_en_here) > FT(1e-6)
                tke_specific = tke_val / (ρ_c_here * a_en_here)
            else
                tke_specific = FT(0)
            end
            
            # Floor w_conv so SGS term works at initialization
            w_conv_here = max(w_convective[k], FT(0.01))


            # --- C. Timescales ---
            # We calculate this here so we can use it for the Unit Fix in Term Mean
            tau_conv = ℓ_mix / (w_conv_here + FT(0.1)) 

            inv_tau_micro = FT(0)
            f_active = aux_en.frac_supersat[k]
            
            # Liquid
            if (f_active > 0) || (aux_en.q_liq[k] > FT(1e-9))
                τ_liq_here = (edmf.moisture_model isa NonEquilibriumMoisture) ? aux_en.τ_liq[k] : FT(0)
                inv_tau_micro += 1.0 / (τ_liq_here + eps(FT))
            end
            # Ice
            if (T_here < FT(273.15)) 
                if (f_active > 0) || (aux_en.q_ice[k] > FT(1e-9))
                    τ_ice_here = (edmf.moisture_model isa NonEquilibriumMoisture) ? aux_en.τ_ice[k] : FT(0)
                    inv_tau_micro += 1.0 / (τ_ice_here + eps(FT))
                end
            end
            # Precip
            if prog_pr.q_rai[k] > FT(1e-9); inv_tau_micro += 1.0 / (FT(10) + eps(FT)); end
            if prog_pr.q_sno[k] > FT(1e-9); inv_tau_micro += 1.0 / (FT(1e5) + eps(FT)); end

            # Efficiency (Phi)
            if inv_tau_micro > FT(1e-12)
                tau_micro = 1.0 / inv_tau_micro
                phi = tau_conv / (tau_conv + tau_micro)
            else
                phi = FT(0)
            end

            # --- D. Instability Calculation ---
            # 1. Gradients
            Π_here = TD.exner_given_pressure(thermo_params, p_here, TD.PhasePartition(thermo_params, ts_here))
            L_s = TD.latent_heat_sublim(thermo_params, T_here)
            ∂b∂θv  = g * (R_d * ρ_here / p_here) * Π_here
            N2_dry = ∂b∂θv * aux_tc.∂θv∂z[k]
            
            ∂MSE∂z_here = aux_en.∂MSE∂z[k]

            # adjustment for marginal stable regions
            ∂MSE∂z_here -= min(
            hypot(
                c_p * Π_here * sqrt(max(aux_en.Hvar[k],  FT(0))),
                L_s              * sqrt(max(aux_en.QTvar[k], FT(0))),
            ) / max(aux_tc.mixing_length[k], FT(1)),
            FT(1e-3),
        )
           
            N2_moist    = (g / (c_p * T_here)) * ∂MSE∂z_here

            # 2. Mean Latent Heating Term [UNIT FIX]
            latent_heating = aux_en.latent_heating_pos[k] + -aux_en.latent_heating_neg[k]
            latent_heating = sign(latent_heating) * min(abs(latent_heating), FT(3)) 
            
            # Correction: Multiply by tau_conv (l/w) to convert Rate (Jerk) -> Buoyancy (Acceleration)
            term_mean = (g / (c_p * T_here)) * abs(latent_heating) * tau_conv * phi

            # 3. SGS Gradient Term [LOGIC FIX]
            f_effective = f_active * phi
            
            # Correction: Hard switch. Do not let dry stability penalize moist instability.
            if f_effective > FT(1e-4)
                # Active/Cloudy: Use MOIST gradient only.
                term_sgs = -ℓ_mix * N2_moist
            else
                # Inactive/Dry: Use DRY gradient.
                term_sgs = -ℓ_mix * N2_dry
            end

            # 4. Combine
            # Now both terms are accelerations (m/s^2) and can be added.
            total_instability = term_mean + term_sgs

            # --- E. Production & Output ---
            
            K_buoy = edmf.convective_tke_handler.buoyancy_coeff
            
            if total_instability > 0
                # Production = Flux (w * instability) * Efficiency
                production_rate = w_conv_here * total_instability * K_buoy
                
                # Safety cap
                production_rate = min(production_rate, FT(0.5))

                ρatke_convective_production[k] = (ρ_c_here * a_en_here) * production_rate
                aux_en.instability[k] = total_instability
                aux_en.stability[k]   = 0
            else
                ρatke_convective_production[k] = 0
                aux_en.instability[k] = 0
                aux_en.stability[k]   = -total_instability
            end
            
        end
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    end


    if is_tke && (edmf.convective_tke_handler isa ConvectiveTKE)

        # confine tke_max to CFL
        tke_max = aux_tc.temporary_1 # 1 and 2 used above but available again
        Δz = get_Δz(prog_gm.ρ) # This is just a regular array not a CC field so we 
        @inbounds for k in real_center_indices(grid) # because get_Δz returns a regular array not a CC Field, we do this for loop
            if tke_max[k] > 0
                w_max = cfl_limit * Δz[k] / Δt
                tke_max[k] = min(tke_max[k], FT(0.5) * ρ_c[k] * w_max^2) # kinetic energy from max w
            end
        end

        if any(!isfinite, ρatke_convective_production)
            @warn "Non-finite values in convective TKE production computation."
            error("K_buoy = $(K_buoy); a_en = $(a_en); instability = $(parent(instability)); tke = $(parent(aux_en.tke))")
        end
        
        # -- advection/diffusion -- #
        # ∇c_tkec = CCO.DivergenceF2C(; bottom = CCO.SetValue(wvec(FT(0))), top = CCO.SetValue(wvec(FT(0)))) # we have no sure surface value so we could just extrapolate. However since our sfc generates no tke and we leave that to updraft and shear (for now, could revisit if it's problematic), we'll just go w/ 0 on both (technically tke could go out the top but we enforce w=0 so...)
        # advection (is diffusive, like pressure diffusion I guess...) # [[ they use ::         @. nh_press_adv = Int(ᶠinterp_a(a_up) > 0) * ρ_f * ᶠinterp_a(a_up) * α_a * w_up * ∇(wvec(Ifc(w_up)))      ]]
        # We'll have a fast advection part and a slow diffusion part. diffusion is up gradient, advection is with w_convective [[ could in principle go both directions but should be biased up for convection so go up ]]
        ρatke_convective_advection = aux_en.tke_convective_advection # for some reason these have ρa in ρatke but are named this way so we follow
        k_adv::FT = edmf.convective_tke_handler.advection_coeff
        # @. tke_convective_advection = ρ_c * a_en * k_adv * w_convective * ∇c(wvec(Ifw(tke_convective)))


        if edmf.convective_tke_handler.transport_tke_by_advection

            #=
                ======================================================================
                DERIVATION: SKEWED TWO-STREAM TKE ADVECTION
                ======================================================================

                Objective:
                Compute the vertical advection of Convective TKE (e) by resolving the
                asymmetry between updrafts (up) and downdrafts (dn).

                Standard TKE advection assumes a single bulk velocity. However, in moist
                convection, updrafts and downdrafts often have distinct areas and intensities.
                (e.g., Narrow/Strong Updrafts vs. Broad/Weak Subsidence). This skewness
                results in net TKE transport even if the net mass flux is zero.

                1. Governing Constraints
                ------------------------
                (a) Mass Balance: Net mass flux must be zero.
                    ρ * a_up * w_up + ρ * a_dn * w_dn = 0
                    
                    Defining γ (gamma) = a_up / a_dn:
                    => w_dn = -γ * w_up

                (b) Energy Conservation: The sum of kinetic energy in both streams must 
                    equal the grid-mean TKE (e_mean).
                    e_mean = 0.5 * [ a_up * w_up^2 + a_dn * w_dn^2 ]

                2. Geometric Closure (The Physics)
                ------------------------
                We determine the Updraft Area Fraction (a_up) based on the thermodynamic
                driving forces available in the grid cell.
                
                Let P be the "Potential Forcing Power" (Buoyancy * Length):
                - P_up ~ Latent Heat * saturation_fraction + Base_TKE
                - P_dn ~ Evap Cooling * subsaturation_fraction + Base_TKE
                
                Scaling Law: Stronger forcing creates narrower, faster streams.
                We partition area inversely proportional to the square root of potential:
                    a_up = sqrt(P_dn) / (sqrt(P_up) + sqrt(P_dn))
                
                Behaviors:
                - Moist Unstable (P_up >> P_dn): a_up -> 0 (Narrow Updrafts)
                - Dry Precipitating (P_dn >> P_up): a_up -> 1 (Narrow Downdrafts)
                - Neutral (P_up ≈ P_dn): a_up -> 0.5 (Symmetric Turbulence)

                3. Solution for Velocities
                ------------------------
                Substitute w_dn = -γ * w_up into the Energy Conservation equation:
                
                e_mean = 0.5 * w_up^2 * [ a_up + a_dn * γ^2 ]
                
                Solving for w_up:
                w_up = sqrt( 2 * e_mean / (a_up + a_dn * γ^2) )

                4. Final Flux Calculation
                ------------------------
                Flux = Flux_up + Flux_dn
                    = (MassFlux_up * KineticEnergy_up) + (MassFlux_dn * KineticEnergy_dn)
                    = (ρ * a_up * w_up) * (0.5 * w_up^2) + (ρ * a_dn * w_dn) * (0.5 * w_dn^2)
                    = 0.5 * ρ * [ a_up * w_up^3 + a_dn * w_dn^3 ]

                Since w_dn is negative, w_dn^3 is negative, representing downward transport.
                ======================================================================
            =#

            # --- Physical Constants & Parameters ---
            g      = TCP.grav(param_set)
            c_p    = TCP.cp_d(param_set)
            L_v    = TCP.LH_v0(param_set)
            ℓ_mix  = FT(500)       # Mixing length for potential scaling
            ε_area = FT(0.01)      # Hard geometric limit (1% area)
            ε_tiny = FT(1e-10)     # Division safety floor
            
            has_micro = (edmf.moisture_model isa NonEquilibriumMoisture)

            # --- Aliases (Reused memory) ---
            w_up = w_convective # temporary_3
            w_dn = aux_tc.temporary_2 
            a_up = aux_tc.temporary_4 

            # ------------------------------------------------------------------------
            # 1) Partition TKE into Skewed Up/Down Streams
            # ------------------------------------------------------------------------
            @inbounds for k in real_center_indices(grid)
                
                # --- State Variables ---
                # Ensure TKE is non-negative and convert to specific energy (m^2/s^2)
                tke_mean = max(ρatke_convective[k], FT(0)) / ρ_c[k]
                
                T        = aux_en.T[k]
                f_sat    = aux_en.frac_supersat[k] 
                q_L_tot  = aux_en.q_liq[k] + prog_pr.q_rai[k]
                q_I_tot  = aux_en.q_ice[k] + prog_pr.q_sno[k]

                # --- A. Determine Geometry (a_up) ---
                # Thermodynamic coefficient: g * Lv / (Cp * T)
                therm_coeff = g * L_v / (c_p * T)
                
                # Relaxation weighting
                # Use sqrt(2*TKE) as a characteristic velocity for timescale estimation
                w_scale  = sqrt(2 * tke_mean + ε_tiny) 
                tau_conv = ℓ_mix / w_scale
                
                # Microphysics relaxation factors (phi)
                tau_L    = (has_micro && q_L_tot > FT(1e-9)) ? aux_en.τ_liq[k] : FT(0)
                tau_I    = (has_micro && q_I_tot > FT(1e-9)) ? aux_en.τ_ice[k] : FT(0)
                phi_L    = tau_conv / (tau_conv + tau_L + FT(1e-12))
                phi_I    = tau_conv / (tau_conv + tau_I + FT(1e-12))

                # Potential Magnitudes (P)
                # Pot_micro represents the maximum potential buoyant work available from phase change.
                # We add Base_buoy (1.0) to ensure the system defaults to symmetry (a_up=0.5)
                # in the absence of microphysics.
                Pot_micro = therm_coeff * (q_L_tot * phi_L + q_I_tot * phi_I) * ℓ_mix
                Base_buoy = FT(1) 

                # P_up:  Boosted by Latent Heating in saturated fraction
                # P_dn:  Boosted by Evap Cooling in subsaturated fraction
                P_up = (Pot_micro * f_sat) + Base_buoy
                P_dn = (Pot_micro * (FT(1) - f_sat)) + Base_buoy

                # Geometry Logic: Stronger potential -> Narrower stream
                # a_up = sqrt(P_dn) / (sqrt(P_up) + sqrt(P_dn))
                sqrt_P_up = sqrt(P_up)
                sqrt_P_dn = sqrt(P_dn)
                
                val_a_up = sqrt_P_dn / (sqrt_P_up + sqrt_P_dn)
                val_a_up = clamp(val_a_up, ε_area, FT(1) - ε_area)
                val_a_dn = FT(1) - val_a_up
                
                a_up[k] = val_a_up

                # --- B. Velocity Shape from Mass Balance ---
                # w_dn = - gamma * w_up
                gamma = val_a_up / val_a_dn

                # --- C. Velocity Magnitude from TKE Conservation ---
                # Solves: TKE = 0.5 * w_up^2 * (a_up + a_dn * gamma^2)
                shape_factor = val_a_up + val_a_dn * gamma^2
                val_w_up     = sqrt(2 * tke_mean / shape_factor)
                
                # Assign Velocities (w_up > 0, w_dn < 0)
                w_up[k] = val_w_up
                w_dn[k] = -val_w_up * gamma
            end

            # ------------------------------------------------------------------------
            # 2) Gradient Operators
            # ------------------------------------------------------------------------
            ∇f_tke = CCO.GradientC2F(
                bottom = CCO.SetGradient(wvec(FT(0))),
                top    = CCO.SetGradient(wvec(FT(0)))
            )
            LBF = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(FT(0))) # Upwind for w > 0
            RBF = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(0)))   # Upwind for w < 0

            # ------------------------------------------------------------------------
            # 3) Advective Flux Calculation
            # ------------------------------------------------------------------------
            # The flux is split into two streams.
            # Flux_stream = MassFlux_stream * SpecificEnergy_stream
            #             = (ρ * a * w) * (0.5 * w^2)
            #             = 0.5 * ρ * a * w^3
            #
            # Note: w_dn is negative, so w_dn^3 is negative, correctly indicating 
            # downward transport of energy.
            # ------------------------------------------------------------------------
            @. ρatke_convective_advection = -(
                ∇c(
                    wvec(
                        k_adv * 0.5 * ρ_f * (
                            # --- Updraft Contribution (w_up > 0) ---
                            LBF(w_up)^3 * Ifx(a_up)

                            # --- Downdraft Contribution (w_dn < 0) ---
                            + RBF(w_dn)^3 * Ifx(FT(1) - a_up)
                        )
                    )
                    - wvec(
                        # --- Second-order Diffusion Correction ---
                        # Standard numerical diffusion term to stabilize the advection
                        k_adv * 0.5 * ρ_f * ᶠinterp_a(a_en) * Δt *
                        0.5 * (Ifx(w_up)^2 + Ifx(w_dn)^2) *
                        ∇f_tke(ρatke_convective / (ρ_c * a_en))
                    )
                )
            )


        
            # Add some diffusion
            # @. ρatke_convective_advection += -∇c(wvec(
            #     -ρ_f * ᶠinterp_a(a_en) * min(
                    
            #         # --- TOTAL DIFFUSIVITY ---
            #         # We use the boosted velocity we just calculated.
            #         (30 * Ifx(sqrt(max(w_convective, FT(0))))) 
            #     ) *
            #     ∇f_tke(ρatke_convective / (ρ_c * a_en))
            # ))


        else # :: Diffusion (more stable in principle) K reduction in stable areas
            # ========================================================================================================= #
            aux_tc_f = face_aux_turbconv(state)

            # --- 1. Define the mixing length ---
            ℓ_mix = FT(500.0) 
            
            # --- 2. Define the necessary OPERATOR with boundaries ---
            #
            # THIS IS YOUR FIX: Defining ∇f with the correct wvec
            # boundary condition, which solves the MethodError.
            #
            ∇f_tke = CCO.GradientC2F(bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))
            
            # --- 3. Allocate ONLY the necessary temporary array ---
            # We still need w_conv_f because it's used twice
            w_conv_f = aux_tc_f.temporary_f1 

            # --- 4. Pre-calculate the one repeated value ---
            @. w_conv_f = Ifx(w_convective)


            # --- 4. Calculate Final Tendency (Fully Fused) ---
            K_floor = FT(0.1)

            # @. ρatke_convective_advection = -∇c(wvec(
            #     -ρ_f * ᶠinterp_a(a_en) * # Density * Area
            #     min(
            #         # --- TOTAL DIFFUSIVITY ---
            #         # 1. Convective (Goes to 0 cleanly if w=0)
            #         (k_adv * (ℓ_mix * w_conv_f^2) / 
            #         (max(w_conv_f, FT(0)) + ℓ_mix * Ifx(sqrt(max(aux_tc.∂b∂z, zero(FT)))) + eps(FT)))
                    
            #         # 2. Environmental
            #         + Ifx(aux_tc.KM) # prone to striping
                    
            #         # 3. Floor (The new safety net)
            #         + K_floor,

            #         FT(500.0)
            #     ) *
            #     ∇f_tke(ρatke_convective / (ρ_c * a_en))
            # ))


            # ======================================================================== #
            # 1. SETUP & LOCAL CONSTANTS
            # ======================================================================== #
            aux_tc_f = face_aux_turbconv(state)
            
            g_grav  = TCP.grav(param_set)
            c_p     = TCP.cp_d(param_set)
            θ_ref   = FT(300.0)
            tau     = FT(300.0)
            ℓ_mix   = FT(500.0) # Local definition for diffusion length scale

            # --- Scaling Coefficient for Latent Heat ---
            scaling_coeff = (2 * g_grav * ℓ_mix * tau) / (θ_ref * c_p)

            # ======================================================================== #
            # 2. CALCULATE EFFECTIVE VELOCITY (AT CENTERS)
            #    We use 'temporary_2' as a scratchpad to build the full velocity.
            #    We strictly avoid Indices 1 and 3.
            # ======================================================================== #
            
            # A. Start with the "Kick" (Latent Heat contribution)
            #    w_eff = sqrt(Q * C)
            w_eff = aux_tc.temporary_2
            Q_latent = aux_en.latent_heating_pos + -aux_en.latent_heating_neg
            @. w_eff = clamp(sqrt(abs(Q_latent) * scaling_coeff), FT(0), FT(6.0))

            # B. Add the Convective Velocity (Vector Sum)
            #    w_eff = sqrt(w_conv^2 + w_kick^2)
            @. w_eff = sqrt(w_convective^2 + w_eff^2)

            # C. Apply Stability Suppression
            #    This reduces mixing in very stable (high Ri) regions.
            #    Formula: w_eff_final = w_total * (k_adv / sqrt(1 + Ri))
            @. w_eff = w_eff * (
                k_adv / sqrt(1 + ((sqrt(max(aux_tc.∂b∂z, 0)) * ℓ_mix) / max(w_eff, eps(FT)))^2)
            )

            # ======================================================================== #
            # 3. INTERPOLATE TO FACES
            #    Diffusion happens at faces, so we move the result from Step 2.
            # ======================================================================== #
            
            w_eff_f = aux_tc_f.temporary_f1
            @. w_eff_f = Ifx(w_eff)

            # ======================================================================== #
            # 4. CALCULATE FLUX (DIFFUSION)
            #    K = w_eff_f * ℓ_mix
            # ======================================================================== #

            ∇f_tke = CCO.GradientC2F(bottom = CCO.SetGradient(wvec(FT(0))), top = CCO.SetGradient(wvec(FT(0))))
            K_floor = FT(0.1)

            @. ρatke_convective_advection = -∇c(wvec(
                -ρ_f * ᶠinterp_a(a_en) * min(
                    
                    # --- TOTAL DIFFUSIVITY ---
                    # We use the boosted velocity we just calculated.
                    (w_eff_f * ℓ_mix) 
                    
                    # + Environmental contribution (optional, currently 0)
                    + Ifx(aux_tc.KM) * 0 
                    
                    # + Floor
                    + K_floor,

                    FT(500.0) # Max Cap
                ) *
                ∇f_tke(ρatke_convective / (ρ_c * a_en))
            ))
            # # ========================================================================================================= #

        end


        # # === dissipation === #

        # === Dissipation (Split Physical Closures) === #
        ρatke_convective_dissipation = aux_en.tke_convective_dissipation # for some reason these have ρa in ρatke but are named this way so we follow

        # 1. Constants & Coefficients
        ℓ_mix = FT(500.0)
        k_self_diss = edmf.convective_tke_handler.self_dissipation_coeff # Controls Viscous Drag
        k_diss      = edmf.convective_tke_handler.dissipation_coeff      # Controls Stability Damping

        # 2. Calculate Dissipation (Fully Fused)
        #    Inputs: w_convective (√2k), ∂b∂z (N^2)
        #    Formula: ε = (Drag + Buoyancy) * Mass
        #
        #    Term 1 (Drag): k_self_diss * (k^1.5 / ℓ)
        #    Term 2 (Buoy): k_diss * k * N
        #
        #    (Using identity: k = 0.5 * w^2)
        
        # @. ρatke_convective_dissipation = (ρ_c * area_en) * (0.5 * w_convective^2) * (
        #     # A. Viscous Drag (Inertial Cascade)
        #     #    Scales with w / ℓ
        #     (k_self_diss * (w_convective / sqrt(FT(2))) / ℓ_mix) + # this is kolgomorov, too strong
            
        #     # B. Buoyancy Penalty (Gravity Work)
        #     #    Scales with N (independent of ℓ)
        #     (k_diss * sqrt(max(aux_tc.∂b∂z, 0)))
        # )

        τ_entrainment = FT(60) # really could just get rid of it and merge with k_self_diss...
        @. ρatke_convective_dissipation = (
            # A. Entrainment/Detrainment Mixing (Scales with TKE)
            (k_self_diss * ρatke_convective / τ_entrainment) +
            
            # B. Buoyancy Damping (Work against Gravity)
            # (Remains unchanged as it depends on N, not w)
            # (ρatke_convective * k_diss * sqrt(max(aux_tc.∂b∂z, 0))) # this wrongly only relies on dry buoyancy gradient
            ρatke_convective * k_diss * min(sqrt(max(aux_en.∂MSE∂z * g/(c_p * aux_en.T), zero(FT))), sqrt(max(aux_tc.∂b∂z, zero(FT)))) # use moist gradient for buoyancy damping
        )


        if any(!isfinite, ρatke_convective_production)
            @warn "Non-finite values in convective TKE production limiter computation."
            error("KE_current = $(parent(max.(ρatke_convective ./ (ρ_c .* a_en))));
                ΔKE_timestep = $(parent((ρatke_convective_production .- ρatke_convective_dissipation .+ ρatke_convective_advection) .* Δt ./ (ρ_c .* a_en))); 
                tke_max = $(parent(tke_max)); tke = $(parent(aux_en.tke))")
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

    # shear should for non tke vars contain gradient stuff. but doesn't contain convective tke contribution
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
    Ic = is_tke ? CCO.InterpolateF2C() : Base.identity
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



