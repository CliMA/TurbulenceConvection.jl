### 

# I think for our single-moment we still need sedimentation....

# we may only need sedimentation for cloud ice... (e.g. as in Lohmann et al 2018 The Importance)
# These will need to be called on entire profiles in EDMF_Environment.jl, EDMF_Updraft.jl

function calculate_sedimentation_sources(
    param_set::APS,
    grid::Grid,
    ρ_c,  # figure out types for these
    ts; # figure out types for these
    w = 0.0, # background wind... (might still change this to operate inside update_aux or something? idk...)
    area = 1.0,
    integration_method::Symbol = :upwinding,
    liq_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Blk1MVel,
    ice_velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Chen2022Vel,
    grid_mean::Bool = true, # whether or not this is a grid mean calc or if we have subdomains...
    liq_Dmax::FT2 = Inf, # maximum diameter for liquid particles
    ice_Dmax::FT2 = 62.5e-6, # maximum diameter for ice particles
    liq_scaling_factor::FT2 = 1.0, # scaling factor for terminal velocity
    ice_scaling_factor::FT2 = 1.0, # scaling factor for terminal velocity
    # Nl::Union{FT2, AbstractArray{FT2}, Nothing} = nothing, # N for particle size distrobution
    # Ni::Union{FT2, AbstractArray{FT2}, Nothing} = nothing, # N for particle size distrobution
    Nl::Union{FT2, AbstractArray{FT2}} = NaN, # N for particle size distrobution # testing for type stability, use NaN instead of nothing
    Ni::Union{FT2, AbstractArray{FT2}} = NaN, # N for particle size distrobution # testing for type stability, use NaN instead of nothing
) where {FT2}

    #= Note --
    I'm not sure how to deal w/ density here
    For dycore.jl subsidence, we assumed the overall density was constant,
    Here q is specific humidity (kg/kg) so it's the same?

    We make the same assumption that expansion contraction will handle the ρ part, and we jus need q not ρq...


    However, I feel like area should matter right?
    Like going down from an area with 1 area fraction to one with .5 should increase q just to squeeze in?

    That sounds weird though... it feels like it should be only on the grid mean or something...
    Lke if the updraft gets really wide at the top, the sedimentation should fall into the environment too not just the updraft...
    -- maybe ask Anna about that

    Also, I'm not sure we should not be balancing w/ density i'm not sure droplets compress in sedimentation...
    precip doesn't right?
    In EDMF_Precipitation.jl they do weight by ρ_c before doing that calculation... and then divide out... make sense
    =#

    FT = eltype(ρ_c) # don't have state so not sure how to get it lol
    # area = FT.(area) # convert to FT
    thermo_params = TCP.thermodynamics_params(param_set)


    local velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type}
    local Dmax::FT
    local scaling_factor::FT
    # local Nt::Union{FT2, AbstractArray{FT2}, Nothing}
    local Nt::Union{FT2, AbstractArray{FT2}} # testing for type stability, use NaN instead of nothing


    # TS = typeof(ρ_c)
    # create vector of zeros for sedimentation sources
    S_ql = ρ_c .* 0 # I dont know how to create zeros lol
    S_qi = ρ_c .* 0 # I dont know how to create zeros lol

    S_ql_other = ρ_c .* 0 # I dont know how to create zeros lol
    S_qi_other = ρ_c .* 0 # I dont know how to create zeros lol


    if get_isbits_nt(param_set.user_args, :use_sedimentation, false)


        # @info(ts)
        # @info TD.PhasePartition.(Ref(thermo_params), ts)
        q = TD.PhasePartition.(thermo_params, ts)
        q_ice = (q -> q.ice).(q)
        q_liq = (q -> q.liq).(q)

        wvec = CC.Geometry.WVector
        ∇c = CCO.DivergenceF2C() # F2C to come back from C2F

        UBsed = CCO.UpwindBiasedProductC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # upwinding, extrapolate bc we don't know the boa/toa derivatives


        for (q_type, q, S, S_other, velo_scheme, Dmax, scaling_factor, Nt) in zip(
            (liq_type, ice_type),
            (q_liq, q_ice),
            (S_ql, S_qi),
            (S_ql_other, S_qi_other),
            (liq_velo_scheme, ice_velo_scheme),
            (liq_Dmax, ice_Dmax),
            (liq_scaling_factor, ice_scaling_factor),
            (Nl, Ni),
        )

            # w_sed = calculate_sedimentation_velocity_test.(q, q_type) # this will be on centers not faces then... just like subsidence was

            microphys_params = TCP.microphysics_params(param_set)
            w_sed =
                .-calculate_sedimentation_velocity.(
                    microphys_params,
                    q,
                    ρ_c,
                    q_type,
                    Nt;
                    velo_scheme = velo_scheme,
                    Dmax = Dmax,
                ) # this will be on centers not faces then... just like subsidence was [ waiting till we figure out stability problems :) ]
            w_sed .*= scaling_factor # scale terminal velocity



            # @info max(w_sed[:], q_type, scaling_factor, q)
            # @info("w", (w_sed[:][190], scaling_factor, q_type))
            # negative bc terminal_velocity is defined as positive down, still Right Biased but now also -∇, so matches w, not precip terminal_velocity... better for upwinding I think? bc we have z=0 at sfc

            #=
            - In principle we could have add w tendency to w_sed tendency, (and that should work?)
            - However, then I think your exchange(to_other) comes out wrong? bc that should be based on w+w_sed, not w_sed bc area is also advected by w
            - or maye just the exchange alone is based on w-w_sed and then ∇(ρqaw_sed) + ∇(ρqaw) combine as desired in EDMF_Functions?
            - maybe just using w+w_sed is better w/ upwinding than trying to do both separately? to_other has to be upwinded anyway right?
            =#

            # w_sed .+= w # these must combine, though maybe we can just copy EDMF_Functions method? though there's no EDMF_Functions for Environment
            # This, however, creates an w direction problem... (go back to upwinding? idk...)


            # also this is still bad bc we'd need to remove the w_up in EDMF_Functions anyway now that it's included here... to not double count...

            q .*= ρ_c # convert to absolute amount (area and density weight so fluxes are accurate)

            if !grid_mean
                w_sed_to_other = copy(w_sed) # make a copy for the exchange calculation (not sure if this is enough, maybe deepcopy() or w_sed_to_other = w_sed .* 1.0 or something is better?)
                # for k in real_center_indices(grid)
                # w_sed_to_other[k] += w[k.i] # add background wind (can't figure out how to add to F2C direclty? idk...)  w_sed .+= w didn't work..., minus bc w_sed is defined towards surface, but w is up
                # end
            end

            if integration_method == :upwinding
                # =================================================================================== #
                # UPWINDING (NECSSARY HERE bc w_w_sed could have either sign right...) -- seems unstable...
                # get boa and toa values
                kc_surf = kc_surface(grid)
                q_boa = q[kc_surf]
                kc_toa = kc_top_of_atmos(grid)
                q_toa = q[kc_toa]
                # sedimentation should be face valued so we'll need a C2F call
                C2Fsed = CCO.InterpolateC2F(; bottom = bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # just extrapolate, it's just sedimentation velocity
                C2Fq = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(q_boa)), top = CCO.SetValue(FT(q_toa))) # not sure if this should be right biased or not

                # we to C2F, then ∇c goes F2C, then UBsed got C2F, but we wanna end on C [ there's no UpwindBiasedProductF2C and we had UpwindBiasedProductF2C(u[F], x[C])
                F2Csed = CCO.InterpolateF2C(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # We might have put 0 for no penetration on boundary, but that's not exactly true in our dataset...
                CV32FT = x -> x[1] # convert Contravariant3Vector to Float64

                # regular upwinding 
                wpaq = @. UBsed(wvec(C2Fsed(w_sed)), q * area) # w_sed is towards surface, so we want to upwind w_sed * q
                ∇wρaq = @. -∇c(wvec(wpaq)) # maybe this has to be done after CV32FT for the contravariant 3 vector?
                # ∇wρaq = @. UBsed(wvec(C2Fsed(w_sed)), ∇wρaq) # works, but output type is different than tendencies type
                # F2Csed = CCO.InterpolateF2C(; bottom = CCO.Extrapolate(), top = 0) # no input from top, just extrapolate at sfc
                # ∇wρaq = @. F2Csed(∇wρaq) # convert back to C
                # ∇wρaq = @. CV32FT(∇wρaq) # convert back to Float64 from Contravariant3Vector

                # flux from one partition to the other
                if !grid_mean
                    to_other = @. ∇wρaq * 0
                    _area_top = isa(area, Number) ? area : area[kc_toa]
                    _area_bottom = isa(area, Number) ? area : area[kc_surf]
                    C2Fa =
                        CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(_area_bottom)), top = CCO.SetValue(FT(_area_top)))
                    ∇a = @. ∇c(wvec(C2Fa(area))) # seems stable w/ either sign w...? (Left bias this?)

                    prod = @. q * ∇a * (∇a > 0)
                    to_other = @. -UBsed(wvec(C2Fsed(w_sed_to_other)), prod)  # upwinded
                    to_other = @. F2Csed(to_other) # convert back to C
                    to_other = @. CV32FT(to_other) # convert back to Float64 from Contravariant3Vector
                    @. ∇wρaq -= to_other # remove from source


                    # for k in real_center_indices(grid)
                    #     if (∇a[k] > 0) # &&  (w_sed_to_other[k] > 0) # area is decreasing towards surface and we have a positive net w_sed (towards surface), otherwise it's either getting lofted inside it's regime or falling into the same regime
                    #         to_other[k] = -w_sed_to_other[k] * q[k] * ∇a[k]  # flux to other (if ∇a is large for example, send most of the source to the other side):/
                    #         ∇wρaq[k] -= to_other[k] # remove from source
                    #     end
                    # end
                end

            elseif integration_method == :right_biased
                # =================================================================================== #
                # Right Biased (copy from EDMF_Precipitation) -- seems stable
                # get toa values
                kc_toa = kc_top_of_atmos(grid)
                q_toa = q[kc_toa]
                # right biased operators
                RB = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(q_toa)))
                ∇ = CCO.DivergenceF2C(; bottom = CCO.Extrapolate())
                ∇wρaq = @. -∇(wvec(RB(q * w_sed * area))) # seems stable w/ either sign w...?

                # surely this should be upwinded.... bc w_sed_to_other could have either sign right?
                if !grid_mean
                    to_other = ∇wρaq .* 0
                    _area_top = isa(area, Number) ? area : area[kc_toa]
                    RBa = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(_area_top)))
                    ∇a = @. ∇(wvec(RBa(area))) # seems stable w/ either sign w...? (Left bias this?)

                    for k in real_center_indices(grid)
                        if ∇a[k] > 0  # && w_sed_to_other[k] > 0 # area is decreasing towards surface and we have a positive net w_sed (towards surface), otherwise it's either getting lofted inside it's regime or falling into the same regime
                            to_other[k] = -w_sed_to_other[k] * q[k] * ∇a[k]  # flux to other (if ∇a is large for example, send most of the source to the other side):/
                            if w_sed_to_other[k] > 0
                                @info "w_sed_to_other[k] > 0 at k: $k, w_sed_to_other[k]: $(w_sed_to_other[k])"
                            end

                            # if ∇wρaq[k] < 0 && to_other[k] > 0 # this is ok, it just exacerbates the rate of decline.
                            #     @info "∇wρaq[k] < 0 but to_other > 0 at k: $k, (∇wρaq[k], to_other[k]): ($(∇wρaq[k]), $(to_other[k]))"
                            # end
                            if to_other[k] < 0 # this would be bad, can't happen...
                                @info "to_other < 0 at k: $k, to_other[k]: $(to_other[k])"
                            end


                            ∇wρaq[k] -= to_other[k] # remove from source



                        end
                    end
                end
                # =================================================================================== #
            else
                error("Method not implemented")
            end


            for k in real_center_indices(grid)
                # _area = isa(area, Number) ? area : area[k]
                S[k] = ∇wρaq[k] / ρ_c[k] # undo density weighting
                if !grid_mean
                    S_other[k] = to_other[k] / ρ_c[k] # divide source by local area, undo density weighting, negative cause source is negative of the divergence
                end
            end



            # # get index of maximum
            # if q_type == liq_type 
            #     # kmax = argmin([w_sed[k] for k in real_center_indices(grid)])
            #     kmax = argmax([q[k] for k in real_center_indices(grid)]) 
            #     kmax = real_center_indices(grid)[kmax] # convert to real index
            #     @info ("kmax w: $(kmax.i), w sed: $(w_sed[kmax]), q: $(q[kmax]), ρc: $(ρ_c[kmax]), scaling factor: $scaling_factor, ∇wρaq: $(∇wρaq[kmax]), to other: $(to_other[kmax]), area: $(area[kmax])")

            #     kmax = argmax([∇wρaq[k] for k in real_center_indices(grid)])
            #     kmax = real_center_indices(grid)[kmax] # convert to real index
            #     @info ("kmax ∇wρaq: $(kmax.i), w sed: $(w_sed[kmax]), q: $(q[kmax]), ρc: $(ρ_c[kmax]), scaling factor: $scaling_factor, ∇wρaq: $(∇wρaq[kmax]), to other: $(to_other[kmax]), area: $(area[kmax])")

            #     kmin = argmin([∇wρaq[k] for k in real_center_indices(grid)])
            #     kmin = real_center_indices(grid)[kmin] # convert to real index
            #     @info ("kmin ∇wρaq: $(kmin.i), w sed: $(w_sed[kmin]), q: $(q[kmin]), ρc: $(ρ_c[kmin]), scaling factor: $scaling_factor, ∇wρaq: $(∇wρaq[kmin]), to other: $(to_other[kmin]), area: $(area[kmin])")

            #     kmax = argmax([S[k] for k in real_center_indices(grid)])
            #     kmax = real_center_indices(grid)[kmax] # convert to real index
            #     @info ("kmax S: $(kmax.i), w sed: $(w_sed[kmax]), q: $(q[kmax]), ρc: $(ρ_c[kmax]), scaling factor: $scaling_factor, S: $(S[kmax]), S_other: $(S_other[kmax]), area: $(area[kmax])")

            #     kmin = argmin([S[k] for k in real_center_indices(grid)])
            #     kmin = real_center_indices(grid)[kmin] # convert to real index
            #     @info ("kmin S: $(kmin.i), w sed: $(w_sed[kmin]), q: $(q[kmin]), ρc: $(ρ_c[kmin]), scaling factor: $scaling_factor, S: $(S[kmin]), S_other: $(S_other[kmin]), area: $(area[kmin])")

            #     println("------------------------------------------------")
            # end


        end

    end
    return NoneqMoistureSources{FT}.(S_ql, S_qi), NoneqMoistureSources{FT}.(S_ql_other, S_qi_other)

end


"""
Calculate cloud liq/ice sedimentation velocity
Formula from:

Blk1MVel is what clima_1m uses

This might be in terms of dz/dt so w < 0 is down, w > 0 is up... not sure...
If so gotta check the sign on sedimentation lol...

Is also still unstable even w/ upwinding... maybe right biased will be stable?

"""
# function calculate_sedimentation_velocity_test(
#     q::FT,
#     q_type::TD.Phase
#     ) where {FT}

#     # I think put positive numbers bc that's what terminval_velocity returns
#     # S = -∇(wq) but we effectively put -w instead of w so we then also use S = ∇(wq)

#     if q_type == TD.Liquid()
#         w_sed = FT(0.005)
#     elseif q_type == TD.Ice()
#         w_sed = FT(0.05)
#     end
#     return w_sed # test
# end

#=
NOTE: Our signatures for CM1.terminal_velocity are out of date, try to use definitions from https://github.com/CliMA/CloudMicrophysics.jl/blob/v0.14.0/

liq_type, ice_type, rain_type, snow_type are defined in src/TurbulenceConvection.jl as is Blk1MVel
=#

function calculate_sedimentation_velocity(
    microphys_params::CM.Parameters.AbstractCloudMicrophysicsParameters, # Abstract Parameter Set
    q::FT,
    ρ::FT,
    ::CMT.LiquidType,
    # Nt::Union{FT, Nothing}; # N for particle size distrobution;
    Nt::FT; # N for particle size distrobution; # testing for type stability, use NaN instead of nothing
    velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Blk1MVel,
    Dmax::FT = Inf, # maximum diameter for ice particles
) where {FT}

    return FT(0.0) # not implemented in CloudMicrophysics.jl 0.14

    # if velo_scheme == :Blk1MVel
    #     return FT(0.0) # CM1.terminal_velocity(microphys_params, liq_type, Blk1MVel, ρ, q) does not exist
    # else
    #     error("velo_scheme $velo_scheme not implemented")
    # end
end


function calculate_sedimentation_velocity(
    microphys_params::CM.Parameters.AbstractCloudMicrophysicsParameters,
    q::FT,
    ρ::FT,
    ::CMT.IceType,
    # Nt::Union{FT, Nothing}; # N for particle size distrobution;
    Nt::FT; # N for particle size distrobution; # testing for type stability, use NaN instead of nothing
    velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Chen2022Vel,
    Dmax::FT = Inf, # maximum diameter for ice particles
) where {FT}

    if velo_scheme == Chen2022Vel
        # w  = CM1.terminal_velocity(microphys_params, ice_type, Chen2022Vel, ρ, q) # testing
        w = my_terminal_velocity(microphys_params, ice_type, velo_scheme, ρ, q; Dmax = Dmax, Nt = Nt) # testing
    # I think this is correct, this version has no inkling of diameter/radius
    # Newer versions I believe do but I don't have them..., ask Anna
    elseif velo_scheme == Blk1MVel
        w = CM1.terminal_velocity(microphys_params, ice_type, velo_scheme, ρ, q) # does this velo_scheme even exist for ice_type?
    # no Dmax for nonchen afaik
    else
        error("velo_scheme $velo_scheme not implemented")
    end

    return isnan(w) ? FT(0.0) : w
end
