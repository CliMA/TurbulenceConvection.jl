### 

# I think for our single-moment we still need sedimentation....

# we may only need sedimentation for cloud ice... (e.g. as in Lohmann et al 2018 The Importance)
# These will need to be called on entire profiles in EDMF_Environment.jl, EDMF_Updraft.jl

# function swap_eltype(::Type{CC.DataLayouts.VF{OldS, Nv, OldA}}, ::Type{NewS}) where {OldS, Nv, OldA, NewS}
#     # Extract the array wrapper (Array, CuArray, etc.)
#     ArrayWrapper = Base.typename(OldA).wrapper
#     # Reconstruct: Wrapper{NewElement, 2 dims for VF}
#     NewA = ArrayWrapper{NewS, 2}
#     return CC.DataLayouts.VF{NewS, Nv, NewA}
# end
# # ADD THIS: Overload for Fields
# # This catches Field{VF, Space} and swaps the VF inside
# function swap_eltype(::Type{CC.Fields.Field{V, S}}, ::Type{NewS}) where {V, S, NewS}
#     NewV = swap_eltype(V, NewS) # Calls the function above
#     return CC.Fields.Field{NewV, S}
# end

"""
This function operates not element wise but on entire CC Fields because it uses vertical gradients
"""
function calculate_sedimentation_sources!(
    sedimentation::SCF,
    sedimentation_other::SCF,
    param_set::APS,
    # q_type::Union{CMT.LiquidType, CMT.IceType}, # deprecate bc we no longer dispatch on this
    ρ::CC.Fields.Field,  # figure out types for these
    q::CC.Fields.Field, # figure out types for these
    w_sed::CC.Fields.Field,
    w::CC.Fields.Field,
    area::CC.Fields.Field,
    grid::Grid
    ;
    # ρ_q::FT;  # figure out types for these # not used
    # w::Union{CC.Fields.Field, FT} = 0.0, # background wind... (might still change this to operate inside update_aux or something? idk...)
    differencing_scheme::AbstractIntegrationScheme = UpwindDifferencingScheme(),
    # velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Chen2022Vel, # deprecated
    grid_mean::Bool = false, # whether or not this is a grid mean calc or if we have subdomains...
    use_relative_w::Bool = false,
    scratch1::CF,
    scratch2::CF,
    scratch3::CF,
    scratch4::CF,
    scratch1F::FF,
    scratch2F::FF
) where {SCF <: CC.Fields.Field, CF <: CC.Fields.Field, FF <: CC.Fields.Field}
 #::Tuple{swap_eltype(CF, NoneqMoistureSources{eltype(CF)}), swap_eltype(CF, NoneqMoistureSources{eltype(CF)})} where {CF<:CC.Fields.Field, FF<:CC.Fields.Field}


    #= Note --
    I'm not sure how to deal w/ density here
    For dycore.jl subsidence, we assumed the overall density was constant,
    Here q is specific humidity (kg/kg) so it's the same?

    We make the same assumption that expansion contraction will handle the ρ part, and we jus need q not ρq...

    we used to have w_sed = -calculate_sedimentation_velocity(...) but we've flipped now so we use -w_sed where applicable

    However, I feel like area should matter right?
    Like going down from an area with 1 area fraction to one with .5 should increase q just to squeeze in?

    That sounds weird though... it feels like it should be only on the grid mean or something...
    Lke if the updraft gets really wide at the top, the sedimentation should fall into the environment too not just the updraft...
    -- maybe ask Anna about that

    Also, I'm not sure we should not be balancing w/ density i'm not sure droplets compress in sedimentation...
    precip doesn't right?
    In EDMF_Precipitation.jl they do weight by ρ before doing that calculation... and then divide out... make sense
    =#

    # FT = eltype(ρ) # don't have state so not sure how to get it lol
    FT = eltype(param_set) # Do we need this? or just asume liq/ice_Dmax is guaranteed to give us the type?
    # thermo_params::TDPS = TCP.thermodynamics_params(param_set)
    # microphys_params::ACMP = TCP.microphysics_params(param_set)


    # create vector of zeros for sedimentation sources
    # S_q = ρ .* 0 # I dont know how to create zeros lol (similar doesn't fill -- do we need it filled?)
    # S_other = ρ .* 0 # I dont know how to create zeros lol
    # S = similar(ρ) # a copy but don't actually fill in
    # S_other = similar(ρ) # a copy but don't actually fill in

    # S = scratch1 # reuse temporary storage
    # S_other = scratch2 # reuse temporary storage

    S = sedimentation
    S_other = sedimentation_other

    UBsed = CCO.UpwindBiasedProductC2F(; bottom = CCO.Extrapolate(), top = CCO.SetValue(FT(0))) # upwinding, extrapolate bc we don't know the boa/toa derivatives, no flux through the top

    #=
        Actually I think Extrapolate() is bad -- it forces the flux to be the same at cell faces because it interpolates inwards, see https://clima.github.io/ClimaCore.jl/dev/operators/#ClimaCore.Operators.DivergenceF2C
        This means there can never be flux convergence, and the value will never change.
        Because we extrapolate the terminal velocity, it is alr
        I think for the bottom cell, so for the bottom cell, assume the w is the same but 

        We could try to set up our own system to extrapolate q to the lower cell face for example, but that's hard...
        RB doesn't have any such problem... given the boundary has no vertical velocity BC, i think for the lowest point we should use the RB value.

        The flux would just be the flux at the middle, shifted down
    =#

    # negative bc terminal_velocity is defined as positive down, still Right Biased but now also -∇, so matches w, not precip terminal_velocity... better for upwinding I think? bc we have z=0 at sfc

    #=
    - In principle we could have add w tendency to w_sed tendency, (and that should work?)
    - However, then I think your exchange(to_other) comes out wrong? bc that should be based on w+w_sed, not w_sed bc area is also advected by w
    - or maye just the exchange alone is based on w-w_sed and then ∇(ρqaw_sed) + ∇(ρqaw) combine as desired in EDMF_Functions?
    - maybe just using w+w_sed is better w/ upwinding than trying to do both separately? to_other has to be upwinded anyway right?

    What you would really want to do if you wanted to separate things would be to only calculate using w_sed+w here and remove the advection term in EDMF_functions().
        - as written what we have is dispersive (handling up and down separately)... but that's ok I guess... for now. For very strong updrafts that might be bad but then the tendencies should be mismatched anwyay
        - for small updrafts light dispersion is probably ok given the uncertainty in the size distribution and varying terminal velocity with droplet size.
    =#


    # This, however, creates an w direction problem... (go back to upwinding? idk...)

    # also this is still bad bc we'd need to remove the w_up in EDMF_Functions anyway now that it's included here... to not double count...


    if !grid_mean
        # w_sed_to_other = copy(w_sed) # make a copy for the exchange calculation (not sure if this is enough, maybe deepcopy() or w_sed_to_other = w_sed .* 1.0 or something is better?)
        # w_sed_to_other = similar(w_sed) # a copy but don't actually fill in 
        w_sed_to_other = w_sed # rn, we're not actually addting background wind to w_sed, so this is fine we can reuse the same object... This also assumes the area is advecting with w..., which isn't exactly true... e.g. w_en may not see a spike in w_up coming so da/dt is unknowable, and thus the crosover rate. instantatenously, you probably want advection from en and up together 
        # w_sed_to_other[k] += w[k.i] # add background wind (can't figure out how to add to F2C direclty? idk...)  w_sed .+= w didn't work..., minus bc w_sed is defined towards surface, but w is up
        # end
    end


    if differencing_scheme isa UpwindDifferencingScheme
        # =================================================================================== #
        # UPWINDING (NECSSARY HERE bc w_w_sed could have either sign right...) -- seems unstable...
        # get boa and toa values
        kc_surf = kc_surface(grid)
        kf_surf = kf_surface(grid)
        # q_boa = q[kc_surf]
        kc_toa = kc_top_of_atmos(grid)
        # q_toa = q[kc_toa]
        # sedimentation should be face valued so we'll need a C2F call
        C2Fsed = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # just extrapolate, it's just sedimentation velocity
        # C2Fq = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(q_boa)), top = CCO.SetValue(FT(q_toa))) # not sure if this should be right biased or not

        # we to C2F, then ∇c goes F2C, then UBsed got C2F, but we wanna end on C [ there's no UpwindBiasedProductF2C and we had UpwindBiasedProductF2C(u[F], x[C])
        F2Csed = Ic # don't need bcs for F2C
        # CV32FT = x -> x[1] # convert Contravariant3Vector to Float64

        w_sed_f = scratch1F # reuse temporary storage
        if use_relative_w # seems to induce weird behavior around w_sed = 0...
            @. w_sed_f = C2Fsed(w_sed) - w # preserve full information in w
        else
            @. w_sed_f = C2Fsed(w_sed) # convert to face values
        end
        w_sed = w_sed_f # alias

        # regular upwinding 
        # wρaq = scratch2F
        # @. wρaq = UBsed(wvec(-w_sed), ρ * q * area) # -w_sed is towards surface, so we want to upwind w/ -w_sed  and do -∇c after instead of just skipping and cancelling the negatives 


        # replace the lowest cell with the right-biased value [[ it's guaranted w = 0, so the terminal velocity should be all... ]]

        # I deally we would use the boundary condition option in UBSub for this conversion. however, if the upper flux that UBsub pulls is zero, we wont be able to adjust the lower (otherwise we could use a scaling factor on the ratio between the value at k+1 and k and put that as our CCO.SetValue(...))
        # wρaq[kf_surf] = CCG.Contravariant3Vector(ρq[kc_surf] * area[kc_surf] * -w_sed[kf_surf]) # not quite right as far as i can tell... needs wvec first for physical units

        # replace the highest cell with 0
        # wρaq[kf_toa] = FT(0) # i think the bc has this alraedy

        # wvec_wρaq = @. wvec(wρaq) # convert to vector type for divergence [to get back to physical units]
        wvec_wρaq = @. wvec( UBsed(wvec(-w_sed), ρ * q * area) )
        scratch_wvec = wvec_wρaq # alias
        # @warn "typeof(wvec_wρaq) = $(typeof(wvec_wρaq))"
        # wvec_wρaq[kf_surf] = q[kc_surf] * area[kc_surf] * -w_sed[kf_surf]
        wvec_wρaq[kf_surf] = map(x -> ρ[kc_surf] * q[kc_surf] * area[kc_surf] * -w_sed[kf_surf], wvec_wρaq[kf_surf]) # i think this works, following the logic in # /home/jbenjami/.julia/packages/ClimaCore/vJw0m/src/Geometry/axistensors.jl. Also would have worked on the Contravariant3vector but the units would have been wrong (not physical w/o wvec)
        
        ∇wρaq = scratch1 # reuse temporary storage
        @. ∇wρaq = -∇c(wvec_wρaq) # maybe this has to be done after CV32FT for the contravariant 3 vector?

        # flux from one partition to the other
        if !grid_mean
            _area_top = isa(area, FT) ? area : area[kc_toa]
            _area_bottom = isa(area, FT) ? area : area[kc_surf]
            C2Fa = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(_area_bottom)), top = CCO.SetValue(FT(_area_top)))
            ∇a = scratch2
            @. ∇a = ∇c(wvec(C2Fa(area))) # seems stable w/ either sign w...? (Left bias this?)
            prod = scratch2 # in place
            @. prod = ifelse(∇a > 0, ρ * q * ∇a, FT(0)) # only send to other if area is increasing towards surface

            to_other_f = scratch_wvec # reuse wvec vector
            @. to_other_f = -wvec(UBsed(wvec(C2Fsed(-w_sed_to_other)), prod))  # upwinded, outer wvec important to get physical units
            to_other = scratch2 # reuse temporary storage

            # @warn "typeof(F2Csed.(to_other_f) = $(typeof(F2Csed.(to_other_f)))"
            # @. to_other = F2Csed(to_other_f).components.data.:1 # from Dennis # remove from source.  # convert back to C
            @. to_other = toscalar(F2Csed(to_other_f))  # from Dennis # remove from source.  # convert back to C
            @. ∇wρaq -= to_other
        end


    elseif differencing_scheme isa RightBiasedDifferencingScheme
        # =================================================================================== #
        # Right Biased (copy from EDMF_Precipitation) -- seems stable
        # get toa values
        kc_toa = kc_top_of_atmos(grid)
        # right biased operators
        RB = CCO.RightBiasedC2F(; top = CCO.Extrapolate())
        ∇ = CCO.DivergenceF2C(; bottom = CCO.Extrapolate())

        C2Fsed = CCO.InterpolateC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # just extrapolate, it's just sedimentation velocity
        F2Csed = Ic # don't need bcs for F2c

        w_sed = scratch1F # reuse temporary storage
        if use_relative_w
            @. w_sed = F2Csed(C2Fsed(w_sed) - w) # preserve full information in w, however for RB we need to go back to C values
        else
            @. w_sed = C2Fsed(w_sed) # convert to face values
        end

        ∇wρaq = scratch1
        @. ∇wρaq = ∇(wvec(RB(ρ * q * w_sed * area))) # this is the way they did rain and snow...

        # surely this should be upwinded.... bc w_sed_to_other could have either sign right?
        if !grid_mean
            # to_other = similar(∇wρaq) # a copy but don't actually fill in
            to_other = scratch2 # reuse temporary storage

            # _area_top = isa(area, Number) ? area : area[kc_toa]
            _area_top = isa(area, FT) ? area : area[kc_toa]
            RBa = CCO.RightBiasedC2F(; top = CCO.SetValue(FT(_area_top)))
            ∇a = similar(∇wρaq) # a copy but don't actually fill in
            @. ∇a = ∇(wvec(RBa(area))) # seems stable w/ either sign w...? (Left bias this?)

            @inbounds for k in real_center_indices(grid)
                if ∇a[k] > 0  # && w_sed_to_other[k] > 0 # area is decreasing towards surface and we have a positive net w_sed (towards surface), otherwise it's either getting lofted inside it's regime or falling into the same regime
                    to_other[k] = w_sed_to_other[k] * q[k] * ∇a[k]  # flux to other (if ∇a is large for example, send most of the source to the other side):/
                    ∇wρaq[k] -= to_other[k] # remove from source
                end
            end
        end
        # =================================================================================== #
    else
        error("Method not implemented")
    end


    @inbounds for k in real_center_indices(grid)
        # _area = isa(area, Number) ? area : area[k]
        S[k] = ∇wρaq[k] / ρ[k] # undo density weighting
        if !grid_mean
            S_other[k] = to_other[k] / ρ[k] # divide source by local area, undo density weighting, negative cause source is negative of the divergence
        end
    end


    return nothing
    # return NoneqMoistureSource{FT}.(S), NoneqMoistureSource{FT}.(S_other)

end


"""
Calculate cloud liq/ice sedimentation velocity
Formula from:

Blk1MVel is what clima_1m uses

This might be in terms of dz/dt so w < 0 is down, w > 0 is up... not sure...
If so gotta check the sign on sedimentation lol...

Is also still unstable even w/ upwinding... maybe right biased will be stable?

"""

#=
NOTE: Our signatures for CM1.terminal_velocity are out of date, try to use definitions from https://github.com/CliMA/CloudMicrophysics.jl/blob/v0.14.0/

liq_type, ice_type, rain_type, snow_type are defined in src/TurbulenceConvection.jl as is Blk1MVel
=#

function calculate_sedimentation_velocity(
    # microphys_params::CM.Parameters.AbstractCloudMicrophysicsParameters, # Abstract Parameter Set
    param_set::APS,
    q::FT,
    ρ::FT,
    ::CMT.LiquidType,
    # Nt::Union{FT, Nothing}; # N for particle size distrobution;
    Nt::FT; # N for particle size distribution; # testing for type stability, use NaN instead of nothing
    velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Blk1MVel,
    Dmax::FT = FT(Inf), # maximum diameter for ice particles
    rain_type::CMT.RainType = rain_type, # rain_type is available from a a global, we dispatch this way bc liq_type doesn't have terminal velocity defined
) where {FT <: Real}

    Nt = isnan(Nt) ? FT(param_set.user_params.N_CCN) : Nt # default to 250 / cm^-3 as that's a reasonable rain amount, do not let the default value for rain win as you'll get something ridiculously high. Even this might need to be scaled down, idk.
    # Nt = isnan(Nt) ? FT(250*1e6) : Nt # default to 250 / cm^-3 as that's a reasonable rain amount, do not let the default value for rain win as you'll get something ridiculously high. Even this might need to be scaled down, idk.

    if velo_scheme == Chen2022Vel
        w = my_terminal_velocity(param_set, rain_type, velo_scheme, ρ, q; Dmax = Dmax, Nt = Nt) # testing [ this version makes you pass in raintype to avoid creating another object] [ rain_type is avilable from a global]
    elseif velo_scheme == Blk1MVel
        error("velo_scheme $velo_scheme not implemented for liquid_type in current CM1 version. If you want liquid sedimentation, use Chen2022Vel (TODO: turn this into an error at namelist construction)")
        # return FT(0.0) # not implemented in CloudMicrophysics.jl 0.14
    end

    return (isnan(w) || (w < 0)) ? zero(FT) : w # nans are bad, also for some really small numbers can return negative values... so just set to 0

end


function calculate_sedimentation_velocity(
    # microphys_params::CM.Parameters.AbstractCloudMicrophysicsParameters,
    param_set::APS,
    q::FT,
    ρ::FT,
    ice_type::CMT.IceType,
    # Nt::Union{FT, Nothing}; # N for particle size distrobution;
    Nt::FT; # N for particle size distrobution; # testing for type stability, use NaN instead of nothing
    velo_scheme::Union{CMT.Blk1MVelType, CMT.Chen2022Type} = Chen2022Vel,
    Dmax::FT = FT(Inf), # maximum diameter for ice particles
) where {FT <: Real}

    if velo_scheme == Chen2022Vel
        # w  = CM1.terminal_velocity(microphys_params, ice_type, Chen2022Vel, ρ, q) # testing
        w = my_terminal_velocity(param_set, ice_type, velo_scheme, ρ, q; Dmax = Dmax, Nt = Nt) # testing
    # I think this is correct, this version has no inkling of diameter/radius
    # Newer versions I believe do but I don't have them..., ask Anna
    elseif velo_scheme == Blk1MVel
        # w = CM1.terminal_velocity(microphys_params, ice_type, velo_scheme, ρ, q) # does this velo_scheme even exist for ice_type?
        error("velo_scheme $velo_scheme not implemented for ice_type in current CM1 version. If you want ice sedimentation, use Chen2022Vel (TODO: turn this into an error at namelist construction)")
    # no Dmax for nonchen afaik
    else
        error("velo_scheme $velo_scheme not implemented")
    end

    return (isnan(w) || (w < 0)) ? FT(0.0) : w # nans are bad, also for some really small numbers can return negative values... so just set to 0
    # return isnan(w) ? FT(0.0) : w # nans are bad, also for some really small numbers can return negative values... so just set to 0

end

