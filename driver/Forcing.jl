initialize(::ForcingBase, grid, state) = nothing

function initialize(::ForcingBase{ForcingLES}, grid, state, LESDat::LESData)
    aux_gm = TC.center_aux_grid_mean(state)
    nt = NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        zc_les = Array(TC.get_nc_data(data, "zc"))
        getvar(var) = TC.pyinterp(vec(grid.zc.z), zc_les, TC.mean_nc_data(data, "profiles", var, imin, imax))
        get_nudgevar(var) = TC.pyinterp(vec(grid.zc.z), zc_les, TC.init_nc_data(data, "profiles", var))

        dTdt_hadv = getvar("dtdt_hadv")
        H_nudge = get_nudgevar("thetali_mean")
        dTdt_fluc = getvar("dtdt_fluc")
        dqtdt_hadv = getvar("dqtdt_hadv")
        qt_nudge = get_nudgevar("qt_mean")
        dqtdt_fluc = getvar("dqtdt_fluc")
        subsidence = getvar("ls_subsidence")
        u_nudge = get_nudgevar("u_mean")
        v_nudge = get_nudgevar("v_mean")
        (; dTdt_hadv, H_nudge, dTdt_fluc, dqtdt_hadv, qt_nudge, dqtdt_fluc, subsidence, u_nudge, v_nudge)
    end
    for k in TC.real_center_indices(grid)
        aux_gm.dTdt_hadv[k] = nt.dTdt_hadv[k]
        aux_gm.H_nudge[k] = nt.H_nudge[k]
        aux_gm.dTdt_fluc[k] = nt.dTdt_fluc[k]
        aux_gm.dqtdt_hadv[k] = nt.dqtdt_hadv[k]
        aux_gm.qt_nudge[k] = nt.qt_nudge[k]
        aux_gm.dqtdt_fluc[k] = nt.dqtdt_fluc[k]
        aux_gm.subsidence[k] = nt.subsidence[k]
        aux_gm.u_nudge[k] = nt.u_nudge[k]
        aux_gm.v_nudge[k] = nt.v_nudge[k]
    end
end




# For surface params they create the forcing funcs in a call to surface_params which gets stored in a dedicated way (like a surf_params = FixedSurfaceFlux() object) defined in types.jl...
# we would have to hard code in our own storage mechanism then and figure out in main and everywhere else how to pass it around...
# for ARM they just have a function externally stored in APL, maybe that's an easier road... we dont have anything in 
# function create_forcing_funcs(::ForcingBase{ForcingSOCRATES_RF09_obs}, grid, state, param_set; Dat::SOCRATES_RF09_obsData) # THEY HAVE SOMETHING TO DO THIS IN TYPES.JL, also where float_or_func that i used to make TCf in dycore.jl came from... https://github.com/CliMA/TurbulenceConvection.jl/blob/7b4666baca418b00bb60e929de96fcc06100de57/src/types.jl
#     """
#     Returns, for each forcing, vectors in height F[k] where F[k] is a forcing functions f(t)... 
#     ARM does f(t,z) but we already know the output grid so this seems better for now?

#     # returns  named tuple (; dTdt_hadv(t), H_nudge(t), dqtdt_hadv(t), qt_nudge(t), subsidence(t), u_nudge(t), v_nudge(t)) 
#     """

#     thermo_params = TCP.thermodynamics_params(param_set)
#     R_d  = TCP.R_d(param_set) # TD.Parameters.R_d(thermo_params)
#     grav = TCP.grav(param_set)  # TD.Parameters.grav(thermo_params)
#     molmass_ratio = TCP.molmass_ratio(param_set)
#     FT = eltype(param_set)
#     aux_gm = TC.center_aux_grid_mean(state)

#     data_filename_obs  = Dat.data_filename #namelist["meta"]["datafile"]
#     data_filename_ERA5 = replace(data_filename_obs, "obs"=>"ERA5") # still needed to force winds antypd such (time-varying forcing)



#     function interp_along_dim(var, interp_dim, interp_dim_in; interp_dim_out=nothing,data=nothing, base=nothing, data_func = nothing, varg=nothing, interp_dim_in_is_full_array=true, reshape_ground=true, verbose=false)
#         " interpolation data 
#         - if we set interp_dim_out, we want to actually evaluate teh function along interp_dim at those locations
#         - otherwise, we want open ended fcns along that dimension

#         - data_func is a func to be applied to the raw data before it is processed, though perhaps it it most useful if interp_dim_out is unset and we are returning functions...
#         - vectorize_in means your input is an array and you need to loop over it too (as opposed to just being a fixed template vector)

#         To Do : decide types
#         "

#         vardata = isa(var,String) ? data[base][var]  : var # should return an ncdataset still

#         if isa(interp_dim, Number)
#             interp_dim_num = interp_dim
#         elseif isa(interp_dim,String)
#             if isa(vardata, NC.CFVariable)
#                 dimnames = NC.dimnames(vardata)
#                 interp_dim_num = findfirst(x->x==interp_dim,dimnames)
#             elseif isa(vardata, AbstractArray)
#                 error("can only specify numeric interp_dim for numeric vardata")                
#             end
#         end        

#         if !isnothing(varg)
#             vardatag   = isa(varg,String) ? data[base][varg]  : varg
#             if reshape_ground
#                 if isa(vardatag, Number) # allow for you to pass a single number and still concat
#                     sz_vardatag = collect(size(vardata)) # array
#                     sz_vardatag[interp_dim_num] = 1
#                     vardatag = fill(vardatag, sz_vardatag... ) # create array full with just this one value
#                 elseif isa(vardatag, NC.CFVariable) # vardatag doesn't have lev as a dimension so we need to add it in (and I guess check the others are in the same order)
#                     dimnamesg = NC.dimnames(vardatag)
#                     # in same order as vardata just in case
#                     vardatag = reshape(vardatag, (size(vardatag)...,1)) # add trailing singleton for lev
#                     dimnamesg = [dimnamesg..., "lev"]
#                     vardatag = permutedims(vardatag, [findfirst(x->x==dim,dimnames) for dim in dimnamesg] ) # permute into the right order
#                 elseif isa(vardatag, AbstractArray) # already processed in -- i think we can't guarantee then that the lev axis exists you may need to do that yourself
#                     # assume we need to add a new dimension at the same location as in the full array and order otherwise is preserved (quick check looks ok)
#                     sz_vardatag = collect(size(vardata)) # array
#                     insert!(sz_vardatag,interp_dim_num,1) # insert sz 1 at this location 
#                     vardatag =  reshape(vardatag, sz_vardatag...) # reshape
#                 end
#             end
#             vardata = cat(vardata,vardatag;dims=interp_dim_num)
#         end
        
#         if !isnothing(data_func) # apply data_func if we need to
#             vardata  = data_func(vardata)
#         end

#         # mapslices to apply along timedim, see https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices
#         if !interp_dim_in_is_full_array
#             if isnothing(interp_dim_out)
#                 return mapslices(d -> dd -> TC.pyinterp(dd, interp_dim_in, d)      , vardata, dims=[interp_dim_num,]) # wll return a lambda fcn that can be evaluated along that dimensoin
#             else
#                 return mapslices(d -> TC.pyinterp(interp_dim_out, interp_dim_in, d), vardata, dims=[interp_dim_num,]) # lambda fcn will evaluate
#             end
#         else # vectorize over input dim values as well as data (no support for vectorize over output dim yet)
#             # stack on new catd dimension, then split apart inside the fcn call
#             catd = ndims(vardata)+1
#             _input = cat(interp_dim_in, vardata;dims=catd ) # although maybe the input is just a vector in which case this won't work..... we could just pass it in
#             if isnothing(interp_dim_out)
#                 return dropdims( mapslices(d -> dd -> TC.pyinterp(dd , d[:,1], d[:,2])      , _input, dims=[interp_dim_num,catd]); dims=catd) # wll return a lambda fcn that can be evaluated along that dimensoin
#             else
#                 return dropdims( mapslices(d -> TC.pyinterp(interp_dim_out, d[:,1], d[:,2]), _input, dims=[interp_dim_num,catd]); dims=catd) # lambda fcn will evaluate

#             end
#         end
#     end

#     nt = NC.Dataset(data_filename_obs, "r")  do obs_data
#         ERA5_data = NC.Dataset(data_filename_ERA5, "r") 
#         data = Dict("obs" => obs_data, "ERA5" => ERA5_data)
#         t = data["obs"]["tsec"][:]
#         # wish dims could be labeled like xarray (DimensionalData?), rn we have to be careful about dim order -- rn seems to load as (?,?,lev, time)

#         function add_dim_from_template(vardata, dimnum)
#             # both data need to be labeled
#             sz_vardata = collect(size(vardata)) # array
#             insert!(sz_vardata,dimnum,1) # insert sz 1 at this location 
#             return reshape(vardata, sz_vardata...) # reshape
#         end


#         # we need all of these to calculate z in both obs and ERA5 based (θ_liq_ice is technically needed for obs in obs based run but maybe simpler to copy this formulation later for ERA5 based)
#         p, pg, omega, Ptend, T, Tg, q, ts, tsg, ρ, ρg, θ_liq_ice, θ_liq_iceg, pvg, qtg, Tz, qz, pz, T_bar, p_frac, dz, z, z0, tsz, Tvz = (Dict() for _ = 1:25) # assign all at once but maybe this is poorly optimizeable (though it's not in the computational side so idk...) see https://discourse.julialang.org/t/initialize-multiple-variables-in-one-line/28916/12
#         for base = ["obs","ERA5"]
#             T[ base] = data[base]["T"] # Provided data is liquid water tempERA5ture, T_L = T - L q_c/c_p but we also initialize with no condensate
#             dimnames = NC.dimnames(T[base])
#             lev_dim_num = findfirst(x->x=="lev",dimnames)
#             ldn = lev_dim_num
#             Tg[base] = add_dim_from_template(data[base]["Tg", ], ldn)
#             q[ base] = data[base]["q"] # i think this is qt, though they've labeled it water vapor on the netcdf, matches the plots they give...
#             pg[base] = add_dim_from_template(data[base]["Ps"], ldn ) # add lev to pg, is already 3d dataset (lon,lat, time)
#             pg[base] = add_dim_from_template(data[base]["Ps"], ldn ) # add lev to pg, is already 3d dataset (lon,lat, time)
#             omega[base] = data[base]["omega"]
#             Ptend[base] = add_dim_from_template(data[base]["Ptend"], ldn )
#             p[ base] = data[base]["lev"] # is 1D 
#             L = length(p[base])
#             for dim in dimnames # these i believe go left to right so adding like this should work
#                 if dim != "lev"
#                     p[base] = add_dim_from_template(p[base],  findfirst(x->x==dim,dimnames)) # add all dims cause p default only has lev
#                 end
#             end
#             repeat_nums = Int.(size(T[base]) ./ size(p[base])) # for building p into a full larray from lev...
#             p[base] = repeat(p[base], outer = repeat_nums)
#             pvg[base]        = TD.saturation_vapor_pressure.(thermo_params, Tg[base], TD.Liquid())
#             qtg[base]        = (1 / molmass_ratio) .* pvg[base] ./ (pg[base] .- pvg[base]) #Total water mixing ratio at surface   
#             ts[base]         = TD.PhaseEquil_pTq.(thermo_params, p[base], T[base], q[base]) 
#             tsg[base]        = TD.PhaseEquil_pTq.(thermo_params, pg[base], Tg[base], qtg[base]) 
#             ρ[ base]         = TD.air_density.(thermo_params,ts[base])
#             ρg[ base]        = TD.air_density.(thermo_params,tsg[base])
#             θ_liq_ice[base]  = TD.liquid_ice_pottemp.(thermo_params, ts[base]) # 
#             θ_liq_iceg[base] = TD.liquid_ice_pottemp.(thermo_params, tsg[base]) # 

#             # calculate z (will need to do for both obs and ERA5) (could maybe replace some of this with mapslices and a function but...)
#             Tz[ base] = cat( T[base], Tg[base];dims=ldn)
#             qz[ base] = cat( q[base],qtg[base];dims=ldn)
#             pz[ base] = cat( p[base], pg[base];dims=ldn)
#             tsz[base] = cat(ts[base],tsg[base];dims=ldn)
#             Tvz[base] = TD.virtual_temperature.(thermo_params, tsz[base]) # virtual temp, we havent returned these for now...
#             Lz = L+1 # cause we extended it using the ground...
#             T_bar[base]   = Statistics.mean((selectdim(Tvz[base], ldn, 1:Lz-1), selectdim(Tvz[base], lev_dim_num, 2:Lz)))
#             p_frac[base]  =  selectdim( pz[base], lev_dim_num, 2:Lz) ./ selectdim( pz[base], ldn, 1:Lz-1)
#             dz[base] = @. (R_d * T_bar[base] / grav) * log(p_frac[base])
#             z[base]  = reverse(cumsum(reverse(dz[base];dims=ldn);dims=ldn);dims=ldn) #cumsum(dz[base]) # grid is already defined from  (from Grid.jl)
#             s_sz = collect(size(z[base])) # should now be same dims as T
#             s_sz[ldn] = 1 # we just want to add a single slice of zeros for the ground
#             _z0       = zeros(FT,s_sz...)
#             z0[base]  = cat(z[base], _z0; dims=ldn)
#         end

#         function drop_lat_lon(vardata;dims=tuple(collect(findfirst(x->x==dim,NC.dimnames(data["ERA5"]["T"])) for dim in ["lat","lon"] )...),bases=nothing)
#             # gotta do this at the end after we've made full use of our dimnames since our data is unlabeled
#             # in general we could cheat in the future since the socrates order seems to always be lon lat lev time regardless of which dims exist...
#             # i guess i also don't know what happens if we've turned the time dim into just a fcn -- in principle it's last so that shouldnt hurt

#             if isnothing(bases)
#                 vardata = dropdims(vardata, dims = dims)
#             else
#                 for base in bases
#                     vardata[base] = dropdims(vardata[base], dims = dims)
#                 end
#             end
#             return vardata
#         end

#         drop_lat_lon_bases(vardata) = drop_lat_lon(vardata; bases=["obs","ERA5"])

#         dimnames = NC.dimnames(data["ERA5"]["T"])
#         lev_dim_num = findfirst(x->x=="lev",dimnames)
#         ldn = lev_dim_num
#         tdn = findfirst(x->x=="time",dimnames)

#         # convert our output vars to evaluations on new z
#         # for the data version they may be labeled but also might not be if more calculations were done so just use the dimnums from above (none of them as is delete any dims so this is still safe)
#         getvar_z(var,base; kwargs...)                             = interp_along_dim(var    , "lev", reverse(z[ base];dims=ldn); interp_dim_out=vec(grid.zc.z), data=data, base=base, data_func=x->reverse(x;dims=ldn), kwargs...) # reverse data so it's with monotonic increasing z (maybe also pad w/ ground since clima grid might have a lower lower bound... current extrapolation i believe is just to repeat the ends)
#         getvardata_z(vardata,base; kwargs...)                     = interp_along_dim(vardata,  ldn , reverse(z[ base];dims=ldn); interp_dim_out=vec(grid.zc.z),                       data_func=x->reverse(x;dims=ldn), kwargs...) # use ground data to extend interpolation since i think the clima grid has lower lower bound (default was just to repeat value but maybe could also later code pyinterp here to have an extend keyword instead)
#         getvar_z_wg(var, varg, base; kwargs...)                   = interp_along_dim(var    , "lev", reverse(z0[base];dims=ldn); interp_dim_out=vec(grid.zc.z), data=data, base=base, data_func=x->reverse(x;dims=ldn), varg=varg, kwargs...) # reverse data so it's with monotonic increasing z (maybe also pad w/ ground since clima grid might have a lower lower bound... current extrapolation i believe is just to repeat the ends)
#         getvardata_z_wg(vardata, vardatag, base; kwargs...)       = interp_along_dim(vardata,  ldn , reverse(z0[base];dims=ldn); interp_dim_out=vec(grid.zc.z),                       data_func=x->reverse(x;dims=ldn), varg=vardatag, kwargs...) # use ground data to extend interpolation since i think the clima grid has lower lower bound (default was just to repeat value but maybe could also later code pyinterp here to have an extend keyword instead)
       
#         # convert our output vars to functions of t for later evalutation
#         getvar_t(var,base; kwargs...)                             = interp_along_dim(var    , "time", t; data=data, base=base, interp_dim_in_is_full_array=false, kwargs...)
#         getvardata_t(vardata; kwargs...)                          = interp_along_dim(vardata,    tdn, t;                       interp_dim_in_is_full_array=false, kwargs...)
#         # combine
#         getvar_z_tfunc(var,base; kwargs...)                       = getvardata_t(getvar_z(var,base; kwargs...)) 
#         getvardata_z_tfunc(vardata,base; kwargs...)               = getvardata_t(getvardata_z(vardata, base; kwargs...)) 
#         getvar_z_tfunc_wg(var, varg, base; kwargs...)             = getvardata_t(getvar_z_wg(    var    , varg,base; data_func=x->reverse(x;dims=ldn), kwargs...)) 
#         getvardata_z_tfunc_wg(vardata, vardatag, base; kwargs...) = getvardata_t(getvardata_z_wg(vardata, vardatag, base; data_func=x->reverse(x;dims=ldn), kwargs...)) 

#         # we have safety w/ _wg methods to extrapolate towards ground below lowest model level to match prescribed grid... but what are we setting as surface?

#         # Paper explicitly mentions nudging T_L (essentially T, q?), u, (Tg but maybe only for GCM runs?)
#         dTdt_hadv  = drop_lat_lon(getvar_z_tfunc_wg("divT", FT(0), "ERA5"))            # getvar("dtdt_hadv")       <-- references from the LES setup -- not sure what exactly is needed/not needed... some of these we don't have for example... and I just made them 
#         dqtdt_hadv = drop_lat_lon(getvar_z_tfunc_wg("divq", FT(0), "ERA5"))            # getvar("dqtdt_hadv")
        
#         ρ_sub = ρ["ERA5"]
#         # pad wit ground at bottom...
#         sz_ρ_sub = collect(size(ρ_sub)) # array
#         ρ_subg = ρg["ERA5"] # fill(0, sz_ρ_sub... ) # create array full with just this one value ( this one can't be zero so....)
#         ρ_sub = cat(ρ_sub,ρ_subg;dims=ldn)

#         dpdt_g    = Ptend["ERA5"] # should be lon lat lev tkme (gotta expand it to match size of omega + ground... or will broadcasting work...)
#         psub      = cat( p["ERA5"], pg["ERA5"]; dims=ldn )
#         omega_sub = cat( omega["ERA5"], dpdt_g; dims=ldn )

#         c   = 100
#         a   = -2 + exp(250/c)
#         f_p = p -> @. 2(a+1) / (a+exp(p/c)) - 1
#         # stack p and pg
#         f_p = f_p(psub)

#         subsidence = drop_lat_lon(getvar_z_tfunc_wg("omega", "Ptend", "ERA5", data_func = x -> -(x .- omega_sub.*f_p )  ./ (ρ_sub .* grav)))  # getvar("ls_subsidence")   # I believe we want subsidence = -omega/(rho*g) as per https://github.com/CliMA/pycles/blob/6be800b394a8c61e84d1141277cf825569c40c2e/ForcingGCMFixed.pyx#L265

#         # Paper says Both LES experiments are forced by ERA5-derived geostrophic winds and nudged toward the ERA5 horizontal winds with a nudging timescale of 1 hr for the ERA5-based simulation and 20 min for the Obs-based simulation.
#         # so not sure what to do, we don't have a separate forcing and nudging of the mean  (and for the forcing it's one column anyway now) so I just went with the real one not the geostrophic...
#         # in the online files the obs forcings are constant in time while the ERA5 are bit
#         u_nudge    = drop_lat_lon(getvar_z_tfunc_wg("u", FT(0), "ERA5"))               # getvar("u_mean")
#         v_nudge    = drop_lat_lon(getvar_z_tfunc_wg("v", FT(0), "ERA5"))               # getvar("v_mean")
#         # u_force    = drop_lat_lon(getvar_z_tfunc_wg("ug", FT(0), "ERA5"))               # getvar("u_mean") # turned off for now cause confused on using
#         # v_force    = drop_lat_lon(getvar_z_tfunc_wg("vg", FT(0), "ERA5"))               # getvar("v_mean")

#         # since we don't know the parititioning (we only have T_L, q_t, p, we start with no condensate -- but given theta_li conservation this should converge to the right T over time right... (given consistent forcing...))
#         # we could solve T given T = T_L + (L/c_p)(q_t - q_v^*(T)) but we dont have density, so since the obs profiles for q,T are fixed hopefully we can just achieve convergence...

#         # we already padded pg for calculating p earlier for calculating z... so thθ_liq_iceg, qtg doesn't need creating of a lev dim... reshape_ground -> false then
#         H_nudge    = drop_lat_lon(getvardata_z_tfunc_wg(θ_liq_ice["obs"], θ_liq_iceg["obs"], "obs" ; reshape_ground=false))   # getvar("thetali_mean") # this is basically our T_L 
#         qt_nudge   = drop_lat_lon(getvardata_z_tfunc_wg(q["obs"]        ,        qtg["obs"], "obs" ; reshape_ground=false, verbose=false)) # getvar("qt_mean")      # tbis is basically our q              do we want this? should would we relax global q? 

#         # we also have access to surface pressure tendency, not sure how we should use
#         # dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge = drop_lat_lon.(  (dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge) ) # drop lat/lon here now tht everything is done... no bases here... ( we wont do initally but in each calc is fine to do so dont need this)
#         close(ERA5_data) # make sure we close this file since didn't use do block

#         # geostrophic velocity profiles # adapted from dycoms/

#         prog_gm = TC.center_prog_grid_mean(state) # do this way cause the aux_gm vars are fixed to only have this one for geostrophic winds... format copied form dycore.jl
#         ρ_c = prog_gm.ρ
#         lg = CC.Fields.local_geometry_field(axes(ρ_c))
        
#         # uₕ_g = ((u,v,l) -> t -> [CCG.Covariant12Vector(CCG.UVVector(u(t)[1], v(t)[1]), l)]).(u_force,v_force, tuple([lg[k] for k in TC.real_center_indices(grid)]...)) # can you build a point on a function? idk... might need to put in manual evaluation later for this to keep types safe...
#         # ^ this line is very slow for some reason ^, the iteration at the end seemed to work cause braodcasting wouldnt work without it..., not sure why not... (maybe tuple(lg) works? idk...)

#         (; dTdt_hadv, H_nudge, dqtdt_hadv, qt_nudge, subsidence, u_nudge, v_nudge) # removed geostrophic cause it's slow and idk how to use geostrophic forcing with relax to regular... seems it is one or the other no?
#     end
#     # update_forcing() calls this and is in Cases.jl
# end


function initialize(forcing::ForcingBase{T}, grid, state, param_set) where {T<:ForcingSOCRATES}# added param_set so we can calculate stuff....
    FT = eltype(param_set)
    aux_gm = TC.center_aux_grid_mean(state)


    ug_keys = (:ug_nudge, :vg_nudge)
    forcing_funcs = forcing.forcing_funcs[]
    @show(keys(forcing_funcs))

    # set the geostrophic velocity profiles -- need to check if we actually have a fcn that can take in t=0 and return a profile... ( i think ours is one func for each z so whoops... maybe need to reconstruct?)
    g_func = (f) -> f([FT(0)])[1]
    prof_ug = g_func.(forcing_funcs[:ug_nudge]) # map over the forcing funcs to get the profile at t=0
    prof_vg = g_func.(forcing_funcs[:vg_nudge])
    # it wants a fcn out, could edit src/Fields.jl I guess to add another method but maybe it needs to face/center points idk...
    prof_ug = Dierckx.Spline1D(vec(grid.zc.z), vec(prof_ug); k = 1)
    prof_vg = Dierckx.Spline1D(vec(grid.zc.z), vec(prof_vg); k = 1)

    # prof_ug = map.(g_func, collect(forcing_funcs[:ug_nudge])) # map over the forcing funcs to get the profile at t=0
    # prof_ug = map.(g_func, collect(forcing_funcs[:vg_nudge]))

    # prof_ug = forcing_funcs[:ug_nudge]([FT(0)])[1] # full profile;  add this (also needs its own forcing update function... cause we gotta do the set_z thing...)
    # prof_vg = forcing_funcs[:vg_nudge]([FT(0)])[1] # 
    aux_gm_uₕ_g = TC.grid_mean_uₕ_g(state)
    TC.set_z!(aux_gm_uₕ_g, prof_ug, prof_vg) 
    



    forcing_funcs = forcing_funcs[setdiff(keys(forcing_funcs), ug_keys)] # keys we don't need.
    for (name,funcs) in zip(keys(forcing_funcs),forcing_funcs) # iterate over the named tuple of our forcings...
        for k in TC.real_center_indices(grid)
            func = funcs[k]
            getproperty(aux_gm, name)[k] = func([FT(0)])[1] # apply to time = 0 and apply to aux_gm, turn to vec cause needs to be cast as in https://github.com/CliMA/TurbulenceConvection.jl/blob/a9ebce1f5f15f049fc3719a013ddbc4a9662943a/src/utility_functions.jl#L48
        end
    end
    return nothing
end
