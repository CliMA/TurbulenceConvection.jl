update_radiation(self::RadiationBase, state, t::Real, param_set) = nothing
initialize(self::RadiationBase{RadiationNone}, state) = nothing

"""
see eq. 3 in Stevens et. al. 2005 DYCOMS paper
"""
function update_radiation(self::RadiationBase{RadiationDYCOMS_RF01}, state, t::Real, param_set)
    grid = TC.Grid(state)
    cp_d = TCP.cp_d(param_set)
    aux_gm = TC.center_aux_grid_mean(state)
    aux_gm_f = TC.face_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    q_tot_f = TC.face_aux_turbconv(state).temporary_1
    ρ_f = aux_gm_f.ρ
    ρ_c = prog_gm.ρ
    # find zi (level of 8.0 g/kg isoline of qt)
    # TODO: report bug: zi and ρ_i are not initialized
    FT = TC.float_type(state)
    zi = FT(0)
    ρ_i = FT(0)
    kc_surf = TC.kc_surface(grid)
    q_tot_surf = aux_gm.q_tot[kc_surf]
    If = CCO.InterpolateC2F(; bottom = CCO.SetValue(q_tot_surf), top = CCO.Extrapolate())
    @. q_tot_f .= If(aux_gm.q_tot)
    @inbounds for k in TC.real_face_indices(grid)
        if (q_tot_f[k] < 8.0 / 1000)
            idx_zi = k
            # will be used at cell faces
            zi = grid.zf[k].z
            ρ_i = ρ_f[k]
            break
        end
    end

    ρ_z = Dierckx.Spline1D(vec(grid.zc.z), vec(ρ_c); k = 1)
    q_liq_z = Dierckx.Spline1D(vec(grid.zc.z), vec(aux_gm.q_liq); k = 1)

    integrand(ρq_l, params, z) = params.κ * ρ_z(z) * q_liq_z(z)
    rintegrand(ρq_l, params, z) = -integrand(ρq_l, params, z)

    z_span = (grid.zmin, grid.zmax)
    rz_span = (grid.zmax, grid.zmin)
    params = (; κ = self.kappa)

    Δz = TC.get_Δz(prog_gm.ρ)[1]
    rprob = ODE.ODEProblem(rintegrand, 0.0, rz_span, params; dt = Δz)
    rsol = ODE.solve(rprob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_0 = rsol.(vec(grid.zf.z))

    prob = ODE.ODEProblem(integrand, 0.0, z_span, params; dt = Δz)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    q_1 = sol.(vec(grid.zf.z))
    parent(aux_gm_f.f_rad) .= self.F0 .* exp.(-q_0)
    parent(aux_gm_f.f_rad) .+= self.F1 .* exp.(-q_1)

    # cooling in free troposphere
    @inbounds for k in TC.real_face_indices(grid)
        if grid.zf[k].z > zi
            cbrt_z = cbrt(grid.zf[k].z - zi)
            aux_gm_f.f_rad[k] += ρ_i * cp_d * self.divergence * self.alpha_z * (cbrt_z^4 / 4 + zi * cbrt_z)
        end
    end

    ∇c = CCO.DivergenceF2C()
    wvec = CC.Geometry.WVector
    @. aux_gm.dTdt_rad = -∇c(wvec(aux_gm_f.f_rad)) / ρ_c / cp_d

    return
end

function initialize(self::RadiationBase{RadiationLES}, state, LESDat::LESData)
    grid = TC.Grid(state)
    # load from LES
    aux_gm = TC.center_aux_grid_mean(state)
    dTdt = NC.Dataset(LESDat.les_filename, "r") do data
        imin = LESDat.imin
        imax = LESDat.imax

        # interpolate here
        zc_les = Array(TC.get_nc_data(data, "zc"))
        meandata = TC.mean_nc_data(data, "profiles", "dtdt_rad", imin, imax)
        pyinterp(vec(grid.zc.z), zc_les, meandata)
    end
    @inbounds for k in TC.real_center_indices(grid)
        aux_gm.dTdt_rad[k] = dTdt[k]
    end
    return
end


function initialize(self::RadiationBase{RadiationTRMM_LBA}, state)
    grid = TC.Grid(state)
    aux_gm = TC.center_aux_grid_mean(state)
    rad = APL.TRMM_LBA_radiation(eltype(grid))
    @inbounds for k in TC.real_center_indices(grid)
        aux_gm.dTdt_rad[k] = rad(0, grid.zc[k].z)
    end
    return nothing
end

function update_radiation(self::RadiationBase{RadiationTRMM_LBA}, state, t::Real, param_set)
    grid = TC.Grid(state)
    aux_gm = TC.center_aux_grid_mean(state)
    rad = APL.TRMM_LBA_radiation(eltype(grid))
    @inbounds for k in TC.real_center_indices(grid)
        aux_gm.dTdt_rad[k] = rad(t, grid.zc[k].z)
    end
    return nothing
end


# technically the code as written already sets this cause it's in the SSCF output...
function initialize(radiation::RadiationBase{RadiationSOCRATES}, state) #where {T <: RadiationSOCRATES}
    FT = TC.float_type(state)
    grid = TC.Grid(state)
    aux_gm = TC.center_aux_grid_mean(state)
    radiation_funcs = radiation.radiation_funcs[]

    @inbounds for k in TC.real_center_indices(grid)
        aux_gm.dTdt_rad[k] = radiation_funcs[:dTdt_rad][k]([FT(0)])[1] # apply to time = 0 and apply to aux_gm, turn to vec cause needs to be cast as in https://github.com/CliMA/TurbulenceConvection.jl/blob/a9ebce1f5f15f049fc3719a013ddbc4a9662943a/src/utility_functions.jl#L48
    end
    return nothing
end

function update_radiation(radiation::RadiationBase{RadiationSOCRATES}, state, t::Real, param_set) # where {T <: RadiationSOCRATES}
    grid = TC.Grid(state)
    FT = TC.float_type(state)
    aux_gm = TC.center_aux_grid_mean(state)
    radiation_funcs = radiation.radiation_funcs[]

    @inbounds for k in TC.real_center_indices(grid)
        aux_gm.dTdt_rad[k] = radiation_funcs[:dTdt_rad][k]([FT(t)])[1] # apply to time = 0 and apply to aux_gm, turn to vec cause needs to be cast as in
    end
    return nothing
end
    