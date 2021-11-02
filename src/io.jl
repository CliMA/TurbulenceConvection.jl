#= This file is being deprecated in favor of using diagnostics.jl =#

initialize_io(self::ForcingBase, Stats) = nothing

initialize_io(self::RadiationBase, Stats::NetCDFIO_Stats) = nothing
io(self::RadiationBase, grid, state, Stats::NetCDFIO_Stats) = nothing

function initialize_io(self::RadiationBase{RadiationLES}, Stats::NetCDFIO_Stats)
    add_profile(Stats, "rad_dTdt")
    add_profile(Stats, "rad_flux")
    return
end

function io(self::RadiationBase{RadiationLES}, Stats::NetCDFIO_Stats)
    write_profile(Stats, "rad_dTdt", self.dTdt)
    write_profile(Stats, "rad_flux", self.f_rad)
    return
end

function initialize_io(self::RadiationBase{RadiationDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    add_profile(Stats, "rad_dTdt")
    add_profile(Stats, "rad_flux")
    return
end

function io(self::RadiationBase{RadiationDYCOMS_RF01}, Stats::NetCDFIO_Stats)
    write_profile(Stats, "rad_dTdt", self.dTdt)
    write_profile(Stats, "rad_flux", self.f_rad)
    return
end

function initialize_io(en::EnvironmentVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "env_cloud_base")
    add_ts(Stats, "env_cloud_top")
    add_ts(Stats, "env_cloud_cover")
    add_ts(Stats, "env_lwp")
    add_ts(Stats, "env_iwp")
    return
end

function io(en::EnvironmentVariables, grid, state, Stats::NetCDFIO_Stats)
    # Assuming amximum overlap in environmental clouds
    write_ts(Stats, "env_cloud_cover", en.cloud_cover)
    write_ts(Stats, "env_cloud_base", en.cloud_base)
    write_ts(Stats, "env_cloud_top", en.cloud_top)
    write_ts(Stats, "env_lwp", en.lwp)
    write_ts(Stats, "env_iwp", en.iwp)
    return
end

function initialize_io(precip::PrecipVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "rwp_mean")
    add_ts(Stats, "cutoff_precipitation_rate")
    return
end

function io(precip::PrecipVariables, grid, state, Stats::NetCDFIO_Stats)
    write_ts(Stats, "rwp_mean", precip.mean_rwp)

    #TODO - change to rain rate that depends on rain model choice
    write_ts(Stats, "cutoff_precipitation_rate", precip.cutoff_precipitation_rate)
    return
end

function initialize_io(up::UpdraftVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "updraft_cloud_fraction")

    add_ts(Stats, "updraft_cloud_cover")
    add_ts(Stats, "updraft_cloud_base")
    add_ts(Stats, "updraft_cloud_top")
    add_ts(Stats, "updraft_lwp")
    add_ts(Stats, "updraft_iwp")
    return
end

function io(up::UpdraftVariables, grid, state, Stats::NetCDFIO_Stats)
    write_profile(Stats, "updraft_cloud_fraction", up.cloud_fraction)
    # Note definition of cloud cover : each updraft is associated with a cloud cover equal to the maximum
    # area fraction of the updraft where ql > 0. Each updraft is assumed to have maximum overlap with respect to
    # itup (i.e. no consideration of tilting due to shear) while the updraft classes are assumed to have no overlap
    # at all. Thus total updraft cover is the sum of each updraft"s cover
    write_ts(Stats, "updraft_cloud_cover", sum(up.cloud_cover))
    write_ts(Stats, "updraft_cloud_base", minimum(abs.(up.cloud_base)))
    write_ts(Stats, "updraft_cloud_top", maximum(abs.(up.cloud_top)))
    write_ts(Stats, "updraft_lwp", up.lwp)
    write_ts(Stats, "updraft_iwp", up.iwp)
    return
end

function initialize_io(gm::GridMeanVariables, Stats::NetCDFIO_Stats)
    add_ts(Stats, "lwp_mean")
    add_ts(Stats, "iwp_mean")
    add_ts(Stats, "cloud_base_mean")
    add_ts(Stats, "cloud_top_mean")
    add_ts(Stats, "cloud_cover_mean")
    return
end

function io(gm::GridMeanVariables, grid, state, Stats::NetCDFIO_Stats)
    write_ts(Stats, "cloud_cover_mean", gm.cloud_cover)

    write_ts(Stats, "lwp_mean", gm.lwp)
    write_ts(Stats, "iwp_mean", gm.iwp)
    write_ts(Stats, "cloud_base_mean", gm.cloud_base)
    write_ts(Stats, "cloud_top_mean", gm.cloud_top)
    return
end

# Initialize the IO pertaining to this class
function initialize_io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)

    initialize_io(edmf.UpdVar, Stats)
    initialize_io(edmf.EnvVar, Stats)
    initialize_io(edmf.Precip, Stats)

    add_profile(Stats, "entrainment_sc")
    add_profile(Stats, "detrainment_sc")
    add_profile(Stats, "nh_pressure")
    add_profile(Stats, "nh_pressure_adv")
    add_profile(Stats, "nh_pressure_drag")
    add_profile(Stats, "nh_pressure_b")
    add_profile(Stats, "asp_ratio")
    add_profile(Stats, "horiz_K_eddy")
    add_profile(Stats, "sorting_function")
    add_profile(Stats, "b_mix")
    add_ts(Stats, "rd")
    add_profile(Stats, "turbulent_entrainment")
    add_profile(Stats, "turbulent_entrainment_full")
    add_profile(Stats, "turbulent_entrainment_W")
    add_profile(Stats, "turbulent_entrainment_H")
    add_profile(Stats, "turbulent_entrainment_QT")
    add_profile(Stats, "massflux")
    add_profile(Stats, "massflux_h")
    add_profile(Stats, "massflux_qt")
    add_profile(Stats, "massflux_tendency_h")
    add_profile(Stats, "massflux_tendency_qt")
    add_profile(Stats, "diffusive_flux_h")
    add_profile(Stats, "diffusive_flux_u")
    add_profile(Stats, "diffusive_flux_v")
    add_profile(Stats, "diffusive_flux_qt")
    add_profile(Stats, "diffusive_tendency_h")
    add_profile(Stats, "diffusive_tendency_qt")
    add_profile(Stats, "total_flux_h")
    add_profile(Stats, "total_flux_qt")
    add_profile(Stats, "updraft_qt_precip")
    add_profile(Stats, "updraft_thetal_precip")
    # Diff mixing lengths: Ignacio
    add_profile(Stats, "ed_length_scheme")
    add_profile(Stats, "mixing_length_ratio")
    add_profile(Stats, "entdet_balance_length")
    return
end

function io(edmf::EDMF_PrognosticTKE, grid, state, Stats::NetCDFIO_Stats, TS::TimeStepping, param_set)

    mean_nh_pressure = face_field(grid)
    mean_nh_pressure_adv = face_field(grid)
    mean_nh_pressure_drag = face_field(grid)
    mean_nh_pressure_b = face_field(grid)

    mean_asp_ratio = center_field(grid)
    mean_entr_sc = center_field(grid)
    mean_detr_sc = center_field(grid)
    massflux = center_field(grid)
    mean_frac_turb_entr = center_field(grid)
    mean_horiz_K_eddy = center_field(grid)
    mean_sorting_function = center_field(grid)
    mean_b_mix = center_field(grid)

    io(edmf.UpdVar, grid, state, Stats)
    io(edmf.EnvVar, grid, state, Stats)
    io(edmf.Precip, grid, state, Stats)

    aux_up = center_aux_updrafts(state)
    a_up_bulk = center_aux_turbconv(state).bulk.area

    write_ts(Stats, "rd", StatsBase.mean(edmf.pressure_plume_spacing))

    @inbounds for k in real_center_indices(grid)
        if a_up_bulk[k] > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                massflux[k] += interpf2c(edmf.m, grid, k, i)
                mean_entr_sc[k] += aux_up[i].area[k] * edmf.entr_sc[i, k] / a_up_bulk[k]
                mean_detr_sc[k] += aux_up[i].area[k] * edmf.detr_sc[i, k] / a_up_bulk[k]
                mean_asp_ratio[k] += aux_up[i].area[k] * edmf.asp_ratio[i, k] / a_up_bulk[k]
                mean_frac_turb_entr[k] += aux_up[i].area[k] * edmf.frac_turb_entr[i, k] / a_up_bulk[k]
                mean_horiz_K_eddy[k] += aux_up[i].area[k] * edmf.horiz_K_eddy[i, k] / a_up_bulk[k]
                mean_sorting_function[k] += aux_up[i].area[k] * edmf.sorting_function[i, k] / a_up_bulk[k]
                mean_b_mix[k] += aux_up[i].area[k] * edmf.b_mix[i, k] / a_up_bulk[k]
            end
        end
    end

    @inbounds for k in real_face_indices(grid)
        a_up_bulk_f =
            interpc2f(a_up_bulk, grid, k; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        if a_up_bulk_f > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                a_up_f = interpc2f(
                    aux_up[i].area,
                    grid,
                    k;
                    bottom = SetValue(edmf.area_surface_bc[i]),
                    top = SetZeroGradient(),
                )
                mean_nh_pressure[k] += a_up_f * edmf.nh_pressure[i, k] / a_up_bulk_f
                mean_nh_pressure_b[k] += a_up_f * edmf.nh_pressure_b[i, k] / a_up_bulk_f
                mean_nh_pressure_adv[k] += a_up_f * edmf.nh_pressure_adv[i, k] / a_up_bulk_f
                mean_nh_pressure_drag[k] += a_up_f * edmf.nh_pressure_drag[i, k] / a_up_bulk_f
            end
        end
    end

    write_profile(Stats, "turbulent_entrainment", mean_frac_turb_entr)
    write_profile(Stats, "horiz_K_eddy", mean_horiz_K_eddy)
    write_profile(Stats, "entrainment_sc", mean_entr_sc)
    write_profile(Stats, "detrainment_sc", mean_detr_sc)
    write_profile(Stats, "nh_pressure", mean_nh_pressure)
    write_profile(Stats, "nh_pressure_adv", mean_nh_pressure_adv)
    write_profile(Stats, "nh_pressure_drag", mean_nh_pressure_drag)
    write_profile(Stats, "nh_pressure_b", mean_nh_pressure_b)
    write_profile(Stats, "asp_ratio", mean_asp_ratio)
    write_profile(Stats, "massflux", massflux)
    write_profile(Stats, "massflux_h", edmf.massflux_h)
    write_profile(Stats, "massflux_qt", edmf.massflux_qt)
    write_profile(Stats, "massflux_tendency_h", edmf.massflux_tendency_h)
    write_profile(Stats, "massflux_tendency_qt", edmf.massflux_tendency_qt)
    write_profile(Stats, "diffusive_flux_h", edmf.diffusive_flux_h)
    write_profile(Stats, "diffusive_flux_qt", edmf.diffusive_flux_qt)
    write_profile(Stats, "diffusive_flux_u", edmf.diffusive_flux_u)
    write_profile(Stats, "diffusive_flux_v", edmf.diffusive_flux_v)
    write_profile(Stats, "diffusive_tendency_h", edmf.diffusive_tendency_h)
    write_profile(Stats, "diffusive_tendency_qt", edmf.diffusive_tendency_qt)
    write_profile(Stats, "total_flux_h", edmf.massflux_h .+ edmf.diffusive_flux_h)
    write_profile(Stats, "total_flux_qt", edmf.massflux_h .+ edmf.diffusive_flux_qt)
    write_profile(Stats, "updraft_qt_precip", edmf.UpdThermo.qt_tendency_rain_formation_tot)
    write_profile(Stats, "updraft_thetal_precip", edmf.UpdThermo.θ_liq_ice_tendency_rain_formation_tot)

    #Different mixing lengths : Ignacio
    write_profile(Stats, "ed_length_scheme", edmf.mls)
    write_profile(Stats, "mixing_length_ratio", edmf.ml_ratio)
    write_profile(Stats, "entdet_balance_length", edmf.l_entdet)
    return
end