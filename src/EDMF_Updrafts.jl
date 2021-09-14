function initialize(tptke, self::UpdraftVariables, GMV::GridMeanVariables)
    kc_surf = kc_surface(self.grid)

    self.W.values .= 0
    self.B.values .= 0
    @inbounds for i in xrange(self.n_updrafts)
        @inbounds for k in real_center_indices(self.grid)
            # Simple treatment for now, revise when multiple updraft closures
            # become more well defined
            if self.prognostic
                self.Area.values[i, k] = 0.0 #self.updraft_fraction/self.n_updrafts
            else
                self.Area.values[i, k] = self.updraft_fraction / self.n_updrafts
            end
            self.QT.values[i, k] = GMV.QT.values[k]
            self.QL.values[i, k] = GMV.QL.values[k]
            self.H.values[i, k] = GMV.H.values[k]
            self.T.values[i, k] = GMV.T.values[k]
        end

        self.Area.values[i, kc_surf] = self.updraft_fraction / self.n_updrafts
    end

    set_means(tptke, self, GMV)

    return
end

function initialize_DryBubble(tptke, self::UpdraftVariables, GMV::GridMeanVariables, ref_state::ReferenceState)
    dz = self.grid.Δz

    # criterion 2: b>1e-4
    #! format: off
    z_in = [
          75.,  125.,  175.,  225.,  275.,  325.,  375.,  425.,  475.,
         525.,  575.,  625.,  675.,  725.,  775.,  825.,  875.,  925.,
         975., 1025., 1075., 1125., 1175., 1225., 1275., 1325., 1375.,
        1425., 1475., 1525., 1575., 1625., 1675., 1725., 1775., 1825.,
        1875., 1925., 1975., 2025., 2075., 2125., 2175., 2225., 2275.,
        2325., 2375., 2425., 2475., 2525., 2575., 2625., 2675., 2725.,
        2775., 2825., 2875., 2925., 2975., 3025., 3075., 3125., 3175.,
        3225., 3275., 3325., 3375., 3425., 3475., 3525., 3575., 3625.,
        3675., 3725., 3775., 3825., 3875., 3925.]

    thetal_in = [
        299.9882, 299.996 , 300.0063, 300.0205, 300.04  , 300.0594,
        300.0848, 300.1131, 300.1438, 300.1766, 300.2198, 300.2567,
        300.2946, 300.3452, 300.3849, 300.4245, 300.4791, 300.5182,
        300.574 , 300.6305, 300.6668, 300.7222, 300.7771, 300.8074,
        300.8591, 300.9092, 300.9574, 300.9758, 301.0182, 301.0579,
        301.0944, 301.1276, 301.1572, 301.1515, 301.1729, 301.1902,
        301.2033, 301.2122, 301.2167, 301.2169, 301.2127, 301.2041,
        301.1913, 301.1743, 301.1533, 301.1593, 301.1299, 301.097 ,
        301.0606, 301.0212, 300.9788, 300.9607, 300.9125, 300.8625,
        300.8108, 300.7806, 300.7256, 300.6701, 300.6338, 300.5772,
        300.5212, 300.482 , 300.4272, 300.3875, 300.3354, 300.2968,
        300.2587, 300.2216, 300.1782, 300.1452, 300.1143, 300.0859,
        300.0603, 300.0408, 300.0211, 300.0067, 299.9963, 299.9884    ]

    Area_in = [
        0.04 , 0.055, 0.07 , 0.08 , 0.085, 0.095, 0.1  , 0.105, 0.11 ,
        0.115, 0.115, 0.12 , 0.125, 0.125, 0.13 , 0.135, 0.135, 0.14 ,
        0.14 , 0.14 , 0.145, 0.145, 0.145, 0.15 , 0.15 , 0.15 , 0.15 ,
        0.155, 0.155, 0.155, 0.155, 0.155, 0.155, 0.16 , 0.16 , 0.16 ,
        0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 ,
        0.155, 0.155, 0.155, 0.155, 0.155, 0.155, 0.15 , 0.15 , 0.15 ,
        0.15 , 0.145, 0.145, 0.145, 0.14 , 0.14 , 0.14 , 0.135, 0.135,
        0.13 , 0.13 , 0.125, 0.12 , 0.115, 0.115, 0.11 , 0.105, 0.1  ,
        0.095, 0.085, 0.08 , 0.07 , 0.055, 0.04    ]

    W_in = [
        0.017 , 0.0266, 0.0344, 0.0417, 0.0495, 0.0546, 0.061 , 0.0668,
        0.0721, 0.0768, 0.0849, 0.0887, 0.092 , 0.0996, 0.1019, 0.1037,
        0.1106, 0.1114, 0.1179, 0.1243, 0.1238, 0.1297, 0.1355, 0.1335,
        0.1387, 0.1437, 0.1485, 0.1448, 0.1489, 0.1527, 0.1564, 0.1597,
        0.1628, 0.1565, 0.1588, 0.1609, 0.1626, 0.1641, 0.1652, 0.166 ,
        0.1665, 0.1667, 0.1666, 0.1662, 0.1655, 0.1736, 0.1722, 0.1706,
        0.1686, 0.1664, 0.1639, 0.1698, 0.1667, 0.1634, 0.1599, 0.1641,
        0.1601, 0.1559, 0.1589, 0.1543, 0.1496, 0.1514, 0.1464, 0.1475,
        0.1422, 0.1425, 0.1424, 0.1419, 0.1361, 0.135 , 0.1335, 0.1316,
        0.1294, 0.1302, 0.1271, 0.1264, 0.1269, 0.1256    ]

    T_in = [
        299.2557, 298.775 , 298.2969, 297.8227, 297.3536, 296.8843,
        296.421 , 295.9603, 295.502 , 295.0456, 294.5994, 294.1468,
        293.6951, 293.2556, 292.8054, 292.3549, 291.9188, 291.4677,
        291.0325, 290.5978, 290.1434, 289.7073, 289.2706, 288.81  ,
        288.3698, 287.928 , 287.4842, 287.0118, 286.5622, 286.1099,
        285.6544, 285.1957, 284.7335, 284.2379, 283.7677, 283.2937,
        282.8157, 282.3337, 281.8476, 281.3574, 280.8631, 280.3649,
        279.8626, 279.3565, 278.8467, 278.362 , 277.8447, 277.3241,
        276.8006, 276.2742, 275.7454, 275.2388, 274.705 , 274.1694,
        273.6327, 273.1155, 272.576 , 272.0363, 271.514 , 270.9736,
        270.4339, 269.9094, 269.3711, 268.8465, 268.311 , 267.7877,
        267.2649, 266.7432, 266.2159, 265.698 , 265.1821, 264.6685,
        264.1574, 263.6518, 263.1461, 262.6451, 262.1476, 261.6524]
    #! format: on

    Area_in = pyinterp(self.grid.zc, z_in, Area_in)
    thetal_in = pyinterp(self.grid.zc, z_in, thetal_in)
    T_in = pyinterp(self.grid.zc, z_in, T_in)
    @inbounds for i in xrange(self.n_updrafts)
        @inbounds for k in real_face_indices(self.grid)
            if minimum(z_in) <= self.grid.zf[k] <= maximum(z_in)
                self.W.values[i, k] = 0.0
            end
        end

        @inbounds for k in real_center_indices(self.grid)
            if minimum(z_in) <= self.grid.zc[k] <= maximum(z_in)
                self.Area.values[i, k] = Area_in[k] #self.updraft_fraction/self.n_updrafts
                self.H.values[i, k] = thetal_in[k]
                self.QT.values[i, k] = 0.0
                self.QL.values[i, k] = 0.0

                # for now temperature is provided as diagnostics from LES
                self.T.values[i, k] = T_in[k]
            else
                self.Area.values[i, k] = 0.0 #self.updraft_fraction/self.n_updrafts
                self.H.values[i, k] = GMV.H.values[k]
                self.T.values[i, k] = GMV.T.values[k]
            end
        end
    end

    set_means(tptke, self, GMV)

    return
end

function initialize_io(self::UpdraftVariables, Stats::NetCDFIO_Stats)
    add_profile(Stats, "updraft_area")
    add_profile(Stats, "updraft_w")
    add_profile(Stats, "updraft_qt")
    add_profile(Stats, "updraft_ql")
    add_profile(Stats, "updraft_RH")
    add_profile(Stats, "updraft_thetal")
    add_profile(Stats, "updraft_temperature")
    add_profile(Stats, "updraft_buoyancy")
    add_profile(Stats, "updraft_cloud_fraction")

    add_ts(Stats, "updraft_cloud_cover")
    add_ts(Stats, "updraft_cloud_base")
    add_ts(Stats, "updraft_cloud_top")
    add_ts(Stats, "updraft_lwp")
    return
end

function set_means(tptke, self::UpdraftVariables, GMV::GridMeanVariables)

    self.Area.bulkvalues .= up_sum(self.Area.values)
    self.W.bulkvalues .= 0.0
    self.QT.bulkvalues .= 0.0
    self.QL.bulkvalues .= 0.0
    self.H.bulkvalues .= 0.0
    self.T.bulkvalues .= 0.0
    self.B.bulkvalues .= 0.0
    self.RH.bulkvalues .= 0.0
    grid = get_grid(GMV)

    @inbounds for k in real_center_indices(self.grid)
        if self.Area.bulkvalues[k] > 1.0e-20
            @inbounds for i in xrange(self.n_updrafts)
                self.QT.bulkvalues[k] += self.Area.values[i, k] * self.QT.values[i, k] / self.Area.bulkvalues[k]
                self.QL.bulkvalues[k] += self.Area.values[i, k] * self.QL.values[i, k] / self.Area.bulkvalues[k]
                self.H.bulkvalues[k] += self.Area.values[i, k] * self.H.values[i, k] / self.Area.bulkvalues[k]
                self.T.bulkvalues[k] += self.Area.values[i, k] * self.T.values[i, k] / self.Area.bulkvalues[k]
                self.RH.bulkvalues[k] += self.Area.values[i, k] * self.RH.values[i, k] / self.Area.bulkvalues[k]
                self.B.bulkvalues[k] += self.Area.values[i, k] * self.B.values[i, k] / self.Area.bulkvalues[k]
            end
        else
            self.QT.bulkvalues[k] = GMV.QT.values[k]
            self.QL.bulkvalues[k] = 0.0
            self.H.bulkvalues[k] = GMV.H.values[k]
            self.RH.bulkvalues[k] = GMV.RH.values[k]
            self.T.bulkvalues[k] = GMV.T.values[k]
            self.B.bulkvalues[k] = 0.0
        end
        if self.QL.bulkvalues[k] > 1e-8 && self.Area.bulkvalues[k] > 1e-3
            self.cloud_fraction[k] = 1.0
        else
            self.cloud_fraction[k] = 0.0
        end
    end

    @inbounds for k in real_face_indices(self.grid)
        a_bulk_bcs = (; bottom = SetValue(sum(tptke.area_surface_bc)), top = SetZeroGradient())
        a_bulk_f = interpc2f(self.Area.bulkvalues, grid, k; a_bulk_bcs...)
        if a_bulk_f > 1.0e-20
            @inbounds for i in xrange(self.n_updrafts)
                a_up_bcs = (; bottom = SetValue(tptke.area_surface_bc[i]), top = SetZeroGradient())
                a_up_f = interpc2f(self.Area.values, grid, k, i; a_up_bcs...)
                self.W.bulkvalues[k] += a_up_f * self.W.values[i, k] / a_bulk_f
            end
        else
            self.W.bulkvalues[k] = 0.0
        end
    end
    return
end

# quick utility to set "new" arrays with values in the "values" arrays
function set_new_with_values(self::UpdraftVariables)
    @inbounds for i in xrange(self.n_updrafts)
        @inbounds for k in real_face_indices(self.grid)
            self.W.new[i, k] = self.W.values[i, k]
        end

        @inbounds for k in real_center_indices(self.grid)
            self.Area.new[i, k] = self.Area.values[i, k]
            self.QT.new[i, k] = self.QT.values[i, k]
            self.H.new[i, k] = self.H.values[i, k]
        end
    end
    return
end

# quick utility to set "new" arrays with values in the "values" arrays
function set_old_with_values(self::UpdraftVariables)
    @inbounds for i in xrange(self.n_updrafts)
        @inbounds for k in real_center_indices(self.grid)
            self.Area.old[i, k] = self.Area.values[i, k]
        end
    end
    return
end

# quick utility to set "tmp" arrays with values in the "new" arrays
function set_values_with_new(self::UpdraftVariables)
    @inbounds for i in xrange(self.n_updrafts)
        @inbounds for k in real_face_indices(self.grid)
            self.W.values[i, k] = self.W.new[i, k]
        end

        @inbounds for k in real_center_indices(self.grid)
            self.W.values[i, k] = self.W.new[i, k]
            self.Area.values[i, k] = self.Area.new[i, k]
            self.QT.values[i, k] = self.QT.new[i, k]
            self.H.values[i, k] = self.H.new[i, k]
        end
    end
    return
end

function io(self::UpdraftVariables, Stats::NetCDFIO_Stats, ref_state::ReferenceState)
    write_profile(Stats, "updraft_area", self.Area.bulkvalues)
    write_profile(Stats, "updraft_w", self.W.bulkvalues)
    write_profile(Stats, "updraft_qt", self.QT.bulkvalues)
    write_profile(Stats, "updraft_ql", self.QL.bulkvalues)
    write_profile(Stats, "updraft_RH", self.RH.bulkvalues)
    write_profile(Stats, "updraft_thetal", self.H.bulkvalues)
    write_profile(Stats, "updraft_temperature", self.T.bulkvalues)
    write_profile(Stats, "updraft_buoyancy", self.B.bulkvalues)

    upd_cloud_diagnostics(self, ref_state)
    write_profile(Stats, "updraft_cloud_fraction", self.cloud_fraction)
    # Note definition of cloud cover : each updraft is associated with a cloud cover equal to the maximum
    # area fraction of the updraft where ql > 0. Each updraft is assumed to have maximum overlap with respect to
    # itself (i.e. no consideration of tilting due to shear) while the updraft classes are assumed to have no overlap
    # at all. Thus total updraft cover is the sum of each updraft"s cover
    write_ts(Stats, "updraft_cloud_cover", sum(self.cloud_cover))
    write_ts(Stats, "updraft_cloud_base", minimum(abs.(self.cloud_base)))
    write_ts(Stats, "updraft_cloud_top", maximum(abs.(self.cloud_top)))
    write_ts(Stats, "updraft_lwp", self.lwp)
    return
end

function upd_cloud_diagnostics(self::UpdraftVariables, ref_state::ReferenceState)
    self.lwp = 0.0

    @inbounds for i in xrange(self.n_updrafts)
        self.cloud_base[i] = zc_toa(self.grid)
        self.cloud_top[i] = 0.0
        self.updraft_top[i] = 0.0
        self.cloud_cover[i] = 0.0

        @inbounds for k in real_center_indices(self.grid)
            if self.Area.values[i, k] > 1e-3
                self.updraft_top[i] = max(self.updraft_top[i], self.grid.zc[k])
                self.lwp += ref_state.rho0_half[k] * self.QL.values[i, k] * self.Area.values[i, k] * self.grid.Δz

                if self.QL.values[i, k] > 1e-8
                    self.cloud_base[i] = min(self.cloud_base[i], self.grid.zc[k])
                    self.cloud_top[i] = max(self.cloud_top[i], self.grid.zc[k])
                    self.cloud_cover[i] = max(self.cloud_cover[i], self.Area.values[i, k])
                end
            end
        end
    end

    return
end

"""
clear precipitation source terms for QT and H from each updraft
"""
function clear_precip_sources(self::UpdraftThermodynamics)
    self.prec_source_qt .= 0
    self.prec_source_h .= 0
    return
end

"""
sum precipitation source terms for QT and H from all sub-timesteps
"""
function update_total_precip_sources(self::UpdraftThermodynamics)
    self.prec_source_h_tot .= up_sum(self.prec_source_h)
    self.prec_source_qt_tot .= up_sum(self.prec_source_qt)
    return
end

function buoyancy(
    self::UpdraftThermodynamics,
    UpdVar::UpdraftVariables,
    EnvVar::EnvironmentVariables,
    GMV::GridMeanVariables,
    extrap::Bool,
)

    grid = self.grid
    kc_surf = kc_surface(grid)
    param_set = parameter_set(GMV)

    UpdVar.Area.bulkvalues .= up_sum(UpdVar.Area.values)

    if !extrap
        @inbounds for i in xrange(self.n_updraft)
            @inbounds for k in real_center_indices(grid)
                if UpdVar.Area.values[i, k] > 0.0
                    ts = TD.PhaseEquil_pTq(
                        param_set,
                        self.ref_state.p0_half[k],
                        UpdVar.T.values[i, k],
                        UpdVar.QT.values[i, k],
                    )
                    rho = TD.air_density(ts)
                    UpdVar.B.values[i, k] = buoyancy_c(param_set, self.ref_state.rho0_half[k], rho)
                else
                    UpdVar.B.values[i, k] = EnvVar.B.values[k]
                end
                ts = TD.PhaseEquil_pθq(
                    param_set,
                    self.ref_state.p0_half[k],
                    UpdVar.H.values[i, k],
                    UpdVar.QT.values[i, k],
                )
                UpdVar.RH.values[i, k] = TD.relative_humidity(ts)
            end
        end
    else
        @inbounds for i in xrange(self.n_updraft)
            @inbounds for k in real_center_indices(grid)
                if UpdVar.Area.values[i, k] > 0.0
                    qt = UpdVar.QT.values[i, k]
                    h = UpdVar.H.values[i, k]
                    ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[k], h, qt)
                    rho = TD.air_density(ts)
                    UpdVar.B.values[i, k] = buoyancy_c(param_set, self.ref_state.rho0_half[k], rho)
                    UpdVar.RH.values[i, k] = TD.relative_humidity(ts)
                elseif k > kc_surf
                    if UpdVar.Area.values[i, k - 1] > 0.0
                        qt = UpdVar.QT.values[i, k - 1]
                        h = UpdVar.H.values[i, k - 1]
                        ts = TD.PhaseEquil_pθq(param_set, self.ref_state.p0_half[k], h, qt)
                        t = TD.air_temperature(ts)
                        rho = TD.air_density(ts)
                        UpdVar.B.values[i, k] = buoyancy_c(param_set, self.ref_state.rho0_half[k], rho)
                        UpdVar.RH.values[i, k] = TD.relative_humidity(ts)
                    else
                        UpdVar.B.values[i, k] = EnvVar.B.values[k]
                        UpdVar.RH.values[i, k] = EnvVar.RH.values[k]
                    end
                else
                    UpdVar.B.values[i, k] = EnvVar.B.values[k]
                    UpdVar.RH.values[i, k] = EnvVar.RH.values[k]
                end
            end
        end
    end


    @inbounds for k in real_center_indices(self.grid)
        GMV.B.values[k] = (1.0 - UpdVar.Area.bulkvalues[k]) * EnvVar.B.values[k]
        @inbounds for i in xrange(self.n_updraft)
            GMV.B.values[k] += UpdVar.Area.values[i, k] * UpdVar.B.values[i, k]
        end
        @inbounds for i in xrange(self.n_updraft)
            UpdVar.B.values[i, k] -= GMV.B.values[k]
        end
        EnvVar.B.values[k] -= GMV.B.values[k]
    end

    return
end

"""
compute precipitation source terms
"""
function microphysics(self::UpdraftThermodynamics, UpdVar::UpdraftVariables, Rain::RainVariables, dt)
    rst = rain_struct()
    mph = mph_struct()
    param_set = parameter_set(Rain)

    @inbounds for i in xrange(self.n_updraft)
        @inbounds for k in real_center_indices(self.grid)

            # autoconversion and accretion
            mph = microphysics_rain_src(
                param_set,
                Rain.rain_model,
                UpdVar.QT.values[i, k],
                UpdVar.QL.values[i, k],
                Rain.Upd_QR.values[k],
                UpdVar.Area.values[i, k],
                UpdVar.T.values[i, k],
                self.ref_state.p0_half[k],
                self.ref_state.rho0_half[k],
                dt,
            )

            # update rain sources of state variables
            self.prec_source_qt[i, k] -= mph.qr_src * UpdVar.Area.values[i, k]
            self.prec_source_h[i, k] += mph.thl_rain_src * UpdVar.Area.values[i, k]
        end
    end
    return
end
