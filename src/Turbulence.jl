
function ParameterizationFactory(namelist, Gr::Grid, Ref::ReferenceState, param_set::PS) where {PS}
    scheme = namelist["turbulence"]["scheme"]
    if scheme == "EDMF_PrognosticTKE"
        return EDMF_PrognosticTKE(namelist, Gr, Ref, param_set)
    elseif scheme == "SimilarityED"
        return SimilarityED(namelist, Gr, Ref, param_set)
    else
        error("Did not recognize parameterization $scheme")
    end
end

# Calculate the tendency of the grid mean variables due to turbulence as
# the difference between the values at the beginning and  end of all substeps taken
function update(self::ParameterizationBase, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    @inbounds for k in real_center_indicies(self.Gr)
        GMV.H.tendencies[k] += (GMV.H.new[k] - GMV.H.values[k]) * TS.dti
        GMV.QT.tendencies[k] += (GMV.QT.new[k] - GMV.QT.values[k]) * TS.dti
        GMV.U.tendencies[k] += (GMV.U.new[k] - GMV.U.values[k]) * TS.dti
        GMV.V.tendencies[k] += (GMV.V.new[k] - GMV.V.values[k]) * TS.dti
    end

    return
end

# Update the diagnosis of the inversion height, using the maximum temperature gradient method
function update_inversion(self::ParameterizationBase, GMV::GridMeanVariables, option)
    theta_rho = center_field(self.Gr)
    ∇θ_liq_max = 0.0
    k_fi = first_center(self.Gr)

    @inbounds for k in real_center_indicies(self.Gr)
        qv = GMV.QT.values[k] - GMV.QL.values[k]
        theta_rho[k] = theta_rho_c(self.Ref.p0_half[k], GMV.T.values[k], GMV.QT.values[k], qv)
    end


    if option == "theta_rho"
        @inbounds for k in real_center_indicies(self.Gr)
            if theta_rho[k] > theta_rho[k_fi]
                self.zi = self.Gr.z_half[k]
                break
            end
        end
    elseif option == "thetal_maxgrad"

        @inbounds for k in real_center_indicies(self.Gr)
            ∇θ_liq = ∇_upwind(GMV.THL.values, self.Gr, k)
            if ∇θ_liq > ∇θ_liq_max
                ∇θ_liq_max = ∇θ_liq
                self.zi = self.Gr.z[k]
            end
        end
    elseif option == "critical_Ri"
        self.zi = get_inversion(theta_rho, GMV.U.values, GMV.V.values, self.Gr, Ri_bulk_crit(self))

    else
        error("INVERSION HEIGHT OPTION NOT RECOGNIZED")
    end

    return
end



# Compute eddy diffusivities from similarity theory (Siebesma 2007)
function compute_eddy_diffusivities_similarity(self::ParameterizationBase, GMV::GridMeanVariables, Case::CasesBase)
    update_inversion(self, GMV, Case.inversion_option)
    self.wstar = get_wstar(Case.Sur.bflux, self.zi)

    ustar = Case.Sur.ustar
    KH = diffusivity_h(self)
    KM = diffusivity_m(self)
    @inbounds for k in real_center_indicies(self.Gr)
        zzi = self.Gr.z_half[k] / self.zi
        if zzi <= 1.0
            if self.wstar < 1e-6
                KH.values[k] = 0.0
                KM.values[k] = 0.0
            else
                KH.values[k] =
                    vkb *
                    ((ustar / self.wstar)^3 + 39.0 * vkb * zzi)^(1.0 / 3.0) *
                    zzi *
                    (1.0 - zzi) *
                    (1.0 - zzi) *
                    self.wstar *
                    self.zi
                KM.values[k] = KH.values[k] * prandtl_number(self)
            end
        else
            KH.values[k] = 0.0
            KM.values[k] = 0.0
        end
    end


    # Set the boundary points at top and bottom of domain
    set_bcs(KH, self.Gr)
    set_bcs(KM, self.Gr)
    return
end

function update_GMV_diagnostics(self, GMV)
    return
end


#####################################################################################################################

function SimilarityED(namelist, Gr::Grid, Ref::ReferenceState, param_set::PS) where {PS}
    extrapolate_buoyancy = false
    base = ParameterizationBase(namelist, Gr, Ref)
    return SimilarityED{PS}(param_set, base, extrapolate_buoyancy)
end

initialize(self::SimilarityED, Case::CasesBase, GMV::GridMeanVariables, Ref::ReferenceState) = nothing

function initialize_io(self::SimilarityED, Stats::NetCDFIO_Stats)
    add_profile(Stats, "eddy_viscosity")
    add_profile(Stats, "eddy_diffusivity")
    return
end

function io(self::SimilarityED, Stats::NetCDFIO_Stats, TS::TimeStepping)
    write_profile(Stats, "eddy_viscosity", diffusivity_m(self).values[(self.Gr.gw):(self.Gr.nzg - self.Gr.gw)])
    write_profile(Stats, "eddy_diffusivity", diffusivity_h(self).values[(self.Gr.gw):(self.Gr.nzg - self.Gr.gw)])
    return
end

function update(self::SimilarityED, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    set_bcs(GMV.H, self.Gr)
    set_bcs(GMV.QT, self.Gr)
    set_bcs(GMV.THL, self.Gr)

    compute_eddy_diffusivities_similarity(self.base, GMV, Case)

    gw = self.Gr.gw
    nzg = self.Gr.nzg
    nz = self.Gr.nz
    a = pyzeros(nz)
    b = pyzeros(nz)
    c = pyzeros(nz)
    x = pyzeros(nz)
    dummy_ae = pyones(nzg)
    rho_K_m = pyzeros(nzg)
    KH = diffusivity_h(self).values

    @inbounds for k in face_indicies(self.Gr)
        rho_K_m[k] = 0.5 * (KH[k] + KH[k + 1]) * self.Ref.rho0[k]
    end


    # Matrix is the same for all variables that use the same eddy diffusivity
    construct_tridiag_diffusion(self.Gr, TS.dt, rho_K_m, self.Ref.rho0_half, dummy_ae, a, b, c)

    # Solve QT

    @inbounds for k in xrange(nz)
        x[k] = GMV.QT.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_qtflux * self.Gr.dzi * self.Ref.alpha0_half[gw]

    x .= tridiag_solve(x, a, b, c)

    @inbounds for k in xrange(nz)
        GMV.QT.new[k + gw] = x[k]
    end

    # Solve H
    @inbounds for k in xrange(nz)
        x[k] = GMV.H.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_hflux * self.Gr.dzi * self.Ref.alpha0_half[gw]

    x .= tridiag_solve(x, a, b, c)
    @inbounds for k in xrange(nz)
        GMV.H.new[k + gw] = x[k]
    end


    # Solve U
    @inbounds for k in xrange(nz)
        x[k] = GMV.U.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_uflux * self.Gr.dzi * self.Ref.alpha0_half[gw]

    x .= tridiag_solve(x, a, b, c)

    @inbounds for k in xrange(nz)
        GMV.U.new[k + gw] = x[k]
    end

    # Solve V
    @inbounds for k in xrange(nz)
        x[k] = GMV.V.values[k + gw]
    end
    x[0] = x[0] + TS.dt * Case.Sur.rho_vflux * self.Gr.dzi * self.Ref.alpha0_half[gw]

    x .= tridiag_solve(x, a, b, c)
    @inbounds for k in xrange(nz)
        GMV.V.new[k + gw] = x[k]
    end

    update_GMV_diagnostics(self, GMV)
    update(self.base, GMV, Case, TS)

    return
end

update_inversion(self::SimilarityED, GMV, option) = update_inversion(self.base, GMV, option)

function update_GMV_diagnostics(self::SimilarityED, GMV)
    # Ideally would write this to be able to use an SGS condensation closure, but unless the need arises,
    # we will just do an all-or-nothing treatment as a placeholder
    GMV.satadjust()
    return
end
