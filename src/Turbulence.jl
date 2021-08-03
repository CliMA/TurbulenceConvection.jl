
function ParameterizationFactory(namelist, Gr::Grid, Ref::ReferenceState, param_set::PS) where {PS}
    scheme = namelist["turbulence"]["scheme"]
    return EDMF_PrognosticTKE(namelist, Gr, Ref, param_set)
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


function update_GMV_diagnostics(self, GMV)
    return
end
