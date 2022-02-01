import TurbulenceConvection
const TC = TurbulenceConvection
import Thermodynamics
const TD = Thermodynamics

function initialize_turbconv(
    edmf::TC.EDMF_PrognosticTKE,
    grid::TC.Grid,
    state::TC.State,
    case,
    gm::TC.GridMeanVariables,
    t::Real,
)
    surf_params = case.surf_params
    param_set = TC.parameter_set(gm)
    aux_tc = TC.center_aux_turbconv(state)
    prog_gm = TC.center_prog_grid_mean(state)
    p0_c = TC.center_ref_state(state).p0
    parent(aux_tc.prandtl_nvec) .= edmf.prandtl_number
    @inbounds for k in TC.real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        aux_tc.θ_virt[k] = TD.virtual_pottemp(ts)
    end
    surf = get_surface(surf_params, grid, state, gm, t, param_set)
    TC.initialize_edmf(edmf, grid, state, case, gm, surf)
    return
end
