using Pkg
using Plots
Pkg.activate(expanduser("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test"))
import SOCRATESSingleColumnForcings as SSCF
import CLIMAParameters as CP # use CLIMAParameters = "0.7, 0.8, 0.9, 0.10"
# import ClimaParams as CPP # would using this trouble w/ TC.jl? it's a different uuid technically..., use ClimaParams = "0.10"
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP


FT = Float64
ylims = FT[1500, 3000]


toml_dict = CP.create_toml_dict(FT; dict_type = "alias") # CP 0.7 and below, Thermodynamics 0.11 and above
aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)

data = SSCF.open_atlas_les_input(10)
new_zc = data[:grid_data]

case_data = SSCF.process_case(
    10;
    forcing_type = :obs_data,
    new_z = (;
        dTdt_hadv = new_zc,
        H_nudge = new_zc,
        dqtdt_hadv = new_zc,
        qt_nudge = new_zc,
        subsidence = new_zc,
        u_nudge = new_zc,
        v_nudge = new_zc,
        ug_nudge = new_zc,
        vg_nudge = new_zc,
        dTdt_rad = new_zc,
    ),
    initial_condition = false,
    thermo_params = thermo_params,
    conservative_interp = true,
)


subsidence_init =  (f->f([0])[1]).(case_data[:subsidence])
using Plots
ENV["GKSwstype"]="nul"


plot(subsidence_init,
    new_zc,
    label = "subsidence",
    xlabel = "subsidence",
    ylabel = "z",
    dpi = 600,
    title = "subsidence",
    lw = 2,
    size = (800, 800),
    legend = :topleft,
    color = :black,
    fmt = :png,
    grid = true,
)

t = 12*3600
subsidence_final =  (f->f([t])[1]).(case_data[:subsidence])

plot!(subsidence_final, new_zc, label = "subsidence final", lw = 2, color = :red)
# plot!(ylims = ylims)
savefig("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test/subsidence_init.png")




using NCDatasets
test_data_file = "/home/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES_postprocess_runs_storage/subexperiments/SOCRATES_Base/Calibrate_and_Run/tau_autoconv_noneq/adapt_dt__dt_min_2.0__dt_max_4.0/iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean/postprocessing/output/Atlas_LES/RFAll_obs/best_particle/data/Output.SOCRATES_RF10_obs_data.1_1/stats/Stats.SOCRATES_RF10_obs_data.nc"
test_data = NCDatasets.Dataset(test_data_file, "r") do data
    # @info "data = $data"
    getvar(var; group = "profiles") = data.group[group][var][:]
    zc_les = getvar("zc"; group = "reference")
    zf_les = getvar("zf"; group = "reference")
    subsidence = getvar("subsidence")
    (subsidence, zc_les)
    qi_mean = getvar("qi_mean")

    qi_mean_ls_vert_adv = getvar("qi_mean_ls_vert_adv")

    # get last index
    t_ind = size(qi_mean_ls_vert_adv, 2)
    qi_mean_ls_vert_adv = qi_mean_ls_vert_adv[:, t_ind]
    qi_mean = qi_mean[:, t_ind]
    subsidence = subsidence[:, t_ind]
    return (;subsidence, zc_les, zf_les, qi_mean, qi_mean_ls_vert_adv)
end
subsidence, zc_les, zf_les, qi_mean, qi_mean_ls_vert_adv = test_data

diff(x) = x[2:end] .- x[1:end-1]

subsidence = SSCF.Dierckx.Spline1D(zc_les, subsidence; k = 1)
qi_mean = SSCF.Dierckx.Spline1D(zc_les, qi_mean; k = 1)
 
plot(qi_mean_ls_vert_adv,
    zc_les,
    label = "qi_mean_ls_vert_adv",
    xlabel = "qi_mean_ls_vert_adv",
    ylabel = "z",
    dpi = 600,
    title = "qi_mean_ls_vert_adv",
    lw = 2,
    size = (800, 800),
    legend = :topleft,
    color = :black,
    fmt = :png,
    grid = true,
    marker = :circle,
    markerstrokewidth = 0
)

derivative = SSCF.Dierckx.derivative
qi_mean_ls_vert_adv_my_calc1 = -derivative(qi_mean, zc_les) .* subsidence(zc_les)
qi_mean_ls_vert_adv_my_calc2 = -derivative(subsidence, zc_les) .* qi_mean(zc_les)
qi_mean_ls_vert_adv_my_calc3 = -derivative(SSCF.Dierckx.Spline1D(zc_les, qi_mean(zc_les) .* subsidence(zc_les); k = 1), zc_les)
plot!(qi_mean_ls_vert_adv_my_calc1, zc_les, label = "qi_mean_ls_vert_adv_my_calc advection", lw = 2, color = :red, marker = :circle, markerstrokewidth = 0)
plot!(qi_mean_ls_vert_adv_my_calc2, zc_les, label = "qi_mean_ls_vert_adv_my_calc convergence", lw = 2, color = :blue, marker = :circle, markerstrokewidth = 0)
plot!(qi_mean_ls_vert_adv_my_calc3, zc_les, label = "qi_mean_ls_vert_adv_my_calc full", lw = 2, color = :green, marker = :circle, markerstrokewidth = 0)

# qi_mean_ls_vert_adv_my_calc4 = -derivative(qi_mean, zc_les) # .* subsidence(zc_les).^3
# qi_mean_ls_vert_adv_my_calc4 =  -derivative(subsidence, zc_les)

v = subsidence(zc_les)
a⁺ = [derivative(qi_mean, zc_les[1:end-1])..., FT(0)]
a⁻ = [FT(0), derivative(qi_mean, zc_les[2:end])...]
qi_mean_ls_vert_adv_my_calc4 = -(( (v + abs.(v)) .* a⁻ .+ (v - abs.(v)) .* a⁺) ) ./ 2

plot!(qi_mean_ls_vert_adv_my_calc4, zc_les, label = "qi_mean_ls_vert_adv_my_calc full", lw = 2, color = :purple, marker = :circle, markerstrokewidth = 0)

savefig("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test/qi_mean_ls_vert_adv.png")


# ==================================================================================================== #

Pkg.activate(expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl/integration_tests"))
import ClimaCore as CC
import ClimaCore.Operators as CCO


# new_z = zc_les
# z_mesh = CC.Geometry.ZPoint{FT}.(new_z) # added 0 to beginning? copy from the file #Array(TC.get_nc_data(data, "zc")) also idk what to do about paths like this
# # z_mesh = CC.Geometry.ZPoint{FT}.([FT(0), new_z...]) # added 0 to beginning? copy from the file #Array(TC.get_nc_data(data, "zc")) also idk what to do about paths like this (not needed if we do zc to zf conversion)
# nz = length(z_mesh)
# z₀, z₁ = z_mesh[1], z_mesh[end]
# zmax = z_mesh[end]
# domain = CC.Domains.IntervalDomain(
#     CC.Geometry.ZPoint{FT}(z₀),
#     CC.Geometry.ZPoint{FT}(z₁),
#     boundary_tags = (:bottom, :top),
# )
# z_mesh = CC.Meshes.IntervalMesh(domain, z_mesh)

# using ClimaComms

# device = ClimaComms.device()
# CC.Fields.coordinate_field(z_mesh)

import TurbulenceConvection as TC
TCP = TC.Parameters
import JSON
namelist_file = "/central/groups/esm/jbenjami/Research_Schneider/CliMA/CalibrateEDMF.jl/experiments/SOCRATES_postprocess_runs_storage/subexperiments/SOCRATES_Base/Calibrate_and_Run/tau_autoconv_noneq/adapt_dt__dt_min_2.0__dt_max_4.0/iwp_mean__lwp_mean__qi_mean__qip_mean__ql_mean__qr_mean/postprocessing/output/Atlas_LES/RFAll_obs/best_particle/data/Output.SOCRATES_RF10_obs_data.1_1/namelist_SOCRATES.in"
namelist = JSON.parsefile(namelist_file, dicttype = Dict, inttype = Int, null = FT(NaN))
namelist["meta"]["casename"]

namelist["meta"]["forcing_type"] = Symbol(namelist["meta"]["forcing_type"])

toml_dict = CP.create_toml_dict(FT; dict_type = "alias")

tc_dir = expanduser("~/Research_Schneider/CliMA/TurbulenceConvection.jl")

using ForwardDiff # for netcdf.io
include(joinpath(tc_dir, "driver", "NetCDFIO.jl"))
include(joinpath(tc_dir, "driver", "generate_namelist.jl"))
include(joinpath(tc_dir, "driver", "Cases.jl"))
include(joinpath(tc_dir, "driver", "parameter_set.jl"))
include(joinpath(tc_dir, "driver", "dycore.jl"))
param_set = create_parameter_set(namelist, toml_dict, FT) # this creates an override file in  a directory we don' need...
thermo_params = TCP.thermodynamics_params(param_set)

include(joinpath(tc_dir, "driver", "common_spaces.jl"))
spaces = get_spaces(namelist, param_set, FT)

@info "spaces = $spaces"

space_c = spaces[1]
space_f = spaces[2]

coords = CC.Fields.coordinate_field.(spaces)
coords_c = coords[1]
coords_f = coords[2]


subsidence_c = similar(coords_c.z)
subsidence_c_parent = parent(subsidence_c)
@. subsidence_c_parent = subsidence(zc_les)

qi_c = similar(coords_c.z)
qi_c_parent = parent(qi_c)
@. qi_c_parent = qi_mean(zc_les)


C123 = CCG.Covariant123Vector
C12 = CCG.Contravariant12Vector
q_liq_gm_toa = FT(0.0)
q_ice_gm_toa = FT(0.0)
q_liq_boa = FT(0.0)
q_ice_boa = FT(0.0)
# C2Fq_liq = CCO.InterpolateC2F(bottom = CCO.SetValue(FT(q_liq_boa)), top = CCO.SetValue(FT(q_liq_gm_toa))) # not sure if this should be right biased or not
# C2Fq_ice = CCO.InterpolateC2F(bottom = CCO.SetValue(FT(q_ice_boa)), top = CCO.SetValue(FT(q_ice_gm_toa))) # not sure if this should be right biased or not
C2Fq_liq = CCO.InterpolateC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
C2Fq_ice = CCO.InterpolateC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
F2Cq_liq = CCO.InterpolateF2C(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
F2Cq_ice = CCO.InterpolateF2C(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
wvec = CC.Geometry.WVector
∇c = CCO.DivergenceF2C() # F2C to come back from C2F

RBq_ice = CCO.RightBiasedC2F(; top = CCO.SetValue(q_ice_gm_toa))
LBq_ice = CCO.LeftBiasedC2F(; bottom = CCO.SetValue(q_ice_boa))
C2Fsub = CCO.InterpolateC2F(; bottom = CCO.SetValue(FT(0)), top = CCO.Extrapolate()) # not sure but maybe we should use this for stability as it's not biased for your upwind algorithm... no penetration. # I think this hsould be the more stable option? Though according to cgarkue github issue, it doesn't always work (see https://github.com/CliMA/TurbulenceConvection.jl/pull/1146)

# F2Csub = CCO.InterpolateF2C(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # We might have put 0 for no penetration on boundary, but that's not exactly true in our dataset...
F2Csub = Ic # don't need bcs for F2c
CV32FT = x -> x[1] # convert Contravariant3Vector to Float64

# Note C2F is lossy, peaks will be cut off and diffused out. So it's not exactly a harmless opearation like RB/LB. Just doing the full sol'n (see _5) is conservative.
∇q_ice_gm_c = similar(coords_c.z)
@. ∇q_ice_gm_c = ∇c(wvec(C2Fq_ice(qi_c)))
∇q_ice_gm_c_parent = parent(∇q_ice_gm_c)

UBsub = CCO.UpwindBiasedProductC2F(; bottom = CCO.Extrapolate(), top = CCO.Extrapolate()) # upwinding, extrapolate bc we don't know the boa/toa derivatives

# Base upwind 
# w∇q_ice_gm_c = @. UBsub(wvec(C2Fsub(subsidence_c)), ∇q_ice_gm_c) # works, but output type is different than tendencies type
w∇q_ice_gm_c = @. wvec(UBsub(wvec(C2Fsub(subsidence_c)), ∇q_ice_gm_c)) # works, but output type is different than tendencies type
w∇q_ice_gm_c = @. F2Csub(w∇q_ice_gm_c) # convert back to C
w∇q_ice_gm_c = @. CV32FT(w∇q_ice_gm_c) # convert back to Float64 from Contravariant3Vector
w∇q_ice_gm_c_parent = parent(w∇q_ice_gm_c)

# 2 | flip which is vector
# w∇q_ice_gm_c_2 = @. UBsub(wvec(C2Fsub(∇q_ice_gm_c)), subsidence_c) # works, but output type is different than tendencies type
w∇q_ice_gm_c_2 = @. wvec(UBsub(wvec(C2Fq_ice(∇q_ice_gm_c)), subsidence_c)) # works, but output type is different than tendencies type
w∇q_ice_gm_c_2 = @. F2Csub(w∇q_ice_gm_c_2) # convert back to C
w∇q_ice_gm_c_2 = @. CV32FT(w∇q_ice_gm_c_2) # convert back to Float64 from Contravariant3Vector
w∇q_ice_gm_c_2_parent = parent(w∇q_ice_gm_c_2)

# right biased
∇q_ice_gm_c_RB = @. ∇c(wvec(RBq_ice(qi_c)))
∇q_ice_gm_c_RB_parent = parent(∇q_ice_gm_c_RB)
w∇q_ice_gm_c_RB = @. ∇q_ice_gm_c_RB * subsidence_c
w∇q_ice_gm_c_RB_parent = parent(w∇q_ice_gm_c_RB)

# left biased
∇q_ice_gm_c_LB = @. ∇c(wvec(LBq_ice(qi_c)))
∇q_ice_gm_c_LB_parent = parent(∇q_ice_gm_c_LB)
w∇q_ice_gm_c_LB = @. ∇q_ice_gm_c_LB * subsidence_c
w∇q_ice_gm_c_LB_parent = parent(w∇q_ice_gm_c_LB)

# 3 | switch order of wvec and F2Csub
w∇q_ice_gm_c_3 = @. UBsub(wvec(C2Fq_ice(∇q_ice_gm_c)), subsidence_c) # works, but output type is different than tendencies type
w∇q_ice_gm_c_3 = @. wvec(F2Csub(w∇q_ice_gm_c_3)) # convert back to C
w∇q_ice_gm_c_3 = @. CV32FT(w∇q_ice_gm_c_3) # convert back to Float64 from Contravariant3Vector
w∇q_ice_gm_c_3_parent = parent(w∇q_ice_gm_c_3)

# 4 | No bias
w∇q_ice_gm_c_4 = @. F2Csub(wvec(C2Fq_ice(∇q_ice_gm_c))) * subsidence_c
w∇q_ice_gm_c_4_parent = parent(w∇q_ice_gm_c_4)

#5 | upwind biased product then gradient (so full contribution)
# w∇q_ice_gm_c_5 = @. UBsub(wvec(C2Fsub(subsidence_c)), qi_c) # works, but output type is different than tendencies type
w∇q_ice_gm_c_5 = @. wvec(UBsub(wvec(C2Fsub(subsidence_c)), qi_c)) # works, but output type is different than tendencies type
w∇q_ice_gm_c_5 = @. ∇c(w∇q_ice_gm_c_5)
w∇q_ice_gm_c_5_parent = parent(w∇q_ice_gm_c_5)

# 6 | total contribution constructed from base and flipped
# w∇q_ice_gm_c_6 = @. w∇q_ice_gm_c + w∇q_ice_gm_c_2
# w∇q_ice_gm_c_6_parent = parent(w∇q_ice_gm_c_6)

# Note things may not be exactly conservative at the edges, these are not conservative round trip interpolations, you inherently are losing information from/cutting off peaks with each F2C/C2F.

plot(qi_mean_ls_vert_adv_my_calc3[:],
    zc_les,
    label = "qi_mean_ls_vert_adv",
    xlabel = "w∇q_ice_gm",
    ylabel = "z",
    dpi = 600,
    title = "w∇q_ice_gm",
    lw = 5,
    size = (800, 800),
    legend = :topleft,
    color = :black,
    fmt = :png,
    grid = true,
    marker = :circle,
    markerstrokewidth = 0
)
plot!(-w∇q_ice_gm_c_parent[:], zc_les, label = "w∇q_ice_gm_c", lw = 2, color = :red, marker = :circle, markerstrokewidth = 0)
plot!(-w∇q_ice_gm_c_RB_parent[:], zc_les, label = "w∇q_ice_gm_c_RB", lw = 2, color = :pink, marker = :circle, markerstrokewidth = 0, linestyle = :dash)
plot!(-w∇q_ice_gm_c_LB_parent[:], zc_les, label = "w∇q_ice_gm_c_LB", lw = 2, color = :cyan, marker = :circle, markerstrokewidth = 0, linestyle = :dot)
plot!(-w∇q_ice_gm_c_2_parent[:], zc_les, label = "w∇q_ice_gm_c_2", lw = 2, color = :blue, marker = :circle, markerstrokewidth = 0, linsestyle = :dash)
plot!(-w∇q_ice_gm_c_3_parent[:], zc_les, label = "w∇q_ice_gm_c_3", lw = 2, color = :green, marker = :circle, markerstrokewidth = 0, linestyle=:dot)
plot!(-w∇q_ice_gm_c_4_parent[:], zc_les, label = "w∇q_ice_gm_c_4", lw = 2, color = :purple, marker = :circle, markerstrokewidth = 0)
plot!(-w∇q_ice_gm_c_5_parent[:], zc_les, label = "w∇q_ice_gm_c_5", lw = 2, color = :orange, marker = :circle, markerstrokewidth = 0)
# plot!(-w∇q_ice_gm_c_6_parent[:], zc_les, label = "w∇q_ice_gm_c_6", lw = 2, color = :brown, marker = :circle, markerstrokewidth = 0)
plot!(ylims = ylims)
savefig("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test/w∇q_ice_gm.png")


plot(∇q_ice_gm_c_parent[:],
    zc_les,
    label = "∇q_ice_gm_c",
    xlabel = "∇q_ice_gm",
    ylabel = "z",
    dpi = 600,
    title = "∇q_ice_gm",
    lw = 2,
    size = (800, 800),
    legend = :topleft,
    color = :black,
    fmt = :png,
    grid = true,
    marker = :circle,
    markerstrokewidth = 0
)
plot!(∇q_ice_gm_c_RB_parent[:], zc_les, label = "∇q_ice_gm_c_RB", lw = 2, color = :red, marker = :circle, markerstrokewidth = 0)
plot!(ylims = ylims)
savefig("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test/∇q_ice_gm.png")



q_ice_gm_c = @. C2Fq_ice(qi_c)
q_ice_gm_c_parent = parent(q_ice_gm_c)
q_ice_gm_c_LB = @. LBq_ice(qi_c)
q_ice_gm_c_LB_parent = parent(q_ice_gm_c_LB)
q_ice_gm_c_RB = @. RBq_ice(qi_c)
q_ice_gm_c_RB_parent = parent(q_ice_gm_c_RB)

q_ice_gm_c_loop = @. F2Cq_ice(C2Fq_ice(qi_c))
q_ice_gm_c_loop_parent = parent(q_ice_gm_c_loop)

plot(q_ice_gm_c_parent[:],
    zf_les,
    label = "q_ice_gm_c",
    xlabel = "q_ice_gm",
    ylabel = "z",
    dpi = 600,
    title = "q_ice_gm",
    lw = 2,
    size = (800, 800),
    legend = :topleft,
    color = :black,
    fmt = :png,
    grid = true,
    marker = :circle,
    markerstrokewidth = 0
)
plot!(q_ice_gm_c_LB_parent[:], zf_les, label = "q_ice_gm_c_LB", lw = 2, color = :red, marker = :circle, markerstrokewidth = 0)
plot!(q_ice_gm_c_RB_parent[:], zf_les, label = "q_ice_gm_c_RB", lw = 2, color = :blue, marker = :circle, markerstrokewidth = 0)

plot!(q_ice_gm_c_loop_parent[:], zc_les, label = "q_ice_gm_c_loop (note not conserved round trip bc peaks were cut)", lw = 2, color = :green, marker = :circle, markerstrokewidth = 0, linestyle = :dash)
plot!(ylims = ylims)
savefig("~/Research_Schneider/CliMA/SOCRATESSingleColumnForcings.jl/test/q_ice_gm_c.png")