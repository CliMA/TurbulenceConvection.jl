if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end
using Test
using LinearAlgebra


include(joinpath(@__DIR__, "..", "driver", "parameter_set.jl"))

import TurbulenceConvection
const TC = TurbulenceConvection

import ClimaCore
const CC = ClimaCore
const CCG = CC.Geometry
const CCO = CC.Operators

import UnPack

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

# This experiment tests
#     1) hydrostatic and geostrophic balance;
#     2) linear instability.
# - "baroclinic_wave": the defaul simulation, following https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2241.

# Parameters
const R = 6.371229e6 # radius
const grav = 9.80616 # gravitational constant
const Ω = 7.29212e-5 # Earth rotation (radians / sec)
const R_d = 287.0 # R dry (gas constant / mol mass dry air)
const κ = 2 / 7 # kappa
const γ = 1.4 # heat capacity ratio
const cp_d = R_d / κ # heat capacity at constant pressure
const cv_d = cp_d - R_d # heat capacity at constant volume
const p_0 = 1.0e5 # reference pressure
const k = 3
const T_e = 310 # temperature at the equator
const T_p = 240 # temperature at the pole
const T_0 = 0.5 * (T_e + T_p)
const T_tri = 273.16 # triple point temperature
const Γ = 0.005
const A = 1 / Γ
const B = (T_0 - T_p) / T_0 / T_p
const C = 0.5 * (k + 2) * (T_e - T_p) / T_e / T_p
const b = 2
const H = R_d * T_0 / grav
const z_t = 15.0e3
const λ_c = 20.0
const ϕ_c = 40.0
const d_0 = R / 6
const V_p = 1.0

τ_z_1(z) = exp(Γ * z / T_0)
τ_z_2(z) = 1 - 2 * (z / b / H)^2
τ_z_3(z) = exp(-(z / b / H)^2)
τ_1(z) = 1 / T_0 * τ_z_1(z) + B * τ_z_2(z) * τ_z_3(z)
τ_2(z) = C * τ_z_2(z) * τ_z_3(z)
τ_int_1(z) = A * (τ_z_1(z) - 1) + B * z * τ_z_3(z)
τ_int_2(z) = C * z * τ_z_3(z)
F_z(z) = (1 - 3 * (z / z_t)^2 + 2 * (z / z_t)^3) * (z ≤ z_t)
I_T(ϕ) = cosd(ϕ)^k - k / (k + 2) * (cosd(ϕ))^(k + 2)
T(ϕ, z) = (τ_1(z) - τ_2(z) * I_T(ϕ))^(-1)
p(ϕ, z) = p_0 * exp(-grav / R_d * (τ_int_1(z) - τ_int_2(z) * I_T(ϕ)))
r(λ, ϕ) = R * acos(sind(ϕ_c) * sind(ϕ) + cosd(ϕ_c) * cosd(ϕ) * cosd(λ - λ_c))
U(ϕ, z) = grav * k / R * τ_int_2(z) * T(ϕ, z) * (cosd(ϕ)^(k - 1) - cosd(ϕ)^(k + 1))
u(ϕ, z) = -Ω * R * cosd(ϕ) + sqrt((Ω * R * cosd(ϕ))^2 + R * cosd(ϕ) * U(ϕ, z))
v(ϕ, z) = 0.0
c3(λ, ϕ) = cos(π * r(λ, ϕ) / 2 / d_0)^3
s1(λ, ϕ) = sin(π * r(λ, ϕ) / 2 / d_0)
cond(λ, ϕ) = (0 < r(λ, ϕ) < d_0) * (r(λ, ϕ) != R * pi)
δu(λ, ϕ, z) =
    -16 * V_p / 3 / sqrt(3) *
    F_z(z) *
    c3(λ, ϕ) *
    s1(λ, ϕ) *
    (-sind(ϕ_c) * cosd(ϕ) + cosd(ϕ_c) * sind(ϕ) * cosd(λ - λ_c)) / sin(r(λ, ϕ) / R) * cond(λ, ϕ)
δv(λ, ϕ, z) =
    16 * V_p / 3 / sqrt(3) * F_z(z) * c3(λ, ϕ) * s1(λ, ϕ) * cosd(ϕ_c) * sind(λ - λ_c) / sin(r(λ, ϕ) / R) * cond(λ, ϕ)
const κ₄ = 1.0e17 # m^4/s
uu(λ, ϕ, z) = u(ϕ, z) + δu(λ, ϕ, z)
uv(λ, ϕ, z) = v(ϕ, z) + δv(λ, ϕ, z)

# set up function space
function sphere_3D(::Type{FT}; R = 6.37122e6, helem = 4, npoly = 4, namelist) where {FT}
    zelem = namelist["grid"]["nz"]
    Δz = namelist["grid"]["dz"]
    zlim = (0, zelem * Δz)
    vertdomain =
        CC.Domains.IntervalDomain(CCG.ZPoint{FT}(zlim[1]), CCG.ZPoint{FT}(zlim[2]); boundary_names = (:bottom, :top))
    vertmesh = CC.Meshes.IntervalMesh(vertdomain, nelems = zelem)
    grid = TC.Grid(vertmesh)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = CC.Domains.SphereDomain(R)
    horzmesh = CC.Meshes.EquiangularCubedSphere(horzdomain, helem)
    horztopology = CC.Topologies.Topology2D(horzmesh)
    quad = CC.Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = CC.Spaces.SpectralElementSpace2D(horztopology, quad)

    hv_center_space = CC.Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = CC.Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space, grid)
end

Φ(z) = grav * z

function pressure(ρ, e, normuvw, z)
    I = e - Φ(z) - normuvw^2 / 2
    T = I / cv_d + T_tri
    return ρ * R_d * T
end

# set up 3D domain - spherical shell
include(joinpath(@__DIR__, "..", "driver", "generate_namelist.jl"))
import .NameList
namelist = NameList.default_namelist("Bomex")
hv_center_space, hv_face_space, grid = sphere_3D(Float64; R, helem = 4, npoly = 4, namelist)

# initial conditions
coords = CC.Fields.coordinate_field(hv_center_space)
local_geometries = CC.Fields.local_geometry_field(hv_center_space)
face_coords = CC.Fields.coordinate_field(hv_face_space)

function initial_condition(ϕ, λ, z, edmf_vars)
    ρ = p(ϕ, z) / R_d / T(ϕ, z)
    e = cv_d * (T(ϕ, z) - T_tri) + Φ(z) + (uu(λ, ϕ, z)^2 + uv(λ, ϕ, z)^2) / 2
    ρe = ρ * e

    return (; ρ = ρ, ρe = ρe, edmf_vars...)
end

function initial_condition_velocity(local_geometry)
    coord = local_geometry.coordinates
    ϕ = coord.lat
    λ = coord.long
    z = coord.z
    return CCG.transform(CCG.Covariant12Axis(), CCG.UVVector(uu(λ, ϕ, z), uv(λ, ϕ, z)), local_geometry)
end
# Coriolis
const f = @. CCG.Contravariant3Vector(CCG.WVector(2 * Ω * sind(coords.lat)))
const If2c = CCO.InterpolateF2C()
const Ic2f = CCO.InterpolateC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())

function init_state(edmf, grid, coords, face_coords, local_geometries)
    FT = eltype(edmf)
    edmf_vars = cent_prognostic_vars(FT, edmf)
    Yc = map(coords) do coord
        initial_condition(coord.lat, coord.long, coord.z, edmf_vars)
    end
    uₕ = map(local_geometries) do local_geometry
        initial_condition_velocity(local_geometry)
    end
    w = map(_ -> CCG.Covariant3Vector(0.0), face_coords)

    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)
    cent_prog_fields() = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars(FT, edmf))
    face_prog_fields() = TC.FieldFromNamedTuple(fspace, face_prognostic_vars(FT, edmf))
    Y = CC.Fields.FieldVector(Yc = Yc, uₕ = uₕ, w = w, cent = cent_prog_fields(), face = face_prog_fields())
    return Y
end

include(joinpath(@__DIR__, "..", "driver", "Cases.jl"))
import .Cases

function get_aux(edmf, grid)
    FT = eltype(grid)
    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT, edmf))
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT, edmf))
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    return aux
end

function get_edmf_cache(grid, namelist)
    param_set = create_parameter_set(namelist)
    Ri_bulk_crit = namelist["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"]
    case_type = Cases.get_case(namelist)
    Fo = TC.ForcingBase(case_type, param_set)
    Rad = TC.RadiationBase(case_type)
    surf_ref_state = Cases.surface_ref_state(case_type, param_set, namelist)
    surf_params = Cases.surface_params(case_type, grid, surf_ref_state, param_set; Ri_bulk_crit)
    inversion_type = Cases.inversion_type(case_type)
    case = Cases.CasesBase(case_type; inversion_type, surf_params, Fo, Rad)
    precip_name = TC.parse_namelist(
        namelist,
        "microphysics",
        "precipitation_model";
        default = "None",
        valid_options = ["None", "cutoff", "clima_1m"],
    )
    precip_model = if precip_name == "None"
        TC.NoPrecipitation()
    elseif precip_name == "cutoff"
        TC.CutoffPrecipitation()
    elseif precip_name == "clima_1m"
        TC.Clima1M()
    else
        error("Invalid precip_name $(precip_name)")
    end
    edmf = TC.EDMFModel(namelist, precip_model)
    return (; edmf, case, grid, param_set, aux = get_aux(edmf, grid))
end

function get_gm_cache(Y, coords)
    cuvw = CCG.Covariant123Vector.(Y.uₕ)
    cω³ = CCO.Curl().(Y.uₕ)
    fω¹² = CCO.Curl().(Y.w)
    fu¹² = CCG.Contravariant12Vector.(CCG.Covariant123Vector.(Ic2f.(Y.uₕ)))
    fu³ = CCG.Contravariant3Vector.(CCG.Covariant123Vector.(Y.w))
    cp = @. pressure(Y.Yc.ρ, Y.Yc.ρe / Y.Yc.ρ, norm(cuvw), coords.z)
    cE = @. (norm(cuvw)^2) / 2 + Φ(coords.z)
    return (; cuvw, cω³, fω¹², fu¹², fu³, cp, cE, coords)
end

function ∑tendencies_bw!(tendencies::FV, prog::FV, cache, t::Real) where {FV <: CC.Fields.FieldVector}
    UnPack.@unpack edmf_cache, Δt = cache
    UnPack.@unpack edmf, grid, param_set, aux, case = edmf_cache

    state = TC.State(prog, aux, tendencies)

    surf = get_surface(case.surf_params, grid, state, t, param_set)
    force = case.Fo
    radiation = case.Rad

    TC.affect_filter!(edmf, grid, state, param_set, surf, case.casename, t)

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.
    Cases.update_forcing(case, grid, state, t, param_set)
    Cases.update_radiation(case.Rad, grid, state, param_set)

    TC.update_aux!(edmf, grid, state, surf, param_set, t, Δt)

    tends_face = tendencies.face
    tends_cent = tendencies.cent
    parent(tends_face) .= 0
    parent(tends_cent) .= 0
    # compute tendencies
    TC.compute_turbconv_tendencies!(edmf, grid, state, param_set, surf, Δt)
    ∑tendencies_gm!(tendencies, prog, cache, t)

    return nothing
end
function ∑tendencies_gm!(dY, Y, cache, t)
    UnPack.@unpack cuvw, cω³, fω¹², fu¹², fu³, cp, cE, coords = cache
    cρ = Y.Yc.ρ # scalar on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant12Vector on centers
    cρe = Y.Yc.ρe # scalar on centers

    dρ = dY.Yc.ρ
    dw = dY.w
    duₕ = dY.uₕ
    dρe = dY.Yc.ρe
    z = coords.z

    # 0) update w at the bottom
    # fw = -g^31 cuₕ/ g^33

    hdiv = CCO.Divergence()
    hwdiv = CCO.WeakDivergence()
    hgrad = CCO.Gradient()
    hwgrad = CCO.WeakGradient()
    hcurl = CCO.Curl()
    hwcurl = CCO.WeakCurl()

    @. dρ = 0 * cρ

    ### HYPERVISCOSITY
    # 1) compute hyperviscosity coefficients

    @. dρe = hwdiv(hgrad(cρe / cρ))
    χe = dρe
    @. duₕ = hwgrad(hdiv(cuₕ)) - CCG.Covariant12Vector(hwcurl(CCG.Covariant3Vector(hcurl(cuₕ))),)
    χuₕ = duₕ

    CC.Spaces.weighted_dss!(dρe)
    CC.Spaces.weighted_dss!(duₕ)

    @. dρe = -κ₄ * hwdiv(cρ * hgrad(χe))
    @. duₕ = -κ₄ * (hwgrad(hdiv(χuₕ)) - CCG.Covariant12Vector(hwcurl(CCG.Covariant3Vector(hcurl(χuₕ))),))

    # 1) Mass conservation
    @. cuvw = CCG.Covariant123Vector(cuₕ) + CCG.Covariant123Vector(If2c(fw))

    @. dw = fw * 0

    # 1.a) horizontal divergence
    @. dρ -= hdiv(cρ * (cuvw))

    # 1.b) vertical divergence
    vdivf2c = CCO.DivergenceF2C(
        top = CCO.SetValue(CCG.Contravariant3Vector(0.0)),
        bottom = CCO.SetValue(CCG.Contravariant3Vector(0.0)),
    )
    # we want the total u³ at the boundary to be zero: we can either constrain
    # both to be zero, or allow one to be non-zero and set the other to be its
    # negation

    # explicit part
    @. dρ -= vdivf2c(Ic2f(cρ * cuₕ))
    # implicit part
    @. dρ -= vdivf2c(Ic2f(cρ) * fw)

    # 2) Momentum equation

    # curl term
    # effectively a homogeneous Dirichlet condition on u₁ at the boundary
    vcurlc2f = CCO.CurlC2F(
        bottom = CCO.SetCurl(CCG.Contravariant12Vector(0.0, 0.0)),
        top = CCO.SetCurl(CCG.Contravariant12Vector(0.0, 0.0)),
    )
    @. cω³ = hcurl(cuₕ) # Contravariant3Vector
    @. fω¹² = hcurl(fw) # Contravariant12Vector
    @. fω¹² += vcurlc2f(cuₕ) # Contravariant12Vector

    # cross product
    # convert to contravariant
    # these will need to be modified with topography
    @. fu¹² = CCG.Contravariant12Vector(CCG.Covariant123Vector(Ic2f(cuₕ))) # Contravariant12Vector in 3D
    @. fu³ = CCG.Contravariant3Vector(CCG.Covariant123Vector(fw))
    @. dw -= fω¹² × fu¹² # Covariant3Vector on faces
    @. duₕ -= If2c(fω¹² × fu³)

    # Needed for 3D:
    @. duₕ -= (f + cω³) × CCG.Contravariant12Vector(CCG.Covariant123Vector(cuₕ))

    @. cp = pressure(cρ, cρe / cρ, norm(cuvw), z)

    @. duₕ -= hgrad(cp) / cρ
    vgradc2f = CCO.GradientC2F(
        bottom = CCO.SetGradient(CCG.Covariant3Vector(0.0)),
        top = CCO.SetGradient(CCG.Covariant3Vector(0.0)),
    )
    @. dw -= vgradc2f(cp) / Ic2f(cρ)

    @. cE = (norm(cuvw)^2) / 2 + Φ(z)
    @. duₕ -= hgrad(cE)
    @. dw -= vgradc2f(cE)

    # 3) total energy

    @. dρe -= hdiv(cuvw * (cρe + cp))
    @. dρe -= vdivf2c(fw * Ic2f(cρe + cp))
    @. dρe -= vdivf2c(Ic2f(cuₕ * (cρe + cp)))

    CC.Spaces.weighted_dss!(dY.Yc)
    CC.Spaces.weighted_dss!(dY.uₕ)
    CC.Spaces.weighted_dss!(dY.w)

    return dY
end

include(joinpath(@__DIR__, "..", "driver", "dycore.jl"))
include(joinpath(@__DIR__, "..", "driver", "Surface.jl"))

edmf_cache = get_edmf_cache(grid, namelist)
Y = init_state(edmf_cache.edmf, grid, coords, face_coords, local_geometries)
# Solve the ODE
gm_cache = get_gm_cache(Y, coords)

dt = 5
cache = (; gm_cache..., edmf_cache, Δt = dt)
time_end = 10.0
prob = ODE.ODEProblem(∑tendencies_bw!, Y, (0.0, time_end), cache)

integrator = ODE.init(
    prob,
    ODE.SSPRK33(),
    dt = dt,
    saveat = 200 * dt,
    progress = true,
    adaptive = false,
    progress_message = (dt, u, p, t) -> t,
)

if haskey(ENV, "CI_PERF_SKIP_RUN") # for performance analysis
    throw(:exit_profile)
end

sol = @timev ODE.solve!(integrator)

# # The latest ClimaCore has compat issues with ClimaCorePlots
# # so we've temporarily remove our dependence on ClimaCorePlots
# # so that we can update the latest ClimaCore version.

# import ClimaCorePlots, Plots
# ENV["GKSwstype"] = "nul"
# # visualization artifacts

# @info "Solution L₂ norm at time t = 0: ", norm(Y.Yc.ρe)
# @info "Solution L₂ norm at time t = $(time_end): ", norm(sol.u[end].Yc.ρe)

# anim = Plots.@animate for sol1 in sol.u
#     uₕ = sol1.uₕ
#     uₕ_phy = CCG.transform.(Ref(CCG.UVAxis()), uₕ)
#     v = uₕ_phy.components.data.:2
#     Plots.plot(v, level = 3, clim = (-6, 6))
# end

# dir = "baroclinic_wave"
# path = joinpath(@__DIR__, "output", dir)
# mkpath(path)
# Plots.mp4(anim, joinpath(path, "v.mp4"), fps = 5)
