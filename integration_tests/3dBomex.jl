if !haskey(ENV, "BUILDKITE")
    import Pkg
    Pkg.develop(Pkg.PackageSpec(; path = dirname(@__DIR__)))
end
using Test
using LinearAlgebra, StaticArrays

include(joinpath(@__DIR__, "..", "driver", "parameter_set.jl"))

import TurbulenceConvection
const TC = TurbulenceConvection

import ClimaCore
const CC = ClimaCore
const CCG = CC.Geometry
const CCO = CC.Operators
import ClimaCore.Geometry: ⊗

import UnPack

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

# set up function space

function hvspace_3D(namelist; xlim = (-5000, 5000), ylim = (-5000, 5000), xelem = 2, yelem = 2, npoly = 4)
    FT = Float64

    zelem = namelist["grid"]["nz"]
    Δz = namelist["grid"]["dz"]
    zlim = (0, zelem * Δz)

    xdomain = CC.Domains.IntervalDomain(CCG.XPoint{FT}(xlim[1]), CCG.XPoint{FT}(xlim[2]), periodic = true)
    ydomain = CC.Domains.IntervalDomain(CCG.YPoint{FT}(ylim[1]), CCG.YPoint{FT}(ylim[2]), periodic = true)

    horzdomain = CC.Domains.RectangleDomain(xdomain, ydomain)
    horzmesh = CC.Meshes.RectilinearMesh(horzdomain, xelem, yelem)
    horztopology = CC.Topologies.Topology2D(horzmesh)

    zdomain =
        CC.Domains.IntervalDomain(CCG.ZPoint{FT}(zlim[1]), CCG.ZPoint{FT}(zlim[2]); boundary_names = (:bottom, :top))
    vertmesh = CC.Meshes.IntervalMesh(zdomain, nelems = zelem)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(vertmesh)
    grid = TC.Grid(vertmesh)

    quad = CC.Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = CC.Spaces.SpectralElementSpace2D(horztopology, quad)

    hv_center_space = CC.Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = CC.Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space, grid)
end


# set up 3D domain - doubly periodic box
include(joinpath(@__DIR__, "..", "driver", "dycore.jl"))
include(joinpath(@__DIR__, "..", "driver", "Surface.jl"))
include(joinpath(@__DIR__, "..", "driver", "generate_namelist.jl"))
import .NameList
namelist = NameList.default_namelist("Bomex")
hv_center_space, hv_face_space, grid = hvspace_3D(namelist)

const MSLP = 1e5 # mean sea level pressure
const grav = 9.8 # gravitational constant
const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const γ = 1.4 # heat capacity ratio
const C_p = R_d * γ / (γ - 1) # heat capacity at constant pressure
const C_v = R_d / (γ - 1) # heat capacity at constant volume
const R_m = R_d # moist R, assumed to be dry

function pressure(ρθ)
    if ρθ >= 0
        return MSLP * (R_d * ρθ / MSLP)^γ
    else
        return NaN
    end
end

Φ(z) = grav * z

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section 5a
function init_dry_rising_bubble_3d(x, y, z)
    x_c = 0.0
    y_c = 0.0
    z_c = 350.0
    r_c = 250.0
    θ_b = 300.0
    θ_c = 0.5
    cp_d = C_p
    cv_d = C_v
    p_0 = MSLP
    g = grav

    # auxiliary quantities
    r = sqrt((x - x_c)^2 + (z - z_c)^2)
    θ_p = r < r_c ? 0.5 * θ_c * (1.0 + cospi(r / r_c)) : 0.0 # potential temperature perturbation

    θ = θ_b + θ_p # potential temperature
    π_exn = 1.0 - g * z / cp_d / θ # exner function
    T = π_exn * θ # temperature
    p = p_0 * π_exn^(cp_d / R_d) # pressure
    ρ = p / R_d / T # density
    ρθ = ρ * θ # potential temperature density

    # Horizontal momentum defined on cell centers
    return (ρ = ρ, ρθ = ρθ, ρuₕ = ρ * CCG.UVVector(0.0, 0.0))
end

# initial conditions
coords = CC.Fields.coordinate_field(hv_center_space)
face_coords = CC.Fields.coordinate_field(hv_face_space)

function init_state(grid, coords, face_coords)
    Yc = map(coords) do coord
        bubble = init_dry_rising_bubble_3d(coord.x, coord.y, coord.z)
        bubble
    end
    # Vertical momentum defined on cell faces
    ρw = map(face_coords) do coord
        CCG.WVector(0.0)
    end

    FT = Float64
    N_up = 1
    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)
    cent_prog_fields() = TC.FieldFromNamedTuple(cspace, cent_prognostic_vars(FT, N_up))
    face_prog_fields() = TC.FieldFromNamedTuple(fspace, face_prognostic_vars(FT, N_up))
    Y = CC.Fields.FieldVector(Yc = Yc, ρw = ρw, cent = cent_prog_fields(), face = face_prog_fields())
    return Y
end

function energy(Yc, ρu, z)
    ρ = Yc.ρ
    ρθ = Yc.ρθ
    u = ρu / ρ
    kinetic = ρ * norm(u)^2 / 2
    potential = z * grav * ρ
    internal = C_v * pressure(ρθ) / R_d
    return kinetic + potential + internal
end
function combine_momentum(ρuₕ, ρw)
    CCG.transform(CCG.UVWAxis(), ρuₕ) + CCG.transform(CCG.UVWAxis(), ρw)
end
function center_momentum(Y)
    If2c = CCO.InterpolateF2C()
    combine_momentum.(Y.Yc.ρuₕ, If2c.(Y.ρw))
end
function total_energy(Y)
    ρ = Y.Yc.ρ
    ρu = center_momentum(Y)
    ρθ = Y.Yc.ρθ
    z = CC.Fields.coordinate_field(axes(ρ)).z
    sum(energy.(Y.Yc, ρu, z))
end

include(joinpath(@__DIR__, "..", "driver", "Cases.jl"))
import .Cases

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
    edmf = TC.EDMFModel(namelist)
    return (; edmf, case, grid, param_set, aux = get_aux(grid))
end

function get_gm_cache(Y, coords)
    return (; coords)
end

function ∑tendencies_3d_bomex!(tendencies, prog, cache, t)
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

    ∑tendencies_3d_bomex_gm!(tendencies, prog, cache, t)
    return nothing
end

function ∑tendencies_3d_bomex_gm!(dY, Y, cache, t)
    UnPack.@unpack coords = cache
    ρw = Y.ρw
    Yc = Y.Yc
    dYc = dY.Yc
    dρw = dY.ρw
    ρ = Yc.ρ
    ρuₕ = Yc.ρuₕ
    ρθ = Yc.ρθ
    dρθ = dYc.ρθ
    dρuₕ = dYc.ρuₕ
    dρ = dYc.ρ

    # spectral horizontal operators
    hdiv = CCO.Divergence()
    hgrad = CCO.Gradient()
    hwdiv = CCO.WeakDivergence()
    hwgrad = CCO.WeakGradient()

    # vertical FD operators with BC's
    vdivf2c = CCO.DivergenceF2C(bottom = CCO.SetValue(CCG.WVector(0.0)), top = CCO.SetValue(CCG.WVector(0.0)))
    vvdivc2f =
        CCO.DivergenceC2F(bottom = CCO.SetDivergence(CCG.WVector(0.0)), top = CCO.SetDivergence(CCG.WVector(0.0)))
    uvdivf2c = CCO.DivergenceF2C(
        bottom = CCO.SetValue(CCG.WVector(0.0) ⊗ CCG.UVVector(0.0, 0.0)),
        top = CCO.SetValue(CCG.WVector(0.0) ⊗ CCG.UVVector(0.0, 0.0)),
    )
    If = CCO.InterpolateC2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    Ic = CCO.InterpolateF2C()
    ∂ = CCO.DivergenceF2C(bottom = CCO.SetValue(CCG.WVector(0.0)), top = CCO.SetValue(CCG.WVector(0.0)))
    ∂f = CCO.GradientC2F()
    ∂c = CCO.GradientF2C()
    B = CCO.SetBoundaryOperator(bottom = CCO.SetValue(CCG.WVector(0.0)), top = CCO.SetValue(CCG.WVector(0.0)))

    fcc = CCO.FluxCorrectionC2C(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())
    fcf = CCO.FluxCorrectionF2F(bottom = CCO.Extrapolate(), top = CCO.Extrapolate())

    uₕ = @. ρuₕ / ρ
    w = @. ρw / If(ρ)
    wc = @. Ic(ρw) / ρ
    p = @. pressure(ρθ)
    θ = @. ρθ / ρ
    Yfρ = @. If(ρ)

    ### HYPERVISCOSITY
    # 1) compute hyperviscosity coefficients
    @. dρθ = hdiv(hgrad(θ))
    @. dρuₕ = hdiv(hgrad(uₕ))
    @. dρw = hdiv(hgrad(w))
    CC.Spaces.weighted_dss!(dYc)
    CC.Spaces.weighted_dss!(dρw)

    κ₄ = 100.0 # m^4/s
    @. dρθ = -κ₄ * hdiv(ρ * hgrad(dρθ))
    @. dρuₕ = -κ₄ * hdiv(ρ * hgrad(dρuₕ))
    @. dρw = -κ₄ * hdiv(Yfρ * hgrad(dρw))

    # density
    @. dρ = -∂(ρw)
    @. dρ -= hdiv(ρuₕ)

    # potential temperature
    @. dρθ += -(∂(ρw * If(ρθ / ρ)))
    @. dρθ -= hdiv(uₕ * ρθ)

    # Horizontal momentum
    @. dρuₕ += -uvdivf2c(ρw ⊗ If(uₕ))
    Ih = Ref(CCG.Axis2Tensor((CCG.UVAxis(), CCG.UVAxis()), @SMatrix [1.0 0.0; 0.0 1.0]),)
    @. dρuₕ -= hdiv(ρuₕ ⊗ uₕ + p * Ih)

    # vertical momentum
    z = coords.z
    @. dρw += B(CCG.transform(CCG.WAxis(), -(∂f(p)) - If(ρ) * ∂f(Φ(z))) - vvdivc2f(Ic(ρw ⊗ w)))
    uₕf = @. If(ρuₕ / ρ) # requires boundary conditions
    @. dρw -= hdiv(uₕf ⊗ ρw)

    ### UPWIND FLUX CORRECTION
    upwind_correction = true
    if upwind_correction
        @. dρ += fcc(w, ρ)
        @. dρθ += fcc(w, ρθ)
        @. dρuₕ += fcc(w, ρuₕ)
        @. dρw += fcf(wc, ρw)
    end

    ### DIFFUSION
    κ₂ = 0.0 # m^2/s
    #  1a) horizontal div of horizontal grad of horiz momentun
    @. dρuₕ += hdiv(κ₂ * (ρ * hgrad(ρuₕ / ρ)))
    #  1b) vertical div of vertical grad of horiz momentun
    @. dρuₕ += uvdivf2c(κ₂ * (Yfρ * ∂f(ρuₕ / ρ)))

    #  1c) horizontal div of horizontal grad of vert momentum
    @. dρw += hdiv(κ₂ * (Yfρ * hgrad(ρw / Yfρ)))
    #  1d) vertical div of vertical grad of vert momentun
    @. dρw += vvdivc2f(κ₂ * (Yc.ρ * ∂c(ρw / Yfρ)))

    #  2a) horizontal div of horizontal grad of potential temperature
    @. dρθ += hdiv(κ₂ * (Yc.ρ * hgrad(ρθ / ρ)))
    #  2b) vertical div of vertial grad of potential temperature
    @. dρθ += ∂(κ₂ * (Yfρ * ∂f(ρθ / ρ)))

    CC.Spaces.weighted_dss!(dYc)
    CC.Spaces.weighted_dss!(dρw)
    return dY
end

Y = init_state(grid, coords, face_coords)
energy_0 = total_energy(Y)
mass_0 = sum(Y.Yc.ρ) # Computes ∫ρ∂Ω such that quadrature weighting is accounted for.
theta_0 = sum(Y.Yc.ρθ)

# Solve the ODE
gm_cache = get_gm_cache(Y, coords)
function get_aux(grid; N_up = 1)
    FT = eltype(grid)
    cspace = TC.center_space(grid)
    fspace = TC.face_space(grid)
    aux_cent_fields = TC.FieldFromNamedTuple(cspace, cent_aux_vars(FT, N_up))
    aux_face_fields = TC.FieldFromNamedTuple(fspace, face_aux_vars(FT, N_up))
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    return aux
end
edmf_cache = get_edmf_cache(grid, namelist)
Δt = 0.05
time_end = 2 * Δt
cache = (; gm_cache..., edmf_cache, Δt)

prob = ODE.ODEProblem(∑tendencies_3d_bomex!, Y, (0.0, time_end), cache)
integrator =
    ODE.init(prob, ODE.SSPRK33(), dt = Δt, saveat = 50.0, progress = true, progress_message = (dt, u, p, t) -> t);

if haskey(ENV, "CI_PERF_SKIP_RUN") # for performance analysis
    throw(:exit_profile)
end

sol = @timev OrdinaryDiffEq.solve!(integrator)

ENV["GKSwstype"] = "nul"
import Plots
Plots.GRBackend()

dir = "3dBomex"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

# post-processing
Es = [total_energy(u) for u in sol.u]
Mass = [sum(u.Yc.ρ) for u in sol.u]
Theta = [sum(u.Yc.ρθ) for u in sol.u]

Plots.png(Plots.plot((Es .- energy_0) ./ energy_0), joinpath(path, "energy.png"))
Plots.png(Plots.plot((Mass .- mass_0) ./ mass_0), joinpath(path, "mass.png"))
Plots.png(Plots.plot((Theta .- theta_0) ./ theta_0), joinpath(path, "rho_theta.png"))

function linkfig(figpath, alt = "")
    # buildkite-agent upload figpath
    # link figure in logs if we are running on CI
    if get(ENV, "BUILDKITE", "") == "true"
        artifact_url = "artifact://$figpath"
        print("\033]1338;url='$(artifact_url)';alt='$(alt)'\a\n")
    end
end

linkfig(relpath(joinpath(path, "energy.png"), joinpath(@__DIR__, "../..")), "Total Energy")
linkfig(relpath(joinpath(path, "rho_theta.png"), joinpath(@__DIR__, "../..")), "Potential Temperature")
linkfig(relpath(joinpath(path, "mass.png"), joinpath(@__DIR__, "../..")), "Mass")
