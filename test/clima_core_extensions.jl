import TurbulenceConvection
const TC = TurbulenceConvection

import Dierckx

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

import ClimaCore
const CC = ClimaCore

prognostic_vars(FT, n_up) = (; prognostic_vars_gm(FT)..., prognostic_vars_edmf(FT, n_up)...)
prognostic_vars_gm(FT) = (; U = FT(0), V = FT(0), H = FT(0), QT = FT(0))
prognostic_vars_up(FT) = (; Area = FT(0), H = FT(0), QT = FT(0))
prognostic_vars_en(FT) = (; TKE = FT(0), Hvar = FT(0), QTvar = FT(0), HQTcov = FT(0))
prognostic_vars_edmf(FT, n_up) =
    (; turbconv = (; en = prognostic_vars_en(FT), up = ntuple(i -> prognostic_vars_up(FT), n_up)))
center_space(grid::TC.Grid) = grid.cs
face_space(grid::TC.Grid) = grid.fs

FT = Float64
n_updrafts = 2
n_cells = 100
Δz = 1 / FT(n_cells)
grid = TC.Grid(Δz, n_cells)
cent_state = TC.FieldFromNamedTuple(center_space(grid), prognostic_vars(FT, n_updrafts))
face_state = TC.FieldFromNamedTuple(face_space(grid), prognostic_vars(FT, n_updrafts))
state = CC.Fields.FieldVector(cent = cent_state, face = face_state)

using Test
for k in TC.real_center_indices(grid)
    for i in 1:n_updrafts
        x = state.cent.turbconv.up[i].Area[k]
        state.cent.turbconv.up[i].Area[k] = 2

        @test_throws ErrorException x = state.face.turbconv.up[i].Area[k]
        @test_throws ErrorException state.face.turbconv.up[i].Area[k] = 2
    end
end

for k in TC.real_face_indices(grid)
    for i in 1:n_updrafts
        x = state.face.turbconv.up[i].Area[k]
        state.face.turbconv.up[i].Area[k] = 2

        @test_throws ErrorException x = state.cent.turbconv.up[i].Area[k]
        @test_throws ErrorException state.cent.turbconv.up[i].Area[k] = 2
    end
end

function ∫field(field::CC.Fields.FiniteDifferenceField, field_zmin = 0)
    space = axes(field)
    # TODO: use API and less internals
    topology = space.topology
    domain = topology.mesh.domain
    zmin = domain.coord_min
    zmax = domain.coord_max
    z_span = (zmin.z, zmax.z)
    ld = CC.Spaces.local_geometry_data(space)
    z = ld.coordinates
    field_z = Dierckx.Spline1D(vec(parent(z)), vec(field); k = 1)
    integrand(x, params, z) = field_z(z)
    # Assumes uniform grid spacing:
    prob = ODE.ODEProblem(integrand, field_zmin, z_span; dt = grid.Δz)
    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-12, abstol = 1e-12)
    return sol
end


@testset "Vertical integrals" begin
    cent_state = TC.FieldFromNamedTuple(center_space(grid), (; y = FT(0)))
    face_state = TC.FieldFromNamedTuple(face_space(grid), (; y = FT(0)))
    state = CC.Fields.FieldVector(cent = cent_state, face = face_state)
    yc = state.cent.y
    yf = state.face.y
    cents = CC.Fields.coordinate_field(center_space(grid))
    faces = CC.Fields.coordinate_field(face_space(grid))

    yc .= 1 .+ sin.(cents.z)
    yf .= 1 .+ sin.(faces.z)

    ∫yc = ∫field(yc)
    ∫yf = ∫field(yf)

    boundary_fact(k) = k.i == 1 || k.i == n_cells + 1 ? 0.5 : 1
    ∫yc_simple = sum(k -> yc[k] * grid.Δz, TC.real_center_indices(grid))
    ∫yf_simple = sum(k -> boundary_fact(k) * yf[k] * grid.Δz, TC.real_face_indices(grid))

    ∫yc_new = sum(yc)
    ∫yf_new = sum(yf)

    ∫yc_analytic = 1 - cos(1) - (0 - cos(0))
    ∫yf_analytic = 1 - cos(1) - (0 - cos(0))

    @testset "Vertical definite integrals" begin
        @test ∫yc_simple ≈ ∫yc_analytic atol = grid.Δz
        @test ∫yf_simple ≈ ∫yf_analytic atol = grid.Δz

        @test ∫yc[end] ≈ ∫yc_analytic atol = grid.Δz
        @test ∫yf[end] ≈ ∫yf_analytic atol = grid.Δz

        @test ∫yc_new ≈ ∫yc_analytic atol = grid.Δz
        @test ∫yf_new ≈ ∫yf_analytic atol = grid.Δz
    end

    ∫yc_simple = copy(yc)
    ∫yf_simple = copy(yf)
    ∫yc_simple .= parent(yc)[1]
    ∫yf_simple .= parent(yf)[1]
    ∫yc_simple .= 0
    ∫yf_simple .= 0

    for k in TC.real_center_indices(grid)
        ∫yc_simple[k] += yc[k] * grid.Δz
        if k.i > 1
            ∫yc_simple[k] += ∫yc_simple[k - 1]
        end
    end
    for k in TC.real_face_indices(grid)
        ∫yf_simple[k] += yf[k] * grid.Δz * boundary_fact(k)
        if k.i > 1
            ∫yf_simple[k] += ∫yf_simple[k - 1]
        end
    end

    ∫yc_analytic = cents.z .- cos.(cents.z) .- (0 .- cos(0))
    ∫yf_analytic = faces.z .- cos.(faces.z) .- (0 .- cos(0))

    ∫yc_new = cumsum(yc)
    ∫yf_new = cumsum(yf)

    @testset "Vertical indefinite integrals" begin
        @test all(isapprox.(parent(∫yc_simple), parent(∫yc_analytic), atol = grid.Δz))
        @test all(isapprox.(parent(∫yf_simple), parent(∫yf_analytic), atol = grid.Δz))

        @test all(isapprox.(parent(∫yc.(cents.z)), parent(∫yc_analytic), atol = grid.Δz))
        @test all(isapprox.(parent(∫yf.(faces.z)), parent(∫yf_analytic), atol = grid.Δz))

        @test all(isapprox.(parent(∫yc_new), parent(∫yc_analytic), atol = grid.Δz))
        @test all(isapprox.(parent(∫yf_new), parent(∫yf_analytic), atol = grid.Δz))
    end
end

nothing
