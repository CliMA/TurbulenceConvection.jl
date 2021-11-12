import TurbulenceConvection
const TC = TurbulenceConvection

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
grid = TC.Grid(FT(0.1), 10)
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

