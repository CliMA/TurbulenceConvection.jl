import JSON
import ArgParse
import TurbulenceConvection
import ClimaCore
const CC = ClimaCore

include("parameter_set.jl")
include("Cases.jl")
import .Cases

struct State{P, A, T, D}
    prog::P
    aux::A
    tendencies::T
    diagnostics::D
end

struct Simulation1d
    io_nt::NamedTuple
    grid
    state
    ref_state
    GMV
    Case
    Turb
    TS
    Stats
    param_set
end

# TODO: not quite sure, yet, where this all should live.

#####
##### Auxiliary fields
#####
# Face & Center
aux_vars_ref_state(FT) = (; ref_state = (ρ0 = FT(0), α0 = FT(0), p0 = FT(0)))

# Center only
cent_aux_vars_gm(FT) = ()
cent_aux_vars_edmf(FT, n_up) = (; turbconv = (; bulk = (; area = FT(0), θ_liq_ice = FT(0), q_tot = FT(0))))
cent_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., cent_aux_vars_gm(FT)..., cent_aux_vars_edmf(FT, n_up)...)

# Face only
face_aux_vars_gm(FT) = ()
face_aux_vars_edmf(FT, n_up) = ()
face_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., face_aux_vars_gm(FT)..., face_aux_vars_edmf(FT, n_up)...)

#####
##### Diagnostic fields
#####

# Center only
cent_diagnostic_vars_gm(FT) = ()
cent_diagnostic_vars_edmf(FT, n_up) = ()
cent_diagnostic_vars(FT, n_up) = (; cent_diagnostic_vars_gm(FT)..., cent_diagnostic_vars_edmf(FT, n_up)...)

# Face only
face_diagnostic_vars_gm(FT) = ()
face_diagnostic_vars_edmf(FT, n_up) = ()
face_diagnostic_vars(FT, n_up) = (; face_diagnostic_vars_gm(FT)..., face_diagnostic_vars_edmf(FT, n_up)...)

#####
##### Prognostic fields
#####

# Center only
cent_prognostic_vars(FT, n_up) = (; cent_prognostic_vars_gm(FT)..., cent_prognostic_vars_edmf(FT, n_up)...)
cent_prognostic_vars_gm(FT) = (; U = FT(0), V = FT(0), H = FT(0), QT = FT(0))
cent_prognostic_vars_up(FT) = (; area = FT(0), H = FT(0), QT = FT(0))
cent_prognostic_vars_en(FT) = (; TKE = FT(0), Hvar = FT(0), QTvar = FT(0), HQTcov = FT(0))
cent_prognostic_vars_edmf(FT, n_up) =
    (; turbconv = (; en = cent_prognostic_vars_en(FT), up = ntuple(i -> cent_prognostic_vars_up(FT), n_up)))
# cent_prognostic_vars_edmf(FT, n_up) = (;) # could also use this for empty model

# Face only
face_prognostic_vars(FT, n_up) = (; face_prognostic_vars_edmf(FT, n_up)...)
face_prognostic_vars_up(FT) = (; W = FT(0))
face_prognostic_vars_edmf(FT, n_up) = (; turbconv = (; up = ntuple(i -> face_prognostic_vars_up(FT), n_up)))
# face_prognostic_vars_edmf(FT, n_up) = (;) # could also use this for empty model

function Simulation1d(namelist)
    TC = TurbulenceConvection
    param_set = create_parameter_set(namelist)

    FT = Float64
    grid = TC.Grid(FT(namelist["grid"]["dz"]), namelist["grid"]["nz"])
    Stats = TC.NetCDFIO_Stats(namelist, grid)
    case = Cases.get_case(namelist)
    ref_params = Cases.reference_params(case, grid, param_set, namelist)
    ref_state = TC.ReferenceState(grid, param_set, Stats; ref_params...)

    GMV = TC.GridMeanVariables(namelist, grid, param_set)
    Sur = TC.SurfaceBase(Cases.get_surface_type(case); namelist, ref_params)
    Fo = TC.ForcingBase{Cases.get_forcing_type(case)}()
    Rad = TC.RadiationBase{Cases.get_radiation_type(case)}()

    Case = Cases.CasesBase(case, namelist, grid, param_set, Sur, Fo, Rad)
    Turb = TC.EDMF_PrognosticTKE(namelist, grid, param_set)
    TS = TC.TimeStepping(namelist)

    n_updrafts = Turb.n_updrafts

    cent_prog_fields = TC.FieldFromNamedTuple(TC.center_space(grid), cent_prognostic_vars(FT, n_updrafts))
    face_prog_fields = TC.FieldFromNamedTuple(TC.face_space(grid), face_prognostic_vars(FT, n_updrafts))
    aux_cent_fields = TC.FieldFromNamedTuple(TC.center_space(grid), cent_aux_vars(FT, n_updrafts))
    aux_face_fields = TC.FieldFromNamedTuple(TC.face_space(grid), face_aux_vars(FT, n_updrafts))
    diagnostic_cent_fields = TC.FieldFromNamedTuple(TC.center_space(grid), cent_diagnostic_vars(FT, n_updrafts))
    diagnostic_face_fields = TC.FieldFromNamedTuple(TC.face_space(grid), face_diagnostic_vars(FT, n_updrafts))

    parent(aux_face_fields.ref_state.p0) .= ref_state.p0
    parent(aux_face_fields.ref_state.ρ0) .= ref_state.rho0
    parent(aux_face_fields.ref_state.α0) .= ref_state.alpha0
    parent(aux_cent_fields.ref_state.p0) .= ref_state.p0_half
    parent(aux_cent_fields.ref_state.ρ0) .= ref_state.rho0_half
    parent(aux_cent_fields.ref_state.α0) .= ref_state.alpha0_half

    prog = CC.Fields.FieldVector(cent = cent_prog_fields, face = face_prog_fields)
    tendencies = CC.Fields.FieldVector(cent = cent_prog_fields, face = face_prog_fields)
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    diagnostics = CC.Fields.FieldVector(cent = diagnostic_cent_fields, face = diagnostic_face_fields)

    state = State(prog, aux, tendencies, diagnostics)

    io_nt = (;
        ref_state = TC.io_dictionary_ref_state(state),
        aux = TC.io_dictionary_aux(state),
        diagnostics = TC.io_dictionary_diagnostics(state),
        prog = TC.io_dictionary_state(state),
        tendencies = TC.io_dictionary_tendencies(state),
    )

    return Simulation1d(io_nt, grid, state, ref_state, GMV, Case, Turb, TS, Stats, param_set)
end

function TurbulenceConvection.initialize(sim::Simulation1d, namelist)
    TC = TurbulenceConvection
    state = sim.state
    Cases.initialize_profiles(sim.Case, sim.grid, sim.GMV, TC.center_ref_state(state))
    TC.satadjust(sim.GMV, sim.grid, sim.state)

    Cases.initialize_surface(sim.Case, sim.grid, state, sim.param_set)
    Cases.initialize_forcing(sim.Case, sim.grid, state, sim.GMV, sim.param_set)
    Cases.initialize_radiation(sim.Case, sim.grid, state, sim.GMV, sim.param_set)

    TC.initialize(sim.Turb, sim.grid, state, sim.Case, sim.GMV, sim.TS)

    TC.initialize_io(sim)
    TC.io(sim)
    return
end

function run(sim::Simulation1d)
    TC = TurbulenceConvection
    iter = 0
    grid = sim.grid
    state = sim.state
    TC.open_files(sim.Stats) # #removeVarsHack
    while sim.TS.t <= sim.TS.t_max
        TC.update(sim.Turb, grid, state, sim.GMV, sim.Case, sim.TS)
        TC.update(sim.TS)

        if mod(iter, 100) == 0
            progress = sim.TS.t / sim.TS.t_max
            @show progress
        end

        if mod(sim.TS.t, sim.Stats.frequency) == 0
            # TODO: is this the best location to call diagnostics?
            TC.compute_diagnostics!(sim.Turb, sim.GMV, grid, state, sim.Case, sim.TS)

            # TODO: remove `vars` hack that avoids
            # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
            # opening/closing files every step should be okay. #removeVarsHack
            # TurbulenceConvection.io(sim) # #removeVarsHack
            TC.write_simulation_time(sim.Stats, sim.TS.t) # #removeVarsHack

            TC.io(sim.io_nt.aux, sim.Stats)
            TC.io(sim.io_nt.diagnostics, sim.Stats)
            TC.io(sim.io_nt.prog, sim.Stats)
            TC.io(sim.io_nt.tendencies, sim.Stats)

            TC.io(sim.GMV, grid, state, sim.Stats) # #removeVarsHack
            TC.io(sim.Case, grid, state, sim.Stats) # #removeVarsHack
            TC.io(sim.Turb, grid, state, sim.Stats, sim.TS, sim.param_set) # #removeVarsHack
        end
        iter += 1
    end
    TC.close_files(sim.Stats) # #removeVarsHack
    return
end

function TurbulenceConvection.initialize_io(sim::Simulation1d)
    TC = TurbulenceConvection
    TC.initialize_io(sim.io_nt.ref_state, sim.Stats)
    TC.io(sim.io_nt.ref_state, sim.Stats) # since the reference prog is static

    TC.initialize_io(sim.io_nt.aux, sim.Stats)
    TC.initialize_io(sim.io_nt.diagnostics, sim.Stats)
    TC.initialize_io(sim.io_nt.prog, sim.Stats)
    TC.initialize_io(sim.io_nt.tendencies, sim.Stats)

    # TODO: depricate
    TC.initialize_io(sim.GMV, sim.Stats)
    TC.initialize_io(sim.Case, sim.Stats)
    TC.initialize_io(sim.Turb, sim.Stats)
    return
end

function TurbulenceConvection.io(sim::Simulation1d)
    TC = TurbulenceConvection
    TC.open_files(sim.Stats)
    TC.write_simulation_time(sim.Stats, sim.TS.t)

    TC.io(sim.io_nt.aux, sim.Stats)
    TC.io(sim.io_nt.diagnostics, sim.Stats)
    TC.io(sim.io_nt.prog, sim.Stats)
    TC.io(sim.io_nt.tendencies, sim.Stats)

    # TODO: depricate
    TC.io(sim.GMV, sim.grid, sim.state, sim.Stats)
    TC.io(sim.Case, sim.grid, sim.state, sim.Stats)
    TC.io(sim.Turb, sim.grid, sim.state, sim.Stats, sim.TS, sim.param_set)
    TC.close_files(sim.Stats)
    return
end

function main(namelist; kwargs...)
    main1d(namelist; kwargs...)
end

function main1d(namelist; time_run = false)
    Simulation = Simulation1d(namelist)
    TurbulenceConvection.initialize(Simulation, namelist)
    if time_run
        @time run(Simulation)
    else
        run(Simulation)
    end
    println("The simulation has completed.")
    return Simulation.Stats.path_plus_file
end


function parse_commandline()
    s = ArgParse.ArgParseSettings(; description = "Run case input")

    ArgParse.@add_arg_table! s begin
        "case_name"
        help = "The case name"
        arg_type = String
        required = true
    end

    return ArgParse.parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__

    args = parse_commandline()
    case_name = args["case_name"]

    namelist = open("namelist_" * "$case_name.in", "r") do io
        JSON.parse(io; dicttype = Dict, inttype = Int64)
    end
    main(namelist)
end
