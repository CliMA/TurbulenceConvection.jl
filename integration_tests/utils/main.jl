import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

import UnPack
import JSON
import ArgParse
import TurbulenceConvection

import ClimaCore
const CC = ClimaCore

import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

include("initial_conditions.jl")
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
cent_aux_vars_gm(FT) = (;
    tke = FT(0),
    Hvar = FT(0),
    QTvar = FT(0),
    HQTcov = FT(0),
    q_liq = FT(0),
    q_ice = FT(0),
    RH = FT(0),
    T = FT(0),
    buoy = FT(0),
    cloud_fraction = FT(0),
    H_third_m = FT(0),
    W_third_m = FT(0),
    QT_third_m = FT(0),
)
cent_aux_vars_en_2m(FT) = (;
    dissipation = FT(0),
    shear = FT(0),
    entr_gain = FT(0),
    detr_loss = FT(0),
    press = FT(0),
    buoy = FT(0),
    interdomain = FT(0),
    rain_src = FT(0),
)
cent_aux_vars_up(FT) = (;
    q_liq = FT(0), q_ice = FT(0), T = FT(0), RH = FT(0), buoy = FT(0), area = FT(0), q_tot = FT(0), θ_liq_ice = FT(0),
)
cent_aux_vars_edmf(FT, n_up) = (;
    turbconv = (;
        bulk = (;
            area = FT(0),
            θ_liq_ice = FT(0),
            RH = FT(0),
            buoy = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            T = FT(0),
        ),
        up = ntuple(i -> cent_aux_vars_up(FT), n_up),
        en = (;
            area = FT(0),
            q_tot = FT(0),
            q_liq = FT(0),
            q_ice = FT(0),
            θ_liq_ice = FT(0),
            RH = FT(0),
            T = FT(0),
            buoy = FT(0),
            cloud_fraction = FT(0),
        ),
        en_2m = (;
            tke = cent_aux_vars_en_2m(FT),
            Hvar = cent_aux_vars_en_2m(FT),
            QTvar = cent_aux_vars_en_2m(FT),
            HQTcov = cent_aux_vars_en_2m(FT),
        ),
        KM = FT(0),
        KH = FT(0),
    ),
)
cent_aux_vars(FT, n_up) = (; aux_vars_ref_state(FT)..., cent_aux_vars_gm(FT)..., cent_aux_vars_edmf(FT, n_up)...)

# Face only
face_aux_vars_gm(FT) = ()
face_aux_vars_up(FT) = (; w = FT(0))

face_aux_vars_edmf(FT, n_up) = (;
    turbconv = (;
        bulk = (; w = FT(0)),
        ρ_ae_KM = FT(0),
        ρ_ae_KH = FT(0),
        en = (; w = FT(0)),
        up = ntuple(i -> face_aux_vars_up(FT), n_up),
    ),
)
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
cent_prognostic_vars_gm(FT) = (; u = FT(0), v = FT(0), θ_liq_ice = FT(0), q_tot = FT(0))
cent_prognostic_vars_up(FT) = (; ρarea = FT(0), ρaθ_liq_ice = FT(0), ρaq_tot = FT(0))
cent_prognostic_vars_en(FT) = (; tke = FT(0), Hvar = FT(0), QTvar = FT(0), HQTcov = FT(0))
cent_prognostic_vars_edmf(FT, n_up) = (;
    turbconv = (;
        en = cent_prognostic_vars_en(FT), up = ntuple(i -> cent_prognostic_vars_up(FT), n_up), ra = (; qr = FT(0)),
    ),
)
# cent_prognostic_vars_edmf(FT, n_up) = (;) # could also use this for empty model

# Face only
face_prognostic_vars(FT, n_up) = (; w = FT(0), face_prognostic_vars_edmf(FT, n_up)...)
face_prognostic_vars_up(FT) = (; ρaw = FT(0))
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

    prog = CC.Fields.FieldVector(cent = cent_prog_fields, face = face_prog_fields)
    tendencies = CC.Fields.FieldVector(cent = deepcopy(cent_prog_fields), face = deepcopy(face_prog_fields))
    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    diagnostics = CC.Fields.FieldVector(cent = diagnostic_cent_fields, face = diagnostic_face_fields)

    state = State(prog, aux, tendencies, diagnostics)

    TC.compute_ref_state!(state, grid, param_set, Stats; ref_params...)

    io_nt = (;
        ref_state = TC.io_dictionary_ref_state(state),
        aux = TC.io_dictionary_aux(state),
        diagnostics = TC.io_dictionary_diagnostics(state),
        prog = TC.io_dictionary_state(state),
        tendencies = TC.io_dictionary_tendencies(state),
    )

    return Simulation1d(io_nt, grid, state, GMV, Case, Turb, TS, Stats, param_set)
end

function TurbulenceConvection.initialize(sim::Simulation1d, namelist)
    TC = TurbulenceConvection
    state = sim.state
    Cases.initialize_profiles(sim.Case, sim.grid, sim.GMV, state)
    TC.satadjust(sim.GMV, sim.grid, sim.state)

    Cases.initialize_surface(sim.Case, sim.grid, state, sim.param_set)
    Cases.initialize_forcing(sim.Case, sim.grid, state, sim.GMV, sim.param_set)
    Cases.initialize_radiation(sim.Case, sim.grid, state, sim.GMV, sim.param_set)

    initialize_edmf(sim.Turb, sim.grid, state, sim.Case, sim.GMV, sim.TS)

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

    t_span = (0.0, sim.TS.t_max)
    params = (;
        edmf = sim.Turb,
        grid = grid,
        state = state,
        gm = sim.GMV,
        Case = sim.Case,
        TS = sim.TS,
        Stats = sim.Stats,
        io_nt = sim.io_nt,
    )

    condition_io(u, t, integrator) = mod(round(Int, t), round(Int, sim.Stats.frequency)) == 0
    condition_every_iter(u, t, integrator) = true

    function affect_io!(integrator)
        UnPack.@unpack edmf, io_nt, grid, state, Case, gm, TS, Stats = integrator.p
        parent(state.prog) .= integrator.u
        param_set = TC.parameter_set(gm)
        # TODO: is this the best location to call diagnostics?
        TC.compute_diagnostics!(edmf, gm, grid, state, Case, TS)

        # TODO: remove `vars` hack that avoids
        # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
        # opening/closing files every step should be okay. #removeVarsHack
        # TurbulenceConvection.io(sim) # #removeVarsHack
        TC.write_simulation_time(Stats, TS.t) # #removeVarsHack

        TC.io(io_nt.aux, Stats)
        TC.io(io_nt.diagnostics, Stats)
        TC.io(io_nt.prog, Stats)
        TC.io(io_nt.tendencies, Stats)

        TC.io(gm, grid, state, Stats) # #removeVarsHack
        TC.io(Case, grid, state, Stats) # #removeVarsHack
        TC.io(edmf, grid, state, Stats, TS, param_set) # #removeVarsHack
    end

    function affect_filter!(integrator)
        UnPack.@unpack edmf, grid, state, gm, Case, TS = integrator.p
        parent(state.prog) .= parent(integrator.u)
        prog_en = TC.center_prog_environment(state)
        up = edmf.UpdVar
        ###
        ### Filters
        ###
        TC.set_edmf_surface_bc(edmf, grid, state, up, Case.Sur)
        TC.filter_updraft_vars(edmf, grid, state, gm)

        @inbounds for k in TC.real_center_indices(grid)
            prog_en.tke[k] = max(prog_en.tke[k], 0.0)
            prog_en.Hvar[k] = max(prog_en.Hvar[k], 0.0)
            prog_en.QTvar[k] = max(prog_en.QTvar[k], 0.0)
            prog_en.HQTcov[k] = max(prog_en.HQTcov[k], -sqrt(prog_en.Hvar[k] * prog_en.QTvar[k]))
            prog_en.HQTcov[k] = min(prog_en.HQTcov[k], sqrt(prog_en.Hvar[k] * prog_en.QTvar[k]))
        end
        parent(integrator.u) .= parent(state.prog)
    end

    # TODO: this should not be handled via callbacks
    function affect_implicit!(integrator)
        UnPack.@unpack edmf, grid, state, gm, TS = integrator.p
        parent(state.prog) .= parent(integrator.u)
        implicit_eqs = edmf.implicit_eqs
        # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
        param_set = TC.parameter_set(gm)
        up = edmf.UpdVar
        prog_en = TC.center_prog_environment(state)

        common_args = (
            grid,
            param_set,
            state,
            TS,
            up.n_updrafts,
            edmf.minimum_area,
            edmf.pressure_plume_spacing,
            edmf.frac_turb_entr,
            edmf.entr_sc,
            edmf.mixing_length,
        )

        implicit_eqs.A_TKE .= TC.construct_tridiag_diffusion_en(common_args..., true)
        implicit_eqs.A_Hvar .= TC.construct_tridiag_diffusion_en(common_args..., false)
        implicit_eqs.A_QTvar .= TC.construct_tridiag_diffusion_en(common_args..., false)
        implicit_eqs.A_HQTcov .= TC.construct_tridiag_diffusion_en(common_args..., false)

        implicit_eqs.b_TKE .= TC.en_diffusion_tendencies(grid, state, TS, :tke, up.n_updrafts)
        implicit_eqs.b_Hvar .= TC.en_diffusion_tendencies(grid, state, TS, :Hvar, up.n_updrafts)
        implicit_eqs.b_QTvar .= TC.en_diffusion_tendencies(grid, state, TS, :QTvar, up.n_updrafts)
        implicit_eqs.b_HQTcov .= TC.en_diffusion_tendencies(grid, state, TS, :HQTcov, up.n_updrafts)
        parent(prog_en.tke) .= implicit_eqs.A_TKE \ implicit_eqs.b_TKE
        parent(prog_en.Hvar) .= implicit_eqs.A_Hvar \ implicit_eqs.b_Hvar
        parent(prog_en.QTvar) .= implicit_eqs.A_QTvar \ implicit_eqs.b_QTvar
        parent(prog_en.HQTcov) .= implicit_eqs.A_HQTcov \ implicit_eqs.b_HQTcov

        parent(integrator.u) .= parent(state.prog)

    end
    callback_io = ODE.DiscreteCallback(condition_io, affect_io!)
    callback_filters = ODE.DiscreteCallback(condition_every_iter, affect_filter!)
    callback_implicit_solve = ODE.DiscreteCallback(condition_every_iter, affect_implicit!)

    prob = ODE.ODEProblem(TC.step!, state.prog, t_span, params; dt = sim.TS.dt)
    sol = ODE.solve(
        prob,
        ODE.Euler(),
        reltol = 1e-12,
        abstol = 1e-12,
        callback = ODE.CallbackSet(callback_implicit_solve, callback_filters, callback_io),
        progress = true,
        progress_message = (dt, u, p, t) -> t,
    )

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
