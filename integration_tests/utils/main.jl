import JSON
using ArgParse
import TurbulenceConvection
include("Cases.jl")
import .Cases

mutable struct Simulation1d
    Gr
    Ref
    GMV
    Case
    Turb
    TS
    Stats
end

function Simulation1d(namelist, paramlist, inpath=pwd())
    Gr = TurbulenceConvection.Grid(namelist)
    Ref = TurbulenceConvection.ReferenceState(Gr)
    GMV = TurbulenceConvection.GridMeanVariables(namelist, Gr, Ref)
    Case = Cases.CasesFactory(namelist, paramlist, Gr, Ref)
    Turb = TurbulenceConvection.ParameterizationFactory(namelist,paramlist, Gr, Ref)
    TS = TurbulenceConvection.TimeStepping(namelist)
    Stats = TurbulenceConvection.NetCDFIO_Stats(namelist, paramlist, Gr, inpath)
    return Simulation1d(Gr,Ref,GMV,Case,Turb,TS,Stats)
end

function TurbulenceConvection.initialize(self::Simulation1d, namelist)
    Cases.initialize_reference(self.Case, self.Gr, self.Ref, self.Stats)
    Cases.initialize_profiles(self.Case, self.Gr, self.GMV, self.Ref)
    Cases.initialize_surface(self.Case, self.Gr, self.Ref )
    Cases.initialize_forcing(self.Case, self.Gr, self.Ref, self.GMV)
    Cases.initialize_radiation(self.Case, self.Gr, self.Ref, self.GMV)
    TurbulenceConvection.initialize(self.Turb, self.Case, self.GMV, self.Ref)
    TurbulenceConvection.initialize_io(self)
    TurbulenceConvection.io(self)

    return
end

function run(self::Simulation1d)
    iter = 0
    TurbulenceConvection.open_files(self.Stats) # #removeVarsHack
    while self.TS.t <= self.TS.t_max
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "1")
        TurbulenceConvection.zero_tendencies(self.GMV)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "2")
        Cases.update_surface(self.Case, self.GMV, self.TS)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "3")
        Cases.update_forcing(self.Case, self.GMV, self.TS)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "4")
        Cases.update_radiation(self.Case, self.GMV, self.TS)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "5")
        TurbulenceConvection.update(self.Turb, self.GMV, self.Case, self.TS)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "6")
        TurbulenceConvection.update(self.TS)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "7")
        # Apply the tendencies, also update the BCs and diagnostic thermodynamics
        TurbulenceConvection.update(self.GMV, self.TS)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "8")
        TurbulenceConvection.update_GMV_diagnostics(self.Turb, self.GMV)
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "9")

        if mod(iter, 100) == 0
            progress = self.TS.t/self.TS.t_max
            @show progress
        end
        if mod(self.TS.t, self.Stats.frequency) == 0
            # TODO: remove `vars` hack that avoids
            # https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
            # opening/closing files every step should be okay. #removeVarsHack
            # TurbulenceConvection.io(self) # #removeVarsHack
            TurbulenceConvection.write_simulation_time(self.Stats, self.TS.t) # #removeVarsHack
            TurbulenceConvection.io(self.GMV, self.Stats) # #removeVarsHack
            TurbulenceConvection.io(self.Case, self.Stats) # #removeVarsHack
            TurbulenceConvection.io(self.Turb, self.Stats, self.TS) # #removeVarsHack
        end
        TurbulenceConvection.check_nans(self.GMV, self.Turb, "9")
        iter+=1
    end
    TurbulenceConvection.close_files(self.Stats) # #removeVarsHack
    return
end

function TurbulenceConvection.initialize_io(self::Simulation1d)

    TurbulenceConvection.initialize_io(self.GMV, self.Stats)
    TurbulenceConvection.initialize_io(self.Case, self.Stats)
    TurbulenceConvection.initialize_io(self.Turb, self.Stats)
    return
end

function TurbulenceConvection.io(self::Simulation1d)
    TurbulenceConvection.open_files(self.Stats)
    TurbulenceConvection.write_simulation_time(self.Stats, self.TS.t)
    TurbulenceConvection.io(self.GMV, self.Stats)
    TurbulenceConvection.io(self.Case, self.Stats)
    TurbulenceConvection.io(self.Turb, self.Stats, self.TS)
    TurbulenceConvection.close_files(self.Stats)
    return
end

function force_io(self::Simulation1d)
    return
end

function main(namelist, paramlist, inpath=pwd(); kwargs...)
    main1d(namelist, paramlist, inpath; kwargs...)
end

function main1d(namelist, paramlist, inpath=pwd(); time_run = false)
    Simulation = Simulation1d(namelist, paramlist, inpath)
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
    s = ArgParseSettings(;description="Run case input")

    @add_arg_table! s begin
        "case_name"
            help = "The case name"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__

    args = parse_commandline()
    case_name = args["case_name"]

    namelist = open("$case_name.in", "r") do io
       JSON.parse(io; dicttype=Dict, inttype=Int64)
    end
    paramlist = open("paramlist_"*"$case_name.in", "r") do io
       JSON.parse(io; dicttype=Dict, inttype=Int64)
    end
    main(namelist, paramlist)
end
