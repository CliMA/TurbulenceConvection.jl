module NameList
#Adapated from PyCLES: https://github.com/pressel/pycles
export default_namelist
function default_namelist(case_name)
    namelist_defaults = Dict()
    namelist_defaults["grid"] = Dict()
    namelist_defaults["grid"]["dims"] = 1
    namelist_defaults["grid"]["gw"] = 2

    namelist_defaults["thermodynamics"] = Dict()
    namelist_defaults["thermodynamics"]["thermal_variable"] = "thetal"
    namelist_defaults["thermodynamics"]["sgs"] = "quadrature"
    namelist_defaults["thermodynamics"]["quadrature_order"] = 3
    namelist_defaults["thermodynamics"]["quadrature_type"] = "log-normal" #"gaussian" or "log-normal"

    namelist_defaults["time_stepping"] = Dict()

    namelist_defaults["microphysics"] = Dict()
    namelist_defaults["microphysics"]["rain_model"] = "None"

    namelist_defaults["turbulence"] = Dict()
    namelist_defaults["turbulence"]["scheme"] = "EDMF_PrognosticTKE"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"] = Dict()
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = 1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "moisture_deficit"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["extrapolate_buoyancy"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_steady_updrafts"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_local_micro"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_constant_plume_spacing"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["use_similarity_diffusivity"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["constant_area"] = false
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["calculate_tke"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["calc_scalar_var"] = true
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["mixing_length"] = "sbtd_eq"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_asp_label"] = "const"

    namelist_defaults["output"] = Dict()
    namelist_defaults["output"]["output_root"] = "./"

    namelist_defaults["stats_io"] = Dict()
    namelist_defaults["stats_io"]["stats_dir"] = "stats"
    namelist_defaults["stats_io"]["frequency"] = 60.0

    namelist_defaults["meta"] = Dict()

    if case_name == "Bomex"
        namelist = Bomex(namelist_defaults)
    elseif case_name == "Nieuwstadt"
        namelist = Nieuwstadt(namelist_defaults)
    elseif case_name == "life_cycle_Tan2018"
        namelist = life_cycle_Tan2018(namelist_defaults)
    elseif case_name == "Soares"
        namelist = Soares(namelist_defaults)
    elseif case_name == "Rico"
        namelist = Rico(namelist_defaults)
    elseif case_name == "TRMM_LBA"
        namelist = TRMM_LBA(namelist_defaults)
    elseif case_name == "ARM_SGP"
        namelist = ARM_SGP(namelist_defaults)
    elseif case_name == "GATE_III"
        namelist = GATE_III(namelist_defaults)
    elseif case_name == "DYCOMS_RF01"
        namelist = DYCOMS_RF01(namelist_defaults)
    elseif case_name == "GABLS"
        namelist = GABLS(namelist_defaults)
    elseif case_name == "SP"
        namelist = SP(namelist_defaults)
    elseif case_name == "DryBubble"
        namelist = DryBubble(namelist_defaults)
    else
        error("Not a valid case name")
    end
    # write_file(namelist)
end


function Soares(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["dt"] = 30.0
    namelist["time_stepping"]["t_max"] = 8 * 3600.0

    namelist["meta"]["simname"] = "Soares"
    namelist["meta"]["casename"] = "Soares"

    return namelist
end

function Nieuwstadt(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["dt"] = 10.0
    namelist["time_stepping"]["t_max"] = 8 * 3600.0

    namelist["meta"]["simname"] = "Nieuwstadt"
    namelist["meta"]["casename"] = "Nieuwstadt"

    return namelist
end

function Bomex(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 60
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["dt"] = 20.0
    namelist["time_stepping"]["t_max"] = 21600.0

    namelist["meta"]["simname"] = "Bomex"
    namelist["meta"]["casename"] = "Bomex"

    return namelist
end

function life_cycle_Tan2018(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 75
    namelist["grid"]["dz"] = 40.0

    namelist["time_stepping"]["dt"] = 30.0
    namelist["time_stepping"]["t_max"] = 6*3600.0
    namelist["meta"]["simname"] = "life_cycle_Tan2018"
    namelist["meta"]["casename"] = "life_cycle_Tan2018"

    return namelist
end

function Rico(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 120
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["dt"] = 20.0
    namelist["time_stepping"]["t_max"] = 86400.0

    # namelist["microphysics"]["rain_model"] = "cutoff"
    namelist["microphysics"]["rain_model"] = "clima_1m"

    namelist["meta"]["simname"] = "Rico"
    namelist["meta"]["casename"] = "Rico"

    return namelist
end

function TRMM_LBA(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 320
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["dt"] = 30.0
    namelist["time_stepping"]["t_max"] = 21600.0

    namelist["microphysics"]["rain_model"] = "cutoff"
    # namelist["microphysics"]["rain_model"] = "clima_1m"

    namelist["meta"]["simname"] = "TRMM_LBA"
    namelist["meta"]["casename"] = "TRMM_LBA"

    return namelist
end

function ARM_SGP(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 88
    namelist["grid"]["dz"] =  50.0

    namelist["time_stepping"]["dt"] = 10.0
    namelist["time_stepping"]["t_max"] = 3600.0 * 14.5
    namelist["meta"]["simname"] = "ARM_SGP"
    namelist["meta"]["casename"] = "ARM_SGP"

    return namelist
end

function GATE_III(namelist_defaults)

    # adopted from: "Large eddy simulation of Maritime Deep Tropical Convection",
    # By Khairoutdinov et al (2009)  JAMES, vol. 1, article #15
    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 1700
    namelist["grid"]["dz"] = 10

    namelist["time_stepping"]["dt"] = 5.0
    namelist["time_stepping"]["t_max"] = 3600.0 * 24.0
    namelist["meta"]["simname"] = "GATE_III"
    namelist["meta"]["casename"] = "GATE_III"

    return namelist
end

function DYCOMS_RF01(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 30
    namelist["grid"]["dz"] = 50

    namelist["time_stepping"]["dt"] = 10.0
    # namelist["time_stepping"]["t_max"] = 60 * 60 * 16.
    namelist["time_stepping"]["t_max"] = 30.0
    namelist["meta"]["simname"] = "DYCOMS_RF01"
    namelist["meta"]["casename"] = "DYCOMS_RF01"

    return namelist
end

function GABLS(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 8
    namelist["grid"]["dz"] = 50.0

    namelist["time_stepping"]["dt"] = 1.0
    namelist["time_stepping"]["t_max"] = 9 * 3600.0
    namelist["meta"]["simname"] = "GABLS"
    namelist["meta"]["casename"] = "GABLS"

    return namelist
end

# Sullivan Patton not fully implemented - Ignacio
function SP(namelist_defaults)

    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 256
    namelist["grid"]["dz"] = 8

    namelist["time_stepping"]["dt"] = 5.0
    namelist["time_stepping"]["t_max"] = 7200.0
    namelist["meta"]["simname"] = "SP"
    namelist["meta"]["casename"] = "SP"

    return namelist
end

function DryBubble(namelist_defaults)
    namelist = deepcopy(namelist_defaults)

    namelist["grid"]["nz"] = 200
    namelist["grid"]["dz"] = 50.0

    namelist["stats_io"]["frequency"] = 10.0
    namelist["time_stepping"]["dt"] = 10.0
    namelist["time_stepping"]["t_max"] = 1000.0
    namelist["meta"]["simname"] = "DryBubble"
    namelist["meta"]["casename"] = "DryBubble"

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "moisture_deficit_div"

    return namelist
end

# function write_file(namelist)

#     try
#         type(namelist["meta"]["simname"])
#     catch
#         print("Casename not specified in namelist dictionary!")
#         print("FatalError")
#         exit()
#     end

#     namelist["meta"]["uuid"] = str(uuid.uuid4())

#     fh = open(namelist["meta"]["simname"] + ".in", "w")
#     #pprint.pprint(namelist)
#     json.dump(namelist, fh, sort_keys=true, indent=4)
#     fh.close()

#     return
# end

# if __name__ == "__main__"
#     main()

end
