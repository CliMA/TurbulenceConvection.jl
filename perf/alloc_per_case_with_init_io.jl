# Launch with `julia --project --track-allocation=user`
include(joinpath(@__DIR__, "common.jl"))
import Profile

case_name = ENV["ALLOCATION_CASE_NAME"]
@info "Recording allocations for `$case_name`"
namelist = NameList.default_namelist(case_name)
namelist["time_stepping"]["t_max"] = 10800.0 # run for shorter time

main(namelist) # compile first
Profile.clear_malloc_data()
main(namelist)
