#=
A simple script for updating the manifest
files in all of our environments.

For example,
```
julia --project=.dev .dev/rm_add_up_deps.jl --pkg SurfaceFluxes --ver 0.4.2
```
=#

import ArgParse
function parse_commandline()
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        "--pkgname"
        help = "Package to rm/add"
        arg_type = String
        required = true
        default = ""
        "--ver"
        help = "Version to add"
        arg_type = String
        required = true
        default = ""
    end
    parsed_args = ArgParse.parse_args(ARGS, s)
    return (s, parsed_args)
end
(s, parsed_args) = parse_commandline()
pkg = parsed_args["pkgname"]
ver = parsed_args["ver"]
@assert !isnothing(ver)
@assert !isnothing(pkg)
up_deps_only = isnothing(pkg) && isnothing(ver)
root = dirname(@__DIR__)
dirs =
    (root, joinpath(root, "test"), joinpath(root, "perf"), joinpath(root, "docs"), joinpath(root, "integration_tests"))

cd(root) do
    for dir in dirs
        @info "Pkg.up for environment $dir"
        pkgname = "\"$pkg\""
        vernum = "\"$ver\""
        cmd = `$(Base.julia_cmd()) --project=$dir -e """import Pkg; Pkg.rm($(pkgname)); Pkg.add(Pkg.PackageSpec(name=$(pkgname),version=$(vernum)))"""`
        run(cmd)
        @warn "Package $pkgname was removed from the Project toml from directory $dir. Please add it back in!"
        # run(`git checkout $dir/Project.toml`) # for convenience.
    end
end

# https://github.com/JuliaLang/Pkg.jl/issues/3014
for dir in dirs
    cd(dir) do
        rm("LocalPreferences.toml"; force = true)
    end
end
