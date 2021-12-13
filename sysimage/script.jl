#!/usr/bin/env julia
#
# Called with no arguments will create the system image
#     TurbulenceConvection.so
# in the `@__DIR__` directory.
#
# Called with a single argument the system image will be placed in the path
# specified by the argument (relative to the callers path)
#
# Called with a specified systemimg path and `true`, the system image will
# compile the TurbulenceConvection package module (useful for CI)

sysimage_path =
    isempty(ARGS) ? joinpath(@__DIR__, "TurbulenceConvection.so") : abspath(ARGS[1])

turbulenceconvection_pkg = get(ARGS, 2, "false") == "true"

@info "Creating system image object file at: '$(sysimage_path)'"
@info "Building TurbulenceConvection into system image: $(turbulenceconvection_pkg)"

start_time = time()

using Pkg
Pkg.add("PackageCompiler")

Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate(verbose = true)

pkgs = Symbol[]
if turbulenceconvection_pkg
    push!(pkgs, :TurbulenceConvection)
else
    append!(
        pkgs,
        [Symbol(v.name) for v in values(Pkg.dependencies()) if v.is_direct_dep],
    )
end

# use package compiler
using PackageCompiler
PackageCompiler.create_sysimage(
    pkgs,
    sysimage_path = sysimage_path,
    precompile_execution_file = joinpath(
        @__DIR__,
        "..",
        "integration_tests",
        "Bomex.jl",
    ),
)

tot_secs = Int(floor(time() - start_time))
@info "Created system image object file at: $(sysimage_path)"
@info "System image build time: $tot_secs sec"
