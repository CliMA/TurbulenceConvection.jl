#=
A simple script for updating the manifest
files in all of our environments.
=#

root = dirname(@__DIR__)
dirs =
    (root, joinpath(root, "test"), joinpath(root, "perf"), joinpath(root, "docs"), joinpath(root, "integration_tests"))

for dir in dirs
    cmd = `$(Base.julia_cmd()) --project=$dir 'import Pkg; Pkg.up()'`
    run(cmd)
end

# https://github.com/JuliaLang/Pkg.jl/issues/3014
for dir in dirs
    cd(dir) do
        rm("LocalPreferences.toml"; force = true)
    end
end
