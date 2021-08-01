# Name changes
# z => zf
import Glob

replacements = Dict(
    "delta" => "δ",
    "kappa" => "κ",
    "beta" => "β",
    "gamma" => "γ",
    "qtg" => "q_tot_g",
    "Gr" => "grid",
    "zeta" => "ζ",
    "psi" => "ψ",
    "temperature_half" => "Tc",
    "temperature" => "T",
    "qt" => "q_tot",
    "ql" => "q_liq",
    "qi" => "q_ice",
    "qr" => "q_rai",
    "qv" => "q_vap",
    "thetal" => "θ_liq_ice",
    "theta" => "θ",
    "THL" => "θ_liq_ice",
    "alpha0_half" => "α0_c",
    # "alpha" => "α0_f", # may need to be done manually
    "rho0_half" => "ρ0_c",
    # "rho" => "ρ0_f", # may need to be done manually
    "p0_half" => "p0_c",
    # "p" => "p0_f", # may need to be done manually
    )
# for filename in readdir(".")
all_files = [
    joinpath(root, f)
    for (root, dirs, files) in Base.Filesystem.walkdir(".")
    for f in files
]
all_files = filter(x->!occursin(".git", x), all_files)
all_files = filter(x->!occursin("Output.", x), all_files)
all_files = filter(x->!occursin(".dev", x), all_files)
all_files = filter(x->occursin(".jl", x), all_files)
all_files = filter(x->!occursin("rename_variables.jl", x), all_files)
all_files = filter(x->!occursin("viz", x), all_files)
all_files = filter(x->!occursin("/test/", x), all_files)
all_files = filter(x->!occursin("/docs/", x), all_files)
for filename in all_files
    @show filename
    contents = open(filename, "r") do io
        read(io, String)
    end
    for p in replacements
        contents = replace(contents, p)
    end
    open(filename, "w") do io
        write(io, contents)
    end

end

