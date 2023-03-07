using Glob

if isempty(Glob.glob("TurbulenceConvection.so"))

    using Pkg
    using PackageCompiler
    Pkg.develop(;path="..")
    using TurbulenceConvection

    tc_dir = pkgdir(TurbulenceConvection)
    pkgs = [:TurbulenceConvection]
    append!(pkgs, [Symbol(v.name) for v in values(Pkg.dependencies()) if v.is_direct_dep])

    # Avoid unnecessary pkgs, and packages that are in julia's Base.loaded_modules
    do_not_compile_pkgs = [
        :CairoMakie,
        :Makie,
        :ForwardDiff,
        :PackageCompiler,
        :NPZ,
        :Test,
        :Dates,
        :LinearAlgebra,
        :Statistics,
        :Random,
        :Logging,
        :SparseArrays,
        :TerminalLoggers,
        :OrdinaryDiffEq,
        :StochasticDiffEq,
        :DiffEqBase,
        :SciMLBase,
    ]
    filter!(pkg -> pkg ∉ do_not_compile_pkgs, pkgs)

    create_sysimage(
        pkgs;
        sysimage_path = "TurbulenceConvection.so",
        # Caltech Central CPU architecture, `native` leads to issues as well. # see julia -C help for full list of options..., https://ask.cyberinfrastructure.org/t/how-do-i-get-the-list-of-features-and-resources-of-each-node-in-slurm/201/2 for sinfo help
        cpu_target = "generic",
        # cpu_target = "skylake-avx512", # This one works most of the time, but needs a failsafe mechanism.
        # cpu_target = "icelake; skylake-avx512; broadwell; cascadelake", # sinfo -o "%20N %25f  %10G "
        # cpu_target = "generic;skylake-avx512,clone_all;znver2,clone_all", # from https://github.com/JuliaLang/julia/issues/48579 (still doesn't work...)
        # precompile_execution_file = joinpath(tc_dir, "test", "runtests.jl"),
    )

    # Other cpu_target options
    # cpu_target = "skylake-avx512",
    #cpu_target = "native", # Default 
    #cpu_target = PackageCompiler.default_app_cpu_target(),
    #cpu_target = "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)", # Same as Julia Base
    #cpu_target = "generic",

end