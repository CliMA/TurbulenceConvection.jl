using TurbulenceConvection, Documenter
using DocumenterCitations
import Glob

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

#! format: off
pages = Any[
    "Home" => "index.md",
    "Formulation" => "Formulation.md",
    "Reference states" => "ReferenceStates.md",
    "EDMF equations" => "EDMFEquations.md",
    "API" => "API.md",
    "Developer docs" => "dev.md",
    "References" => "References.md",
]

mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)
#! format: on

makedocs(
    bib,
    sitename = "TurbulenceConvection.jl",
    strict = true,
    format = format,
    checkdocs = :exports,
    clean = true,
    doctest = true,
    modules = [TurbulenceConvection],
    pages = pages,
)

# Clean up
build_dir = joinpath(@__DIR__, "build")
if isdir(build_dir)
    cd(build_dir) do
        for filename in Glob.glob("Output.*")
            rm(filename; force = true, recursive = true)
        end
        for filename in Glob.glob("*.in")
            rm(filename; force = true)
        end
    end
end

deploydocs(
    repo = "github.com/CliMA/TurbulenceConvection.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
