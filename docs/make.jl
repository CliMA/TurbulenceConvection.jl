using TurbulenceConvection, Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

#! format: off
pages = Any[
    "Home" => "index.md",
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

deploydocs(
    repo = "github.com/CliMA/TurbulenceConvection.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
