push!(LOAD_PATH,"./src/")

using Documenter, InstrumentOperator

pages = Any[
        "Home"                  => "index.md"
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

makedocs(
    sitename = "Instrument Operator",
    format = format,
    clean = false,
    modules = [InstrumentOperator],
    pages = pages
)

deploydocs(
    repo = "github.com/RadiativeTransfer/InstrumentOperator.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true
)
