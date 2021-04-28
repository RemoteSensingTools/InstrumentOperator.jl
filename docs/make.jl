push!(LOAD_PATH,joinpath(@__DIR__, "../src"))

using Documenter, InstrumentOperator, Literate, Plots

function build()

    tutorials = ["CrossSection_convolution.jl"] # , 
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials"))
    end
    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
            "Home"                  => "index.md",
            "Tutorials"             => tutorials_md
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
end

build()

deploydocs(
    repo = "github.com/RadiativeTransfer/InstrumentOperator.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true
)
