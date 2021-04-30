push!(LOAD_PATH,joinpath(@__DIR__, "../src"))

using Documenter, InstrumentOperator, Literate, Plots
using DocumenterTools: Themes


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
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
            ],
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = mathengine,
        collapselevel = 1,
        )
    # download the themes
    for file in ("radiativetransfer-lightdefs.scss", "radiativetransfer-darkdefs.scss", "radiativetransfer-style.scss")
        download("https://raw.githubusercontent.com/RadiativeTransfer/doctheme/main/$file", joinpath(@__DIR__, file))
    end
    # create the themes
    for w in ("light", "dark")
        header = read(joinpath(@__DIR__, "radiativetransfer-style.scss"), String)
        theme = read(joinpath(@__DIR__, "radiativetransfer-$(w)defs.scss"), String)
        write(joinpath(@__DIR__, "radiativetransfer-$(w).scss"), header*"\n"*theme)
    end
    # compile the themes
    Themes.compile(joinpath(@__DIR__, "radiativetransfer-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
    Themes.compile(joinpath(@__DIR__, "radiativetransfer-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

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
