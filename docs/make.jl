push!(LOAD_PATH,"./src/")

using Documenter, InstrumentOperator

makedocs(sitename="InstrumentOperator Documentation")

deploydocs(
    repo = "github.com/RadiativeTransfer/InstrumentOperator.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true,
)