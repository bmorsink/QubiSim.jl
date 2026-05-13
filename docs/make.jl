using Documenter
using QubiSim

makedocs(
    sitename = "QubiSim",
    format = Documenter.HTML(size_threshold = nothing, edit_link = "main"),
    modules = [QubiSim],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",

        "Tutorials" => [
            "Bell state" => "tutorials_bell_state.md",
            "Examples" => "tutorials_examples.md",
        ],

        "Guide" => [
            "Overview" => "guide_overview.md",
            "Quantum Circuits" => "guide_circuits.md",
            "Compilation" => "guide_compilation.md",
            "Execution" => "guide_execution.md",
            "Results" => "guide_results.md",
        ],

        "Internals" => "internals.md",

        "API Reference" => "api.md",
    ],

)

deploydocs(
    repo = "github.com/bmorsink/QubiSim.jl.git",
    devbranch = "main"
)
