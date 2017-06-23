using Documenter, Dolo

makedocs(
    modules=[Dolo],
    authors="Spencer Lyon, Pablo Winant, and contributors",
    clean=false,
    format=:html,
    sitename="Dolo.jl",
    doctest=true,
    pages=[
        "Guide" => [
            "index.md",
        ],
    ],
)

deploydocs(
    repo   = "github.com/EconForge/Dolo.jl.git",
    julia  = "0.6",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
