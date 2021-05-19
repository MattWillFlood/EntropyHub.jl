using EntropyHub
using Documenter

DocMeta.setdocmeta!(EntropyHub, :DocTestSetup, :(using EntropyHub); recursive=true)

makedocs(;
    modules=[EntropyHub],
    authors="Matthew W. Flood <entropyhubproject@gmail.com>",
    repo="https://github.com/MattWillFlood/EntropyHub.jl/blob/{commit}{path}#{line}",
    sitename="EntropyHub.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MattWillFlood.github.io/EntropyHub.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MattWillFlood/EntropyHub.jl",
)
