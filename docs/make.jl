using EntropyHub
using Documenter, DocumenterTools 
#OutdatedWarning.generate("src")

DocMeta.setdocmeta!(EntropyHub, :DocTestSetup, :(using EntropyHub); recursive=true)

makedocs(
    source="src",
    modules=[EntropyHub],
    authors="Matthew W. Flood <info@entropyhub.xyz>",
    #repo="https://github.com/MattWillFlood/EntropyHub.jl/blob/{commit}{path}#{line}",
    #repo = Remotes.repourl("https://github.com/MattWillFlood/EntropyHub.jl"),
    sitename="EntropyHub.jl",
    doctest=false,
    draft=false,
    clean=true,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://mattwillflood.github.io/EntropyHub.jl/",
        assets = ["assets/favicon.ico"],
        collapselevel = 1,
        footer=" [www.EntropyHub.xyz](https://www.EntropyHub.xyz) ",
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => ["Base Entropies" => "Guide/Base_Entropies.md",
                "Cross-Entropies" => "Guide/Cross_Entropies.md",
                "Multivariate Entropies" => "Guide/Multivariate_Entropies.md",
                "Bidimensional Entropies" => "Guide/Bidimensional_Entropies.md",
                "Multiscale Entropies" => "Guide/Multiscale_Entropies.md",
                "Multiscale Cross-Entropies" => "Guide/Multiscale_Cross_Entropies.md",
                "Multivariate Multiscale Entropies" => "Guide/Multivariate_Multiscale_Entropies.md",
                "Other Functions" => "Guide/Other.md",
                ],
        "Examples" => ["Notes on Examples" => "Examples/Examples.md",
                "Ex.1: Sample Entropy" =>  "Examples/Example1.md",
                "Ex.2: Permutation Entropy" =>  "Examples/Example2.md",
                "Ex.3: Phase Entropy" =>  "Examples/Example3.md",
                "Ex.4: Cross-Distribution Entropy" =>  "Examples/Example4.md",
                "Ex.5: Multiscale Entropy Object" =>  "Examples/Example5.md",
                "Ex.6: Multiscale [Increment] Entropy" =>  "Examples/Example6.md",
                "Ex.7: Refined Multiscale [Sample] Entropy" =>  "Examples/Example7.md",
                "Ex.8: Composite Multiscale Cross-Approximate Entropy" =>  "Examples/Example8.md",
                "Ex.9: Hierarchical Multiscale corrected Cross-Conditional Entropy" =>  "Examples/Example9.md",
                "Ex.10: Bidimensional Fuzzy Entropy" =>  "Examples/Example10.md",
                "Ex.11: Multivariate Dispersion Entropy" =>  "Examples/Example11.md",
                "Ex.12: [Generalized] Refined-composite Multivariate Multiscale Fuzzy Entropy" =>  "Examples/Example12.md",
                "Ex.13: Window Data Tool" =>  "Examples/Example13.md",
                ],
    ],
)

deploydocs(
    repo="github.com/MattWillFlood/EntropyHub.jl.git",
    versions = ["stable" => "v^", "v#.#"],
    branch = "gh-pages",

    #versions = nothing,
)
#=
deploydocs(;
    repo="github.com/MattWillFlood/EntropyHub.jl",
    devbranch = "master",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#.#", devurl =>devurl],
)
"""
=#