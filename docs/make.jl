using Documenter
using AbsSmoothFrankWolfe


makedocs(;
    sitename = "AbsSmoothFrankWolfe.jl",
    modules = [AbsSmoothFrankWolfe],
    format = Documenter.HTML(),
    pages = [
            "Home" => "index.md",
            "Examples" => "examples.md",
            "References" => "references.md"
   ],
)

deploydocs(; repo = "github.com/ZIB-IOL/AbsSmoothFrankWolfe.jl.git", devbranch="main", push_preview=true)

