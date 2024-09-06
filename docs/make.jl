using Documenter
using AbsSmoothFW

push!(LOAD_PATH,"../src/")
makedocs(;
    sitename = "AbsSmoothFW.jl",
    modules = [AbsSmoothFW],
    format = Documenter.HTML(),
    pages = [
            "Home" => "index.md",
   ],
)

deploydocs(; repo = "github.com/shtadinada/AbsSmoothFW.jl.git", devbranch = "main")


