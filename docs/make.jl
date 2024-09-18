using Documenter
using AbsSmoothFrankWolfe

generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/ZIB-IOL/AbsSmoothFrankWolfe.jl/tree/main/"
isdir(generated_path) || mkdir(generated_path)

open(joinpath(generated_path, "index.md"), "w") do io
    # Point to source license file
    println(
        io,
        """
        ```@meta
        EditURL = "$(base_url)README.md"
        ```
        """,
    )
    # Write the contents out below the meta block
    for line in eachline(joinpath(dirname(@__DIR__), "README.md"))
        println(io, line)
    end
end

makedocs(;
    sitename = "AbsSmoothFrankWolfe.jl",
    modules = [AbsSmoothFrankWolfe],
    format = Documenter.HTML(),
    pages = [
            "Home" => "index.md"
   ],
)

deploydocs(; repo = "github.com/ZIB-IOL/AbsSmoothFrankWolfe.jl/tree/main", push_preview=true)

