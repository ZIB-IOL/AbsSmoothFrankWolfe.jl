using Documenter
using AbsSmoothFW

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
    sitename = "AbsSmoothFW.jl",
    modules = [AbsSmoothFW],
    format = Documenter.HTML(),
    pages = [
            "Home" => "index.md",
   ],
)

deploydocs(; repo = "github.com/ZIB-IOL/AbsSmoothFW.jl", push_preview=true)

