using Documenter
using AbsSmoothFW

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "AbsSmoothFW.jl Documentation",
    pages = [
            "Index" => "index.md",
    format = Documenter.HTML(prettyurls = false)
 
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/shtadinada/AbsSmoothFW.jl.git",
    devbranch = "main"
)


