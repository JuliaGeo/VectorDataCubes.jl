using VectorDataCubes
using Documenter
using DocumenterVitepress
using Literate

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")
const GENERATED_DIR = joinpath(@__DIR__, "src", "examples")

example_scripts = sort(filter(f -> endswith(f, ".jl"), readdir(EXAMPLES_DIR)))

mkpath(GENERATED_DIR)
example_pages = map(example_scripts) do script
    src = joinpath(EXAMPLES_DIR, script)
    Literate.markdown(src, GENERATED_DIR)
    joinpath("examples", splitext(script)[1] * ".md")
end

makedocs(;
    modules = [VectorDataCubes],
    authors = "Anshul Singhvi <anshulsinghvi@gmail.com> and contributors",
    sitename = "VectorDataCubes.jl",
    repo = "https://github.com/asinghvi17/VectorDataCubes.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/asinghvi17/VectorDataCubes.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => example_pages,
        "API reference" => "api.md",
    ],
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/asinghvi17/VectorDataCubes.jl",
    devbranch = "main",
    push_preview = true,
)
