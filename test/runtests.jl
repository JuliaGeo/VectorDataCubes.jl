using SafeTestsets

# Each file runs in its own module (`@safetestset`), so every test file must be
# self-contained: it imports VectorDataCubes and everything else it needs. Nothing
# leaks in from here.
@safetestset "construction" begin include("construction.jl") end
@safetestset "selectors" begin include("selectors.jl") end
@safetestset "zonal" begin include("zonal.jl") end
@safetestset "extract" begin include("extract.jl") end
@safetestset "tables" begin include("tables.jl") end
@safetestset "makie" begin include("makie.jl") end
@safetestset "basics" begin include("basics.jl") end
