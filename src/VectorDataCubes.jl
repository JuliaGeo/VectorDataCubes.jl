module VectorDataCubes

include("geometry_lookup.jl")
include("selectors.jl")
include("zonal.jl")
include("extract.jl")
include("tables.jl")

export GeometryLookup
export Geometry
export vectordatacube, vectordatacubetable
# `zonal` and `extract` are deliberately not exported: Rasters exports both names,
# so use `VectorDataCubes.zonal` / `using VectorDataCubes: zonal`.

end # module VectorDataCubes
