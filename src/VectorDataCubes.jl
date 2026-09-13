module VectorDataCubes

include("geometry_lookup.jl")
include("selectors.jl")
include("zonal.jl")
include("tables.jl")

export GeometryLookup
export Geometry
export vectordatacube, vectordatacubetable
# `zonal` is deliberately not exported: Rasters exports a `zonal` too, so use
# `VectorDataCubes.zonal` or `using VectorDataCubes: zonal`.

end # module VectorDataCubes
