using Test
using VectorDataCubes

using Rasters, DimensionalData
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO

@test isnothing(Base.get_extension(VectorDataCubes, :VectorDataCubesProjExt))

ring(points) = GI.Polygon([GI.LinearRing([points..., first(points)])])
dateline = ring([(170.0, -10.0), (-170.0, -10.0), (-170.0, 10.0), (170.0, 10.0)])
asraster(lookup) = Raster([1], (Geometry(lookup),); name=:value)

default = GeometryLookup([dateline]; manifold=GO.Spherical())
table = vectordatacubetable(asraster(default))
@test DD.lookup(vectordatacube(table), Geometry).manifold == GO.Spherical()

custom = GeometryLookup([dateline]; manifold=GO.Spherical(radius=1.0))
@test_throws "without a CRS" vectordatacubetable(asraster(custom))

withcrs = GeometryLookup([dateline]; manifold=GO.Spherical(), crs=EPSG(4326))
@test_throws "requires Proj.jl" vectordatacubetable(asraster(withcrs))
@test_throws "requires Proj.jl" Rasters.setcrs(default, EPSG(4326))
