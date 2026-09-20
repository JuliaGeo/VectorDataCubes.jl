using Test
using VectorDataCubes

using Rasters, DimensionalData
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO

ring(points) = GI.Polygon([GI.LinearRing([points..., first(points)])])
dateline = ring([(170.0, -10.0), (-170.0, -10.0), (-170.0, 10.0), (170.0, 10.0)])
asraster(lookup) = Raster([1], (Geometry(lookup),); name=:value)

@testset "without Proj" begin
    testproject = dirname(Base.active_project())
    script = joinpath(@__DIR__, "spherical_datum_no_proj.jl")
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$testproject $script`
    @test success(cmd)
end

import Proj

@testset "datum-derived spherical radius" begin
    @test !isnothing(Base.get_extension(VectorDataCubes, :VectorDataCubesProjExt))

    wgs = GeometryLookup([dateline]; manifold=GO.Spherical(), crs=EPSG(4326))
    table = vectordatacubetable(asraster(wgs))
    back = DD.lookup(vectordatacube(table), Geometry)
    expected = (2 * 6378137.0 + 6356752.314245179) / 3
    @test back.manifold.radius ≈ expected atol=1e-6
    @test GI.crs(back) == EPSG(4326)
    @test parent(back) == parent(wgs)

    spherecrs = ProjString("+proj=longlat +R=1234567.89 +type=crs")
    custom = GeometryLookup(
        [dateline]; manifold=GO.Spherical(radius=1234567.89), crs=spherecrs
    )
    customtable = vectordatacubetable(asraster(custom))
    customback = DD.lookup(vectordatacube(customtable), Geometry)
    @test customback.manifold.radius == 1234567.89
    @test GI.crs(customback) == spherecrs

    # A table may leave its CRS unset while its geometry values carry it.
    point = GI.Wrappers.Point(1.0, 2.0; crs=spherecrs)
    pointcube = asraster(GeometryLookup([point]; manifold=GO.Spherical()))
    geometrycrstable = VectorDataCubes.VectorDataCubeTable(
        DD.DimTable(pointcube), (:Geometry,), nothing,
        (Geometry=(; edges="spherical"),),
    )
    fromgeometry = DD.lookup(vectordatacube(geometrycrstable), Geometry)
    @test fromgeometry.manifold.radius == 1234567.89
    @test GI.crs(fromgeometry) == spherecrs

    mismatch = GeometryLookup(
        [dateline]; manifold=GO.Spherical(radius=1234568.0), crs=spherecrs
    )
    @test_throws "disagrees" vectordatacubetable(asraster(mismatch))

    projected = GeometryLookup(
        [dateline]; manifold=GO.Spherical(), crs=EPSG(3857)
    )
    @test_throws "projected" vectordatacubetable(asraster(projected))
    @test_throws "projected" Rasters.setcrs(wgs, EPSG(3857))

    relabeled = Rasters.setcrs(custom, spherecrs)
    @test relabeled.manifold === custom.manifold
    @test parent(relabeled) === parent(custom)
end
