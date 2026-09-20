using Test
using VectorDataCubes
using VectorDataCubes: extract

using Rasters, DimensionalData
using Rasters.Lookups
import DE9IM
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO

const USP = GO.UnitSpherical.UnitSphericalPoint
usp(point) = USP(point)

@testset "UnitSpherical input is stored as longitude/latitude" begin
    sphere = GO.Spherical()
    points = usp.([(179.0, 0.0), (-179.0, 5.0)])
    lookup = GeometryLookup(points; manifold=sphere)

    @test all(p -> !(p isa USP), parent(lookup))
    @test all(p -> length(p) == 2, parent(lookup))
    @test all(isapprox.(parent(lookup)[1], (179.0, 0.0)))
    @test GI.crs(lookup) === nothing

    explicit = GeometryLookup(points; manifold=sphere, crs=EPSG(4326))
    @test GI.crs(explicit) == EPSG(4326)
    @test all(p -> !(p isa USP), parent(explicit))

    rebuilt = DD.rebuild(lookup; data=reverse(points))
    @test all(p -> !(p isa USP), parent(rebuilt))
    @test all(isapprox.(parent(rebuilt)[1], (-179.0, 5.0)))
end

@testset "whole geometries and spherical queries normalize USP coordinates" begin
    sphere = GO.Spherical()
    lonlat = [(170.0, -10.0), (-170.0, -10.0), (-170.0, 10.0),
              (170.0, 10.0), (170.0, -10.0)]
    polygon = GI.Polygon([GI.LinearRing(usp.(lonlat))])
    lookup = GeometryLookup([polygon]; manifold=sphere)

    @test !VectorDataCubes._hasusp(only(parent(lookup)))
    @test Lookups.selectindices(lookup, Contains(usp((179.0, 0.0)))) == [1]

    point_lookup = GeometryLookup(usp.([(179.0, 0.0), (-179.0, 0.0)]); manifold=sphere)
    @test Lookups.selectindices(point_lookup, At(usp((179.0, 0.0)))) == 1
    @test Lookups.selectindices(point_lookup, Near(usp((-178.0, 0.0)))) == 2
    @test Lookups.selectindices(point_lookup,
        At(usp.([(-179.0, 0.0), (179.0, 0.0)]))) == [2, 1]
    @test Lookups.selectindices(lookup, DE9IM.Covers(usp((179.0, 0.0)))) == [1]
end

@testset "tables and extract emit longitude/latitude" begin
    sphere = GO.Spherical()
    points = usp.([(0.5, 0.5), (1.5, 1.5)])
    table = (geometry=points, label=["a", "b"])
    cube = vectordatacube(table; manifold=sphere, crs=EPSG(4326))
    @test all(p -> !(p isa USP), parent(DD.lookup(cube, Geometry)))

    raster = Raster(reshape(1.0:4.0, 2, 2),
        (X(Sampled(0.5:1.0:1.5; sampling=Intervals(Center()))),
         Y(Sampled(0.5:1.0:1.5; sampling=Intervals(Center())))))
    sampled = extract(raster, points; manifold=sphere, crs=EPSG(4326))
    @test collect(sampled) == [1.0, 4.0]
    @test DD.lookup(sampled, Geometry).manifold == sphere
    @test all(p -> !(p isa USP), parent(DD.lookup(sampled, Geometry)))

    emitted = vectordatacubetable(cube)
    @test all(p -> !(p isa USP), emitted.Geometry)
    @test GI.crs(emitted) == EPSG(4326)
end
