using Test
using VectorDataCubes

using Rasters, DimensionalData
using Rasters.Lookups
import DE9IM
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO
import Tables
import Proj
import DataAPI

const USP = GO.UnitSpherical.UnitSphericalPoint
usp(point) = USP(point)

struct USPTable{G,L,C}
    geometry::G
    label::L
    crs::C
end
Tables.istable(::Type{<:USPTable}) = true
Tables.columnaccess(::Type{<:USPTable}) = true
Tables.columns(table::USPTable) = table
Tables.columnnames(::USPTable) = (:geometry, :label)
Tables.getcolumn(table::USPTable, name::Symbol) = getfield(table, name)
Tables.getcolumn(table::USPTable, i::Int) = Tables.getcolumn(table, Tables.columnnames(table)[i])
GI.geometrycolumns(::USPTable) = (:geometry,)
GI.crs(table::USPTable) = table.crs
DataAPI.colmetadatasupport(::Type{<:USPTable}) = (read=true, write=false)
DataAPI.colmetadata(::USPTable, col, key, default) =
    col == :geometry && key == "edges" ? "spherical" : default

@testset "UnitSpherical input is stored as longitude/latitude" begin
    sphere = GO.Spherical()
    points = usp.([(179.0, 0.0), (-179.0, 5.0)])
    lookup = GeometryLookup(points; manifold=sphere)

    @test all(p -> !(p isa USP), parent(lookup))
    @test all(p -> length(p) == 2, parent(lookup))
    @test all(isapprox.(parent(lookup)[1], (179.0, 0.0)))
    @test GI.crs(lookup) === nothing
    bare = GeometryLookup(first(points); manifold=sphere)
    @test GI.crs(bare) === nothing
    @test parent(bare) == parent(lookup)[1:1]

    explicit = GeometryLookup(points; manifold=sphere, crs=EPSG(4326))
    @test GI.crs(explicit) == EPSG(4326)
    @test all(p -> !(p isa USP), parent(explicit))

    rebuilt = DD.rebuild(lookup; data=reverse(points))
    @test all(p -> !(p isa USP), parent(rebuilt))
    @test all(isapprox.(parent(rebuilt)[1], (-179.0, 5.0)))

    xyz = GI.Point(1.0, 2.0, 3.0)
    xyz_lookup = GeometryLookup([xyz]; manifold=sphere)
    @test parent(xyz_lookup)[1] === xyz
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
    at = At(usp((179.0, 0.0)))
    near = Near(usp((-178.0, 0.0)))
    contains = Contains(usp((179.0, 0.0)))
    @test Lookups.selectindices(point_lookup, at) == 1
    @test Lookups.selectindices(point_lookup, near) == 2
    @test Lookups.selectindices(point_lookup, contains) == [1]
    @test Lookups.hasselection(point_lookup, at)
    @test Lookups.hasselection(point_lookup, near)
    @test Lookups.hasselection(point_lookup, contains)
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
    sampled = VectorDataCubes.extract(raster, points; manifold=sphere, crs=EPSG(4326))
    @test collect(sampled) == [1.0, 4.0]
    @test DD.lookup(sampled, Geometry).manifold == sphere
    @test all(p -> !(p isa USP), parent(DD.lookup(sampled, Geometry)))

    emitted = vectordatacubetable(cube)
    @test all(p -> !(p isa USP), emitted.Geometry)
    @test GI.crs(emitted) == EPSG(4326)

    cartesian = USPTable(points, ["a", "b"], GI.crs(first(points)))
    err = try
        vectordatacube(cartesian; manifold=sphere)
    catch exception
        exception
    end
    @test err isa ArgumentError
    @test occursin("geographic CRS", err.msg)
    corrected = vectordatacube(cartesian; manifold=sphere, crs=EPSG(4326))
    @test GI.crs(DD.lookup(corrected, Geometry)) == EPSG(4326)
end

@testset "USP normalization precedes datum inference" begin
    points = usp.([(1.0, 2.0), (3.0, 4.0)])
    customcrs = ProjString("+proj=longlat +R=1234567.89 +type=crs")
    table = USPTable(points, [1, 2], customcrs)
    for lookup in (GeometryLookup(table), DD.lookup(vectordatacube(table), Geometry))
        @test lookup.manifold.radius == 1234567.89
        @test GI.crs(lookup) == customcrs
        @test all(p -> !(p isa USP), parent(lookup))
        emitted = vectordatacubetable(Raster([1, 2], Geometry(lookup); name=:value))
        @test GI.crs(emitted) == customcrs
        @test DD.lookup(vectordatacube(emitted), Geometry).manifold == lookup.manifold
    end

    plain = (geometry=points, value=[1, 2])
    lookup = DD.lookup(vectordatacube(plain; manifold=GO.Spherical()), Geometry)
    @test GI.crs(lookup) === nothing
    @test all(p -> !(p isa USP), parent(lookup))

    inferred = DD.lookup(vectordatacube(USPTable(points, [1, 2], nothing)), Geometry)
    @test GI.crs(inferred) === nothing
    @test inferred.manifold == GO.Spherical()
end
