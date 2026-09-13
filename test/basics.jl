using NaturalEarth
using Test
using VectorDataCubes

using Rasters, DimensionalData
using Rasters.Lookups
import GeometryOps as GO, GeoInterface as GI
using Extents

import DimensionalData as DD

@testset "construction" begin
    # fetch land polygons from Natural Earth
    country_polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
    # create a DimVector from the land polygons
    gl = GeometryLookup(country_polygons)
    @test crs(gl) == EPSG(4326)
    @test all(splat(GO.equals), zip(val(gl), country_polygons))
end


@testset "Interacting with the geometry lookup" begin
    # fetch land polygons from Natural Earth
    country_fc = NaturalEarth.naturalearth("admin_0_countries", 110)
    country_polygons = country_fc.geometry
    # create a DimVector from the land polygons
    gl = GeometryLookup(country_polygons)
    ras = rand(Geometry(gl))
    @testset "indexing" begin
        # TODO: fix this to only return a single index, if necessary...
        @test only(ras[Geometry = (X(At(-103)), Y(At(44)))]) == ras[findfirst(==("USA"), country_fc.ADM0_A3)]
    end

    # Every `Where(GO.pred(g))` form is covered in test/selectors.jl.
    @test ras[Geometry=Where(GO.disjoint(gl[1]))] == ras[Geometry=2:DD.End()]
end