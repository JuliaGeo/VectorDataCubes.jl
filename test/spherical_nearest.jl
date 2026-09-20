using Test, Random, VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import GeometryOps as GO, GeoInterface as GI, Extents
using GeometryOps.FlexibleRTrees: STR, HPR, Unsorted

const unit = GO.Spherical(radius=1.0)
coords(p) = (Float64(GI.x(p)), Float64(GI.y(p)))
reference(points, query) = argmin(p -> GO.distance(unit, coords(query), coords(p)), points)

@testset "indexed spherical nearest agrees with exact scan" begin
    rng = MersenneTwister(27)
    points = [(360rand(rng) - 180, 180rand(rng) - 90) for _ in 1:256]
    append!(points, [(179.0, 0.0), (-179.0, 0.0), (0.0, 90.0), (0.0, -90.0)])
    queries = [(360rand(rng) - 180, 180rand(rng) - 90) for _ in 1:60]
    append!(queries, [(180.0, 0.0), (-180.0, 0.0), (45.0, 90.0), (45.0, -90.0)])
    for data in (points, map(p -> Float32.(p), points))
        scan = GeometryLookup(data; manifold=unit, tree=nothing)
        for algorithm in (STR(), HPR(), Unsorted())
            lookup = GeometryLookup(data; manifold=unit, tree=algorithm)
            @test isnothing(VectorDataCubes._builttree(lookup))
            for query in queries
                index = Lookups.selectindices(lookup, Near(query))
                @test index == Lookups.selectindices(scan, Near(query))
                @test data[index] == reference(data, query)
            end
            @test !isnothing(VectorDataCubes._builttree(lookup))
        end
    end
end

@testset "spherical nearest bounds and ties" begin
    points = [(179.0, 0.0), (0.0, 89.0), (179.0, 0.0)]
    for radius in (1.0, GO.Spherical().radius), algorithm in (STR(), HPR(), Unsorted(), nothing)
        lookup = GeometryLookup(points; manifold=GO.Spherical(; radius), tree=algorithm)
        @test Lookups.selectindices(lookup, Near((180.0, 0.0))) == 1
        @test Lookups.selectindices(lookup, Near((170.0, 89.0))) == 2
        @test Lookups.selectindices(lookup, (X(Near(180.0)), Y(Near(0.0)))) == 1
    end
    for data in ([(0.0, 0.0)], [(179.0, 0.0), (-179.0, 0.0)], [(0.0, 89.0), (90.0, 89.0)])
        tree = VectorDataCubes.spatialtree(GeometryLookup(data; manifold=unit))
        for query in [(180.0, 0.0), (180.0 - 1e-10, 0.0), (0.0, -90.0), (0.0, 90.0)]
            lower = VectorDataCubes._extent_distance(VectorDataCubes._unitpoint(query), Extents.extent(tree))
            @test lower <= minimum(p -> GO.distance(unit, query, p), data)
        end
    end
end
