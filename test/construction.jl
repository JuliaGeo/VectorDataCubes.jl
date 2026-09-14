using Test
using VectorDataCubes
using VectorDataCubes: spatialtree

using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
using GeometryOps.FlexibleRTrees: RTree, STR, HPR, Unsorted
import Extents

_csquare(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

csq1 = _csquare(0.0, 0.0, 1.0, 1.0)
csq2 = _csquare(1.0, 0.0, 2.0, 1.0)
csq3 = _csquare(10.0, 10.0, 11.0, 11.0)
csquares = [csq1, csq2, csq3]

# The tree holder is internal; tests only need to know whether a tree has been built yet.
_built(gl) = !isnothing(VectorDataCubes._builttree(gl))

@testset "constructor variants" begin
    @testset "from a vector of geometries" begin
        gl = GeometryLookup(csquares)
        @test all(splat(GO.equals), zip(val(gl), csquares))
        @test spatialtree(gl) isa RTree{STR}
        @test gl.manifold == GO.Planar()
    end

    @testset "from a table with geometrycolumn" begin
        tbl = (; geom = csquares, value = 1:3)
        gl = GeometryLookup(tbl; geometrycolumn = :geom)
        @test all(splat(GO.equals), zip(val(gl), csquares))
    end

    @testset "Union{Missing} eltype is narrowed" begin
        geoms = Union{Missing, eltype(csquares)}[csquares...]
        gl = GeometryLookup(geoms)
        @test !(Missing <: eltype(val(gl)))
        @test length(gl) == 3
    end

    @testset "empty input" begin
        gl = GeometryLookup(empty(csquares))
        @test isempty(gl)
        @test spatialtree(gl) === nothing
        @test Lookups.bounds(gl) == ((nothing, nothing), (nothing, nothing))
        dv = rand(Geometry(gl))
        @test isempty(dv)
        @test Extents.extent(dv) == Extents.Extent(X = (nothing, nothing), Y = (nothing, nothing))
        @test isempty(dv[X(0 .. 1), Y(0 .. 1)])
        @test isempty(Rasters.crop(Raster(dv); to = Extents.Extent(X = (0.0, 1.0), Y = (0.0, 1.0))))
    end

    @testset "internal dims" begin
        @test DD.name.(DD.dims(GeometryLookup(csquares))) == (:X, :Y)
        @test DD.name.(DD.dims(GeometryLookup(csquares, (Y(), X())))) == (:Y, :X)
        @test DD.name.(DD.dims(GeometryLookup(csquares, (X, Y)))) == (:X, :Y)
        @test_throws ArgumentError GeometryLookup(csquares, (Dim{:a}(), Dim{:b}()))
        @test_throws ArgumentError GeometryLookup(csquares, (X(), X()))
        @test_throws ArgumentError GeometryLookup(csquares, (X(),))
        @test_throws ArgumentError GeometryLookup(csquares, (X(), Y(), Ti()))
    end

    @testset "error paths" begin
        @test_throws ArgumentError GeometryLookup([1, 2, 3])
        @test_throws ArgumentError GeometryLookup(csquares; tree = 42)
        withmissing = Union{Missing, eltype(csquares)}[csq1, missing, csq3]
        err = try
            GeometryLookup(withmissing)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("[2]", err.msg)
    end
end

@testset "tree keyword" begin
    @testset "default tree is built lazily, once" begin
        gl = GeometryLookup(csquares)
        @test !_built(gl)
        tree = spatialtree(gl)
        @test tree isa RTree{STR}
        @test tree.data === parent(gl)
        @test spatialtree(gl) === tree
    end

    @testset "the tree type is known from the lookup type" begin
        gl = GeometryLookup(csquares)
        @test only(Base.return_types(spatialtree, (typeof(gl),))) == Union{Nothing, typeof(spatialtree(gl))}
        @test isconcretetype(typeof(spatialtree(gl)))
        # a view keeps its own, still concrete, tree type
        v = view(gl, 1:2)
        @test only(Base.return_types(spatialtree, (typeof(v),))) == Union{Nothing, typeof(spatialtree(v))}
        # so does a geometry vector whose eltype tells the compiler nothing
        anygl = DD.rebuild(gl; data = Any[csquares...])
        @test only(Base.return_types(spatialtree, (typeof(anygl),))) == Union{Nothing, typeof(spatialtree(anygl))}
        @test isconcretetype(typeof(spatialtree(anygl)))
        @test spatialtree(GeometryLookup(csquares; tree = nothing)) === nothing
    end

    @testset "slicing, view, reverse and rebuild do not build" begin
        gl = GeometryLookup(csquares)
        spatialtree(gl)
        dv = rand(Geometry(gl))
        @test !_built(DD.lookup(dv[1:2], Geometry))
        @test !_built(DD.lookup(view(dv, 1:2), Geometry))
        @test !_built(reverse(gl))
        @test !_built(DD.rebuild(gl; data = csquares[1:2]))
        # identical data keeps the built tree
        @test DD.rebuild(gl; data = parent(gl)).tree === gl.tree
        @test spatialtree(DD.lookup(dv[1:2], Geometry)) isa RTree{STR}
    end

    @testset "the lookup type does not depend on the geometry count" begin
        gl = GeometryLookup(csquares)
        @test typeof(gl[1:1]) == typeof(gl[1:3])
        @test typeof(DD.rebuild(gl; data = empty(csquares))) == typeof(gl)
    end

    @testset "tree = nothing disables the accelerator and survives" begin
        gl = GeometryLookup(csquares; tree = nothing)
        @test spatialtree(gl) === nothing
        @test spatialtree(gl[1:2]) === nothing
        @test spatialtree(view(gl, 1:2)) === nothing
        @test spatialtree(reverse(gl)) === nothing
        @test spatialtree(DD.rebuild(gl; data = csquares[1:2])) === nothing
        dv = rand(Geometry(gl))
        @test spatialtree(DD.lookup(dv[Geometry = Contains((0.5, 0.5))], Geometry)) === nothing
    end

    @testset "concurrent queries build one tree" begin
        gl = GeometryLookup(csquares)
        trees = Vector{Any}(undef, 8)
        Threads.@sync for i in eachindex(trees)
            Threads.@spawn trees[i] = spatialtree(gl)
        end
        @test all(t -> t === first(trees), trees)
        @test first(trees) === spatialtree(gl)
    end

    @testset "tree as a bulk-load algorithm" begin
        for algorithm in (STR(), HPR(), Unsorted())
            gl = GeometryLookup(csquares; tree = algorithm)
            @test !_built(gl)
            @test spatialtree(gl) isa RTree{typeof(algorithm)}
            @test spatialtree(gl[1:2]) isa RTree{typeof(algorithm)}
            # every algorithm gives the one tree type the lookup type promises
            @test typeof(spatialtree(gl)) == VectorDataCubes.XYRTree{typeof(algorithm), typeof(csquares)}
            dv = DimArray([1, 2, 3], Geometry(gl))
            @test dv[Geometry = Contains((10.5, 10.5))] == [3]
        end
    end

    @testset "tree as a prebuilt instance" begin
        tree = RTree(HPR(), csquares)
        gl = GeometryLookup(csquares; tree)
        @test spatialtree(gl) === tree
        @test spatialtree(gl[1:2]) isa RTree{HPR}
        # a tree over another vector cannot index this lookup
        @test_throws ArgumentError GeometryLookup(csquares; tree = RTree(STR(), copy(csquares)))
    end

    @testset "rebuild with an explicit tree" begin
        gl = GeometryLookup(csquares)
        @test spatialtree(DD.rebuild(gl; tree = nothing)) === nothing
        rb = DD.rebuild(gl; data = csquares[1:2], tree = Unsorted())
        @test spatialtree(rb) isa RTree{Unsorted}
    end
end

@testset "crs" begin
    @testset "no crs anywhere gives nothing" begin
        gl = GeometryLookup(csquares)
        @test GI.crs(gl) === nothing
        @test crs(gl) === nothing
    end

    @testset "explicit crs keyword" begin
        gl = GeometryLookup(csquares; crs = EPSG(4326))
        @test crs(gl) == EPSG(4326)
        @test crs(Geometry(gl)) == EPSG(4326)
    end

    @testset "setcrs on the lookup" begin
        gl = GeometryLookup(csquares; crs = EPSG(4326))
        gl2 = Rasters.setcrs(gl, EPSG(3857))
        @test crs(gl2) == EPSG(3857)
        # only the crs changed - data and tree holder are reused
        @test parent(gl2) === parent(gl)
        @test gl2.tree === gl.tree
    end

    @testset "setcrs on dimensions and cubes" begin
        gl = GeometryLookup(csquares; crs = EPSG(4326))
        dv = rand(Geometry(gl))
        @test crs(Rasters.setcrs(dims(dv, Geometry), EPSG(3857))) == EPSG(3857)
        @test crs(dims(dv, Geometry)) == EPSG(4326)
        r = Rasters.setcrs(Raster(dv), EPSG(3857))
        @test crs(dims(r, Geometry)) == EPSG(3857)
        @test parent(r) == parent(dv)
        cube = Rasters.setcrs(Raster(rand(Geometry(gl), Ti(1:2))), EPSG(3857))
        @test crs(dims(cube, Geometry)) == EPSG(3857)
    end

    @testset "reproject" begin
        gl = GeometryLookup(csquares; crs = EPSG(4326))
        # Another test file may already have loaded Proj into this session.
        if isnothing(Base.get_extension(GO, :GeometryOpsProjExt))
            @test_throws MethodError reproject(EPSG(3857), gl)
        else
            @test crs(reproject(EPSG(3857), gl)) == EPSG(3857)
        end
        @test_throws ArgumentError reproject(EPSG(3857), GeometryLookup(csquares))
    end
end

@testset "equality and hashing" begin
    gl = GeometryLookup(csquares)
    same = GeometryLookup(csquares)
    @test gl == gl
    @test gl == same
    @test isequal(gl, same)
    @test hash(gl) == hash(same)
    # the tree does not take part
    spatialtree(gl)
    @test gl == same
    @test hash(gl) == hash(same)
    @test gl == GeometryLookup(csquares; tree = nothing)

    withcrs = GeometryLookup(csquares; crs = EPSG(4326))
    @test gl != withcrs
    @test !isequal(gl, withcrs)
    @test hash(gl) != hash(withcrs)
    @test withcrs == GeometryLookup(csquares; crs = EPSG(4326))
    @test hash(withcrs) == hash(GeometryLookup(csquares; crs = EPSG(4326)))
    @test withcrs != GeometryLookup(csquares; crs = EPSG(3857))

    @test gl != GeometryLookup(csquares[1:2])
    @test gl != GeometryLookup(csquares, (Y(), X()))
    @test gl == GeometryLookup(csquares[1:3])
    @test hash(gl) == hash(GeometryLookup(csquares[1:3]))
end

@testset "DimensionalData interface" begin
    gl = GeometryLookup(csquares)

    @testset "dims, order, parent" begin
        @test DD.name.(DD.dims(gl)) == (:X, :Y)
        @test DD.name.(DD.dims(Geometry(gl))) == (:X, :Y)
        @test Lookups.hasinternaldimensions(gl)
        @test DD.order(gl) == Lookups.Unordered()
        @test parent(gl) === val(gl)
        @test length(gl) == 3
    end

    @testset "Geometry dimension" begin
        @test Geometry <: DD.Dimension
        @test DD.name(Geometry) == :Geometry
        dv = rand(Geometry(gl))
        @test DD.dims(dv, Geometry) isa Geometry
        @test val(DD.lookup(dv, Geometry)) == val(gl)
    end

    @testset "rebuild" begin
        rb_new = DD.rebuild(gl; data = csquares[1:2])
        @test length(rb_new) == 2
        @test rb_new.tree !== gl.tree
        rb_empty = DD.rebuild(gl; data = empty(csquares))
        @test isempty(rb_empty)
        @test spatialtree(rb_empty) === nothing
        @test crs(DD.rebuild(gl; crs = EPSG(4326))) == EPSG(4326)
        @test DD.name.(DD.dims(DD.rebuild(gl; dims = (Y(), X())))) == (:Y, :X)
        @test_throws ArgumentError DD.rebuild(gl; dims = (Dim{:a}(), Dim{:b}()))
        @test_throws ArgumentError DD.rebuild(gl; dims = (X(),))
    end

    @testset "reverse" begin
        rgl = reverse(gl)
        @test rgl isa GeometryLookup
        @test parent(rgl) == reverse(csquares)
        dv = DimArray([1, 2, 3], Geometry(gl))
        rv = reverse(dv; dims = Geometry)
        @test DD.lookup(rv, Geometry) isa GeometryLookup
        @test rv[Geometry = Contains((10.5, 10.5))] == [3]
        @test rv[Geometry = Contains((0.5, 0.5))] == [1]
    end

    @testset "set" begin
        dv = rand(Geometry(gl))
        new = GeometryLookup(csquares; crs = EPSG(4326))
        @test DD.lookup(set(dv, Geometry => new), Geometry) === new
        @test DD.lookup(set(dv, Geometry(new)), Geometry) === new

        # setting the values keeps the lookup, and checks them like the constructor
        moved = [_csquare(0.0, 0.0, 0.5, 0.5), _csquare(1.0, 0.0, 1.5, 0.5), _csquare(2.0, 2.0, 3.0, 3.0)]
        withcrs = rand(Geometry(GeometryLookup(csquares; crs = EPSG(4326))))
        l = DD.lookup(set(withcrs, Geometry => moved), Geometry)
        @test l isa GeometryLookup
        @test parent(l) === moved
        @test crs(l) == EPSG(4326)
        @test !_built(l)
        @test_throws ArgumentError set(dv, Geometry => [1, 2, 3])
        @test_throws ArgumentError set(dv, Geometry => Union{Missing, eltype(csquares)}[csq1, missing, csq3])
    end

    @testset "metadata" begin
        @test DD.metadata(gl) == Lookups.NoMetadata()
        md = Dict(:source => "test")
        glm = GeometryLookup(csquares; metadata = md)
        @test DD.metadata(glm) == md
        @test DD.metadata(DD.rebuild(glm; data = csquares[1:2])) == md
        @test DD.metadata(DD.rebuild(glm; metadata = Lookups.NoMetadata())) == Lookups.NoMetadata()
        @test DD.metadata(glm[1:2]) == md
        dv = rand(Geometry(glm))
        @test DD.metadata(dims(dv, Geometry)) == md
    end

    @testset "bounds and extent" begin
        @test Lookups.bounds(gl) == ((0.0, 11.0), (0.0, 11.0))
        @test Lookups.bounds(Geometry(gl)) == ((0.0, 11.0), (0.0, 11.0))
        @test Lookups.bounds(GeometryLookup(csquares, (Y(), X()))) == ((0.0, 11.0), (0.0, 11.0))
        @test Lookups.bounds(gl[1:2]) == ((0.0, 2.0), (0.0, 1.0))
        @test Lookups.bounds(GeometryLookup([csq3])) == ((10.0, 11.0), (10.0, 11.0))
        @test Lookups.bounds(GeometryLookup(csquares; tree = nothing)) == ((0.0, 11.0), (0.0, 11.0))
        # the same answer before and after the tree exists
        spatialtree(gl)
        @test Lookups.bounds(gl) == ((0.0, 11.0), (0.0, 11.0))
        # including for coordinates the tree widens to Float64
        gl32 = GeometryLookup([_csquare(0f0, 0f0, 1f0, 1f0), _csquare(5f0, 5f0, 6f0, 6f0)])
        before = Lookups.bounds(gl32)
        spatialtree(gl32)
        @test Lookups.bounds(gl32) === before === ((0.0, 6.0), (0.0, 6.0))
        dv = rand(Geometry(gl))
        @test Extents.extent(dv) == Extents.Extent(X = (0.0, 11.0), Y = (0.0, 11.0))
        @test Extents.extent(dv[1:2]) == Extents.Extent(X = (0.0, 2.0), Y = (0.0, 1.0))
    end

    @testset "unwrapped X/Y selectors and crop" begin
        dv = DimArray([1, 2, 3], Geometry(gl))
        @test dv[X(-0.1 .. 1.1), Y(-0.1 .. 1.1)] == [1]
        @test dv[X(-0.1 .. 2.1), Y(-0.1 .. 1.1)] == [1, 2]
        @test dv[X(At(0.5)), Y(At(0.5))] == 1
        @test dv[X(At(10.5)), Y(At(10.5))] == 3
        @test dv[X(Contains(1.5)), Y(Contains(0.5))] == [2]
        cropped = Rasters.crop(Raster(dv); to = Extents.Extent(X = (-0.1, 1.1), Y = (-0.1, 1.1)))
        @test parent(cropped) == [1]
        @test DD.lookup(cropped, Geometry) isa GeometryLookup
    end

    @testset "show" begin
        dv = rand(Geometry(gl))
        str = sprint(show, MIME"text/plain"(), dv)
        @test occursin("Geometry", str)
        @test occursin("GeometryLookup{Polygon}", str)
        @test sprint(Lookups.show_compact, MIME"text/plain"(), gl) == "GeometryLookup{Polygon}"
        @test occursin("GeometryLookup", sprint(show, MIME"text/plain"(), gl))
        # mixed geometries: the element type is named, never spelled out in full
        point = GI.Point((0.0, 0.0))
        mixed = GeometryLookup(Union{typeof(csq1), typeof(point)}[csq1, point])
        @test sprint(Lookups.show_compact, MIME"text/plain"(), mixed) ==
            "GeometryLookup{Union{Point, Polygon}}"
        anyeltype = GeometryLookup(Any[csq1, point])
        @test sprint(Lookups.show_compact, MIME"text/plain"(), anyeltype) == "GeometryLookup{Any}"
    end
end
