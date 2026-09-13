using Test
using VectorDataCubes

using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import DE9IM
using Extents
using Random

# Four hand-made squares with exactly known spatial relations:
# sq2 touches sq1 along the edge x = 1, sq3 lies strictly within sq1,
# and sq4 is disjoint from everything else.
_square(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

sq1 = _square(0.0, 0.0, 1.0, 1.0)
sq2 = _square(1.0, 0.0, 2.0, 1.0)
sq3 = _square(0.25, 0.25, 0.75, 0.75)
sq4 = _square(10.0, 10.0, 11.0, 11.0)
squares = [sq1, sq2, sq3, sq4]

# Sorting keeps multi-index comparisons independent of the tree's traversal order.
selinds(gl, sel) = sort(Lookups.selectindices(gl, sel))

# `@inferred` needs concrete argument types, which globals in a test file do not have.
inferred_selectindices(gl, sel) = @inferred Lookups.selectindices(gl, sel)

@testset "selectors on hand-made squares (tree = $treedesc)" for (treedesc, treekw) in
    (("STRtree", (;)), ("nothing", (; tree = nothing)))

    gl = GeometryLookup(squares; treekw...)
    dv = rand(Geometry(gl))

    @testset "standard indices pass through" begin
        @test dv[Geometry = 2] == dv[2]
        @test dv[Geometry = 1:2] == dv[1:2]
        @test Lookups.selectindices(gl, 3) == 3
        @test Lookups.selectindices(gl, 1:2) == 1:2
    end

    @testset "Contains(point)" begin
        @test selinds(gl, Contains((0.5, 0.5))) == [1, 3]
        @test selinds(gl, Contains((1.5, 0.5))) == [2]
        @test selinds(gl, Contains((50.0, 50.0))) == Int[]
        # point-in-polygon is closed: a point on the shared border belongs to both squares
        @test selinds(gl, Contains((1.0, 0.5))) == [1, 2]
        @test selinds(gl, Contains(GI.Point(1.0, 0.5))) == [1, 2]
        @test dv[Geometry(Contains((1.5, 0.5)))] == dv[[2]]
        @test dv[Geometry = Contains((1.5, 0.5))] == dv[[2]]
        @test_throws ArgumentError Lookups.selectindices(gl, Contains(1.0))
        @test_throws ArgumentError Lookups.selectindices(gl, Contains(0 .. 1))
        # a vector of points is the sorted union
        @test Lookups.selectindices(gl, Contains([(1.5, 0.5), (0.5, 0.5)])) == [1, 2, 3]
        @test Lookups.selectindices(gl, Contains([(50.0, 50.0)])) == Int[]
        @test dv[Geometry = Contains([(1.5, 0.5), (0.5, 0.5)])] == dv[1:3]
    end

    @testset "At(geometry)" begin
        @test Lookups.selectindices(gl, At(sq2)) == 2
        # equality is geometric (GO.equals), not object identity
        @test Lookups.selectindices(gl, At(_square(1.0, 0.0, 2.0, 1.0))) == 2
        @test dv[Geometry(At(sq3))] == dv[3]
        @test_throws ArgumentError Lookups.selectindices(gl, At(_square(5.0, 5.0, 6.0, 6.0)))
        @test_throws ArgumentError Lookups.selectindices(gl, At(1.0))
        @test_throws ArgumentError Lookups.selectindices(gl, At(GI.extent(sq1)))
        @test Lookups.selectindices(gl, At([sq4, sq2])) == [4, 2]
        @test dv[Geometry = At([sq4, sq2])] == dv[[4, 2]]
    end

    @testset "Near(point)" begin
        @test Lookups.selectindices(gl, Near((10.4, 10.5))) == 4
        @test Lookups.selectindices(gl, Near((2.5, 0.5))) == 2
        @test Lookups.selectindices(gl, Near(GI.Point(2.5, 0.5))) == 2
        @test dv[Geometry(Near((2.5, 0.5)))] == dv[2]
        @test Lookups.selectindices(gl, Near([(2.5, 0.5), (10.4, 10.5)])) == [2, 4]
        @test dv[Geometry = Near([(2.5, 0.5), (10.4, 10.5)])] == dv[[2, 4]]
        # geometries covering the point are at distance zero, and ties go to the lowest
        # index, so the tree and the linear scan cannot disagree
        @test Lookups.selectindices(gl, Near((0.5, 0.5))) == 1
        @test Lookups.selectindices(gl, Near((0.3, 0.3))) == 1
        @test Lookups.selectindices(gl, Near((1.0, 0.5))) == 1
        @test Lookups.selectindices(GeometryLookup([sq3, sq1]; treekw...), Near((0.5, 0.5))) == 1
        # only point geometries are supported for Near
        @test_throws ArgumentError Lookups.selectindices(gl, Near(sq1))
        @test_throws ArgumentError Lookups.selectindices(gl, Near(1.0))
        @test_throws ArgumentError Lookups.selectindices(gl, Near((NaN, 0.5)))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(Near(Inf)), Y(Near(0.5))))
        @test_throws ArgumentError Lookups.selectindices(GeometryLookup(empty(squares); treekw...), Near((0.0, 0.0)))
    end

    @testset "(X(At), Y(At)) point lookup" begin
        @test Lookups.selectindices(gl, (X(At(1.5)), Y(At(0.5)))) == 2
        @test dv[Geometry = (X(At(1.5)), Y(At(0.5)))] == dv[2]
        @test dv[X(At(1.5)), Y(At(0.5))] == dv[2]
        # closed semantics: a border point matches (the first covering square)
        @test Lookups.selectindices(gl, (X(At(1.0)), Y(At(0.5)))) in (1, 2)
        # At requires an exact match, so points in no geometry are errors,
        # whether inside the overall extent or not
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(5.0)), Y(At(5.0))))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(50.0)), Y(At(50.0))))
    end

    @testset "(X(Contains), Y(Contains)) point lookup" begin
        @test selinds(gl, (X(Contains(0.5)), Y(Contains(0.5)))) == [1, 3]
        @test selinds(gl, (X(Contains(1.0)), Y(Contains(0.5)))) == [1, 2]
        @test selinds(gl, (X(Contains(50.0)), Y(Contains(50.0)))) == Int[]
        @test length(dv[Geometry = (X(Contains(0.5)), Y(Contains(0.5)))]) == 2
    end

    @testset "(X(Near), Y(Near)) point lookup" begin
        @test Lookups.selectindices(gl, (X(Near(2.5)), Y(Near(0.5)))) == 2
        @test dv[X(Near(10.4)), Y(Near(10.5))] == dv[4]
    end

    @testset "(X(Touches), Y(Touches)) extent lookup" begin
        @test selinds(gl, (X(Touches(0.9, 1.1)), Y(Touches(0.4, 0.6)))) == [1, 2]
        @test selinds(gl, (X(Touches(9.0, 12.0)), Y(Touches(9.0, 12.0)))) == [4]
        @test selinds(gl, (X(Touches(50.0, 60.0)), Y(Touches(50.0, 60.0)))) == Int[]
        # one-sided: the other axis is unbounded
        @test selinds(gl, (X(Touches(0.9, 1.1)),)) == [1, 2]
        @test sort(Lookups.selectindices(gl, (Touches(0.9, 1.1), nothing))) == [1, 2]
        @test dv[X(Touches(9.0, 12.0))] == dv[[4]]
    end

    @testset "(X(a .. b), Y(a .. b)) interval covers lookup" begin
        @test selinds(gl, (X(-0.1 .. 2.1), Y(-0.1 .. 1.1))) == [1, 2, 3]
        @test selinds(gl, (X(-0.1 .. 1.1), Y(-0.1 .. 1.1))) == [1, 3]
        # the interval only intersecting a geometry is not enough - it must cover it
        @test selinds(gl, (X(0.4 .. 0.6), Y(0.4 .. 0.6))) == Int[]
        @test dv[Geometry = (X(-0.1 .. 2.1), Y(-0.1 .. 1.1))] == dv[1:3]
        @test dv[X(-0.1 .. 2.1), Y(-0.1 .. 1.1)] == dv[1:3]
        # one-sided: the other axis is unbounded
        @test selinds(gl, (X(-0.1 .. 1.1),)) == [1, 3]
        @test dv[X(-0.1 .. 1.1)] == dv[[1, 3]]
        @test dv[Y(9.0 .. 12.0)] == dv[[4]]
    end

    @testset "unsupported pairs are ArgumentErrors" begin
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(0.5)), Y(Contains(0.5))))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(Contains(0.5)), Y(At(0.5))))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(0.5)),))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(Contains(0.5)),))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(Near(0.5)),))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(0.5)), Y(Touches(0.4, 0.6))))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(0 .. 1), Y(Touches(0.4, 0.6))))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(Where(x -> true)), Y(Where(x -> true))))
        @test_throws ArgumentError Lookups.selectindices(gl, (At(0.5), At(0.5), At(0.5)))
        @test_throws ArgumentError Lookups.selectindices(gl, (nothing, nothing))
        @test_throws ArgumentError dv[X(At(0.5)), Y(Contains(0.5))]
    end

    @testset "Touches(extent)" begin
        @test selinds(gl, Touches(GI.extent(sq1))) == [1, 2, 3]
        @test selinds(gl, Touches(GI.extent(sq4))) == [4]
        @test selinds(gl, Touches(Extent(X = (50.0, 60.0), Y = (50.0, 60.0)))) == Int[]
        @test dv[Geometry = Touches(GI.extent(sq4))] == dv[[4]]
        # the interval form is per axis only, and a 2-tuple is not read as a point
        @test_throws ArgumentError Lookups.selectindices(gl, Touches(0.4, 0.6))
        @test_throws ArgumentError Lookups.selectindices(gl, Touches((0.9, 1.1), (0.4, 0.6)))
        @test_throws ArgumentError Lookups.selectindices(gl, Touches([(0.4, 0.6)]))
        @test_throws ArgumentError Lookups.selectindices(gl, Touches(nothing))
    end

    @testset "Where with GeometryOps predicates" begin
        @test selinds(gl, Where(GO.intersects(sq1))) == [1, 2, 3]
        # GO.equals has no curried form, so build the Fix2 directly
        @test selinds(gl, Where(Base.Fix2(GO.equals, sq2))) == [2]
        @test selinds(gl, Where(GO.contains(sq3))) == [1, 3]
        @test selinds(gl, Where(GO.within(sq1))) == [1, 3]
        @test selinds(gl, Where(GO.covers(sq3))) == [1, 3]
        @test selinds(gl, Where(GO.coveredby(sq1))) == [1, 3]
        @test selinds(gl, Where(GO.touches(sq1))) == [2]
        @test selinds(gl, Where(GO.disjoint(sq1))) == [4]
        # a line entering sq1 below sq3, ending inside sq1: crosses sq1 only
        @test selinds(gl, Where(GO.crosses(GI.LineString([(-1.0, 0.1), (0.5, 0.1)])))) == [1]
        # a square offset by 0.5 partly overlaps sq1, sq2 and sq3
        @test selinds(gl, Where(GO.overlaps(_square(0.5, 0.5, 1.5, 1.5)))) == [1, 2, 3]
        @test dv[Geometry = Where(GO.disjoint(sq1))] == dv[[4]]
    end

    @testset "Where(GO.pred(g)) dispatches to the tree-narrowed fast path" begin
        for f in (GO.intersects, GO.contains, GO.within, GO.covers, GO.coveredby,
                  GO.touches, GO.crosses, GO.overlaps, GO.disjoint)
            sel = Where(f(sq1))
            m = which(Lookups.selectindices, Tuple{typeof(gl), typeof(sel)})
            @test m.module == VectorDataCubes
            @test occursin("Fix2", string(m.sig))
        end
        m = which(Lookups.selectindices, Tuple{typeof(gl), typeof(Where(Base.Fix2(GO.equals, sq1)))})
        @test occursin("Fix2", string(m.sig))
    end

    @testset "generic Where" begin
        @test Lookups.selectindices(gl, Where(g -> GO.area(g) < 0.5)) == [3]
        @test Lookups.selectindices(gl, Where(x -> false)) == Int[]
        @test Lookups.selectindices(gl, Where(x -> false)) isa Vector{Int}
        empty_dv = dv[Geometry = Where(x -> false)]
        @test isempty(empty_dv)
        @test empty_dv isa DD.AbstractDimVector
    end

    @testset "DE9IM predicates" begin
        @test selinds(gl, DE9IM.Intersects(sq1)) == [1, 2, 3]
        @test selinds(gl, DE9IM.Disjoint(sq1)) == [4]
        @test selinds(gl, DE9IM.Contains(sq3)) == [1, 3]
        @test selinds(gl, DE9IM.Within(sq1)) == [1, 3]
        @test selinds(gl, DE9IM.Covers(sq3)) == [1, 3]
        @test selinds(gl, DE9IM.CoveredBy(sq1)) == [1, 3]
        @test selinds(gl, DE9IM.Touches(sq1)) == [2]
        @test selinds(gl, DE9IM.Crosses(GI.LineString([(-1.0, 0.1), (0.5, 0.1)]))) == [1]
        @test selinds(gl, DE9IM.Overlaps(_square(0.5, 0.5, 1.5, 1.5))) == [1, 2, 3]
        @test selinds(gl, DE9IM.Equals(sq2)) == [2]
        @test selinds(gl, DE9IM.Equals(_square(5.0, 5.0, 6.0, 6.0))) == Int[]
        @test dv[Geometry = DE9IM.Disjoint(sq1)] == dv[[4]]
        @test dv[Geometry(DE9IM.Touches(sq1))] == dv[[2]]
        # vector-valued: sorted union
        @test Lookups.selectindices(gl, DE9IM.Intersects([sq4, sq2])) == [1, 2, 4]
        @test Lookups.selectindices(gl, DE9IM.Equals([sq4, sq2])) == [2, 4]
        @test dv[Geometry = DE9IM.Equals([sq4, sq2])] == dv[[2, 4]]
        # the wrapped value must be a geometry (or a vector of them), never an extent,
        # nothing or a keyword-carrying predicate
        @test_throws ArgumentError Lookups.selectindices(gl, DE9IM.Intersects(sq1; foo = 1))
        @test_throws ArgumentError Lookups.selectindices(gl, DE9IM.Intersects())
        @test_throws ArgumentError Lookups.selectindices(gl, DE9IM.Intersects(1.0))
        @test_throws ArgumentError Lookups.selectindices(gl, DE9IM.Intersects(GI.extent(sq1)))
        @test_throws ArgumentError Lookups.selectindices(gl, DE9IM.Intersects(Any[sq1, 1.0]))
    end

    @testset "hasselection" begin
        @test Lookups.hasselection(gl, At(sq2))
        @test !Lookups.hasselection(gl, At(_square(5.0, 5.0, 6.0, 6.0)))
        @test Lookups.hasselection(gl, Contains((0.5, 0.5)))
        @test !Lookups.hasselection(gl, Contains((50.0, 50.0)))
        @test Lookups.hasselection(gl, Near((50.0, 50.0)))
        @test !Lookups.hasselection(gl, Near(sq1))
        @test !Lookups.hasselection(GeometryLookup(empty(squares); treekw...), Near((0.0, 0.0)))
        @test Lookups.hasselection(gl, Touches(GI.extent(sq4)))
        @test !Lookups.hasselection(gl, Touches(Extent(X = (50.0, 60.0), Y = (50.0, 60.0))))
        @test Lookups.hasselection(gl, Where(GO.disjoint(sq1)))
        @test Lookups.hasselection(gl, Where(g -> GO.area(g) < 0.5))
        @test !Lookups.hasselection(gl, Where(x -> false))
        @test Lookups.hasselection(gl, DE9IM.Touches(sq1))
        @test !Lookups.hasselection(gl, DE9IM.Touches(sq4))
        # vector-valued selectors: every value must select
        @test Lookups.hasselection(gl, At([sq1, sq2]))
        @test !Lookups.hasselection(gl, At([sq1, _square(5.0, 5.0, 6.0, 6.0)]))
        @test Lookups.hasselection(gl, Near([(2.5, 0.5)]))
        @test !Lookups.hasselection(gl, Near([sq1]))
        @test Lookups.hasselection(gl, Contains([(0.5, 0.5)]))
        @test !Lookups.hasselection(gl, Contains([(50.0, 50.0)]))
        # a value this lookup cannot use is `false`, not an error: indexing would fail
        @test !Lookups.hasselection(gl, At(1.0))
        @test !Lookups.hasselection(gl, Contains(1.0))
        @test !Lookups.hasselection(gl, Touches(0.4, 0.6))
        @test !Lookups.hasselection(gl, DE9IM.Intersects(1.0))
        @test !Lookups.hasselection(gl, DE9IM.Intersects(sq1; foo = 1))
        # through the cube and its dimension
        @test Lookups.hasselection(dv, Geometry(Contains((0.5, 0.5))))
        @test !Lookups.hasselection(dv, Geometry(Contains((50.0, 50.0))))
        @test Lookups.hasselection(dv, Geometry(DE9IM.Touches(sq1)))
        @test !Lookups.hasselection(dv, Geometry(DE9IM.Touches(sq4)))
    end

    @testset "selectindices is type stable" begin
        @test inferred_selectindices(gl, Contains((0.5, 0.5))) == [1, 3]
        @test inferred_selectindices(gl, Contains([(0.5, 0.5)])) == [1, 3]
        @test inferred_selectindices(gl, At(sq2)) == 2
        @test inferred_selectindices(gl, At([sq2])) == [2]
        @test inferred_selectindices(gl, Near((2.5, 0.5))) == 2
        @test inferred_selectindices(gl, Near([(2.5, 0.5)])) == [2]
        @test inferred_selectindices(gl, Touches(GI.extent(sq4))) == [4]
        @test inferred_selectindices(gl, Where(GO.disjoint(sq1))) == [4]
        @test inferred_selectindices(gl, Where(g -> GO.area(g) < 0.5)) == [3]
        @test inferred_selectindices(gl, DE9IM.Touches(sq1)) == [2]
        @test inferred_selectindices(gl, DE9IM.Touches([sq1])) == [2]
        @test inferred_selectindices(gl, (X(At(1.5)), Y(At(0.5)))) == 2
        @test inferred_selectindices(gl, (X(Near(2.5)), Y(Near(0.5)))) == 2
        @test sort(inferred_selectindices(gl, (X(Contains(0.5)), Y(Contains(0.5))))) == [1, 3]
        @test sort(inferred_selectindices(gl, (X(-0.1 .. 1.1), Y(-0.1 .. 1.1)))) == [1, 3]
        @test sort(inferred_selectindices(gl, (X(Touches(0.9, 1.1)), Y(Touches(0.4, 0.6))))) == [1, 2]
    end

    @testset "mask" begin
        @test VectorDataCubes.mask(gl, Contains((0.5, 0.5))) == [true, false, true, false]
        @test VectorDataCubes.mask(gl, At(sq2)) == [false, true, false, false]
        @test VectorDataCubes.mask(gl, Near((2.5, 0.5))) == [false, true, false, false]
        @test VectorDataCubes.mask(gl, DE9IM.Disjoint(sq1)) == [false, false, false, true]
        @test VectorDataCubes.mask(gl, (X(-0.1 .. 1.1), Y(-0.1 .. 1.1))) == [true, false, true, false]
        @test VectorDataCubes.mask(gl, Where(x -> false)) == falses(4)
        m = VectorDataCubes.mask(dv, Contains((0.5, 0.5)))
        @test m isa DD.DimArray{Bool}
        @test DD.dims(m) == DD.dims(dv)
        @test parent(m) == [true, false, true, false]
        @test dv[Geometry = m] == dv[[1, 3]]
        cube = rand(Geometry(gl), Ti(1:3))
        @test parent(VectorDataCubes.mask(cube, At(sq4))) == [false, false, false, true]
        @test parent(VectorDataCubes.mask(DimStack((a = cube,)), At(sq4))) == [false, false, false, true]
        @test_throws ArgumentError VectorDataCubes.mask(rand(X(1:3), Y(1:3)), At(sq4))
        @test_throws ArgumentError VectorDataCubes.mask(rand(Geometry(gl), Dim{:Origin}(gl)), At(sq4))
    end

    @testset "empty selections produce empty arrays" begin
        empty_dv = dv[Geometry(Contains((50.0, 50.0)))]
        @test isempty(empty_dv)
        @test empty_dv isa DD.AbstractDimVector
    end
end

@testset "(Y(), X()) lookup answers like (X(), Y())" begin
    xy = GeometryLookup(squares, (X(), Y()))
    yx = GeometryLookup(squares, (Y(), X()))
    pairs = (
        (X(At(1.5)), Y(At(0.5))),
        (X(Contains(0.5)), Y(Contains(0.5))),
        (X(Near(2.5)), Y(Near(0.5))),
        (X(-0.1 .. 2.1), Y(-0.1 .. 1.1)),
        (X(Touches(0.9, 1.1)), Y(Touches(0.4, 0.6))),
        (X(-0.1 .. 1.1),),
        (Y(9.0 .. 12.0),),
        (X(Touches(0.9, 1.1)),),
    )
    for sel in pairs
        @test sort(vec([Lookups.selectindices(xy, sel);])) == sort(vec([Lookups.selectindices(yx, sel);]))
    end
    # the pair is matched by dimension type, not by position
    @test Lookups.selectindices(yx, (X(At(1.5)), Y(At(0.5)))) == 2
    @test_throws ArgumentError Lookups.selectindices(yx, (X(At(0.5)), Y(At(1.5))))
    @test sort(Lookups.selectindices(yx, (X(Touches(9.0, 12.0)), Y(Touches(-1.0, 2.0))))) == Int[]
    @test sort(Lookups.selectindices(yx, (X(-1.0 .. 3.0), Y(-1.0 .. 2.0)))) == [1, 2, 3]
    dv = rand(Geometry(yx))
    @test dv[X(At(1.5)), Y(At(0.5))] == dv[2]
    @test dv[X(-0.1 .. 1.1), Y(-0.1 .. 1.1)] == dv[[1, 3]]
    @test dv[Geometry = (X(-0.1 .. 1.1), Y(-0.1 .. 1.1))] == dv[[1, 3]]
    @test dv[X(-0.1 .. 1.1)] == dv[[1, 3]]
    @test dv[Y(9.0 .. 12.0)] == dv[[4]]
    @test dv[X(Touches(9.0, 12.0))] == dv[[4]]
end

@testset "Near: tree search agrees with the linear scan" begin
    rng = MersenneTwister(1234)
    # Squares large enough to overlap heavily, so most query points sit inside several
    # of them and the zero-distance ties are exercised.
    polys = map(1:500) do _
        x, y = 30 .* rand(rng, 2)
        s = 1.0 + 9.0 * rand(rng)
        _square(x, y, x + s, y + s)
    end
    withtree = GeometryLookup(polys)
    notree = GeometryLookup(polys; tree = nothing)
    for _ in 1:50
        p = (-5.0 + 45.0 * rand(rng), -5.0 + 45.0 * rand(rng))
        linear = argmin(i -> GO.distance(p, polys[i]), eachindex(polys))
        @test Lookups.selectindices(notree, Near(p)) == linear
        @test Lookups.selectindices(withtree, Near(p)) == linear
    end
end
