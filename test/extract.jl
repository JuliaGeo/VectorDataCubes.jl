using Test
using VectorDataCubes

using Rasters, DimensionalData
using Rasters.Lookups
# Explicitly bind VectorDataCubes' own `extract` (Rasters exports one too).
using VectorDataCubes: extract
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
using Dates

# The 10x10 raster of test/zonal.jl (value = x index), but with `Intervals`
# sampling so a point anywhere inside a cell selects it.
ex = X(Sampled(0.5:1.0:9.5; sampling=Intervals(Center())))
ey = Y(Sampled(0.5:1.0:9.5; sampling=Intervals(Center())))
eti = Ti(DateTime(2020, 1, 1):Month(1):DateTime(2020, 3, 1))
ras2d = Raster([Float64(xi) for xi in 1:10, yi in 1:10], (ex, ey); name=:vals)
ras3d = Raster([Float64(xi * t) for xi in 1:10, yi in 1:10, t in 1:3], (ex, ey, eti); name=:cube)

# x cells 2 and 6, and a point outside the raster
pts = [(1.2, 1.3), (5.5, 2.5), (20.0, 20.0)]

@testset "extract at points" begin
    @testset "2D raster -> vector over Geometry" begin
        res = extract(ras2d, pts)
        @test res isa Raster
        @test DD.name(res) == :vals
        @test DD.lookup(res, Geometry) isa GeometryLookup
        @test parent(DD.lookup(res, Geometry)) == pts
        @test res[1] == 2.0
        @test res[2] == 6.0
        @test ismissing(res[3])
        @test res[Geometry(Near((5.4, 2.4)))] == 6.0
    end

    @testset "3D raster -> cube over (Ti, Geometry)" begin
        res = extract(ras3d, pts)
        @test size(res) == (3, 3)
        @test map(DD.name, DD.dims(res)) == (:Ti, :Geometry)
        @test DD.dims(res, Ti) == DD.dims(ras3d, Ti)
        @test res[Ti=2, Geometry=1] == 4.0
        @test res[Ti=3, Geometry=2] == 18.0
        @test all(ismissing, res[Geometry=3])
        # every point outside the raster still keeps the cube shape
        off = extract(ras3d, [(20.0, 20.0), (30.0, 30.0)])
        @test size(off) == (3, 2)
        @test DD.dims(off, Ti) == DD.dims(ras3d, Ti)
        @test all(ismissing, off)
    end

    @testset "RasterStack -> stack of cubes" begin
        flat = rebuild(ras2d; metadata=Dict{Symbol,Any}(:units => "K"))
        st = RasterStack((flat=flat, cube=ras3d); metadata=Dict{Symbol,Any}(:source => "test"))
        res = extract(st, pts)
        @test res isa RasterStack
        @test DD.metadata(res)[:source] == "test"
        @test DD.metadata(res[:flat])[:units] == "K"
        @test size(res[:flat]) == (3,)
        @test size(res[:cube]) == (3, 3)
        @test res[:flat][2] == 6.0
        @test res[:cube][Ti=1, Geometry=2] == 6.0
        @test ismissing(res[:flat][3])
    end

    @testset "skipmissing drops points outside the raster" begin
        res = extract(ras3d, pts; skipmissing=true)
        @test size(res) == (3, 2)
        @test parent(DD.lookup(res, Geometry)) == pts[1:2]
        @test !any(ismissing, res)
        res2d = extract(ras2d, pts; skipmissing=true)
        @test collect(res2d) == [2.0, 6.0]
        # a stack drops a point when any layer is missing there
        st = RasterStack((flat=ras2d, cube=ras3d))
        rst = extract(st, pts; skipmissing=true)
        @test size(rst[:cube]) == (3, 2)
        gappy = [Float64(yi) for xi in 1:10, yi in 1:10]
        gappy[2, 2] = -9999.0  # the cell of pts[1], in one layer only
        st2 = RasterStack((a=ras2d, b=Raster(gappy, (ex, ey); name=:b, missingval=-9999.0)))
        rst2 = extract(st2, pts; skipmissing=true)
        @test parent(DD.lookup(rst2[:a], Geometry)) == pts[2:2]
        @test collect(rst2[:a]) == [6.0]
        # nothing left keeps the cube's shape
        @test size(extract(ras3d, [(20.0, 20.0)]; skipmissing=true)) == (3, 0)
        # a sentinel missingval counts as missing
        data = [Float64(xi) for xi in 1:10, yi in 1:10]
        data[1, 1] = -9999.0
        rasmv = Raster(data, (ex, ey); name=:vals, missingval=-9999.0)
        @test collect(extract(rasmv, [(0.5, 0.5), (1.5, 1.5)]; skipmissing=true)) == [2.0]
        @test isequal(collect(extract(rasmv, [(0.5, 0.5), (1.5, 1.5)])), [-9999.0, 2.0])
    end

    @testset "skipmissing drops points in an all-missing cell, missing or sentinel" begin
        # pts[1] sits in cell (2, 2), pts[2] in cell (6, 3)
        withmissing = Array{Union{Missing,Float64}}(parent(ras2d))
        withmissing[2, 2] = missing
        sentinel = copy(parent(ras2d))
        sentinel[2, 2] = -9999.0
        for r in (Raster(withmissing, (ex, ey)), Raster(sentinel, (ex, ey); missingval=-9999.0))
            res = extract(r, pts[1:2]; skipmissing=true)
            @test collect(res) == [6.0]
            @test parent(DD.lookup(res, Geometry)) == pts[2:2]
        end
        # on a 3-D raster the cell is a slice, and only an all-missing slice drops the point
        withmissing3 = Array{Union{Missing,Float64}}(parent(ras3d))
        withmissing3[2, 2, :] .= missing
        withmissing3[6, 3, 1] = missing
        resm = extract(Raster(withmissing3, (ex, ey, eti)), pts[1:2]; skipmissing=true)
        @test parent(DD.lookup(resm, Geometry)) == pts[2:2]
        @test isequal(vec(collect(resm)), [missing, 12.0, 18.0])
        sentinel3 = copy(parent(ras3d))
        sentinel3[2, 2, :] .= -9999.0
        sentinel3[6, 3, 1] = -9999.0
        ress = extract(Raster(sentinel3, (ex, ey, eti); missingval=-9999.0), pts[1:2]; skipmissing=true)
        @test parent(DD.lookup(ress, Geometry)) == pts[2:2]
        @test vec(collect(ress)) == [-9999.0, 12.0, 18.0]
    end

    @testset "a reverse-ordered Y axis gives the same answer" begin
        # value = y index, so a y axis read the wrong way round would show
        yfwd = Raster([Float64(yi) for xi in 1:10, yi in 1:10], (ex, ey); name=:yv)
        yrev = Raster([Float64(yi) for xi in 1:10, yi in 10:-1:1],
            (ex, Y(Sampled(9.5:-1.0:0.5; sampling=Intervals(Center())))); name=:yv)
        @test isequal(collect(extract(yfwd, pts)), collect(extract(yrev, pts)))
        @test collect(extract(yrev, pts[1:2])) == [2.0, 3.0]
    end

    @testset "table input and crs" begin
        tbl = (; geometry=pts, name=["a", "b", "c"])
        res = extract(ras2d, tbl; crs=EPSG(4326))
        @test crs(DD.lookup(res, Geometry)) == EPSG(4326)
        @test collect(res[1:2]) == [2.0, 6.0]
        # the lookup's crs is kept, and the keyword overrides it
        gl = GeometryLookup(pts; crs=EPSG(4326))
        @test crs(DD.lookup(extract(ras2d, gl), Geometry)) == EPSG(4326)
        @test crs(DD.lookup(extract(ras2d, gl; crs=EPSG(3857)), Geometry)) == EPSG(3857)
        # mismatched crs warns, like zonal
        ras3857 = Rasters.setcrs(ras2d, EPSG(3857))
        @test_logs (:warn, r"crs") match_mode=:any extract(ras3857, gl)
    end

    @testset "dimension input keeps its name" begin
        res = extract(ras3d, Dim{:Station}(GeometryLookup(pts)))
        @test map(DD.name, DD.dims(res)) == (:Ti, :Station)
        @test res[Ti=1, Station=2] == 6.0
    end

    @testset "points on a cell edge and just outside" begin
        # Intervals(Center) cells run [0, 1), [1, 2), ...: an edge point is the upper cell's
        @test collect(extract(ras2d, [(1.0, 1.0), (0.0, 0.0), (9.999, 9.999)])) == [2.0, 1.0, 10.0]
        @test all(ismissing, extract(ras2d, [(10.0, 5.0), (-0.001, 5.0), (5.0, 10.0)]))
    end

    @testset "Points sampling uses At with atol" begin
        rasp = Raster(parent(ras2d), (X(0.5:1.0:9.5), Y(0.5:1.0:9.5)); name=:vals)
        @test isequal(collect(extract(rasp, [(1.5, 1.5), (1.51, 1.5)])), [2.0, missing])
        @test collect(extract(rasp, [(1.5, 1.5), (1.51, 1.5)]; atol=0.1)) == [2.0, 2.0]
    end

    @testset "non-point geometries and empty input throw" begin
        square = GI.Polygon([GI.LinearRing([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)])])
        @test_throws ArgumentError extract(ras2d, [square])
        @test_throws ArgumentError extract(ras2d, GeometryLookup([square]))
        @test_throws ArgumentError extract(ras2d, DD.rebuild(GeometryLookup(pts); data=empty(pts)))
        # the lookup's dims must be on the raster
        @test_throws ArgumentError extract(Raster(rand(10, 3), (ex, eti)), pts)
    end
end
