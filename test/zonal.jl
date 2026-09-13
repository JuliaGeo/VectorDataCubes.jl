using Test
using VectorDataCubes

using Rasters, DimensionalData
using Rasters.Lookups
# Explicitly bind VectorDataCubes' own `zonal` (Rasters exports one too).
using VectorDataCubes: zonal
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import Proj # the "reproject" testset needs the Proj extension
using Statistics: mean
using Dates

_zsquare(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

# A 10x10 raster on [0, 10]^2 with cell centers at 0.5, 1.5, ..., 9.5,
# where each cell's value is its x index. Zone values are exactly known:
# zoneA covers x cells 1:2, zoneB covers x cells 5:7, zoneC is off-raster.
zx = X(0.5:1.0:9.5)
zy = Y(0.5:1.0:9.5)
zti = Ti(DateTime(2020, 1, 1):Month(1):DateTime(2020, 3, 1))
ras2d = Raster([Float64(xi) for xi in 1:10, yi in 1:10], (zx, zy); name=:vals)
# 3D: value = x index * time index
ras3d = Raster([Float64(xi * t) for xi in 1:10, yi in 1:10, t in 1:3], (zx, zy, zti); name=:cube)

zoneA = _zsquare(0.0, 0.0, 2.0, 2.0)    # 2x2 cells, x values {1, 2}
zoneB = _zsquare(4.0, 2.0, 7.0, 5.0)    # 3x3 cells, x values {5, 6, 7}
zoneC = _zsquare(20.0, 20.0, 22.0, 22.0) # entirely outside the raster
zgl = GeometryLookup([zoneA, zoneB, zoneC])

@testset "zonal with a GeometryLookup" begin
    @testset "2D raster -> vector over Geometry" begin
        res = zonal(sum, ras2d; of=zgl, progress=false)
        @test res isa Raster
        @test DD.dims(res, Geometry) isa Geometry
        @test val(DD.lookup(res, Geometry)) == parent(zgl)
        @test res[1] == (1 + 2) * 2        # {1,2} over 2 y-cells
        @test res[2] == (5 + 6 + 7) * 3    # {5,6,7} over 3 y-cells
        @test ismissing(res[3])            # off-raster geometry
        # the result is a real vector data cube: spatial selectors work on it
        @test res[Geometry(Contains((1.0, 1.0)))] == res[[1]]
    end

    @testset "of = Geometry(lookup) behaves the same" begin
        res_lookup = zonal(sum, ras2d; of=zgl, progress=false)
        res_dim = zonal(sum, ras2d; of=Geometry(zgl), progress=false)
        @test isequal(res_lookup, res_dim)
    end

    @testset "non-lookup `of` forwards to Rasters.zonal" begin
        res = zonal(sum, ras2d; of=[zoneA, zoneB], progress=false)
        @test res isa Vector
        @test res == [(1 + 2) * 2, (5 + 6 + 7) * 3]
    end

    @testset "3D raster -> cube over (Ti, Geometry)" begin
        res = zonal(mean, ras3d; of=zgl, progress=false)
        @test size(res) == (3, 3)
        @test DD.dims(res, Ti) == DD.dims(ras3d, Ti)
        @test DD.dims(res, Geometry) isa Geometry
        @test res[Ti=1, Geometry=1] ≈ 1.5    # mean({1,2} * 1)
        @test res[Ti=2, Geometry=1] ≈ 3.0    # mean({1,2} * 2)
        @test res[Ti=3, Geometry=2] ≈ 18.0   # mean({5,6,7} * 3)
        @test all(ismissing, res[Geometry=3])
        @test !any(ismissing, res[Geometry=1:2])
    end

    @testset "a reverse-ordered Y axis gives the same answer" begin
        # value = y index, so a y axis read the wrong way round would show
        yfwd = Raster([Float64(yi) for xi in 1:10, yi in 1:10], (zx, zy); name=:yv)
        yrev = Raster([Float64(yi) for xi in 1:10, yi in 10:-1:1], (zx, Y(9.5:-1.0:0.5)); name=:yv)
        @test isequal(zonal(sum, yfwd; of=zgl, progress=false), zonal(sum, yrev; of=zgl, progress=false))
        @test zonal(sum, yrev; of=zgl, progress=false)[2] == (3 + 4 + 5) * 3
    end

    @testset "spatialslices = false reduces over all dims" begin
        res = zonal(mean, ras3d; of=zgl, spatialslices=false, progress=false)
        @test size(res) == (3,)
        @test res[1] ≈ mean([xi * t for xi in 1:2, _ in 1:2, t in 1:3]) # 3.0
        @test ismissing(res[3])
    end

    @testset "RasterStack -> stack of cubes" begin
        st = RasterStack((flat=ras2d, cube=ras3d))
        res = zonal(mean, st; of=zgl, progress=false)
        @test res isa RasterStack
        @test size(res[:flat]) == (3,)
        @test size(res[:cube]) == (3, 3)
        @test res[:flat][1] ≈ 1.5
        @test res[:cube][Ti=2, Geometry=2] ≈ 12.0  # mean({5,6,7} * 2)
        @test ismissing(res[:flat][3])
        @test all(ismissing, res[:cube][Geometry=3])
    end

    @testset "emptyval fills empty slices" begin
        # data missing at t = 2 in zoneD's cells: that slice is empty under
        # skipmissing, so it gets emptyval while other slices are computed
        data = Array{Union{Missing,Float64}}([Float64(xi * t) for xi in 1:10, yi in 1:10, t in 1:3])
        data[9:10, 9:10, 2] .= missing
        rasm = Raster(data, (zx, zy, zti); name=:gappy)
        zoneD = _zsquare(8.0, 8.0, 10.0, 10.0) # x cells 9:10, y cells 9:10
        # -1 rather than NaN: `mean` of an empty iterator is NaN either way
        res = zonal(mean, rasm; of=GeometryLookup([zoneD]), emptyval=-1.0, progress=false)
        @test res[Ti=2, Geometry=1] == -1.0
        @test res[Ti=1, Geometry=1] ≈ 9.5
        @test res[Ti=3, Geometry=1] ≈ 28.5
    end

    @testset "emptyval with a masked-out and an off-raster geometry" begin
        # zoneE's bounding box holds the cell center (0.5, 0.5) but the triangle does not:
        # the crop is non-empty yet the mask removes every cell, so each slice -> emptyval.
        # zoneC crops to nothing -> missing. Eltypes then differ between geometries.
        zoneE = GI.Polygon([GI.LinearRing([(0.0, 0.0), (0.9, 0.0), (0.0, 0.9), (0.0, 0.0)])])
        gl = GeometryLookup([zoneA, zoneE, zoneC])
        res = zonal(mean, ras3d; of=gl, emptyval=-1.0, progress=false)
        @test size(res) == (3, 3)
        @test !any(ismissing, res[Geometry=1])
        @test res[Ti=1, Geometry=1] ≈ 1.5
        @test all(==(-1.0), res[Geometry=2])   # masked out -> emptyval per slice
        @test all(ismissing, res[Geometry=3])  # off-raster -> missing
        resm = zonal(mean, ras3d; of=gl, emptyval=missing, progress=false)
        @test all(ismissing, resm[Geometry=2:3])
        @test resm[Ti=2, Geometry=1] ≈ 3.0
    end

    @testset "sentinel missingval" begin
        data2d = [Float64(xi) for xi in 1:10, yi in 1:10]
        data2d[1, 1] = -9999.0
        rasmv = Raster(data2d, (zx, zy); name=:vals, missingval=-9999.0)
        res = zonal(sum, rasmv; of=zgl, progress=false)
        @test res[1] == (1 + 2) * 2 - 1
        @test res[2] == (5 + 6 + 7) * 3
        @test ismissing(res[3])
        data3d = [Float64(xi * t) for xi in 1:10, yi in 1:10, t in 1:3]
        data3d[1:2, 1:2, 2] .= -9999.0
        rasmv3 = Raster(data3d, (zx, zy, zti); name=:cube, missingval=-9999.0)
        res3 = zonal(sum, rasmv3; of=zgl, emptyval=NaN, progress=false)
        @test res3[Ti=1, Geometry=1] == 6.0
        @test isnan(res3[Ti=2, Geometry=1])
        @test res3[Ti=3, Geometry=1] == 18.0
        @test res3[Ti=2, Geometry=2] == 108.0
    end

    @testset "the geometry dimension keeps its name" begin
        res = zonal(sum, ras2d; of=Dim{:Origin}(zgl), progress=false)
        @test DD.hasdim(res, Dim{:Origin})
        @test !DD.hasdim(res, Geometry)
        @test DD.lookup(res, Dim{:Origin}) isa GeometryLookup
        @test res[Origin=1] == (1 + 2) * 2
        res3 = zonal(sum, ras3d; of=Dim{:Origin}(zgl), progress=false)
        @test map(DD.name, DD.dims(res3)) == (:Ti, :Origin)
    end

    @testset "spatialslices as a tuple" begin
        res = zonal(mean, ras3d; of=zgl, spatialslices=(X, Y), progress=false)
        @test isequal(res, zonal(mean, ras3d; of=zgl, progress=false))
        res_all = zonal(mean, ras3d; of=zgl, spatialslices=(X, Y, Ti), progress=false)
        @test isequal(res_all, zonal(mean, ras3d; of=zgl, spatialslices=false, progress=false))
        @test res_all[1] ≈ 3.0
        # per-geometry crops differ in spatial size, so the lookup's dims can never be left out
        @test_throws ArgumentError zonal(mean, ras3d; of=zgl, spatialslices=(Ti,), progress=false)
        @test_throws ArgumentError zonal(mean, ras3d; of=zgl, spatialslices=(X, Ti), progress=false)
        @test_throws ArgumentError zonal(mean, ras3d; of=zgl, spatialslices=(X, Y, Band), progress=false)
        @test_throws ArgumentError zonal(mean, ras3d; of=zgl, spatialslices=:all, progress=false)
    end

    @testset "raster without the lookup's dims throws" begin
        rasxt = Raster(rand(10, 3), (zx, zti); name=:noy)
        @test_throws ArgumentError zonal(sum, rasxt; of=zgl, progress=false)
    end

    @testset "of as a cube, a stack or a dim tuple" begin
        cube = zonal(sum, ras2d; of=zgl, progress=false)
        for of in (cube, RasterStack((a=cube,)), (Geometry(zgl),), (zti, Geometry(zgl)))
            res = zonal(sum, ras2d; of, progress=false)
            @test isequal(res, cube)
        end
        # the geometry dim's name is taken from `of`
        res = zonal(sum, ras2d; of=(zti, Dim{:Origin}(zgl)), progress=false)
        @test DD.hasdim(res, Dim{:Origin})
        # a raster without a geometry dim is an extent for Rasters.zonal
        @test zonal(sum, ras2d; of=ras2d) == sum(ras2d)
        # two geometry dims is ambiguous
        od = Raster(rand(3, 3), (Geometry(zgl), Dim{:Destination}(zgl)))
        @test_throws ArgumentError zonal(sum, ras2d; of=od, progress=false)
    end

    @testset "crs mismatch warns" begin
        ras3857 = Rasters.setcrs(ras2d, EPSG(3857))
        gl4326 = GeometryLookup([zoneA, zoneB]; crs=EPSG(4326))
        gl3857 = GeometryLookup([zoneA, zoneB]; crs=EPSG(3857))
        @test_logs (:warn, r"crs") match_mode=:any zonal(sum, ras3857; of=gl4326, progress=false)
        @test_logs min_level=Base.CoreLogging.Warn zonal(sum, ras3857; of=gl3857, progress=false)
        @test_logs min_level=Base.CoreLogging.Warn zonal(sum, ras2d; of=gl4326, progress=false)
        # one warning per call, not one per layer
        st3857 = Rasters.setcrs(RasterStack((a=ras2d, b=ras2d, c=ras2d)), EPSG(3857))
        @test_logs (:warn, r"crs") zonal(sum, st3857; of=gl4326, progress=false)
    end

    @testset "name and metadata pass through" begin
        rasmeta = Raster(parent(ras2d), (zx, zy); name=:vals, metadata=Dict{Symbol,Any}(:units => "K"))
        res = zonal(sum, rasmeta; of=zgl, progress=false)
        @test DD.name(res) == :vals
        @test DD.metadata(res)[:units] == "K"
        st = RasterStack((flat=ras2d, cube=ras3d); metadata=Dict{Symbol,Any}(:source => "test"))
        @test DD.metadata(zonal(sum, st; of=zgl, progress=false))[:source] == "test"
    end

    @testset "all geometries off-raster keeps the cube shape" begin
        zoneC2 = _zsquare(30.0, 30.0, 32.0, 32.0)
        off_gl = GeometryLookup([zoneC, zoneC2])
        res = zonal(mean, ras3d; of=off_gl, progress=false)
        @test size(res) == (3, 2)  # (Ti, Geometry), not collapsed to (2,)
        @test DD.dims(res, Ti) == DD.dims(ras3d, Ti)
        @test all(ismissing, res)
        res2d = zonal(mean, ras2d; of=off_gl, progress=false)
        @test size(res2d) == (2,)
        @test all(ismissing, res2d)
    end

    @testset "empty lookup throws" begin
        empty_gl = DD.rebuild(zgl; data=empty(parent(zgl)))
        @test_throws ArgumentError zonal(sum, ras2d; of=empty_gl, progress=false)
    end
end

@testset "reproject" begin
    gl = GeometryLookup([zoneA, zoneB]; crs=EPSG(4326))
    gl3857 = reproject(EPSG(3857), gl)
    @test crs(gl3857) == EPSG(3857)
    # Web Mercator coordinates are in meters, so points are far from the origin
    @test GI.x(GI.getpoint(GI.getexterior(gl3857[1]), 2)) ≈ 2.0 * 20037508.34 / 180 rtol = 1e-3
    # round trip back to EPSG:4326 recovers the original coordinates (up to float error)
    glback = reproject(EPSG(4326), gl3857)
    @test all(zip(val(glback), val(gl))) do (a, b)
        all(zip(GI.getpoint(a), GI.getpoint(b))) do (p, q)
            isapprox(GI.x(p), GI.x(q); atol=1e-6) && isapprox(GI.y(p), GI.y(q); atol=1e-6)
        end
    end
    # reprojecting a whole vector data cube reprojects the lookup
    dv = rand(Geometry(gl))
    dv3857 = reproject(EPSG(3857), dv)
    @test crs(DD.lookup(dv3857, Geometry)) == EPSG(3857)
    # no crs is an error
    @test_throws ArgumentError reproject(EPSG(3857), GeometryLookup([zoneA]))
end
