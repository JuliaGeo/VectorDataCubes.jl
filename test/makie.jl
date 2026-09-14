using Test
using VectorDataCubes

using Makie
import Makie.GeometryBasics as GB
import GeoInterface as GI
import DimensionalData as DD

const EXT = Base.get_extension(VectorDataCubes, :VectorDataCubesMakieExt)

_msquare(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

msq1 = _msquare(0.0, 0.0, 1.0, 1.0)
msq2 = _msquare(2.0, 0.0, 3.0, 1.0)
msq3 = _msquare(10.0, 10.0, 11.0, 11.0)
msquares = [msq1, msq2, msq3]
multi = GI.MultiPolygon([msq1, msq3])
line1 = GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])
line2 = GI.LineString([(3.0, 3.0), (4.0, 4.0)])
multiline = GI.MultiLineString([line1, line2])
ring = GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])
points = [GI.Point(0.0, 0.0), GI.Point(1.0, 2.0), GI.Point(3.0, 1.0)]
multipoint = GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)])

# The extension's methods name the lookup type, which every geometry package's
# `AbstractArray{<:ItsGeometry}` method loses to, so one conversion covers every
# element type: a concrete GeoInterface wrapper as well as `Any`.
anyvec(v) = Any[v...]

@testset "the lookup's own methods win" begin
    for gl in (GeometryLookup(msquares), GeometryLookup(anyvec(msquares)), GeometryLookup(points))
        @test Base.which(Makie.convert_arguments, Tuple{Type{Makie.Poly},typeof(gl)}).module === EXT
        @test Base.which(Makie.convert_arguments, Tuple{Type{Makie.Lines},typeof(gl)}).module === EXT
        @test Base.which(Makie.convert_arguments, Tuple{Makie.PointBased,typeof(gl)}).module === EXT
        @test Base.which(Makie.plottype, Tuple{typeof(gl)}).module === EXT
    end
    dim = Geometry(GeometryLookup(msquares))
    @test Base.which(Makie.convert_arguments, Tuple{Type{Makie.Poly},typeof(dim)}).module === EXT
    @test Base.which(Makie.plottype, Tuple{typeof(dim)}).module === EXT
end

@testset "convert_arguments" begin
    @testset "polygons" begin
        for gl in (GeometryLookup(msquares), GeometryLookup(anyvec(msquares)))
            (polys,) = Makie.convert_arguments(Makie.Poly, gl)
            @test polys isa AbstractVector{<:GB.Polygon}
            @test length(polys) == 3
            @test polys[1] == GI.convert(GB, msq1)
        end
    end

    @testset "mixed polygons and multipolygons lift to multipolygons" begin
        for gl in (GeometryLookup([msq1, multi, msq2]), GeometryLookup(Any[msq1, multi, msq2]))
            (polys,) = Makie.convert_arguments(Makie.Poly, gl)
            @test polys isa AbstractVector{<:GB.MultiPolygon}
            @test length(polys) == 3
            @test length(polys[1].polygons) == 1
            @test length(polys[2].polygons) == 2
        end
    end

    @testset "linestrings" begin
        for gl in (GeometryLookup([line1, line2]), GeometryLookup(Any[line1, line2]))
            (pts,) = Makie.convert_arguments(Makie.Lines, gl)
            @test pts isa AbstractVector{<:GB.Point{2}}
            # two linestrings joined by one NaN separator
            @test length(pts) == 6
            @test count(p -> isnan(p[1]), pts) == 1
        end
        for gl in (GeometryLookup([line1, multiline]), GeometryLookup(Any[line1, multiline]))
            (pts,) = Makie.convert_arguments(Makie.Lines, gl)
            @test pts isa AbstractVector{<:GB.Point{2}}
            @test count(p -> isnan(p[1]), pts) == 2
        end
    end

    @testset "points" begin
        for gl in (GeometryLookup(points), GeometryLookup(anyvec(points)))
            (pts,) = Makie.convert_arguments(Makie.Scatter, gl)
            @test pts isa AbstractVector{<:GB.Point{2}}
            @test length(pts) == 3
            @test pts[2] == GB.Point2(1.0, 2.0)
        end
    end

    @testset "linear rings become polygons" begin
        for gl in (GeometryLookup([ring, ring]), GeometryLookup(Any[ring, ring]))
            (polys,) = Makie.convert_arguments(Makie.Poly, gl)
            @test polys isa AbstractVector{<:GB.Polygon}
            @test length(polys) == 2
            @test GI.npoint(polys[1]) == GI.npoint(ring)
        end
    end

    @testset "multipoints" begin
        gl = GeometryLookup([multipoint, multipoint])
        (pts,) = Makie.convert_arguments(Makie.Scatter, gl)
        @test pts isa AbstractVector{<:GB.Point{2}}
        @test length(pts) == 4
    end

    @testset "through a Geometry dimension" begin
        gl = GeometryLookup(msquares)
        dim = Geometry(gl)
        @test Makie.convert_arguments(Makie.Poly, dim) == Makie.convert_arguments(Makie.Poly, gl)
        # `isequal`, since the separators between geometries are `NaN`s
        @test isequal(Makie.convert_arguments(Makie.Lines, dim), Makie.convert_arguments(Makie.Lines, gl))
        @test Makie.plottype(dim) == Makie.Poly

        # any dimension wrapping the lookup, whatever its name
        origin = DD.Dim{:Origin}(gl)
        @test Makie.convert_arguments(Makie.Poly, origin) == Makie.convert_arguments(Makie.Poly, gl)
        @test Makie.plottype(origin) == Makie.Poly
    end

    @testset "unsupported mixtures" begin
        @test_throws ArgumentError Makie.convert_arguments(Makie.Poly, GeometryLookup(Any[msq1, line1]))
        @test_throws ArgumentError Makie.convert_arguments(Makie.Poly, Geometry(GeometryLookup(Any[msq1, points[1]])))
    end

    @testset "empty lookups" begin
        for gl in (GeometryLookup(GI.Polygon[]), GeometryLookup(Any[]))
            @test_throws ArgumentError Makie.convert_arguments(Makie.Poly, gl)
            @test_throws ArgumentError Makie.plottype(gl)
        end
        @test_throws ArgumentError Makie.plottype(Geometry(GeometryLookup(GI.Polygon[])))
    end
end

@testset "plottype" begin
    @test Makie.plottype(GeometryLookup(msquares)) == Makie.Poly
    @test Makie.plottype(GeometryLookup([multi])) == Makie.Poly
    @test Makie.plottype(GeometryLookup([ring])) == Makie.Poly
    @test Makie.plottype(GeometryLookup([line1, line2])) == Makie.Lines
    @test Makie.plottype(GeometryLookup([multiline])) == Makie.Lines
    @test Makie.plottype(GeometryLookup(points)) == Makie.Scatter
    @test Makie.plottype(GeometryLookup([multipoint])) == Makie.Scatter
    @test Makie.plottype(GeometryLookup(anyvec(msquares))) == Makie.Poly
    @test Makie.plottype(GeometryLookup(Any[line1, multiline])) == Makie.Lines
    @test Makie.plottype(GeometryLookup(anyvec(points))) == Makie.Scatter
    @test Makie.plottype(Geometry(GeometryLookup(anyvec(points)))) == Makie.Scatter
end

@testset "plots construct without a backend" begin
    gl = GeometryLookup(msquares)
    fig, ax, plt = poly(gl; color = [1.0, 2.0, 3.0])
    @test plt isa Makie.Poly
    @test length(plt[1][]) == 3
    @test plt.color[] == [1.0, 2.0, 3.0]

    fig, ax, plt = poly(Geometry(gl); color = [1.0, 2.0, 3.0])
    @test plt isa Makie.Poly
    @test length(plt[1][]) == 3

    fig, ax, plt = plot(gl)
    @test plt isa Makie.Poly

    fig, ax, plt = plot(GeometryLookup(anyvec(msquares)))
    @test plt isa Makie.Poly

    fig, ax, plt = plot(GeometryLookup([ring, ring]))
    @test plt isa Makie.Poly

    fig, ax, plt = scatter(GeometryLookup(points))
    @test plt isa Makie.Scatter

    fig, ax, plt = lines(GeometryLookup([line1, line2]))
    @test plt isa Makie.Lines

    fig, ax, plt = plot(Geometry(GeometryLookup(points)))
    @test plt isa Makie.Scatter
end

@testset "a choropleth from a cube" begin
    cube = rand(Geometry(GeometryLookup(msquares)), DD.Ti(1:2))
    fig, ax, plt = poly(DD.dims(cube, Geometry); color = cube[DD.Ti(1)])
    @test plt isa Makie.Poly
    @test length(plt[1][]) == 3
    @test plt.color[] == cube[DD.Ti(1)]
end
