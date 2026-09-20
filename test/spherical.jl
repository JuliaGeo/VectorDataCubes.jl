using Test
using VectorDataCubes
using VectorDataCubes: spatialtree, zonal

using Rasters, DimensionalData
using Rasters.Lookups
import DataAPI
import DE9IM
import DimensionalData as DD
import Extents
import GeoInterface as GI
import GeometryOps as GO

ring(points) = GI.Polygon([GI.LinearRing([points..., first(points)])])

# With unoriented spherical edges this is the small rectangle straddling the dateline,
# rather than the large planar rectangle between -170 and 170 degrees.
dateline = ring([(170.0, -10.0), (-170.0, -10.0), (-170.0, 10.0), (170.0, 10.0)])
inside = (179.0, 0.0)
overlap = ring([(175.0, -5.0), (-175.0, -5.0), (-175.0, 15.0), (175.0, 15.0)])
crossing = GI.LineString([(160.0, 0.0), (-160.0, 0.0)])
far = ring([(20.0, -5.0), (30.0, -5.0), (30.0, 5.0), (20.0, 5.0)])

@testset "spherical lookup lifecycle" begin
    planar = GeometryLookup([dateline])
    sphere = GO.Spherical()
    spherical = GeometryLookup([dateline]; manifold=sphere)

    @test planar.manifold == GO.Planar()
    @test spherical.manifold === sphere
    @test DD.rebuild(spherical; data=parent(spherical)).manifold === sphere
    @test view(spherical, 1:1).manifold === sphere
    @test reverse(spherical).manifold === sphere
    @test spherical != planar
    @test !isequal(spherical, planar)
    @test hash(spherical) != hash(planar)
    usp = GO.UnitSpherical.UnitSphericalPoint((179.0, 0.0))
    @test_throws ArgumentError GeometryLookup([usp]; manifold=sphere)
    @test_throws ArgumentError DD.rebuild(spherical; data=[usp])

    # The public bounds remain longitude/latitude bounds even though the spherical tree is XYZ.
    bounds_before = Lookups.bounds(spherical)
    @test bounds_before == Lookups.bounds(planar)
    @test hasproperty(Extents.extent(spatialtree(spherical)), :Z)
    @test !hasproperty(Extents.extent(spatialtree(planar)), :Z)
    @test Lookups.bounds(spherical) == bounds_before

    # A cached index is tied to the manifold, including spherical orientation.
    planar_tree = spatialtree(planar)
    rebuilt_sphere = DD.rebuild(planar; manifold=GO.Spherical())
    @test isnothing(VectorDataCubes._builttree(rebuilt_sphere))
    @test spatialtree(rebuilt_sphere) !== planar_tree
    rebuilt_oriented = DD.rebuild(rebuilt_sphere; manifold=GO.Spherical(oriented=true))
    @test isnothing(VectorDataCubes._builttree(rebuilt_oriented))
    @test spatialtree(rebuilt_oriented) !== spatialtree(rebuilt_sphere)
end

@testset "spherical selectors" begin
    geoms = [dateline, overlap, far]
    sphere = GO.Spherical()
    indexed = GeometryLookup(geoms; manifold=sphere)
    scanned = GeometryLookup(geoms; manifold=sphere, tree=nothing)
    inds(gl, sel) = sort(Lookups.selectindices(gl, sel))

    @test inds(indexed, Contains(inside)) == [1, 2]
    @test inds(indexed, Contains(inside)) == inds(scanned, Contains(inside))
    @test Lookups.selectindices(indexed, At(dateline)) == 1

    predicate_cases = (
        Where(GO.intersects(overlap)),
        Where(GO.crosses(crossing)),
        Where(GO.overlaps(overlap)),
        DE9IM.Intersects(overlap),
        DE9IM.Covers(inside),
    )
    for selector in predicate_cases
        @test inds(indexed, selector) == inds(scanned, selector)
    end
    @test inds(indexed, Where(GO.crosses(crossing))) == [1, 2]
    @test inds(indexed, Where(GO.overlaps(overlap))) == [1]

    # Finite coordinate boxes become polygons whose four edges are geodesics.
    interval = (X(175.0 .. 179.0), Y(-5.0 .. 5.0))
    extent = Extents.Extent(X=(175.0, 179.0), Y=(-5.0, 5.0))
    @test inds(indexed, interval) == inds(scanned, interval)
    @test inds(indexed, Touches(extent)) == inds(scanned, Touches(extent))

    for invalid in (
        (X(-Inf .. 10.0), Y(-5.0 .. 5.0)),
        (X(10.0 .. 190.0), Y(-5.0 .. 5.0)),
        (X(10.0 .. 5.0), Y(-5.0 .. 5.0)),
        (X(10.0 .. 20.0), Y(-90.0 .. 5.0)),
    )
        @test_throws ArgumentError Lookups.selectindices(indexed, invalid)
    end
    @test_throws ArgumentError Lookups.selectindices(
        indexed, Touches(Extents.Extent(X=(-Inf, 10.0), Y=(-5.0, 5.0)))
    )

    # The northern edge is a great-circle arc, which bows north of 60 degrees.
    curved_box_points = GeometryLookup(
        [GI.Point(0.0, 65.0), GI.Point(0.0, 70.0)]; manifold=sphere
    )
    curved_box = (X(-45.0 .. 45.0), Y(0.0 .. 60.0))
    @test Lookups.selectindices(curved_box_points, curved_box) == [1]

    polar_cap = ring([(-135.0, 80.0), (-45.0, 80.0), (45.0, 80.0), (135.0, 80.0)])
    polar_lookup = GeometryLookup([polar_cap]; manifold=sphere)
    @test Lookups.selectindices(polar_lookup, Contains((0.0, 90.0))) == [1]
end

@testset "Float32 spherical index extents" begin
    sphere = GO.Spherical()
    polygon32 = ring([(170.0f0, -10.0f0), (-170.0f0, -10.0f0),
                      (-170.0f0, 10.0f0), (170.0f0, 10.0f0)])
    indexed32 = GeometryLookup([polygon32]; manifold=sphere)
    scanned32 = GeometryLookup([polygon32]; manifold=sphere, tree=nothing)
    point32 = GI.Point(179.0f0, 0.0f0)

    @test Lookups.selectindices(indexed32, Contains(point32)) == [1]
    @test Lookups.selectindices(indexed32, Contains(point32)) ==
        Lookups.selectindices(scanned32, Contains(point32))
    @test hasproperty(Extents.extent(spatialtree(indexed32)), :Z)
end

@testset "spherical Near is exact for points" begin
    sphere = GO.Spherical()
    points = [GI.Point(-179.0, 0.0), GI.Point(150.0, 0.0), GI.Point(-179.0, 0.0)]
    for gl in (GeometryLookup(points; manifold=sphere),
               GeometryLookup(points; manifold=sphere, tree=nothing))
        @test Lookups.selectindices(gl, Near((179.0, 0.0))) == 1
        @test Lookups.selectindices(gl, Near((-179.0, 0.0))) == 1
    end
    @test_throws ArgumentError Lookups.selectindices(
        GeometryLookup([crossing]; manifold=sphere), Near((179.0, 0.0))
    )
    @test_throws ArgumentError Lookups.selectindices(
        GeometryLookup([dateline]; manifold=sphere), Near((179.0, 0.0))
    )
    @test !Lookups.hasselection(GeometryLookup([dateline]; manifold=sphere), Near((179.0, 0.0)))
end

@testset "operations that require another coordinate space" begin
    spherical = GeometryLookup([dateline]; manifold=GO.Spherical(), crs=EPSG(4326))
    @test_throws ArgumentError Rasters.reproject(EPSG(3857), spherical)

    raster = Raster(ones(4, 4), (X(172.5:5.0:187.5), Y(-7.5:5.0:7.5)))
    @test_throws ArgumentError zonal(sum, raster; of=spherical, progress=false)
end

@testset "table edge metadata" begin
    import Proj
    planar = GeometryLookup([far]; crs=EPSG(4326))
    spherical = GeometryLookup([dateline]; manifold=GO.Spherical(), crs=EPSG(4326))
    oriented = GeometryLookup([dateline]; manifold=GO.Spherical(oriented=true), crs=EPSG(4326))

    mixed = Raster(reshape(1:1, 1, 1),
        (Dim{:PlanarGeometry}(planar), Dim{:SphericalGeometry}(oriented)); name=:value)
    table = vectordatacubetable(mixed)

    @test GI.geometrycolumns(table) == (:PlanarGeometry, :SphericalGeometry)
    @test GI.crs(table) == EPSG(4326)
    @test DataAPI.colmetadata(table, :PlanarGeometry, "edges") == "planar"
    @test DataAPI.colmetadata(table, :SphericalGeometry, "edges") == "spherical"
    @test DataAPI.colmetadata(table, :SphericalGeometry, "orientation") == "counterclockwise"
    @test DataAPI.colmetadata(table, :PlanarGeometry, "orientation", nothing) === nothing
    @test DataAPI.colmetadata(table, "SphericalGeometry", "edges") == "spherical"
    @test DataAPI.colmetadata(table, 2, "edges") == "spherical"
    @test DataAPI.colmetadata(table, 2, "orientation"; style=true) ==
        ("counterclockwise", :note)
    @test Set(DataAPI.metadatakeys(table)) ==
        Set((GI.GEOINTERFACE_GEOMETRYCOLUMNS_KEY, GI.GEOINTERFACE_CRS_KEY))

    inherited = vectordatacube(table; geometrycolumn=:SphericalGeometry, layers=:value)
    inherited_lookup = DD.lookup(inherited, Geometry)
    @test inherited_lookup.manifold isa GO.Spherical
    @test inherited_lookup.manifold.oriented
    @test parent(inherited_lookup) == table.SphericalGeometry

    overridden = vectordatacube(table; geometrycolumn=:SphericalGeometry, layers=:value,
        manifold=GO.Planar())
    @test DD.lookup(overridden, Geometry).manifold == GO.Planar()

    roundtrip = vectordatacubetable(Raster([1], (Geometry(spherical),); name=:value))
    back = vectordatacube(roundtrip)
    @test DD.lookup(back, Geometry).manifold isa GO.Spherical
    @test parent(DD.lookup(back, Geometry)) == parent(spherical)

    custom_radius = GeometryLookup([dateline]; manifold=GO.Spherical(radius=1.0))
    @test_throws ArgumentError vectordatacubetable(
        Raster([1], (Geometry(custom_radius),); name=:value)
    )
end
