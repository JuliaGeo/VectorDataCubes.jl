using Test
using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import Tables, DataAPI

# A small axis-aligned unit square at (x, y), like the idiom in test/selectors.jl.
# Coordinates are always Float64: mixed Int/Float64 tuples break GeometryOps predicates.
function square(x, y; s=1.0)
    x, y = float(x), float(y)
    GI.Polygon([GI.LinearRing([(x, y), (x + s, y), (x + s, y + s), (x, y + s), (x, y)])])
end

@testset "Tables.jl integration" begin
    geoms = [square(0, 0), square(2, 0), square(0, 2), square(2, 2)]
    gl = GeometryLookup(geoms; crs=EPSG(4326))

    @testset "1-D cube (Geometry only)" begin
        A = rand(Geometry(gl))
        for tbl in (DD.DimTable(A), vectordatacubetable(A))
            ct = Tables.columntable(tbl)
            @test :Geometry in keys(ct)
            @test :value in keys(ct)
            @test length(ct.Geometry) == length(geoms)
            # The Geometry column holds the *actual* geometry objects.
            @test all(GI.isgeometry, ct.Geometry)
            @test all(splat(GO.equals), zip(ct.Geometry, geoms))

            rt = Tables.rowtable(tbl)
            @test length(rt) == length(geoms)
            @test GO.equals(rt[1].Geometry, geoms[1])
        end
    end

    @testset "2-D cube (Geometry × Ti)" begin
        A = rand(Geometry(gl), Ti(1:3))
        for tbl in (DD.DimTable(A), vectordatacubetable(A))
            ct = Tables.columntable(tbl)
            @test Set(keys(ct)) == Set((:Geometry, :Ti, :value))
            @test length(ct.Geometry) == length(geoms) * 3
            @test all(GI.isgeometry, ct.Geometry)
            @test Set(ct.Ti) == Set(1:3)
            # Every (geometry, time) combination appears exactly once.
            @test length(unique(zip(map(GO.centroid, ct.Geometry), ct.Ti))) == length(geoms) * 3
        end
    end

    @testset "DimStack cube (multiple value columns)" begin
        st = DimStack((a=rand(Geometry(gl)), b=rand(Geometry(gl))))
        for tbl in (DD.DimTable(st), vectordatacubetable(st))
            ct = Tables.columntable(tbl)
            @test Set(keys(ct)) == Set((:Geometry, :a, :b))
            @test length(ct.Geometry) == length(geoms)
            @test all(GI.isgeometry, ct.Geometry)
            @test eltype(ct.a) == Float64
            @test eltype(ct.b) == Float64
        end
    end

    @testset "the table forwards the Tables.jl interface to its DimTable" begin
        A = rand(Geometry(gl), Ti(1:2); name=:v)
        tbl = vectordatacubetable(A)
        @test tbl isa VectorDataCubes.VectorDataCubeTable
        @test parent(tbl) isa DD.DimTable
        @test Tables.istable(tbl)
        @test Tables.columnaccess(typeof(tbl))
        @test Tables.columns(tbl) === tbl
        @test Tables.columnnames(tbl) == (:Geometry, :Ti, :v)
        @test Tables.schema(tbl).names == (:Geometry, :Ti, :v)
        @test Tables.getcolumn(tbl, :Ti) == [1, 1, 1, 1, 2, 2, 2, 2]
        @test Tables.getcolumn(tbl, 2) == Tables.getcolumn(tbl, Ti)
        @test Tables.getcolumn(tbl, Float64, 3, :v) == Tables.getcolumn(tbl, :v)
        @test tbl.v == vec(parent(A))
        @test (DataAPI.nrow(tbl), DataAPI.ncol(tbl)) == (8, 3)
        rows = Tables.rows(tbl)
        @test length(rows) == 8
        @test Tables.getcolumn(first(rows), :Ti) == 1
    end

    @testset "a sliced cube keeps its refdims as columns" begin
        A = rand(Geometry(gl), Ti(1:2); name=:v)
        tbl = vectordatacubetable(A[Ti=1])
        @test Tables.columnnames(tbl) == (:Geometry, :Ti, :v)
        @test GI.geometrycolumns(tbl) == (:Geometry,)
        @test GI.crs(tbl) == EPSG(4326)
        # The geometry dimension itself can be the refdim.
        onegeom = vectordatacubetable(A[Geometry=1])
        @test GI.geometrycolumns(onegeom) == (:Geometry,)
        @test length(Tables.getcolumn(onegeom, :Geometry)) == 2
    end

    @testset "crs and geometry columns are DataAPI metadata" begin
        A = rand(Geometry(gl), Ti(1:2))
        tbl = vectordatacubetable(A)
        @test GI.crs(tbl) == EPSG(4326)
        @test GI.geometrycolumns(tbl) == (:Geometry,)
        # What DataFrames copies when it builds a DataFrame from the table.
        @test DataAPI.metadatasupport(typeof(tbl)) == (read=true, write=false)
        @test Set(DataAPI.metadatakeys(tbl)) == Set(("GEOINTERFACE:geometrycolumns", "GEOINTERFACE:crs"))
        @test DataAPI.metadata(tbl, "GEOINTERFACE:crs") == EPSG(4326)
        @test DataAPI.metadata(tbl, "GEOINTERFACE:crs"; style=true) == (EPSG(4326), :note)
        @test DataAPI.metadata(tbl, "GEOINTERFACE:geometrycolumns"; style=true) == ((:Geometry,), :note)
        @test DataAPI.metadata(tbl, "nope", 1) == 1
        @test DataAPI.metadata(tbl, "nope", 1; style=true) == (1, :default)
        @test_throws ArgumentError DataAPI.metadata(tbl, "nope")
        @test DataAPI.metadata(tbl) == Dict(
            "GEOINTERFACE:geometrycolumns" => (:Geometry,), "GEOINTERFACE:crs" => EPSG(4326),
        )
    end

    @testset "show names the geometry columns and the crs" begin
        A = rand(Geometry(gl), Ti(1:2); name=:v)
        str = sprint(show, vectordatacubetable(A))
        @test startswith(str, "VectorDataCubeTable with 8 rows, 3 columns, geometry column :Geometry, crs ")
        @test occursin(":v", str)
        # No crs, no crs clause.
        nocrs = sprint(show, vectordatacubetable(rand(Geometry(GeometryLookup(geoms)))))
        @test startswith(nocrs, "VectorDataCubeTable with 4 rows, 2 columns, geometry column :Geometry, and schema:")
        # A long crs (a WKT string, a proj string) is truncated, not splashed over the header.
        proj = ProjString("+proj=longlat +datum=WGS84 +no_defs +note=" * repeat("x", 200))
        long = sprint(show, vectordatacubetable(rand(Geometry(GeometryLookup(geoms; crs=proj)))))
        @test length(first(split(long, '\n'))) < 160
        @test occursin("…", long)
    end

    @testset "no-crs cube has no crs metadata" begin
        gl2 = GeometryLookup(geoms)  # hand-made polygons have no crs
        A = rand(Geometry(gl2))
        tbl = vectordatacubetable(A)
        @test GI.crs(tbl) === nothing
        @test GI.geometrycolumns(tbl) == (:Geometry,)
        @test collect(DataAPI.metadatakeys(tbl)) == ["GEOINTERFACE:geometrycolumns"]
        @test all(GI.isgeometry, Tables.columntable(tbl).Geometry)
    end

    @testset "several geometry dimensions" begin
        Origin = Dim{:Origin}(GeometryLookup(geoms; crs=EPSG(4326)))
        Destination = Dim{:Destination}(GeometryLookup(geoms; crs=EPSG(4326)))
        od = Raster(reshape(1:16, 4, 4), (Origin, Destination); name=:trips)
        tbl = vectordatacubetable(od)
        @test Tables.columnnames(tbl) == (:Origin, :Destination, :trips)
        @test GI.geometrycolumns(tbl) == (:Origin, :Destination)
        @test GI.crs(tbl) == EPSG(4326)
        ct = Tables.columntable(tbl)
        @test all(GI.isgeometry, ct.Origin)
        @test all(GI.isgeometry, ct.Destination)
        @test ct.trips == 1:16
        @test occursin("geometry columns :Origin, :Destination", sprint(show, tbl))
        # Which of the two indexes a cube is the user's call, not ours.
        err = try vectordatacube(tbl) catch e e end
        @test err isa ArgumentError
        @test occursin(":Origin, :Destination", err.msg)
        @test keys(vectordatacube(tbl; geometrycolumn=:Origin)) == (:Destination, :trips)

        # A lookup without a crs does not count as disagreement.
        NoCRS = Dim{:Destination}(GeometryLookup(geoms))
        @test GI.crs(vectordatacubetable(Raster(zeros(4, 4), (Origin, NoCRS)))) == EPSG(4326)
        # Two different crs do.
        Other = Dim{:Destination}(GeometryLookup(geoms; crs=EPSG(3857)))
        @test_throws ArgumentError vectordatacubetable(Raster(zeros(4, 4), (Origin, Other)))
    end

    @testset "errors on non-vector cube" begin
        ras = rand(X(1:3), Y(1:3))
        @test_throws ArgumentError vectordatacubetable(ras)
    end
end

@testset "vectordatacube: table -> cube" begin
    geoms = [square(0, 0), square(2, 0), square(0, 2), square(2, 2)]
    tbl = (geometry=geoms, name=["a", "b", "c", "d"], pop=[10, 20, 30, 40])

    @testset "attribute columns become layers over Geometry" begin
        cube = vectordatacube(tbl; crs=EPSG(4326))
        @test cube isa DD.DimStack
        @test keys(cube) == (:name, :pop)
        @test DD.lookup(cube, Geometry) isa GeometryLookup
        @test crs(DD.lookup(cube, Geometry)) == EPSG(4326)
        @test all(splat(GO.equals), zip(val(DD.lookup(cube, Geometry)), geoms))
        @test parent(cube[:name]) == tbl.name
        @test parent(cube[:pop]) == tbl.pop
    end

    @testset "attributes stay aligned under subsetting" begin
        cube = vectordatacube(tbl)
        sub = cube[Geometry=2:3]
        @test parent(sub[:name]) == ["b", "c"]
        @test length(DD.lookup(sub, Geometry)) == 2
        # spatial selectors too
        hit = cube[Geometry(Contains((2.5, 0.5)))]
        @test parent(hit[:name]) == ["b"]
    end

    @testset "layers keyword" begin
        @test keys(vectordatacube(tbl; layers=(:pop,))) == (:pop,)
        @test keys(vectordatacube(tbl; layers=:pop)) == (:pop,)
        @test keys(vectordatacube(tbl; layers="pop")) == (:pop,)
        @test keys(vectordatacube(tbl; layers=[:pop, :name])) == (:pop, :name)
        @test keys(vectordatacube(tbl; layers=["name"])) == (:name,)
        @test_throws ArgumentError vectordatacube(tbl; layers=:nope)
        @test_throws ArgumentError vectordatacube(tbl; layers=(:pop, :nope))
        @test_throws ArgumentError vectordatacube(tbl; layers=(:pop, :pop))
        @test_throws ArgumentError vectordatacube(tbl; layers=())
        @test keys(vectordatacube(tbl; layers=(n for n in (:name, :pop)))) == (:name, :pop)
        @test occursin(":nope", sprint(showerror, try vectordatacube(tbl; layers=:nope) catch e e end))
        # A wide table lists a few of its columns, not all of them.
        wide = merge((geometry=geoms,), NamedTuple{ntuple(i -> Symbol(:c, i), 12)}(ntuple(i -> 1:4, 12)))
        msg = sprint(showerror, try vectordatacube(wide; layers=:nope) catch e e end)
        @test occursin("c1, c2", msg) && !occursin("c12", msg)
    end

    @testset "geometrycolumn keyword" begin
        renamed = (geom=geoms, name=tbl.name)
        @test_throws ArgumentError vectordatacube(renamed)
        cube2 = vectordatacube(renamed; geometrycolumn=:geom)
        @test parent(cube2[:name]) == tbl.name
        # A String names a column just as well as a Symbol.
        @test parent(vectordatacube(renamed; geometrycolumn="geom")[:name]) == tbl.name
    end

    @testset "missing geometries are an error naming the rows" begin
        holes = (geometry=[geoms[1], missing, geoms[3], missing], name=tbl.name)
        err = try vectordatacube(holes) catch e e end
        @test err isa ArgumentError
        @test occursin("rows 2, 4", err.msg)
        # A `Union{Missing, T}` column without actual missings is fine, and the lookup
        # that comes out of it holds a concrete geometry type.
        nomissing = (geometry=Union{Missing,eltype(geoms)}[geoms...], name=tbl.name)
        cube = vectordatacube(nomissing)
        @test parent(cube[:name]) == tbl.name
        @test eltype(DD.lookup(cube, Geometry)) == eltype(geoms)
    end

    @testset "round trip through vectordatacubetable" begin
        cube = vectordatacube(tbl; crs=EPSG(4326))
        t = vectordatacubetable(cube)
        ct = Tables.columntable(t)
        @test Set(keys(ct)) == Set((:Geometry, :name, :pop))
        @test all(splat(GO.equals), zip(ct.Geometry, geoms))
        @test ct.name == tbl.name
        @test ct.pop == tbl.pop
        # Back to a cube: the geometry column and the crs come from the table's metadata.
        back = vectordatacube(t)
        @test keys(back) == (:name, :pop)
        @test crs(DD.lookup(back, Geometry)) == EPSG(4326)
        @test all(splat(GO.equals), zip(val(DD.lookup(back, Geometry)), geoms))
        @test parent(back[:name]) == tbl.name
        @test parent(back[:pop]) == tbl.pop

        # Without a crs the round trip stays without one.
        plain = vectordatacube(vectordatacubetable(vectordatacube(tbl)))
        @test crs(DD.lookup(plain, Geometry)) === nothing
        @test keys(plain) == (:name, :pop)
    end

    @testset "a table with no rows makes an empty cube" begin
        empty = (geometry=eltype(geoms)[], name=String[])
        cube = vectordatacube(empty)
        @test keys(cube) == (:name,)
        @test isempty(DD.lookup(cube, Geometry))
    end

    @testset "errors" begin
        @test_throws ArgumentError vectordatacube(geoms)             # not a table
        @test_throws ArgumentError vectordatacube((geometry=geoms,)) # no attribute columns
    end
end
