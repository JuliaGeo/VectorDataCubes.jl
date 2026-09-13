#=
# Two geometry dimensions: NYC taxi origin × destination cube

A port of xvec's "Indexing" tutorial (https://xvec.readthedocs.io/en/stable/indexing.html):
arrange NYC yellow-taxi trip counts as a cube whose *two* dimensions are both indexed by
taxi-zone polygons (origins and destinations), then slice it spatially on either axis.

Nothing in `GeometryLookup` is tied to the `Geometry` dimension: any dimension
can carry one, so a cube can have several geometry-indexed axes at once.

Data (both from the NYC Taxi & Limousine Commission's open data, no auth):
- yellow-taxi trip records, January 2022 (~37 MB Parquet, ~2.5M trips),
- the taxi-zone polygons (~1 MB zipped shapefile, EPSG:2263).
=#

using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import GeometryOps as GO, GeoInterface as GI
import GeoDataFrames, Parquet2, ZipFile
import Tables
using Downloads: download

datadir(args...) = joinpath(@__DIR__, "data", args...)
mkpath(datadir())

zonefile = datadir("taxi_zones.shp")
if !isfile(zonefile)
    zippath = download("https://d37ci6vzurychx.cloudfront.net/misc/taxi_zones.zip")
    archive = ZipFile.Reader(zippath)
    for f in archive.files
        endswith(f.name, "/") && continue # skip directory entries
        write(datadir(basename(f.name)), read(f))
    end
    close(archive)
end

tripfile = datadir("yellow_tripdata_2022-01.parquet")
isfile(tripfile) || download(
    "https://d37ci6vzurychx.cloudfront.net/trip-data/yellow_tripdata_2022-01.parquet",
    tripfile,
)

#=
## The geometry dimensions

The taxi zones ship in EPSG:2263 (NY Long Island, US feet); we reproject them to lon/lat once
so the spatial selectors below can use familiar coordinates. Origins and destinations share
the same zones, so both lookups are built from the same geometry vector.
=#

zones = GeoDataFrames.read(zonefile)
geoms = GO.reproject(
    GO.get_geometries(zones);
    source_crs=EPSG(2263), target_crs=EPSG(4326), always_xy=true,
)

Origin = Dim{:Origin}(GeometryLookup(geoms; crs=EPSG(4326)))
Destination = Dim{:Destination}(GeometryLookup(geoms; crs=EPSG(4326)))

#=
## The cube

Count trips per (pickup zone, dropoff zone) pair. Trip records reference
zones by `LocationID`; IDs not in the shapefile (264/265, "unknown") are
dropped.
=#

trips = Parquet2.Dataset(tripfile)
pickup_ids = Tables.getcolumn(trips, :PULocationID)
dropoff_ids = Tables.getcolumn(trips, :DOLocationID)

rowof = Dict(id => i for (i, id) in enumerate(zones.LocationID))
counts = zeros(Int, length(geoms), length(geoms))
for (pu, dropoff) in zip(pickup_ids, dropoff_ids)
    i = get(rowof, pu, 0)
    j = get(rowof, dropoff, 0)
    (i == 0 || j == 0) && continue
    counts[i, j] += 1
end

od = Raster(counts, (Origin, Destination); name=:trips)
println(sum(od), " of ", length(pickup_ids), " trips between known zones")

#=
## Spatial selectors on either axis

`Contains` finds the zone containing a point — on whichever geometry
dimension it is applied to. Where do trips from the zone around Times Square
go?
=#

times_square = (-73.9857, 40.7580)
jfk = (-73.7781, 40.6413)

from_tsq = vec(parent(od[Origin=Contains(times_square)]))
top5 = sortperm(from_tsq; rev=true)[1:5]
println("top destinations from Times Square:")
for j in top5
    println("  ", rpad(zones.zone[j], 28), from_tsq[j], " trips")
end

# ...and on the destination axis: where do trips *into* JFK come from?
# (The busiest "origin" is JFK itself — intra-zone trips — so show the top 3.)
into_jfk = vec(parent(od[Destination=Contains(jfk)]))
println("top origins for JFK dropoffs:")
for i in sortperm(into_jfk; rev=true)[1:3]
    println("  ", rpad(zones.zone[i], 28), into_jfk[i], " trips")
end

# Both axes at once — how many trips ran Times Square -> JFK that month?
tsq_to_jfk = only(od[Origin=Contains(times_square), Destination=Contains(jfk)])
println("Times Square -> JFK: ", tsq_to_jfk, " trips")
tsq_to_jfk == from_tsq[only(Lookups.selectindices(val(Destination), Contains(jfk)))]

#=
## To a table

[`vectordatacubetable`](@ref) flattens the origin-destination cube to one row per zone pair
with real polygons in both geometry columns: every `GeometryLookup`-backed dimension becomes
one, recorded with the shared crs as table metadata that `GeoInterface` and `DataFrame` read.
=#

tbl = vectordatacubetable(od)

# The columns are the two geometry dimensions plus the value layer, and both
# dimensions are geometry columns:
GI.geometrycolumns(tbl), GI.crs(tbl)
