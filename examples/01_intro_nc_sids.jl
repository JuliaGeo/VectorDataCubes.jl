#=
# Building a vector data cube: North Carolina SIDS

A port of the classic North Carolina sudden-infant-death-syndrome (SIDS) example that
introduces vector data cubes in both R's `stars` (https://r-spatial.org/r/2022/09/12/vdc.html)
and Python's `xvec` (https://xvec.readthedocs.io/en/stable/intro.html).

The dataset is 100 county polygons with birth and SIDS counts for two periods (1974 and
1979). We arrange them as a *vector data cube*: a `DimStack` whose first dimension is indexed
by the county geometries (via [`GeometryLookup`](@ref)) and whose second is the year.

The data ships with R's `sf` package; we download the shapefile straight from
its GitHub repository (~150 kB).
=#

using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import GeometryOps as GO, GeoInterface as GI
import GeoDataFrames
using DataFrames
using Downloads: download

datadir(args...) = joinpath(@__DIR__, "data", args...)
mkpath(datadir())

for ext in ("shp", "dbf", "shx", "prj")
    file = datadir("nc.$ext")
    isfile(file) || download(
        "https://raw.githubusercontent.com/r-spatial/sf/main/inst/shape/nc.$ext", file
    )
end

counties = GeoDataFrames.read(datadir("nc.shp"))

#=
## Constructing the cube

The geometry dimension is a `Geometry` dimension wrapping a `GeometryLookup` of the county
polygons. The lookup builds a spatial tree over the geometries, so spatial selectors on the
cube are fast, and it carries the CRS (the data is in NAD27, EPSG:4267).
=#

gl = GeometryLookup(counties)
# You could also write `gl = GeometryLookup(GO.get_geometries(counties); crs=GI.crs(counties))` if you want to do this manually.

# Let's now create a time dimension, called `Year`.
years = Dim{:Year}([1974, 1979])
# Finally, we can construct a DimStack with the data, which is a set of dimensional arrays
# which share some common axes.
cube = DimStack((;
    births = DimArray([counties.BIR74 counties.BIR79], (Geometry(gl), years)),
    sids = DimArray([counties.SID74 counties.SID79], (Geometry(gl), years)),
    nonwhite_births = DimArray([counties.NWBIR74 counties.NWBIR79], (Geometry(gl), years)),
))

#=
## Spatial selectors

The cube is indexed by real geometries, so we can ask spatial questions
directly. Which county contains Raleigh?
=#

raleigh = (-78.6382, 35.7796)
wake = cube[Geometry = Contains(raleigh)]
# This returns a DimStack that has certain elements.  We can select a layer of the stack:
wake.births
# and check it has the correct things:
size(wake[:births]) == (1, 2)

# The selector machinery also works on the lookup itself, which lets us recover
# the county's attributes from the original table:
wake_idx = only(Lookups.selectindices(gl, Contains(raleigh)))
println("Raleigh is in $(counties.NAME[wake_idx]) county")

#=
## Derived layers

Cube arithmetic works as usual — dimensions (including the geometry lookup) are
carried through broadcasting. Here is the SIDS rate per 1000 births:
=#

rate = @d cube.sids ./ cube.births .* 1000 name=:rate

worst_rate, worst_idx = findmax(rate[Year=At(1974)])

#=
## To a table

[`vectordatacubetable`](@ref) flattens the cube to a table: one row per county × year, the
county polygons in a `Geometry` column, and the crs as table metadata. Any Tables.jl
constructor accepts it, and `DataFrame` keeps the metadata with the data:
=#

tbl = DataFrame(vectordatacubetable(rate))
GI.crs(tbl)
# We can sort the table by the value, and get the top 5 counties by rate:
first(sort(tbl, :rate; rev=true), 5)
# and write it to a shapefile — `GeoDataFrames.write` finds the geometry column
# and the crs in the table's metadata, so the `.prj` records NAD27:
GeoDataFrames.write(datadir("nc_rates.shp"), tbl)
GI.crs(GeoDataFrames.read(datadir("nc_rates.shp")))