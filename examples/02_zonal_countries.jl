#=
# Zonal statistics over countries: NCEP air temperature

A port of xvec's "Zonal statistics" tutorial (https://xvec.readthedocs.io/en/stable/zonal_stats.html)
and the R `stars` `aggregate(raster, by = polygons)` example (https://r-spatial.github.io/stars/):
aggregate a gridded temperature time series over country polygons into a `(Ti, Geometry)` cube.

Data:
- `air_temperature.nc`, xarray's tutorial dataset: NCEP reanalysis 2 m air temperature over
  North America, 6-hourly 2013–2014, ~7 MB NetCDF, from the `pydata/xarray-data` repository.
- Country polygons from Natural Earth (1:110m admin-0), via NaturalEarth.jl.
=#

using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
## Rasters exports a `zonal` too, so bind VectorDataCubes' explicitly: it
## returns vector data cubes for geometry-lookup `of`s and forwards any other
## `of` to `Rasters.zonal`.
using VectorDataCubes: zonal
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import NCDatasets # activates Rasters' NetCDF backend
import Proj # activates reprojection of geometry lookups
using NaturalEarth
using DataFrames
using Statistics: mean
using Dates
using Downloads: download

datadir(args...) = joinpath(@__DIR__, "data", args...)
mkpath(datadir())
airfile = datadir("air_temperature.nc")
isfile(airfile) || download(
    "https://raw.githubusercontent.com/pydata/xarray-data/master/air_temperature.nc",
    airfile,
)

#=
## The raster cube

The NetCDF stores temperature in Kelvin on 0–360° longitudes; we shift the longitudes to the
-180–180° convention Natural Earth uses, convert to °C, and reduce the 6-hourly series to a
per-cell monthly climatology with `groupby`/`combine`: an `(X, Y, Ti)` cube of 12 time steps.
=#

air = RasterStack(airfile)[:air]
## Shift with a range (not a plain vector) so the lookup keeps its regular
## step — Rasters' masking machinery requires a regular grid.
xs = DD.lookup(air, X)
air = set(air, X => range(first(xs) - 360, last(xs) - 360; length=length(xs)))
air = air .- 273.15

monthly = DD.combine(mean, groupby(air, Ti => month); dims=Ti)

#=
## The geometry dimension: countries as a vector data cube

Natural Earth returns a GeoJSON FeatureCollection. `vectordatacube` lifts it into a `DimStack`
over a `Geometry` dimension: the polygons become a `GeometryLookup` and the attribute columns
become layers over it, so subsetting the cube keeps attributes and geometries aligned.

We keep three attribute layers and select the North American countries.
=#

countries = vectordatacube(
    naturalearth("admin_0_countries", 110);
    layers=(:NAME, :CONTINENT, :POP_EST), crs=EPSG(4326),
)
northam = countries[Geometry=findall(==("North America"), countries[:CONTINENT])]
geodim = DD.dims(northam, Geometry)

#=
## Zonal aggregation -> vector data cube

Passing the geometry dimension (or its lookup) as `of` makes `zonal` return a cube over
`(Ti, Geometry)`: the mean is computed per country *per month* (each spatial slice
separately), and the result keeps the geometry lookup, so spatial selectors work on it.

Two coverage caveats:
- values are means over the raster's coverage (lon 160°W–30°W, lat 15°N–75°N); a country
  partly outside that window is averaged over the overlapping part only;
- at 2.5° resolution an island nation smaller than a grid cell can cover *no* cell center;
  `emptyval = missing` makes those `missing` rather than `mean` of an empty slice (`NaN`).
=#

temps = zonal(mean, monthly; of=geodim, emptyval=missing, progress=false)
size(temps) == (12, length(geodim))

# Seasonal cycle of the country containing a point in the US Great Plains:
usa = temps[Geometry(Contains((-100.0, 40.0)))]
round.(vec(parent(usa)); digits=1)

#=
Countries entirely outside the raster's coverage (fully south of 15°N, like Panama and Costa
Rica) come back as `missing`, so reductions over the geometry axis should `skipmissing`.

Which country is warmest in July, and coldest in January? `temps` and `northam[:NAME]` share
the geometry dimension, so a position found in one indexes the other.
=#
july = temps[Ti=At(7)]
january = temps[Ti=At(1)]
println("warmest in July:    ", northam[:NAME][argmax(skipmissing(parent(july)))])
println("coldest in January: ", northam[:NAME][argmin(skipmissing(parent(january)))])

#=
## Reprojection

The geometry lookup carries its CRS, so the whole cube can be reprojected;
only the geometry dimension is affected.
=#

temps_3857 = reproject(EPSG(3857), temps)
crs(DD.lookup(temps_3857, Geometry)) == EPSG(3857)

#=
## To a table

One row per (month, country); the `:Geometry` column holds the country
polygons.
=#

df = DataFrame(vectordatacubetable(temps))
july_df = dropmissing(subset(df, :Ti => ByRow(==(7))))
first(sort(july_df, :value; rev=true), 5)
