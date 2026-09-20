# VectorDataCubes.jl

Work with multidimensional data over points, lines, and polygons.

A vector data cube associates array values with geometries and other dimensions,
such as time or measurement type. VectorDataCubes connects those geometries to
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) and
[Rasters.jl](https://github.com/rafaqz/Rasters.jl) through the `Geometry` dimension.
You can select regions spatially while keeping values, attributes, and other
dimensions aligned.

The package lets you:

- Build arrays and stacks with a spatially indexed geometry axis.
- Select geometries by point, bounding box, or geometric relationship.
- Aggregate rasters over geometries, retaining time and other dimensions.
- Convert between tables with geometry columns and vector data cubes.

## Quick start

In the Julia REPL, install the package from this repository and the packages
used in the examples:

```julia
using Pkg
Pkg.add("VectorDataCubes")
Pkg.add(["DimensionalData", "GeoInterface", "Rasters", "Statistics", "Tables"]) # needed to use the package
```

This example stores two years of observations over two square regions:

```julia
using VectorDataCubes
using DimensionalData
import GeoInterface as GI

square(x, y) = GI.Polygon([GI.LinearRing([
    (x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1), (x, y),
])])

geometries = [square(0.0, 0.0), square(2.0, 0.0)]
regions = GeometryLookup(geometries)
cube = DimArray([10 12; 20 24], (Geometry(regions), Ti([2020, 2021]));
                name = :observations)

size(cube)                                      # (2, 2)
cube[Geometry(Contains((0.5, 0.5))), Ti(At(2021))] # one region, value 12
cube[Geometry = (X(-0.1 .. 1.1), Y(-0.1 .. 1.1))] # first region, both years
```

`Contains(point)` returns all geometries containing the point, so the geometry
axis remains even when only one region matches. The `X` and `Y` intervals select
geometries fully covered by the box. Keep these coordinate selectors together
inside `Geometry`. Use `X(Touches(a, b))` and `Y(Touches(c, d))` there to include
geometries that intersect the box instead.

The lookup accepts GeoInterface-compatible geometries, including those read
from GeoJSON and Shapefile tables. Supply `crs` when the input does not carry its
coordinate reference system. Spatial queries use the geometries' coordinates;
transform query geometries to the same CRS before selecting.

## Spherical geometries

Select longitude/latitude geometries with great-circle edges explicitly:

```julia
import GeometryOps as GO
regions = GeometryLookup(geometries; manifold=GO.Spherical())
attributes = vectordatacube((geometry=geometries, value=[10, 20]); manifold=GO.Spherical())
```

The lazy index uses unit-sphere XYZ bounds while stored and emitted geometries
remain longitude/latitude. Finite `X`/`Y` interval boxes become polygons with
great-circle edges. This first cut accepts longitude widths below 180 degrees
and latitude bounds strictly between the poles; pass a geometry for other regions.
Spherical `Near` supports point lookups, using the XYZ tree when enabled and an
exhaustive scan with `tree=nothing`. Spherical zonal statistics and reprojection
remain follow-ups. When Proj is loaded, spherical table input derives its radius
from the CRS datum: a declared sphere keeps its radius and an ellipsoid uses
`(2a+b)/3`. The `edges` metadata still selects spherical behavior independently.

Table conversion preserves `GEOINTERFACE:crs` and `GEOINTERFACE:geometrycolumns`,
plus per-column DataAPI `edges` and `orientation` metadata using GeoParquet
semantics. Absent `edges` means planar. A custom radius can round-trip only with a
CRS that describes the same sphere; VectorDataCubes never invents or relabels one.
An explicit lookup manifold remains authoritative in memory; table conversion
and `setcrs` require Proj to validate a supplied CRS against its radius.
See the [spherical implementation plan](docs/src/spherical-plan.md).

## Tables and attributes

`vectordatacube` turns each attribute column into a layer over a shared geometry
axis.  You can go the other direction using `vectordatacubetable`.

Using the geometries above:

```julia
table = (geometry = geometries, region = ["West", "East"], population = [100, 200])
attributes = vectordatacube(table)

selected = attributes[Geometry(Contains((2.5, 0.5)))]
only(selected[:region])      # "East"
only(selected[:population])  # 200

using DataFrames
columns = DataFrame(vectordatacubetable(cube))
length(columns.Geometry)    # 4: one row per region and year
```

Use `geometrycolumn = :geom` for a differently named geometry column, and
`layers = (:population,)` to keep only selected attributes. 

Table output contains the actual geometry objects in a `Geometry` column. 
`vectordatacubetable` also records the geometry CRS, when present, in 
its parent cube's metadata.

## Zonal statistics

Aggregate a raster over the same regions with `VectorDataCubes.zonal`:

```julia
using Rasters
using Statistics: mean

raster = Raster(
    [Float64(x * t) for x in 1:6, y in 1:2, t in 1:2],
    (X(0.25:0.5:2.75), Y(0.25:0.5:0.75), Ti([2020, 2021]));
    name = :temperature,
)

averages = VectorDataCubes.zonal(mean, raster; of = regions, progress = false)
size(averages)                              # (2, 2): time × geometry
averages[Ti(At(2021)), Geometry(At(geometries[1]))] # 3.0
```

The result retains non-spatial dimensions and adds the geometry axis, so spatial
selectors work on the aggregated values too. A `RasterStack` produces a stack of
cubes. Set `spatialslices = false` to reduce all dimensions to one value per
geometry.

Cropping, masking, and missing-value handling use Rasters' zonal machinery.
Geometries entirely outside the raster produce `missing`; pass
`emptyval = missing` to also handle spatial slices with no valid cells.
Call this function qualified: Rasters exports its own `zonal` function.

## API

| API | Purpose |
| --- | --- |
| `GeometryLookup(data; crs, geometrycolumn, ...)` | Collect geometries and build a spatial index. |
| `Geometry(lookup)` | Attach the lookup to an array dimension. |
| `vectordatacube(table; geometrycolumn, layers, crs)` | Build a stack with one layer per attribute column. |
| `vectordatacubetable(cube)` | Flatten a cube to a Tables.jl-compatible table. |
| `VectorDataCubes.zonal(f, raster; of, ...)` | Aggregate over a lookup and return a vector data cube. |
| `Rasters.reproject(target_crs, lookup)` | Transform geometries and rebuild their spatial index. |

See the [API reference](docs/src/api.md) and the docstrings for
[geometry lookups](src/geometry_lookup.jl), [table conversion](src/tables.jl), and
[zonal statistics](src/zonal.jl) for further options.

## How it works

A vector data cube is a dimensional array or stack whose geometry axis carries
a `GeometryLookup`. The lookup stores a geometry vector and an STRtree spatial
index. Spatial containment and intersection queries first narrow candidates by
extent, then use GeometryOps predicates to check the actual geometries.

The lookup associates `(X(), Y())` coordinates with a single geometry axis.
Geometry selectors and coordinate selector tuples resolve to indices on that
axis. Subsetting keeps the geometries aligned with the array values and rebuilds
the spatial index when the geometry collection changes. Geometry calculations
currently use planar coordinates.

## Going further

The [examples](examples/) cover county observations over time, country-level
zonal statistics, sampling rasters at points, and taxi trips with separate
origin and destination geometry axes. They adapt vector data cube tutorials
from the Python `xvec` and R `stars` ecosystems.

These are also rendered in the docs, so please check those out!

From a checkout, run an example with the documentation environment:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs examples/01_intro_nc_sids.jl
```

Each example downloads its datasets into `examples/data/` on first use.

## AI disclosure

This package was written with the help of generative AI, including Claude and Codex.
