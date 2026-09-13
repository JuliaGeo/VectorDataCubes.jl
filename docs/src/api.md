# API reference

## Index

```@index
```

## Building cubes

```@docs
GeometryLookup
Geometry
vectordatacube
vectordatacubetable
VectorDataCubes.VectorDataCubeTable
VectorDataCubes.spatialtree
```

## Selectors

A [`GeometryLookup`](@ref) resolves DimensionalData's selectors spatially: each one narrows
the candidates with the spatial tree, then refines them with an exact GeometryOps predicate.
[`VectorDataCubes.mask`](@ref) turns any of them into a `Bool` mask.

On the `Geometry` axis:

- `Contains(point)`: every geometry covering the point; closed, so a point on a shared
  border belongs to each geometry it lies on.
- `At(geom)`: the geometry equal to `geom`; `Near(point)`: the nearest geometry, ties to
  the lowest index.
- `Touches(extent)`: the geometries intersecting the extent.
- `Where(f)`: the geometries `f` holds for; a curried GeometryOps predicate such as
  `GO.intersects(geom)` is narrowed by the tree first.
- A DE9IM.jl predicate, `DE9IM.Covers(geom)` and the like.

On the internal `X`/`Y` axes, pairs: `At`/`At`, `Contains`/`Contains`, `Near`/`Near`, two
intervals (geometries covered by the box) and `Touches`/`Touches` (geometries intersecting
it), as `cube[X(At(x)), Y(At(y))]` or `cube[X(a .. b), Y(c .. d)]`.

```@docs
VectorDataCubes.mask
```

## Coordinate reference systems

A lookup carries its crs (`GeoInterface.crs`, set with `Rasters.setcrs`), and
`Rasters.reproject` reprojects it through `GeometryOps.reproject`, which needs Proj.jl loaded.
Reprojecting a lookup with no crs is an `ArgumentError`; set one first.

## Zonal statistics

`zonal` is deliberately **not** exported (Rasters exports a `zonal` too); call it
qualified as `VectorDataCubes.zonal`, or bind it with `using VectorDataCubes: zonal`.

```@docs
VectorDataCubes.zonal
```

## Point extraction

`extract` is likewise **not** exported (Rasters exports an `extract` too); call it
qualified as `VectorDataCubes.extract`, or bind it with `using VectorDataCubes: extract`.

```@docs
VectorDataCubes.extract
```

## Plotting

With Makie loaded, a `GeometryLookup` — and any dimension wrapping one, such as the
`Geometry` dimension of a cube — is a plottable object: `poly`, `lines` and `scatter`
convert its geometries to GeometryBasics and hand them to Makie's own recipes, and
`plot` picks the recipe from the kind of the first geometry:

- polygons, multipolygons and linear rings to `poly`,
- linestrings and multilinestrings to `lines`,
- points and multipoints to `scatter`.

Pass one value per geometry as `color` for a choropleth:

```julia
using Makie
poly(dims(cube, Geometry); color = cube[Ti = 1])
```

A lookup mixing single and multi geometries of one family (polygon/multipolygon,
linestring/multilinestring, point/multipoint) is lifted to the multi kind before
plotting. Any other mixture, and an empty lookup, is an `ArgumentError`.

A cube passed to `plot` goes to DimensionalData's own Makie extension, which has no
notion of a geometry axis and fails to convert it. Plot the `Geometry` dimension and
pass the cube's values as `color`, as above.
