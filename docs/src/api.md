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

A `GeometryLookup` resolves DimensionalData's selectors spatially: each one narrows the
candidates with the spatial tree, then refines them with an exact GeometryOps predicate.
On the `Geometry` axis, `Contains`, `At`, `Near`, `Touches`, `Where` and the predicates
of DE9IM.jl are accepted; on the internal `X`/`Y` axes, pairs such as
`cube[X(At(x)), Y(At(y))]` or `cube[X(a .. b), Y(c .. d)]`.

```@meta
# The signatures below name `Lookups`, `Extents` and `DE9IM`, which are in scope here.
CurrentModule = VectorDataCubes
```

```@docs
Lookups.selectindices(::GeometryLookup, ::Lookups.Contains)
Lookups.selectindices(::GeometryLookup, ::Lookups.At)
Lookups.selectindices(::GeometryLookup, ::Lookups.Near)
Lookups.selectindices(::GeometryLookup, ::Lookups.Touches{<:Extents.Extent})
Lookups.selectindices(::GeometryLookup, ::Lookups.Where)
Lookups.selectindices(::GeometryLookup, ::DE9IM.DE9IMPredicate)
Lookups.selectindices(::GeometryLookup, ::Tuple)
VectorDataCubes.mask
```

## Coordinate reference systems

A lookup carries its crs (`GeoInterface.crs`, set with `Rasters.setcrs`); with Proj
loaded it can be reprojected.

```@docs
RA.reproject(::RA.GeoFormat, ::GeometryLookup)
```

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
