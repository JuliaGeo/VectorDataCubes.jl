```@raw html
---
layout: home

hero:
  name: VectorDataCubes.jl
  text: Spatial indexing for vector data cubes
  tagline: A DimensionalData / Rasters array whose Geometry dimension is backed by a spatial index — so cube[Geometry(Contains(point))] just works.
  actions:
    - theme: brand
      text: Get started
      link: /examples/01_intro_nc_sids
    - theme: alt
      text: API reference
      link: /api
    - theme: alt
      text: View on GitHub
      link: https://github.com/asinghvi17/VectorDataCubes.jl

features:
  - title: Spatial selectors
    details: Index a cube by geometry — cube[Geometry(Contains(point))], or plain X(a..b), Y(c..d) on the same axis — narrowed by an R-tree and refined by exact GeometryOps predicates.
  - title: Zonal statistics → cubes
    details: VectorDataCubes.zonal aggregates a raster over geometries into a (Ti, Geometry) cube you can keep slicing spatially.
  - title: Tables ↔ cubes
    details: vectordatacube lifts a GeoJSON / Shapefile / DataFrame into a cube; vectordatacubetable flattens it back, with real geometries and the CRS preserved.
---
```

## What is a vector data cube?

A **vector data cube** is a [DimensionalData](https://github.com/rafaqz/DimensionalData.jl) /
[Rasters](https://github.com/rafaqz/Rasters.jl) array or stack whose `Geometry`
dimension is backed by a [`GeometryLookup`](@ref) — a vector of geometries plus a
packed R-tree spatial index (`GeometryOps.FlexibleRTrees`). You get everything
DimensionalData and Rasters already give you (named dimensions, selectors,
broadcasting, table conversion), **plus** spatial indexing on the geometry axis:

```julia
cube[Geometry(Contains(point))]   # the geometry that contains a point
cube[X(a .. b), Y(c .. d)]        # the same axis, selected by bounding box
```

The geometry lookup spans both the `Geometry` dimension *and* the `(X(), Y())` it
wraps, which is why geometry selectors and `X`/`Y` selectors both resolve against the
one axis.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/asinghvi17/VectorDataCubes.jl")
```

## Where to go next

- The [examples](examples/01_intro_nc_sids.md) walk through building cubes, spatial
  selectors, zonal statistics, point extraction, and cubes with two geometry axes —
  ported from the Python (`xvec`) and R (`stars`) tutorials.
- The [API reference](api.md) documents the public surface: [`GeometryLookup`](@ref),
  [`vectordatacube`](@ref), [`vectordatacubetable`](@ref), and
  [`VectorDataCubes.zonal`](@ref).
