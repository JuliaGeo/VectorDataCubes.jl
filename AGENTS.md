# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

VectorDataCubes.jl makes a **vector data cube**: a `DimensionalData`/`Rasters` array or
stack whose `Geometry` dimension is backed by a `GeometryLookup` (a geometry vector + a
packed R-tree). The payoff is spatial indexing — `cube[Geometry(Contains(point))]`,
`cube[X(a..b), Y(c..d)]` — on top of everything DimensionalData/Rasters already gives you.

Five `src/` files, one concern each:

- `geometry_lookup.jl` — the `Geometry` dimension and `GeometryLookup` (the spatial-indexing core:
  struct, constructor, lazily built tree via `spatialtree`, DD interface, crs, `reproject`).
- `selectors.jl` — every `Lookups.selectindices` method for a `GeometryLookup`.
- `zonal.jl` — `VectorDataCubes.zonal`, aggregating a raster over geometries into a cube.
- `extract.jl` — `VectorDataCubes.extract`, sampling a raster at points into a cube.
- `tables.jl` — `vectordatacube` / `vectordatacubetable`, round-tripping a table ↔ cube;
  `VectorDataCubeTable` (the returned table) carries geometry columns and crs as DataAPI metadata.

Plus one package extension, `ext/VectorDataCubesMakieExt.jl` (loaded with Makie), which
makes a `GeometryLookup`, and any dimension wrapping one, plottable by converting its
geometries to GeometryBasics and forwarding to Makie's `poly`/`lines`/`scatter` recipes.
Its `convert_arguments`/`plottype` methods name the lookup type, which is more specific
than the `AbstractArray{<:SomeGeometry}` methods the geometry packages install through
`GeoInterface.@enable_makie`, so one conversion covers every element type. Nothing there
touches `DimArray`s — DimensionalData's own Makie extension owns those. Reprojection needs
no extension: `Rasters.reproject` on a `GeometryLookup` calls `GeometryOps.reproject`, which
needs Proj.jl loaded and says so itself (a `MethodError` carrying GeometryOps' hint).

## Commands

Julia workspace package: root `Project.toml` has `[workspace] projects = ["docs", "test"]`; `test/` and `docs/` carry their own `Project.toml`.  When running code that wants to use this package, use the `docs/` environment, since that contains the ecosystem packages too.

```sh
# Full suite
julia --project=. -e 'using Pkg; Pkg.test()'
# Faster iteration (no sandbox build)
julia --project=test test/runtests.jl
# A single test file — each is self-contained (does its own imports), so just include it.
julia --project=test -e 'include("test/zonal.jl")'
# An example (each self-downloads its data to examples/data/)
julia --project=docs examples/02_zonal_countries.jl
```

`runtests.jl` includes each test file via `SafeTestsets.@safetestset`, so every file runs
in its own module — each must do its own `using VectorDataCubes` (and other imports);
nothing leaks in from `runtests.jl`. A new test file must follow suit.

`test/basics.jl` and the examples need network access. CI runs Julia 1.10 and `1`.

For iterative work, prefer the persistent Julia REPL via the `mcp__julia__*` tools — it
avoids re-paying Julia's per-process compile latency on every run.

## Architecture

- **`Lookups.selectindices` (in `selectors.jl`) is the heart.** Every spatial
  selector resolves to indices there: narrow with an R-tree extent query
  (`_maybe_get_candidates`), then refine with an exact GeometryOps predicate
  (`_select_predicate`). A new selector means a new method here. Supported today:
  `Contains(point)` (closed, `GO.covers`), `At(geom)`, `Near(point)` (branch-and-bound
  over the tree, ties to the lowest index so a lookup with a tree and one without answer
  the same), `Touches(extent)` (`GO.intersects`), `Where(f)` with tree-narrowed fast
  paths for every curried `GO.pred(g)`, the DE9IM.jl predicates (`DE9IM.Covers(g)`, …;
  never re-exported, their names clash with DD's selectors), and on the internal `X`/`Y`
  dims the pairs `At`/`At`, `Contains`/`Contains`, `Near`/`Near`, interval/interval
  (covered by the box) and `Touches`/`Touches` (intersecting the box), matched to
  coordinates by dimension type. `VectorDataCubes.mask(lookup_or_cube, sel)` turns any of
  them into a `Bool` mask.
- **The lookup spans `(X(), Y())` *and* the `Geometry` dim wrapping it** — that's why both
  `Geometry(...)` and `X()/Y()` selectors work on one axis.
- **The tree is lazy.** `spatialtree(lookup)` builds it on the first spatial query and
  caches it; slicing, `view`, `reverse` and `DD.rebuild` with new data hand back an
  unbuilt index, so the array type never depends on the number of geometries. Nothing
  outside `geometry_lookup.jl` touches `lookup.tree` or `lookup.data`. It always indexes
  `Float64` `X`/`Y` extents (`_xyextent`) with a `Vector{Int}` of leaf indices, so its
  type — `XYRTree{algorithm, geometryvector}` — follows from the lookup's own type,
  whatever the geometries' coordinate type or dimensionality.
- **`zonal` and `extract` are package-owned, not methods of the Rasters ones** (Rasters
  can't dispatch on `zonal`'s `of`, and its `extract` returns rows), and not exported
  (call them qualified). A `GeometryLookup` `of` yields a cube; anything else forwards to
  `Rasters.zonal`. Both assemble their result with `_geometry_cube`, which keeps the name
  and metadata of `x` and gives a cube over `(otherdims..., geometry dim)`.

## Conventions

- **Import aliases**, consistent everywhere: `DD`, `GO`, `GOCore`, `GI`, `RA`, plus
  `Extents`, `Missings`.
- **`nokw` / `isnokw`** (from Rasters) is the "keyword not supplied" sentinel, distinct
  from a meaningful `nothing` (e.g. no CRS / no tree).
