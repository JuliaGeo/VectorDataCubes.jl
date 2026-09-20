# Spherical geometry plan

VectorDataCubes will use a lookup's manifold consistently for its spatial index
and exact predicates. Public geometries remain in their input coordinates;
unit-sphere XYZ bounds and prepared geometries are internal accelerators.

## Agreed design

- Expose the existing typed `GeometryLookup.manifold`. Planar remains the default;
  spherical input uses longitude/latitude in degrees and great-circle edges.
- Keep lazy R-trees: planar lookups index XY bounds, spherical lookups index
  conservative XYZ bounds from `GeometryOps.extent`. Rebuild the cache when data
  or manifold changes. Public dimension bounds remain longitude/latitude.
- Use manifold-aware predicates, including RelateNG for spherical `crosses` and
  `overlaps`. Keep `At(geometry)` coordinate equality, boundary-inclusive
  `Contains`, sorted selections, and lowest-index nearest ties.
- Interpret finite interval boxes as polygons with great-circle edges in the
  first cut. Coordinate windows with constant-latitude boundaries and cross-space
  operations are later work. Reject unbounded or ambiguous spherical boxes.
- Follow [GeoParquet's column metadata semantics](https://geoparquet.org/releases/v1.1.0/):
  `edges` is `planar` or `spherical`; `orientation`, when asserted, is
  `counterclockwise`. Missing `edges` means planar. Preserve these per geometry
  column through DataAPI column metadata.
- Continue emitting the existing DataAPI table keys `GEOINTERFACE:crs` and
  `GEOINTERFACE:geometrycolumns`. Datum and sphere/ellipsoid parameters belong to
  the CRS. Edge and ring-interior interpretation are separate geometry semantics.
- Emit longitude/latitude geometry, never internal unit-sphere coordinates under
  a geographic CRS. File writers remain responsible for mapping DataAPI metadata
  into their supported formats; this package does not become a GeoParquet writer.

## Small first cut

1. Add explicit manifold construction and preservation, XY/XYZ lazy indexing,
   manifold-aware predicates, and finite geodesic interval selection.
2. Implement spherical nearest-point selection with an exhaustive scan. General
   spherical point-to-line/polygon distance and tree pruning remain follow-ups.
3. Round-trip edge/orientation semantics through table column metadata. Preserve
   supplied CRS unchanged. Until datum-to-radius resolution is implemented, reject
   table export of a nondefault spherical radius rather than silently losing it.
4. Preserve the manifold through point extraction. Reject spherical zonal burning
   and reprojection until their space-conversion semantics are implemented.

The implementation should use released GeometryOps APIs and explicit manifold
bounds, without depending on an unmerged tree-query API or adding CRS inference.

## Verification

Compare indexed selection against exhaustive spherical predicates for antimeridian
and polar geometries, boundaries, holes, crossings, overlaps, Float32 input, and
geodesic boxes. Check lazy cache invalidation, slicing/view/reverse, unchanged
public bounds after index construction, nearest-point ties, mixed-manifold table
columns, metadata round-trips, and unchanged longitude/latitude output. Run the
existing planar test suite.

## Follow-ups

- Resolve sphere parameters from the CRS datum, and preserve custom datums through
  interoperable output. Do not relabel WGS84 coordinates merely to select a
  spherical computational approximation.
- General spherical distance and a proven XYZ-box lower bound for nearest queries.
- Explicit reprojection and spherical rasterization/zonal semantics.
- USP input conversion at the public boundary and format-specific metadata bridges.

Relevant upstream work: [GeoDataFrames #167](https://github.com/evetion/GeoDataFrames.jl/pull/167),
[GeometryOps #506](https://github.com/JuliaGeo/GeometryOps.jl/pull/506), and
[GeometryOps #504](https://github.com/JuliaGeo/GeometryOps.jl/pull/504).
