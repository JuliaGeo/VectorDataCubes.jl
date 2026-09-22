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
2. Implement spherical nearest-point selection with an exhaustive scan, then
   accelerate it using XYZ chord bounds converted to angular distance. General
   spherical point-to-line/polygon distance remains a follow-up.
3. Round-trip edge/orientation semantics through table column metadata. Preserve
   the supplied CRS unchanged. With Proj loaded, a spherical CRS keeps its declared
   radius exactly; an ellipsoidal datum uses the arithmetic mean radius `(2a+b)/3`.
   A custom radius must agree with the CRS, and without a CRS only the GeometryOps
   default radius can round-trip.
4. Preserve the manifold through point extraction. Reject spherical zonal burning
   and reprojection until their space-conversion semantics are implemented.

The implementation uses released GeometryOps APIs and explicit manifold bounds,
without depending on an unmerged tree-query API. Proj is an optional extension:
GeoParquet `edges` still chooses the manifold independently, while the CRS supplies
and validates its physical radius.

## Verification

Compare indexed selection against exhaustive spherical predicates for antimeridian
and polar geometries, boundaries, holes, crossings, overlaps, Float32 input, and
geodesic boxes. Check lazy cache invalidation, slicing/view/reverse, unchanged
public bounds after index construction, nearest-point ties, mixed-manifold table
columns, metadata round-trips, and unchanged longitude/latitude output. Run the
existing planar test suite.

## Follow-ups

- General spherical distance; point lookup tree pruning is implemented using an
  XYZ-box chord lower bound and exhaustive-scan equivalence tests.
- Explicit reprojection and spherical rasterization/zonal semantics.
- USP input conversion at the public boundary and format-specific metadata bridges.

Relevant upstream work: [GeoDataFrames #167](https://github.com/evetion/GeoDataFrames.jl/pull/167),
[GeometryOps #506](https://github.com/JuliaGeo/GeometryOps.jl/pull/506), and
[GeometryOps #504](https://github.com/JuliaGeo/GeometryOps.jl/pull/504).
