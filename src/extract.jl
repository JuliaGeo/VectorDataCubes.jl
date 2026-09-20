# `VectorDataCubes.extract` is package-owned like `zonal`: Rasters' `extract` returns rows,
# cannot dispatch on a lookup of points, and takes exactly two spatial dimensions (an
# `(X, Y, Ti)` raster is a `MethodError` in Rasters 0.15). Cells come from Rasters' selectors.

"""
    VectorDataCubes.extract(x, points; geometrycolumn = nothing, crs = nokw, manifold = nokw, skipmissing = false, atol = nothing)

Point sampling as a vector data cube: the value of the cell of `x` containing each point, as a
`Raster` (or a `RasterStack`, one layer per layer of `x`) over the non-spatial dimensions of
`x` plus a geometry dimension holding a [`GeometryLookup`](@ref) of the points.

The point counterpart of [`zonal`](@ref VectorDataCubes.zonal): a 2-D `x` gives a vector over
the points, an `(X, Y, Ti)` one a `(Ti, Geometry)` cube.

# Arguments

- `x`: a `Raster` or `RasterStack` with the lookup's dimensions (`X` and `Y` by default).
- `points`, in one of these forms:
  - a vector of point geometries;
  - a table or feature collection with a point geometry column;
  - a `GeometryLookup` of points, giving a result over `Geometry`;
  - a dimension wrapping one (`Dim{:Station}(lookup)`), whose name the result keeps.
  Non-point geometries are an error; [`zonal`](@ref VectorDataCubes.zonal) handles those.

# Keywords

- `geometrycolumn`: the geometry column of a table `points`, when it is not the one
  the table declares.
- `crs`: the CRS of the points. When not given, it is taken from the lookup, the table
  or the geometries, in that order.
- `manifold`: overrides the output geometry lookup's manifold. Defaults to the
  input lookup or table column metadata, otherwise planar. Sampling still uses
  the raster's coordinate axes.
- `skipmissing`: `true` drops points outside `x` or in an all-missing cell (any layer, for a
  stack) from the result and its lookup; `false` (the default) keeps them as `missing`.
- `atol`: the tolerance for matching a point to a cell centre on `Points` lookups; `Intervals`
  lookups use the cell containing the point and ignore it.

The result keeps the name and metadata of `x`. When `x` and the points both
carry a CRS of the same kind and they differ, a warning is emitted; nothing is
reprojected.

Unexported, since Rasters exports an `extract` too: call it qualified, or bind
it with `using VectorDataCubes: extract`.
"""
extract(x::Union{RA.AbstractRaster,RA.AbstractRasterStack}, points;
    geometrycolumn=nothing, crs=nokw, manifold=nokw, skipmissing=false, atol=nothing
) = _extract(x, _pointdim(points, geometrycolumn, crs, manifold); skipmissing, atol)

function _pointdim(points::DD.Dimension{<:GeometryLookup}, geometrycolumn, crs, manifold)
    geomdim = isnokw(crs) ? points : DD.rebuild(points, RA.setcrs(val(points), crs))
    isnokw(manifold) || (geomdim = DD.rebuild(geomdim, DD.rebuild(val(geomdim); manifold)))
    return _checkpoints(geomdim)
end
_pointdim(lookup::GeometryLookup, geometrycolumn, crs, manifold) = _pointdim(Geometry(lookup), geometrycolumn, crs, manifold)
_pointdim(points, geometrycolumn, crs, manifold) =
    _checkpoints(Geometry(GeometryLookup(points; geometrycolumn, crs, manifold)))

function _checkpoints(geomdim)
    geoms = parent(val(geomdim))
    isempty(geoms) && throw(ArgumentError("Cannot extract at an empty `GeometryLookup`."))
    i = findfirst(g -> !(GI.trait(g) isa GI.AbstractPointTrait), geoms)
    isnothing(i) || throw(ArgumentError(
        "`extract` samples `x` at points, but geometry $i has trait " *
        "`$(nameof(typeof(GI.trait(geoms[i]))))`; use `VectorDataCubes.zonal` to " *
        "aggregate over lines and polygons."
    ))
    return geomdim
end

function _extract(x, geomdim; skipmissing, atol)
    lookup = val(geomdim)
    xydims = _lookupdims(x, lookup)
    _warn_crs_mismatch(x, lookup)
    cells = map(p -> _cell(xydims, p, atol), parent(lookup))
    if skipmissing
        keep = map(c -> !isnothing(c) && !_missingcell(x, xydims, c), cells)
        cells = cells[keep]
        geomdim = DD.rebuild(geomdim, lookup[keep])
    end
    return _cellcube(x, xydims, cells, geomdim)
end

# The index of the cell containing `point`, or `nothing` when it lies outside
# `xydims`. Rasters' own point sampling: `Contains` for `Intervals` lookups,
# `At` within `atol` for `Points` lookups.
function _cell(xydims, point, atol)
    selectors = map(d -> RA._at_or_contains(d, RA._dimcoord(d, point), atol), xydims)
    DD.hasselection(xydims, selectors) || return nothing
    return DD.dims2indices(xydims, selectors)
end

_celldims(xydims, cell) = map(DD.rebuild, xydims, cell)

# Rasters' `skipmissing` skips `missing` and the raster's `missingval` alike.
_missingcell(x::RA.AbstractRaster, xydims, cell) =
    isempty(skipmissing(view(x, _celldims(xydims, cell)...)))
_missingcell(st::RA.AbstractRasterStack, xydims, cell) =
    any(layer -> _missingcell(layer, xydims, cell), DD.layers(st))

function _cellcube(x::RA.AbstractRaster, xydims, cells, geomdim)
    zs = map(c -> isnothing(c) ? missing : x[_celldims(xydims, c)...], cells)
    return _geometry_cube(x, zs, geomdim, DD.otherdims(x, xydims))
end
_cellcube(st::RA.AbstractRasterStack, xydims, cells, geomdim) =
    DD.maplayers(A -> _cellcube(A, xydims, cells, geomdim), st)
