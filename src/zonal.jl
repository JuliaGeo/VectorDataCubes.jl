# `VectorDataCubes.zonal` is package-owned, not a method of `Rasters.zonal`,
# which cannot dispatch on its `of` keyword: geometry-lookup `of`s are handled
# here, everything else forwards to `Rasters.zonal`.

"""
    VectorDataCubes.zonal(f, x; of, spatialslices = true, kw...)

Zonal statistics as a vector data cube. Like `Rasters.zonal`, `f` reduces the cells of `x`
covered by each geometry; with a [`GeometryLookup`](@ref) `of`, the result is a `Raster` (or
a `RasterStack`, one layer per layer of `x`) over that lookup, ready for spatial selectors.

# Arguments

- `f`: a function reducing an iterable to one value, such as `sum` or `Statistics.mean`.
- `x`: a `Raster` or `RasterStack` with the lookup's dimensions (`X` and `Y` by default).

# Keywords

- `of`: the geometries, in one of these forms:
  - a `GeometryLookup`, giving a result over `Geometry`;
  - a dimension wrapping one (`Geometry(lookup)`, `Dim{:Origin}(lookup)`), whose name
    the result keeps;
  - a `DimArray`, `DimStack`, `Raster` or `RasterStack` with exactly one such dimension,
    or a tuple of dimensions containing one;
  - anything else, forwarded to `Rasters.zonal` unchanged.
- `spatialslices`: the dimensions `f` reduces over when `x` has more dimensions than
  the lookup spans (`Ti` or `Band` on top of `X` and `Y`, say):
  - `true` (the default): the lookup's dimensions, so `f` sees one spatial slice at a time
    (`mapslices(f, masked; dims = (X, Y))`) and the result spans the other dimensions too;
  - `false`: every dimension, so `f` sees the whole masked raster per geometry and the
    result is a vector over the geometry dimension;
  - a tuple of dimensions containing the lookup's, such as `(X, Y, Ti)`: those, with the result
    over the rest. The lookup's dimensions are required: each crop has its own spatial size.
- `emptyval`: the value for a geometry (or slice, when slicing) covering no non-missing cell
  under `skipmissing`; without it, `f` is called on an empty iterator.
- `skipmissing`, `progress`, `threaded`, `boundary`, `shape`: as in `Rasters.zonal`.

A geometry entirely outside `x` gives `missing` (a `missing`-filled slice when slicing). The
result keeps the name and metadata of `x`. A CRS mismatch between `x` and the lookup (same
CRS kind, different value) only warns; nothing is reprojected.

Unexported, since Rasters exports a `zonal` too: call it qualified, or bind it
with `using VectorDataCubes: zonal`.
"""
zonal(f, x; of, kw...) = _zonal(f, x, of; kw...)

_zonal(f, x, of; kw...) = RA.zonal(f, x; of, kw...)
_zonal(f, x, lookup::GeometryLookup; kw...) = _zonal(f, x, Geometry(lookup); kw...)
_zonal(f, x, of::Union{DD.AbstractDimArray,DD.AbstractDimStack,DD.DimTuple}; kw...) =
    _zonal(f, x, _geometrydim(of), of; kw...)
_zonal(f, x, ::Nothing, of; kw...) = RA.zonal(f, x; of, kw...)
_zonal(f, x, geomdim::DD.Dimension, of; kw...) = _zonal(f, x, geomdim; kw...)

function _zonal(f, x::Union{RA.AbstractRaster,RA.AbstractRasterStack},
    geomdim::DD.Dimension{<:GeometryLookup}; kw...
)
    lookup = val(geomdim)
    isempty(parent(lookup)) &&
        throw(ArgumentError("Cannot compute zonal statistics with an empty `GeometryLookup`."))
    _warn_crs_mismatch(x, lookup)
    return _zonal_geometries(f, x, geomdim; kw...)
end

# Stacks fan out by layer, so layers with different dimensions each produce
# a cube of the right shape.
function _zonal_geometries(f, st::RA.AbstractRasterStack, geomdim; kw...)
    K = keys(st)
    layers = map(K) do k
        _zonal_geometries(f, st[k], geomdim; kw...)
    end
    return RA.RasterStack(NamedTuple{K}(layers); metadata=DD.metadata(st))
end
function _zonal_geometries(f, x::RA.AbstractRaster, geomdim;
    spatialslices=true, skipmissing=true, emptyval=nokw, progress=true, threaded=true, kw...
)
    lookup = val(geomdim)
    geoms = parent(lookup)
    slicedims = _zonal_slicedims(spatialslices, x, _lookupdims(x, lookup))
    return Base.open(x) do o
        xp = RA._prepare_for_burning(o)
        zs = if isnothing(slicedims)
            RA._zonal(f, xp, nothing, geoms; skipmissing, emptyval, progress, threaded, kw...)
        else
            # `emptyval` applies per slice in the wrapper, so an all-empty geometry still yields
            # a slice-shaped result. Slice eltypes can differ between geometries (all-`emptyval`
            # for a sub-cell one), so collect untyped; Rasters' loop types from the first.
            inner = _SpatialSliceify(f, DD.dims(xp, slicedims), emptyval)
            desc = "Applying $f to each geometry..."
            _zonal_eachgeom(inner, xp, geoms, desc; skipmissing, progress, threaded, kw...)
        end
        otherdims = isnothing(slicedims) ? () : DD.otherdims(xp, slicedims)
        _geometry_cube(xp, zs, geomdim, otherdims)
    end
end

# The single dimension of `x` backed by a `GeometryLookup`, or `nothing`.
function _geometrydim(x)
    geomdims = filter(d -> val(d) isa GeometryLookup, DD.dims(x))
    length(geomdims) <= 1 || throw(ArgumentError(
        "`of` has $(length(geomdims)) dimensions backed by a `GeometryLookup` " *
        "($(join(map(DD.name, geomdims), ", "))); pass the one to use, e.g. `of = dims(of, Geometry)`."
    ))
    return isempty(geomdims) ? nothing : only(geomdims)
end

# The dimensions of `x` the lookup spans, in the lookup's order.
function _lookupdims(x, lookup)
    xydims = DD.dims(x, DD.dims(lookup))
    length(xydims) == 2 || throw(ArgumentError(
        "The `GeometryLookup` spans dimensions $(map(DD.name, DD.dims(lookup))) " *
        "but `x` has dimensions $(map(DD.name, DD.dims(x))); both of the lookup's dimensions must be present in `x`."
    ))
    return xydims
end

function _warn_crs_mismatch(x, lookup)
    xcrs, lcrs = RA.crs(x), GI.crs(lookup)
    (isnothing(xcrs) || isnothing(lcrs)) && return nothing
    DD.basetypeof(xcrs) === DD.basetypeof(lcrs) || return nothing
    xcrs == lcrs && return nothing
    @warn "`x` has crs $(xcrs) but the `GeometryLookup` has crs $(lcrs); values are computed as if both were the same. Reproject one of them first."
    return nothing
end

_zonal_slicedims(spatialslices::Bool, x, xydims) = spatialslices ? xydims : nothing
function _zonal_slicedims(spatialslices::Tuple, x, xydims)
    slicedims = DD.dims(x, spatialslices)
    length(slicedims) == length(spatialslices) || throw(ArgumentError(
        "`spatialslices = $spatialslices` names dimensions `x` does not have; `x` has $(map(DD.name, DD.dims(x)))."
    ))
    length(DD.dims(slicedims, xydims)) == 2 || throw(ArgumentError(
        "`spatialslices = $spatialslices` must contain the lookup's dimensions $(map(DD.name, xydims)): " *
        "each geometry's crop has its own spatial size, so per-geometry results cannot be stacked along them."
    ))
    return slicedims
end
_zonal_slicedims(spatialslices, x, xydims) = throw(ArgumentError(
    "`spatialslices` must be `true`, `false` or a tuple of dimensions containing the lookup's; got `$spatialslices`."
))

# Like Rasters' `_zonal(f, x, ::Nothing, geoms)`, reusing its per-geometry
# crop/mask path and `_run` threading/progress, but collecting into an
# untyped vector that is narrowed afterwards.
function _zonal_eachgeom(f, x, geoms, desc; skipmissing, progress, threaded, kw...)
    zs = Vector{Any}(undef, length(geoms))
    RA._run(eachindex(zs), threaded, progress, desc) do i
        zs[i] = RA._zonal(f, x, geoms[i]; skipmissing, emptyval=nokw, kw...)
    end
    return map(identity, zs)
end

# Wraps `f` to reduce each spatial slice, returning a `Raster` over the remaining dims. With
# `skipmissing=true` Rasters passes the wrapper `skipmissing(masked)`; that is unwrapped and
# `skipmissing` re-applied per slice.
struct _SpatialSliceify{F,D,E}
    f::F
    dims::D
    emptyval::E
end

(s::_SpatialSliceify)(x::DD.AbstractDimArray) =
    _mapspatialslices(_empty_aware(s.f, s.emptyval), x, s.dims)
(s::_SpatialSliceify)(sm::Base.SkipMissing) =
    _mapspatialslices(_empty_aware(s.f, s.emptyval) ∘ Base.skipmissing, sm.x, s.dims)
(s::_SpatialSliceify)(sm::RA.SkipMissingVal) =
    _mapspatialslices(_empty_aware(s.f, s.emptyval) ∘ Base.skipmissing, sm.x, s.dims)

# If `emptyval` was passed, return it for empty (e.g. fully-masked) slices
# instead of calling `f` on an empty iterator.
function _empty_aware(f, emptyval)
    isnokw(emptyval) && return f
    return el -> isempty(el) ? emptyval : f(el)
end

function _mapspatialslices(g, x::DD.AbstractDimArray, slicedims)
    otherdims = DD.otherdims(x, slicedims)
    isempty(otherdims) && return g(x)
    slices = eachslice(x; dims=otherdims)
    return DD.rebuild(x; data=[g(slice) for slice in slices], dims=DD.dims(slices), refdims=())
end

# Assemble the per-geometry results (scalars, `Raster`s when slicing, or
# `missing` for geometries entirely outside the raster) into a vector data
# cube along `geomdim`, keeping the name and metadata of `x`.
function _geometry_cube(x::RA.AbstractRaster, zs::AbstractVector, geomdim::DD.Dimension, otherdims::Tuple)
    name, metadata = DD.name(x), DD.metadata(x)
    if isempty(zs)
        data = Array{eltype(x)}(undef, length.(otherdims)..., 0)
        return RA.Raster(data, (otherdims..., geomdim); name, metadata)
    end
    i = findfirst(z -> z isa DD.AbstractDimArray, zs)
    if isnothing(i)
        # Scalar results: a vector over the geometry dimension only...
        isempty(otherdims) && return RA.Raster(zs, (geomdim,); name, metadata)
        # ...unless slicing over `otherdims` was requested and every geometry
        # was outside the raster - then keep the cube shape, so the output
        # dimensionality doesn't depend on data coverage.
        data = Base.stack(map(z -> fill(z, length.(otherdims)), zs))
        return RA.Raster(data, (otherdims..., geomdim); name, metadata)
    end
    # Geometries entirely outside the raster came back as `missing` and are
    # expanded to missing-filled slices.
    template = zs[i]
    arrays = map(zs) do z
        z isa DD.AbstractDimArray ? parent(z) : fill(z, size(template))
    end
    data = Base.stack(arrays)
    return RA.Raster(data, (DD.dims(template)..., geomdim); name, metadata)
end
