import DE9IM

const STI = GO.SpatialTreeInterface

# Every selector on a `GeometryLookup` resolves here: narrow candidates with an extent query
# on `spatialtree(lookup)` (all indices without a tree), then refine with an exact GeometryOps
# predicate. Multi-index results are `Vector{Int}`; `At`, `Near` and `(At, At)` return an `Int`.

# Dimension-wrapped selectors: sort into the lookup's internal dim order, then pair up.
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple)
    Lookups.selectindices(lookup, map(_val_or_nothing, DD.sortdims(sel, DD.dims(lookup))))
end
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    Lookups.selectindices(lookup, map(DD.rebuild, map(DD.name2dim, K), values(sel)))
end
_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)

function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    length(sel) == 2 || _pair_error(sel)
    return _select_pair(lookup, _xy(lookup, sel)...)
end
# The pair follows the lookup's own dim order, so `dims = (Y(), X())` pairs up too.
function _xy(lookup::GeometryLookup, sel::Tuple)
    first(DD.dims(lookup)) isa DD.XDim ? (sel[1], sel[2]) : (sel[2], sel[1])
end
_pair_error(sel) = throw(ArgumentError("""
    Unsupported selector pair `$sel` on a `GeometryLookup`. Supported: `(X(At(x)), Y(At(y)))`,
    `(X(Contains(x)), Y(Contains(y)))`, `(X(Near(x)), Y(Near(y)))`, `(X(a .. b), Y(c .. d))`,
    `(X(Touches(a, b)), Y(Touches(c, d)))`, and an interval or `Touches` on one axis alone.
    """))

const _PairInterval = Union{DD.IntervalSets.Interval,Nothing}
const _PairTouches = Union{Lookups.Touches{<:Tuple{Real,Real}},Nothing}

_select_pair(::GeometryLookup, x, y) = _pair_error((x, y))
_select_pair(::GeometryLookup, ::Nothing, ::Nothing) = _pair_error((nothing, nothing))
function _select_pair(lookup::GeometryLookup, x::Lookups.At, y::Lookups.At)
    point = (val(x), val(y))
    i = _first_covering(lookup, point)
    isnothing(i) && throw(ArgumentError("No geometry in the lookup covers the point $point."))
    return i
end
_select_pair(lookup::GeometryLookup, x::Lookups.Contains, y::Lookups.Contains) =
    _select_predicate(lookup, GO.covers, (val(x), val(y)))
_select_pair(lookup::GeometryLookup, x::Lookups.Near, y::Lookups.Near) =
    _nearest(lookup, (val(x), val(y)))
function _select_pair(lookup::GeometryLookup, x::_PairInterval, y::_PairInterval)
    box = Extents.Extent(X=_axisbounds(x), Y=_axisbounds(y))
    return _select_predicate(lookup, GO.coveredby, box)
end
function _select_pair(lookup::GeometryLookup, x::_PairTouches, y::_PairTouches)
    box = Extents.Extent(X=_axisbounds(x), Y=_axisbounds(y))
    return _select_predicate(lookup, GO.intersects, box)
end
_axisbounds(::Nothing) = (-Inf, Inf)
_axisbounds(i::DD.IntervalSets.Interval) = extrema(i)
_axisbounds(t::Lookups.Touches) = extrema(val(t))

# Single selectors on the `Geometry` axis

# Closed point-in-polygon (`GO.covers`): a point on a shared border selects every geometry
# whose boundary it lies on.
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Contains) =
    _select_predicate(lookup, GO.covers, _checked_geometry(sel))
function Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Contains{<:AbstractVector})
    return _union(lookup, val(sel)) do v
        Lookups.selectindices(lookup, DD.rebuild(sel; val=v))
    end
end
Lookups.selectindices(lookup::GeometryLookup,
    sel::Lookups.Contains{<:GO.UnitSpherical.UnitSphericalPoint}) =
    _select_predicate(lookup, GO.covers, _checked_geometry(sel))

Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.At) = _select_at(lookup, sel)
Lookups.selectindices(lookup::GeometryLookup,
    sel::Lookups.At{<:GO.UnitSpherical.UnitSphericalPoint}) = _select_at(lookup, sel)
function _select_at(lookup, sel)
    geom = _checked_geometry(sel)
    i = _at(lookup, geom)
    isnothing(i) && throw(ArgumentError("No geometry in the lookup equals `$(sel)`."))
    return i
end
function _at(lookup::GeometryLookup, geom)
    geom = _querygeometry(lookup.manifold, geom)
    candidates = _maybe_get_candidates(lookup, geom)
    geoms = parent(lookup)
    k = findfirst(i -> GO.equals(geoms[i], geom), candidates)
    return isnothing(k) ? nothing : candidates[k]
end

Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Near) = _nearest(lookup, _checked_point(sel))
Lookups.selectindices(lookup::GeometryLookup,
    sel::Lookups.Near{<:GO.UnitSpherical.UnitSphericalPoint}) =
    _nearest(lookup, _checked_point(sel))
# Same body as DimensionalData's `selectindices(::Lookup, ::Selector{<:AbstractVector})`, but
# needed: the `At`/`Near` methods above and that generic one are ambiguous for a vector value.
Lookups.selectindices(lookup::GeometryLookup, sel::Union{Lookups.At{<:AbstractVector},Lookups.Near{<:AbstractVector}}) =
    Int[Lookups.selectindices(lookup, DD.rebuild(sel; val=v)) for v in val(sel)]

# DimensionalData's `Touches` is the loose relation (`GO.intersects`); the strict DE-9IM one
# is `DE9IM.Touches(geom)`.
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Touches{<:Extents.Extent}) =
    _select_predicate(lookup, GO.intersects, val(sel))
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Touches) = throw(ArgumentError("""
    `Touches` on a `GeometryLookup` wraps an `Extents.Extent`; got `$sel`. For the
    bounds form use it per axis: `(X(Touches(a, b)), Y(Touches(c, d)))`.
    """))

Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Where) = findall(val(sel), parent(lookup))

for pred in (:equals, :intersects, :contains, :within, :covers, :coveredby, :touches,
             :crosses, :overlaps, :disjoint)
    @eval Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Where{<:Base.Fix2{typeof(GO.$pred)}}) =
        _select_predicate(lookup, GO.$pred, val(sel).x)
end

function Lookups.selectindices(lookup::GeometryLookup, sel::DE9IM.DE9IMPredicate)
    geom = _checked_geometry(sel)
    pred = _predicate(sel)
    _isgeometry(geom) && return _select_predicate(lookup, pred, geom)
    return _union(g -> _select_predicate(lookup, pred, g), lookup, geom)
end
_predicate(::DE9IM.Intersects) = GO.intersects
_predicate(::DE9IM.Disjoint) = GO.disjoint
_predicate(::DE9IM.Contains) = GO.contains
_predicate(::DE9IM.Within) = GO.within
_predicate(::DE9IM.Covers) = GO.covers
_predicate(::DE9IM.CoveredBy) = GO.coveredby
_predicate(::DE9IM.Touches) = GO.touches
_predicate(::DE9IM.Crosses) = GO.crosses
_predicate(::DE9IM.Overlaps) = GO.overlaps
_predicate(::DE9IM.Equals) = GO.equals

# `hasselection` answers whether indexing would succeed: the selector must be one this
# lookup accepts and must select at least one geometry. `At`, `Contains` and `Near` need
# a method each to beat DimensionalData's per-selector fallbacks.
Lookups.hasselection(lookup::GeometryLookup, sel::Lookups.At) = _selects(lookup, sel)
Lookups.hasselection(lookup::GeometryLookup, sel::Lookups.Contains) = _selects(lookup, sel)
Lookups.hasselection(lookup::GeometryLookup, sel::Lookups.Near) = _selects(lookup, sel)
Lookups.hasselection(lookup::GeometryLookup, sel::Union{Lookups.Touches,Lookups.Where}) = _selects(lookup, sel)
Lookups.hasselection(lookup::GeometryLookup, sel::DE9IM.DE9IMPredicate) = _selects(lookup, sel)
# A DE9IM predicate is not a `Selector`, so it needs the step from dimension to lookup
# that DimensionalData makes for its own selectors.
Lookups.hasselection(dim::DD.Dimension{<:GeometryLookup}, sel::DE9IM.DE9IMPredicate) =
    Lookups.hasselection(val(dim), sel)

_selects(lookup::GeometryLookup, sel::Lookups.At) =
    _isgeometry(val(sel)) && !isnothing(_at(lookup, val(sel)))
_selects(lookup::GeometryLookup, sel::Lookups.Near) =
    _ispoint(val(sel)) && !isempty(parent(lookup)) &&
    (lookup.manifold isa GO.Planar || all(_ispoint, parent(lookup)))
_selects(lookup::GeometryLookup, sel::Lookups.At{<:GO.UnitSpherical.UnitSphericalPoint}) =
    !isnothing(_at(lookup, val(sel)))
_selects(lookup::GeometryLookup, sel::Lookups.Near{<:GO.UnitSpherical.UnitSphericalPoint}) =
    !isempty(parent(lookup)) && (lookup.manifold isa GO.Planar || all(_ispoint, parent(lookup)))
_selects(lookup::GeometryLookup, sel::Union{Lookups.At{<:AbstractVector},Lookups.Near{<:AbstractVector}}) =
    all(v -> _selects(lookup, DD.rebuild(sel; val=v)), val(sel))
_selects(lookup::GeometryLookup, sel::Lookups.Contains) =
    _isgeometry(val(sel)) && !isempty(Lookups.selectindices(lookup, sel))
_selects(lookup::GeometryLookup, sel::Lookups.Contains{<:GO.UnitSpherical.UnitSphericalPoint}) =
    !isempty(Lookups.selectindices(lookup, sel))
_selects(lookup::GeometryLookup, sel::Lookups.Contains{<:AbstractVector}) =
    all(_isgeometry, val(sel)) && !isempty(Lookups.selectindices(lookup, sel))
_selects(::GeometryLookup, ::Lookups.Touches) = false
_selects(lookup::GeometryLookup, sel::Lookups.Touches{<:Extents.Extent}) =
    !isempty(Lookups.selectindices(lookup, sel))
_selects(lookup::GeometryLookup, sel::Lookups.Where) = any(val(sel), parent(lookup))
# A curried predicate has a tree-narrowed `selectindices`, worth the allocation.
_selects(lookup::GeometryLookup, sel::Lookups.Where{<:Base.Fix2}) =
    !isempty(Lookups.selectindices(lookup, sel))
_selects(lookup::GeometryLookup, sel::DE9IM.DE9IMPredicate) =
    _wrapsgeometry(sel) && !isempty(Lookups.selectindices(lookup, sel))

# The shared "narrow by extent, refine by predicate" core. `pred(A, geom)` is tested
# for every candidate `A`.

function _select_predicate(lookup::GeometryLookup, pred, geom)
    geom = _querygeometry(lookup.manifold, geom)
    candidates = _maybe_get_candidates(lookup, geom)
    geoms = parent(lookup)
    return filter(i -> _predicate(lookup.manifold, pred, geoms[i], geom), candidates)
end
# Everything outside the candidate set is disjoint, so mark the intersecting candidates
# and keep the rest.
function _select_predicate(lookup::GeometryLookup, ::typeof(GO.disjoint), geom)
    geom = _querygeometry(lookup.manifold, geom)
    keep = trues(length(lookup))
    geoms = parent(lookup)
    for i in _maybe_get_candidates(lookup, geom)
        keep[i] = !_predicate(lookup.manifold, GO.intersects, geoms[i], geom)
    end
    return findall(keep)
end

_predicate(::GO.Planar, pred, a, b) = pred(a, b)
_predicate(m::GO.Spherical, pred, a, b) = pred(m, a, b)
_predicate(::GO.Spherical, ::typeof(GO.equals), a, b) = GO.equals(a, b)
_predicate(m::GO.Spherical, pred::Union{typeof(GO.crosses),typeof(GO.overlaps)}, a, b) =
    pred(GO.RelateNG(m), a, b)

_querygeometry(m, geom) = geom
function _querygeometry(::GO.Spherical, geom)
    return _normalizespherical(geom)
end
function _querygeometry(::GO.Spherical, box::Extents.Extent)
    xmin, xmax = box.X
    ymin, ymax = box.Y
    all(isfinite, (xmin, xmax, ymin, ymax)) && 0 < xmax - xmin < 180 &&
        -90 < ymin < ymax < 90 || throw(ArgumentError(
            "Spherical interval boxes need finite X/Y bounds, longitude width between 0 and 180 degrees, " *
            "and latitude bounds strictly between -90 and 90 degrees. Pass a geometry for other regions."
        ))
    return GI.Polygon([[(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]])
end

# Only extent relations that hold for a node whenever they hold for one of its children can
# narrow a tree query, which is why the narrowing is always `query`'s extent intersection:
# every predicate here needs at least one shared point. Candidates come back sorted.
function _maybe_get_candidates(lookup::GeometryLookup, geom)
    tree = spatialtree(lookup)
    isnothing(tree) && return 1:length(lookup)
    selector_extent = _indexextent(lookup.manifold, geom)
    isnothing(selector_extent) && return 1:length(lookup)
    Extents.disjoint(Extents.extent(tree), selector_extent) && return Int[]
    return GO.FlexibleRTrees.query(tree, selector_extent)
end

function _union(f, lookup::GeometryLookup, values)
    selected = falses(length(lookup))
    for v in values
        _mark!(selected, f(v))
    end
    return findall(selected)
end
_mark!(selected::AbstractVector{Bool}, i::Int) = (selected[i] = true; selected)
_mark!(selected::AbstractVector{Bool}, is) = (selected[is] .= true; selected)

function _first_covering(lookup::GeometryLookup, point)
    candidates = _maybe_get_candidates(lookup, point)
    geoms = parent(lookup)
    k = findfirst(i -> _predicate(lookup.manifold, GO.covers, geoms[i], point), candidates)
    return isnothing(k) ? nothing : candidates[k]
end

# An `Extents.Extent` satisfies `GI.isgeometry` (as a `RectangleTrait`), but it is a bounding
# box, not a geometry: `Touches(extent)` and the interval pairs take extents, the rest geometries.
_isgeometry(x) = GI.isgeometry(x) && !(GI.geomtrait(x) isa GI.RectangleTrait)
_ispoint(x) = _isgeometry(x) && GI.geomtrait(x) isa GI.AbstractPointTrait

function _checked_geometry(sel::Lookups.Selector)
    x = val(sel)
    _isgeometry(x) || throw(ArgumentError("""
        `$(nameof(typeof(sel)))` on a `GeometryLookup` wraps a GeoInterface geometry (a
        2-tuple of reals is a point); got `$sel`. For a selection per axis use
        `(X($(nameof(typeof(sel)))(x)), Y($(nameof(typeof(sel)))(y)))`.
        """))
    return x
end
function _checked_point(sel::Lookups.Near)
    x = val(sel)
    _ispoint(x) || throw(ArgumentError("""
        `Near` on a `GeometryLookup` wraps a point (a 2-tuple of reals or a GeoInterface
        point); got `$sel`.
        """))
    return x
end
function _checked_geometry(sel::DE9IM.DE9IMPredicate)
    _wrapsgeometry(sel) || throw(ArgumentError("""
        A DE9IM predicate on a `GeometryLookup` wraps a GeoInterface geometry or a vector
        of them, and cannot carry keyword arguments; got `$sel`. For an extent query use
        `Touches(extent)`.
        """))
    return parent(sel)
end
function _wrapsgeometry(sel::DE9IM.DE9IMPredicate)
    isempty(DE9IM.keywords(sel)) || return false
    geom = parent(sel)
    return _isgeometry(geom) || (geom isa AbstractVector && all(_isgeometry, geom))
end

# Nearest geometry to a point: branch and bound over the tree, visiting children by
# point-to-extent distance and pruning subtrees that cannot beat the best exact distance.
# Nodes at exactly that distance are still visited: ties go to the lowest index, like the scan.

function _nearest(lookup::GeometryLookup, point)
    point = _querygeometry(lookup.manifold, point)
    geoms = parent(lookup)
    isempty(geoms) && throw(ArgumentError("`Near` on an empty `GeometryLookup` has no nearest geometry."))
    # A non-finite coordinate prunes every node of the tree and leaves no nearest index.
    (isfinite(GI.x(point)) && isfinite(GI.y(point))) ||
        throw(ArgumentError("`Near` needs a point with finite coordinates; got `$point`."))
    if lookup.manifold isa GO.Spherical
        all(_ispoint, geoms) || throw(ArgumentError(
            "Spherical Near currently supports point lookups only; GeometryOps needs general spherical distance first."
        ))
        point = _unitpoint(point)
    end
    tree = spatialtree(lookup)
    isnothing(tree) && return last(findmin(g -> _nearest_distance(point, g), geoms))
    return first(_nearest(tree, point, geoms, 0, Inf))::Int
end
function _nearest(node, point, geoms, best_i, best_d)
    if STI.isleaf(node)
        for (i, ext) in STI.child_indices_extents(node)
            _extent_distance(point, ext) <= best_d || continue
            d = _nearest_distance(point, geoms[i])
            if d < best_d || (d == best_d && i < best_i)
                best_i, best_d = i, d
            end
        end
    else
        children = collect(STI.getchild(node))
        distances = map(c -> _extent_distance(point, STI.node_extent(c)), children)
        for k in sortperm(distances)
            distances[k] <= best_d || break
            best_i, best_d = _nearest(children[k], point, geoms, best_i, best_d)
        end
    end
    return best_i, best_d
end

_nearest_distance(point, geom) = GO.distance(point, geom)
_unitpoint(point) = GO.UnitSpherical.UnitSphericalPoint((Float64(GI.x(point)), Float64(GI.y(point))))
# Radius is common to the lookup, so angular distance preserves nearest ordering.
_nearest_distance(point::GO.UnitSpherical.UnitSphericalPoint, geom) =
    GO.UnitSpherical.spherical_distance(point, _unitpoint(geom))

function _extent_distance(point, ext::Extents.Extent)
    x, y = GI.x(point), GI.y(point)
    (xmin, xmax), (ymin, ymax) = ext.X, ext.Y
    dx = max(xmin - x, zero(x), x - xmax)
    dy = max(ymin - y, zero(y), y - ymax)
    return hypot(dx, dy)
end
function _extent_distance(point::GO.UnitSpherical.UnitSphericalPoint, ext::Extents.Extent)
    dx = max(ext.X[1] - point.x, 0.0, point.x - ext.X[2])
    dy = max(ext.Y[1] - point.y, 0.0, point.y - ext.Y[2])
    dz = max(ext.Z[1] - point.z, 0.0, point.z - ext.Z[2])
    # Unit XYZ boxes give a chord lower bound; round down before converting to angle.
    chord = max(0.0, hypot(dx, dy, dz) - 16eps(Float64))
    return 2asin(min(chord / 2, 1.0))
end

"""
    mask(lookup::GeometryLookup, sel) -> Vector{Bool}
    mask(A::Union{AbstractDimArray,AbstractDimStack}, sel) -> DimArray{Bool}

A boolean mask over the geometries of a lookup: `true` where `sel` selects. `sel` is any
selector `Lookups.selectindices` accepts on a [`GeometryLookup`](@ref): `Contains`, `At`,
`Near`, `Touches`, `Where`, a DE9IM predicate or a tuple of `X`/`Y` selectors.

On a cube, the mask is a `DimArray` over the dimension carrying the `GeometryLookup`;
a cube with no such dimension, or more than one, is an `ArgumentError` (pass the
lookup itself instead).

```julia
m = VectorDataCubes.mask(cube, Contains((0.5, 0.5)))
cube[Geometry = m]
```
"""
mask(lookup::GeometryLookup, sel) = _mark!(falses(length(lookup)), Lookups.selectindices(lookup, sel))
function mask(A::Union{DD.AbstractDimArray,DD.AbstractDimStack}, sel)
    geometry_dims = filter(d -> DD.lookup(d) isa GeometryLookup, DD.dims(A))
    length(geometry_dims) == 1 || throw(ArgumentError("""
        `mask` needs exactly one dimension with a `GeometryLookup`, but the object has
        $(length(geometry_dims)): `$(DD.name(geometry_dims))`. Pass the lookup directly:
        `mask(lookup(A, Geometry), sel)`.
        """))
    d = only(geometry_dims)
    return DD.DimArray(mask(DD.lookup(d), sel), (d,))
end
