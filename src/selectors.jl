import DE9IM

const STI = GO.SpatialTreeInterface

# Every selector on a `GeometryLookup` resolves here: narrow the candidates with an
# extent query on `spatialtree(lookup)` (all indices when there is no tree), then
# refine with an exact GeometryOps predicate. Multi-index results are `Vector{Int}`,
# single-match selectors (`At`, `Near`, the `(At, At)` pair) return an `Int`.

Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel

# Dimension-wrapped selectors: sort into the lookup's internal dim order, then pair up.
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple)
    Lookups.selectindices(lookup, map(_val_or_nothing, DD.sortdims(sel, DD.dims(lookup))))
end
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    Lookups.selectindices(lookup, map(DD.rebuild, map(DD.name2dim, K), values(sel)))
end
_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)

"""
    Lookups.selectindices(lookup::GeometryLookup, (xsel, ysel)::Tuple)

Select geometries by a pair of selectors on the internal `X` and `Y` dimensions,
written on a cube as `cube[X(xsel), Y(ysel)]` or `cube[Geometry = (X(xsel), Y(ysel))]`.
Each selector is matched to a coordinate by its dimension type, so the pair works the
same on a lookup with `dims = (Y(), X())`.

- `(At(x), At(y))`: the geometry covering the point, the lowest matching index as an
  `Int`; none is an `ArgumentError`.
- `(Contains(x), Contains(y))`: every geometry covering the point (closed: a point on
  a shared border belongs to each geometry it lies on), a `Vector{Int}`.
- `(Near(x), Near(y))`: the geometry nearest to the point, an `Int`.
- `(a .. b, c .. d)`: geometries fully covered by the box.
- `(Touches(a, b), Touches(c, d))`: geometries intersecting the box.
- `X(a .. b)` or `X(Touches(a, b))` alone: the other axis is unbounded.

Any other combination is an `ArgumentError`.
"""
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    length(sel) == 2 || _pair_error(sel)
    return _select_pair(lookup, _xy(lookup, sel)...)
end
function _xy(lookup::GeometryLookup, sel::Tuple)
    first(DD.dims(lookup)) isa DD.XDim ? (sel[1], sel[2]) : (sel[2], sel[1])
end
_pair_error(sel) = throw(ArgumentError("""
    Unsupported selector pair `$sel` on a `GeometryLookup`. The supported pairs are
    `(X(At(x)), Y(At(y)))`, `(X(Contains(x)), Y(Contains(y)))`, `(X(Near(x)), Y(Near(y)))`,
    `(X(a .. b), Y(c .. d))` and `(X(Touches(a, b)), Y(Touches(c, d)))`; an interval or
    `Touches` may also be given on one axis alone.
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
    return _select_predicate(lookup, _covered_by_box, box)
end
_covered_by_box(geom, box) = GO.covers(box, geom)
function _select_pair(lookup::GeometryLookup, x::_PairTouches, y::_PairTouches)
    box = Extents.Extent(X=_axisbounds(x), Y=_axisbounds(y))
    return _select_predicate(lookup, GO.intersects, box)
end
_axisbounds(::Nothing) = (-Inf, Inf)
_axisbounds(i::DD.IntervalSets.Interval) = extrema(i)
_axisbounds(t::Lookups.Touches) = extrema(val(t))

# Single selectors on the `Geometry` axis

"""
    Lookups.selectindices(lookup::GeometryLookup, sel::Contains)

The indices of the geometries covering the point (or geometry) `Contains` wraps,
as a `Vector{Int}`: `cube[Geometry = Contains((x, y))]`.

Point-in-polygon is closed (`GeometryOps.covers`): a point on a shared border belongs
to every geometry whose boundary it lies on. A vector of points selects the sorted
union of the matches for each.

The value must be a GeoInterface geometry; a 2-tuple of reals is a point.
"""
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Contains) =
    _select_predicate(lookup, GO.covers, _checked_geometry(sel))
function Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Contains{<:AbstractVector})
    return _union(lookup, val(sel)) do v
        Lookups.selectindices(lookup, DD.rebuild(sel; val=v))
    end
end

"""
    Lookups.selectindices(lookup::GeometryLookup, sel::At)

The index of the geometry equal (`GeometryOps.equals`) to the geometry `At` wraps;
no match is an `ArgumentError`. Several equal geometries give the lowest index. A
vector of geometries gives a `Vector{Int}`, one index per value.
"""
function Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.At)
    geom = _checked_geometry(sel)
    i = _at(lookup, geom)
    isnothing(i) && throw(ArgumentError("No geometry in the lookup equals `$(sel)`."))
    return i
end
function _at(lookup::GeometryLookup, geom)
    candidates = _maybe_get_candidates(lookup, GI.extent(geom))
    geoms = parent(lookup)
    k = findfirst(i -> GO.equals(geoms[i], geom), candidates)
    return isnothing(k) ? nothing : candidates[k]
end

"""
    Lookups.selectindices(lookup::GeometryLookup, sel::Near)

The index of the geometry nearest (`GeometryOps.distance`) to the point `Near` wraps;
geometries containing the point are at distance zero, and ties go to the lowest index.
The search is a branch-and-bound over [`spatialtree`](@ref), or a linear scan when
the lookup has no tree. Only points are supported; an empty lookup is an
`ArgumentError`. A vector of points gives a `Vector{Int}`, one index per point.
"""
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Near) = _nearest(lookup, _checked_point(sel))
Lookups.selectindices(lookup::GeometryLookup, sel::Union{Lookups.At{<:AbstractVector},Lookups.Near{<:AbstractVector}}) =
    Int[Lookups.selectindices(lookup, DD.rebuild(sel; val=v)) for v in val(sel)]

"""
    Lookups.selectindices(lookup::GeometryLookup, sel::Touches{<:Extents.Extent})

The indices of the geometries intersecting (`GeometryOps.intersects`) the
`Extents.Extent` `Touches` wraps, as a `Vector{Int}`:
`cube[Geometry = Touches(GI.extent(geom))]`. This is the loose "touches" of
DimensionalData's `Touches`; the strict DE-9IM relation (boundaries meet, interiors do
not) is `DE9IM.Touches(geom)`.

`Touches` accepts an extent or a pair of bounds only — DimensionalData's type
constrains its value — so intersection with a geometry is `Where(GO.intersects(geom))`
or `DE9IM.Intersects(geom)`, and the bounds form is per axis:
`cube[X(Touches(a, b)), Y(Touches(c, d))]`.
"""
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Touches{<:Extents.Extent}) =
    _select_predicate(lookup, GO.intersects, val(sel))
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Touches) = throw(ArgumentError("""
    `Touches` on a `GeometryLookup` wraps an `Extents.Extent`; got `$sel`. For the
    bounds form use it per axis: `(X(Touches(a, b)), Y(Touches(c, d)))`.
    """))

"""
    Lookups.selectindices(lookup::GeometryLookup, sel::Where)

The indices `i` for which `f(geometry_i)` holds, as a `Vector{Int}`.

A curried GeometryOps predicate — `Where(GO.intersects(geom))`, and likewise
`equals`, `contains`, `within`, `covers`, `coveredby`, `touches`, `crosses`,
`overlaps`, `disjoint` — is narrowed with the spatial tree first; any other function
is a linear `findall`.
"""
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Where) = findall(val(sel), parent(lookup))

for pred in (:equals, :intersects, :contains, :within, :covers, :coveredby, :touches,
             :crosses, :overlaps, :disjoint)
    @eval Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.Where{<:Base.Fix2{typeof(GO.$pred)}}) =
        _select_predicate(lookup, GO.$pred, val(sel).x)
end

"""
    Lookups.selectindices(lookup::GeometryLookup, sel::DE9IM.DE9IMPredicate)

The indices of the geometries `A` for which `pred(A, geom)` holds, for a DE-9IM
predicate `pred(geom)` from DE9IM.jl: `Intersects`, `Disjoint`, `Contains`, `Within`,
`Covers`, `CoveredBy`, `Touches`, `Crosses`, `Overlaps`, `Equals`. The wrapped geometry
is the second argument, as in `cube[Geometry = DE9IM.Covers(geom)]`.

A vector of geometries selects the sorted union of the matches for each. Keyword
arguments on the predicate are not supported.
"""
function Lookups.selectindices(lookup::GeometryLookup, sel::DE9IM.DE9IMPredicate)
    geom = _checked_geometry(sel)
    pred = _predicate(sel)
    geom isa AbstractVector || return _select_predicate(lookup, pred, geom)
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
    _ispoint(val(sel)) && !isempty(parent(lookup))
_selects(lookup::GeometryLookup, sel::Union{Lookups.At{<:AbstractVector},Lookups.Near{<:AbstractVector}}) =
    all(v -> _selects(lookup, DD.rebuild(sel; val=v)), val(sel))
_selects(lookup::GeometryLookup, sel::Lookups.Contains) =
    _isgeometry(val(sel)) && !isempty(Lookups.selectindices(lookup, sel))
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
    candidates = _maybe_get_candidates(lookup, GI.extent(geom))
    geoms = parent(lookup)
    return filter(i -> pred(geoms[i], geom), candidates)
end
# Everything outside the candidate set is disjoint, so mark the intersecting candidates
# and keep the rest.
function _select_predicate(lookup::GeometryLookup, ::typeof(GO.disjoint), geom)
    keep = trues(length(lookup))
    geoms = parent(lookup)
    for i in _maybe_get_candidates(lookup, GI.extent(geom))
        keep[i] = !GO.intersects(geoms[i], geom)
    end
    return findall(keep)
end

# Only extent relations that hold for a node whenever they hold for one of its children
# can narrow a tree query, which is why this is always `Extents.intersects`: every
# predicate here needs at least one shared point.
function _maybe_get_candidates(lookup::GeometryLookup, selector_extent)
    tree = spatialtree(lookup)
    (isnothing(tree) || isnothing(selector_extent)) && return 1:length(lookup)
    Extents.disjoint(GI.extent(tree), selector_extent) && return Int[]
    return STI.query(tree, Base.Fix1(Extents.intersects, selector_extent))::Vector{Int}
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
    x, y = point
    candidates = _maybe_get_candidates(lookup, Extents.Extent(X=(x, x), Y=(y, y)))
    geoms = parent(lookup)
    k = findfirst(i -> GO.covers(geoms[i], point), candidates)
    return isnothing(k) ? nothing : candidates[k]
end

# An `Extents.Extent` satisfies `GI.isgeometry`, but it is a bounding box rather than a
# geometry: `Touches(extent)` and the interval pairs take extents, the rest take geometries.
_isgeometry(x) = GI.isgeometry(x) && !(x isa Extents.Extent)
_ispoint(x) = _isgeometry(x) && GI.trait(x) isa GI.PointTrait

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
    return geom isa AbstractVector ? all(_isgeometry, geom) : _isgeometry(geom)
end

# Nearest geometry to a point: branch and bound over the spatial tree, visiting children
# in order of point-to-extent distance and pruning subtrees that cannot beat the best
# exact distance so far. Subtrees at exactly that distance are still visited, so ties
# resolve to the lowest index, as in the linear scan.

function _nearest(lookup::GeometryLookup, point)
    geoms = parent(lookup)
    isempty(geoms) && throw(ArgumentError("`Near` on an empty `GeometryLookup` has no nearest geometry."))
    # A non-finite coordinate prunes every node of the tree and leaves no nearest index.
    (isfinite(GI.x(point)) && isfinite(GI.y(point))) ||
        throw(ArgumentError("`Near` needs a point with finite coordinates; got `$point`."))
    tree = spatialtree(lookup)
    isnothing(tree) && return last(findmin(g -> GO.distance(point, g), geoms))
    return first(_nearest(tree, point, geoms, 0, Inf))::Int
end
function _nearest(node, point, geoms, best_i, best_d)
    if STI.isleaf(node)
        for (i, ext) in STI.child_indices_extents(node)
            _extent_distance(point, ext) <= best_d || continue
            d = GO.distance(point, geoms[i])
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
function _extent_distance(point, ext::Extents.Extent)
    x, y = GI.x(point), GI.y(point)
    (xmin, xmax), (ymin, ymax) = ext.X, ext.Y
    dx = max(xmin - x, zero(x), x - xmax)
    dy = max(ymin - y, zero(y), y - ymax)
    return hypot(dx, dy)
end

"""
    mask(lookup::GeometryLookup, sel) -> Vector{Bool}
    mask(A::Union{AbstractDimArray,AbstractDimStack}, sel) -> DimArray{Bool}

A boolean mask over the geometries of a lookup: `true` where the selector `sel`
selects. `sel` is any selector `Lookups.selectindices` accepts on a
[`GeometryLookup`](@ref) — `Contains`, `At`, `Near`, `Touches`, `Where`, a DE9IM
predicate or a tuple of `X`/`Y` selectors.

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
