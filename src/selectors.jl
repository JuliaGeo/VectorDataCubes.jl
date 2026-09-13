# Return an `Int` or Vector{Bool}
# Base case: got a standard index that can go into getindex on a base Array
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel
# other cases: 
# - decompose selectors
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple)
    DD.selectindices(lookup, map(_val_or_nothing, DD.sortdims(sel, DD.dims(lookup))))
end
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    dimsel = map(DD.rebuild, map(DD.name2dim, K), DD.values(sel))
    DD.selectindices(lookup, dimsel) 
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    if (length(sel) == length(DD.dims(lookup))) && all(map(s -> s isa At, sel))
        i = findfirst(x -> all(map(DD.Dimensions._matches, sel, x)), lookup)
        isnothing(i) && _coord_not_found_error(sel)
        return i
    else
        return [DD.Dimensions._matches(sel, x) for x in lookup]
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Contains)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    filter(potential_candidates) do candidate
        GO.contains(lookup.data[candidate], val(sel))
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.At)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    candidates = _maybe_get_candidates(lookup, GI.extent(geom))
    x = findfirst(candidates) do candidate
        GO.equals(geom, lookup.data[candidate])
    end
    if isnothing(x)
        throw(ArgumentError("$sel not found in lookup"))
    else
        return candidates[x]
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Near)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    # TODO: temporary
    @assert GI.trait(geom) isa GI.PointTrait "Only point geometries are supported for the near lookup at this point!  We will add more geometry support in the future."

    # Get the nearest geometry
    # TODO: this sucks!  Use some branch and bound algorithm
    # on the spatial tree instead.
    # if pointtrait
    return findmin(x -> GO.distance(geom, x), lookup.data)[2]
    # else
    #     findmin(x -> GO.distance(GO.GEOS(), geom, x), lookup.data)[2]
    # end 
    # this depends on LibGEOS being installed.

end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Touches)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], val(sel))
    end
end
function Lookups.selectindices(
    lookup::GeometryLookup, 
    (xs, ys)::Tuple{Union{<:DD.Lookups.Touches}, Union{<:DD.Lookups.Touches}}
)
    target_ext = Extents.Extent(X = (first(xs), last(xs)), Y = (first(ys), last(ys)))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], target_ext)
    end
end
function Lookups.selectindices(
    lookup::GeometryLookup, 
    (xs, ys)::Tuple{Union{<:DD.IntervalSets.ClosedInterval},Union{<:DD.IntervalSets.ClosedInterval}}
)
    target_ext = Extents.Extent(X = extrema(xs), Y = extrema(ys))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    filter(potential_candidates) do candidate
        GO.covers(target_ext, lookup.data[candidate])
    end
end
function Lookups.selectindices(
    lookup::GeometryLookup, 
    (x, y)::Tuple{Union{<:DD.Lookups.At,<:DD.Lookups.Contains}, Union{<:DD.Lookups.At,<:DD.Lookups.Contains}}
)
    xval, yval = val(x), val(y)
    # `At` requires an exact match, so a point matching no geometry is an error;
    # `Contains` just returns the (possibly empty) set of matches.
    is_at = x isa DD.Lookups.At && y isa DD.Lookups.At
    potential_candidates = if isnothing(lookup.tree)
        collect(eachindex(lookup.data))
    else
        # The tree query already rejects points outside the root node's
        # extent, so no separate extent check is needed here.
        GO.SpatialTreeInterface.query(lookup.tree, (xval, yval))
    end
    if is_at
        # If both selectors are At(), return a single index (first match) for speed and clarity
        for candidate in potential_candidates
            if GO.contains(lookup.data[candidate], (xval, yval))
                return candidate
            end
        end
        throw(ArgumentError("Point ($xval, $yval) not found in lookup"))
    else
        # For Contains selectors, return all matching indices
        filter(potential_candidates) do candidate
            GO.contains(lookup.data[candidate], (xval, yval))
        end
    end
end

for fname in (:equals, :intersects, 
                :contains, :within, :covers, 
                :coveredby, :touches)
    @eval begin
        function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Where{Base.Fix2{typeof(GO.$fname)}})
            sel_ext = GI.extent(val(sel).x)
            potential_candidates = _maybe_get_candidates(lookup, sel_ext)
            f = val(sel)
            return filter(potential_candidates) do idx
                f(lookup.data[idx])
            end
        end
    end
end
# Disjoint needs a specialized implementation, which will look at intersects instead.
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Where{Base.Fix2{typeof(GO.disjoint)}})
    sel_ext = GI.extent(val(sel).x)
    potential_candidates = _maybe_get_candidates(lookup, sel_ext, Extents.intersects)
    f = val(sel)
    actual_intersections = filter(potential_candidates) do idx
        f(lookup.data[idx])
    end
    return setdiff(1:length(lookup.data), actual_intersections)
end
# Local functions
_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)

# Get the candidates for the selector extent.  
# If the selector extent is disjoint from the tree rootnode extent,
# you should raise an error.  We should have an error type that can be
# plotted etc. to allow debugging and understanding.
function _maybe_get_candidates(lookup::GeometryLookup, selector_extent, operation::O) where O
    tree = lookup.tree
    isnothing(tree) && return 1:length(lookup)
    Extents.disjoint(GI.extent(tree), selector_extent) && return Int[]
    potential_candidates = GO.SpatialTreeInterface.query(
        tree,
        Base.Fix1(operation, selector_extent)
    )
    isempty(potential_candidates) && return Int[]
    return potential_candidates
end

_maybe_get_candidates(lookup::GeometryLookup, selector_extent) = _maybe_get_candidates(lookup, selector_extent, Extents.intersects)