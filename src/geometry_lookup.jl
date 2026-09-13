import DimensionalData as DD
import GeometryOps as GO, GeometryOpsCore as GOCore
import GeoInterface as GI
import Rasters as RA
import Extents
import Missings

using Rasters: isnokw, nokw, Lookups, val
using GeometryOps.FlexibleRTrees: RTree, BulkLoadAlgorithm, STR

"""
    Geometry

A dimension meant to be used with a [`GeometryLookup`](@ref).
"""
DD.@dim Geometry "Geometry"

# The tree indexes each geometry's `X`/`Y` extent as `Float64`, so its type is known
# from the geometry vector type alone, whatever the geometries' coordinate type or
# dimensionality.
const XYExtent = Extents.Extent{(:X, :Y),Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}
const XYRTree{A,D} = RTree{A,XYExtent,D,Vector{Int}}

function _xyextent(geometry)
    ext = GI.extent(geometry)
    return Extents.Extent(X=Float64.(ext.X), Y=Float64.(ext.Y))
end

# The lazily built spatial accelerator of a `GeometryLookup`. `spatialtree` builds the
# tree on its first call and caches it here; slicing or rebuilding a lookup with new
# geometries starts from a fresh, unbuilt index with the same algorithm. The slot is
# atomic because the build is double-checked: readers take no lock once a tree is there.
mutable struct SpatialIndex{A<:BulkLoadAlgorithm,Tr<:RTree}
    const algorithm::A
    @atomic tree::Union{Nothing,Tr}
    const lock::ReentrantLock
end
SpatialIndex(algorithm::BulkLoadAlgorithm, ::Type{D}) where {D<:AbstractVector} =
    SpatialIndex{typeof(algorithm),XYRTree{typeof(algorithm),D}}(algorithm, nothing, ReentrantLock())
SpatialIndex(tree::RTree) = SpatialIndex{typeof(tree.algorithm),typeof(tree)}(tree.algorithm, tree, ReentrantLock())

# `indices` as a `Vector` fixes the leaf index type for every algorithm (`Unsorted` would
# otherwise keep a `Base.OneTo`); the typed comprehension fixes the extent type for every
# geometry vector, where `map` would follow its eltype and its array type.
_buildtree(algorithm::BulkLoadAlgorithm, geometries::AbstractVector) =
    RTree(algorithm, geometries;
        indices=collect(eachindex(geometries)), extents=XYExtent[_xyextent(g) for g in geometries])

"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn, crs, tree, metadata)

A DimensionalData `Lookup` over geometries, with spatial indexing.

`GeometryLookup` is the lookup of the [`Geometry`](@ref) dimension of a
vector data cube. It holds a vector of geometries and a lazily built spatial
tree, so selectors such as `Contains(point)`, `Touches(extent)`,
`Where(GO.intersects(geom))` and the DE9IM.jl predicates resolve to indices with
a tree query followed by an exact `GeometryOps` predicate.

The lookup spans the two internal dimensions in `dims` as well as the
dimension it is wrapped in, so on a cube built as
`DimArray(values, Geometry(GeometryLookup(geoms)))` both
`cube[Geometry = Contains(point)]` and `cube[X(a..b), Y(c..d)]` work.

# Arguments

- `data`: a vector of geometries, or any table / collection that
  `GeometryOpsCore.get_geometries` understands. `missing` elements are an error.
- `dims`: one `X` and one `Y` dimension, in either order. Geometry coordinates
  are always `(x, y)`; `dims` only names the internal dimensions.

# Keywords

- `geometrycolumn`: the geometry column to read when `data` is a table.
- `crs`: the coordinate reference system. Defaults to `GeoInterface.crs` of
  `data`, then of its first geometry, then `nothing`.
- `tree`: the spatial accelerator, a `GeometryOps.FlexibleRTrees.RTree`. One of
  - not given: a sort-tile-recursive tree, built lazily on the first spatial query;
  - `nothing`: no accelerator, every query scans all geometries;
  - a bulk-load algorithm (`STR()`, `HPR()`, `Unsorted()`): built lazily with it;
  - a prebuilt `RTree` over the lookup's own geometry vector: stored as is.
- `metadata`: dimension metadata, `DimensionalData.NoMetadata()` by default.

# Examples

```julia
using Rasters

using NaturalEarth
import GeometryOps as GO

# construct the polygon lookup
polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
polygon_lookup = GeometryLookup(polygons, (X(), Y()))

# create a DimArray with the polygon lookup
dv = rand(Geometry(polygon_lookup))

# select the polygon with the centroid of the 88th polygon
only(dv[Geometry(Contains(GO.centroid(polygons[88])))]) == dv[88] # true
```
"""
struct GeometryLookup{T,A<:AbstractVector{T},D,M<:GO.Manifold,Tree<:Union{Nothing,SpatialIndex},CRS,Me} <: DD.Dimensions.MultiDimensionalLookup{T}
    manifold::M
    data::A
    tree::Tree
    dims::D
    crs::CRS
    metadata::Me
end
function GeometryLookup(
        data, dims=(DD.X(), DD.Y());
        geometrycolumn=nothing, crs=nokw, tree=nokw, metadata=Lookups.NoMetadata()
    )
    geometries = _checked_geometries(GOCore.get_geometries(data; geometrycolumn))
    if isnokw(crs)
        crs = GI.crs(data)
        if isnothing(crs) && !isempty(geometries)
            crs = GI.crs(first(geometries))
        end
    end
    return GeometryLookup(GO.Planar(), geometries, _spatialindex(tree, geometries), _checked_dims(dims), crs, metadata)
end

function _checked_geometries(geometries)
    if Missing <: eltype(geometries)
        any(ismissing, geometries) && _missing_geometries_error(geometries)
        geometries = Missings.disallowmissing(geometries)
    end
    all(GI.isgeometry, geometries) || _not_geometries_error(geometries)
    return geometries
end

@noinline function _missing_geometries_error(geometries)
    at = findall(ismissing, geometries)
    throw(ArgumentError("""
        `GeometryLookup` cannot hold `missing` geometries, but the collection has them
        at indices $(first(at, 5))$(length(at) > 5 ? ", …" : ""). Drop or fill those rows first.
        """))
end
@noinline function _not_geometries_error(geometries)
    at = findall(!GI.isgeometry, geometries)
    throw(ArgumentError("""
        Every element of a `GeometryLookup` must be a GeoInterface geometry
        (`GeoInterface.isgeometry(x) == true`), but the elements at indices
        $(first(at, 5))$(length(at) > 5 ? ", …" : "") are not.
        """))
end

function _checked_dims(dims)
    based = dims isa Tuple ? DD.basedims(dims) : dims
    ok = based isa Tuple && length(based) == 2 &&
        count(d -> d isa DD.XDim, based) == 1 && count(d -> d isa DD.YDim, based) == 1
    ok || throw(ArgumentError("""
        The `dims` of a `GeometryLookup` must be one `X` and one `Y` dimension,
        like `(X(), Y())` or `(Y(), X())`; got `$dims`.
        """))
    return based
end

_spatialindex(tree, geometries) = isnokw(tree) ? SpatialIndex(STR(), typeof(geometries)) : _spatialindex_error(tree)
_spatialindex(::Nothing, geometries) = nothing
_spatialindex(algorithm::BulkLoadAlgorithm, geometries) = SpatialIndex(algorithm, typeof(geometries))
function _spatialindex(tree::RTree, geometries)
    tree.data === geometries || throw(ArgumentError("""
        A prebuilt `tree` must index the lookup's own geometry vector — for a table, its
        geometry column — but `tree.data` is a different object. Build the tree over that
        vector, or pass a bulk-load algorithm (`STR()`, `HPR()`, `Unsorted()`) instead.
        """))
    return SpatialIndex(tree)
end
@noinline _spatialindex_error(tree) = throw(ArgumentError("""
    `tree` must be one of: not given (a lazily built `STR()` tree), `nothing` (no
    accelerator), a `GeometryOps.FlexibleRTrees` bulk-load algorithm (`STR()`, `HPR()`,
    `Unsorted()`), or a prebuilt `GeometryOps.FlexibleRTrees.RTree`; got a `$(typeof(tree))`.
    """))

_fresh(::Nothing, geometries) = nothing
_fresh(index::SpatialIndex, geometries) = SpatialIndex(index.algorithm, typeof(geometries))

"""
    spatialtree(l::GeometryLookup)

The spatial tree accelerating queries on `l`, or `nothing` when `l` was
constructed with `tree = nothing` or is empty.

The tree is built on the first call and cached in the lookup; slicing,
`view`, `reverse` and `rebuild` with new geometries never build one.
"""
spatialtree(l::GeometryLookup) = _spatialtree(l.tree, parent(l))
_spatialtree(::Nothing, geometries) = nothing
function _spatialtree(index::SpatialIndex, geometries)
    isempty(geometries) && return nothing
    tree = _builttree(index)
    isnothing(tree) || return tree
    return lock(index.lock) do
        built = _builttree(index)
        isnothing(built) || return built
        new_tree = _buildtree(index.algorithm, geometries)
        @atomic :release index.tree = new_tree
        return new_tree
    end
end

# The tree of an index that has already been built, without ever building one.
_builttree(l::GeometryLookup) = _builttree(l.tree)
_builttree(::Nothing) = nothing
_builttree(index::SpatialIndex) = @atomic :acquire index.tree

GI.crs(l::GeometryLookup) = l.crs
# Rasters reaches a lookup through `setcrs(dim::Dimension, crs)`, which passes the
# dimension it came from as a keyword.
RA.setcrs(l::GeometryLookup, crs; dim=nothing) = DD.rebuild(l; crs)

"""
    Rasters.reproject(target, l::GeometryLookup)

Reproject every geometry of `l` from its crs to `target`, returning a new lookup.

Needs Proj.jl: run `import Proj` first. A lookup without a crs cannot be
reprojected; set one with `Rasters.setcrs`.
"""
function RA.reproject(target::RA.GeoFormat, l::GeometryLookup)
    isnothing(GI.crs(l)) && throw(ArgumentError(
        "Cannot reproject a `GeometryLookup` with no crs. Set one first with `Rasters.setcrs`."
    ))
    return DD.rebuild(l; data=_reproject(target, l), crs=target)
end
# `ext/VectorDataCubesProjExt.jl` adds the method that does the work.
_reproject(target, ::GeometryLookup) = throw(ArgumentError(
    "Reprojecting a `GeometryLookup` needs Proj.jl: run `import Proj` and try again."
))

# DimensionalData interface

DD.dims(l::GeometryLookup) = l.dims
DD.dims(d::DD.Dimension{<:GeometryLookup}) = DD.dims(DD.val(d))
DD.order(::GeometryLookup) = Lookups.Unordered()
DD.parent(l::GeometryLookup) = l.data
Lookups.metadata(l::GeometryLookup) = l.metadata
# `format` rebuilds a lookup from its values, which cannot recover the internal dims
# or the spatial index; a `GeometryLookup` is complete as constructed.
DD.Dimensions.format(l::GeometryLookup, ::Type, values, axis::AbstractRange) = l

function DD.rebuild(
        l::GeometryLookup;
        data=l.data, tree=nokw, dims=l.dims, crs=nokw, manifold=l.manifold, metadata=l.metadata
    )
    index = if isnokw(tree)
        data === l.data ? l.tree : _fresh(l.tree, data)
    else
        _spatialindex(tree, data)
    end
    new_crs = if isnokw(crs)
        data_crs = GI.crs(data)
        isnothing(data_crs) ? l.crs : data_crs
    else
        crs
    end
    new_dims = dims === l.dims ? dims : _checked_dims(dims)
    return GeometryLookup(manifold, data, index, new_dims, new_crs, metadata)
end

Base.reverse(l::GeometryLookup) = DD.rebuild(l; data=reverse(parent(l)))

Lookups._set_lookup(::Lookups.Safety, ::Lookups.Lookup, new::GeometryLookup) = new
# `set(cube, Geometry => geometries)` replaces the lookup values; check them as the
# constructor does, so a lookup can never come to hold `missing` or a non-geometry.
Lookups._set_lookup_parent(::Lookups.Safe, l::GeometryLookup, values::AbstractVector) =
    DD.rebuild(l; data=_checked_geometries(values))
Lookups._set_lookup_parent(::Lookups.Safe, l::GeometryLookup, ::Lookups.AutoValues) = l

function Lookups.bounds(l::GeometryLookup)
    ext = _extent(l)
    return map(d -> _dimbounds(ext, d), DD.dims(l))
end
_dimbounds(::Nothing, ::DD.Dimension) = (nothing, nothing)
_dimbounds(ext::Extents.Extent, ::DD.XDim) = ext.X
_dimbounds(ext::Extents.Extent, ::DD.YDim) = ext.Y

function _extent(l::GeometryLookup)
    geometries = parent(l)
    isempty(geometries) && return nothing
    tree = _builttree(l)
    isnothing(tree) || return Extents.extent(tree)
    # `_xyextent`, like the tree's own extent, so bounds do not change once it is built.
    return mapreduce(_xyextent, Extents.union, geometries)
end

function Base.:(==)(a::GeometryLookup, b::GeometryLookup)
    a === b && return true
    return DD.name(DD.dims(a)) == DD.name(DD.dims(b)) && GI.crs(a) == GI.crs(b) &&
        (parent(a) === parent(b) || parent(a) == parent(b))
end
function Base.isequal(a::GeometryLookup, b::GeometryLookup)
    a === b && return true
    return isequal(DD.name(DD.dims(a)), DD.name(DD.dims(b))) && isequal(GI.crs(a), GI.crs(b)) &&
        (parent(a) === parent(b) || isequal(parent(a), parent(b)))
end
Base.hash(l::GeometryLookup, h::UInt) =
    hash(parent(l), hash(GI.crs(l), hash(DD.name(DD.dims(l)), hash(:GeometryLookup, h))))

@inline Lookups.reducelookup(::GeometryLookup) = Lookups.NoLookup(Base.OneTo(1))

function Lookups.show_compact(io::IO, mime, l::GeometryLookup)
    print(io, "GeometryLookup{", _elname(eltype(l)), "}")
end
# Sorted, so the header does not depend on the order Julia happens to store a union in.
_elname(T::Union) = string("Union{", join(sort!(map(string ∘ _elname, Base.uniontypes(T))), ", "), "}")
_elname(T::Type) = nameof(T)

function Lookups.show_properties(io::IO, mime, l::GeometryLookup)
    print(io, " ")
    show(IOContext(io, :inset => "", :dimcolor => 244), mime, DD.basedims(l))
end
