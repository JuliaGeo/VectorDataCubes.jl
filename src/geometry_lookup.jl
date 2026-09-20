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

# Fixed extent types keep lazy lookups concrete even when empty or sliced.
const XYExtent = Extents.Extent{(:X, :Y),Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}
const XYZExtent = Extents.Extent{(:X, :Y, :Z),NTuple{3,Tuple{Float64,Float64}}}
_extenttype(::GO.Planar) = XYExtent
_extenttype(::GO.Spherical) = XYZExtent

function _xyextent(geometry)
    ext = GI.extent(geometry)
    return Extents.Extent(X=Float64.(ext.X), Y=Float64.(ext.Y))
end

_indexextent(::GO.Planar, geometry) = _xyextent(geometry)
function _indexextent(m::GO.Spherical, geometry)
    # GO 0.1.47's spherical extent conversion does not promote Float32 input.
    geom = GO.apply(GI.PointTrait(), geometry) do p
        (Float64(GI.x(p)), Float64(GI.y(p)))
    end
    return GO.extent(m, geom)
end

# The lazily built spatial accelerator of a `GeometryLookup`. `spatialtree` builds and caches
# the tree on first use; a rebuild with new geometries gets a fresh index, same algorithm. The
# slot is atomic so readers take no lock once a tree is there (the build is double-checked).
mutable struct SpatialIndex{A<:BulkLoadAlgorithm,Tr<:RTree}
    const algorithm::A
    @atomic tree::Union{Nothing,Tr}
    const lock::ReentrantLock
end
function SpatialIndex(m, algorithm::BulkLoadAlgorithm, ::Type{D}) where {D<:AbstractVector}
    # RTree's type parameter layout can change independently of its constructor API.
    Tr = Base.promote_op(_buildtree, typeof(m), typeof(algorithm), D)
    return SpatialIndex{typeof(algorithm),Tr}(algorithm, nothing, ReentrantLock())
end
SpatialIndex(tree::RTree) = SpatialIndex{typeof(tree.algorithm),typeof(tree)}(tree.algorithm, tree, ReentrantLock())

# `indices` as a `Vector` fixes the leaf index type for every algorithm (`Unsorted` would
# otherwise keep a `Base.OneTo`); the typed comprehension fixes the extent type for every
# geometry vector, where `map` would follow its eltype and its array type.
function _buildtree(m, algorithm::BulkLoadAlgorithm, geometries::AbstractVector)
    E = _extenttype(m)
    return RTree(algorithm, geometries;
        indices=collect(eachindex(geometries)), extents=E[_indexextent(m, g) for g in geometries])
end

"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn, crs, manifold, tree, metadata)

A DimensionalData `Lookup` over geometries, with spatial indexing.

The lookup of the [`Geometry`](@ref) dimension of a vector data cube: a vector of geometries
plus a lazily built spatial tree. Selectors — `Contains(point)`, `Touches(extent)`,
`Where(GO.intersects(geom))`, DE9IM.jl predicates — narrow by the tree, then test exactly.

The lookup spans its internal `dims` as well as the dimension wrapping it, so on
`DimArray(values, Geometry(GeometryLookup(geoms)))` both `cube[Geometry = Contains(point)]`
and `cube[X(a..b), Y(c..d)]` work.

# Arguments

- `data`: a vector of geometries, or any table / collection that
  `GeometryOpsCore.get_geometries` understands. `missing` elements are an error.
- `dims`: one `X` and one `Y` dimension, in either order. Geometry coordinates
  are always `(x, y)`; `dims` only names the internal dimensions.

# Keywords

- `geometrycolumn`: the geometry column to read when `data` is a table.
- `crs`: the coordinate reference system. Defaults to `GeoInterface.crs` of
  `data`, then of its first geometry, then `nothing`.
- `manifold`: `GeometryOps.Planar()` or `GeometryOps.Spherical()`. Defaults to
  the selected table column's `edges`/`orientation` metadata, otherwise planar.
  Spherical geometries use longitude/latitude in degrees and great-circle edges.
  `UnitSphericalPoint` coordinates are converted to longitude/latitude at construction.
  Finite interval boxes become great-circle polygons; `Near` currently supports
  spherical point lookups only. `Spherical(oriented=true)` asserts that ring
  interiors lie to the left of the stored vertex order. The CRS is preserved
  independently. With Proj loaded, spherical table input derives and validates the
  physical radius from the CRS datum.
- `tree`: the spatial accelerator, a `GeometryOps.FlexibleRTrees.RTree`. One of
  - not given: a sort-tile-recursive tree, built lazily on the first spatial query;
  - `nothing`: no accelerator, every query scans all geometries;
  - a bulk-load algorithm (`STR()`, `HPR()`, `Unsorted()`): built lazily with it;
  - a prebuilt `RTree` over the lookup's own geometry vector: stored as is for
    planar lookups. Spherical lookups currently require an algorithm or `nothing`.
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
        geometrycolumn=nothing, crs=nokw, manifold=nokw, tree=nokw, metadata=Lookups.NoMetadata()
    )
    geometries = _checked_geometries(GOCore.get_geometries(data; geometrycolumn))
    infer_manifold = isnokw(manifold)
    infer_manifold && (manifold = _inputmanifold(data, geometrycolumn))
    _checkmanifold(manifold)
    isnokw(crs) && (crs = _inputcrs(data, geometries, manifold))
    geometries, normalized = _normalizecoordinates(manifold, geometries)
    if normalized && crs == _unitsphericalcrs()
        throw(ArgumentError(
            "UnitSphericalPoint input was converted to longitude/latitude, but its CRS is the " *
            "internal unit-sphere Cartesian CRS. Pass the geographic CRS explicitly with `crs`."
        ))
    end
    infer_manifold && (manifold = _inputmanifold(data, geometrycolumn, crs))
    return GeometryLookup(manifold, geometries, _spatialindex(manifold, tree, geometries), _checked_dims(dims), crs, metadata)
end

function _inputcrs(data, geometries, manifold)
    crs = GI.isgeometry(data) ? _geometrycrs(manifold, data) : GI.crs(data)
    if isnothing(crs) && !isempty(geometries)
        crs = _geometrycrs(manifold, first(geometries))
    end
    return crs
end
_geometrycrs(::GO.Planar, geometry) = GI.crs(geometry)
function _geometrycrs(::GO.Spherical, geometry)
    crs = GI.crs(geometry)
    # USP's Cartesian CRS describes its storage, not the normalized public geometry.
    return _hasusp(geometry) && crs == _unitsphericalcrs() ? nothing : crs
end

_checkmanifold(::GO.Planar) = nothing
function _checkmanifold(m::GO.Spherical)
    isfinite(m.radius) && m.radius > 0 || throw(ArgumentError("A spherical radius must be finite and positive."))
end
_checkmanifold(m) = throw(ArgumentError("GeometryLookup supports Planar() or Spherical(); got $m."))

_normalizecoordinates(::GO.Planar, geometries) = (geometries, false)
function _normalizecoordinates(::GO.Spherical, geometries)
    any(_hasusp, geometries) || return (geometries, false)
    return map(_normalizespherical, geometries), true
end

_hasusp(geometry) = GO.applyreduce(|, GI.PointTrait(), geometry; init=false) do point
    point isa GO.UnitSpherical.UnitSphericalPoint
end

_unitsphericalcrs() = GI.crs(GO.UnitSpherical.UnitSphericalPoint((0.0, 0.0)))

function _normalizespherical(geometry)
    _hasusp(geometry) || return geometry
    inverse = GO.UnitSpherical.GeographicFromUnitSphere()
    return GO.apply(GI.PointTrait(), geometry; crs=nothing) do point
        point isa GO.UnitSpherical.UnitSphericalPoint ? inverse(point) : point
    end
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

_spatialindex(m, tree, geometries) = isnokw(tree) ? SpatialIndex(m, STR(), typeof(geometries)) : _spatialindex_error(tree)
_spatialindex(m, ::Nothing, geometries) = nothing
_spatialindex(m, algorithm::BulkLoadAlgorithm, geometries) = SpatialIndex(m, algorithm, typeof(geometries))
function _spatialindex(m, tree::RTree, geometries)
    tree.data === geometries || throw(ArgumentError("""
        A prebuilt `tree` must index the lookup's own geometry vector — for a table, its
        geometry column — but `tree.data` is a different object. Build the tree over that
        vector, or pass a bulk-load algorithm (`STR()`, `HPR()`, `Unsorted()`) instead.
    """))
    # Released RTree has no manifold field to validate, including ring orientation.
    m isa GO.Planar || throw(ArgumentError(
        "Prebuilt trees are supported for planar lookups only; pass a bulk-load algorithm for a spherical lookup."
    ))
    return SpatialIndex(tree)
end
@noinline _spatialindex_error(tree) = throw(ArgumentError("""
    `tree` must be one of: not given (a lazily built `STR()` tree), `nothing` (no
    accelerator), a `GeometryOps.FlexibleRTrees` bulk-load algorithm (`STR()`, `HPR()`,
    `Unsorted()`), or a prebuilt `GeometryOps.FlexibleRTrees.RTree`; got a `$(typeof(tree))`.
    """))

_fresh(m, ::Nothing, geometries) = nothing
_fresh(m, index::SpatialIndex, geometries) = SpatialIndex(m, index.algorithm, typeof(geometries))

"""
    spatialtree(l::GeometryLookup)

The spatial tree accelerating queries on `l`, or `nothing` when `l` was
constructed with `tree = nothing` or is empty.

The tree is built on the first call and cached in the lookup; slicing,
`view`, `reverse` and `rebuild` with new geometries never build one.
"""
spatialtree(l::GeometryLookup) = _spatialtree(l.manifold, l.tree, parent(l))
_spatialtree(m, ::Nothing, geometries) = nothing
function _spatialtree(m, index::SpatialIndex, geometries)
    isempty(geometries) && return nothing
    tree = _builttree(index)
    isnothing(tree) || return tree
    return lock(index.lock) do
        built = _builttree(index)
        isnothing(built) || return built
        new_tree = _buildtree(m, index.algorithm, geometries)
        @atomic :release index.tree = new_tree
        return new_tree
    end
end

# The tree of an index that has already been built, without ever building one.
_builttree(l::GeometryLookup) = _builttree(l.tree)
_builttree(::Nothing) = nothing
_builttree(index::SpatialIndex) = @atomic :acquire index.tree

GI.crs(l::GeometryLookup) = l.crs
GOCore.manifold(l::GeometryLookup) = l.manifold
# Rasters reaches a lookup through `setcrs(dim::Dimension, crs)`, which passes the
# dimension it came from as a keyword.
function RA.setcrs(l::GeometryLookup, crs; dim=nothing)
    _validate_manifold_crs(l.manifold, crs)
    return DD.rebuild(l; crs)
end

# Needs Proj.jl loaded, like `GeometryOps.reproject` itself.
function RA.reproject(target::RA.GeoFormat, l::GeometryLookup)
    l.manifold isa GO.Planar || throw(ArgumentError(
        "Reprojecting a spherical GeometryLookup requires an explicit edge-conversion policy and is not supported yet."
    ))
    isnothing(GI.crs(l)) && throw(ArgumentError(
        "Cannot reproject a `GeometryLookup` with no crs. Set one first with `Rasters.setcrs`."
    ))
    geometries = GO.reproject(parent(l); source_crs=GI.crs(l), target_crs=target)
    return DD.rebuild(l; data=geometries, crs=target)
end

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
    _checkmanifold(manifold)
    if data !== l.data || manifold != l.manifold
        data = first(_normalizecoordinates(manifold, _checked_geometries(data)))
    end
    index = if isnokw(tree)
        data === l.data && manifold == l.manifold ? l.tree : _fresh(manifold, l.tree, data)
    else
        _spatialindex(manifold, tree, data)
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
    l.manifold isa GO.Planar && !isnothing(tree) && return Extents.extent(tree)
    # `_xyextent`, like the tree's own extent, so bounds do not change once it is built.
    return mapreduce(_xyextent, Extents.union, geometries)
end

function Base.:(==)(a::GeometryLookup, b::GeometryLookup)
    a === b && return true
    return a.manifold == b.manifold && DD.name(DD.dims(a)) == DD.name(DD.dims(b)) && GI.crs(a) == GI.crs(b) &&
        (parent(a) === parent(b) || parent(a) == parent(b))
end
function Base.isequal(a::GeometryLookup, b::GeometryLookup)
    a === b && return true
    return isequal(a.manifold, b.manifold) && isequal(DD.name(DD.dims(a)), DD.name(DD.dims(b))) && isequal(GI.crs(a), GI.crs(b)) &&
        (parent(a) === parent(b) || isequal(parent(a), parent(b)))
end
Base.hash(l::GeometryLookup, h::UInt) =
    hash(parent(l), hash(GI.crs(l), hash(DD.name(DD.dims(l)), hash(l.manifold, hash(:GeometryLookup, h)))))

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
    l.manifold isa GO.Spherical && print(io, " ", l.manifold)
end
