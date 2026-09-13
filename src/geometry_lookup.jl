import DimensionalData as DD
import GeometryOps as GO, GeometryOpsCore as GOCore
import GeoInterface as GI
import Rasters as RA
import Extents
import Missings
import SortTileRecursiveTree
import Proj # to load GeometryOps' Proj extension for `reproject`

using Rasters: isnokw, nokw, Lookups, val, Dimensions, dims

"""
    Geometry

A dimension meant to be used with a [`GeometryLookup`](@ref).
"""
DD.@dim Geometry "Geometry"

"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing)

A lookup type for geometry dimensions in vector data cubes.

`GeometryLookup` provides efficient spatial indexing and lookup for 
geometries using an STRtree (Sort-Tile-Recursive tree). 

It is used as the lookup type for geometry dimensions in vector 
data cubes, enabling fast spatial queries and operations.

It spans the dimensions given to it in `dims`, as well as the dimension
 it's wrapped in - you would construct a DimArray with a GeometryLookup
like `DimArray(data, Geometry(GeometryLookup(data, dims)))`.
Here, `Geometry` is a dimension - but selectors in X and Y will also work!

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
struct GeometryLookup{T,A<:AbstractVector{T},D,M<:GO.Manifold,Tree,CRS} <: DD.Dimensions.Lookup{T, 1}
    manifold::M
    data::A
    tree::Tree
    dims::D
    crs::CRS
end
function GeometryLookup(data, dims=(DD.X(), DD.Y()); geometrycolumn=nothing, crs=nokw, tree=nokw)
    # First, retrieve the geometries - from a table, vector of geometries, etc.
    geometries = GOCore.get_geometries(data; geometrycolumn)
    geometries = Missings.disallowmissing(geometries)

    if isnokw(crs)
        crs = GI.crs(data)
        if isnothing(crs)
            crs = GI.crs(first(geometries))
        end
    end
    
    # Check that the geometries are actually geometries
    if any(!GI.isgeometry, geometries)
        throw(ArgumentError("""
        The collection passed in to `GeometryLookup` has some elements that are not geometries 
        (`GI.isgeometry(x) == false` for some `x` in `data`).
        """))
    end
    # Make sure there are only two dimensions
    if length(dims) != 2
        throw(ArgumentError("""
        The `dims` argument to `GeometryLookup` must have two dimensions, but it has $(length(dims)) dimensions (`$(dims)`).
        Please make sure that it has only two dimensions, like `(X(), Y())`.
        """))
    end
    # Build the lookup accelerator tree
    tree = if isnokw(tree)
        SortTileRecursiveTree.STRtree(geometries)
    elseif GO.SpatialTreeInterface.isspatialtree(tree)
        if tree isa Type
            tree(geometries)
        else
            tree
        end
    elseif isnothing(tree)
        nothing
    else
        throw(ArgumentError("""
        Got an argument for `tree` which is not a valid spatial tree (according to `GeometryOps.SpatialTreeInterface`)
        nor `nokw` or `nothing`

        Type is $(typeof(tree))
        """))
    end
    # TODO: auto manifold detection and best tree type for that manifold
    GeometryLookup(GO.Planar(), geometries, tree, dims, crs)
end

GI.crs(l::GeometryLookup) = l.crs
RA.setcrs(l::GeometryLookup, crs) = DD.rebuild(l; crs)

# To reproject a GeometryLookup, you need to reproject the underlying geometry.
function RA.reproject(target::RA.GeoFormat, l::GeometryLookup)
    isnothing(l.crs) && throw(ArgumentError("Cannot reproject a `GeometryLookup` with no crs. Set one first with `setcrs`."))
    new_data = GO.reproject(l.data; source_crs=l.crs, target_crs=target, always_xy=true)
    return DD.rebuild(l; data=new_data, crs=target)
end

#=

## DD methods for the lookup

Here we define DimensionalData's methods for the lookup.
This is broadly standard except for the `rebuild` method, which is used to update the tree accelerator when the data changes.
=#

DD.dims(l::GeometryLookup) = l.dims
# This has to return itself
# DD.dims(d::DD.Dimension{<:GeometryLookup}) = dims(val(d))
DD.order(::GeometryLookup) = Lookups.Unordered()
DD.parent(lookup::GeometryLookup) = lookup.data
# TODO: format for geometry lookup
DD.Dimensions.format(l::GeometryLookup, D::Type, values, axis::AbstractRange) = l

# Make sure that the tree is rebuilt if the data changes
function DD.rebuild(
        lookup::GeometryLookup; 
        data=lookup.data, tree=nokw, 
        dims=lookup.dims, crs=nokw, 
        manifold=nokw, metadata=nokw
    )
    # TODO: metadata support for geometry lookup
    new_tree = if isnokw(tree)
        if data == lookup.data
            lookup.tree
        elseif isempty(data)
            nothing
        else
            SortTileRecursiveTree.STRtree(data)
        end
    elseif GO.SpatialTreeInterface.isspatialtree(tree)
        if tree isa Type
            tree(data)
        else
            tree
        end
    else
        SortTileRecursiveTree.STRtree(data)
    end
    new_crs = if isnokw(crs)
        data_crs = GI.crs(data)
        if isnothing(data_crs)
            lookup.crs
        else
            data_crs
        end
    else
        crs
    end

    new_manifold = isnokw(manifold) ? lookup.manifold : manifold

    return GeometryLookup(new_manifold, Missings.disallowmissing(data), new_tree, dims, new_crs)
end

# # Bounds - get the bounds of the lookup
# function Lookups.bounds(lookup::GeometryLookup)
#     if isempty(lookup.data)
#         Extents.Extent(NamedTuple{DD.name.(lookup.dims)}(ntuple(2) do i; (nothing, nothing); end))
#     else
#         if isnothing(lookup.tree)
#             mapreduce(GI.extent, Extents.union, lookup.data)
#         else
#             GI.extent(lookup.tree)
#         end
#     end
# end

@inline Lookups.reducelookup(l::GeometryLookup) = Lookups.NoLookup(Base.OneTo(1))

function Lookups.show_properties(io::IO, mime, lookup::GeometryLookup)
    print(io, " ")
    show(IOContext(io, :inset => "", :dimcolor => 244), mime, DD.basedims(lookup))
end

# Dimension methods

@inline _reducedims(lookup::GeometryLookup, dim::DD.Dimension) =
    DD.rebuild(dim, [map(x -> zero(x), dim.val[1])])

function DD.format(dim::DD.Dimension{<:GeometryLookup}, axis::AbstractRange)
    # checkaxis(dim, axis)
    return dim
end

