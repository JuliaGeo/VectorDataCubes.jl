module VectorDataCubesMakieExt

using VectorDataCubes
using Makie
import Makie.GeometryBasics as GB
import DimensionalData as DD
import GeoInterface as GI

const GeometryDim = DD.Dimension{<:GeometryLookup}
const Plottable = Union{GeometryLookup,GeometryDim}

_geometries(l::GeometryLookup) = parent(l)
_geometries(d::GeometryDim) = parent(DD.lookup(d))

# A `GeometryLookup` is an `AbstractVector` of geometries, so the Makie methods every
# geometry package installs for its own type (`GeoInterface.@enable_makie` on
# `AbstractArray{<:WrapperGeometry}`, and Shapefile's, GeoJSON's, ArchGDAL's and
# LibGEOS's equivalents) match a lookup of those geometries as well. Naming the lookup
# type is more specific than any `AbstractArray{<:SomeGeometry}`, so these three methods
# convert every lookup the same way, whatever its element type.
Makie.convert_arguments(t::Type{<:Makie.Poly}, x::Plottable) =
    Makie.convert_arguments(t, _geometrybasics(x))
# `Lines` needs a method of its own because those packages define
# `convert_arguments(::Type{<:Lines}, …)` next to their `PointBased` one; `Scatter` and
# `LineSegments` reach the `PointBased` method below.
Makie.convert_arguments(t::Type{<:Makie.Lines}, x::Plottable) =
    Makie.convert_arguments(t, _geometrybasics(x))
Makie.convert_arguments(t::Makie.PointBased, x::Plottable) =
    Makie.convert_arguments(t, _geometrybasics(x))

# The recipe `plot(x)` reaches for, taken from the kind of the first geometry.
function Makie.plottype(x::Plottable)
    geoms = _geometries(x)
    isempty(geoms) && throw(ArgumentError(
        "cannot choose a plot type for an empty `GeometryLookup`; call `poly`, `lines` or `scatter` directly"
    ))
    return _plottype(GI.trait(first(geoms)))
end

_plottype(::Union{GI.PolygonTrait,GI.MultiPolygonTrait,GI.LinearRingTrait}) = Makie.Poly
_plottype(::Union{GI.LineStringTrait,GI.MultiLineStringTrait}) = Makie.Lines
_plottype(::Union{GI.PointTrait,GI.MultiPointTrait}) = Makie.Scatter
_plottype(trait) = throw(ArgumentError(
    "no default plot type for geometries with trait $trait; a `GeometryLookup` plots as \
    `poly` (polygons, multipolygons, linear rings), `lines` (linestrings, multilinestrings) \
    or `scatter` (points, multipoints)"
))

"""
    _geometrybasics(x)

A vector of GeometryBasics geometries, one per geometry of `x` (a `GeometryLookup` or a
dimension wrapping one), sharing a single element type so Makie's own
`convert_arguments` methods apply.

A lookup mixing single and multi geometries of one family (polygon/multipolygon,
linestring/multilinestring, point/multipoint) is lifted to the multi kind. Any other
mixture, and an empty lookup, is an `ArgumentError`.
"""
function _geometrybasics(x::Plottable)
    geoms = _geometries(x)
    isempty(geoms) && throw(ArgumentError("cannot plot an empty `GeometryLookup`"))
    kind = GI.trait(first(geoms))
    all(g -> GI.trait(g) === kind, geoms) && return [_geometrybasic(kind, g) for g in geoms]
    multi = _multitrait(kind)
    if isnothing(multi) || any(g -> _multitrait(GI.trait(g)) !== multi, geoms)
        throw(ArgumentError(
            "cannot plot a `GeometryLookup` mixing $(join(unique(map(GI.trait, geoms)), ", ")); \
            only polygon/multipolygon, linestring/multilinestring and point/multipoint \
            mixtures are lifted to one element type"
        ))
    end
    return [_geometrybasic(multi, GI.trait(g), g) for g in geoms]
end

# A linear ring bounds an area, and `_plottype` sends it to `poly`, but `GI.convert`
# gives it a `LineString` that `poly` has no method for.
_geometrybasic(::GI.LinearRingTrait, geom) = GB.Polygon(GI.convert(GB, geom))
_geometrybasic(::GI.AbstractTrait, geom) = GI.convert(GB, geom)

_geometrybasic(::T, kind::T, geom) where {T<:GI.AbstractTrait} = _geometrybasic(kind, geom)
_geometrybasic(::GI.MultiPolygonTrait, ::GI.PolygonTrait, geom) = GB.MultiPolygon([GI.convert(GB, geom)])
_geometrybasic(::GI.MultiLineStringTrait, ::GI.LineStringTrait, geom) = GB.MultiLineString([GI.convert(GB, geom)])
_geometrybasic(::GI.MultiPointTrait, ::GI.PointTrait, geom) = GB.MultiPoint([GI.convert(GB, geom)])

_multitrait(::Union{GI.PolygonTrait,GI.MultiPolygonTrait}) = GI.MultiPolygonTrait()
_multitrait(::Union{GI.LineStringTrait,GI.MultiLineStringTrait}) = GI.MultiLineStringTrait()
_multitrait(::Union{GI.PointTrait,GI.MultiPointTrait}) = GI.MultiPointTrait()
_multitrait(trait) = nothing

end # module
