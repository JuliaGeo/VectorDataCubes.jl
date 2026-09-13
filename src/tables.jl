# Tables.jl integration in both directions: [`vectordatacube`](@ref) lifts a
# flat table (or feature collection) into a `DimStack` over a `Geometry`
# dimension, one layer per attribute column; [`vectordatacubetable`](@ref)
# flattens a cube back to a table, one row per geometry × other-dim coordinate,
# with real geometry objects in the geometry columns and the lookup's crs
# carried as DataAPI table metadata.

import Tables
import DataAPI

"""
    VectorDataCubeTable <: Tables.AbstractColumns

The table [`vectordatacubetable`](@ref) returns: a column table wrapping the
`DimensionalData.DimTable` of a vector data cube, which additionally carries
the cube's geometry columns and crs as DataAPI.jl table metadata under
GeoInterface's keys `"GEOINTERFACE:geometrycolumns"` and `"GEOINTERFACE:crs"`.
Any consumer that reads that metadata — `GeoInterface.geometrycolumns`,
`GeoInterface.crs`, `DataFrame`, `GeoDataFrames.write`, ... — therefore sees
the geometry columns and the crs without being told about them.

`parent(tbl)` is the wrapped `DimTable`; the Tables.jl columns interface
forwards to it unchanged.
"""
struct VectorDataCubeTable{T<:DD.DimTable,C} <: Tables.AbstractColumns
    table::T
    geometrycolumns::Tuple{Vararg{Symbol}}
    crs::C
end

Base.parent(t::VectorDataCubeTable) = getfield(t, :table)
GI.geometrycolumns(t::VectorDataCubeTable) = getfield(t, :geometrycolumns)
GI.crs(t::VectorDataCubeTable) = getfield(t, :crs)

Tables.istable(::Type{<:VectorDataCubeTable}) = true
Tables.columnaccess(::Type{<:VectorDataCubeTable}) = true
Tables.columns(t::VectorDataCubeTable) = t
Tables.columnnames(t::VectorDataCubeTable) = Tables.columnnames(parent(t))
Tables.schema(t::VectorDataCubeTable) = Tables.schema(parent(t))
Tables.getcolumn(t::VectorDataCubeTable, i::Int) = Tables.getcolumn(parent(t), i)
Tables.getcolumn(t::VectorDataCubeTable, key::Symbol) = Tables.getcolumn(parent(t), key)
Tables.getcolumn(t::VectorDataCubeTable, ::Type{T}, i::Int, key::Symbol) where {T} =
    Tables.getcolumn(parent(t), T, i, key)
Tables.getcolumn(t::VectorDataCubeTable, dim::Union{DD.Dimension,Type{<:DD.Dimension}}) =
    Tables.getcolumn(parent(t), dim)

function Base.show(io::IO, t::VectorDataCubeTable)
    geomcols = GI.geometrycolumns(t)
    print(io, "VectorDataCubeTable with ", DataAPI.nrow(t), " rows, ", DataAPI.ncol(t),
        length(geomcols) == 1 ? " columns, geometry column " : " columns, geometry columns ",
        join((":$c" for c in geomcols), ", "))
    tablecrs = GI.crs(t)
    isnothing(tablecrs) || print(io, ", crs ", _brief(tablecrs))
    print(io, ", and schema:\n")
    show(IOContext(io, :print_schema_header => false), Tables.schema(t))
end

# A crs is a whole WKT string as often as it is an `EPSG` code.
function _brief(crs; width=60)
    str = sprint(show, crs)
    return length(str) <= width ? str : first(str, width - 1) * "…"
end

DataAPI.metadatasupport(::Type{<:VectorDataCubeTable}) = (read=true, write=false)
function DataAPI.metadatakeys(t::VectorDataCubeTable)
    geomkey = (GI.GEOINTERFACE_GEOMETRYCOLUMNS_KEY,)
    return isnothing(GI.crs(t)) ? geomkey : (geomkey..., GI.GEOINTERFACE_CRS_KEY)
end
function DataAPI.metadata(t::VectorDataCubeTable, key::AbstractString; style::Bool=false)
    valid = DataAPI.metadatakeys(t)
    key in valid || throw(ArgumentError("""
    $(repr(key)) is not table metadata of this `VectorDataCubeTable`, which carries \
    $(join(map(repr, valid), " and ")).
    Pass a default — `DataAPI.metadata(table, key, default)` — to read a key that may be absent.
    """))
    value = key == GI.GEOINTERFACE_CRS_KEY ? GI.crs(t) : GI.geometrycolumns(t)
    return style ? (value, :note) : value
end
function DataAPI.metadata(t::VectorDataCubeTable, key::AbstractString, default; style::Bool=false)
    key in DataAPI.metadatakeys(t) || return style ? (default, :default) : default
    return DataAPI.metadata(t, key; style)
end

"""
    vectordatacube(table; geometrycolumn=nothing, layers=nothing, crs=nokw)

Convert a table with a geometry column (a GeoJSON `FeatureCollection`, a
`Shapefile.Table`, a `DataFrame`, ...) to a vector data cube: a `DimStack`
over a `Geometry` dimension carrying a [`GeometryLookup`](@ref) of the
geometries, with one layer per remaining column.

Because the attributes are layers over the same `Geometry` dimension,
subsetting the cube (by index or spatial selector) keeps them aligned with the
geometries — there is no separate attribute table to keep in sync.

# Keywords

- `geometrycolumn`: the name of the geometry column, a `Symbol` or a `String`.
  Defaults to the table's own metadata (`GeoInterface.geometrycolumns`), which
  is `:geometry` for most formats and `:Geometry` for a
  [`vectordatacubetable`](@ref); a table declaring several geometry columns
  must name one here. Other geometry-typed columns are kept as ordinary layers.
- `layers`: the column names to keep as layers — a `Symbol`, a `String`, or
  any iterable of them. Defaults to every column except the geometry column.
- `crs`: the coordinate reference system of the geometries. Defaults to the
  crs of the table or its geometries, if they carry one.

A `missing` geometry is an `ArgumentError` naming the offending rows.

For a cube whose only dimension is `Geometry`, [`vectordatacubetable`](@ref)
is the inverse: `vectordatacube(vectordatacubetable(cube))` recovers the
layers and the crs.

# Example

The country containing a point, from the 110 m Natural Earth countries:

```@example vectordatacube
using VectorDataCubes, NaturalEarth
countries = vectordatacube(naturalearth("admin_0_countries", 110))
countries[Geometry(Contains((9.0, 50.0)))][:NAME]
```
"""
function vectordatacube(table; geometrycolumn=nothing, layers=nothing, crs=nokw)
    Tables.istable(table) || throw(ArgumentError("""
    `vectordatacube` requires a Tables.jl-compatible table with a geometry column,
    but `Tables.istable` is false for the input ($(typeof(table))).
    To build a cube from a plain geometry vector, use `Geometry(GeometryLookup(geoms))` directly.
    """))
    cols = Tables.columns(table)
    colnames = Tables.columnnames(cols)
    geomcol = isnothing(geometrycolumn) ? _geometrycolumn(table) : Symbol(geometrycolumn)
    geomcol in colnames || throw(ArgumentError("""
    No geometry column :$geomcol found in the table (columns: $(_first_few(colnames; n=10))).
    Pass the right column name with the `geometrycolumn` keyword.
    """))
    geometries = collect(Tables.getcolumn(cols, geomcol))
    if Missing <: eltype(geometries)
        missingrows = findall(ismissing, geometries)
        isempty(missingrows) || throw(ArgumentError("""
        `missing` geometries cannot index a cube, but the geometry column :$geomcol has \
        $(length(missingrows)) of them (rows $(_first_few(missingrows))).
        Filter those rows out first, e.g. with `Tables.subset` or `filter`.
        """))
    end
    if isnokw(crs)
        table_crs = GI.crs(table)
        isnothing(table_crs) || (crs = table_crs)
    end
    gl = GeometryLookup(geometries; crs)
    layernames = _layernames(layers, colnames, geomcol)
    gdim = Geometry(gl)
    return DD.DimStack(NamedTuple{layernames}(map(layernames) do name
        DD.DimArray(collect(Tables.getcolumn(cols, name)), gdim; name)
    end))
end

_first_few(xs; n=5) =
    length(xs) <= n ? join(xs, ", ") : join(first(xs, n), ", ") * ", ..."

# The geometry column detection of `GeometryOpsCore.get_geometries`, except that a table
# declaring several geometry columns (a cube with several geometry dimensions, a
# GeoParquet file, ...) has to be told which one indexes the cube.
function _geometrycolumn(table)
    geomcols = GI.geometrycolumns(table)
    (isnothing(geomcols) || isempty(geomcols)) && return :geometry
    length(geomcols) == 1 || throw(ArgumentError("""
    The table declares $(length(geomcols)) geometry columns \
    ($(join((":$c" for c in geomcols), ", "))), but a cube is indexed by one of them.
    Name it with the `geometrycolumn` keyword; the others stay ordinary layers.
    """))
    return only(geomcols)
end

function _layernames(::Nothing, colnames, geomcol)
    names = Tuple(n for n in colnames if n != geomcol)
    isempty(names) && throw(ArgumentError("""
    The table has no columns other than the geometry column :$geomcol, so there are
    no layers to build. Use `Geometry(GeometryLookup(geoms))` directly instead.
    """))
    return names
end
_layernames(layers::Union{Symbol,AbstractString}, colnames, geomcol) =
    _layernames((layers,), colnames, geomcol)
function _layernames(layers, colnames, geomcol)
    names = map(Symbol, Tuple(layers))
    isempty(names) && throw(ArgumentError("`layers` must name at least one column."))
    unknown = filter(!in(colnames), names)
    isempty(unknown) || throw(ArgumentError("""
    `layers` names columns that are not in the table: $(join((":$n" for n in unknown), ", ")) \
    (columns: $(_first_few(colnames; n=10))).
    """))
    duplicates = unique(filter(n -> count(==(n), names) > 1, names))
    isempty(duplicates) || throw(ArgumentError("""
    `layers` names columns more than once: $(join((":$n" for n in duplicates), ", ")).
    """))
    return names
end

"""
    vectordatacubetable(cube)

Convert a vector data cube (a `DimArray`/`Raster`/`DimStack`/`RasterStack` with
at least one dimension whose lookup is a [`GeometryLookup`](@ref)) to a
[`VectorDataCubeTable`](@ref): a column table with one row per
combination of dimension coordinates — one column per dimension, holding the
actual geometry objects for the geometry dimensions, and one value column per
layer — that also carries the geometry column names and the crs as DataAPI
table metadata, so `GeoInterface.geometrycolumns`, `GeoInterface.crs`,
`DataFrame` and `GeoDataFrames.write` all see them.

Every geometry dimension (`Geometry`, `Dim{:Origin}`, ...) becomes a geometry
column named after the dimension; their lookups must agree on the crs
(lookups without a crs are ignored), otherwise this is an `ArgumentError`.
"""
function vectordatacubetable(cube::Union{DD.AbstractDimArray,DD.AbstractDimStack})
    geomdims = filter(d -> DD.lookup(d) isa GeometryLookup, (DD.dims(cube)..., DD.refdims(cube)...))
    isempty(geomdims) && throw(ArgumentError("""
    `vectordatacubetable` requires a vector data cube with a dimension whose lookup is a
    `GeometryLookup`, but the input has dimensions $(DD.basedims(cube)).
    Wrap your geometries in a `Geometry(GeometryLookup(geoms))` axis first.
    """))
    return VectorDataCubeTable(DD.DimTable(cube), map(DD.name, geomdims), _shared_crs(geomdims))
end

function _shared_crs(geomdims)
    withcrs = filter(d -> !isnothing(GI.crs(DD.lookup(d))), geomdims)
    isempty(withcrs) && return nothing
    crs = GI.crs(DD.lookup(first(withcrs)))
    all(d -> GI.crs(DD.lookup(d)) == crs, withcrs) || throw(ArgumentError("""
    The geometry dimensions of the cube disagree on the crs, so the table cannot carry one:
    $(join(("$(DD.name(d)) => $(GI.crs(DD.lookup(d)))" for d in withcrs), ", ")).
    Reproject the lookups to a common crs (or drop the crs with `setcrs(lookup, nothing)`).
    """))
    return crs
end
