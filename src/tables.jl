# Tables.jl integration in both directions: [`vectordatacube`](@ref) lifts a table into a
# `DimStack` over `Geometry`, one layer per attribute column; [`vectordatacubetable`](@ref)
# flattens a cube to one row per coordinate, with geometry columns and crs as table metadata.

import Tables
import DataAPI

"""
    VectorDataCubeTable <: Tables.AbstractColumns

The table [`vectordatacubetable`](@ref) returns: a column table wrapping a vector data cube's
`DimensionalData.DimTable`, carrying the geometry columns and crs as DataAPI.jl table metadata
that `GeoInterface.crs`, `DataFrame`, `GeoDataFrames.write` and other metadata readers see.

The metadata keys are GeoInterface's `"GEOINTERFACE:geometrycolumns"` and `"GEOINTERFACE:crs"`.
Geometry columns also carry DataAPI `"edges"` and, when asserted, `"orientation"`
metadata with GeoParquet semantics. These are column metadata, not new GeoInterface keys.
`parent(tbl)` is the wrapped `DimTable`; the Tables.jl columns interface forwards to it.
"""
struct VectorDataCubeTable{T<:DD.DimTable,C,M} <: Tables.AbstractColumns
    table::T
    geometrycolumns::Tuple{Vararg{Symbol}}
    crs::C
    columnmetadata::M
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

DataAPI.colmetadatasupport(::Type{<:VectorDataCubeTable}) = (read=true, write=false)
DataAPI.colmetadatakeys(t::VectorDataCubeTable) =
    Dict(col => DataAPI.colmetadatakeys(t, col) for col in GI.geometrycolumns(t))
DataAPI.colmetadatakeys(t::VectorDataCubeTable, col) = String.(keys(_columnmetadata(t, col)))
function _columnmetadata(t::VectorDataCubeTable, col)
    name = col isa Integer ? Tables.columnnames(t)[col] : Symbol(col)
    name in Tables.columnnames(t) || throw(ArgumentError("Table has no column $col."))
    return get(getfield(t, :columnmetadata), name, (;))
end
function DataAPI.colmetadata(t::VectorDataCubeTable, col, key::AbstractString; style::Bool=false)
    md = _columnmetadata(t, col)
    haskey(md, Symbol(key)) || throw(ArgumentError("Column $col has no metadata key $(repr(key))."))
    value = md[Symbol(key)]
    return style ? (value, :note) : value
end
function DataAPI.colmetadata(t::VectorDataCubeTable, col, key::AbstractString, default; style::Bool=false)
    key in DataAPI.colmetadatakeys(t, col) || return style ? (default, :default) : default
    return DataAPI.colmetadata(t, col, key; style)
end

_inputmanifold(l::GeometryLookup, geometrycolumn, crs=nothing) = l.manifold
function _inputmanifold(table, geometrycolumn, crs=nothing)
    DataAPI.colmetadatasupport(typeof(table)).read || return GO.Planar()
    col = isnothing(geometrycolumn) ? _geometrycolumn(table) : Symbol(geometrycolumn)
    edges = DataAPI.colmetadata(table, col, "edges", "planar")
    orientation = DataAPI.colmetadata(table, col, "orientation", nothing)
    edges in ("planar", "spherical") || throw(ArgumentError("Unsupported edges metadata $(repr(edges)) on column $col."))
    orientation in (nothing, "counterclockwise") || throw(ArgumentError(
        "Unsupported orientation metadata $(repr(orientation)) on column $col."
    ))
    edges == "planar" && return GO.Planar()
    radius = isnothing(crs) ? GO.Spherical().radius : _crsdatum(crs).radius
    return GO.Spherical(; radius, oriented=orientation == "counterclockwise")
end

_crsdatum(crs) = throw(ArgumentError(
    "Resolving the datum of CRS $(repr(crs)) requires Proj.jl. Load Proj before " *
    "importing, exporting, or assigning a CRS to spherical geometries."
))

_validate_manifold_crs(::GO.Planar, crs) = nothing
function _validate_manifold_crs(m::GO.Spherical, crs)
    isnothing(crs) && return nothing
    (; radius, sphere) = _crsdatum(crs)
    # Ellipsoids allow the 5 cm rounding in GeometryOps' default WGS84 mean radius.
    atol = sphere ? 8eps(max(abs(m.radius), abs(radius))) : 0.05
    isapprox(m.radius, radius; rtol=0, atol) || throw(ArgumentError(
        "The spherical radius $(m.radius) disagrees with the radius $radius " *
        "derived from CRS $(repr(crs)). Use a CRS whose datum describes this sphere."
    ))
    return nothing
end

_geometrymetadata(::GO.Planar, crs) = (; edges="planar")
function _geometrymetadata(m::GO.Spherical, crs)
    if isnothing(crs)
        m.radius == GO.Spherical().radius || throw(ArgumentError(
            "A custom spherical radius cannot be preserved in table metadata without a CRS. " *
            "Supply a geographic CRS whose datum describes the sphere."
        ))
    else
        _validate_manifold_crs(m, crs)
    end
    return m.oriented ? (; edges="spherical", orientation="counterclockwise") : (; edges="spherical")
end

"""
    vectordatacube(table; geometrycolumn=nothing, layers=nothing, crs=nokw, manifold=nokw)

Convert a table with a geometry column (a GeoJSON `FeatureCollection`, a `Shapefile.Table`,
a `DataFrame`, ...) to a vector data cube: a `DimStack` over a `Geometry` dimension carrying
a [`GeometryLookup`](@ref) of the geometries, with one layer per remaining column.

Because the attributes are layers over the same `Geometry` dimension,
subsetting the cube (by index or spatial selector) keeps them aligned with the
geometries — there is no separate attribute table to keep in sync.

# Keywords

- `geometrycolumn`: the geometry column, a `Symbol` or a `String`; other geometry-typed
  columns stay ordinary layers. Defaults to the table's `GeoInterface.geometrycolumns`:
  - `:geometry` for most formats, `:Geometry` for a [`vectordatacubetable`](@ref);
  - a table declaring several geometry columns must name one here.
- `layers`: the column names to keep as layers — a `Symbol`, a `String`, or
  any iterable of them. Defaults to every column except the geometry column.
- `crs`: the coordinate reference system of the geometries. Defaults to the
  crs of the table or its geometries, if they carry one.
- `manifold`: an explicit `GeometryOps.Planar()` or `GeometryOps.Spherical()`.
  Defaults to the selected column's DataAPI `"edges"` and `"orientation"` metadata;
  absent `"edges"` means planar, following GeoParquet.

A `missing` geometry is an `ArgumentError` naming the offending rows.

For a cube whose only dimension is `Geometry`, [`vectordatacubetable`](@ref)
is the inverse: `vectordatacube(vectordatacubetable(cube))` recovers the
layers, crs, and edge/orientation semantics.

# Example

The country containing a point, from the 110 m Natural Earth countries:

```julia
using VectorDataCubes, NaturalEarth
countries = vectordatacube(naturalearth("admin_0_countries", 110))
countries[Geometry(Contains((9.0, 50.0)))][:NAME]
```
"""
function vectordatacube(table; geometrycolumn=nothing, layers=nothing, crs=nokw, manifold=nokw)
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
    infer_manifold = isnokw(manifold)
    infer_manifold && (manifold = _inputmanifold(table, geomcol))
    _checkmanifold(manifold)
    isnokw(crs) && (crs = _inputcrs(table, geometries, manifold))
    infer_manifold && (manifold = _inputmanifold(table, geomcol, crs))
    gl = GeometryLookup(geometries; crs, manifold)
    _validate_manifold_crs(manifold, crs)
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

Convert a vector data cube (a `DimArray`/`Raster`/`DimStack`/`RasterStack` with at least one
dimension backed by a [`GeometryLookup`](@ref)) to a [`VectorDataCubeTable`](@ref): a column
table with one row per combination of dimension coordinates, and these columns:

- one per dimension, named after it, holding the actual geometry objects for a geometry
  dimension (`Geometry`, `Dim{:Origin}`, ...);
- one value column per layer.

The geometry column names and the crs travel as DataAPI table metadata, which
`GeoInterface.geometrycolumns`, `GeoInterface.crs`, `DataFrame` and `GeoDataFrames.write` read.
Geometry lookups carrying a crs must agree on it; disagreement is an `ArgumentError`.

Each geometry column carries DataAPI `"edges"` (`"planar"` or `"spherical"`) metadata.
`Spherical(oriented=true)` also emits `"orientation" => "counterclockwise"`, asserting
the supplied ring convention; other lookups omit it. Coordinates remain unchanged.
File writers must explicitly translate these keys to their own format metadata.
A custom spherical radius is exported only when the supplied CRS describes the
same sphere; its physical parameters belong in the CRS.
"""
function vectordatacubetable(cube::Union{DD.AbstractDimArray,DD.AbstractDimStack})
    geomdims = filter(d -> DD.lookup(d) isa GeometryLookup, (DD.dims(cube)..., DD.refdims(cube)...))
    isempty(geomdims) && throw(ArgumentError("""
    `vectordatacubetable` requires a vector data cube with a dimension whose lookup is a
    `GeometryLookup`, but the input has dimensions $(DD.basedims(cube)).
    Wrap your geometries in a `Geometry(GeometryLookup(geoms))` axis first.
    """))
    cols = map(DD.name, geomdims)
    crs = _shared_crs(geomdims)
    columnmetadata = NamedTuple{cols}(map(
        d -> _geometrymetadata(DD.lookup(d).manifold, crs), geomdims
    ))
    return VectorDataCubeTable(DD.DimTable(cube), cols, crs, columnmetadata)
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
