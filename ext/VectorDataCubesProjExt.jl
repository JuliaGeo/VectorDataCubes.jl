module VectorDataCubesProjExt

using VectorDataCubes
import Proj

function _crsdatum(crs)
    projcrs = try
        convert(Proj.CRS, crs)
    catch err
        throw(ArgumentError("Proj could not parse the CRS $(repr(crs)): $(sprint(showerror, err))"))
    end
    Proj.is_projected(projcrs) && return (; kind=:projected, radius=nothing, sphere=false)
    Proj.is_geographic(projcrs) || throw(ArgumentError(
        "Spherical edges require a geographic CRS; $(repr(crs)) is neither geographic nor projected."
    ))

    ellipsoid = Proj.proj_get_ellipsoid(projcrs)
    ellipsoid == C_NULL && throw(ArgumentError(
        "Proj could not obtain an ellipsoid from the geographic CRS $(repr(crs))."
    ))
    a = Ref{Cdouble}()
    b = Ref{Cdouble}()
    computed = Ref{Cint}()
    invf = Ref{Cdouble}()
    ok = try
        Proj.proj_ellipsoid_get_parameters(ellipsoid, a, b, computed, invf)
    finally
        Proj.proj_destroy(ellipsoid)
    end
    ok == 1 || throw(ArgumentError(
        "Proj could not read ellipsoid parameters from the geographic CRS $(repr(crs))."
    ))
    major, minor = Float64(a[]), Float64(b[])
    valid = isfinite(major) && isfinite(minor) && major > 0 && minor > 0
    valid || throw(ArgumentError("The CRS $(repr(crs)) has invalid ellipsoid axes ($major, $minor)."))

    # A true sphere keeps its declared radius exactly. For an ellipsoid use the
    # IUGG arithmetic mean radius, (2a + b) / 3.
    sphere = major == minor
    radius = sphere ? major : (2major + minor) / 3
    return (; kind=:geographic, radius, sphere)
end

VectorDataCubes._crsdatum(crs::Proj.CRS) = _crsdatum(crs)
VectorDataCubes._crsdatum(crs::AbstractString) = _crsdatum(crs)
VectorDataCubes._crsdatum(crs::Proj.GFT.CoordinateReferenceSystemFormat) = _crsdatum(crs)
VectorDataCubes._crsdatum(crs::Proj.GFT.MixedFormat) = _crsdatum(crs)

end
