module VectorDataCubesProjExt

import Proj
import VectorDataCubes
import VectorDataCubes: GeometryLookup
import GeometryOps as GO
import GeoInterface as GI
import Rasters as RA

VectorDataCubes._reproject(target::RA.GeoFormat, l::GeometryLookup) =
    GO.reproject(parent(l); source_crs=GI.crs(l), target_crs=target, always_xy=true)

end
