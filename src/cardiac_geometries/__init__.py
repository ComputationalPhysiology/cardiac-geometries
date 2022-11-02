from . import calculus
from ._import_checks import has_dolfin
from ._import_checks import has_gmsh
from ._import_checks import has_ldrb
from ._import_checks import has_meshio
from ._import_checks import has_mshr


if has_gmsh():
    from . import _gmsh as gmsh
else:
    gmsh = None  # type: ignore


if has_dolfin():
    from ._dolfin_utils import gmsh2dolfin
    from . import _slab_fibers as slab_fibers
    from . import _lv_ellipsoid_fibers as lv_ellipsoid_fibers
else:
    gmsh2dolfin = None  # type: ignore
    slab_fibers = None  # type: ignore
    lv_ellipsoid_fibers = None  # type:ignore

if has_mshr():
    from . import _mshr as mshr
else:
    mshr = None  # type: ignore

__all__ = [
    "has_dolfin",
    "has_gmsh",
    "has_ldrb",
    "has_meshio",
    "has_mshr",
    "calculus",
    "mshr",
    "gmsh",
    "gmsh2dolfin",
    "slab_fibers",
    "lv_ellipsoid_fibers",
]
