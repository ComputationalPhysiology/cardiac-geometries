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
    from .dolfin_utils import gmsh2dolfin
    from .fibers import _slab as slab_fibers
    from .fibers import _lv_ellipsoid as lv_ellipsoid_fibers
    from .fibers import _biv_ellipsoid as biv_ellipsoid_fibers
    from . import _mesh as mesh
    from ._mesh import (
        create_biv_ellipsoid,
        create_lv_ellipsoid,
        create_slab,
        create_biv_ellipsoid_torso,
    )
else:
    gmsh2dolfin = None  # type: ignore
    slab_fibers = None  # type: ignore
    lv_ellipsoid_fibers = None  # type:ignore
    biv_ellipsoid_fibers = None  # type: ignore
    create_biv_ellipsoid = None  # type: ignore
    create_lv_ellipsoid = None  # type: ignore
    create_slab = None  # type: ignore
    create_biv_ellipsoid_torso = None  # type: ignore

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
    "biv_ellipsoid_fibers",
    "create_biv_ellipsoid",
    "create_biv_ellipsoid_torso",
    "create_lv_ellipsoid",
    "create_slab",
    "mesh",
]
