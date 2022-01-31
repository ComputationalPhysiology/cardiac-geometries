from .import_checks import has_dolfin
from .import_checks import has_gmsh
from .import_checks import has_ldrb
from .import_checks import has_meshio
from .import_checks import has_mshr


__all__ = [
    "has_dolfin",
    "has_gmsh",
    "has_ldrb",
    "has_meshio",
    "has_mshr",
]

if has_gmsh():
    from . import _gmsh as gmsh  # noqa: F401

    __all__.append("gmsh")

if has_dolfin():
    from ._dolfin_utils import gmsh2dolfin  # noqa: F401

    __all__.append("gmsh2dolfin")
