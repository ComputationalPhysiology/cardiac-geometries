from . import calculus
from ._import_checks import has_dolfin
from ._import_checks import has_gmsh
from ._import_checks import has_ldrb
from ._import_checks import has_meshio
from ._import_checks import has_mshr

__all__ = ["has_dolfin", "has_gmsh", "has_ldrb", "has_meshio", "has_mshr", "calculus"]

if has_gmsh():
    from . import _gmsh as gmsh  # noqa: F401

    __all__.append("gmsh")

if has_dolfin():
    from ._dolfin_utils import gmsh2dolfin  # noqa: F401

    __all__.append("gmsh2dolfin")
