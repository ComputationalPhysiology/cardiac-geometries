try:
    import meshio  # noqa: F401

    _has_meshio = True
except ImportError:
    _has_meshio = False


try:
    import dolfin  # noqa: F401

    _has_dolfin = True
except ImportError:
    _has_dolfin = False


try:
    import gmsh  # noqa: F401

    _has_gmsh = True
except (ImportError, OSError):
    _has_gmsh = False

try:
    import mshr  # noqa: F401

    _has_mshr = True
except (ImportError, OSError):
    _has_mshr = False

try:
    import ldrb  # noqa: F401

    _has_ldrb = True
except ImportError:
    _has_ldrb = False


def has_meshio() -> bool:
    return _has_meshio


def has_dolfin() -> bool:
    return _has_dolfin


def has_gmsh() -> bool:
    return _has_gmsh


def has_ldrb() -> bool:
    return _has_ldrb


def has_mshr() -> bool:
    return _has_mshr
