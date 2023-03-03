from typing import Dict
from typing import Optional
from typing import Tuple

import dolfin
import mshr

from ..dolfin_utils import Geometry
from ..dolfin_utils import MarkerFunctions
from ._utils import default_markers


def create_lv_mesh(
    N: int = 13,
    a_endo: float = 1.5,
    b_endo: float = 0.5,
    c_endo: float = 0.5,
    a_epi: float = 2.0,
    b_epi: float = 1.0,
    c_epi: float = 1.0,
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    base_x: float = 0.0,
    markers: Optional[Dict[str, Tuple[int, int]]] = None,
) -> Geometry:
    r"""
    Create an lv-ellipsoidal mesh.

    An ellipsoid is given by the equation

    .. math::

        \frac{x^2}{a} + \frac{y^2}{b} + \frac{z^2}{c} = 1

    We create two ellipsoids, one for the endocardium and one
    for the epicardium and subtract them and then cut the base.
    For simplicity we assume that the longitudinal axis is in
    in :math:`x`-direction and as default the base is located
    at the :math:`x=0` plane.
    """
    dolfin.info("Creating LV mesh. This could take some time...")
    # LV
    # The center of the LV ellipsoid
    center = dolfin.Point(*center)

    # Markers
    if markers is None:
        markers = default_markers()
        markers.pop("rv", None)

    class Endo(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return (x[0] - center.x()) ** 2 / a_endo**2 + (
                x[1] - center.y()
            ) ** 2 / b_endo**2 + (
                x[2] - center.z()
            ) ** 2 / c_endo**2 - 1.1 < dolfin.DOLFIN_EPS and on_boundary

    class Base(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return x[0] - base_x < dolfin.DOLFIN_EPS and on_boundary

    class Epi(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return (x[0] - center.x()) ** 2 / a_epi**2 + (
                x[1] - center.y()
            ) ** 2 / b_epi**2 + (
                x[2] - center.z()
            ) ** 2 / c_epi**2 - 0.9 > dolfin.DOLFIN_EPS and on_boundary

    # The plane cutting the base
    diam = -2 * a_epi
    box = mshr.Box(dolfin.Point(base_x, -diam, -diam), dolfin.Point(diam, diam, diam))

    # LV epicardium
    el_lv = mshr.Ellipsoid(center, a_epi, b_epi, c_epi)
    # LV endocardium
    el_lv_endo = mshr.Ellipsoid(center, a_endo, b_endo, c_endo)

    # LV geometry (subtract the smallest ellipsoid)
    lv = el_lv - el_lv_endo

    # LV geometry
    m = lv - box

    # Create mesh
    print("Generate mesh. This can take some time...")
    mesh = mshr.generate_mesh(m, N)

    ffun = dolfin.MeshFunction("size_t", mesh, 2)
    ffun.set_all(0)

    endo = Endo()
    endo.mark(ffun, markers["lv"])

    base = Base()
    base.mark(ffun, markers["base"])

    epi = Epi()
    epi.mark(ffun, markers["epi"])

    marker_functions = MarkerFunctions(ffun=ffun, vfun=None, efun=None, cfun=None)
    return Geometry(mesh=mesh, markers=markers, marker_functions=marker_functions)
