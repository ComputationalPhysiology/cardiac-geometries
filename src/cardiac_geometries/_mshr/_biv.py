from typing import Dict
from typing import Optional
from typing import Tuple

import dolfin
import mshr

from .._dolfin_utils import Geometry
from .._dolfin_utils import MarkerFunctions
from ._utils import default_markers


def create_biv_mesh(
    N: int = 13,
    a_endo_lv: float = 1.5,
    b_endo_lv: float = 0.5,
    c_endo_lv: float = 0.5,
    a_epi_lv: float = 2.0,
    b_epi_lv: float = 1.0,
    c_epi_lv: float = 1.0,
    center_lv: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    a_endo_rv: float = 1.45,
    b_endo_rv: float = 1.25,
    c_endo_rv: float = 0.75,
    a_epi_rv: float = 1.75,
    b_epi_rv: float = 1.5,
    c_epi_rv: float = 1.0,
    center_rv: Tuple[float, float, float] = (0.0, 0.5, 0.0),
    base_x: float = 0.0,
    markers: Optional[Dict[str, Tuple[int, int]]] = None,
) -> Geometry:
    r"""
    Create an biv-ellipsoidal mesh.

    An ellipsoid is given by the equation

    .. math::

        \frac{x^2}{a} + \frac{y^2}{b} + \frac{z^2}{c} = 1

    We create three ellipsoids, one for the LV and RV endocardium
    and one for the epicardium and subtract them and then cut the base.
    For simplicity we assume that the longitudinal axis is in
    in :math:`x`-direction and as default the base is located
    at the :math:`x=0` plane.

    """
    dolfin.info("Creating BiV mesh. This could take some time...")

    # The center of the LV ellipsoid
    center_lv = dolfin.Point(*center_lv)
    # The center of the RV ellipsoid (slightly translated)
    center_rv = dolfin.Point(*center_rv)

    # Markers
    if markers is None:
        markers = default_markers()

    class EndoLV(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return (x[0] - center_lv.x()) ** 2 / a_endo_lv**2 + (
                x[1] - center_lv.y()
            ) ** 2 / b_endo_lv**2 + (
                x[2] - center_lv.z()
            ) ** 2 / c_endo_lv**2 - 1 < dolfin.DOLFIN_EPS and on_boundary

    class Base(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return x[0] - base_x < dolfin.DOLFIN_EPS and on_boundary

    class EndoRV(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return (
                (x[0] - center_rv.x()) ** 2 / a_endo_rv**2
                + (x[1] - center_rv.y()) ** 2 / b_endo_rv**2
                + (x[2] - center_rv.z()) ** 2 / c_endo_rv**2
                - 1
                < dolfin.DOLFIN_EPS
                and (x[0] - center_lv.x()) ** 2 / a_epi_lv**2
                + (x[1] - center_lv.y()) ** 2 / b_epi_lv**2
                + (x[2] - center_lv.z()) ** 2 / c_epi_lv**2
                - 0.9
                > dolfin.DOLFIN_EPS
            ) and on_boundary

    class Epi(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return (
                (x[0] - center_rv.x()) ** 2 / a_epi_rv**2
                + (x[1] - center_rv.y()) ** 2 / b_epi_rv**2
                + (x[2] - center_rv.z()) ** 2 / c_epi_rv**2
                - 0.9
                > dolfin.DOLFIN_EPS
                and (x[0] - center_lv.x()) ** 2 / a_epi_lv**2
                + (x[1] - center_lv.y()) ** 2 / b_epi_lv**2
                + (x[2] - center_lv.z()) ** 2 / c_epi_lv**2
                - 0.9
                > dolfin.DOLFIN_EPS
                and on_boundary
            )

    # The plane cutting the base
    a_epi = max(a_epi_lv, a_epi_rv)
    diam = -5 * a_epi
    box = mshr.Box(dolfin.Point(base_x, a_epi, a_epi), dolfin.Point(diam, diam, diam))
    # Generate mesh

    # LV epicardium
    el_lv = mshr.Ellipsoid(center_lv, a_epi_lv, b_epi_lv, c_epi_lv)
    # LV endocardium
    el_lv_endo = mshr.Ellipsoid(center_lv, a_endo_lv, b_endo_lv, c_endo_lv)

    # LV geometry (subtract the smallest ellipsoid)
    lv = el_lv - el_lv_endo

    # LV epicardium
    el_rv = mshr.Ellipsoid(center_rv, a_epi_rv, b_epi_rv, c_epi_rv)
    # LV endocardium
    el_rv_endo = mshr.Ellipsoid(center_rv, a_endo_rv, b_endo_rv, c_endo_rv)

    # RV geometry (subtract the smallest ellipsoid)
    rv = el_rv - el_rv_endo - el_lv

    # BiV geometry
    m = lv + rv - box

    # Create mesh
    mesh = mshr.generate_mesh(m, N)

    ffun = dolfin.MeshFunction("size_t", mesh, 2)
    ffun.set_all(0)

    endolv = EndoLV()
    endolv.mark(ffun, markers["lv"])
    base = Base()
    base.mark(ffun, markers["base"])
    endorv = EndoRV()
    endorv.mark(ffun, markers["rv"])
    epi = Epi()
    epi.mark(ffun, markers["epi"])

    marker_functions = MarkerFunctions(ffun=ffun, vfun=None, efun=None, cfun=None)
    return Geometry(mesh=mesh, markers=markers, marker_functions=marker_functions)


def mark_biv_mesh(
    mesh: dolfin.Mesh,
    ffun: Optional[dolfin.MeshFunction] = None,
    markers: Optional[Dict[str, int]] = None,
    tol: float = 0.01,
    values: Dict[str, int] = {"lv": 0, "septum": 1, "rv": 2},
) -> dolfin.MeshFunction:

    from .ldrb import scalar_laplacians

    scalars = scalar_laplacians(mesh=mesh, ffun=ffun, markers=markers)

    for cell in dolfin.cells(mesh):

        lv = scalars["lv"](cell.midpoint())
        rv = scalars["rv"](cell.midpoint())
        epi = scalars["epi"](cell.midpoint())

        print(cell.index(), "lv = {}, rv = {}".format(lv, rv))

        if (lv > tol or epi > 1 - tol) and rv < tol:
            print("LV")
            value = values["lv"]
            if lv < tol and rv > lv:
                value = values["rv"]
        elif (rv > tol or epi > 1 - tol) and lv < tol:
            print("RV")
            value = values["rv"]
        else:
            print("SEPTUM")
            value = values["septum"]

        mesh.domains().set_marker((cell.index(), value), 3)

    sfun = dolfin.MeshFunction("size_t", mesh, 3, mesh.domains())
    return sfun
