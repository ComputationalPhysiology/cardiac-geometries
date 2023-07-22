import datetime
import json
import math
import tempfile
from importlib.metadata import metadata
from pathlib import Path
from typing import Optional
from typing import Union

from ._import_checks import has_dolfin
from ._import_checks import has_gmsh
from .geometry import Geometry
from .geometry import MeshTypes
from .utils import json_serial

meta = metadata("cardiac_geometries")
__version__ = meta["Version"]


def create_biv_ellipsoid(
    outdir: Union[str, Path, None] = None,
    char_length: float = 0.5,
    center_lv_y: float = 0.0,
    center_lv_z: float = 0.0,
    a_endo_lv: float = 2.5,
    b_endo_lv: float = 1.0,
    c_endo_lv: float = 1.0,
    a_epi_lv: float = 3.0,
    b_epi_lv: float = 1.5,
    c_epi_lv: float = 1.5,
    center_rv_y: float = 0.5,
    center_rv_z: float = 0.0,
    a_endo_rv: float = 3.0,
    b_endo_rv: float = 1.5,
    c_endo_rv: float = 1.5,
    a_epi_rv: float = 4.0,
    b_epi_rv: float = 2.5,
    c_epi_rv: float = 2.0,
    create_fibers: bool = False,
    fiber_angle_endo: float = -60,
    fiber_angle_epi: float = +60,
    fiber_space: str = "P_1",
) -> Optional[Geometry]:
    """Create BiV ellipsoidal geometry

    Parameters
    ----------
    outdir : Union[str, Path, None], optional
        Directory where to save the results. If not provided a temporary
        directory will be created, by default None, by default None
    char_length : float, optional
        Characteristic length of mesh, by default 0.5
    center_lv_y : float, optional
        Y-coordinate for the center of the lv, by default 0.0
    center_lv_z : float, optional
        Z-coordinate for the center of the lv, by default 0.0
    a_endo_lv : float, optional
        Dilation of lv endo ellipsoid in the x-direction, by default 2.5
    b_endo_lv : float, optional
       Dilation of lv endo ellipsoid in the y-direction, by default 1.0
    c_endo_lv : float, optional
       Dilation of lv endo ellipsoid in the z-direction, by default 1.0
    a_epi_lv : float, optional
        Dilation of lv epi ellipsoid in the x-direction, by default 3.0
    b_epi_lv : float, optional
        Dilation of lv epi ellipsoid in the y-direction, by default 1.5
    c_epi_lv : float, optional
        Dilation of lv epi ellipsoid in the z-direction, by default 1.5
    center_rv_y : float, optional
        Y-coordinate for the center of the rv, by default 0.5
    center_rv_z : float, optional
        Z-coordinate for the center of the rv, by default 0.0
    a_endo_rv : float, optional
       Dilation of rv endo ellipsoid in the x-direction, by default 3.0
    b_endo_rv : float, optional
       Dilation of rv endo ellipsoid in the y-direction, by default 1.5
    c_endo_rv : float, optional
       Dilation of rv endo ellipsoid in the z-direction, by default 1.5
    a_epi_rv : float, optional
        Dilation of rv epi ellipsoid in the x-direction, by default 4.0
    b_epi_rv : float, optional
        Dilation of rv epi ellipsoid in the y-direction, by default 2.5
    c_epi_rv : float, optional
        Dilation of rv epi ellipsoid in the z-direction, by default 2.0
    create_fibers : bool, optional
        If True create analytic fibers, by default False
    fiber_angle_endo : float, optional
        Angle for the endocardium, by default -60
    fiber_angle_epi : float, optional
        Angle for the epicardium, by default +60
    fiber_space : str, optional
        Function space for fibers of the form family_degree, by default "P_1"

    Returns
    -------
    Optional[Geometry]
        A Geometry with the mesh, markers, markers functions and fibers.
        Returns None if dolfin is not installed.

    Raises
    ------
    ImportError
        If gmsh is not installed
    """
    if not has_gmsh():
        raise ImportError("Cannot create BiV ellipsoid. Gmsh is not installed")

    _tmpfile = None
    if outdir is None:
        _tmpfile = tempfile.TemporaryDirectory()
        outdir = _tmpfile.__enter__()

    from ._gmsh import biv_ellipsoid

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    with open(outdir / "info.json", "w") as f:
        json.dump(
            {
                "char_length": char_length,
                "center_lv": (0.0, center_lv_y, center_lv_z),
                "a_endo_lv": a_endo_lv,
                "b_endo_lv": b_endo_lv,
                "c_endo_lv": c_endo_lv,
                "a_epi_lv": a_epi_lv,
                "b_epi_lv": b_epi_lv,
                "c_epi_lv": c_epi_lv,
                "center_rv": (0.0, center_rv_y, center_rv_z),
                "a_endo_rv": a_endo_rv,
                "b_endo_rv": b_endo_rv,
                "c_endo_rv": c_endo_rv,
                "a_epi_rv": a_epi_rv,
                "b_epi_rv": b_epi_rv,
                "c_epi_rv": c_epi_rv,
                "create_fibers": create_fibers,
                "fibers_angle_endo": fiber_angle_endo,
                "fibers_angle_epi": fiber_angle_epi,
                "fiber_space": fiber_space,
                "mesh_type": MeshTypes.biv_ellipsoid.value,
                "cardiac_geometry_version": __version__,
                "timestamp": datetime.datetime.now().isoformat(),
            },
            f,
            indent=2,
            default=json_serial,
        )

    mesh_name = outdir / "biv_ellipsoid.msh"

    biv_ellipsoid(
        mesh_name=mesh_name.as_posix(),
        char_length=char_length,
        center_lv=(0.0, center_lv_y, center_lv_z),
        a_endo_lv=a_endo_lv,
        b_endo_lv=b_endo_lv,
        c_endo_lv=c_endo_lv,
        a_epi_lv=a_epi_lv,
        b_epi_lv=b_epi_lv,
        c_epi_lv=c_epi_lv,
        center_rv=(0.0, center_rv_y, center_rv_z),
        a_endo_rv=a_endo_rv,
        b_endo_rv=b_endo_rv,
        c_endo_rv=c_endo_rv,
        a_epi_rv=a_epi_rv,
        b_epi_rv=b_epi_rv,
        c_epi_rv=c_epi_rv,
    )

    if not has_dolfin():
        return None

    from .dolfin_utils import gmsh2dolfin

    geometry = gmsh2dolfin(mesh_name, unlink=False)

    with open(outdir / "markers.json", "w") as f:
        json.dump(geometry.markers, f, default=json_serial)

    if create_fibers:
        from .fibers._biv_ellipsoid import create_biv_fibers

        create_biv_fibers(
            mesh=geometry.mesh,
            ffun=geometry.marker_functions.ffun,
            markers=geometry.markers,
            fiber_space=fiber_space,
            alpha_endo=fiber_angle_endo,
            alpha_epi=fiber_angle_epi,
            outdir=outdir,
        )

    geo = Geometry.from_folder(outdir)

    if _tmpfile is not None:
        _tmpfile.__exit__(None, None, None)

    return geo


def create_biv_ellipsoid_torso(
    outdir: Union[str, Path, None] = None,
    char_length: float = 0.5,
    heart_as_surface: bool = True,
    torso_length: float = 20.0,
    torso_width: float = 20.0,
    torso_height: float = 20.0,
    rotation_angle: float = math.pi / 6,
    center_lv_y: float = 0.0,
    center_lv_z: float = 0.0,
    a_endo_lv: float = 2.5,
    b_endo_lv: float = 1.0,
    c_endo_lv: float = 1.0,
    a_epi_lv: float = 3.0,
    b_epi_lv: float = 1.5,
    c_epi_lv: float = 1.5,
    center_rv_y: float = 0.5,
    center_rv_z: float = 0.0,
    a_endo_rv: float = 3.0,
    b_endo_rv: float = 1.5,
    c_endo_rv: float = 1.5,
    a_epi_rv: float = 4.0,
    b_epi_rv: float = 2.5,
    c_epi_rv: float = 2.0,
    create_fibers: bool = False,
    fiber_angle_endo: float = -60,
    fiber_angle_epi: float = +60,
    fiber_space: str = "P_1",
) -> Optional[Geometry]:
    """Create BiV ellipsoidal geometry

    Parameters
    ----------
    outdir : Union[str, Path, None], optional
        Directory where to save the results. If not provided a temporary
        directory will be created, by default None, by default None
    char_length : float, optional
        Characteristic length of mesh, by default 0.5
    heart_as_surface: bool
        If true, create the heart as a a surface inside the torso,
        otherwise let the heart be a volume, by default True.
    torso_length : float, optional
        Length of torso in the x-direction, by default 20.0
    torso_width : float, optional
        Length of torso in the y-direction, by default 20.0
    torso_height : float, optional
        Length of torso in the z-direction, by default 20.0
    rotation_angle: float, optional
        Angle to rotate the torso in order to object realistic position of
        the heart in a torso, by default pi / 6
    center_lv_y : float, optional
        Y-coordinate for the center of the lv, by default 0.0
    center_lv_z : float, optional
        Z-coordinate for the center of the lv, by default 0.0
    a_endo_lv : float, optional
        Dilation of lv endo ellipsoid in the x-direction, by default 2.5
    b_endo_lv : float, optional
       Dilation of lv endo ellipsoid in the y-direction, by default 1.0
    c_endo_lv : float, optional
       Dilation of lv endo ellipsoid in the z-direction, by default 1.0
    a_epi_lv : float, optional
        Dilation of lv epi ellipsoid in the x-direction, by default 3.0
    b_epi_lv : float, optional
        Dilation of lv epi ellipsoid in the y-direction, by default 1.5
    c_epi_lv : float, optional
        Dilation of lv epi ellipsoid in the z-direction, by default 1.5
    center_rv_y : float, optional
        Y-coordinate for the center of the rv, by default 0.5
    center_rv_z : float, optional
        Z-coordinate for the center of the rv, by default 0.0
    a_endo_rv : float, optional
       Dilation of rv endo ellipsoid in the x-direction, by default 3.0
    b_endo_rv : float, optional
       Dilation of rv endo ellipsoid in the y-direction, by default 1.5
    c_endo_rv : float, optional
       Dilation of rv endo ellipsoid in the z-direction, by default 1.5
    a_epi_rv : float, optional
        Dilation of rv epi ellipsoid in the x-direction, by default 4.0
    b_epi_rv : float, optional
        Dilation of rv epi ellipsoid in the y-direction, by default 2.5
    c_epi_rv : float, optional
        Dilation of rv epi ellipsoid in the z-direction, by default 2.0
    create_fibers : bool, optional
        If True create analytic fibers, by default False
    fiber_angle_endo : float, optional
        Angle for the endocardium, by default -60
    fiber_angle_epi : float, optional
        Angle for the epicardium, by default +60
    fiber_space : str, optional
        Function space for fibers of the form family_degree, by default "P_1"

    Returns
    -------
    Optional[Geometry]
        A Geometry with the mesh, markers, markers functions and fibers.
        Returns None if dolfin is not installed.

    Raises
    ------
    ImportError
        If gmsh is not installed
    """
    if not has_gmsh():
        raise ImportError("Cannot create BiV ellipsoid. Gmsh is not installed")

    _tmpfile = None
    if outdir is None:
        _tmpfile = tempfile.TemporaryDirectory()
        outdir = _tmpfile.__enter__()

    from ._gmsh import biv_ellipsoid_torso

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    with open(outdir / "info.json", "w") as f:
        json.dump(
            {
                "char_length": char_length,
                "heart_as_surface": heart_as_surface,
                "torso_length": torso_length,
                "torso_width": torso_width,
                "torso_height": torso_height,
                "rotation_angle": rotation_angle,
                "center_lv": (0.0, center_lv_y, center_lv_z),
                "a_endo_lv": a_endo_lv,
                "b_endo_lv": b_endo_lv,
                "c_endo_lv": c_endo_lv,
                "a_epi_lv": a_epi_lv,
                "b_epi_lv": b_epi_lv,
                "c_epi_lv": c_epi_lv,
                "center_rv": (0.0, center_rv_y, center_rv_z),
                "a_endo_rv": a_endo_rv,
                "b_endo_rv": b_endo_rv,
                "c_endo_rv": c_endo_rv,
                "a_epi_rv": a_epi_rv,
                "b_epi_rv": b_epi_rv,
                "c_epi_rv": c_epi_rv,
                "create_fibers": create_fibers,
                "fibers_angle_endo": fiber_angle_endo,
                "fibers_angle_epi": fiber_angle_epi,
                "fiber_space": fiber_space,
                "mesh_type": MeshTypes.biv_ellipsoid.value,
                "cardiac_geometry_version": __version__,
                "timestamp": datetime.datetime.now().isoformat(),
            },
            f,
            indent=2,
            default=json_serial,
        )

    mesh_name = outdir / "biv_ellipsoid_torso.msh"

    biv_ellipsoid_torso(
        mesh_name=mesh_name.as_posix(),
        char_length=char_length,
        heart_as_surface=heart_as_surface,
        torso_length=torso_length,
        torso_height=torso_height,
        torso_width=torso_width,
        rotation_angle=rotation_angle,
        center_lv=(0.0, center_lv_y, center_lv_z),
        a_endo_lv=a_endo_lv,
        b_endo_lv=b_endo_lv,
        c_endo_lv=c_endo_lv,
        a_epi_lv=a_epi_lv,
        b_epi_lv=b_epi_lv,
        c_epi_lv=c_epi_lv,
        center_rv=(0.0, center_rv_y, center_rv_z),
        a_endo_rv=a_endo_rv,
        b_endo_rv=b_endo_rv,
        c_endo_rv=c_endo_rv,
        a_epi_rv=a_epi_rv,
        b_epi_rv=b_epi_rv,
        c_epi_rv=c_epi_rv,
    )

    if not has_dolfin():
        return None

    from .dolfin_utils import gmsh2dolfin

    geometry = gmsh2dolfin(mesh_name, unlink=False)

    with open(outdir / "markers.json", "w") as f:
        json.dump(geometry.markers, f, default=json_serial)

    geo = Geometry.from_folder(outdir)

    if _tmpfile is not None:
        _tmpfile.__exit__(None, None, None)

    return geo


def create_lv_ellipsoid(
    outdir: Optional[Union[Path, str]] = None,
    r_short_endo: float = 7.0,
    r_short_epi: float = 10.0,
    r_long_endo: float = 17.0,
    r_long_epi: float = 20.0,
    psize_ref: float = 3,
    mu_apex_endo: float = -math.pi,
    mu_base_endo: float = -math.acos(5 / 17),
    mu_apex_epi: float = -math.pi,
    mu_base_epi: float = -math.acos(5 / 20),
    create_fibers: bool = False,
    fiber_angle_endo: float = -60,
    fiber_angle_epi: float = +60,
    fiber_space: str = "P_1",
) -> Optional[Geometry]:
    """Create an LV ellipsoidal geometry

    Parameters
    ----------
    outdir : Optional[Path], optional
        Directory where to save the results. If not provided a temporary
        directory will be created, by default None
    r_short_endo : float, optional
        Shortest radius on the endocardium layer, by default 7.0
    r_short_epi : float, optional
       Shortest radius on the epicardium layer, by default 10.0
    r_long_endo : float, optional
        Longest radius on the endocardium layer, by default 17.0
    r_long_epi : float, optional
        Longest radius on the epicardium layer, by default 20.0
    psize_ref : float, optional
        The reference point size (smaller values yield as finer mesh, by default 3
    mu_apex_endo : float, optional
        Angle for the endocardial apex, by default -math.pi
    mu_base_endo : float, optional
        Angle for the endocardial base, by default -math.acos(5 / 17)
    mu_apex_epi : float, optional
        Angle for the epicardial apex, by default -math.pi
    mu_base_epi : float, optional
        Angle for the epicardial apex, by default -math.acos(5 / 20)
    create_fibers : bool, optional
        If True create analytic fibers, by default False
    fiber_angle_endo : float, optional
        Angle for the endocardium, by default -60
    fiber_angle_epi : float, optional
        Angle for the epicardium, by default +60
    fiber_space : str, optional
        Function space for fibers of the form family_degree, by default "P_1"

    Returns
    -------
    Optional[Geometry]
        A Geometry with the mesh, markers, markers functions and fibers.
        Returns None if dolfin is not installed.

    Raises
    ------
    ImportError
        If gmsh is not installed
    """

    if not has_gmsh():
        raise ImportError("Cannot create BiV ellipsoid. Gmsh is not installed")

    _tmpfile = None
    if outdir is None:
        _tmpfile = tempfile.TemporaryDirectory()
        outdir = _tmpfile.__enter__()

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    with open(outdir / "info.json", "w") as f:
        json.dump(
            {
                "r_short_endo": r_short_endo,
                "r_short_epi": r_short_epi,
                "r_long_endo": r_long_endo,
                "r_long_epi": r_long_epi,
                "psize_ref": psize_ref,
                "mu_apex_endo": mu_apex_endo,
                "mu_base_endo": mu_base_endo,
                "mu_apex_epi": mu_apex_epi,
                "mu_base_epi": mu_base_epi,
                "create_fibers": create_fibers,
                "fibers_angle_endo": fiber_angle_endo,
                "fibers_angle_epi": fiber_angle_epi,
                "fiber_space": fiber_space,
                "mesh_type": MeshTypes.lv_ellipsoid.value,
                "cardiac_geometry_version": __version__,
                "timestamp": datetime.datetime.now().isoformat(),
            },
            f,
            indent=2,
            default=json_serial,
        )

    from ._gmsh import lv_ellipsoid

    mesh_name = outdir / "lv_ellipsoid.msh"
    lv_ellipsoid(
        mesh_name=mesh_name.as_posix(),
        r_short_endo=r_short_endo,
        r_short_epi=r_short_epi,
        r_long_endo=r_long_endo,
        r_long_epi=r_long_epi,
        mu_base_endo=mu_base_endo,
        mu_base_epi=mu_base_epi,
        mu_apex_endo=mu_apex_endo,
        mu_apex_epi=mu_apex_epi,
        psize_ref=psize_ref,
    )

    if not has_dolfin():
        return None

    from .dolfin_utils import gmsh2dolfin

    geometry = gmsh2dolfin(mesh_name, unlink=False)

    with open(outdir / "markers.json", "w") as f:
        json.dump(geometry.markers, f, default=json_serial)

    if create_fibers:
        from .fibers._lv_ellipsoid import create_microstructure

        create_microstructure(
            mesh=geometry.mesh,
            ffun=geometry.marker_functions.ffun,
            markers=geometry.markers,
            function_space=fiber_space,
            r_short_endo=r_short_endo,
            r_short_epi=r_short_epi,
            r_long_endo=r_long_endo,
            r_long_epi=r_long_epi,
            alpha_endo=fiber_angle_endo,
            alpha_epi=fiber_angle_epi,
            outdir=outdir,
        )

    geo = Geometry.from_folder(outdir)

    if _tmpfile is not None:
        _tmpfile.__exit__(None, None, None)

    return geo


def create_slab(
    outdir: Optional[Union[Path, str]] = None,
    lx: float = 20.0,
    ly: float = 7.0,
    lz: float = 3.0,
    dx: float = 1.0,
    create_fibers: bool = True,
    fiber_angle_endo: float = -60,
    fiber_angle_epi: float = +60,
    fiber_space: str = "P_1",
) -> Optional[Geometry]:
    """Create slab geometry

    Parameters
    ----------
    outdir : Optional[Path], optional
        Directory where to save the results. If not provided a temporary
        directory will be created, by default None
    lx : float, optional
        Length of slab the x-direction, by default 20.0
    ly : float, optional
        Length of slab the x-direction, by default 7.0
    lz : float, optional
        Length of slab the z-direction, by default 3.0
    dx : float, optional
        Element size, by default 1.0
    create_fibers : bool, optional
        If True create analytic fibers, by default True
    fiber_angle_endo : float, optional
        Angle for the endocardium, by default -60
    fiber_angle_epi : float, optional
        Angle for the epicardium, by default +60
    fiber_space : str, optional
        Function space for fibers of the form family_degree, by default "P_1"

    Returns
    -------
    Optional[Geometry]
        A Geometry with the mesh, markers, markers functions and fibers.
        Returns None if dolfin is not installed.

    Raises
    ------
    ImportError
        If gmsh is not installed
    """

    if not has_gmsh():
        raise ImportError("Cannot create BiV ellipsoid. Gmsh is not installed")

    _tmpfile = None
    if outdir is None:
        _tmpfile = tempfile.TemporaryDirectory()
        outdir = _tmpfile.__enter__()

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    with open(outdir / "info.json", "w") as f:
        json.dump(
            {
                "lx": lx,
                "ly": ly,
                "lz": lz,
                "dx": dx,
                "create_fibers": create_fibers,
                "fibers_angle_endo": fiber_angle_endo,
                "fibers_angle_epi": fiber_angle_epi,
                "fiber_space": fiber_space,
                "mesh_type": MeshTypes.slab.value,
                "cardiac_geometry_version": __version__,
                "timestamp": datetime.datetime.now().isoformat(),
            },
            f,
            indent=2,
            default=json_serial,
        )

    from ._gmsh import slab

    mesh_name = outdir / "slab.msh"
    slab(mesh_name=mesh_name.as_posix(), lx=lx, ly=ly, lz=lz, dx=dx)

    if not has_dolfin():
        return None

    from .dolfin_utils import gmsh2dolfin

    geometry = gmsh2dolfin(mesh_name, unlink=False)

    with open(outdir / "markers.json", "w") as f:
        json.dump(geometry.markers, f, default=json_serial)

    if create_fibers:
        from .fibers._slab import create_microstructure

        create_microstructure(
            mesh=geometry.mesh,
            ffun=geometry.marker_functions.ffun,
            markers=geometry.markers,
            function_space=fiber_space,
            alpha_endo=fiber_angle_endo,
            alpha_epi=fiber_angle_epi,
            outdir=outdir,
        )

    geo = Geometry.from_folder(outdir)

    if _tmpfile is not None:
        _tmpfile.__exit__(None, None, None)

    return geo
