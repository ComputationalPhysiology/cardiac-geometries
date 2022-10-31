import datetime
import json
import math
from importlib.metadata import metadata
from pathlib import Path
from typing import Optional
from typing import Union

import numpy as np
import rich_click as click
from cardiac_geometries.geometry import Geometry

from . import has_dolfin


def json_serial(obj):
    if isinstance(obj, (np.ndarray)):
        return obj.tolist()
    else:
        try:
            return str(obj)
        except Exception:
            raise TypeError("Type %s not serializable" % type(obj))


meta = metadata("cardiac_geometries")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]


@click.group()
@click.version_option(__version__, prog_name="cardiac_geometries")
def app():
    """
    Cardiac Geometries - A library for creating meshes of
    cardiac geometries
    """
    pass


@click.command()
@click.argument(
    "outdir",
    required=True,
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "--r-short-endo",
    default=7.0,
    type=float,
    help="Shortest radius on the endocardium layer",
    show_default=True,
)
@click.option(
    "--r-short-epi",
    default=10.0,
    type=float,
    help="Shortest radius on the epicardium layer",
    show_default=True,
)
@click.option(
    "--r-long-endo",
    default=17.0,
    type=float,
    help="Longest radius on the endocardium layer",
    show_default=True,
)
@click.option(
    "--r-long-epi",
    default=20.0,
    type=float,
    help="Longest radius on the epicardium layer",
    show_default=True,
)
@click.option(
    "--psize-ref",
    default=3.0,
    type=float,
    help="The reference point size (smaller values yield as finer mesh",
    show_default=True,
)
@click.option(
    "--mu-apex-endo",
    default=-math.pi,
    type=float,
    help="Angle for the endocardial apex",
    show_default=True,
)
@click.option(
    "--mu-base-endo",
    default=-math.acos(5 / 17),
    type=float,
    help="Angle for the endocardial base",
    show_default=True,
)
@click.option(
    "--mu-apex-epi",
    default=-math.pi,
    type=float,
    help="Angle for the epicardial apex",
    show_default=True,
)
@click.option(
    "--mu-base-epi",
    default=-math.acos(5 / 20),
    type=float,
    help="Angle for the epicardial base",
    show_default=True,
)
@click.option(
    "--create-fibers",
    default=False,
    is_flag=True,
    type=bool,
    help="If True create analytic fibers",
    show_default=True,
)
@click.option(
    "--fiber-angle-endo",
    default=-60,
    type=float,
    help="Angle for the endocardium",
    show_default=True,
)
@click.option(
    "--fiber-angle-epi",
    default=+60,
    type=float,
    help="Angle for the epicardium",
    show_default=True,
)
@click.option(
    "--fiber-space",
    default="P_1",
    type=str,
    help="Function space for fibers of the form family_degree",
    show_default=True,
)
def create_lv_ellipsoid(
    outdir: Path,
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
):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
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
                "mesh_type": "lv_ellipsoid",
                "cardiac_geometry_version": __version__,
                "timestamp": datetime.datetime.now().isoformat(),
            },
            f,
            indent=2,
            default=json_serial,
        )
    if not has_dolfin():
        return 0

    from ._dolfin_utils import gmsh2dolfin

    geometry = gmsh2dolfin(mesh_name, unlink=False)

    with open(outdir / "markers.json", "w") as f:
        json.dump(geometry.markers, f, default=json_serial)

    if create_fibers:
        from ._lv_ellipsoid_fibers import create_microstructure

        f0, s0, n0 = create_microstructure(
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
        )
        import dolfin

        path = outdir / "microstructure.h5"
        with dolfin.HDF5File(geometry.mesh.mpi_comm(), path.as_posix(), "w") as h5file:
            h5file.write(f0, "f0")
            h5file.write(s0, "s0")
            h5file.write(n0, "n0")

    from .geometry import Geometry

    geo = Geometry.from_folder(outdir)
    geo.save(outdir / "lv_ellipsoid.h5")


@click.command()
@click.argument(
    "outdir",
    required=True,
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "--lx",
    default=20.0,
    type=float,
    help="Length of slab in the x-direction",
    show_default=True,
)
@click.option(
    "--ly",
    default=7.0,
    type=float,
    help="Length of slab in the y-direction",
    show_default=True,
)
@click.option(
    "--lz",
    default=1.0,
    type=float,
    help="Length of slab in the z-direction",
    show_default=True,
)
@click.option(
    "--dx",
    default=1.0,
    type=float,
    help="Element size",
    show_default=True,
)
@click.option(
    "--create-fibers",
    default=False,
    is_flag=True,
    type=bool,
    help="If True create analytic fibers",
    show_default=True,
)
@click.option(
    "--fiber-angle-endo",
    default=-60,
    type=float,
    help="Angle for the endocardium",
    show_default=True,
)
@click.option(
    "--fiber-angle-epi",
    default=+60,
    type=float,
    help="Angle for the epicardium",
    show_default=True,
)
@click.option(
    "--fiber-space",
    default="P_1",
    type=str,
    help="Function space for fibers of the form family_degree",
    show_default=True,
)
def create_slab(
    outdir: Path,
    lx: float = 20.0,
    ly: float = 7.0,
    lz: float = 3.0,
    dx: float = 1.0,
    create_fibers: bool = True,
    fiber_angle_endo: float = -60,
    fiber_angle_epi: float = +60,
    fiber_space: str = "P_1",
):

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    from ._gmsh import slab

    mesh_name = outdir / "slab.msh"
    slab(mesh_name=mesh_name.as_posix(), lx=lx, ly=ly, lz=lz, dx=dx)
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
                "mesh_type": "slab",
                "cardiac_geometry_version": __version__,
                "timestamp": datetime.datetime.now().isoformat(),
            },
            f,
            indent=2,
            default=json_serial,
        )

    if not has_dolfin():
        return 0

    from ._dolfin_utils import gmsh2dolfin

    geometry = gmsh2dolfin(mesh_name, unlink=False)

    with open(outdir / "markers.json", "w") as f:
        json.dump(geometry.markers, f, default=json_serial)

    if create_fibers:
        from ._slab_fibers import create_microstructure

        f0, s0, n0 = create_microstructure(
            mesh=geometry.mesh,
            ffun=geometry.marker_functions.ffun,
            markers=geometry.markers,
            function_space=fiber_space,
            alpha_endo=fiber_angle_endo,
            alpha_epi=fiber_angle_epi,
        )
        import dolfin

        path = outdir / "microstructure.h5"
        with dolfin.HDF5File(geometry.mesh.mpi_comm(), path.as_posix(), "w") as h5file:
            h5file.write(f0, "f0")
            h5file.write(s0, "s0")
            h5file.write(n0, "n0")

    from .geometry import Geometry

    geo = Geometry.from_folder(outdir)
    geo.save(outdir / "slab.h5")


@click.command(
    help=(
        "Convert microstructure.h5 into separate .xdmf-files for f0, s0, and n0. "
        "Assumes a folder containing a 'mesh.xdmf' and 'microstructure.h5'."
    ),
)
@click.argument(
    "folder",
    required=True,
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "--base-direction",
    default="+z",
    type=str,
    help="Direction pointing towards the base. One of +x,-x,+y,-y,+z,-z",
    show_default=True,
)
def fibers_to_xdmf(folder, base_direction: str):

    outdir = Path(folder)
    tetra_mesh_name = outdir / "mesh.xdmf"
    microstructure_path = outdir / "microstructure.h5"

    import dolfin
    from dolfin import FiniteElement  # noqa: F401
    from dolfin import tetrahedron  # noqa: F401
    from dolfin import VectorElement  # noqa: F401
    from .viz import fiber_to_xdmf

    mesh = dolfin.Mesh()

    with dolfin.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(mesh)

    from .geometry import load_microstructure

    f0, s0, n0 = load_microstructure(mesh=mesh, microstructure_path=microstructure_path)

    fiber_to_xdmf(fun=f0, fname=outdir / "f0", base_direction=base_direction)
    fiber_to_xdmf(fun=s0, fname=outdir / "s0", base_direction=base_direction)
    fiber_to_xdmf(fun=n0, fname=outdir / "n0", base_direction=base_direction)


@click.command(help="Convert folder with geometry files to a single file")
@click.argument(
    "folder",
    required=True,
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        writable=True,
        readable=True,
        resolve_path=True,
    ),
)
@click.option("--outfile", type=str, default=None)
def folder2h5(folder: str, outfile: Optional[Union[str, Path]]):
    try:
        import h5py  # noqa: F401
    except ImportError as e:
        msg = "Please install h5py: 'python -m pip install h5py --no-binary=h5py'"
        raise ImportError(msg) from e

    folder_ = Path(folder)
    from .geometry import Geometry

    geo = Geometry.from_folder(folder_)

    if outfile is None:
        outfile = folder_.with_suffix(".h5")
    else:
        outfile = Path(outfile)

    geo.save(outfile)
    print(f"Saved to {outfile}.")


@click.command(help="List information about geometry file")
@click.argument(
    "mesh_path",
    required=True,
    type=click.Path(
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
)
@click.option(
    "--schema-path",
    default=None,
    type=click.Path(
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
)
def info(
    mesh_path: Union[str, Path],
    schema_path: Optional[Union[str, Path]] = None,
) -> None:
    geo = Geometry.from_file(mesh_path, schema_path=schema_path)

    info = getattr(geo, "info", {})
    from rich.console import Console
    from rich.table import Table

    table = Table(title=f"Info about {mesh_path}")

    table.add_column("Key", style="cyan", no_wrap=True)
    table.add_column("Value", justify="right", style="green")
    for k, v in info.items():
        table.add_row(k, str(v))

    console = Console()
    console.print(table)


app.add_command(create_lv_ellipsoid)
app.add_command(create_slab)
app.add_command(fibers_to_xdmf)
app.add_command(folder2h5)
app.add_command(info)
